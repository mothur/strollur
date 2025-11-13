#' @title import_dataset
#' @description
#' The import_dataset function will create a \link{dataset} object from the
#' exported table of a \link{dataset} object.
#' @param table a table containing the data from a \link{dataset} object. You
#' can create the table using 'export(dataset)'.
#' @param tags a vector of strings containing the items you wish to export.
#' Options are 'sequence_data', 'bin_data', 'metadata',
#' 'references', 'sequence_tree', 'sample_tree', and 'reports'.
#' By default, everything is imported.
#' @examples
#'
#' miseq <- miseq_sop_example()
#' data <- import_dataset(export_dataset(miseq))
#' data
#'
#' @return a \link{dataset} object
#' @seealso [dataset$export()]
#' @export
import_dataset <- function(table, tags = NULL) {
  table_version <- attributes(table)$rdataset_version

  # check attributes for valid version
  if (utils::packageVersion("rdataset") != table_version) {
    message <- paste0(
      "[ERROR]: Unable to create 'dataset' object. ",
      "The table was created with rdataset version ",
      table_version, " and you are running rdataset version ",
      utils::packageVersion("rdataset"), "."
    )
    cli::cli_abort(message)
  }

  ht <- !(is.null(tags))

  data <- new_dataset(attributes(table)$dataset_name)

  names <- names(table)
  has_sequence_data <- FALSE

  # extract bin assignment types, ie 'otu', 'asv'
  bin_data_names <- names[grepl("_bin_data", names)]

  # its in the table and you requested it to be imported
  if (("sequence_data" %in% names) && (!ht || ("sequence_data" %in% tags))) {
    has_sequence_data <- TRUE

    # filter sequence_data
    table$sequence_data <- table$sequence_data[
      table$sequence_data$include_sequence,
    ]

    matched_indices <- match(
      table$sequence_abundance_table$sequence_ids,
      table$sequence_data$sequence_ids
    )

    # filter sequence_abundance_table
    table$sequence_abundance_table$sequence_ids <- table$sequence_data$
      sequence_ids[matched_indices]
    table$sequence_abundance_table <- na.omit(table$sequence_abundance_table)
  } else if (("sequence_data" %in% tags) && !("sequence_data" %in% names)) {
    .abort_missing_tag_alert("sequence_data")
  }

  if (length(bin_data_names) != 0) {
    if (!ht || ("bin_data" %in% tags)) {
      # "otu", "asv" or whatever the user specified
      bin_types <- sub("_bin_data", "", bin_data_names)

      for (type in bin_types) {
        bt <- paste0(type, "_bin_data")

        # filter bins from _bin_data
        table[[bt]] <- table[[bt]][table[[bt]]$include_bin, ]

        if (has_sequence_data) {
          b <- paste0(type, "_sequence_bin_assignments")

          # filter seqs from _sequence_bin_assignments
          matched_indices <- match(
            table[[b]]$sequence_ids,
            table$sequence_data$sequence_ids
          )

          table[[b]]$sequence_ids <- table$sequence_data$sequence_ids[
            matched_indices
          ]

          # filter bins from _sequence_bin_assignments
          matched_indices <- match(
            table[[b]]$bin_ids,
            table[[bt]]$bin_ids
          )

          table[[b]]$bin_ids <- table[[bt]]$bin_ids[matched_indices]
          table[[b]] <- na.omit(table[[b]])

          obrs <- paste0(type, "_bin_representative_sequences")
          if (obrs %in% names) {
            # filter seqs from _bin_representative_sequences
            matched_indices <- match(
              table[[obrs]]$sequence_ids,
              table$sequence_data$sequence_ids
            )

            table[[obrs]]$sequence_ids <- table$sequence_data$sequence_ids[
              matched_indices
            ]

            # filter bins from _bin_representative_sequences
            matched_indices <- match(
              table[[obrs]]$bin_ids,
              table[[bt]]$bin_ids
            )

            table[[obrs]]$bin_ids <- table[[bt]]$bin_ids[matched_indices]
            table[[obrs]] <- na.omit(table[[obrs]])
          }
        } else {
          b <- paste0(type, "_bin_abundance_table")

          if (b %in% names) {
            matched_indices <- match(
              table[[b]]$bin_ids,
              table[[bt]]$bin_ids
            )

            # filter bins from _bin_abundance_table
            table[[b]]$bin_ids <- table[[bt]]$bin_ids[matched_indices]
            table[[b]] <- na.omit(table[[b]])
          }
        }
      }
    }
  } else if ("bin_data" %in% tags) {
    .abort_missing_tag_alert("bin_data")
  }

  # look at sequence_data
  if (has_sequence_data) {
    sequence_data_names <- names(table$sequence_data)

    # "sequence_names", "sequences", "comments"
    add_sequences(data, table$sequence_data)

    # "sequence_names", "taxonomies"
    if ("taxonomies" %in% sequence_data_names) {
      assign_sequence_taxonomy(data, table$sequence_data)
    }

    # look at sequence_abundance_table
    if ("sequence_abundance_table" %in% names) {
      # create sequence_names from sequence_ids
      matched_indices <- match(
        table$sequence_abundance_table$sequence_ids,
        table$sequence_data$sequence_ids
      )

      table$sequence_abundance_table$sequence_names <-
        table$sequence_data$sequence_names[matched_indices]


      # "sequence_names", "abundances", "samples", "treatments"
      assign_sequence_abundance(data, table$sequence_abundance_table)
    }
  }


  # if we have bin assignments
  if (length(bin_data_names) != 0) {
    if (!ht || ("bin_data" %in% tags)) {
      # "otu", "asv" or whatever the user specified
      bin_types <- sub("_bin_data", "", bin_data_names)

      for (type in bin_types) {
        bin_type <- paste0(type, "_bin_data")
        bin_data_names <- names(table[[bin_type]])

        # import bin data, label+"_bin_data",
        # with sequence assignments         or  without sequence assignments
        # label+"_sequence_bin_assignments" or label+"_bin_abundance_table"

        sequence_bin_assignments <- paste0(type, "_sequence_bin_assignments")
        otu_bin_abund_table <- paste0(type, "_bin_abundance_table")

        # only abundance data
        if (otu_bin_abund_table %in% names) {
          otu_bin_abund_names <- names(table[[otu_bin_abund_table]])

          # create bin_names from bin_ids
          m_indices <- match(
            table[[otu_bin_abund_table]]$bin_ids,
            table[[bin_type]]$bin_ids
          )

          table[[otu_bin_abund_table]]$bin_names <-
            table[[bin_type]]$bin_names[m_indices]

          # bin_id, abund, sample(optional), treatment(optional)
          assign_bins(data, table[[otu_bin_abund_table]], type = type)

          if ("treatments" %in% otu_bin_abund_names) {
            assign_treatments(data, table[[otu_bin_abund_table]])
          }
        } else if ((sequence_bin_assignments %in% names) && has_sequence_data) {
          # requested sequence_data and has sequence_data

          # create bin_names from bin_ids
          m_indices <- match(
            table[[sequence_bin_assignments]]$bin_ids,
            table[[bin_type]]$bin_ids
          )

          table[[sequence_bin_assignments]]$bin_names <-
            table[[bin_type]]$bin_names[m_indices]

          # create sequence_names from sequence_ids
          m_indices <- match(
            table[[sequence_bin_assignments]]$sequence_ids,
            table$sequence_data$sequence_ids
          )

          table[[sequence_bin_assignments]]$sequence_names <-
            table$sequence_data$sequence_names[m_indices]

          assign_bins(data, table[[sequence_bin_assignments]], type = type)
        } else {
          # does not want sequence data, just bins and abundances
          assign_bins(data, table[[bin_type]], type = type)
        }

        otu_bin_rep_table <- paste0(type, "_bin_representative_sequences")
        if ((otu_bin_rep_table %in% names) && has_sequence_data) {
          # create bin_names from bin_ids
          m_indices <- match(
            table[[otu_bin_rep_table]]$bin_ids,
            table[[bin_type]]$bin_ids
          )

          table[[otu_bin_rep_table]]$bin_names <-
            table[[bin_type]]$bin_names[m_indices]

          # create sequence_names from sequence_ids
          m_indices <- match(
            table[[otu_bin_rep_table]]$sequence_ids,
            table$sequence_data$sequence_ids
          )

          table[[otu_bin_rep_table]]$sequence_names <-
            table$sequence_data$sequence_names[m_indices]

          assign_bin_representative_sequences(
            data,
            table[[otu_bin_rep_table]],
            type = type
          )
        }

        if ("taxonomies" %in% bin_data_names) {
          assign_bin_taxonomy(data, table[[bin_type]], type)
        }
      }
    }
  }

  # add metadata
  if (("metadata" %in% names) && (!ht || ("metadata" %in% tags))) {
    add_metadata(data, table$metadata)
  } else if (("metadata" %in% tags) && !("metadata" %in% names)) {
    .abort_missing_tag_alert("metadata")
  }

  # add references
  if (("references" %in% names) && (!ht || ("references" %in% tags))) {
    add_references(data, table$references)
  } else if (("references" %in% tags) && !("references" %in% names)) {
    .abort_missing_tag_alert("references")
  }

  if (("reports" %in% names) && (!ht || ("reports" %in% tags))) {
    report_names <- names(table$reports)

    for (name in report_names) {
      name_col <- attr(table$reports[[name]], "sequence_name")
      add_report(data, table$reports[[name]], name, name_col)
    }
  } else if (("reports" %in% tags) && !("reports" %in% names)) {
    .abort_missing_tag_alert("reports")
  }

  # add sequence_tree
  if (("sequence_tree" %in% names) && (!ht || ("sequence_tree" %in% tags))) {
    data$add_sequence_tree(table$sequence_tree)
  } else if (("sequence_tree" %in% tags) && !("sequence_tree" %in% names)) {
    .abort_missing_tag_alert("sequence_tree")
  }

  # add sample_tree
  if (("sample_tree" %in% names) && (!ht || ("sample_tree" %in% tags))) {
    data$add_sample_tree(table$sample_tree)
  } else if (("sample_tree" %in% tags) && !("sample_tree" %in% names)) {
    .abort_missing_tag_alert("sample_tree")
  }

  data
}
