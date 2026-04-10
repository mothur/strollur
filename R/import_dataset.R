#' @title import_dataset
#' @description The import_dataset function will create a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#'   from the exported table of a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur}
#'   object.
#'
#' @param table a table containing the data from a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur}
#'   object. You can create the table using 'export(data)'.
#'
#' @examples
#'
#' miseq <- miseq_sop_example()
#' data <- import_dataset(export_dataset(miseq))
#' data
#'
#' @return a
#'   \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
#' @seealso [dataset$export()]
#' @export
import_dataset <- function(table) {
  table_version <- attributes(table)$strollur_version

  # check attributes for valid version
  if (utils::packageVersion("strollur") != table_version) {
    message <- paste0(
      "[ERROR]: Unable to create 'strollur' object. ",
      "The table was created with strollur version ",
      table_version, " and you are running strollur version ",
      utils::packageVersion("strollur"), "."
    )
    cli::cli_abort(message)
  }

  data <- new_dataset(attributes(table)$dataset_name)

  names <- names(table)
  has_sequence_data <- FALSE

  # extract bin assignment types, ie 'otu', 'asv'
  bin_data_names <- names[grepl("_bin_data", names)]

  # its in the table
  if ("sequence_data" %in% names) {
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
  }

  if (length(bin_data_names) != 0) {
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
      }
    }
  }

  # look at sequence_data
  if (has_sequence_data) {
    sequence_data_names <- names(table$sequence_data)

    # "sequence_names", "sequences", "comments"
    add(data = data, table = table$sequence_data, type = "sequences")

    # "sequence_names", "taxonomies"
    if ("taxonomies" %in% sequence_data_names) {
      xdev_assign_sequence_taxonomy(
        data = data, table = table$sequence_data
      )
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
      xdev_assign_sequence_abundance(data, table$sequence_abundance_table)
    }
  }

  # if we have bin assignments
  if (length(bin_data_names) != 0) {
    # "otu", "asv" or whatever the user specified
    bin_types <- sub("_bin_data", "", bin_data_names)

    for (type in bin_types) {
      bin_type <- paste0(type, "_bin_data")
      bin_data_names <- names(table[[bin_type]])

      # import bin data, label+"_bin_data",
      # with sequence assignments         or  without sequence assignments
      # label+"_sequence_bin_assignments" or label+"_bin_abundance_table"

      sequence_bin_assignments <- paste0(type, "_sequence_bin_assignments")

      if ((sequence_bin_assignments %in% names) && has_sequence_data) {
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

        xdev_assign_bins(
          data = data, table = table[[sequence_bin_assignments]],
          bin_type = type
        )
      } else {
        # does not want sequence data, just bins and abundances
        xdev_assign_bins(
          data = data, table = table[[bin_type]],
          bin_type = type
        )
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

        xdev_assign_bin_representative_sequences(
          data,
          table[[otu_bin_rep_table]],
          bin_type = type
        )
      }

      if ("taxonomies" %in% bin_data_names) {
        xdev_assign_bin_taxonomy(
          data = data, table = table[[bin_type]], bin_type = type
        )
      }
    }
  }

  # add metadata
  if ("metadata" %in% names) {
    add(data = data, table = table$metadata, type = "metadata")
  }

  # add references
  if ("references" %in% names) {
    add(data = data, table = table$references, type = "references")
  }

  # list all report names
  report_names <- names[grepl("_report", names)]
  report_names <- report_names[!report_names %in% "sequence_report"]

  if (length(report_names) != 0) {
    for (name in report_names) {
      name_col <- attr(table[[name]], "sequence_name")
      xdev_add_report(
        data = data, table = table[[name]],
        type = name, sequence_name = name_col
      )
    }
  }

  # add sequence_tree
  if ("sequence_tree" %in% names) {
    data$add_sequence_tree(table$sequence_tree)
  }

  # add sample_tree
  if ("sample_tree" %in% names) {
    data$add_sample_tree(table$sample_tree)
  }

  data
}
