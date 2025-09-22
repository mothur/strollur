#' @title import
#' @description
#' The import function will create a 'sequence_data' object from the exported
#' table of a 'sequence_data' object.
#' @param table a table containing the data from a 'sequence_data' object. You
#' can create the table using 'sequence_data$export()'.
#' @examples
#'
#' miseq <- miseq_sop_example()
#' dataset <- import(miseq$export())
#' dataset
#'
#' @return A 'sequence_data' object
#' @seealso [sequence_data$export()]
#' @export
import <- function(table) {
  table_version <- attributes(table)$rdataset_version

  # check attributes for valid version
  if (utils::packageVersion("rdataset") != table_version) {
    message <- paste0(
      "[ERROR]: Unable to create 'sequence_data' object. ",
      "The table was created with rdataset version ",
      attributes(table)$rdataset_version,
      " and you are running rdataset version ",
      utils::packageVersion("rdataset"), "."
    )
    cli::abort(message)
  }

  dataset <- sequence_data$new(name = attributes(table)$dataset_name)

  table_names <- names(table)
  has_sequence_data <- FALSE

  # extract bin assigment types
  bin_data_names <- table_names[grepl("_bin_data", table_names)]

  if ("sequence_data" %in% table_names) {
    has_sequence_data <- TRUE

    # filter sequence_data
    table$sequence_data <- table$sequence_data[
      table$sequence_data$include_sequence,
    ]

    matched_indices <- match(
      table$sequence_report$sequence_ids,
      table$sequence_data$sequence_ids
    )

    # filter sequence_report
    table$sequence_report$sequence_ids <- table$sequence_data$sequence_ids[
      matched_indices
    ]
    table$sequence_report <- na.omit(table$sequence_report)

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
      } else {
        b <- paste0(type, "_bin_abundance_table")

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


  # look at sequence_data
  if ("sequence_data" %in% table_names) {
    has_sequence_data <- TRUE
    sequence_data_names <- names(table$sequence_data)

    # "sequence_names", "sequences", "comments"
    dataset$add_sequences(table$sequence_data)

    # "sequence_names", "taxonomies"
    if ("taxonomies" %in% sequence_data_names) {
      dataset$assign_sequence_taxonomy(table$sequence_data)
    }

    # look at sequence_abundance_table
    if ("sequence_abundance_table" %in% table_names) {
      # create sequence_names from sequence_ids
      matched_indices <- match(
        table$sequence_abundance_table$sequence_ids,
        table$sequence_data$sequence_ids
      )

      sequence_names <- table$sequence_data$sequence_names[matched_indices]

      # "sequence_names", "abundances", "samples", "treatments"
      dataset$assign_sequence_abundance(
        data = NULL, sequence_names,
        table$sequence_abundance_table$abundances,
        table$sequence_abundance_table$samples,
        table$sequence_abundance_table$treatments
      )
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

      if (has_sequence_data) {
        sequence_bin_assignments <- paste0(
          type,
          "_sequence_bin_assignments"
        )

        # create bin_names from bin_ids
        m_indices <- match(
          table[[sequence_bin_assignments]]$bin_ids,
          table[[bin_type]]$bin_ids
        )

        bin_names <- table[[bin_type]]$bin_names[m_indices]

        # create sequence_names from sequence_ids
        m_indices <- match(
          table[[sequence_bin_assignments]]$sequence_ids,
          table$sequence_data$sequence_ids
        )

        sequence_names <- table$sequence_data$sequence_names[m_indices]

        dataset$assign_bins(
          bin_names = bin_names,
          sequence_names = sequence_names,
          type = type
        )
      } else {
        otu_bin_abund_table <- paste0(type, "_bin_abundance_table")

        otu_bin_abund_names <- names(table[[otu_bin_abund_table]])

        # create bin_names from bin_ids
        m_indices <- match(
          table[[otu_bin_abund_table]]$bin_ids,
          table[[bin_type]]$bin_ids
        )

        bin_names <- table[[bin_type]]$bin_names[m_indices]

        # bin_id, abund, sample(optional), treatment(optional)
        if ("samples" %in% otu_bin_abund_names) {
          dataset$assign_bins(
            bin_names = bin_names,
            abundances = table[[otu_bin_abund_table]]$abundances,
            samples = table[[otu_bin_abund_table]]$samples,
            type = type
          )
        } else {
          dataset$assign_bins(
            bin_names = bin_names,
            abundances = table[[otu_bin_abund_table]]$abundances,
            type = type
          )
        }

        if ("treatments" %in% otu_bin_abund_names) {
          dataset$assign_treatments(table[[otu_bin_abund_table]])
        }
      }

      if ("taxonomies" %in% bin_data_names) {
        dataset$assign_bin_taxonomy(table[[bin_type]], type = type)
      }
    }
  }

  # add metadata
  if ("metadata" %in% table_names) {
    dataset$add_metadata(table$metadata)
  }

  # add references
  if ("references" %in% table_names) {
    dataset$add_references(table$references)
  }

  # add alignment report
  if ("alignment_report" %in% table_names) {
    dataset$add_alignment_report(
      table$alignment_report,
      attributes(table$alignment_report)$sequence_name_column
    )
  }

  # add contigs report
  if ("contigs_assembly_report" %in% table_names) {
    dataset$add_contigs_assembly_report(
      table$contigs_assembly_report,
      attributes(table$contigs_assembly_report)$sequence_name_column
    )
  }

  # add sequence_tree
  if ("sequence_tree" %in% table_names) {
    dataset$add_sequence_tree(table$sequence_tree)
  }

  # add sample_tree
  if ("sample_tree" %in% table_names) {
    dataset$add_sample_tree(table$sample_tree)
  }

  dataset
}
