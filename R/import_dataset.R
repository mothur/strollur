#' @title Import strollur object from exported data.frame.
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
#' @seealso [export_dataset()]
#' @export
import_dataset <- function(table) {
  table_version <- attributes(table)$strollur_version
  current_version <- new_dataset()$get_version()

  # check attributes for valid version
  if (current_version != table_version) {
    if (.import_compatible(table_version)) {
      message <- paste0(
        "The table was created with strollur version ", table_version,
        ", which is compatible with the current strollur version: ",
        current_version, ". Converting and importing."
      )
      cli::cli_alert(message)
    } else {
      message <- paste0(
        "Unable to create 'strollur' object. ",
        "The table was created with strollur version ", table_version,
        ", which is not compatible with the current strollur version: ",
        current_version, "."
      )
      cli::cli_abort(message)
    }
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
      table$sequence_abundance_table$sequence_id,
      table$sequence_data$sequence_id
    )

    # filter sequence_abundance_table
    table$sequence_abundance_table$sequence_id <- table$sequence_data$
      sequence_id[matched_indices]
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
        b <- paste0(type, "_sequence_bin_assignment")

        # filter seqs from _sequence_bin_assignments
        matched_indices <- match(
          table[[b]]$sequence_id,
          table$sequence_data$sequence_id
        )

        table[[b]]$sequence_id <- table$sequence_data$sequence_id[
          matched_indices
        ]

        # filter bins from _sequence_bin_assignments
        matched_indices <- match(
          table[[b]]$bin_id,
          table[[bt]]$bin_id
        )

        table[[b]]$bin_id <- table[[bt]]$bin_id[matched_indices]
        table[[b]] <- na.omit(table[[b]])

        obrs <- paste0(type, "_bin_representative_sequence")
        if (obrs %in% names) {
          # filter seqs from _bin_representative_sequences
          matched_indices <- match(
            table[[obrs]]$sequence_id,
            table$sequence_data$sequence_id
          )

          table[[obrs]]$sequence_id <- table$sequence_data$sequence_id[
            matched_indices
          ]

          # filter bins from _bin_representative_sequences
          matched_indices <- match(
            table[[obrs]]$bin_id,
            table[[bt]]$bin_id
          )

          table[[obrs]]$bin_id <- table[[bt]]$bin_id[matched_indices]
          table[[obrs]] <- na.omit(table[[obrs]])
        }
      }
    }
  }

  # look at sequence_data
  if (has_sequence_data) {
    sequence_data_names <- names(table$sequence_data)

    # "sequence_name", "sequence", "comment"
    add(data = data, table = table$sequence_data, type = "sequence")

    # "sequence_name", "taxonomy"
    if ("taxonomy" %in% sequence_data_names) {
      xdev_assign_sequence_taxonomy(
        data = data, table = table$sequence_data
      )
    }

    # look at sequence_abundance_table
    if ("sequence_abundance_table" %in% names) {
      # create sequence_names from sequence_ids
      matched_indices <- match(
        table$sequence_abundance_table$sequence_id,
        table$sequence_data$sequence_id
      )

      table$sequence_abundance_table$sequence_name <-
        table$sequence_data$sequence_name[matched_indices]


      # "sequence_name", "abundance", "sample", "treatment"
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
      # label+"_sequence_bin_assignment" or label+"_bin_abundance_table"

      sequence_bin_assignment <- paste0(type, "_sequence_bin_assignment")

      if ((sequence_bin_assignment %in% names) && has_sequence_data) {
        # requested sequence_data and has sequence_data

        # create bin_names from bin_ids
        m_indices <- match(
          table[[sequence_bin_assignment]]$bin_id,
          table[[bin_type]]$bin_id
        )

        table[[sequence_bin_assignment]]$bin_name <-
          table[[bin_type]]$bin_name[m_indices]

        # create sequence_names from sequence_ids
        m_indices <- match(
          table[[sequence_bin_assignment]]$sequence_id,
          table$sequence_data$sequence_id
        )

        table[[sequence_bin_assignment]]$sequence_name <-
          table$sequence_data$sequence_name[m_indices]

        xdev_assign_bins(
          data = data, table = table[[sequence_bin_assignment]],
          bin_type = type
        )
      }

      otu_bin_rep_table <- paste0(type, "_bin_representative_sequence")
      if ((otu_bin_rep_table %in% names) && has_sequence_data) {
        # create bin_names from bin_ids
        m_indices <- match(
          table[[otu_bin_rep_table]]$bin_id,
          table[[bin_type]]$bin_id
        )

        table[[otu_bin_rep_table]]$bin_name <-
          table[[bin_type]]$bin_name[m_indices]

        # create sequence_names from sequence_ids
        m_indices <- match(
          table[[otu_bin_rep_table]]$sequence_id,
          table$sequence_data$sequence_id
        )

        table[[otu_bin_rep_table]]$sequence_name <-
          table$sequence_data$sequence_name[m_indices]

        xdev_assign_bin_representative_sequences(
          data,
          table[[otu_bin_rep_table]],
          bin_type = type
        )
      }

      if ("taxonomy" %in% bin_data_names) {
        xdev_assign_bin_taxonomy(
          data = data, table = table[[bin_type]], bin_type = type
        )
      }
    }
  }

  # add references
  if ("resource_reference" %in% names) {
    xdev_add_references(data = data, table = table$resource_reference)
  }

  # list all report names
  table_names <- c(
      "sequence_data", "sequence_report",
      "sequence_abundance_table", "otu_bin_data",
      "otu_sequence_bin_assignment", "otu_bin_representative_sequence",
      "asv_bin_data",
      "asv_sequence_bin_assignment", "phylotype_bin_data",
      "phylotype_sequence_bin_assignment", "resource_reference",
      "sequence_tree", "sample_tree"
  )
  report_names <- names[!names %in% table_names]

  if (length(report_names) != 0) {
    for (name in report_names) {
      name_col <- attr(table[[name]], "sequence_name")
      if (is.null(name_col)) {
          name_col <- "none"
      }
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
