#' @title write_phyloseq
#' @description
#' The `write_phyloseq()` function will take any strollur object and
#' return it as a "phyloseq" object.
#' @param data the strollur object you created using one of the many
#' read functions in this package.
#' @return returns a "phyloseq" object.
#' @examples
#' miseq <- miseq_sop_example()
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'   phylo_obj <- write_phyloseq(miseq)
#' } else {
#'   message(paste(
#'     "To use this functionality you have to install the",
#'     "phyloseq package."
#'   ))
#' }
#'
#' @export
write_phyloseq <- function(data) {
  if (!require_namespace("phyloseq")) {
    stop("To use this functionality you have to install the phyloseq package.")
  }

  if (!inherits(data, "strollur")) {
    stop("The data parameter must be an object of type `strollur`.")
  }
  phyloseq_parameter_list <- vector("list", 4)
  if (nrow(abundance(data = data, type = "sequence")) > 0) {
    abundances <- abundance(data = data, type = "sequence", by_sample = TRUE)
    treatments <- NULL
    if (any(colnames(abundances) == "treatment")) {
      treatments <- abundances$treatment
      abundances <- abundances[, -which(colnames(abundances) == "treatment")]
    }
    sequence_names <- names(data = data, type = "sample")
    abundances <- reshape(abundances,
      direction = "wide",
      timevar = "sample",
      idvar = "sequence_name",
    )
    abundances[is.na(abundances)] <- 0
    colnames(abundances) <- sub("[^.]+\\.", "", colnames(abundances))
    otu_table <- as.matrix(abundances[, 2:ncol(abundances)])
    rownames(otu_table) <- abundances$sequence_name
    phyloseq_parameter_list[[1]] <- phyloseq::otu_table(otu_table, TRUE)
  }


  # taxonomies
  if (nrow(report(data = data, type = "sequence_taxonomy")) > 0) {
    df <- report(data = data, type = "sequence_taxonomy")

    if (any(colnames(df) == "confidence")) {
      df$taxonomy <- paste0(df$taxonomy, "(", df$confidence, ")")
      df$confidence <- NULL
    }
    taxas <-
      reshape(df,
        direction = "wide",
        timevar = "level",
        idvar = "sequence_name",
      )

    colnames(taxas) <- c("id", paste("level_", seq(1, ncol(taxas) - 1),
      sep = ""
    ))
    sequence_names <- taxas$id
    taxas <- as.matrix(taxas[, colnames(taxas) != "id"])
    taxas[which(taxas == "NA")] <- NA
    rownames(taxas) <- sequence_names
    phyloseq_parameter_list[[2]] <- phyloseq::tax_table(taxas)
  }

  if (!is.null(data$get_sequence_tree())) {
    phyloseq_parameter_list[[3]] <-
      phyloseq::phy_tree(data$get_sequence_tree())
  }


  metadata <- report(data, "metadata")
  if (!is.null(metadata) && nrow(metadata) > 0) {
    sample_assignments <- report(data, "sample_assignment")
    if (!is.null(sample_assignments) && nrow(sample_assignments) > 0) {
      colnames(sample_assignments)[1] <- "sample"
      metadata <- merge(metadata, sample_assignments, by = "sample")
    }
    rownames(metadata) <- metadata$sample
    metadata$sample <- NULL
    phyloseq_parameter_list[[4]] <- phyloseq::sample_data(metadata)
  }


  indexes <- which(!sapply(phyloseq_parameter_list, is.null))
  if (length(indexes) <= 0) {
    stop("You have an empty object that cannot become a phyloseq object.")
  }

  do.call(
    phyloseq::phyloseq,
    phyloseq_parameter_list[indexes]
  )
}
