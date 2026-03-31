#' @importFrom stats reshape
#' @title read_phyloseq
#' @description
#' The `read_phyloseq()` function reads phyloseq objects created from
#' the phyloseq package
#' (https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
#' and converts it into a strollur object.
#' @param phyloseq_object the phyloseq object that is returned when using
#' any read function in the phyloseq package. It has to be of type "phyloseq"
#' @param treatment_column_name the column name inside your phyloseq object
#' within your sample data that is used to descrbe treatments. It must
#' be a character. Defaults to NULL.
#' @param dataset_name A string containing a name for your dataset.
#' @return a strollur object.
#' @examples
#' miseq <- miseq_sop_example()
#' phylo_obj <- write_phyloseq(miseq)
#' miseq_re_read <- read_phyloseq(phylo_obj)
#' @export
read_phyloseq <- function(phyloseq_object, treatment_column_name = NULL,
                          dataset_name = "") {
  if (!require_namespace("phyloseq")) {
    stop("To use this functionality you have to install the phyloseq package.")
  }

  if (!inherits(phyloseq_object, "phyloseq")) {
    stop("phyloseq_object has to an object created using the phyloseq package.")
  }

  rdaset_object <- strollur$new(dataset_name)

  # samples
  if (!is.null(phyloseq::get_sample(phyloseq_object))) {
    sample_df <- phyloseq::get_sample(phyloseq_object)
    xdev_add_sequences(rdaset_object, data.frame(
      sequence_names =
        unique(rownames(sample_df))
    ))
  }

  # otu table
  if (!is.null(phyloseq_object@otu_table)) {
    shaped <- phyloseq::otu_table(phyloseq_object)@.Data
    names <- rownames(shaped)
    shaped <- cbind(names, shaped)
    colnames(shaped)[1] <- "sequence_names"
    shared_table <- reshape(
      data.frame(shaped),
      varying = as.integer(2:ncol(shaped)),
      v.names = "abundances",
      direction = "long",
      times = colnames(shaped)[-1],
      timevar = "samples",
    )
    shared_table$abundances <- as.numeric(shared_table$abundances)
    xdev_assign_sequence_abundance(
      rdaset_object,
      shared_table[, c("sequence_names", "abundances", "samples")]
    )
  }


  # taxonomy table
  if (!is.null(phyloseq_object@tax_table)) {
    tax_table_phylo <- phyloseq::tax_table(phyloseq_object)@.Data
    taxas <- apply(tax_table_phylo, 1, paste, collapse = ";")
    taxas <- data.frame(sequence_names = names(taxas), taxonomies = taxas)
    rownames(taxas) <- NULL
    xdev_assign_sequence_taxonomy(rdaset_object, taxas)
  }

  # sample data
  if (!is.null(phyloseq_object@sam_data)) {
    df <- data.frame(phyloseq::sample_data(phyloseq_object)@.Data)
    df <- data.frame(apply(df, 2, as.character))
    colnames(df) <- phyloseq::sample_data(phyloseq_object)@names
    df <- cbind(sample = phyloseq::sample_names(phyloseq_object), df)
    if (!is.null(treatment_column_name)) {
      df <- df[, which(colnames(df) != treatment_column_name)]
    }
    add(
      data = rdaset_object,
      table = df,
      type = "metadata"
    )
  }

  # phy tree
  if (!is.null(phyloseq_object@phy_tree)) {
    rdaset_object$add_sequence_tree(phyloseq::phy_tree(phyloseq_object))
  }

  # Treatments
  if (!is.null(treatment_column_name)) {
    tcn <- treatment_column_name
    if (any(phyloseq::sample_variables(phyloseq_object) == tcn)) {
      treatment_df <- data.frame(
        cbind(
          phyloseq::sample_names(phyloseq_object),
          as.character(
            phyloseq::get_variable(phyloseq_object)[[treatment_column_name]]
          )
        )
      )
      colnames(treatment_df) <- c("samples", "treatments")
      assign(data = rdaset_object, table = treatment_df, type = "treatments")
    }
  }
  rdaset_object
}
