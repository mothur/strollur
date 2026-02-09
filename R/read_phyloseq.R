read_phyloseq <- function(phyloseq_object, treatment_column_name = NULL, dataset_name = "") {
  if(!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("To use this functionality you have to install the phyloseq package")
  }
  
  rdaset_object <- dataset$new("dataset_name")

  #samples
  if(!is.null(phyloseq::get_sample(phyloseq_object))) {
    sample_df <- phyloseq::get_sample(phyloseq_object)
    xdev_add_sequences(rdaset_object, data.frame(
    sequence_names =
      unique(rownames(sample_df))
    ))
  }

  #otu table
  if(!is.null(phyloseq_object@otu_table)) {
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
        shared_table[ ,c("sequence_names", "abundances", "samples")]
    )
  }


  # taxonomy table
  if(!is.null(phyloseq_object@tax_table)) {
    tax_table_phylo <- phyloseq::tax_table(phyloseq_object)@.Data
    taxas <- apply(tax_table_phylo, 1, paste, collapse = ";")
    taxas <- data.frame(sequence_names = names(taxas), taxonomies = taxas)
    rownames(taxas) <- NULL
    xdev_assign_sequence_taxonomy(rdaset_object, taxas)
  }

  #sample data
  if(!is.null(phyloseq_object@sam_data)) {
    df <- data.frame(phyloseq::sample_data(phyloseq_object)@.Data)
    df <- data.frame(apply(df, 2, as.character)) # only works as a character (report.cpp line 36)
    colnames(df) <- phyloseq::sample_data(phyloseq_object)@names
    if(!is.null(treatment_column_name)) {
      df <- df[, which(colnames(df) != treatment_column_name)]
    }
    add(
        data = rdaset_object,
        table = df,
        type = "metadata"
      )
  }
  
  #phy tree
  if(!is.null(phyloseq_object@phy_tree)) {
    rdaset_object$add_sequence_tree(phyloseq::phy_tree(phyloseq_object))
  }

  #Treatments
  if(!is.null(treatment_column_name)) {
    if(any(phyloseq::sample_variables(phyloseq_object) == treatment_column_name)) {
      treatment_df <- data.frame(
        cbind(phyloseq::sample_names(phyloseq_object),
              as.character(phyloseq::get_variable(phyloseq_object)[[treatment_column_name]]))
        )
      colnames(treatment_df) <- c("samples", treatment_column_name)
      xdev_assign_treatments(data = rdaset_object, table = treatment_df,
         treatment = treatment_column_name)
    }
  }
  rdaset_object
}
