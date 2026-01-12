read_phyloseq <- function(phyloseq_object, dataset_name = "") {
  rdaset_object <- dataset$new(dataset_name)
  if(!is.null(phyloseq_object@otu_table)) {
    count <- get_sample(phyloseq_object)
    sums <- rowSums(count)
    count_table <- data.frame(Representative_Sequence = rownames(count), total = sums,
                        count)
    rownames(count_table) <- NULL
    xdev_add_sequences(rdaset_object, data.frame(
    sequence_names =
      unique(count_table$Representative_Sequence)
    ))
  }

  # taxonomy table
  if(!is.null(phyloseq_object@tax_table)) {
    tax_table_phylo <- tax_table(phyloseq_object)@.Data
    taxas <- apply(tax_table_phylo, 1, paste, collapse = ";")
    taxas <- data.frame(sequence_names = names(taxas), taxonomies = taxas)
    rownames(taxas) <- NULL
    xdev_assign_sequence_taxonomy(rdaset_object, taxas)
  }
  
  if(!is.null(phyloseq_object@phy_tree)) {
    rdaset_object$add_sequence_tree(phy_tree(phyloseq_object))
  }
  rdaset_object
}
