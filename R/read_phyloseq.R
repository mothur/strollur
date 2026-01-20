read_phyloseq <- function(phyloseq_object, dataset_name = "") {
  if(!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("To use this functionality you have to install the phyloseq package")
  }

  rdaset_object <- dataset$new(dataset_name)
  if(!is.null(phyloseq_object@otu_table)) {
    shared_table <- data.frame(create_shared_table(phyloseq_object))
    sample_df <- phyloseq::get_sample(phyloseq_object)
    xdev_add_sequences(rdaset_object, data.frame(
    sequence_names =
      unique(rownames(sample_df))
    ))
    xdev_assign_sequence_abundance(
        rdaset_object,
        shared_table[ ,c("sequence_names", "abundances", "samples")]
    )
    # assign(data = rdaset_object, table = shared_table, type = "bins")
  }

  # taxonomy table
  if(!is.null(phyloseq_object@tax_table)) {
    tax_table_phylo <- phyloseq::tax_table(phyloseq_object)@.Data
    taxas <- apply(tax_table_phylo, 1, paste, collapse = ";")
    taxas <- data.frame(sequence_names = names(taxas), taxonomies = taxas)
    rownames(taxas) <- NULL
    xdev_assign_sequence_taxonomy(rdaset_object, taxas)
  }
  
  if(!is.null(phyloseq_object@phy_tree)) {
    rdaset_object$add_sequence_tree(phyloseq::phy_tree(phyloseq_object))
  }
  rdaset_object
}


create_shared_table <- function(phyloseq_object){
  otu_data <- phyloseq::otu_table(phyloseq_object)
  sample_size <- ncol(otu_data)
  size <- sample_size * nrow(otu_data)
  shared_data <- vector("numeric", size)
  shared_sequence_names <- vector("numeric", size)
  shared_bin_names <- vector("character", size)
  shared_sample_names <- rep(colnames(otu_data), nrow(otu_data))
  rep(rownames(otu_data)[1:2])
  i <- 1
  j <- 1
  while(i < size) {
    shared_data[i:(j * sample_size)] <- otu_data[j, ]
    shared_bin_names[i:(j * sample_size)] <- paste("Otu", j, sep ="")
    shared_sequence_names <- 
    i <- i + sample_size
    j <- j + 1
  }
  return(list(
    sequence_names = unlist(lapply(rownames(otu_data), rep, times = sample_size)),
    bin_names = shared_bin_names,
    abundances = shared_data,
    samples = shared_sample_names
  ))
}
