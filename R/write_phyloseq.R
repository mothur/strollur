write_phyloseq <- function(dataset) {
  dataset <- eso

  ab <- abundance(data = data1, type = "sequences")
  if(nrow(abundance(data = dataset, type = "sequences")) > 0) {
    abundances <- abundance(data = dataset, type = "sequences", by_sample = TRUE) 
    sequence_names <- names(data = dataset, type = "samples")
    abundances <- reshape(abundances,
            direction = "wide",
            timevar = "samples",
            idvar = "sequence_names",
    )
    colnames(abundances) <- sub(".*\\.", "", colnames(abundances))
    otu_table <- as.matrix(abundances[,2:ncol(abundances)])
    rownames(otu_table) <- abundances$sequence_names
  }
  

  # taxonomies
  
  taxas <- 
    reshape(report(data = dataset, type = "sequence_taxonomy"),
            direction = "wide",
            timevar = "level",
            idvar = "id",
    )

  colnames(taxas) <- c("id", paste("level_", seq(1, ncol(taxas) - 1), sep=""))
  rownames(taxas) <- taxas$id  
  taxas <- as.matrix(taxas[, colnames(taxas) != "id"])
  return(phyloseq(otu_table(otu_table, T), tax_table(taxas), phy_tree(dataset$get_sequence_tree())))

  }

# dataset <- read_phyloseq(GlobalPatterns)
# phylo_test <- write_phyloseq(dataset)
# # data(GlobalPatterns)
# library(phyloseq)
# a <- read_phyloseq(GlobalPatterns)  

# data(esophagus)

# obj <- c()

# s <- sample_data(GlobalPatterns)
