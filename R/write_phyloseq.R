write_phyloseq <- function(data) {
  if(!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("To use this functionality you have to install the phyloseq package")
  }

  phyloseq_parameter_list <- vector("list", 4)
  if(nrow(abundance(data = data, type = "sequences")) > 0) {
    abundances <- abundance(data = data, type = "sequences", by_sample = TRUE) 
    treatments <- NULL
    if(any(colnames(abundances) == "treatments")) {
      treatments <- abundances$treatments
      abundances <- abundances[, -which(colnames(abundances) == "treatments")]
    }
    sequence_names <- names(data = data, type = "samples")
    abundances <- reshape(abundances,
            direction = "wide",
            timevar = "samples",
            idvar = "sequence_names",
    )
    abundances[is.na(abundances)] <- 0
    colnames(abundances) <- sub(".*\\.", "", colnames(abundances))
    otu_table <- as.matrix(abundances[,2:ncol(abundances)])
    rownames(otu_table) <- abundances$sequence_names
    phyloseq_parameter_list[[1]] <- phyloseq::otu_table(otu_table, TRUE)
  }
  

  # taxonomies
  if(nrow(report(data = data, type = "sequence_taxonomy")) > 0) { 
    df <- report(data = data, type = "sequence_taxonomy")
    
    if(any(colnames(df) == "confidence")) {
      df$taxon <- paste0(df$taxon, "(", df$confidence, ")")
      df$confidence <- NULL
    }
    taxas <- 
      reshape(df,
              direction = "wide",
              timevar = "level",
              idvar = "id",
      )

    colnames(taxas) <- c("id", paste("level_", seq(1, ncol(taxas) - 1), sep=""))
    rownames(taxas) <- taxas$id  
    taxas <- as.matrix(taxas[, colnames(taxas) != "id"])
    phyloseq_parameter_list[[2]] <- phyloseq::tax_table(taxas)
  }
  
  if(!is.null(data$get_sequence_tree())) {
     phyloseq_parameter_list[[3]] <- phyloseq::phy_tree(data$get_sequence_tree())
  }

  if(!is.null(report(data, "metadata"))) {
    df <- report(data, "metadata")
    rownames(df) <- df$rownames
    df$rownames <- NULL
    phyloseq_parameter_list[[4]] <- phyloseq::sample_data(df)
  } 

  indexes <- which(!sapply(phyloseq_parameter_list, is.null))
  if(length(indexes) <= 0) {
    stop("You have an empty object that cannot become a phyloseq object.")
  }

  return(do.call(phyloseq::phyloseq,
     phyloseq_parameter_list[indexes]))  
}
# library(phyloseq)
# data(GlobalPatterns)

# dat <- read_phyloseq(GlobalPatterns)
# phylo <- write_phyloseq(dat)
# # x <- function() {
# # #   browser()
# # #   tax_table(df)
# # # }

# # # sample_data(df)
# # # x()
# # names(data = dat, type = "samples")
# sam <- report(dat, "metadata")
# # phylo <- 
#  df <- report(dat, "metadata")
#     rownames(df) <- df$rownames
#     df$rownames <- NULL
# d <- phyloseq::sample_data(df)
