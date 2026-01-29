test_that("write phyloseq creates phyloseq object", {

  phyloseq_object <- import_mothur(mothur_list_file = rdataset_example("final.asv.list"),
                                  mothur_shared_file = rdataset_example("final.opti_mcc.shared"), 
                                  mothur_constaxonomy_file = rdataset_example("final.cons.taxonomy"))

  data <- read_phyloseq(phyloseq_object)

  recreated_phylo_object <- write_phyloseq(data)
  expect_true(inherits(recreated_phylo_object, "phyloseq"))
  expect_true(all(tax_table(phyloseq_object) == tax_table(recreated_phylo_object)))
  expect_true(all(otu_table(phyloseq_object) == otu_table(recreated_phylo_object)))

  dat <- readRDS(rdataset_example("GlobalPatterns.RDS"))
  profvis::profvis({read_phyloseq(dat)})
  phyloseq_object <- read_phyloseq(dat)
  recreated_phylo_object <- write_phyloseq(phyloseq_object)
  expect_true(all(tax_table(dat) == tax_table(recreated_phylo_object), na.rm = TRUE))
  expect_true(all(otu_table(dat) == otu_table(recreated_phylo_object)))
  expect_true(length(phy_tree(dat)$edge) == length(phy_tree(recreated_phylo_object)$edge))
})


test_that("write phyloseq will fail if the dataset is empty", {
  empty_dataset <- dataset$new("")
  expect_error(write_phyloseq(empty_dataset))
})



