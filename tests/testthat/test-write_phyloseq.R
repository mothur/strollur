test_that("write phyloseq creates phyloseq object", {

  phyloseq_object <- import_mothur(mothur_list_file = rdataset_example("final.asv.list"),
                                  mothur_shared_file = rdataset_example("final.opti_mcc.shared"), 
                                  mothur_constaxonomy_file = rdataset_example("final.cons.taxonomy"))
  
  data <- read_phyloseq(phyloseq_object)
  recreated_phylo_object <- write_phyloseq(data)
  expect_true(inherits(recreated_phylo_object, "phyloseq"))
  # expect_identical(phyloseq_object, recreated_phylo_object)
  expect_true(all(tax_table(phyloseq_object) == tax_table(recreated_phylo_object)))
  expect_true(all(otu_table(phyloseq_object) == otu_table(recreated_phylo_object)))

  data(GlobalPatterns)
  phyloseq_object <- read_phyloseq(GlobalPatterns)
  recreated_phylo_object <- write_phyloseq(phyloseq_object)
  expect_true(all(tax_table(GlobalPatterns) == tax_table(recreated_phylo_object)))
  expect_true(all(otu_table(GlobalPatterns) == otu_table(recreated_phylo_object)))
  expect_true(length(phy_tree(GlobalPatterns)$edge) == length(phy_tree(recreated_phylo_object)$edge))
})


test_that("write phyloseq will fail if the dataset is empty", {

  empty_dataset <- dataset$new("")
  expect_error(write_phyloseq(empty_dataset))
})



