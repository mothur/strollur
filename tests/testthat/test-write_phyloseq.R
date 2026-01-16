test_that("write phyloseq creates phyloseq object", {

  phyloseq_object <- import_mothur(mothur_list_file = rdataset_example("final.asv.list"),
                                  mothur_shared_file = rdataset_example("final.opti_mcc.shared"), 
                                  mothur_constaxonomy_file = rdataset_example("final.cons.taxonomy"))
  
  data <- read_phyloseq(phyloseq_object)
  write_mothur_count(data2, get_full_name("data2.taxonomy"))
  recreated_phylo_object <- write_phyloseq(data)
  expect_identical(phyloseq_object, recreated_phylo_object)
  expect_true(tax_table(phyloseq_object) == tax_table(recreated_phylo_object))
  expect_true(otu_table(phyloseq_object) == otu_table(recreated_phylo_object))
})
