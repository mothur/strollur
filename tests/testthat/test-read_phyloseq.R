test_that("test read_phyloseq", {
  phyloseq_object <- phyloseq::import_mothur(rdataset_example("esophagus.fn.list"), rdataset_example("esophagus.good.groups"), 
                          rdataset_example("esophagus.tree"))
  data <- read_phyloseq(phyloseq_object)
  expect_true(count(data) == 675)

  phyloseq_object_w_tax <- phyloseq::import_mothur(mothur_list_file = rdataset_example("final.asv.list"),
                                         mothur_shared_file = rdataset_example("final.opti_mcc.shared"), 
                                         mothur_constaxonomy_file = rdataset_example("final.cons.taxonomy"))
  data <- read_phyloseq(phyloseq_object_w_tax)
  row_data <- phyloseq::tax_table(phyloseq_object_w_tax)[9, ]
  rdata_set_data <- report(data, type = "sequence_taxonomy")
  expect_true(abundance(data)[1, ]$abundances == 12288)
  expect_true(abundance(data)[2, ]$abundances == 8892)
  expect_true(abundance(data)[3, ]$abundances == 7794)

  expect_true(nrow(rdata_set_data) == 3186)
  expect_true(rdata_set_data[49, 3] == row_data[[1]])
  expect_true(rdata_set_data[50, 3] == row_data[[2]])
  expect_true(rdata_set_data[51, 3] == row_data[[3]])
  expect_true(rdata_set_data[52, 3] == row_data[[4]])
  expect_true(rdata_set_data[53, 3] == row_data[[5]])
  expect_true(rdata_set_data[54, 3] == row_data[[6]])

})

test_that("read phyloseq fails if not given an phyloseq object", {
  empty <- c()
  expect_error(read_phyloseq(empty), 
               "phyloseq_object has to an object created using the phyloseq package.")
})

