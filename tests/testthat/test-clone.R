# tests clone of sequence_data object

test_that("clone - deep copy of sequence_data object", {
  temp <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    phylo_list = rdataset_example("final.tx.list"),
    asv_list = rdataset_example("final.asv.list"),
    dataset_name = "miseq_sop"
  )

  dataset <- clone(temp)

  expect_equal(dataset$get_dataset_name(), "miseq_sop")
  expect_equal(dataset$get_num_sequences(TRUE), 2425)
  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_bins("otu"), 531)
  expect_equal(dataset$get_num_bins("phylotype"), 63)
  expect_equal(dataset$get_num_bins("asv"), 2425)
})
