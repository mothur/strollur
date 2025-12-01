# tests clone of dataset object

test_that("clone - deep copy of dataset object", {
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

  data <- clone(temp)

  expect_equal(name(data, "dataset"), "miseq_sop")
  expect_equal(num(data, distinct = TRUE), 2425)
  expect_equal(num(data), 113963)
  expect_equal(num(data, "treatments"), 2)
  expect_equal(num(data, "samples"), 19)
  expect_equal(num(data, "bins", "otu"), 531)
  expect_equal(num(data, "bins", "phylotype"), 63)
  expect_equal(num(data, "bins", "asv"), 2425)
})
