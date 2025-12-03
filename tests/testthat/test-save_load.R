# tests save of dataset object

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

  file_name <- paste0(
    normalizePath(test_path()), .Platform$file.sep,
    "test.rds"
  )

  save_dataset(temp, file_name)
  dataset <- load_dataset(file_name)
  remove_file(file_name)

  expect_equal(names(dataset, "dataset"), names(temp, "dataset"))
  expect_equal(
    num(dataset, "sequences", distinct = TRUE),
    num(temp, "sequences", distinct = TRUE)
  )
  expect_equal(
    num(dataset, "sequences"),
    num(temp, "sequences")
  )
  expect_equal(
    num(dataset, "treatments"),
    num(temp, "treatments")
  )
  expect_equal(
    num(dataset, "samples"),
    num(temp, "samples")
  )
  expect_equal(
    num(dataset, "bins", "otu"),
    num(temp, "bins", "otu")
  )
  expect_equal(
    num(dataset, "bins", "phylotype"),
    num(temp, "bins", "phylotype")
  )
  expect_equal(
    num(dataset, "bins", "asv"),
    num(temp, "bins", "asv")
  )
  expect_equal(
    totals(dataset),
    totals(temp)
  )
  expect_equal(
    totals(dataset, "treatments"),
    totals(temp, "treatments")
  )

  expect_error(load_dataset("non_existant_file.rds"))
  expect_error(save_dataset(data.frame(), "test.rds"))
})
