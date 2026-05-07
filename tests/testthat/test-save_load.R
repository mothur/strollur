# tests save of dataset object

test_that("clone - deep copy of dataset object", {
  temp <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    taxonomy = strollur_example("final.taxonomy.gz"),
    design = strollur_example("mouse.time.design"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    phylo_list = strollur_example("final.tx.list.gz"),
    asv_list = strollur_example("final.asv.list.gz"),
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
    count(dataset, "sequence", distinct = TRUE),
    count(temp, "sequence", distinct = TRUE)
  )
  expect_equal(
    count(dataset, "sequence"),
    count(temp, "sequence")
  )
  expect_equal(
    count(dataset, "treatment"),
    count(temp, "treatment")
  )
  expect_equal(
    count(dataset, "sample"),
    count(temp, "sample")
  )
  expect_equal(
    count(dataset, "bin", "otu"),
    count(temp, "bin", "otu")
  )
  expect_equal(
    count(dataset, "bin", "phylotype"),
    count(temp, "bin", "phylotype")
  )
  expect_equal(
    count(dataset, "bin", "asv"),
    count(temp, "bin", "asv")
  )
  expect_equal(
    abundance(dataset, "sample"),
    abundance(temp, "sample")
  )
  expect_equal(
    abundance(dataset, "treatment"),
    abundance(temp, "treatment")
  )

  expect_error(load_dataset("non_existant_file.rds"))
  expect_error(save_dataset(data.frame(), "test.rds"))
})
