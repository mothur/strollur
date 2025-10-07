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

  save(temp, file_name)
  dataset <- load(file_name)
  remove_file(file_name)

  expect_equal(dataset$get_dataset_name(), temp$get_dataset_name())
  expect_equal(dataset$get_num_sequences(TRUE), temp$get_num_sequences(TRUE))
  expect_equal(dataset$get_num_sequences(), temp$get_num_sequences())
  expect_equal(dataset$get_num_treatments(), temp$get_num_treatments())
  expect_equal(dataset$get_num_samples(), temp$get_num_samples())
  expect_equal(dataset$get_num_bins("otu"), temp$get_num_bins("otu"))
  expect_equal(
    dataset$get_num_bins("phylotype"),
    temp$get_num_bins("phylotype")
  )
  expect_equal(dataset$get_num_bins("asv"), temp$get_num_bins("asv"))
  expect_equal(
    get_sample_totals(dataset$data),
    get_sample_totals(temp$data)
  )
  expect_equal(
    get_treatment_totals(dataset$data),
    get_treatment_totals(temp$data)
  )

  expect_error(load("non_existant_file.rds"))
  expect_error(save(data.frame(), "test.rds"))
})
