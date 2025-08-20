# test read_mothur_shared

test_that("test read_mothur", {
  # test complete dataset
  dataset <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  expect_equal(dataset$get_dataset_name(), "miseq_sop")
  expect_equal(dataset$get_num_sequences(TRUE), 2425)
  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_bins("otu"), 531)

  # test count table only, no samples
  dataset <- read_mothur(count = rdataset_example(
    "test_nogroups.count_table"
  ))

  expect_equal(dataset$get_num_sequences(TRUE), 6)
  expect_equal(dataset$get_num_sequences(), 295696)
  expect_equal(dataset$get_num_treatments(), 0)
  expect_equal(dataset$get_num_samples(), 0)
  expect_equal(dataset$get_num_bins(), 0)

  # test shared file and cons_tax
  dataset <- read_mothur(
    otu_shared = rdataset_example("final.opti_mcc.shared"),
    cons_taxonomy = rdataset_example(
      "final.cons.taxonomy"
    )
  )

  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_bins("otu"), 531)

  report <- dataset$get_bin_taxonomy_report()

  expect_equal(report[2707, 1], "Otu452")
  expect_equal(report[2707, 2], 1)
  expect_equal(report[2707, 3], "Bacteria")
  expect_equal(report[2707, 4], 100)
  expect_equal(report[2708, 1], "Otu452")
  expect_equal(report[2708, 2], 2)
  expect_equal(report[2708, 3], "\"Bacteroidetes\"")
  expect_equal(report[2708, 4], 100)

  # error
  expect_error(read_mothur(rabund = "bad_parameter"))
  expect_error(read_mothur(list = "non_existant_filename"))
})
