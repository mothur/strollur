# tests import of sequence_data object

test_that("import ", {
  # full dataset
  miseq <- miseq_sop_example()

  expect_equal(miseq$get_num_bins("phylotype"), 63)
  expect_equal(miseq$get_num_sequences(TRUE), 2425)

  phylo_bins_to_remove <- c("Phylo01", "Phylo02")
  reasons_to_remove <- c("testing", "testing")

  # remove some bins to allow for filtering
  remove_bins(
    miseq$data, phylo_bins_to_remove,
    reasons_to_remove, "phylotype"
  )

  expect_equal(miseq$get_num_bins("phylotype"), 61)
  expect_equal(miseq$get_num_sequences(TRUE), 825)

  exported_miseq <- miseq$export()

  expect_equal(sum(exported_miseq$sequence_data$include_sequence), 825)

  dataset <- import(exported_miseq)

  expect_equal(dataset$get_dataset_name(), miseq$get_dataset_name())
  expect_equal(dataset$get_num_sequences(TRUE), miseq$get_num_sequences(TRUE))
  expect_equal(dataset$get_num_sequences(), miseq$get_num_sequences())
  expect_equal(dataset$get_num_treatments(), miseq$get_num_treatments())
  expect_equal(dataset$get_num_samples(), miseq$get_num_samples())
  expect_equal(dataset$get_num_bins("otu"), miseq$get_num_bins("otu"))
  expect_equal(
    dataset$get_num_bins("phylotype"),
    miseq$get_num_bins("phylotype")
  )
  expect_equal(dataset$get_num_bins("asv"), miseq$get_num_bins("asv"))
  expect_equal(
    get_sample_totals(dataset$data),
    get_sample_totals(miseq$data)
  )
  expect_equal(
    get_treatment_totals(dataset$data),
    get_treatment_totals(miseq$data)
  )
})
