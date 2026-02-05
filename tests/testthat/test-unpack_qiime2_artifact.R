# test "unpack_qiime2_artifact"

test_that("unpack_qiime2_artifact - errors", {
  expect_error(unpack_qiime2_artifact("non_existant_filename"))
})

test_that("unpack_qiime2_artifact", {
  artifact <- unpack_qiime2_artifact(
    strollur_example("table.qza"),
    file_root <- get_full_name("test-qiime2")
  )

  expect_true(grepl("BIOMV", artifact$format))
  expect_equal(artifact$type, "FeatureTable[Frequency]")
  expect_equal(nrow(artifact$contents), 16)
  expect_equal(
    "6a560288-898e-4c1d-92ac-dd8d7822dcc9/data/feature-table.biom",
    artifact$contents[3, 1]
  )
  unlink(get_full_name("test-qiime2"), recursive = TRUE)
})
