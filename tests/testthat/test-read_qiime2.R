# test read_qiime2_taxonomy

test_that("test read_qiime2 - errors", {
  expect_error(read_qiime2("non_existant_filename"))
  expect_error(read_qiime2(
    strollur_example("table.qza"),
    file_root <- get_full_name("test-qiime2")
  ))
  unlink(get_full_name("test-qiime2"), recursive = TRUE)
})

test_that("test read_qiime2", {
  qza_files <- c(
    strollur_example("rep_seqs.qza"),
    strollur_example("table.qza"),
    strollur_example("taxonomy.qza"),
    strollur_example("rooted-tree.qza")
  )

  data <- read_qiime2(
    qza = qza_files,
    metadata = strollur_example("sample_metadata.tsv"),
    dataset_name = "qiime_moving_pictures"
  )

  expect_equal(names(data, type = "dataset"), "qiime_moving_pictures")
  expect_equal(count(data), 157298)
  expect_equal(count(data, distinct = TRUE), 759)
  expect_equal(count(data, type = "bins", bin_type = "asv"), 759)
  expect_equal(count(data, type = "bins", bin_type = "otu"), 0)

  sample_totals <- abundance(data, type = "samples")

  sample_abunds <- c(
    7865, 7245, 8270, 6486, 6755, 8756, 7922, 7068, 4112,
    4545, 3340, 3485, 5146, 1549, 2526, 4166, 917, 1313,
    1191, 1109, 1130, 1279, 8575, 9961, 10095, 2253, 1827,
    1969, 2132, 2555, 1817, 6892, 6022, 7025
  )

  sample_names <- c(
    "L1S105", "L1S140", "L1S208", "L1S257", "L1S281",
    "L1S57", "L1S76", "L1S8", "L2S155", "L2S175", "L2S204",
    "L2S222", "L2S240", "L2S309", "L2S357", "L2S382",
    "L3S242", "L3S294", "L3S313", "L3S341", "L3S360",
    "L3S378", "L4S112", "L4S137", "L4S63", "L5S104",
    "L5S155", "L5S174", "L5S203", "L5S222", "L5S240",
    "L6S20", "L6S68", "L6S93"
  )

  expect_equal(sample_totals$abundances, sample_abunds)
  expect_equal(sample_totals$samples, sample_names)
})
