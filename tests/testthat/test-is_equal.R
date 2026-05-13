# tests is_equal function

test_that("test equal strollur objects", {
  miseq <- miseq_sop_example()
  data <- copy_dataset(miseq)

  expect_true(miseq$is_equal(data))
  expect_error(miseq$is_equal("not a strollur object"))
  expect_error(is_equal(data, "not a strollur object"))
  expect_error(is_equal("not a strollur object", data))
})

test_that("test strollur objects with different sequence_tree / sample_tree", {
  sample_tree <- strollur_example("final.opti_mcc.jclass.ave.tre")
  sequence_tree <- strollur_example("final.phylip.tre.gz")

  data <- new_dataset()
  seq_tree <- ape::read.tree(sequence_tree)
  data$add_sequence_tree(seq_tree)

  data2 <- new_dataset()
  samp_tree <- ape::read.tree(sample_tree)
  data2$add_sequence_tree(samp_tree)

  # mismatched sequence trees
  actual_message <- evaluate_promise(is_equal(data, data2))$messages
  expect_true(grepl("'sequence_tree' are not equivalent", actual_message))
  expect_false(is_equal(data, data2))

  # null tree to filled tree
  data2$sequence_tree <- NULL
  actual_message <- evaluate_promise(is_equal(data, data2))$messages
  expect_true(grepl("'sequence_tree' are not equivalent", actual_message))
  expect_false(is_equal(data, data2))

  data$clear()
  data2$clear()

  expect_true(is_equal(data, data2))

  df <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))
  xdev_assign_bins(data = data, table = df, bin_type = "otu")
  xdev_assign_bins(data = data2, table = df, bin_type = "otu")

  data$add_sample_tree(samp_tree)
  data2$add_sample_tree(samp_tree)

  # remove sample from data to force tree diff
  xdev_remove_samples(data, c("F3D0"))

  # mismatched sample trees
  actual_message <- evaluate_promise(is_equal(data, data2))$messages
  expect_true(grepl("'sample_tree' are not equivalent", actual_message))
  expect_false(is_equal(data, data2))

  # null tree to filled tree
  data2$sample_tree <- NULL
  actual_message <- evaluate_promise(is_equal(data, data2))$messages
  expect_true(grepl("'sample_tree' are not equivalent", actual_message))
  expect_false(is_equal(data, data2))
})

test_that("test strollur objects with different processors / versions ", {
  data <- strollur$new()
  data2 <- strollur$new()

  data2[[".__enclos_env__"]]$private$version <- "not a matching version"
  data2[[".__enclos_env__"]]$private$processors <- -1

  data_private <- data[[".__enclos_env__"]]$private
  data2_private <- data2[[".__enclos_env__"]]$private

  expect_false(data_private$version == data2_private$version)
  expect_false(data_private$processors == data2_private$processors)
})

test_that("test strollur objects with different raw items ", {
  data <- strollur$new()
  data2 <- strollur$new()

  # set data's raw member
  xint_serialize_dobject(data)

  actual_message <- evaluate_promise(is_equal(data, data2))$messages
  expect_true(grepl("field 'raw' are not  equivalent", actual_message))
  expect_false(is_equal(data, data2))
})

test_that("test strollur objects with different c++ back end ", {
  data <- strollur$new()
  data2 <- strollur$new()

  # add sequences to data
  fasta_data <- read_fasta(strollur_example("final.fasta.gz"))
  xdev_add_sequences(data = data, table = fasta_data)

  actual_message <- evaluate_promise(is_equal(data, data2))$messages
  expect_true(grepl("field 'data' are not  equivalent", actual_message))
  expect_false(is_equal(data, data2))
})
