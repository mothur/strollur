# tests save of dataset object

test_that("load - version 0.1.1", {
    # current version
    miseq <- miseq_sop_example()

    file_name <- file.path(tempdir(), "test.rds")
    save_dataset(miseq, file_name)

    data <- load_dataset(file_name)
    remove_file(file_name)

    expect_equal(count(data), 113963)

    # 0.1.1 not compatible, dataset c++ class members changed.
    # metadata is no longer a member of dataset, causes a serialization issue
    # import is possible with table form
    expect_error(load_dataset(strollur_example("miseq_sop.010.rds")))
})

test_that("load - version 0.1.0", {
  # set version to 0.0.0 - forced compatible
  miseq <- miseq_sop_example()

  miseq[[".__enclos_env__"]]$private$version <- "0.0.0"

  file_name <- file.path(tempdir(), "test.rds")
  save_dataset(miseq, file_name)

  message <- capture_messages(load_dataset(file_name))
  remove_file(file_name)

  expect_true(grepl("Converting and importing", message))

  # set for 1.0.0 - not compatible
  miseq[[".__enclos_env__"]]$private$version <- "1.0.0"

  save_dataset(miseq, file_name)

  expect_error(load_dataset(file_name))
})

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

  file_name <- file.path(tempdir(), "test.rds")
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

  file_name <- file.path(tempdir(), "test.rds")
  expect_error(load_dataset("non_existant_file.rds"))
  expect_error(save_dataset(data.frame(), file_name))
})
