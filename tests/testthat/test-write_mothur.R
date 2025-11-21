test_that("write_mothur - errors", {
  expect_error(write_mothur("Bad_type"))
})

test_that("write_mothur - miseq", {
  miseq <- miseq_sop_example()

  outputs <- write_mothur(miseq,
    dir_path = get_full_name("tmp"),
    compress = TRUE
  )
  zip_file <- get_full_name("tmp.zip")

  unzip(zip_file, exdir = get_full_name("tmp"), junkpaths = TRUE)

  data <- read_mothur(
    fasta = outputs[1],
    count = outputs[3],
    taxonomy = outputs[2],
    design = outputs[4],
    otu_list = outputs[5],
    asv_list = outputs[6],
    phylo_list = outputs[7],
    sample_tree = outputs[17],
    sequence_tree = outputs[18],
    dataset_name = "miseq_sop"
  )

  add_references(data, readr::read_tsv(outputs[16],
    col_names = TRUE,
    show_col_types = FALSE
  ))

  add_report(
    data,
    readr::read_tsv(outputs[14],
      col_names = TRUE,
      show_col_types = FALSE
    ),
    "metadata"
  )

  for (output in outputs) {
    remove_file(output)
  }
  remove_file(zip_file)
  unlink(get_full_name("tmp/"), recursive = TRUE)

  expect_equal(data$get_dataset_name(), "miseq_sop")
  expect_equal(data$get_num_sequences(TRUE), 2425)
  expect_equal(data$get_num_sequences(), 113963)
  expect_equal(data$get_num_treatments(), 2)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_bins("otu"), 531)
  expect_equal(data$get_num_bins("asv"), 2425)
  expect_equal(data$get_num_bins("phylotype"), 63)
  expect_equal(data$get_metadata(), miseq$get_metadata())
  expect_equal(report(data, "references"), report(miseq, "references"))
})

test_that("write_mothur - bin_data only", {
  miseq <- miseq_sop_example()

  outputs <- write_mothur(miseq,
    dir_path = get_full_name("tmp"),
    compress = TRUE,
    tags = c("bin_data")
  )
  zip_file <- get_full_name("tmp.zip")

  unzip(zip_file, exdir = get_full_name("tmp"), junkpaths = TRUE)

  data <- read_mothur(
    count = outputs[4],
    design = outputs[11],
    otu_list = outputs[1],
    asv_list = outputs[2],
    phylo_list = outputs[3],
    dataset_name = "miseq_sop"
  )

  for (output in outputs) {
    remove_file(output)
  }
  remove_file(zip_file)
  unlink(get_full_name("tmp/"), recursive = TRUE)

  expect_equal(data$get_dataset_name(), "miseq_sop")
  expect_equal(data$get_num_sequences(TRUE), 2425)
  expect_equal(data$get_num_sequences(), 113963)
  expect_equal(data$get_num_treatments(), 2)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_bins("otu"), 531)
  expect_equal(data$get_num_bins("asv"), 2425)
  expect_equal(data$get_num_bins("phylotype"), 63)
})

test_that("write_mothur - report_data only", {
  data <- dataset$new()

  add_report(data, readr::read_tsv(
    rdataset_example("alignment_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  ), "alignment_report", "QueryName")

  add_report(data, readr::read_tsv(
    rdataset_example("contigs_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  ), "contigs_report", "Name")

  outputs <- write_mothur(data,
    dir_path = get_full_name("tmp"),
    compress = TRUE,
    tags = c("reports")
  )
  zip_file <- get_full_name("tmp.zip")

  unzip(zip_file, exdir = get_full_name("tmp"), junkpaths = TRUE)

  data2 <- dataset$new()

  seq_col <- list("QueryName", "Name")
  names(seq_col) <- outputs

  for (output in outputs) {
    extension <- tail(strsplit(output, "\\.")[[1]], 1)

    add_report(data2, readr::read_tsv(
      output,
      col_names = TRUE, show_col_types = FALSE
    ), extension, seq_col[[output]])
  }

  for (output in outputs) {
    remove_file(output)
  }
  remove_file(zip_file)
  unlink(get_full_name("tmp/"), recursive = TRUE)

  expect_equal(data$get_dataset_name(), "")
  expect_equal(data$get_num_sequences(TRUE), 5)
  expect_equal(data$get_num_sequences(), 5)
  expect_equal(
    report(data, "contigs_report"),
    report(data2, "contigs_report")
  )
  expect_equal(
    report(data, "alignment_report"),
    report(data2, "alignment_report")
  )

  data <- dataset$new()

  add_report(data, readr::read_tsv(
    rdataset_example("chimera_report.tsv"),
    col_names = TRUE, show_col_types = FALSE
  ), "chimera_report", "Query")

  outputs <- write_mothur(data,
    dir_path = get_full_name("tmp"),
    compress = TRUE
  )

  zip_file <- get_full_name("tmp.zip")

  unzip(zip_file, exdir = get_full_name("tmp"), junkpaths = TRUE)

  data2 <- dataset$new()

  add_report(data2, readr::read_tsv(
    outputs[2],
    col_names = TRUE, show_col_types = FALSE
  ), "chimera_report", "Query")

  for (output in outputs) {
    remove_file(output)
  }
  remove_file(zip_file)
  unlink(get_full_name("tmp/"), recursive = TRUE)

  expect_equal(data$get_dataset_name(), "")
  expect_equal(data$get_num_sequences(TRUE), 71)
  expect_equal(data$get_num_sequences(), 71)
  expect_equal(
    report(data, "chimera_report"),
    report(data2, "chimera_report")
  )
})
