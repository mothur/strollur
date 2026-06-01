# tests is_equal function

test_that("test equal strollur objects", {
  miseq <- miseq_sop_example()
  data <- copy_dataset(miseq)

  expect_true(miseq$is_equal(data))

  table <- export_dataset(miseq)

  data <- import_dataset(table)

  expect_true(miseq$is_equal(data))

  clear(miseq)
  clear(data)

  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")
  comments <- c("ddd", "ftf", "efr", "ssd")

  ids <- c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2",
    "seq3", "seq3",
    "seq4", "seq4"
  )
  samples <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3",
    "sample2", "sample4"
  )
  abundances <- c(
    250, 400, 500,
    25, 40, 50,
    25, 25,
    1, 4
  )
  treatments <- c(
    "early", "early", "late",
    "early", "early", "late",
    "early", "early",
    "early", "late"
  )

  data <- new_dataset("test")

  # include reference
  url <- "https://mothur.org/wiki/silva_reference_files/"
  xdev_add_sequences(
    data,
    data.frame(sequence_name = names, sequence = seqs),
    new_reference(
      name = "silva.bacteria.fasta", version = "1.38.1",
      usage = "alignment by mothur2 v1.0", documentation_url = url
    )
  )

  data$assign(
    table = data.frame(
      sequence_name = ids, abundance = abundances,
      sample = samples, treatment = treatments
    ), type = "sequence_abundance"
  )

  # assign bins
  bins <- c("bin1", "bin2", "bin1", "bin2")
  assign(
    data = data,
    table = data.frame(bin_name = bins, sequence_name = names),
    type = "bin"
  )

  table <- export_dataset(data)

  data2 <- import_dataset(table)
  table2 <- export_dataset(data2)

  expect_equal(table$sequence_data, table2$sequence_data)
  expect_equal(table$sequence_report, table2$sequence_report)
  expect_equal(table$sequence_abundance_table, table2$sequence_abundance_table)
  expect_equal(table$otu_bin_data, table2$otu_bin_data)
  expect_equal(
    table$otu_sequence_bin_assignment,
    table2$otu_sequence_bin_assignment
  )
  expect_equal(table$resource_reference, table2$resource_reference)
  expect_true(data2$is_equal(data))

  expect_error(data2$is_equal("not a strollur object"))
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

test_that("test strollur objects with different versions ", {
  data <- strollur$new()
  data2 <- strollur$new()

  data2[[".__enclos_env__"]]$private$version <- "not a matching version"

  data_private <- data[[".__enclos_env__"]]$private
  data2_private <- data2[[".__enclos_env__"]]$private

  expect_false(data_private$version == data2_private$version)
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
  data2 <- strollur$new("different name")

  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("have different dataset names", message))
  expect_false(is_equal(data, data2))

  xdev_set_dataset_name(data2, "")
  expect_true(is_equal(data, data2))

  # add sequences to data
  fasta_data <- read_fasta(strollur_example("final.fasta.gz"))
  xdev_add_sequences(data = data, table = fasta_data)

  # 'The strollur objects have different sequence data statuses.'
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different sequence data statuses", message))
  expect_false(is_equal(data, data2))

  xdev_add_sequences(data = data2, table = fasta_data)
  expect_true(is_equal(data, data2))

  # add taxonomy for sequences
  taxonomy <- readRDS(strollur_example("miseq_tidy_taxonomy.rds"))
  xdev_assign_sequence_taxonomy_tidy(data, taxonomy)

  # The strollur objects have different classification statuses.
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different classification statuses", message))
  expect_false(is_equal(data, data2))

  xdev_assign_sequence_taxonomy_tidy(data2, taxonomy)
  expect_true(is_equal(data, data2))

  # remove sequence from data
  name_to_remove <- "M00967_43_000000000-A3JHG_1_2101_16474_12783"
  xdev_remove_sequences(data,
    sequence_names = c(name_to_remove),
    c("test")
  )

  # The strollur objects have different number of unique sequences
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different number of unique sequences", message))
  expect_false(is_equal(data, data2))

  xdev_remove_sequences(data2,
    sequence_names = c(name_to_remove),
    c("test")
  )

  expect_true(is_equal(data, data2))

  clear(data)
  clear(data2)

  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")
  comments <- c("ddd", "ftf", "efr", "ssd")
  seq_taxes <- c("A;B;", "A;D;", "A;B;", "A;B;C;")

  xdev_add_sequences(
    data,
    data.frame(sequence_name = names)
  )

  xdev_add_sequences(
    data2,
    data.frame(sequence_name = c("seq10", "seq2", "seq3", "seq4"))
  )

  # The strollur objects include different sequence names
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different sequence names", message))
  expect_false(is_equal(data, data2))

  clear(data2)
  clear(data)

  xdev_add_sequences(
    data,
    data.frame(sequence_name = names, sequence = seqs)
  )

  xdev_add_sequences(
    data2,
    data.frame(
      sequence_name = names,
      sequence = c("AGC", "ATTG", "ATGC", "ATTC")
    )
  )

  # The strollur objects include different sequence strings
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different sequence strings", message))
  expect_false(is_equal(data, data2))

  clear(data2)
  clear(data)

  xdev_add_sequences(
    data,
    data.frame(
      sequence_name = names,
      sequence = seqs,
      comment = comments
    )
  )

  xdev_add_sequences(
    data2,
    data.frame(
      sequence_name = names,
      sequence = seqs
    )
  )

  # The strollur objects include different sequence comments
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different sequence comments", message))
  expect_false(is_equal(data, data2))

  clear(data2)
  clear(data)

  xdev_assign_sequence_taxonomy(
    data,
    data.frame(
      sequence_name = names,
      taxonomy = seq_taxes
    )
  )

  xdev_assign_sequence_taxonomy(
    data2,
    data.frame(
      sequence_name = names,
      taxonomy = c(
        "A;D;",
        "A;D;",
        "A;B;", "A;B;C;"
      )
    )
  )

  # The strollur objects include different sequence taxonomies
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different sequence taxonomies", message))
  expect_false(is_equal(data, data2))

  clear(data2)
  clear(data)

  bins <- c("bin1", "bin2", "bin1", "bin2")

  assign(
    data = data,
    table = data.frame(bin_name = bins, sequence_name = names),
    type = "bin"
  )

  xdev_add_sequences(
    data2,
    data.frame(sequence_name = names)
  )

  # The strollur objects include different sequence bin assignment statuses
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("sequence bin assignment statuses", message))
  expect_false(is_equal(data, data2))

  assign(
    data = data2,
    table = data.frame(bin_name = bins, sequence_name = names),
    type = "bin"
  )

  url <- "https://mothur.org/wiki/silva_reference_files/"
  xdev_add_references(
    data,
    data.frame(
      name = "silva.bacteria.fasta", version = "1.38.1",
      usage = "alignment by mothur2 v1.0", documentation_url = url
    )
  )

  # The strollur objects include different resource_references
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different resource_references", message))
  expect_false(is_equal(data, data2))

  ids <- c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2",
    "seq3", "seq3",
    "seq4", "seq4"
  )
  samples <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3",
    "sample2", "sample4"
  )
  abundances <- c(
    250, 400, 500,
    25, 40, 50,
    25, 25,
    1, 4
  )
  treatments <- c(
    "early", "early", "late",
    "early", "early", "late",
    "early", "early",
    "early", "late"
  )

  clear(data2)
  clear(data)

  data$assign(
    table = data.frame(
      sequence_name = ids, abundance = abundances,
      sample = samples, treatment = treatments
    ), type = "sequence_abundance"
  )

  # add same number of sequences to different samples
  data2$assign(
    table = data.frame(
      sequence_name = ids,
      abundance = c(
        4, 400, 500,
        25, 40, 50,
        25, 25,
        1, 250
      ),
      sample = samples, treatment = treatments
    ), type = "sequence_abundance"
  )

  # The strollur objects include different sequence abundance data
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different sequence abundance data", message))
  expect_false(is_equal(data, data2))

  data2$assign(
    table = data.frame(
      sequence_name = ids, abundance = abundances,
      sample = samples, treatment = treatments
    ), type = "sequence_abundance"
  )
  expect_true(is_equal(data, data2))

  assign(
    data = data,
    table = data.frame(bin_name = bins, sequence_name = names),
    type = "bin"
  )

  assign(
    data = data2,
    table = data.frame(
      bin_name = c("bin1", "bin1", "bin1", "bin2"),
      sequence_name = names
    ),
    type = "bin"
  )

  # The strollur objects include different bin data
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different bin data", message))
  expect_false(is_equal(data, data2))

  clear(data2)
  clear(data)

  contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
  xdev_add_report(
    data = data, table = contigs_report, type = "contigs_report",
    sequence_name = "Name"
  )

  xdev_add_sequences(
    data = data2, table = contigs_report["Name"],
    sequence_name = "Name"
  )

  # The strollur objects include different report data.
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different report data", message))
  expect_false(is_equal(data, data2))

  clear(data2)
  clear(data)

  metadata <- readRDS(strollur_example("miseq_metadata.rds"))
  xdev_add_report(data = data, table = metadata, type = "metadata")

  # The strollur objects include different metadata.
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different metadata", message))
  expect_false(is_equal(data, data2))

  clear(data2)
  clear(data)

  xdev_add_sequences(
    data,
    data.frame(sequence_name = names, sequence = seqs)
  )

  xdev_add_sequences(
    data2,
    data.frame(sequence_name = names, sequence = seqs)
  )

  xdev_remove_sequences(data2,
    sequence_names = c("seq1"),
    c("test")
  )

  xdev_remove_sequences(data,
    sequence_names = c("seq1"),
    c("test-different-trash-code")
  )

  # The strollur objects include different sequence trash codes
  message <- capture_output(is_equal(data, data2))

  expect_true(grepl("different sequence trash codes", message))
  expect_false(is_equal(data, data2))
})
