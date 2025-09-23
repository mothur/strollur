# test "rcpp_dataset"

test_that("rcpp_dataset - new_dataset, copy_dataset", {
  dataset <- new_dataset("miseq_sop", 4)
  expect_equal(get_dataset_name(dataset), "miseq_sop")
  expect_equal(get_num_sequences(dataset), 0)
  expect_equal(get_num_processors(dataset), 4)

  dataset2 <- copy_dataset(dataset)
  expect_equal(get_dataset_name(dataset2), "miseq_sop")
  expect_equal(get_num_sequences(dataset2), 0)
  expect_equal(get_num_processors(dataset2), 4)

  add_sequences(dataset, c("seq1"), "", "")
  expect_equal(get_num_sequences(dataset), 1)
  expect_equal(get_num_sequences(dataset2), 0)
})

test_that("rcpp_dataset - add_sequences", {
  dataset <- new_dataset("miseq_sop", 4)

  sequences <- read_fasta(rdataset_example("final.fasta"))
  comments <- rep("my comments", length(sequences$names))

  # add with all 3 parameters
  add_sequences(
    dataset, sequences$sequence_names,
    sequences$sequences, comments
  )
  expect_equal(get_num_sequences(dataset), 2425)
  expect_true(has_sequence_strings(dataset))
  expect_error(set_sequences(
    dataset, sequences$names,
    sequences$sequences,
    c("not_enough_comments", "comment")
  ))

  clear(dataset, "")

  # add with 2 parameters
  add_sequences(dataset, sequences$sequence_names, sequences$sequences, "")
  expect_equal(get_num_sequences(dataset), 2425)
  expect_true(has_sequence_strings(dataset))

  seq_report <- get_sequence_report(dataset)

  first_three_seqs <- c(
    "M00967_43_000000000-A3JHG_1_2101_16474_12783",
    "M00967_43_000000000-A3JHG_1_1113_12711_3318",
    "M00967_43_000000000-A3JHG_1_2108_14707_9807"
  )
  expect_equal((seq_report[[1]][1:3]), first_three_seqs)
  expect_equal((seq_report[[4]][1:3]), c(253, 253, 253))
  expect_equal((seq_report[[5]][10:12]), c(0, 0, 0))
  expect_equal((seq_report[[7]][10:12]), c(0, 0, 0))
  expect_equal((seq_report[[2]][10:12]), c(1, 1, 1))

  clear(dataset, "")

  # add with only names
  add_sequences(dataset, sequences$sequence_names, "", "")
  expect_equal(get_num_sequences(dataset), 2425)
  expect_false(has_sequence_strings(dataset))

  clear(dataset, "")

  expect_error(add_sequences(dataset))

  message <- capture_output(set_sequences(
    dataset,
    c("non_existant_sequence"),
    c("ATGC"), ""
  ))
  expect_true(grepl("non_existant_sequence", message))
})

test_that("rcpp_dataset - assign_bins", {
  dataset <- new_dataset("miseq_sop", 4)

  # bin labels and seq names (list)
  seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  assign_bins(dataset, bin_ids, 0, "", seq_ids)

  expect_equal(get_num_sequences(dataset), 6)
  expect_equal(get_num_bins(dataset), 3)
  expect_false(has_sequence_strings(dataset))

  abundances <- c(110, 525, 80)
  samples <- c("sample1", "sample1", "sample2")
  expect_error(assign_bins(dataset, unique(bin_ids), abundances, "", ""))
  expect_error(assign_bins(dataset, unique(bin_ids), 0, samples, ""))
  expect_error(assign_bins(dataset, unique(bin_ids), abundances, "", ""))

  clear(dataset, "")

  # bin labels and bin abundances (rabund)
  assign_bins(dataset, unique(bin_ids), abundances, "", "")

  expect_equal(get_num_sequences(dataset), 715)
  expect_equal(get_num_bins(dataset), 3)

  clear(dataset, "")

  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  samples <- c(
    "sample1", "sample2", "sample5", "sample1", "sample3",
    "sample1"
  )
  sample_abundances <- c(10, 100, 1, 500, 25, 80)

  expect_error(assign_bins(dataset, "", sample_abundances, samples, ""))

  # bin labels, samples and bin abundances (shared)
  assign_bins(dataset, bin_ids, sample_abundances, samples, "")

  expect_equal(get_num_sequences(dataset), 716)
  expect_equal(get_num_bins(dataset), 3)
  expect_equal(get_num_samples(dataset), 4)
  expect_equal(get_sample_totals(dataset), c(590, 100, 25, 1))
  expect_equal(get_bin(dataset, "bin1"), "")
  expect_equal(get_bin_abundances(dataset, "bin1"), c(10, 100, 0, 1))

  clear(dataset, "")

  bin_ids <- c(
    "bin1", "bin1", "bin1", "bin1", "bin1", "bin1",
    "bin2", "bin2", "bin2", "bin3", "bin3"
  )
  seq_ids <- c(
    "seq1", "seq1", "seq1", "seq2", "seq4", "seq4",
    "seq3", "seq3", "seq6", "seq5", "seq5"
  )
  samples <- c(
    "sample1", "sample2", "sample5",
    "sample1", "sample3", "sample4",
    "sample2", "sample3", "sample1",
    "sample1", "sample6"
  )
  abundances <- c(10, 100, 1, 500, 25, 80, 20, 5, 60, 15, 50)

  # all four
  assign_bins(dataset, bin_ids, abundances, samples, seq_ids)

  expect_equal(get_bin(dataset, "bin1"), "seq1,seq2,seq4")
  expect_equal(get_bin_abundance(dataset, "bin2"), 85)
  expect_equal(get_num_sequences(dataset), 866)
  expect_equal(get_num_bins(dataset), 3)
  expect_equal(get_num_samples(dataset), 6)

  clear(dataset, "")

  # assign sequence abundances then bins
  names <- c(
    "seq1", "seq1", "seq1", "seq2", "seq2",
    "seq2", "seq3", "seq3", "seq4"
  )
  samples <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  )
  abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)

  assign_sequence_abundance(dataset, names, abundances, samples, "")

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_bins(dataset), 0)
  expect_equal(get_num_samples(dataset), 3)

  # missing seq names
  expect_error(assign_bins(dataset, bin_ids, abundances, samples, seq_ids))

  bin_ids <- c(
    "bin1", "bin1", "bin1", "bin1",
    "bin1", "bin1", "bin1", "bin1", "bin2"
  )
  assign_bins(dataset, bin_ids, abundances, samples, names)

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_bins(dataset), 2)
  expect_equal(get_num_samples(dataset), 3)
})

test_that("rcpp_dataset - assign_sequence_abundance", {
  dataset <- new_dataset("miseq_sop", 4)

  # assign sequence abundances then bins
  names <- c(
    "seq1", "seq1", "seq1", "seq2", "seq2",
    "seq2", "seq3", "seq3", "seq4"
  )
  samples <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  )
  abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)

  expect_error(assign_sequence_abundance(dataset, names, 0, samples, ""))

  assign_sequence_abundance(dataset, names, abundances, samples, "")

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_bins(dataset), 0)
  expect_equal(get_num_samples(dataset), 3)

  abunds_by_sample <- get_sequence_abundances_by_sample(dataset)
  expect_equal(abunds_by_sample[[1]], c(250, 400, 500))
  expect_equal(abunds_by_sample[[3]], c(25, 25, 0))

  abund_table <- get_sequence_abundance_table(dataset)

  expect_equal(abund_table$sequence_names, names)
  expect_equal(abund_table$abundances, abundances)
  expect_equal(abund_table$samples, samples)

  message <- capture_output(set_abundances(
    dataset,
    c("non_existant_sequence"),
    list(c(0, 0, 0))
  ))

  expect_true(grepl("non_existant_sequence", message))
  expect_error(set_abundance(dataset, c("seq1"), c(0)))

  clear(dataset, "")

  names <- c("seq1", "seq2")
  abundances <- c(100, 50)

  assign_sequence_abundance(dataset, names, abundances, "", "")
  expect_equal(get_sequence_abundances(dataset), abundances)

  expect_error(set_abundances(dataset, c("seq1"), list(c(0))))

  message <- capture_output(set_abundance(
    dataset,
    c("non_existant_sequence"),
    c(0)
  ))
  expect_true(grepl("non_existant_sequence", message))
})

test_that("rcpp_dataset - assign_sequence_taxonomy", {
  dataset <- new_dataset("miseq_sop", 4)

  # assign sequence abundances then bins
  names <- c("seq1", "seq2", "seq3", "seq4")
  taxonomies <- c(
    "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;"
  )

  dataset <- new_dataset("my_dataset", 4)

  expect_error(assign_sequence_taxonomy(dataset, names, ""))
  assign_sequence_taxonomy(dataset, names, taxonomies)

  report <- get_sequence_taxonomy_report(dataset)

  expect_equal(report[[1]][[1]], "seq1")
  expect_equal(report[[3]][[1]], "Bacteria")
  expect_equal(report[[3]][[2]], "Bacteroidetes")
  expect_equal(report[[3]][[6]], "Proteobacteria")
  expect_equal(get_num_bins(dataset), 0)
  expect_equal(get_num_samples(dataset), 0)

  message <- capture_output(assign_sequence_taxonomy(
    dataset,
    c("non_existant_sequence"),
    c("bad_tax;")
  ))
  expect_true(grepl("non_existant_sequence", message))
})

test_that("rcpp_dataset - assign_treatments", {
  names <- c(
    "seq1", "seq1", "seq1", "seq2", "seq2",
    "seq2", "seq3", "seq3", "seq4"
  )
  samples <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  )
  abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)
  treatments <- c("early", "early", "late")

  dataset <- new_dataset("my_dataset", 4)

  # no samples in dataset
  expect_error(assign_treatments(dataset, unique(samples), treatments))

  # mismatched lengths
  expect_error(assign_treatments(dataset, samples, treatments))

  # add sample data
  assign_sequence_abundance(dataset, names, abundances, samples, "")

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_bins(dataset), 0)
  expect_equal(get_num_samples(dataset), 3)
  expect_equal(get_num_treatments(dataset), 0)

  missing_samples <- c("sample_not_in_dataset", "sample34", "sample_bad")

  expect_error(assign_treatments(dataset, missing_samples, treatments))

  assign_treatments(dataset, unique(samples), treatments)

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_bins(dataset), 0)
  expect_equal(get_num_samples(dataset), 3)
  expect_equal(get_num_treatments(dataset), 2)

  abunds_by_sample <- list(
    c(250, 400, 500), c(25, 40, 50),
    c(25, 25, 0), c(0, 0, 4)
  )
  expect_equal(get_sequence_abundances_by_sample(dataset), abunds_by_sample)
})

test_that("rcpp_dataset - get_names / sequences_by_sample", {
  dataset <- new_dataset("miseq_sop", 4)

  # assign sequence abundances then bins
  names <- c(
    "seq1", "seq1", "seq1", "seq2", "seq2",
    "seq2", "seq3", "seq3", "seq4"
  )
  samples <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  )
  abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)

  assign_sequence_abundance(dataset, names, abundances, samples, "")

  names_by_sample <- get_sequence_names_by_sample(dataset, "")
  expect_equal(length(names_by_sample), 3)

  # seqs in sample2
  expect_equal(names_by_sample[[1]], c("seq1", "seq2", "seq3"))
  # seqs in sample4
  expect_equal(names_by_sample[[3]], c("seq1", "seq2", "seq4"))

  seqs <- c("ATGCCT", "GTGCCT", "CTGCCT", "TTGCCT")
  set_sequences(dataset, unique(names), seqs, "")

  seqs_by_sample <- get_sequences_by_sample(dataset, "")
  expect_equal(length(seqs_by_sample), 3)

  # seqs in sample2
  expect_equal(seqs_by_sample[[1]], c("ATGCCT", "GTGCCT", "CTGCCT"))
  # seqs in sample4
  expect_equal(seqs_by_sample[[3]], c("ATGCCT", "GTGCCT", "TTGCCT"))

  names_by_sample <- get_sequence_names_by_sample(dataset, c("sample3"))

  # seqs in sample3
  expect_equal(names_by_sample[[1]], c("seq1", "seq2", "seq3"))
  expect_equal(length(names_by_sample), 1)

  seqs_by_sample <- get_sequences_by_sample(dataset, c("sample3"))

  # seqs in sample3
  expect_equal(seqs_by_sample[[1]], c("ATGCCT", "GTGCCT", "CTGCCT"))
  expect_equal(length(seqs_by_sample), 1)
})

test_that("rcpp_dataset - merge_bins / merge_seqs", {
  dataset <- new_dataset("miseq_sop", 4)

  # assign sequence abundances then bins
  names <- c(
    "seq1", "seq1", "seq1", "seq2", "seq2",
    "seq2", "seq3", "seq3", "seq4"
  )
  samples <- c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  )
  abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)

  assign_sequence_abundance(dataset, names, abundances, samples, "")

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_sequences(dataset, TRUE), 4)
  expect_equal(get_num_bins(dataset), 0)
  expect_equal(get_num_samples(dataset), 3)

  expect_equal(get_sequence_abundances(dataset), c(1150, 115, 50, 4))

  merge_sequences(dataset, c("seq3", "seq1"))

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_sequences(dataset, TRUE), 3)
  expect_equal(get_num_bins(dataset), 0)
  expect_equal(get_num_samples(dataset), 3)

  assign_bins(
    dataset, c("bin1", "bin1", "bin2"), 0, "",
    c("seq2", "seq3", "seq4")
  )

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_sequences(dataset, TRUE), 3)
  expect_equal(get_num_bins(dataset), 2)
  expect_equal(get_num_samples(dataset), 3)

  merge_bins(dataset, c("bin1", "bin2"))

  expect_equal(get_num_sequences(dataset), 1319)
  expect_equal(get_num_sequences(dataset, TRUE), 3)
  expect_equal(get_num_bins(dataset), 1)
  expect_equal(get_num_samples(dataset), 3)

  expect_error(merge_sequences(dataset, c("seq3", "non_existant_seq")))

  message <- capture_output(remove_sequences(
    dataset,
    c("non_existent_seq"),
    c("trash_tag")
  ))
  expect_true(grepl("non_existent_seq", message))

  message <- capture_output(remove_sequences(
    dataset,
    c("non_existent_seq"),
    c("trash_tag")
  ))
  expect_true(grepl("non_existent_seq", message))
})

test_that("rcpp_dataset - set_bin_abundance / set_bin_abundances, warnings", {
  dataset <- new_dataset("miseq_sop", 4)

  # test with no samples bins
  assign_bins(
    dataset, c("bin1", "bin2", "bin3"),
    c(10, 20, 30), "", ""
  )

  expect_equal(get_bin_abundance(dataset, "bin1"), 10)
  expect_equal(get_num_sequences(dataset), 60)

  message <- capture_output(set_bin_abundance(
    dataset,
    c("bin1", "non_existent_bin"),
    c(100, 20)
  ))
  # set bin1, ignore non_existent_bin
  expect_equal(get_bin_abundance(dataset, "bin1"), 100)
  expect_equal(get_num_sequences(dataset), 150)
  expect_true(grepl("non_existent_bin", message))

  expect_error(set_bin_abundances(
    dataset,
    c("bin1"), list(c(100, 20))
  ))

  clear(dataset, "")

  # test with samples bins
  assign_bins(
    dataset, c("bin1", "bin1", "bin2", "bin3", "bin3"),
    c(10, 20, 30, 40, 50), c(
      "sample1", "sample2", "sample1",
      "sample1", "sample2"
    ), ""
  )

  expect_equal(get_bin_abundances(dataset, "bin1"), c(10, 20))
  expect_equal(get_num_sequences(dataset), 150)
  expect_equal(get_sample_totals(dataset), c(80, 70))

  message <- capture_output(set_bin_abundances(
    dataset,
    c("bin1", "non_existent_bin"),
    list(c(100, 20), c(10, 0))
  ))
  expect_equal(get_num_sequences(dataset), 240)
  expect_equal(get_sample_totals(dataset), c(170, 70))
  expect_true(grepl("non_existent_bin", message))
  expect_error(set_bin_abundance(
    dataset,
    c("bin1"), c(100)
  ))

  # remove samples that will force bin2 removal
  remove_samples(dataset, c("sample1"))

  expect_equal(get_num_bins(dataset), 2)
  expect_equal(get_rabund_vector(dataset), c(20, 50))
  expect_equal(get_num_samples(dataset), 1)
  expect_equal(get_num_sequences(dataset), 70)

  expect_error(remove_bins(dataset, c("bin1"), c("trash_tag", "extra_one")))
})

test_that("rcpp_dataset - misc ", {
  dataset <- new_dataset("miseq_sop", 4)

  # no bin data
  expect_equal(get_bin_taxonomy_report(dataset), data.frame())

  assign_bins(
    dataset, c("bin1", "bin2", "bin3"),
    c(10, 20, 30), "", ""
  )

  # no bin taxonomies
  expect_equal(get_bin_taxonomy_report(dataset), data.frame())
})
