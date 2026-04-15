# test "rcpp_dataset"

test_that("rcpp_dataset - new_dataset, copy_dataset", {
  data <- new_dataset("miseq_sop", 2)
  expect_equal(names(data, "dataset"), "miseq_sop")
  expect_equal(count(data, "sequences"), 0)

  dataset2 <- copy_dataset(data)
  expect_equal(names(dataset2, "dataset"), "miseq_sop")
  expect_equal(count(dataset2, "sequences"), 0)

  xdev_add_sequences(data, data.frame(sequence_names = c("seq1")))
  expect_equal(count(data, "sequences"), 1)
  expect_equal(count(dataset2, "sequences"), 0)
})

test_that("rcpp_dataset - add_sequences", {
  data <- new_dataset("miseq_sop", 2)

  sequences <- read_fasta(strollur_example("final.fasta.gz"))
  sequences$comments <- rep("my comments", length(sequences$sequence_names))

  # add with all 3 parameters
  xdev_add_sequences(data, sequences)

  expect_equal(count(data, "sequences"), 2425)
  expect_true(has_sequence_strings(data))

  clear(data)

  # add with 2 parameters
  xdev_add_sequences(data, sequences)
  expect_equal(count(data, "sequences"), 2425)
  expect_true(has_sequence_strings(data))

  seq_report <- report(data, "sequences")

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

  clear(data)

  # add with only names
  xdev_add_sequences(
    data,
    data.frame(sequence_names = c("seq1", "seq2", "seq3"))
  )
  expect_equal(count(data, "sequences"), 3)
  expect_false(has_sequence_strings(data))

  clear(data)

  expect_error(xdev_add_sequences(data))

  message <- capture_output(xdev_set_sequences(
    data,
    c("non_existant_sequence"),
    c("ATGC"), ""
  ))
  expect_true(grepl("non_existant_sequence", message))
})

test_that("rcpp_dataset - assign_bins", {
  data <- new_dataset("miseq_sop", 2)

  # bin labels and seq names (list)
  seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  xdev_assign_bins(data, data.frame(
    bin_names = bin_ids,
    sequence_names = seq_ids
  ))

  expect_equal(count(data, "sequences"), 6)
  expect_equal(count(data, "bins"), 3)
  expect_false(has_sequence_strings(data))

  abundances <- c(110, 525, 80)
  samples <- c("sample1", "sample1", "sample2")
  expect_error(xdev_assign_bins(data, data.frame(
    bin_names = unique(bin_ids),
    abundances = abundances
  )))
  expect_error(xdev_assign_bins(data, data.frame(
    bin_names = unique(bin_ids),
    samples = samples
  )))
  clear(data)

  # bin labels and bin abundances (rabund)
  assign(data = data, table = data.frame(
    bin_names = unique(bin_ids),
    abundances = abundances
  ), type = "bins", bin_type = "otu")

  expect_equal(count(data, "sequences"), 715)
  expect_equal(count(data, "bins"), 3)

  clear(data)

  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  samples <- c(
    "sample1", "sample2", "sample5", "sample1", "sample3",
    "sample1"
  )
  sample_abundances <- c(10, 100, 1, 500, 25, 80)

  expect_error(xdev_assign_bins(data, data.frame(
    abundances = sample_abundances,
    samples = samples
  )))

  # bin labels, samples and bin abundances (shared)
  xdev_assign_bins(data, data.frame(
    bin_names = bin_ids,
    abundances = sample_abundances,
    samples = samples
  ))

  expect_equal(count(data, "sequences"), 716)
  expect_equal(count(data, "bins"), 3)
  expect_equal(count(data, "samples"), 4)
  expect_equal(abundance(data, "samples")$abundances, c(590, 100, 25, 1))

  df <- abundance(
    data = data, type = "bins",
    bin_type = "otu", by_sample = TRUE
  )
  expect_equal(df[[2]][1:3], c(10, 100, 1))

  clear(data)

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
  xdev_assign_bins(data, data.frame(
    bin_names = bin_ids,
    abundances = abundances,
    samples = samples,
    sequence_names = seq_ids
  ))


  expect_equal(xdev_get_list_vector(data)[1], "seq1,seq2,seq4")
  expect_equal(abundance(
    data = data, type = "bins",
    bin_type = "otu"
  )[[2]][[2]], 85)
  expect_equal(count(data, "sequences"), 866)
  expect_equal(count(data, "bins"), 3)
  expect_equal(count(data, "samples"), 6)

  clear(data)

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

  xdev_assign_sequence_abundance(data, data.frame(
    sequence_names = names,
    abundances = abundances,
    samples = samples
  ))

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, "bins"), 0)
  expect_equal(count(data, "samples"), 3)

  # missing seq names
  expect_error(xdev_assign_bins(data, data.frame(
    bin_names = bin_ids,
    abundances = abundances,
    samples = samples,
    sequence_names = seq_ids
  )))

  bin_ids <- c(
    "bin1", "bin1", "bin1", "bin1",
    "bin1", "bin1", "bin1", "bin1", "bin2"
  )
  xdev_assign_bins(data, data.frame(
    bin_names = bin_ids,
    abundances = abundances,
    samples = samples,
    sequence_names = names
  ))

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, "bins"), 2)
  expect_equal(count(data, "samples"), 3)
})

test_that("rcpp_dataset - assign_sequence_abundance", {
  data <- new_dataset("miseq_sop", 2)

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

  expect_error(xdev_assign_sequence_abundance(
    data,
    data.frame(
      sequence_names = names,
      samples = samples
    )
  ))

  xdev_assign_sequence_abundance(data, data.frame(
    sequence_names = names,
    abundances = abundances,
    samples = samples
  ))

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, "bins"), 0)
  expect_equal(count(data, "samples"), 3)

  abunds_by_sample <- abundance(data, type = "sequences", by_sample = TRUE)
  expect_equal(abunds_by_sample[[2]][1:3], c(250, 400, 500))
  expect_equal(abunds_by_sample[[2]][7:8], c(25, 25))

  abund_table <- abundance(data, type = "sequences", by_sample = TRUE)

  expect_equal(abund_table$sequence_names, names)
  expect_equal(abund_table$abundances, abundances)
  expect_equal(abund_table$samples, samples)

  message <- capture_output(xdev_set_abundances(
    data,
    c("non_existant_sequence"),
    list(c(0, 0, 0))
  ))

  expect_true(grepl("non_existant_sequence", message))
  expect_error(xdev_set_abundance(data, c("seq1"), c(0)))

  clear(data)

  xdev_assign_sequence_abundance(
    data,
    data.frame(
      sequence_names = c("seq1", "seq2"),
      abundances = c(100, 50)
    )
  )

  expect_equal(abundance(data, type = "sequences")[[2]], c(100, 50))

  expect_error(xdev_set_abundances(data, c("seq1"), list(c(0))))

  message <- capture_output(xdev_set_abundance(
    data,
    c("non_existant_sequence"),
    c(0)
  ))
  expect_true(grepl("non_existant_sequence", message))
})

test_that("rcpp_dataset - assign_sequence_taxonomy", {
  data <- new_dataset("miseq_sop", 2)

  # assign sequence abundances then bins
  names <- c("seq1", "seq2", "seq3", "seq4")
  taxonomies <- c(
    "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;"
  )

  expect_error(xdev_assign_sequence_taxonomy(data, names, ""))
  xdev_assign_sequence_taxonomy(data, data.frame(
    sequence_names = names,
    taxonomies = taxonomies
  ))

  report <- report(data, "sequence_taxonomy")

  expect_equal(report[[1]][[1]], "seq1")
  expect_equal(report[[3]][[1]], "Bacteria")
  expect_equal(report[[3]][[2]], "Bacteroidetes")
  expect_equal(report[[3]][[6]], "Proteobacteria")
  expect_equal(count(data, "bins"), 0)
  expect_equal(count(data, "samples"), 0)

  message <- capture_output(xdev_assign_sequence_taxonomy(
    data,
    data.frame(
      sequence_names = c("non_existant_sequence"),
      taxonomies = c("bad_tax;")
    )
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

  data <- new_dataset("my_dataset", 2)

  # no samples in dataset
  expect_error(xdev_assign_treatments(data, data.frame(
    samples = unique(samples),
    treatments = treatments
  )))

  # add sample data
  xdev_assign_sequence_abundance(data, data.frame(
    sequence_names = names,
    abundances = abundances,
    samples = samples
  ))

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, "bins"), 0)
  expect_equal(count(data, "samples"), 3)
  expect_equal(count(data, "treatments"), 0)

  missing_samples <- c("sample_not_in_dataset", "sample34", "sample_bad")

  expect_error(xdev_assign_treatments(data, data.frame(
    samples = missing_samples,
    treatments = treatments
  )))

  assign(data = data, table = data.frame(
    samples = unique(samples),
    treatments = treatments
  ), type = "treatments")

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, "bins"), 0)
  expect_equal(count(data, "samples"), 3)
  expect_equal(count(data, "treatments"), 2)
  expect_equal(
    abundance(data, "sequences", by_sample = TRUE)[[2]],
    c(250, 400, 500, 25, 40, 50, 25, 25, 4)
  )
})

test_that("rcpp_dataset - get_names / sequences_by_sample", {
  data <- new_dataset("miseq_sop", 2)

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

  xdev_assign_sequence_abundance(data, data.frame(
    sequence_names = names,
    abundances = abundances,
    samples = samples
  ))

  names_by_sample <- xdev_get_by_sample(data)
  expect_equal(length(names_by_sample), 3)

  # seqs in sample2
  expect_equal(names_by_sample[[1]], c("seq1", "seq2", "seq3"))
  # seqs in sample4
  expect_equal(names_by_sample[[3]], c("seq1", "seq2", "seq4"))

  seqs <- c("ATGCCT", "GTGCCT", "CTGCCT", "TTGCCT")
  xdev_set_sequences(data, unique(names), seqs)

  seqs_by_sample <- xdev_get_by_sample(data, "sequences")
  expect_equal(length(seqs_by_sample), 3)

  # seqs in sample2
  expect_equal(seqs_by_sample[[1]], c("ATGCCT", "GTGCCT", "CTGCCT"))
  # seqs in sample4
  expect_equal(seqs_by_sample[[3]], c("ATGCCT", "GTGCCT", "TTGCCT"))

  names_by_sample <- xdev_get_by_sample(data, "sequence_names", c("sample3"))

  # seqs in sample3
  expect_equal(names_by_sample[[1]], c("seq1", "seq2", "seq3"))
  expect_equal(length(names_by_sample), 1)

  seqs_by_sample <- xdev_get_by_sample(data, "sequences", c("sample3"))

  # seqs in sample3
  expect_equal(seqs_by_sample[[1]], c("ATGCCT", "GTGCCT", "CTGCCT"))
  expect_equal(length(seqs_by_sample), 1)
})

test_that("rcpp_dataset - xdev_merge_bins / xdev_merge_seqs", {
  data <- new_dataset("miseq_sop", 2)

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

  xdev_assign_sequence_abundance(data, data.frame(
    sequence_names = names,
    abundances = abundances,
    samples = samples
  ))

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, "sequences", "otu", NULL, TRUE), 4)
  expect_equal(count(data, "bins"), 0)
  expect_equal(count(data, "samples"), 3)

  expect_equal(abundance(data, type = "sequences")[[2]], c(1150, 115, 50, 4))

  xdev_merge_sequences(data, c("seq3", "seq1"))

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, type = "sequences", distinct = TRUE), 3)
  expect_equal(count(data, "bins"), 0)
  expect_equal(count(data, "samples"), 3)

  xdev_assign_bins(data, data.frame(
    bin_names = c("bin1", "bin1", "bin2"),
    sequence_names = c("seq2", "seq3", "seq4")
  ))

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, "sequences", distinct = TRUE), 3)
  expect_equal(count(data, "bins"), 2)
  expect_equal(count(data, "samples"), 3)

  xdev_merge_bins(data, c("bin1", "bin2"))

  expect_equal(count(data, "sequences"), 1319)
  expect_equal(count(data, "sequences", distinct = TRUE), 3)
  expect_equal(count(data, "bins"), 1)
  expect_equal(count(data, "samples"), 3)

  expect_error(xdev_merge_sequences(data, c("seq3", "non_existant_seq")))

  message <- capture_output(xdev_remove_sequences(
    data,
    c("non_existent_seq"),
    c("trash_tag")
  ))
  expect_true(grepl("non_existent_seq", message))

  message <- capture_output(xdev_remove_sequences(
    data,
    c("non_existent_seq"),
    c("trash_tag")
  ))
  expect_true(grepl("non_existent_seq", message))
})

test_that("rcpp_dataset - misc ", {
  data <- new_dataset("miseq_sop", 2)

  # no bin data
  expect_equal(report(data, "bin_taxonomy"), data.frame())

  xdev_assign_bins(
    data, data.frame(
      bin_names = c("bin1", "bin2", "bin3"),
      abundances = c(10, 20, 30)
    )
  )

  # no bin taxonomies
  expect_equal(report(data, "bin_taxonomy"), data.frame())
})

test_that("dataset - xdev_set_abundances, xdev_set_sequences", {
  # create dataset sequences and shared data
  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    dataset_name = "miseq_sop"
  )

  expect_equal(count(dataset_t, "sequences"), 113963)

  expect_error(xdev_set_sequences(
    dataset_t$data,
    rep("seq1", 5),
    rep("ATGC", 5),
    c("not", "enough", "sequence", "comments")
  ))

  # abund = 191
  seqs_to_update <- c("M00967_43_000000000-A3JHG_1_1108_14299_17220")
})

test_that("Tests assignSequenceTaxonomy forces reclassify", {
  data <- new_dataset()

  otus <- c("otu1", "otu2", "otu3", "otu3")
  seqs <- c("seq1", "seq2", "seq3", "seq4")
  taxs <- c(
    "Bacteria(90);",
    "Bacteria(100);",
    "Bacteria(100);Bacteroidetes(95);Bacteroidia(90);",
    "Bacteria(100);Bacteroidetes(93);"
  )
  taxs2 <- c(
    "Bacteria(100);Proteobacteria(89);Betaproteobacteria(85);",
    "Bacteria(100);Firmicutes(99);Bacilli(90);",
    "Bacteria(100);Proteobacteria(90);",
    "Bacteria(100);Bacteroidetes(93);"
  )

  assign(
    data,
    data.frame(sequence_names = seqs, taxonomies = taxs),
    "sequence_taxonomy"
  )
  assign(
    data,
    data.frame(sequence_names = seqs, bin_names = otus),
    "bins"
  )

  report <- report(data, "bin_taxonomy")

  expected <- c(
    "Bacteria", "Bacteria_unclassified", "Bacteria_unclassified",
    "Bacteria", "Bacteria_unclassified", "Bacteria_unclassified",
    "Bacteria", "Bacteroidetes", "Bacteroidia"
  )
  expect_equal(report[["taxonomies"]], expected)

  assign(
    data,
    data.frame(sequence_names = seqs, taxonomies = taxs2),
    "sequence_taxonomy"
  )

  report <- report(data, "bin_taxonomy")
  expected <- c(
    "Bacteria", "Proteobacteria", "Betaproteobacteria",
    "Bacteria", "Firmicutes", "Bacilli",
    "Bacteria", "Bacteroidetes", "Bacteroidetes_unclassified"
  )

  expect_equal(report[["taxonomies"]], expected)
})
