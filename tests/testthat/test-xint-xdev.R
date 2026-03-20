# adds coverage for internal (xint) and developer (xdev) Rcpp functions

test_that("xdev_abundance", {
  data <- new_dataset()

  # bad sample
  message <- capture_output(xdev_abundance(data, type = "bad_type"))
  expect_true(grepl(
    "bad_type is not a valid type for the abundance function",
    message
  ))
})

test_that("xdev_assign_bin_taxonomy", {
  data <- new_dataset()

  # can't assign bin taxonomy when you have not assigned bins
  expect_error(xdev_assign_bin_taxonomy(data, data.frame()))
})

test_that("xdev_assign_sequence_abundance", {
  data <- new_dataset()

  add(data, data.frame(sequence_names = c("seq1")))

  table <- data.frame(
    sequence_names = c("seq2", "seq3"),
    abundances = c(100, 200)
  )

  # must assign abundances for all sequences
  expect_error(xdev_assign_sequence_abundance(data, table))
})

test_that("xdev_get_by_sample", {
  data <- new_dataset()

  # can't assign bin taxonomy when you have not assigned bins
  expect_error(xdev_get_by_sample(data, "badType"))
})

test_that("xdev_get_abundances_by_sample - getSequenceAbundanceBySample", {
  data <- new_dataset()

  seq_names <- c(
    "seq1", "seq2", "seq3", "seq3",
    "seq4", "seq4", "seq5", "seq6",
    "seq7", "seq8", "seq9", "seq9",
    "seq10", "seq10", "seq10", "seq10"
  )
  samples <- c(
    "sample1", "sample2", "sample4", "sample5",
    "sample1", "sample2", "sample1", "sample1",
    "sample2", "sample4", "sample4", "sample5",
    "sample1", "sample3", "sample5", "sample6"
  )
  abundances <- c(
    10, 10, 5, 5, 5, 5,
    10, 10, 10, 10, 5, 5,
    1, 2, 3, 4
  )

  assign(
    data = data,
    table = data.frame(
      sequence_names = seq_names,
      samples = samples,
      abundances = abundances
    ),
    type = "sequence_abundance"
  )

  expected <- list(
    c(10, 1, 5, 10, 10), c(10, 5, 10),
    c(2), c(5, 10, 5), c(3, 5, 5), c(4)
  )

  actual <- xdev_get_abundances_by_sample(data)

  expect_equal(actual, expected)
})


test_that("xdev_assign_bins, assign with reference", {
  data <- new_dataset()

  # bad sample
  reference <- new_reference("myReference")

  table <- data.frame(
    bin_names = c("bin1", "bin2", "bin2"),
    sequence_names = c("seq1", "seq2", "seq3")
  )

  xdev_assign_bins(data, table, reference = reference)

  ref_report <- report(data, type = "references")

  expect_equal(nrow(ref_report), 1)

  xdev_assign_bin_representative_sequences(
    data,
    data.frame(
      bin_names = c("bin1", "bin2"),
      sequence_names = c("seq1", "seq2")
    ),
    reference = reference
  )

  ref_report <- report(data, type = "references")

  expect_equal(nrow(ref_report), 2)
})

test_that("xdev_set_abundances", {
  data <- new_dataset()

  otus <- c(
    "otu1", "otu1", "otu1",
    "otu2", "otu2", "otu3",
    "otu4", "otu4", "otu4", "otu4"
  )
  seqs <- paste0("seq", 1:10)
  abunds <- rep(10, 10)

  # assign seqs to bins
  assign(data = data, table = data.frame(
    sequence_names = seqs,
    abundances = abunds,
    bin_names = otus
  ), type = "bins")

  # all seqs
  expect_equal(count(data, type = "sequences"), 100)
  # unique seqs
  expect_equal(count(data, type = "sequences", distinct = TRUE), 10)
  # numBins
  expect_equal(count(data, type = "bins"), 4)

  list <- xdev_get_list_vector(data, "otu")

  expect_equal(list[1], "seq1,seq2,seq3")
  expect_equal(list[2], "seq4,seq5")
  expect_equal(list[3], "seq6")
  expect_equal(list[4], "seq10,seq7,seq8,seq9")

  bin_abundances <- abundance(data, "bins")
  expected <- c(30, 20, 10, 40)

  expect_equal(bin_abundances[[2]], expected)
  expect_equal(length(bin_abundances[[2]]), 4)

  seqs_to_change <- c("seq1", "seq2", "seq7", "seq8", "seq9", "seq10")
  abunds_to_change <- c(40, 20, 0, 0, 0, 0)

  xdev_set_abundance(data, seqs_to_change, abunds_to_change)

  expect_equal(count(data, type = "sequences"), 100)
  expect_equal(count(data, type = "bins"), 3)
  expect_equal(count(data, type = "sequences", distinct = TRUE), 6)

  bin_abundances <- abundance(data, "bins")
  expected <- c(70, 20, 10)

  expect_equal(bin_abundances[[2]], expected)
  expect_equal(length(bin_abundances[[2]]), 3)
})

test_that("xdev_set_abundances", {
  data <- new_dataset()

  # bins -> seq1,seq2,seq4. seq3,seq6  seq5
  # binAbunds -> (111+500+105) + (25+60) + (65)
  # binSamples ->(1,2,3,4,5) + (1,2,3) + (1,6)
  # sampleTotals -> 1(10+500+60+15) + 2(100+20) + 3(25+5) + 4(80) + 5(1) + 6(50)
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

  xdev_assign_bins(data, data.frame(
    bin_names = bin_ids,
    abundances = abundances,
    samples = samples,
    sequence_names = seq_ids
  ))

  # "seq1" "seq2" "seq5" "seq6"
  sample1_names <- xdev_get_by_sample(
    data = data,
    type = "sequence_names",
    samples = c("sample1")
  )
  expect_equal(sample1_names[[1]], c("seq1", "seq2", "seq5", "seq6"))

  # 6
  num_samples <- count(data = data, type = "samples")
  expect_equal(num_samples, 6)

  new_abunds <- rep(list(rep(0, num_samples)), length(sample1_names[[1]]))

  # remove 4 seqs, 1 bin
  xdev_set_abundances(
    data = data,
    sequence_names = sample1_names[[1]],
    abundances = new_abunds
  )

  # bins -> seq4 seq3
  # binAbunds -> (105) + (25)
  # binSamples ->(3,4) + (2,3)
  # sampleTotals -> 1(0) + 2(20) + 3(25+5) + 4(80) + 5(0) + 6(0)
  expect_equal(count(data = data, type = "samples"), 3)
  expect_equal(count(data = data, type = "sequences"), 130)
  expect_equal(count(data = data, type = "sequences", distinct = TRUE), 2)
  expect_equal(abundance(data = data, type = "samples")[[2]], c(20, 30, 80))
})

test_that("Tests removeBins, getScrapReport, getScrapSummary", {
  otu_names <- c(
    "otu1", "otu1", "otu1", "otu1",
    "otu2", "otu2", "otu2", "otu2",
    "otu2", "otu2", "otu2", "otu2",
    "otu3", "otu3", "otu3", "otu3"
  )
  seq_names <- c(
    "seq1", "seq2", "seq3", "seq3",
    "seq4", "seq4", "seq5", "seq6",
    "seq7", "seq8", "seq9", "seq9",
    "seq10", "seq10", "seq10", "seq10"
  )
  samples <- c(
    "sample1", "sample2", "sample4", "sample5",
    "sample1", "sample2", "sample1", "sample1",
    "sample2", "sample4", "sample4", "sample5",
    "sample1", "sample3", "sample5", "sample6"
  )
  abundances <- c(
    10, 10, 5, 5, 5, 5,
    10, 10, 10, 10, 5, 5,
    1, 2, 3, 4
  )
  treatments <- rep("early", 6)

  data <- new_dataset("mydata")

  assign(
    data = data,
    table = data.frame(
      bin_names = otu_names,
      sequence_names = seq_names,
      samples = samples,
      abundances = abundances
    ),
    type = "bins", bin_type = "otu"
  )

  assign(data,
    table = data.frame(
      samples = unique(samples),
      treatments = treatments
    ),
    type = "treatments"
  )

  bin_abundances <- abundance(data, type = "bins", bin_type = "otu")
  sample_totals <- abundance(data, type = "samples")

  expect_equal(bin_abundances[[2]], c(30, 60, 10))
  expect_equal(sample_totals[[2]], c(36, 25, 2, 20, 13, 4))
  expect_equal(count(data, "bins"), 3)
  expect_equal(abundance(data, "treatments")[[2]], c(100))
  expect_equal(count(data, "treatments"), 1)

  # set last 3 samples treatment assignment to late
  treatments[(length(treatments) - 2):length(treatments)] <- "late"

  assign(data,
    table = data.frame(
      samples = sort(unique(samples)),
      treatments = treatments
    ),
    type = "treatments"
  )

  expect_equal(abundance(data, "treatments")[[2]], c(63, 37))
  expect_equal(count(data, "treatments"), 2)

  # attempt bad entry
  unique_samples <- c(unique(samples), "SampleNotInDataset")
  treatments <- c(treatments, "badEntry")

  # should error
  expect_error(assign(data,
    table = data.frame(
      samples = unique_samples,
      treatments = treatments
    ),
    type = "treatments"
  ))

  expect_equal(abundance(data, "treatments")[[2]], c(63, 37))
  expect_equal(count(data, "treatments"), 2)

  expect_equal(count(data), 100)
  expect_equal(count(data, distinct = TRUE), 10)
  expect_equal(sample_totals[[2]], c(36, 25, 2, 20, 13, 4))
  expect_equal(count(data, "bins"), 3)
  expect_equal(count(data, "treatments"), 2)

  expect_equal(xdev_get_list_vector(data), c(
    "seq1,seq2,seq3",
    "seq4,seq5,seq6,seq7,seq8,seq9",
    "seq10"
  ))

  # all otu abunds by sample
  expect_equal(
    xdev_abundance(data, "bins", "otu", TRUE)[[2]],
    c(10, 10, 5, 5, 25, 15, 15, 5, 1, 2, 3, 4)
  )

  # remove otu1 and test bad otu name
  otus_to_remove <- c("otu1", "non_existant_otu")
  reasons <- rep("bad_bin", 2)

  xdev_remove_bins(data, otus_to_remove, reasons, "otu")

  expect_equal(
    abundance(data, type = "bins", bin_type = "otu")[[2]],
    c(60, 10)
  )
  expect_equal(sample_totals[[2]], c(36, 25, 2, 20, 13, 4))
  expect_equal(count(data), 70)
  expect_equal(count(data, distinct = TRUE), 7)
  expect_equal(count(data, "bins"), 2)
  expect_equal(count(data, "treatments"), 2)
  expect_equal(count(data, "samples"), 6)

  # get scrap report
  otu_scrap_report <- report(data = data, type = "bin_scrap")

  expect_equal(otu_scrap_report[[1]], c("otu1"))
  expect_equal(otu_scrap_report[[2]], c("bad_bin"))

  scrap_summary <- summary(data = data, type = "scrap")

  expect_equal(scrap_summary[[1]], c("sequence", "otu"))
  expect_equal(scrap_summary[[2]], c("bad_bin", "bad_bin"))
  expect_equal(scrap_summary[[3]], c(3, 1))
  expect_equal(scrap_summary[[4]], c(30, 30))
})


test_that("Tests assignSequenceTaxonomy, assignBinTaxonomy, removeLineages", {
  names <- c("seq1", "seq2", "seq3", "seq4")
  otus <- c("otu1", "otu1", "otu2", "otu2")
  abunds <- c(100, 10, 10, 5)
  taxonomies <- c(
    "Bacteria(100);Bacteroidetes(95);Bacteroidia(90);",
    "Bacteria(100);Proteobacteria(89);Betaproteobacteria(85);",
    "Bacteria(100);Firmicutes(99);Bacilli(90);",
    "Bacteria(100);Firmicutes(87);Bacilli(85);"
  )

  data <- new_dataset("testdata")

  # assign sequence taxonomy
  assign(
    data = data,
    table = data.frame(
      sequence_names = names,
      taxonomies = taxonomies
    ),
    type = "sequence_taxonomy"
  )

  # assign sequence abundance
  assign(
    data = data,
    table = data.frame(
      sequence_names = names,
      abundances = abunds
    ),
    type = "sequence_abundance"
  )


  # assign otus, and check otu classifications
  assign(
    data = data,
    table = data.frame(
      sequence_names = names,
      bin_names = otus
    ),
    type = "bins"
  )

  sequence_taxonomy_report <- report(data, "sequence_taxonomy")

  expect_equal(sequence_taxonomy_report[[1]], c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2",
    "seq3", "seq3", "seq3",
    "seq4", "seq4", "seq4"
  ))

  expect_equal(sequence_taxonomy_report[[3]], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia", "Bacteria",
    "Proteobacteria",
    "Betaproteobacteria",
    "Bacteria", "Firmicutes",
    "Bacilli", "Bacteria",
    "Firmicutes", "Bacilli"
  ))

  expect_equal(sequence_taxonomy_report[[4]], c(
    100, 95, 90, 100, 89, 85,
    100, 99, 90, 100, 87, 85
  ))


  bin_taxonomy_report <- report(data, "bin_taxonomy")

  expect_equal(bin_taxonomy_report[[1]], c(
    "otu1", "otu1", "otu1",
    "otu2", "otu2", "otu2"
  ))

  expect_equal(bin_taxonomy_report[[3]], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia",
    "Bacteria", "Firmicutes",
    "Bacilli"
  ))

  expect_equal(bin_taxonomy_report[[4]], c(100, 91, 91, 100, 100, 100))

  # remove "Bacteria(100);Firmicutes(90);Bacilli(85);"
  # the confidence thresholds should remove seq4 and leave seq3

  expect_equal(count(data, distinct = TRUE), 4)
  expect_equal(count(data), 125)
  expect_equal(count(data, "bins"), 2)

  contaminants <- c("Bacteria(100);Firmicutes(90);Bacilli(85);")

  xdev_remove_lineages(data, contaminants, "contaminant")

  expect_equal(count(data, distinct = TRUE), 3)
  expect_equal(count(data), 120)
  expect_equal(count(data, "bins"), 2)

  # bin_taxonomy_report should be the same, but sequence_taxonomy_report
  # should reflect seq4 removal
  bin_taxonomy_report <- report(data, "bin_taxonomy")
  sequence_taxonomy_report <- report(data, "sequence_taxonomy")

  expect_equal(bin_taxonomy_report[[1]], c(
    "otu1", "otu1", "otu1",
    "otu2", "otu2", "otu2"
  ))

  expect_equal(bin_taxonomy_report[[3]], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia",
    "Bacteria", "Firmicutes",
    "Bacilli"
  ))

  expect_equal(bin_taxonomy_report[[4]], c(100, 91, 91, 100, 100, 100))


  expect_equal(sequence_taxonomy_report[[1]], c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2",
    "seq3", "seq3", "seq3"
  ))

  expect_equal(sequence_taxonomy_report[[3]], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia", "Bacteria",
    "Proteobacteria",
    "Betaproteobacteria",
    "Bacteria", "Firmicutes",
    "Bacilli"
  ))

  expect_equal(sequence_taxonomy_report[[4]], c(
    100, 95, 90, 100, 89, 85,
    100, 99, 90
  ))

  # will remove seq3 and otu2
  contaminants <- c("Firmicutes")
  xdev_remove_lineages(data, contaminants, "contaminant")

  expect_equal(count(data, distinct = TRUE), 2)
  expect_equal(count(data), 110)
  expect_equal(count(data, "bins"), 1)

  expect_equal(xdev_get_list_vector(data), c("seq1,seq2"))

  bin_taxonomy_report <- report(data, "bin_taxonomy")

  expect_equal(bin_taxonomy_report[[1]], c("otu1", "otu1", "otu1"))

  expect_equal(bin_taxonomy_report[[3]], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia"
  ))

  expect_equal(bin_taxonomy_report[[4]], c(100, 91, 91))

  sequence_taxonomy_report <- report(data, "sequence_taxonomy")

  expect_equal(sequence_taxonomy_report[[1]], c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2"
  ))

  expect_equal(sequence_taxonomy_report[[3]], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia", "Bacteria",
    "Proteobacteria",
    "Betaproteobacteria"
  ))

  expect_equal(sequence_taxonomy_report[[4]], c(100, 95, 90, 100, 89, 85))

  clear(data)

  # tests remove.lineage with only bin tax assignments
  otus <- c("otu1", "otu2")
  abunds <- c(100, 500)
  taxonomies <- c(
    "Bacteria;Bacteroidetes;Bacteroidia;",
    "Bacteria;Proteobacteria;Betaproteobacteria;"
  )

  # assign otus,
  assign(
    data = data,
    table = data.frame(
      abundances = abunds,
      bin_names = otus
    ),
    type = "bins"
  )

  # assign bin taxonomies
  assign(
    data = data,
    table = data.frame(
      taxonomies = taxonomies,
      bin_names = otus
    ),
    type = "bin_taxonomy"
  )

  bin_taxonomy_report <- report(data, "bin_taxonomy")

  expect_equal(count(data), 600)
  expect_equal(count(data, "bins"), 2)

  expect_equal(bin_taxonomy_report[[3]], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia",
    "Bacteria", "Proteobacteria",
    "Betaproteobacteria"
  ))

  contaminants <- c("Proteobacteria")
  xdev_remove_lineages(data, contaminants, "contaminant")

  expect_equal(count(data), 100)
  expect_equal(count(data, "bins"), 1)
})
