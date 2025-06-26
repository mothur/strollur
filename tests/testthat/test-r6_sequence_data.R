# test "sequence_data"

test_that("test R6 sequence_data - intialize", {
  dataset <- sequence_data$new("mydata")

  expect_equal(dataset$get_dataset_name(), "mydata")

  dataset <- sequence_data$new(
    name = "soil",
    fasta = rdataset_example("test.fasta")
  )

  first_ten <- c(
    "M00967_43_000000000-A3JHG_1_1111_8697_7063",
    "M00967_43_000000000-A3JHG_1_1107_13334_19316",
    "M00967_43_000000000-A3JHG_1_2103_15942_24856",
    "M00967_43_000000000-A3JHG_1_2110_12430_18520",
    "M00967_43_000000000-A3JHG_1_1107_22778_21712",
    "M00967_43_000000000-A3JHG_1_1103_18837_8349",
    "M00967_43_000000000-A3JHG_1_1112_11891_19615",
    "M00967_43_000000000-A3JHG_1_1107_14687_22683",
    "M00967_43_000000000-A3JHG_1_1111_18897_11582",
    "M00967_43_000000000-A3JHG_1_1112_12027_22976"
  )

  expect_equal(first_ten, dataset$get_ids()[1:10])

  seq1 <- paste("TACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCGCGCAGGTGGTTAATT",
    "AAGTCTGATGTGAAAGCCCACGGCTTAACCGTGGAGGGTCATTGGAAACTGGTTGACTTGA",
    "GTGCAGAAGAGGGAAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGAGATATGGAGGAA",
    "CACCAGTGGCGAAGGCGGCTTCCTGGTCTGCAACTGACACTGAGGCGCGAAAGCGTGGGGA",
    "GCAAACAGG",
    sep = "", collaspe = ""
  )
  seq2 <- paste("TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGTC",
    "AAGTCAGCGGTAAAATTGTGGAGCTCAACTCCATCGAGCCGTTGAAACTGACGTTCTTGAG",
    "TGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAAC",
    "TCCGATTGCGAAGGCAGCATACCGGCGCCCTACTGACGCTGAAGCACGAAAGCGTGGGTAT",
    "CGAACAGG",
    sep = "", collaspe = ""
  )
  seq3 <- paste("TACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGGACGC",
    "AAGTCAGATGTGAAAACCACGGGCTCAACCTGTGGCCTGCATTTGAAACTGTGTTTCTTGA",
    "GTACTGGAGAGGCAGACGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAA",
    "CACCAGTGGCGAAGGCGGTCTGCTGGACAGCAACTGACGCTGAGGCGCGAAAGCGTGGGGA",
    "GCAAACAGG",
    sep = "", collaspe = ""
  )
  seq4 <- paste("TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGTGATGC",
    "AAGTCTGAAGTGAAAGGCGGGGGCTCAACCCCCGGACTGCTTTGGAAACTGTATGACTGGA",
    "GTGCAGGAGAGGTAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAA",
    "CACCAGTGGCGAAGGCGGCTTACTGGACTGTAACTGACGTTGAGGCTCGAAAGCGTGGGGA",
    "GCAAACAGG",
    sep = "", collaspe = ""
  )
  seq5 <- paste("TACGTAGGTGGCGAGCGTTGTCCGGATTTACTGGGCGTAAAGGGAGCGTAGGCGGACTTTT",
    "AAGTGAGATGTGAAATACTCGGGCTCAACTTGAGTGCTGCATTTCAAACTGGAAGTCTAGA",
    "GTGCAGGAGAGGAGAATGGAATTCCTAGTGTAGCGGTGAAATGCGTAGAGATTAGGAAGAA",
    "CACCAGTGGCGAAGGCGATTCTCTGGACTGTAACTGACGCTGAGGCTCGAAAGCGTGGGGA",
    "GCAAACAGG",
    sep = "", collaspe = ""
  )
  seq6 <- paste("TACGGAGGNTGGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGG",
    "TTAAGTCAGTGGTCAAATTGAGGGGCTCAACCCCTTCCCGCCATTGAAACTGGCGATCTTG",
    "AGTGGAAGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGA",
    "ACCCCGATTGCGAAGGCAGCATGCCGGCTTCCTACTGACGCTGAAGCACGAAAGCGTGGGT",
    "ATCGAACAGG",
    sep = "", collaspe = ""
  )
  seq7 <- paste("TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATC",
    "AAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAG",
    "TGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAAC",
    "TCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTAT",
    "CGAACAGG",
    sep = "", collaspe = ""
  )
  seq8 <- paste("TACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTT",
    "AAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGA",
    "GTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAA",
    "CTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTA",
    "TCAAACAGG",
    sep = "", collaspe = ""
  )
  seq9 <- paste("TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCCAGAC",
    "AAGTCTGAAGTGAAAATCCAGCGCTTAACGTTGGAAGTGCTTTGGAAACTGCCGGGCTAGA",
    "GTGCAGGAGGGGCAGGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAA",
    "CACCAGTGGCGAAGGCGGCCTGCTGGACTGCAACTGACGTTGAGGCTCGAAGGCGTGGGGA",
    "GCAAACAGG",
    sep = "", collaspe = ""
  )
  seq10 <- paste("TACGTAGGTGGCGAGCGTTGTCCGGATTTACTGGGCGTAAAGGGAGCGTAGGCGGACTTT",
    "TAAGTGAGATGTGAAATACTCGGGCTCAACTTGAGTGCTGCATTTCAAACTGGAAGTCTA",
    "GAGTGCAGGAGAGGAGAATGGAATTCCTAGTGTAGCGGTGAAATGCGTAGAGATTAGGAA",
    "GAACACCAGTGGCGAAGGCGATTCTCTGGACTGTAACTGACGCTGAGGCTCGAAAGCGTG",
    "GGGAGCAAACAGG",
    sep = "", collaspe = ""
  )

  first_ten_seqs <- c(
    seq1, seq2, seq3, seq4, seq5,
    seq6, seq7, seq8, seq9, seq10
  )

  expect_equal(first_ten_seqs, dataset$get_sequences()[1:10])
})

test_that("test R6 sequence_data - addSeqs, assign samples", {
  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")

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

  data <- sequence_data$new("mydata")
  data$add_sequences(names, seqs)
  data$assign_sequence_abundance(ids, abundances, samples, treatments)

  # assign otus
  otus <- c("otu1", "otu2", "otu1", "otu2")
  data$assign_otus(otus, seq_ids = names)
  expect_equal(data$get_num_otus(), 2)

  expect_equal(data$get_list()$otu_id, c("otu1", "otu1", "otu2", "otu2"))
  expect_equal(data$get_list()$seq_id, c("seq1", "seq3", "seq2", "seq4"))

  expect_equal(data$get_rabund()$otu_id, c("otu1", "otu2"))
  expect_equal(data$get_rabund()$abundance, c(1200, 120))

  expect_equal(data$get_shared()$otu_id, c("otu1", "otu1", "otu1",
                                           "otu2", "otu2", "otu2"))
  expect_equal(data$get_shared()$sample, c("sample2", "sample3", "sample4",
                                           "sample2", "sample3", "sample4"))
  expect_equal(data$get_shared()$abundance, c(275, 425, 500, 26, 40, 54))

  expect_true(data$is_aligned())
  expect_equal(data$get_ids(), names)
  expect_equal(data$get_sequences(), seqs)
  expect_equal(data$get_ids("sample2"), names)
  expect_equal(data$get_ids("sample3"), c("seq1", "seq2", "seq3"))
  expect_equal(data$get_ids("sample4"), c("seq1", "seq2", "seq4"))
  expect_equal(data$get_sequences("sample2"), seqs)
  expect_equal(data$get_sequences("sample3"), c("ATTGC", "ATTGC", "ATTGC"))
  expect_equal(data$get_sequences("sample4"), c("ATTGC", "ATTGC", "ATTGC"))

  expect_equal(data$get_num_samples(), 3)
  expect_equal(data$get_num_treatments(), 2)
  expect_equal(data$get_treatments(), c("early", "late"))
  expect_equal(data$get_samples(), c("sample2", "sample3", "sample4"))

  sample_summary <- data$get_sample_summary(TRUE)
  expect_equal(sample_summary[[2]]$total, c(766, 554))
  expect_equal(sample_summary[[1]]$total, c(301, 465, 554))

  # total
  expect_equal(data$get_num_sequences(), 1320)
  # unique
  expect_equal(data$get_num_sequences(TRUE), 4)
  # unique and sample
  expect_equal(data$get_num_sequences(TRUE, "sample2"), 4)
  expect_equal(data$get_num_sequences(TRUE, "sample3"), 3)
  # total and sample
  expect_equal(data$get_num_sequences(sample = "sample2"), 301)

  results <- list(sequence_scrap_report = data.frame(),
                  otu_scrap_report = data.frame())
  expect_equal(results, data$get_scrap_report())
})

test_that("test R6 sequence_data - assign_sequence_abundance, remove_sequences", {
  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")

  ids <- c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2",
    "seq3",
    "seq4"
  )
  groups <- c(
    "sample2", "sample3", "sample4",
    "sample3", "sample4",
    "sample3",
    "sample4"
  )
  abundances <- c(
    250, 400, 500,
    40, 50,
    25,
    4
  )
  rabunds <- c(1150, 90, 25, 4)
  treatments <- c(
      "early", "early", "late",
      "late", "late", "late",
      "late"
  )

  seqs_to_remove <- c("seq1", "seq2")
  trash_codes <- c("trashTest", "trashTest2")

  data <- sequence_data$new("mydata")
  data$add_sequences(names, seqs)

  data$assign_sequence_abundance(names, rabunds)

  expect_equal(data$get_num_sequences(), 1269)
  expect_equal(data$get_num_sequences(TRUE), 4)
  expect_equal(data$get_num_samples(), 0)
  expect_equal(data$get_num_treatments(), 0)

  expect_error(data$assign_sequence_abundance(ids, c()))
  missing_ids <- c(
      "seq1", "seq1", "seq1",
      "seq2", "seq2",
      "seq3",
      "seq3"
  )
  expect_error(data$assign_sequence_abundance(missing_ids, abundances))

  data$assign_sequence_abundance(ids, abundances, groups)

  expect_equal(data$get_num_sequences(), 1269)
  expect_equal(data$get_num_sequences(TRUE), 4)
  expect_equal(data$get_num_samples(), 3)
  expect_equal(data$get_num_treatments(), 0)

  data$assign_sequence_abundance(ids, abundances, groups, treatments)
  expect_equal(data$get_num_treatments(), 2)

  data$data$remove_sequences(seqs_to_remove, trash_codes)

  expect_equal(data$get_num_sequences(), 29)
  expect_equal(data$get_num_sequences(TRUE), 2)
  expect_equal(data$get_num_samples(), 2)
  expect_equal(data$get_num_treatments(), 1)
})

test_that("test R6 sequence_data - get_list, get_rabund, get_shared", {
    dataset <- sequence_data$new("my_dataset")
    seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
    otu_ids <- c("otu1", "otu1", "otu1", "otu2", "otu2", "otu3")
    sequence_abundances <- c(10, 100, 1, 500, 25, 80)
    dataset$assign_otus(otu_ids = otu_ids,
                                 abundances = sequence_abundances,
                                 seq_ids = seq_ids)
    # otus would look like:
    # label  otu1             otu2        otu3 ...
    # 0.03   seq1,seq2,seq4   seq3,seq6   seq5 ...
    # 0.03   110              525         80 ...

    list <- dataset$get_list()

    expect_equal(list$otu_id, otu_ids)
    expect_equal(list$seq_id, seq_ids)

    rabund <- dataset$get_rabund()

    abunds <- c(111, 525, 80)
    expect_equal(rabund$otu_id, unique(otu_ids))
    expect_equal(rabund$abundance, abunds)

    dataset <- sequence_data$new("my_dataset")
    otu_ids <- c("otu1", "otu1", "otu1", "otu2", "otu2", "otu3")
    samples <- c("sample1", "sample2", "sample5",
                 "sample1", "sample3", "sample1")
    sample_abundances <- c(10, 100, 1, 500, 25, 80)
    dataset$assign_otus(otu_ids, sample_abundances, samples)

    shared <- dataset$get_shared()
    expect_equal(shared$otu_id, otu_ids)
    expect_equal(shared$abundance, sample_abundances)
    expect_equal(shared$sample, samples)
    expect_equal(dataset$get_num_otus(), 3)

})

test_that("test R6 sequence_data - print", {
  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")

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

  data <- sequence_data$new("mydata")
  data$add_sequences(names, seqs)
  data$assign_sequence_abundance(ids, abundances, samples, treatments)

  expect_snapshot(
    waldo::compare(data$print(), data$print())
  )
})

test_that("test R6 sequence_data - get_align_report, get_contigs_report", {
    names <- c("seq1", "seq2", "seq3", "seq4")
    seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")

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

    data <- sequence_data$new("mydata")
    data$add_sequences(names, seqs)
    data$assign_sequence_abundance(ids, abundances, samples)

    treatments <- c("early", "early", "late")

    data$assign_treatments(unique(samples), treatments)

    expect_equal(data$get_align_report(), data.frame())

    search_scores <- c(77.7, 87.6, 98.5, 75.6)
    sim_scores <- c(97.7, 97.6, 98.5, 95.6)
    longest_inserts <- c(5, 2, 8, 1)

    data$data$add_align_report(unique(ids), search_scores,
                               sim_scores, longest_inserts)


    align_report <- data$get_align_report()
    expect_equal(align_report$id, unique(ids))
    expect_equal(round(align_report$search_score, digits = 2),
                round(search_scores, digits = 2))
    expect_equal(round(align_report$sim_score, digits = 2),
                 round(sim_scores, digits = 2))
    expect_equal(align_report$longest_insert, longest_inserts)

    overlap_lengths <- c(3, 4, 5, 3)
    overlap_starts <- c(2, 1, 1, 1)
    overlap_ends <- c(5, 4, 5, 3)
    mismatches <- c(0, 1, 0, 0)
    expected_errors <- c(5.7, 0.6, 34.5, 3.6)

    data$data$add_contigs_report(unique(ids), overlap_lengths,
                                 overlap_starts, overlap_ends,
                                 mismatches, expected_errors)


    contigs_report <- data$get_contigs_report()
    expect_equal(contigs_report$id, unique(ids))
    expect_equal(contigs_report$length, c(5,5,5,5))
    expect_equal(contigs_report$num_n, c(0,0,0,0))
    expect_equal(contigs_report$overlap_length, overlap_lengths)
    expect_equal(contigs_report$overlap_start, overlap_starts)
    expect_equal(contigs_report$overlap_end, overlap_ends)
    expect_equal(contigs_report$mismatches, mismatches)
    expect_equal(round(contigs_report$ee, digits = 2),
                 round(expected_errors, digits = 2))

    summary <- data$get_sequence_summary()

    expect_equal(length(summary), 3)
    expect_equal(summary$contigs_summary$olengths, c(3, 3, 3, 3, 3, 5, 5, 3))
    expect_equal(summary$align_summary$longest_inserts,
                 c(1,2,5,5,5,8,8,4))

})

