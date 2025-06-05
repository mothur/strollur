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

  expect_equal(data.frame(), data$get_scrap_report())
})

test_that("test R6 sequence_data - print", {
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
  treatments <- c(
    "early", "early", "late",
    "late", "late", "late",
    "late"
  )

  seqs_to_remove <- c("seq1", "seq2")
  trash_codes <- c("trashTest", "trashTest2")

  data <- sequence_data$new("mydata")
  data$add_sequences(names, seqs)
  data$assign_sequence_abundance(ids, abundances, groups, treatments)

  expect_equal(data$get_num_sequences(), 1269)
  expect_equal(data$get_num_sequences(TRUE), 4)
  expect_equal(data$get_num_samples(), 3)
  expect_equal(data$get_num_treatments(), 2)

  data$data$remove_sequences(seqs_to_remove, trash_codes)

  expect_equal(data$get_num_sequences(), 29)
  expect_equal(data$get_num_sequences(TRUE), 2)
  expect_equal(data$get_num_samples(), 2)
  expect_equal(data$get_num_treatments(), 1)
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
