# test read_fasta

test_that("test read_fasta - errors", {
  expect_error(read_fasta("non_existant_filename"))
})

test_that("test read_fasta - comments", {
  name1 <- ">seq1 my very cool comment"
  seq1 <- paste("TACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCGCGCAGGTGGTTAATT",
    "AAGTCTGATGTGAAAGCCCACGGCTTAACCGTGGAGGGTCATTGGAAACTGGTTGACTTGA",
    "GTGCAGAAGAGGGAAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGAGATATGGAGGAA",
    "CACCAGTGGCGAAGGCGGCTTCCTGGTCTGCAACTGACACTGAGGCGCGAAAGCGTGGGGA",
    "GCAAACAGG",
    sep = "", collaspe = ""
  )
  name2 <- ">seq2 my next comment"
  seq2 <- paste("TACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCGCGCAGGTGGTTAATT",
    "AAGTCTGATGTGAAAGCCCACGGCTTAACCGTGGAGGGTCATTGGAAACTGGTTGACTTGA",
    "GTGCAGAAGAGGGAAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGAGATATGGAGGAA",
    "CACCAGTGGCGAAGGCGGCTTCCTGGTCTGCAACTGACACTGAGGCGCGAAAGCGTGGGGA",
    "GCAAACAGG",
    sep = "", collaspe = ""
  )
  name3 <- ">seq3"
  seq3 <- paste("TACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGCGCGCAGGTGGTTAATT",
    "AAGTCTGATGTGAAAGCCCACGGCTTAACCGTGGAGGGTCATTGGAAACTGGTTGACTTGA",
    "GTGCAGAAGAGGGAAGTGGAATTCCATGTGTAGCGGTGAAATGCGTAGAGATATGGAGGAA",
    "CACCAGTGGCGAAGGCGGCTTCCTGGTCTGCAACTGACACTGAGGCGCGAAAGCGTGGGGA",
    "GCAAACAGG",
    sep = "", collaspe = ""
  )

  create_dummy_file(
    "file_with_comments.fasta",
    c(name1, seq1, name2, seq2, name3, seq3)
  )

  # read fasta file with comments
  results <- read_fasta("file_with_comments.fasta")

  remove_file("file_with_comments.fasta")

  expect_equal(names(results), c(
    "sequence_names",
    "sequences", "comments"
  ))
  expect_equal(results$sequence_names, c("seq1", "seq2", "seq3"))
  expect_equal(rep(seq1, 3), results$sequences)
  expect_equal(
    c("my very cool comment", "my next comment", ""),
    results$comments
  )
})
