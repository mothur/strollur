test_that("write_fasta - errors", {
  expect_error(write_fasta("Bad_type"))

  # no file name with nameless dataset
  data <- dataset$new()
  expect_error(write_fasta(data))
})

test_that("write_fasta", {
  miseq <- miseq_sop_example()

  output <- write_fasta(miseq, get_full_name("miseq"))

  df <- read_fasta(output)

  remove_file(output)

  expect_equal(df[[1]], get_sequence_names(miseq))
  seq1 <- paste0(
    "TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGTGAT--TT-A-T-T--GG-GT--TT-A-A",
    "A-GG-GT-GC-G-TA-GGC-G-G-A-CA-G-T-T-AA-G-T-C-A-G-C-G-G--TA-A-AA-TT-G-A",
    "-GA-GG--CT-C-AA-C-C-T-C-T-T-C--CC-G-C-CGTT-GAAAC-TG-A-TTGTC-TTGA-GT-G",
    "G-GC-GA-G-A---AG-T-A-TGTGGAATGCGTGGTGT-AGCGGT-GAAATGCATAG-AT-A-TC-AC-",
    "GC-AG-AACTCCGAT-TGCGAAGGCA------GCATA-CCG-G-CG-CC-C-A-ACTGACG-CTGA-AG",
    "CA-CGAAA-GCG-TGGGT-ATC-GAACAGG"
  )

  seq2 <- paste0(
    "TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CTG-AT--TT-A-T-T--GG-GT--TT-A-A",
    "A-GG-GT-GC-G-TA-GGC-G-G-T-TC-G-A-T-AA-G-T-T-A-G-A-G-G--TG-A-AA-TC-C-C",
    "-GG-GG--CT-C-AA-C-T-C-C-G-G-C-ACT-G-C-CTCT-GATAC-TG-T-CGGGC-TAGA-GT-T",
    "T-AG-TT-G-C---GG-T-A-GGCGGAATGTATGGTGT-AGCGGT-GAAATGCATAG-AG-A-TC-AT-",
    "AC-AG-AACACCGAT-TGCGAAGGCA------GCTTA-CCA-A-AC-TA-C-G-ACTGACG-TTGA-GG",
    "CA-CGAAA-GCG-TGGGG-AGC-AAACAGG"
  )

  seq3 <- paste0(
    "TAC--GT-AT-GGA--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-A",
    "A-GG-GA-GC-G-TA-GGC-G-G-C-AT-A-C-C-AA-G-C-C-T-G-A-T-G--TG-A-AA-AC-C-C",
    "-GG-GG--CC-C-AA-C-C-C-C-G-G-G-AGT-G-C-ATTG-GGAAC-TG-G-CAAGC-TAGA-GT-G",
    "T-CG-GA-G-A---GG-C-A-GGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-",
    "GA-GG-AACACCAGT-GGCGAAGGCG------GCCTG-CTG-G-AC-GA-T-G-ACTGACG-CTGA-GG",
    "CT-CGAAA-GCG-TGGGG-AGC-AAACAGG"
  )

  expect_equal(df[[2]][1], seq1)
  expect_equal(df[[2]][100], seq2)
  expect_equal(df[[2]][50], seq3)

  data <- read_mothur(
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "data"
  )

  expect_equal(write_fasta(data), "no_sequence_data")
})
