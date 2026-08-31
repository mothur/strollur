test_that("read_biom - errors", {
  expect_error(read_biom("Bad_file"))
})

test_that("read_biom", {
  data <- strollur::read_biom(biom = strollur_example("miseq_sop.otu.biom"))

  bin_types <- data$get_bin_types()

  expect_equal(data$get_bin_types(), "asv")

  expect_equal(names(data, type = "dataset"), "miseq_sop")

  expect_equal(count(data, type = "bin", bin_type = "asv"), 531)
  expect_equal(names(data, type = "bin", bin_type = "asv")[1], "Otu001")
  expect_equal(count(data, type = "sample"), 19)
  expect_equal(names(data, type = "sample"), c(
    "F3D0", "F3D1", "F3D141",
    "F3D142", "F3D143", "F3D144",
    "F3D145", "F3D146", "F3D147",
    "F3D148", "F3D149", "F3D150",
    "F3D2", "F3D3", "F3D5", "F3D6",
    "F3D7", "F3D8", "F3D9"
  ))
  table <- report(data, type = "bin_taxonomy", bin_type = "asv")
  expect_equal(length(unique(table$level)), 6)
  expect_equal(table[56, 3], "Firmicutes")

  otu1_seq <- paste0(
    "TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T",
    "--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-C-TG-T-G-C-AA-G",
    "-T-C-T-G-A-A-G--TG-A-AA-TG-C-C-GG-GG--CT-C-AA-C-C-C-C-",
    "G-G-A-ACT-G-C-TTTG-GAAAC-TG-T-ACAGC-TAGA-GT-GC-AG-GA-G",
    "-G---GG-T-G-AGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-",
    "TT-AG-GA-GG-AACACCGGT-GGCGAAGGCG------GCTCA-CTG-G-AC-T",
    "G-T-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG"
  )
  expect_equal(xdev_get_sequences(data)[1], otu1_seq)
})
