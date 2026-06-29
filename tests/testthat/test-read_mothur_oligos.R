test_that("test read_mothur_oligos - errors", {
  # mixed paired and unpaired# paired oligos file
  oligo_data <- c(
    "primer\tCCTACGGGAGGCAGCAG\tATTACCGCGGCTGCTGG\tV3",
    "primer\tATTAGAWACCCBDGTAGTCC\tCCCGTCAATTCMTTTRAGT\tV5",
    "primer\tACTYAAAKGAATTGACGGG\tACRACACGAGCTGACGAC\tV6",
    "barcode\tTAATCG\tACTTGA\tF61",
    "barcode\tATCACG\tTACAGC\tF62",
    "barcode\tGAGATA\tCTAGCT\tF63",
    "barcode\tCGCGGT\tF91",
    "forward\tCCTACGGGAGGCAGCAG"
  )

  file_name <- file.path(tempdir(), "oligos.txt")

  # non existant file
  expect_error(read_mothur_oligos(file_name))

  create_dummy_file(file_name, oligo_data)

  expect_error(
    read_mothur_oligos(file_name),
    regexp = "cannot mix paired and unpaired oligo data"
  )

  # invalid oligos type
  oligo_data <- c(
    "not_a_valid_type\tCCTACGGGAGGCAGCAG\tATTACCGCGGCTGCTGG\tV3"
  )
  create_dummy_file(file_name, oligo_data)

  expect_error(
    read_mothur_oligos(file_name),
    regexp = paste(
      "`not_a_valid_type` is not a valid oligo type. Options",
      "include: linker, spacer, forward, reverse, barcode and",
      "primer."
    )
  )

  # barcodes without barcode names
  oligo_data <- c(
    "barcode\tCCTACGGGAGG"
  )
  create_dummy_file(file_name, oligo_data)

  expect_error(
    read_mothur_oligos(file_name),
    regexp = "barcode names are required"
  )

  # unpaired primers
  oligo_data <- c(
    "primer\tCCTACGGGAGG"
  )
  create_dummy_file(file_name, oligo_data)

  expect_error(
    read_mothur_oligos(file_name),
    regexp = paste(
      "primers must be paired. You can use",
      "'NONE' as a placeholder for paired reads",
      "where only the forward or reverse primer",
      "is present. Use the `forward` or",
      "`reverse` type for unpaired primers."
    )
  )
  remove_file(file_name)
})

test_that("test read_mothur_oligos - paired", {
  oligos <- read_mothur_oligos(strollur_example("paired_read.oligos"))

  expect_equal(oligos$primer_forward, c(
    "CCTACGGGAGGCAGCAG",
    "ATTAGAWACCCBDGTAGTCC",
    "ACTYAAAKGAATTGACGGG", ""
  ))
  expect_equal(oligos$primer_reverse, c(
    "ATTACCGCGGCTGCTGG",
    "CCCGTCAATTCMTTTRAGT",
    "ACRACACGAGCTGACGAC", ""
  ))
  expect_equal(oligos$barcode_forward, c(
    "TAATCG", "ATCACG",
    "GAGATA", "CGCGGT"
  ))
  expect_equal(oligos$barcode_reverse, c(
    "ACTTGA", "TACAGC",
    "CTAGCT", "GAGTGG"
  ))
  expect_equal(oligos$barcode_name, c(
    "F61", "F62",
    "F63", "F91"
  ))
  expect_equal(oligos$primer_name, c("V3", "V5", "", ""))
})

test_that("test read_mothur_oligos - single", {
  # mixed paired and unpaired# paired oligos file
  oligo_data <- c(
    "barcode\tTAATCG\tF61",
    "barcode\tATCACG\tF62",
    "barcode\tCGCGGT\tF91",
    "forward\tCCTACGGGAGGCA",
    "forward\tCCTACGGGAGGCAGCAG\tmy_favorite",
    "reverse\tCCTACGG",
    "linker\tTC",
    "spacer\tTAGC"
  )

  file_name <- file.path(tempdir(), "oligos.txt")
  create_dummy_file(file_name, oligo_data)

  oligos <- read_mothur_oligos(file_name)

  expect_equal(oligos$linker, c("TC", "", ""))
  expect_equal(oligos$spacer, c("TAGC", "", ""))
  expect_equal(oligos$reverse, c("CCTACGG", "", ""))
  expect_equal(oligos$forward, c("CCTACGGGAGGCA", "CCTACGGGAGGCAGCAG", ""))
  expect_equal(oligos$forward_name, c("", "my_favorite", ""))
  expect_equal(oligos$barcode, c("TAATCG", "ATCACG", "CGCGGT"))
  expect_equal(oligos$barcode_name, c("F61", "F62", "F91"))

  remove_file(file_name)
})
