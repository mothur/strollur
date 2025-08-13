# test "sequence_data"

test_that("sequence_data - intialize from read_mothur / print", {
  dataset <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  expect_equal(dataset$get_dataset_name(), "miseq_sop")
  expect_equal(dataset$get_num_sequences(TRUE), 2425)
  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_bins("otu"), 531)

  # add phylotype list
  phylo_list <- read_mothur_list(list = rdataset_example("final.tx.list"))
  dataset$assign_bins(phylo_list$bin_id,
    abundances = NULL,
    samples = NULL, seq_id = phylo_list$seq_id,
    type = "phylotype"
  )

  expect_equal(dataset$get_num_sequences(TRUE), 2425)
  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_bins("phylotype"), 63)

  asv_list <- read_mothur_list(list = rdataset_example("final.asv.list"))
  dataset$assign_bins(asv_list$bin_id,
    abundances = NULL,
    samples = NULL, seq_id = asv_list$seq_id,
    type = "asv"
  )

  expect_equal(dataset$get_num_sequences(TRUE), 2425)
  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_bins("asv"), 2425)

  expect_snapshot(
    waldo::compare(dataset$print(), dataset$print())
  )

  # remove bin from "phylotype" list and confirm that it removes seqs from all
  # from all list types

  expect_equal(dataset$data$get_bin_abundance("Phylo05", "phylotype"), 5337)
  expect_equal(dataset$data$get_bin_abundance("Phylo06", "phylotype"), 715)

  phylo05 <- dataset$data$get_bin("Phylo05", "phylotype")
  expect_equal(length(split_at_char(phylo05)), 54)

  phylo06 <- dataset$data$get_bin("Phylo06", "phylotype")
  expect_equal(length(split_at_char(phylo06)), 47)

  dataset$data$remove_bins(
    c("Phylo05", "Phylo06"),
    c("test", "test"),
    "phylotype"
  )

  expect_equal(dataset$get_num_bins("phylotype"), 61)
  expect_equal(dataset$get_num_bins("otu"), 512)
  expect_equal(dataset$get_num_bins("asv"), 2324)
  expect_equal(dataset$get_num_sequences(), 107911)
  expect_equal(dataset$get_num_sequences(TRUE), 2324)
  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_sequences(sample = "F3D0"), 5977)
  expect_equal(dataset$get_num_sequences(sample = "F3D1"), 4467)
  # note that the number of seqs removed will be less that 297+266 because
  # some seqs are assigned to both samples and some seqs will be present in
  # other samples
  expect_equal(dataset$get_num_sequences(TRUE, "F3D0"), 297)
  expect_equal(dataset$get_num_sequences(TRUE, "F3D1"), 266)

  # remove samples
  dataset$remove_samples(c("F3D0", "F3D1"))
  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 17)
  expect_equal(dataset$get_num_bins("phylotype"), 57)
  expect_equal(dataset$get_num_bins("otu"), 482)
  expect_equal(dataset$get_num_bins("asv"), 2124)
  expect_equal(dataset$get_num_sequences(), 107700)
  expect_equal(dataset$get_num_sequences(TRUE), 2124)
  expect_equal(dataset$get_num_sequences(TRUE, "F3D0"), 0)
  expect_equal(dataset$get_num_sequences(TRUE, "F3D1"), 0)
  expect_equal(dataset$get_num_sequences(sample = "F3D0"), 0)
  expect_equal(dataset$get_num_sequences(sample = "F3D1"), 0)

  # remove things just classified to bacteria, and things classified to
  # Bacteria;"Bacteroidetes"; with confidence less than 95
  dataset$remove_lineages(c(
    "Bacteria;Bacteria_unclassified;",
    "Bacteria(100);\"Bacteroidetes\"(95);"
  ))

  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 17)
  expect_equal(dataset$get_num_bins("phylotype"), 57)
  expect_equal(dataset$get_num_bins("otu"), 475)
  expect_equal(dataset$get_num_bins("asv"), 2086)
  expect_equal(dataset$get_num_sequences(), 107661)
  expect_equal(dataset$get_num_sequences(TRUE), 2086)
})

test_that("sequence_data - intialize from sequence_data object", {
  temp <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  # add phylotype list
  phylo_list <- read_mothur_list(list = rdataset_example("final.tx.list"))
  temp$assign_bins(phylo_list$bin_id,
    abundances = NULL,
    samples = NULL, seq_id = phylo_list$seq_id,
    type = "phylotype"
  )

  asv_list <- read_mothur_list(list = rdataset_example("final.asv.list"))
  temp$assign_bins(asv_list$bin_id,
    abundances = NULL,
    samples = NULL, seq_id = asv_list$seq_id,
    type = "asv"
  )

  dataset <- sequence_data$new(
    name = "clone_of_miseq", dataset = temp,
    processors = 4
  )

  expect_equal(dataset$get_dataset_name(), "clone_of_miseq")
  expect_equal(dataset$get_num_sequences(TRUE), 2425)
  expect_equal(dataset$get_num_sequences(), 113963)
  expect_equal(dataset$get_num_treatments(), 2)
  expect_equal(dataset$get_num_samples(), 19)
  expect_equal(dataset$get_num_bins("otu"), 531)
  expect_equal(dataset$get_num_bins("phylotype"), 63)
  expect_equal(dataset$get_num_bins("asv"), 2425)
  expect_equal(dataset$data$processors, 4)
})

test_that("sequence_data - addSeqs, assign samples", {
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

  # include reference
  url <- "https://mothur.org/wiki/silva_reference_files/"
  data$add_sequences(names, seqs,
    reference_name = "silva.bacteria.fasta",
    reference_note = "alignment by mothur2 v1.0",
    reference_version = "1.38.1", reference_url = url
  )

  references <- data$get_references()

  expect_equal(nrow(references), 1)
  expect_equal(references[[1, "reference_name"]], "silva.bacteria.fasta")
  expect_equal(references[[1, "version"]], "1.38.1")
  expect_equal(references[[1, "usage"]], NA)
  expect_equal(references[[1, "note"]], "alignment by mothur2 v1.0")
  expect_equal(references[[1, "url"]], url)


  data$assign_sequence_abundance(ids, abundances, samples, treatments)

  # assign bins
  bins <- c("bin1", "bin2", "bin1", "bin2")
  data$assign_bins(bins, seq_ids = names)

  expect_equal(data$get_num_bins(), 2)
  expect_equal(data$get_list()$otu_id, c("bin1", "bin1", "bin2", "bin2"))
  expect_equal(data$get_list()$seq_id, c("seq1", "seq3", "seq2", "seq4"))

  # get asv generated by dataset
  expect_equal(data$get_list("asv")$asv_id, c("ASV1", "ASV2", "ASV3", "ASV4"))
  expect_equal(data$get_list("asv")$seq_id, c("seq1", "seq2", "seq3", "seq4"))

  expect_equal(data$get_rabund()$otu_id, c("bin1", "bin2"))
  expect_equal(data$get_rabund()$abundance, c(1200, 120))

  expect_equal(data$get_shared()$otu_id, c(
    "bin1", "bin1", "bin1",
    "bin2", "bin2", "bin2"
  ))
  expect_equal(data$get_shared()$sample, c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  ))
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

  results <- list(
    sequence_scrap_report = data.frame(),
    otu_scrap_report = data.frame(),
    asv_scrap_report = data.frame()
  )

  expect_equal(results, data$get_scrap_report())
})

test_that("sequence_data - assign_sequence_abundance, remove_sequences", {
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

test_that("sequence_data - get_list get_rabund, get_shared", {
  dataset <- sequence_data$new("my_dataset")
  seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  sequence_abundances <- c(10, 100, 1, 500, 25, 80)
  dataset$assign_bins(
    bin_ids = bin_ids,
    abundances = sequence_abundances,
    seq_ids = seq_ids
  )
  # bins would look like:
  # label  bin1             bin2        bin3 ...
  # 0.03   seq1,seq2,seq4   seq3,seq6   seq5 ...
  # 0.03   110              525         80 ...

  list <- dataset$get_list()

  expect_equal(list$otu_id, bin_ids)
  expect_equal(list$seq_id, seq_ids)

  rabund <- dataset$get_rabund()

  abunds <- c(111, 525, 80)
  expect_equal(rabund$otu_id, unique(bin_ids))
  expect_equal(rabund$abundance, abunds)

  dataset <- sequence_data$new("my_dataset")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  samples <- c(
    "sample1", "sample2", "sample5",
    "sample1", "sample3", "sample1"
  )
  sample_abundances <- c(10, 100, 1, 500, 25, 80)
  dataset$assign_bins(bin_ids, sample_abundances, samples)

  shared <- dataset$get_shared()
  expect_equal(shared$otu_id, bin_ids)
  expect_equal(shared$abundance, sample_abundances)
  expect_equal(shared$sample, samples)
  expect_equal(dataset$get_num_bins(), 3)

  dataset <- sequence_data$new("my_dataset")
  bin_ids <- c(
    "bin1", "bin1", "bin1", "bin1",
    "bin2", "bin2", "bin2", "bin2", "bin2", "bin2", "bin2", "bin2",
    "bin3", "bin3", "bin3", "bin3"
  )
  seq_ids <- c(
    "seq1", "seq2", "seq3", "seq3",
    "seq4", "seq4", "seq5", "seq6", "seq7", "seq8", "seq9", "seq9",
    "seq10", "seq10", "seq10", "seq10"
  )
  samples <- c(
    "sample1", "sample2", "sample4", "sample5",
    "sample1", "sample2", "sample1", "sample1",
    "sample2", "sample4", "sample4", "sample5",
    "sample1", "sample3", "sample5", "sample6"
  )
  sample_abundances <- c(
    10, 10, 5, 5,
    5, 5, 10, 10, 10, 10, 5, 5,
    1, 2, 3, 4
  )

  dataset$assign_bins(bin_ids, sample_abundances, samples, seq_ids)

  expect_equal(dataset$get_num_bins(), 3)
  expect_equal(dataset$get_num_samples(), 6)
  expect_equal(
    dataset$get_sample_summary(TRUE)[[1]]$total,
    c(36, 25, 2, 20, 13, 4)
  )
})


test_that("sequence_data - get_align_report, get_contigs_report", {
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

  data$data$add_align_report(
    unique(ids), search_scores,
    sim_scores, longest_inserts
  )


  align_report <- data$get_align_report()
  expect_equal(align_report$id, unique(ids))
  expect_equal(
    round(align_report$search_score, digits = 2),
    round(search_scores, digits = 2)
  )
  expect_equal(
    round(align_report$sim_score, digits = 2),
    round(sim_scores, digits = 2)
  )
  expect_equal(align_report$longest_insert, longest_inserts)

  overlap_lengths <- c(3, 4, 5, 3)
  overlap_starts <- c(2, 1, 1, 1)
  overlap_ends <- c(5, 4, 5, 3)
  mismatches <- c(0, 1, 0, 0)
  expected_errors <- c(5.7, 0.6, 34.5, 3.6)

  data$data$add_contigs_report(
    unique(ids), overlap_lengths,
    overlap_starts, overlap_ends,
    mismatches, expected_errors
  )


  contigs_report <- data$get_contigs_report()
  expect_equal(contigs_report$id, unique(ids))
  expect_equal(contigs_report$length, c(5, 5, 5, 5))
  expect_equal(contigs_report$num_n, c(0, 0, 0, 0))
  expect_equal(contigs_report$overlap_length, overlap_lengths)
  expect_equal(contigs_report$overlap_start, overlap_starts)
  expect_equal(contigs_report$overlap_end, overlap_ends)
  expect_equal(contigs_report$mismatches, mismatches)
  expect_equal(
    round(contigs_report$ee, digits = 2),
    round(expected_errors, digits = 2)
  )

  summary <- data$get_sequence_summary()

  expect_equal(length(summary), 3)
  expect_equal(summary$contigs_summary$olengths, c(3, 3, 3, 3, 3, 5, 5, 3))
  expect_equal(
    summary$align_summary$longest_inserts,
    c(1, 2, 5, 5, 5, 8, 8, 4)
  )
})

# assign_sequence_taxonomy, get_sequence_taxonomy_report
test_that("sequence_data - ", {
  # same length with confidences
  names <- c("seq1", "seq2", "seq3", "seq4")

  tax1 <- "Bacteria(100);Bacteroidetes(95);Bacteroidia(90);"
  tax2 <- "Bacteria(100);Proteobacteria(89);Betaproteobacteria(85);"
  tax3 <- "Bacteria(100);Firmicutes(99);Bacilli(90);"
  tax4 <- "Bacteria(100);Proteobacteria(87);Gammaproteobacteria(82);"
  taxonomies <- c(tax1, tax2, tax3, tax4)

  dataset <- sequence_data$new("my_dataset")

  url <- paste("https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset",
    "9_032012.pds.zip",
    collapse = ""
  )

  # assign taxonomy with reference
  dataset$assign_sequence_taxonomy(
    names = names, taxonomies = taxonomies,
    reference_name = "trainset9_032012.pds.zip",
    reference_note = "classification by mothur2 v1.0 using default options",
    reference_version = "9_032012", reference_url = url
  )

  references <- dataset$get_references()

  note <- "classification by mothur2 v1.0 using default options"
  expect_equal(nrow(references), 1)
  expect_equal(references[[1, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[1, 2]], "9_032012")
  expect_equal(references[[1, 3]], "sequence_classification")
  expect_equal(references[[1, 4]], note)
  expect_equal(references[[1, "url"]], url)

  report <- dataset$get_sequence_taxonomy_report()

  ids <- c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2",
    "seq3", "seq3", "seq3",
    "seq4", "seq4", "seq4"
  )
  expect_equal(report$id, ids)
  expect_equal(report$taxon[1:6], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia", "Bacteria",
    "Proteobacteria", "Betaproteobacteria"
  ))
  expect_equal(report$confidence[7:12], c(100, 99, 90, 100, 87, 82))

  # same length without confidences
  tax1 <- "Bacteria;Bacteroidetes;Bacteroidia;"
  tax2 <- "Bacteria;Proteobacteria;Betaproteobacteria;"
  tax3 <- "Bacteria;Firmicutes;Bacilli;"
  tax4 <- "Bacteria;Proteobacteria;Gammaproteobacteria;"
  taxonomies <- c(tax1, tax2, tax3, tax4)

  dataset <- sequence_data$new("my_dataset")
  dataset$assign_sequence_taxonomy(names, taxonomies)
  report <- dataset$get_sequence_taxonomy_report()

  expect_equal(report$id, ids)
  expect_equal(report$taxon[7:12], c(
    "Bacteria", "Firmicutes",
    "Bacilli", "Bacteria",
    "Proteobacteria", "Gammaproteobacteria"
  ))

  # different lengths without confidences
  tax1 <- "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;"
  tax2 <- "Bacteria;Proteobacteria;Betaproteobacteria;"
  tax3 <- "Bacteria;Firmicutes;Bacilli;Lactobacillales;Streptococcaceae;"
  tax4 <- "Bacteria;Proteobacteria;"
  taxonomies <- c(tax1, tax2, tax3, tax4)

  dataset <- sequence_data$new("my_dataset")
  dataset$assign_sequence_taxonomy(names, taxonomies)
  report <- dataset$get_sequence_taxonomy_report()

  ids <- c(
    "seq1", "seq1", "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2", "seq2", "seq2",
    "seq3", "seq3", "seq3", "seq3", "seq3",
    "seq4", "seq4", "seq4", "seq4", "seq4"
  )
  expect_equal(report$id, ids)
  expect_equal(report$taxon[1:10], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia", "Bacteroidales",
    "Bacteroidales_unclassified",
    "Bacteria", "Proteobacteria",
    "Betaproteobacteria",
    "Betaproteobacteria_unclassified",
    "Betaproteobacteria_unclassified"
  ))
  expect_equal(report$taxon[11:20], c(
    "Bacteria", "Firmicutes",
    "Bacilli", "Lactobacillales",
    "Streptococcaceae",
    "Bacteria", "Proteobacteria",
    "Proteobacteria_unclassified",
    "Proteobacteria_unclassified",
    "Proteobacteria_unclassified"
  ))

  # different lengths with confidences and extra '()'
  tax1 <- "Bacteria(100);Bacteroidetes(95);Bacteroidia(sub_category)(87);"
  tax1 <- paste0(tax1, "Bacteroidales(85);")
  tax2 <- "Bacteria(100);Proteobacteria(95);Betaproteobacteria(87);"
  tax3 <- "Bacteria(100);Firmicutes(95);Bacilli(87);Lactobacillales(87);"
  tax3 <- paste0(tax3, "Streptococcaceae(85)")
  tax4 <- "Bacteria(100);Proteobacteria(95);"
  taxonomies <- c(tax1, tax2, tax3, tax4)

  dataset <- sequence_data$new("my_dataset")
  dataset$assign_sequence_taxonomy(names, taxonomies)
  report <- dataset$get_sequence_taxonomy_report()

  ids <- c(
    "seq1", "seq1", "seq1", "seq1", "seq1",
    "seq2", "seq2", "seq2", "seq2", "seq2",
    "seq3", "seq3", "seq3", "seq3", "seq3",
    "seq4", "seq4", "seq4", "seq4", "seq4"
  )
  expect_equal(report$id, ids)
  expect_equal(report$taxon[1:10], c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia(sub_category)",
    "Bacteroidales",
    "Bacteroidales_unclassified",
    "Bacteria", "Proteobacteria",
    "Betaproteobacteria",
    "Betaproteobacteria_unclassified",
    "Betaproteobacteria_unclassified"
  ))
  expect_equal(report$taxon[11:20], c(
    "Bacteria", "Firmicutes",
    "Bacilli", "Lactobacillales",
    "Streptococcaceae",
    "Bacteria", "Proteobacteria",
    "Proteobacteria_unclassified",
    "Proteobacteria_unclassified",
    "Proteobacteria_unclassified"
  ))
  expect_equal(report$confidence[1:10], c(
    100, 95, 87, 85, 85, 100, 95, 87,
    87, 87
  ))
  expect_equal(report$confidence[11:20], c(
    100, 95, 87, 87, 85, 100, 95, 95,
    95, 95
  ))

  # add bin assignments
  bins <- c("bin1", "bin1", "bin1", "bin2")
  dataset$assign_bins(bins, seq_ids = names)

  report <- dataset$get_bin_taxonomy_report()

  ids <- c(
    "bin1", "bin1", "bin1", "bin1",
    "bin2", "bin2", "bin2", "bin2"
  )

  expect_equal(ids, report$id)
  expect_equal(c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia(sub_category)",
    "Bacteroidales",
    "Bacteria", "Proteobacteria",
    "Proteobacteria_unclassified",
    "Proteobacteria_unclassified"
  ), report$taxon)
  expect_equal(report$confidence, c(100, 34, 34, 34, 100, 100, 100, 100))

  dataset$clear()
  bin_ids <- c("bin1", "bin2", "bin3", "bin4")
  taxonomies <- c(
    "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;"
  )

  expect_error(dataset$assign_bin_taxonomy(bin_ids, taxonomies))
  expect_error(dataset$assign_bins(bin_ids))
  expect_equal(dataset$get_contigs_report(), data.frame())
  expect_equal(dataset$get_sample_summary(), list())
  expect_equal(dataset$get_sequence_summary(), list())
  expect_false(dataset$has_sample("noSample"))

  abunds <- c(200, 40, 100, 5)
  dataset$assign_bins(bin_ids, abunds)
  dataset$assign_bin_taxonomy(bin_ids, taxonomies,
    reference_name = "trainset9_032012.pds.zip",
    reference_note = "classification by mothur2 v1.0",
    reference_version = "9_032012",
    reference_url = url
  )

  references <- dataset$get_references()

  note <- "classification by mothur2 v1.0"
  expect_equal(nrow(references), 1)
  expect_equal(references[[1, "reference_name"]], "trainset9_032012.pds.zip")
  expect_equal(references[[1, "version"]], "9_032012")
  expect_equal(references[[1, "usage"]], "bin_classification")
  expect_equal(references[[1, "note"]], note)
  expect_equal(references[[1, "url"]], url)

  report <- dataset$get_bin_taxonomy_report()

  ids <- c(
    "bin1", "bin1", "bin1", "bin1",
    "bin2", "bin2", "bin2", "bin2",
    "bin3", "bin3", "bin3", "bin3",
    "bin4", "bin4", "bin4", "bin4"
  )

  expect_equal(report$id, ids)
  expect_equal(report$taxon, c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia", "Bacteroidales",
    "Bacteria", "Proteobacteria",
    "Betaproteobacteria", "Neisseriales",
    "Bacteria", "Firmicutes",
    "Bacilli", "Lactobacillales",
    "Bacteria", "Proteobacteria",
    "Gammaproteobacteria", "Pasteurellales"
  ))

  # different lengths and with confidence scores
  taxonomies <- c(
    "Bacteria(100);Bacteroidetes(95);Bacteroidia(94);",
    "Bacteria(100);Proteobacteria(90);Betaproteobacteria(87);Neisseriales(80);",
    "Bacteria(100);Firmicutes(100);Bacilli(100);Lactobacillales(89);",
    "Bacteria(100);Proteobacteria(65);Gammaproteobacteria(60);"
  )

  dataset$assign_bin_taxonomy(bin_ids, taxonomies)

  report <- dataset$get_bin_taxonomy_report()

  expect_equal(report$id, ids)
  expect_equal(report$taxon, c(
    "Bacteria", "Bacteroidetes",
    "Bacteroidia", "Bacteroidia_unclassified",
    "Bacteria", "Proteobacteria",
    "Betaproteobacteria", "Neisseriales",
    "Bacteria", "Firmicutes",
    "Bacilli", "Lactobacillales",
    "Bacteria", "Proteobacteria",
    "Gammaproteobacteria", "Gammaproteobacteria_unclassified"
  ))
  expect_equal(report$confidence, c(
    100, 95, 94, 94, 100, 90, 87, 80,
    100, 100, 100, 89, 100, 65, 60, 60
  ))

  expect_equal(dataset$get_sequence_taxonomy_report(), data.frame())
  expect_equal(dataset$get_shared(), data.frame())
})

test_that("sequence_data - add_metadata, get_metadata", {
  dataset <- sequence_data$new("my_dataset")

  expect_equal(dataset$get_metadata(), data.frame())
  expect_error(dataset$add_metadata(c("bad_type")))

  metadata <- readr::read_tsv(rdataset_example("sample-metadata.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )

  dataset$add_metadata(metadata)
  metadata <- dataset$get_metadata()

  expect_equal(names(metadata), c(
    "sample-id", "barcode-sequence",
    "body-site", "year", "month", "day",
    "subject", "reported-antibiotic-usage",
    "days-since-experiment-start"
  ))
  expect_equal(nrow(metadata), 34)

  # random spot checks
  expect_equal(metadata[[5, 8]], "No")
  expect_equal(metadata[[8, 3]], "left palm")
  expect_equal(metadata[[2, 3]], "gut")
  expect_equal(metadata[[3, 7]], "subject-1")
})

test_that("sequence_data - add_references, get_references", {
  dataset <- sequence_data$new("my_dataset")

  expect_equal(dataset$get_references(), data.frame())
  expect_error(dataset$add_references(reference = c("bad_type")))
  expect_error(dataset$add_references(data.frame()))
  expect_error(dataset$add_references())

  references <- dataset$get_references()
  expect_equal(nrow(references), 0)

  reference <- readr::read_csv(rdataset_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )

  mothur_url <- "https://github.com/mothur/mothur/releases/tag/v1.48.2"

  # add data.frame and single reference at the same time
  dataset$add_references(
    reference = reference,
    reference_name = "mothur software package",
    usage = "analysis of dataset",
    version = "1.48.2", url = mothur_url
  )

  references <- dataset$get_references()

  # random spot checks
  expect_equal(nrow(references), 3)
  expect_equal(references[[1, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[2, 1]], "silva.v4.fasta")
  expect_equal(references[[3, 1]], "mothur software package")
  expect_equal(references[[1, 4]], NA_character_)
  expect_equal(references[[2, 4]], "1.38.1")
  expect_equal(references[[3, 4]], "1.48.2")


  dataset <- sequence_data$new("my_dataset")

  # add single reference then dataframe
  dataset$add_references(
    reference_name = "mothur software package",
    usage = "analysis of dataset",
    version = "1.48.2", url = mothur_url,
    note = "This is my mothur note"
  )

  references <- dataset$get_references()
  expect_equal(nrow(references), 1)

  dataset$add_references(reference = reference)

  references <- dataset$get_references()
  expect_equal(nrow(references), 3)

  expect_equal(references[[2, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[3, 1]], "silva.v4.fasta")
  expect_equal(references[[1, 1]], "mothur software package")
  expect_equal(references[[2, 2]], NA_character_)
  expect_equal(references[[3, 2]], "1.38.1")
  expect_equal(references[[1, 4]], "This is my mothur note")
})
