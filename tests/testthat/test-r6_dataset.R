# test "dataset"

test_that("dataset - intialize from read_mothur / print", {
  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    phylo_list = rdataset_example("final.tx.list"),
    asv_list = rdataset_example("final.asv.list"),
    dataset_name = "miseq_sop"
  )

  expect_equal(dataset_t$get_bin_types(), c("otu", "asv", "phylotype"))
  expect_equal(dataset_t$get_dataset_name(), "miseq_sop")
  expect_equal(dataset_t$get_num_sequences(TRUE), 2425)
  expect_equal(dataset_t$get_num_sequences(), 113963)
  expect_equal(dataset_t$get_num_treatments(), 2)
  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_bins("otu"), 531)
  expect_equal(dataset_t$get_num_bins("phylotype"), 63)
  expect_equal(dataset_t$get_num_bins("asv"), 2425)

  seqs_summary <- dataset_t$get_summary()[["sequence_summary"]]

  expect_equal(seqs_summary$starts[1], 1)
  expect_equal(seqs_summary$ends[1], 375)
  expect_equal(seqs_summary$ambigs[2], 0)
  expect_equal(seqs_summary$numns[4], 0)
  expect_equal(seqs_summary$nbases[1], 249)
  expect_equal(seqs_summary$nbases[7], 256)
  expect_equal(seqs_summary$polymers[1], 3)
  expect_equal(seqs_summary$polymers[7], 6)

  expect_snapshot(
    waldo::compare(dataset_t$print(), dataset_t$print())
  )

  # remove bin from "phylotype" list and confirm that it removes seqs from all
  # from all list types
  expect_equal(get_bin_abundance(dataset_t, "Phylo05", "phylotype"), 5337)
  expect_equal(get_bin_abundance(dataset_t, "Phylo06", "phylotype"), 715)

  phylo05 <- get_bin(dataset_t, "Phylo05", "phylotype")
  expect_equal(length(.split_at_char(phylo05)), 54)

  phylo06 <- get_bin(dataset_t, "Phylo06", "phylotype")
  expect_equal(length(.split_at_char(phylo06)), 47)

  remove_bins(
    dataset_t,
    c("Phylo05", "Phylo06"),
    c("test", "test"),
    "phylotype"
  )

  expect_equal(dataset_t$get_num_bins("phylotype"), 61)
  expect_equal(dataset_t$get_num_bins("otu"), 512)
  expect_equal(dataset_t$get_num_bins("asv"), 2324)
  expect_equal(dataset_t$get_num_sequences(), 107911)
  expect_equal(dataset_t$get_num_sequences(TRUE), 2324)
  expect_equal(dataset_t$get_num_treatments(), 2)
  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_sequences(sample = "F3D0"), 5977)
  expect_equal(dataset_t$get_num_sequences(sample = "F3D1"), 4467)
  # note that the number of seqs removed will be less that 297+266 because
  # some seqs are assigned to both samples and some seqs will be present in
  # other samples
  expect_equal(dataset_t$get_num_sequences(TRUE, "F3D0"), 297)
  expect_equal(dataset_t$get_num_sequences(TRUE, "F3D1"), 266)

  # remove samples
  remove_samples(dataset_t, c("F3D0", "F3D1"))
  expect_equal(dataset_t$get_num_treatments(), 2)
  expect_equal(dataset_t$get_num_samples(), 17)
  expect_equal(dataset_t$get_num_bins("phylotype"), 57)
  expect_equal(dataset_t$get_num_bins("otu"), 482)
  expect_equal(dataset_t$get_num_bins("asv"), 2124)
  expect_equal(dataset_t$get_num_sequences(), 97467)
  expect_equal(dataset_t$get_num_sequences(TRUE), 2124)
  expect_equal(dataset_t$get_num_sequences(TRUE, "F3D0"), 0)
  expect_equal(dataset_t$get_num_sequences(TRUE, "F3D1"), 0)
  expect_equal(dataset_t$get_num_sequences(sample = "F3D0"), 0)
  expect_equal(dataset_t$get_num_sequences(sample = "F3D1"), 0)

  # remove things just classified to bacteria, and things classified to
  # Bacteria;"Bacteroidetes"; with confidence less than 95
  remove_lineages(dataset_t, c(
    "Bacteria;Bacteria_unclassified;",
    "Bacteria(100);\"Bacteroidetes\"(95);"
  ))

  expect_equal(dataset_t$get_num_treatments(), 2)
  expect_equal(dataset_t$get_num_samples(), 17)
  expect_equal(dataset_t$get_num_bins("phylotype"), 57)
  expect_equal(dataset_t$get_num_bins("otu"), 475)
  expect_equal(dataset_t$get_num_bins("asv"), 2086)
  expect_equal(dataset_t$get_num_sequences(), 97428)
  expect_equal(dataset_t$get_num_sequences(TRUE), 2086)
})

test_that("dataset - intialize from dataset object", {
  temp <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    phylo_list = rdataset_example("final.tx.list"),
    asv_list = rdataset_example("final.asv.list"),
    dataset_name = "miseq_sop"
  )

  dataset_t <- dataset$new(
    name = "clone_of_miseq", dataset = temp,
    processors = 4
  )

  expect_equal(get_dataset_name(dataset_t), "clone_of_miseq")
  expect_equal(dataset_t$get_num_sequences(TRUE), 2425)
  expect_equal(dataset_t$get_num_sequences(), 113963)
  expect_equal(dataset_t$get_num_treatments(), 2)
  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_bins("otu"), 531)
  expect_equal(dataset_t$get_num_bins("phylotype"), 63)
  expect_equal(dataset_t$get_num_bins("asv"), 2425)
  expect_equal(get_num_processors(dataset_t), 4)
})

test_that("dataset - addSeqs, assign samples", {
  data <- dataset$new("mydata")

  # missing data and names
  expect_error(data$add_sequences())

  fasta_data <- read_fasta(rdataset_example("final.fasta"))
  names(fasta_data) <- c("myNameTag", "mySeqTag")

  expect_error(add_sequences(data, fasta_data, NULL, "names", "mySeqTag"))
  expect_error(add_sequences(data, "not_a_data.frame"))

  add_sequences(data, fasta_data, NULL, "myNameTag", "mySeqTag")

  expect_equal(data$get_num_sequences(), 2425)

  clear(data)

  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")
  comments <- c("ddd", "ftf", "efr", "ssd")

  fasta_data <- data.frame(names = names, seqs = seqs, comments = comments)

  expect_error(add_sequences(data, fasta_data,
    sequence_names = "names",
    sequences = "sequences"
  ))
  add_sequences(data, fasta_data, NULL, "names", "seqs")

  expect_equal(data$get_num_sequences(), 4)
  expect_equal(get_sequences(data), seqs)
  clear(data)

  expect_error(add_sequences(data, fasta_data,
    sequences = "seqs",
    comments = "comments23"
  ))
  add_sequences(data, fasta_data, NULL, "names", "seqs", "comments")

  expect_equal(data$get_num_sequences(), 4)
  expect_equal(get_sequences(data), seqs)

  clear(data)
  add_sequences(data, fasta_data, NULL, "names", "", "comments")

  expect_equal(data$get_num_sequences(), 4)
  expect_equal(get_sequences(data), rep("", 4))
  clear(data)

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

  # include reference
  url <- "https://mothur.org/wiki/silva_reference_files/"
  add_sequences(
    data,
    data.frame(sequence_names = names, sequences = seqs),
    new_reference(
      "silva.bacteria.fasta", "1.38.1", "",
      "alignment by mothur2 v1.0", url
    )
  )

  references <- get_references(data)

  expect_equal(nrow(references), 1)
  expect_equal(references[[1, "reference_names"]], "silva.bacteria.fasta")
  expect_equal(references[[1, "reference_versions"]], "1.38.1")
  expect_equal(references[[1, "reference_usages"]], "NA")
  expect_equal(references[[1, "reference_notes"]], "alignment by mothur2 v1.0")
  expect_equal(references[[1, "reference_urls"]], url)

  assign_sequence_abundance(
    data, data.frame(
      sequence_names = ids, abundances = abundances,
      samples = samples, treatments = treatments
    )
  )

  # assign bins
  bins <- c("bin1", "bin2", "bin1", "bin2")
  assign_bins(data, data.frame(bin_names = bins, sequence_names = names))

  expect_equal(data$get_num_bins(), 2)
  expect_equal(get_list(data)$otu_id, c("bin1", "bin1", "bin2", "bin2"))
  expect_equal(get_list(data)$seq_id, c("seq1", "seq3", "seq2", "seq4"))

  # get asv generated by dataset
  expect_equal(get_list(data, "asv")$asv_id, c("ASV1", "ASV2", "ASV3", "ASV4"))
  expect_equal(get_list(data, "asv")$seq_id, c("seq1", "seq2", "seq3", "seq4"))

  expect_equal(get_rabund(data)$otu_id, c("bin1", "bin2"))
  expect_equal(get_rabund(data)$abundance, c(1200, 120))

  expect_equal(get_bin_assignments(data)$bin_names, c(
    "bin1", "bin1", "bin1",
    "bin2", "bin2", "bin2"
  ))
  expect_equal(get_bin_assignments(data)$samples, c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  ))
  expect_equal(get_bin_assignments(data)$abundances, c(
    275, 425, 500,
    26, 40, 54
  ))

  expect_true(is_aligned(data))
  expect_equal(get_sequence_names(data), names)
  expect_equal(get_sequences(data), seqs)
  expect_equal(get_sequence_names(data, "sample2"), names)
  expect_equal(get_sequence_names(data, "sample3"), c("seq1", "seq2", "seq3"))
  expect_equal(get_sequence_names(data, "sample4"), c("seq1", "seq2", "seq4"))
  expect_equal(get_sequences(data, "sample2"), seqs)
  expect_equal(get_sequences(data, "sample3"), c("ATTGC", "ATTGC", "ATTGC"))
  expect_equal(get_sequences(data, "sample4"), c("ATTGC", "ATTGC", "ATTGC"))

  expect_equal(data$get_num_samples(), 3)
  expect_equal(data$get_num_treatments(), 2)
  expect_equal(data$get_treatments(), c("early", "late"))
  expect_equal(data$get_samples(), c("sample2", "sample3", "sample4"))

  sample_summary <- data$get_summary(TRUE)$sample_summary
  treatment_summary <- data$get_summary(TRUE)$treatment_summary

  expect_equal(treatment_summary$total, c(766, 554))
  expect_equal(sample_summary$total, c(301, 465, 554))

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

test_that("dataset - assign_sequence_abundance, remove_sequences", {
  data <- dataset$new("mydata")

  # missing data and names
  expect_error(data$assign_sequence_abundance())

  sequence_abundance <- readr::read_tsv(rdataset_example(
    "mothur2_count_table.tsv"
  ), show_col_types = FALSE)


  assign_sequence_abundance(data, sequence_abundance, "names")

  expect_equal(data$get_num_sequences(TRUE), 2425)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_treatments(), 2)

  names(sequence_abundance) <- c("ids", "abunds", "groups", "time")

  expect_error(assign_sequence_abundance(data, sequence_abundance, "ids"))

  # no treatments
  assign_sequence_abundance(data, sequence_abundance, "ids", "abunds", "groups")

  expect_equal(data$get_num_sequences(TRUE), 2425)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_treatments(), 0)

  assign_sequence_abundance(
    data, sequence_abundance, "ids", "abunds", "groups",
    "time"
  )

  expect_equal(data$get_num_sequences(TRUE), 2425)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_treatments(), 2)

  clear(data)

  names <- c("seq1", "seq2", "seq3", "seq4")
  abunds <- c(10, 20, 30)

  expect_error(assign_sequence_abundance(
    data = NULL,
    sequence_names = names,
    abundances = abunds
  ))

  clear(data)

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

  assign_sequence_abundance(data, data.frame(
    sequence_names = names,
    abundances = rabunds
  ))

  expect_equal(data$get_num_sequences(), 1269)
  expect_equal(data$get_num_sequences(TRUE), 4)
  expect_equal(data$get_num_samples(), 0)
  expect_equal(data$get_num_treatments(), 0)

  missing_id <- c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2",
    "seq3",
    "seq3"
  )

  assign_sequence_abundance(data, data.frame(
    sequence_names = ids,
    abundances = abundances,
    samples = groups
  ))

  expect_equal(data$get_num_sequences(), 1269)
  expect_equal(data$get_num_sequences(TRUE), 4)
  expect_equal(data$get_num_samples(), 3)
  expect_equal(data$get_num_treatments(), 0)

  assign_sequence_abundance(
    data, data.frame(
      sequence_names = ids, abundances = abundances,
      samples = groups, treatments = treatments
    )
  )
  expect_equal(data$get_num_treatments(), 2)

  remove_sequences(data, seqs_to_remove, trash_codes)

  expect_equal(data$get_num_sequences(), 29)
  expect_equal(data$get_num_sequences(TRUE), 2)
  expect_equal(data$get_num_samples(), 2)
  expect_equal(data$get_num_treatments(), 1)
})

test_that("dataset - get_list get_rabund, get_bin_assignments", {
  dataset_t <- new_dataset("my_dataset")

  expect_error(assign_bins(dataset_t))

  seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  sequence_abundances <- c(10, 100, 1, 500, 25, 80)

  assign_bins(dataset_t, data.frame(
    bin_names = bin_ids,
    abundances = sequence_abundances,
    sequence_names = seq_ids
  ))
  # bins would look like:
  # label  bin1             bin2        bin3 ...
  # 0.03   seq1,seq2,seq4   seq3,seq6   seq5 ...
  # 0.03   110              525         80 ...

  list <- get_list(dataset_t)

  expect_equal(list$otu_id, bin_ids)
  expect_equal(list$seq_id, seq_ids)
  expect_equal(get_list(dataset_t, "non_existance_bin_type"), data.frame())
  expect_equal(
    get_list_vector(dataset_t, "non_existance_bin_type"),
    character()
  )

  rabund <- get_rabund(dataset_t)

  abunds <- c(111, 525, 80)
  expect_equal(rabund$otu_id, unique(bin_ids))
  expect_equal(rabund$abundance, abunds)

  expect_equal(
    get_rabund_vector(dataset_t, "non_existance_bin_type"),
    numeric()
  )
  expect_equal(
    get_rabund_vector(dataset_t, "otu"),
    abunds
  )

  dataset_t <- dataset$new("my_dataset")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  samples <- c(
    "sample1", "sample2", "sample5",
    "sample1", "sample3", "sample1"
  )
  sample_abundances <- c(10, 100, 1, 500, 25, 80)
  assign_bins(
    dataset_t,
    data.frame(
      bin_names = bin_ids, abundances = sample_abundances,
      samples = samples
    )
  )

  shared <- get_bin_assignments(dataset_t)
  expect_equal(shared$bin_names, bin_ids)
  expect_equal(shared$abundances, sample_abundances)
  expect_equal(shared$samples, samples)
  expect_equal(dataset_t$get_num_bins(), 3)

  expect_equal(
    get_bin_assignments(dataset_t, "non_existance_bin_type"),
    data.frame()
  )
  expect_equal(get_shared_vector(dataset_t, "otu"), list(
    c(10, 100, 0, 1),
    c(500, 0, 25, 0),
    c(80, 0, 0, 0)
  ))
  expect_equal(length(get_shared_vector(dataset_t, "bad_type")), 0)

  expect_true(has_sample(dataset_t, "sample1"))
  expect_false(has_sample(dataset_t, "non_existant_sample"))

  dataset_t <- dataset$new("my_dataset")
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

  assign_bins(
    dataset_t,
    data.frame(
      bin_names = bin_ids, abundances = sample_abundances,
      samples = samples, sequence_names = seq_ids
    )
  )

  expect_equal(dataset_t$get_num_bins(), 3)
  expect_equal(dataset_t$get_num_samples(), 6)
  expect_equal(
    dataset_t$get_summary(TRUE)[["sample_summary"]]$total,
    c(36, 25, 2, 20, 13, 4)
  )

  clear(dataset_t)

  bin_table <- readr::read_tsv(rdataset_example(
    "mothur2_bin_assignments_shared.tsv"
  ), show_col_types = FALSE)

  expect_error(assign_bins(dataset_t, bin_table, "otu", NULL, "id"))

  assign_bins(dataset_t, bin_table)

  expect_equal(dataset_t$get_num_bins(), 531)
  expect_equal(dataset_t$get_num_samples(), 19)

  clear(dataset_t)

  bin_table <- readr::read_tsv(rdataset_example(
    "mothur2_bin_assignments_list.tsv"
  ), show_col_types = FALSE)

  assign_bins(
    dataset_t, bin_table, "otu", NULL,
    "otu_id", "", "", "seq_id"
  )

  expect_equal(dataset_t$get_num_bins(), 531)
  expect_equal(dataset_t$get_num_sequences(), 2425)
})

# assign_sequence_taxonomy, get_sequence_taxonomy_report
test_that("dataset - ", {
  # same length with confidences
  names <- c("seq1", "seq2", "seq3", "seq4")

  tax1 <- "Bacteria(100);Bacteroidetes(95);Bacteroidia(90);"
  tax2 <- "Bacteria(100);Proteobacteria(89);Betaproteobacteria(85);"
  tax3 <- "Bacteria(100);Firmicutes(99);Bacilli(90);"
  tax4 <- "Bacteria(100);Proteobacteria(87);Gammaproteobacteria(82);"
  taxonomies <- c(tax1, tax2, tax3, tax4)

  dataset_t <- dataset$new("my_dataset")

  url <- paste0(
    "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset",
    "9_032012.pds.zip"
  )

  # assign taxonomy with reference
  assign_sequence_taxonomy(
    dataset_t,
    data.frame(sequence_names = names, taxonomies = taxonomies),
    new_reference(
      "trainset9_032012.pds.zip", "9_032012",
      "classification by mothur2 v1.0 using default options", "",
      url
    )
  )

  references <- get_references(dataset_t)

  note <- "classification by mothur2 v1.0 using default options"
  expect_equal(nrow(references), 1)
  expect_equal(references[[1, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[1, 2]], "9_032012")
  expect_equal(references[[1, 3]], note)
  expect_equal(references[[1, 4]], "NA")
  expect_equal(references[[1, 5]], url)

  report <- get_sequence_taxonomy_report(dataset_t)

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

  dataset_t <- dataset$new("my_dataset")
  assign_sequence_taxonomy(dataset_t, data.frame(
    sequence_names = names,
    taxonomies = taxonomies
  ))
  report <- get_sequence_taxonomy_report(dataset_t)

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

  dataset_t <- dataset$new("my_dataset")
  assign_sequence_taxonomy(dataset_t, data.frame(
    sequence_names = names,
    taxonomies = taxonomies
  ))
  report <- get_sequence_taxonomy_report(dataset_t)

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

  dataset_t <- dataset$new("my_dataset")
  assign_sequence_taxonomy(dataset_t, data.frame(
    sequence_names = names,
    taxonomies = taxonomies
  ))
  report <- get_sequence_taxonomy_report(dataset_t)

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
  assign_bins(dataset_t, data.frame(bin_names = bins, sequence_names = names))

  report <- get_bin_taxonomy_report(dataset_t)

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

  clear(dataset_t)
  bin_ids <- c("bin1", "bin2", "bin3", "bin4")
  taxonomies <- c(
    "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;"
  )

  expect_error(assign_bin_taxonomy(dataset_t, data = "not_a_data.frame"))
  expect_equal(get_reports(dataset_t), list())
  expect_equal(dataset_t$get_summary(), list())
  expect_false(dataset_t$has_sample("noSample"))

  abunds <- c(200, 40, 100, 5)
  assign_bins(dataset_t, data.frame(bin_names = bin_ids, abundances = abunds))

  url <- paste0(
    "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset",
    "9_032012.pds.zip"
  )

  assign_bin_taxonomy(
    dataset_t,
    data.frame(bin_names = bin_ids, taxonomies = taxonomies),
    "otu",
    new_reference(
      "trainset9_032012.pds.zip", "9_032012", "",
      "classification by mothur2 v1.0", url
    )
  )

  references <- get_references(dataset_t)

  expect_equal(nrow(references), 1)
  expect_equal(references[[1, "reference_names"]], "trainset9_032012.pds.zip")
  expect_equal(references[[1, "reference_versions"]], "9_032012")
  expect_equal(references[[1, "reference_usages"]], "NA")
  expect_equal(
    references[[1, "reference_notes"]],
    "classification by mothur2 v1.0"
  )
  expect_equal(references[[1, "reference_urls"]], url)

  report <- get_bin_taxonomy_report(dataset_t)

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


  assign_bin_taxonomy(dataset_t, data.frame(
    bin_names = bin_ids,
    taxonomies = taxonomies
  ))

  report <- get_bin_taxonomy_report(dataset_t)

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

  expect_equal(get_sequence_taxonomy_report(dataset_t), data.frame())
  expect_equal(get_bin_assignments(dataset_t), data.frame())

  clear(dataset_t)

  expect_error(dataset_t$assign_bin_taxonomy())

  table <- readr::read_tsv(rdataset_example("final.cons.taxonomy"),
    show_col_types = FALSE
  )
  bin_table <- readr::read_tsv(rdataset_example(
    "mothur2_bin_assignments_list.tsv"
  ), show_col_types = FALSE)

  assign_bins(dataset_t, bin_table, "otu", NULL, "otu_id", "", "", "seq_id")

  assign_bin_taxonomy(dataset_t, table, "otu", NULL, "OTU", "Taxonomy")

  table <- get_bin_taxonomy_report(dataset_t)

  expect_equal(table[2758, 1], "Otu460")
  expect_equal(table[2758, 3], "\"Bacteroidales\"")
  expect_equal(table[2758, 2], 4)

  expect_equal(table[2881, 1], "Otu481")
  expect_equal(table[2881, 3], "Bacteria")
  expect_equal(table[2881, 2], 1)
})

test_that("dataset - add_metadata, get_metadata", {
  dataset_t <- dataset$new("my_dataset")

  expect_equal(get_metadata(dataset_t), data.frame())

  metadata <- readr::read_tsv(rdataset_example("sample-metadata.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )

  add_metadata(dataset_t, metadata)
  metadata <- get_metadata(dataset_t)

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

  clear(dataset_t, "metadata")
  metadata <- get_metadata(dataset_t)
  expect_equal(nrow(metadata), 0)
})

test_that("dataset - add_references, get_references", {
  dataset_t <- dataset$new("my_dataset")

  expect_equal(get_references(dataset_t), data.frame())
  expect_error(add_references(dataset_t, reference = c("bad_type")))
  expect_error(add_references(dataset_t, data.frame()))
  expect_error(add_references(dataset_t))

  references <- get_references(dataset_t)
  expect_equal(nrow(references), 0)

  reference <- readr::read_csv(rdataset_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )

  mothur_url <- "https://github.com/mothur/mothur/releases/tag/v1.48.2"

  # add data.frame and single reference at the same time
  add_references(dataset_t, reference)

  references <- get_references(dataset_t)

  # random spot checks
  expect_equal(nrow(references), 2)
  expect_equal(references[[1, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[2, 1]], "silva.v4.fasta")
  expect_equal(references[[1, 4]], "NA")
  expect_equal(
    references[[2, 4]],
    "custom reference created by trimming silva.bacteria.fasta to the V4 region"
  )

  dataset_t <- dataset$new("my_dataset")

  # add single reference then dataframe
  ref <- data.frame(
    reference_name = "mothur software package",
    reference_version = "1.48.2",
    reference_usage = "analysis of dataset",
    reference_note = "This is my mothur note",
    reference_url = mothur_url
  )

  add_references(
    dataset_t, ref, "reference_name", "reference_version",
    "reference_usage", "reference_note", "reference_url"
  )

  references <- get_references(dataset_t)
  expect_equal(nrow(references), 1)

  add_references(dataset_t, reference)

  references <- get_references(dataset_t)
  expect_equal(nrow(references), 3)

  expect_equal(references[[2, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[3, 1]], "silva.v4.fasta")
  expect_equal(references[[1, 1]], "mothur software package")
  expect_equal(references[[2, 2]], "NA")
  expect_equal(references[[3, 2]], "1.38.1")
  expect_equal(references[[1, 4]], "This is my mothur note")

  dataset_t$clear("bad_type")
  expect_equal(nrow(get_references(dataset_t)), 3)

  clear(dataset_t, c("references"))
  expect_equal(nrow(get_references(dataset_t)), 0)
})

test_that("dataset - add_alignment_report, get_alignment_report", {
  dataset_t <- dataset$new("my_dataset")

  align_report <- readr::read_tsv(rdataset_example("alignment_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )

  expect_error(add_report(dataset_t, align_report, "align_report", "badName"))
  add_report(dataset_t, align_report, "align_report", "QueryName")

  align_report <- get_reports(dataset_t)[["align_report"]]

  # random spot checks
  expect_equal(nrow(align_report), 5)
  expect_equal(align_report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))
  expect_equal(align_report[[2, 2]], 253)
  expect_equal(align_report[[3, 3]], "AF132257.1")
  expect_equal(align_report[[1, 4]], 293)
  expect_equal(align_report[[5, 8]], 1)
  expect_equal(align_report[[4, 4]], 292)

  clear(dataset_t, "reports")
  expect_equal(get_reports(dataset_t), list())

  # no report added because of missing entries
  add_sequences(dataset_t, data.frame(sequence_names = c("seq6", "seq7")))
  add_report(dataset_t, align_report, "align_report", "QueryName")
  expect_equal(get_reports(dataset_t), list())
})

test_that("dataset - add / get _contigs_assembly_report,", {
  dataset_t <- dataset$new("my_dataset")

  report <- readr::read_tsv(rdataset_example("contigs_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  expect_error(add_report(dataset_t, report, "contigs_report", "badName"))

  add_report(dataset_t, report, "contigs_report", "Name")

  report <- get_reports(dataset_t)[["contigs_report"]]

  # random spot checks
  expect_equal(nrow(report), 5)
  expect_equal(report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))
  expect_equal(report[[2, 2]], 252)
  expect_equal(report[[3, 3]], 249)
  expect_equal(report[[1, 4]], 2)
  expect_equal(round(report[[5, 8]], digits = 4), 0.0257)
  expect_equal(report[[4, 4]], 2)

  clear(dataset_t, "reports")
  expect_equal(length(get_reports(dataset_t)), 0)
  expect_equal(dataset_t$get_num_sequences(), 5)
  add_report(dataset_t, report, "contigs_report", "Name")

  report <- get_reports(dataset_t)[["contigs_report"]]

  expect_equal(nrow(report), 5)
  expect_equal(report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))

  # no report added because of missing entries
  clear(dataset_t)
  add_sequences(dataset_t, data.frame(sequence_names = c("seq6", "seq7")))
  add_report(dataset_t, report, "contigs_report", "Name")
  expect_equal(length(get_reports(dataset_t)), 0)
})

test_that("dataset - add / get _chimera_report,", {
  dataset_t <- dataset$new("my_dataset")

  report <- readr::read_tsv(rdataset_example("chimera_report.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  expect_error(add_report(dataset_t, report, "chimera_report", "badName"))

  add_report(dataset_t, report, "chimera_report", "Query")

  report <- get_reports(dataset_t)[["chimera_report"]]

  # random spot checks
  expect_equal(nrow(report), 71)
  expect_equal(report[, 2], get_sequence_names(dataset_t))
  expect_equal(report[[8, 5]], 82.7)
  expect_equal(report[[8, 17]], "N")
  expect_equal(report[[67, 17]], "Y")

  clear(dataset_t, "reports")
  expect_equal(length(get_reports(dataset_t)), 0)
  expect_equal(dataset_t$get_num_sequences(), 71)
  add_report(dataset_t, report, "chimera_report", "Query")

  report <- get_reports(dataset_t)[["chimera_report"]]

  expect_equal(nrow(report), 71)
  expect_equal(report[, 2], get_sequence_names(dataset_t))

  chimera_summary <- dataset_t$get_summary()[["chimera_report"]]

  expect_equal(ncol(chimera_summary), 13)

  # no report added because of missing entries
  clear(dataset_t, "reports")
  add_sequences(dataset_t, data.frame(sequence_names = c("seq6", "seq7")))
  add_report(dataset_t, report, "chimera_report", "Query")
  expect_equal(length(get_reports(dataset_t)), 0)
})

test_that("dataset - get_sequence_summary,", {
  dataset_t <- dataset$new("my_dataset")

  report <- readr::read_tsv(rdataset_example("contigs_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  add_report(dataset_t, report, "contigs_report", "Name")

  report <- readr::read_tsv(rdataset_example("alignment_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  add_report(dataset_t, report, "alignment_report", "QueryName")

  summary <- dataset_t$get_summary()

  expect_equal(summary$contigs_report$MisMatches, c(0, 0, 1, 2, 7, 7, 7, 2))
  expect_equal(summary$contigs_report$Overlap_End, rep(251, 8))
  expect_equal(summary$contigs_report$Length[1], 252)
  expect_equal(summary$contigs_report$Length[7], 253)

  expect_equal(summary$alignment_report$QueryLength, c(
    252, 252, 253, 253,
    253, 253, 253, 252.6
  ))
  expect_equal(summary$alignment_report$GapsInQuery, rep(0, 8))
  expect_equal(summary$alignment_report$SearchScore[1], 57.55)
  expect_equal(summary$alignment_report$SearchScore[7], 82.44)
})

test_that("dataset - add_sequence_tree / get_sequence_tree,", {
  # create tree from sequences
  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ACTGC", "ATTCC", "GTTGC", "ATGGC")
  dataset_t <- dataset$new()

  expect_error(dataset_t$add_sequence_tree(tree = c("bad_type")))

  add_sequences(dataset_t, data.frame(sequence_names = names, sequences = seqs))
  expect_equal(dataset_t$get_sequence_tree(), NULL)

  # add full tree
  dataset_t <- dataset$new()
  add_sequences(dataset_t, data.frame(sequence_names = names))

  l <- lapply(strsplit(seqs, split = ""), "[")
  names(l) <- names

  # should alert that the tree is missing reads, and not save it
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(get_sequence_names(dataset_t)), sort(tree$tip.label))
  expect_equal(tree$edge[, 1], c(5, 5, 5, 6, 6))
  expect_equal(tree$edge[, 2], c(4, 3, 6, 1, 2))
  expect_equal(
    round(tree$edge.length, digits = 2),
    c(0.26, 0.33, 0.07, 0.33, 0.26)
  )

  # remove seq and make sure tree is pruned as well
  remove_sequences(dataset_t, c("seq1"), c("bad"))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(get_sequence_names(dataset_t)), sort(tree$tip.label))
  expect_equal(tree$edge[, 1], c(4, 4, 4))
  expect_equal(tree$edge[, 2], c(3, 2, 1))
  expect_equal(
    round(tree$edge.length, digits = 2),
    c(0.26, 0.33, 0.33)
  )

  # add tree with missing sequences
  l <- lapply(strsplit(seqs[1:3], split = ""), "[")
  names(l) <- names[1:3]

  # should alert that the tree is missing reads, and not save it
  dataset_t <- dataset$new()
  add_sequences(dataset_t, data.frame(sequence_names = names))
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))
  expect_equal(dataset_t$get_sequence_tree(), NULL)

  # add tree from file
  dataset_t <- dataset$new()
  dataset_t$add_sequence_tree(read.tree(rdataset_example("final.phylip.tre")))
  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(get_sequence_names(dataset_t)), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(2426, 2427, 2427, 2426, 2428))
  expect_equal(tree$edge[1:5, 2], c(2427, 1, 2, 2428, 2429))
  expect_equal(
    round(tree$edge.length[1:5], digits = 3),
    c(NaN, 0.004, 0.004, 0.000, 0.002)
  )

  # add tree with extra sequences, forcing prune
  dataset_t <- dataset$new()
  add_sequences(dataset_t, data.frame(sequence_names = names))

  seqs <- c(seqs, "ACTGC")
  names <- c(names, "seq5")

  # add tree with extra sequences, forcing prune
  l <- lapply(strsplit(seqs, split = ""), "[")
  names(l) <- names
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(get_sequence_names(dataset_t)), sort(tree$tip.label))
  expect_equal(tree$edge[, 1], c(5, 5, 5, 6, 6))
  expect_equal(tree$edge[, 2], c(4, 3, 6, 1, 2))
  expect_equal(
    round(tree$edge.length, digits = 2),
    c(0.26, 0.33, 0.07, 0.33, 0.26)
  )
})

test_that("dataset - add_sample_tree / get_sample_tree,", {
  dataset_t <- dataset$new()
  expect_null(dataset_t$get_sample_tree())

  sample_tree <- ape::read.tree(
    rdataset_example("final.opti_mcc.jclass.ave.tre")
  )

  # should report no samples and not save
  dataset_t$add_sample_tree(sample_tree)
  expect_null(dataset_t$get_sample_tree())

  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  expect_error(dataset_t$add_sample_tree(tree = c("bad_type")))

  sequence_tree <- ape::read.tree(rdataset_example("final.phylip.tre"))

  # should report missing samples since this is a sequence tree and not save
  dataset_t$add_sample_tree(sequence_tree)
  expect_null(dataset_t$get_sample_tree())

  dataset_t$add_sample_tree(sample_tree)

  tree <- dataset_t$get_sample_tree()

  expect_equal(sort(dataset_t$get_samples()), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(20, 21, 22, 23, 24))

  remove_samples(dataset_t, c("F3D1", "F3D141"))

  tree <- dataset_t$get_sample_tree()

  expect_equal(sort(dataset_t$get_samples()), sort(tree$tip.label))

  # add tree with all groups, prune tree on add
  dataset_t$add_sample_tree(sample_tree)

  tree <- dataset_t$get_sample_tree()

  # confirm pruning
  expect_equal(sort(dataset_t$get_samples()), sort(tree$tip.label))
})

test_that("dataset - assign_treatments", {
  # create dataset without treatment assignments
  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    phylo_list = rdataset_example("final.tx.list"),
    asv_list = rdataset_example("final.asv.list"),
    dataset_name = "miseq_sop"
  )

  design_table <- readr::read_tsv(
    rdataset_example(
      "mouse.time.design"
    ),
    show_col_types = FALSE
  )

  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_treatments(), 0)

  expect_error(dataset_t$assign_treatments(design_table, samples = "group"))
  expect_error(dataset_t$assign_treatments(design_table, treaments = "time"))
  expect_error(dataset_t$assign_treatments())
  expect_error(dataset_t$assign_treatments(data = NULL, samples = "not_valid"))
  expect_error(dataset_t$assign_treatments("not_a_data.frame"))

  # test with data.frame
  assign_treatments(dataset_t, design_table)

  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_treatments(), 2)

  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    phylo_list = rdataset_example("final.tx.list"),
    asv_list = rdataset_example("final.asv.list"),
    dataset_name = "miseq_sop"
  )

  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_treatments(), 0)

  # test with samples and treatments
  assign_treatments(dataset_t, design_table)

  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_treatments(), 2)
})

test_that("dataset - assign_sequence_taxonomy", {
  # create dataset without taxonomy assignments
  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  tax_table <- read_mothur_taxonomy(rdataset_example("final.taxonomy"))

  # no taxonomies yet
  expect_equal(get_sequence_taxonomy_report(dataset_t), data.frame())

  expect_error(assign_sequence_taxonomy(
    dataset_t, tax_table, NULL,
    "not_valid"
  ))
  expect_error(assign_sequence_taxonomy(
    dataset_t, tax_table, NULL,
    "sequence_names", "not_valid"
  ))
  expect_error(assign_sequence_taxonomy(dataset_t))

  # test with data.frame
  assign_sequence_taxonomy(dataset_t, tax_table)

  report <- get_sequence_taxonomy_report(dataset_t)

  expect_equal(nrow(get_sequence_taxonomy_report(dataset_t)), 14550)

  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  expect_equal(get_sequence_taxonomy_report(dataset_t), data.frame())

  # test with samples and treatments
  assign_sequence_taxonomy(dataset_t, tax_table)

  report <- get_sequence_taxonomy_report(dataset_t)

  expect_equal(nrow(get_sequence_taxonomy_report(dataset_t)), 14550)
})

test_that("dataset - export,", {
  miseq <- miseq_sop_example()

  miseq_table <- export_dataset(miseq)

  table_names <- c(
    "sequence_data", "sequence_report",
    "sequence_abundance_table", "otu_bin_data",
    "otu_sequence_bin_assignments", "otu_bin_representative_sequences",
    "asv_bin_data",
    "asv_sequence_bin_assignments", "phylotype_bin_data",
    "phylotype_sequence_bin_assignments", "references", "metadata",
    "sequence_tree", "sample_tree"
  )

  sequence_data_names <- c(
    "sequence_ids", "sequence_names",
    "sequences", "taxonomies", "include_sequence"
  )

  sequence_report_names <- c(
    "sequence_ids", "starts", "ends", "lengths",
    "ambigs", "longest_homopolymers", "num_ns"
  )

  sequence_at_names <- c(
    "sequence_ids", "abundances",
    "samples", "treatments"
  )

  bin_data_names <- c(
    "bin_ids", "bin_names", "abundances",
    "taxonomies", "include_bin"
  )

  bin_assignment_names <- c("bin_ids", "sequence_ids")

  metadata_names <- c("sample", "days_post_wean")

  references_names <- c(
    "reference_names", "reference_versions",
    "reference_usages", "reference_notes",
    "reference_urls"
  )

  expect_equal(names(miseq_table), table_names)
  expect_equal(nrow(miseq_table$sequence_data), 2425)
  expect_equal(nrow(miseq_table$sequence_report), 2425)
  expect_equal(nrow(miseq_table$sequence_abundance_table), 5539)

  expect_equal(nrow(miseq_table$otu_bin_data), 531)
  expect_equal(nrow(miseq_table$otu_sequence_bin_assignments), 2425)
  expect_equal(nrow(miseq_table$asv_bin_data), 2425)
  expect_equal(nrow(miseq_table$asv_sequence_bin_assignments), 2425)
  expect_equal(nrow(miseq_table$phylotype_bin_data), 63)
  expect_equal(nrow(miseq_table$phylotype_sequence_bin_assignments), 2425)
  expect_equal(nrow(miseq_table$metadata), 19)
  expect_equal(nrow(miseq_table$references), 2)

  expect_equal(names(miseq_table$sequence_data), sequence_data_names)
  expect_equal(names(miseq_table$sequence_report), sequence_report_names)
  expect_equal(names(miseq_table$sequence_abundance_table), sequence_at_names)

  expect_equal(names(miseq_table$otu_bin_data), bin_data_names)
  expect_equal(names(miseq_table$asv_bin_data), bin_data_names)
  expect_equal(names(miseq_table$phylotype_bin_data), bin_data_names)

  expect_equal(
    names(miseq_table$otu_sequence_bin_assignments),
    bin_assignment_names
  )
  expect_equal(
    names(miseq_table$asv_sequence_bin_assignments),
    bin_assignment_names
  )
  expect_equal(
    names(miseq_table$phylotype_sequence_bin_assignments),
    bin_assignment_names
  )

  expect_equal(names(miseq_table$metadata), metadata_names)
  expect_equal(names(miseq_table$references), references_names)

  expect_equal(
    miseq_table$sequence_data$sequence_names,
    get_sequence_names(miseq)
  )
  expect_equal(
    miseq_table$sequence_data$sequences,
    get_sequences(miseq)
  )
})

test_that("dataset - assign_bin_representative_sequences", {
  # create dataset sequences and shared data
  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    otu_shared = rdataset_example("final.opti_mcc.shared"),
    dataset_name = "miseq_sop"
  )

  # select first 531 seqs to be the representatives
  num_bins <- dataset_t$get_num_bins()
  rep_names <- get_sequence_names(dataset_t)[1:num_bins]
  bin_names <- dataset_t$get_bin_names()

  assign_bin_representative_sequences(
    dataset_t,
    data.frame(
      bin_names = bin_names,
      sequence_names = rep_names
    )
  )

  df <- get_bin_representative_sequences(dataset_t)

  expect_equal(df[[1]], bin_names)
  expect_equal(df[[2]], rep_names)
  expect_equal(df[[3]], get_sequences(dataset_t)[1:num_bins])

  # create dataset only shared data, this forces assign_bin_reps to add seqs
  dataset_t <- read_mothur(
    otu_shared = rdataset_example("final.opti_mcc.shared"),
    dataset_name = "miseq_sop"
  )

  assign_bin_representative_sequences(
    dataset_t,
    data.frame(
      bin_names = bin_names,
      sequence_names = rep_names
    )
  )

  df <- get_bin_representative_sequences(dataset_t)
  expect_equal(df[[1]], bin_names)
  expect_equal(df[[2]], rep_names)
  expect_equal(df[, 3], rep("", num_bins))

  expect_error(assign_bin_representative_sequences(
    dataset_t,
    data.frame(
      bin_names = bin_names,
      sequence_names = c("not", "enough", "sequence", "names")
    )
  ))

  d <- dataset$new()
  expect_equal(get_bin_representative_sequences(d), data.frame())
  expect_equal(get_sample_treatment_assignments(d), data.frame())

  dataset_t <- read_mothur(
    count = rdataset_example("final.count_table"),
    dataset_name = "miseq_sop"
  )
  expect_equal(get_sample_treatment_assignments(dataset_t), data.frame())
})
