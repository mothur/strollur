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

  expect_equal(dataset_t$get_dataset_name(), "miseq_sop")
  expect_equal(dataset_t$get_num_sequences(TRUE), 2425)
  expect_equal(dataset_t$get_num_sequences(), 113963)
  expect_equal(dataset_t$get_num_treatments(), 2)
  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_bins("otu"), 531)
  expect_equal(dataset_t$get_num_bins("phylotype"), 63)
  expect_equal(dataset_t$get_num_bins("asv"), 2425)

  seqs_summary <- dataset_t$get_sequence_summary()[["sequence_summary"]]
  print(dataset_t$get_sequence_summary())
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

  expect_equal(get_bin_abundance(dataset_t$data, "Phylo05", "phylotype"), 5337)
  expect_equal(get_bin_abundance(dataset_t$data, "Phylo06", "phylotype"), 715)

  phylo05 <- get_bin(dataset_t$data, "Phylo05", "phylotype")
  expect_equal(length(split_at_char(phylo05)), 54)

  phylo06 <- get_bin(dataset_t$data, "Phylo06", "phylotype")
  expect_equal(length(split_at_char(phylo06)), 47)

  remove_bins(
    dataset_t$data,
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
  dataset_t$remove_samples(c("F3D0", "F3D1"))
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
  dataset_t$remove_lineages(c(
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

  expect_equal(dataset_t$get_dataset_name(), "clone_of_miseq")
  expect_equal(dataset_t$get_num_sequences(TRUE), 2425)
  expect_equal(dataset_t$get_num_sequences(), 113963)
  expect_equal(dataset_t$get_num_treatments(), 2)
  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_bins("otu"), 531)
  expect_equal(dataset_t$get_num_bins("phylotype"), 63)
  expect_equal(dataset_t$get_num_bins("asv"), 2425)
  expect_equal(get_num_processors(dataset_t$data), 4)
})

test_that("dataset - addSeqs, assign samples", {
  data <- dataset$new("mydata")

  # missing data and names
  expect_error(data$add_sequences())

  fasta_data <- read_fasta(rdataset_example("final.fasta"))
  names(fasta_data) <- c("myNameTag", "mySeqTag")

  expect_error(data$add_sequences(fasta_data,
    sequence_names = "names",
    sequences = "mySeqTag"
  ))

  expect_error(data$add_sequences("not_a_data.frame"))

  data$add_sequences(fasta_data,
    sequence_names = "myNameTag",
    sequences = "mySeqTag"
  )

  expect_equal(data$get_num_sequences(), 2425)

  data$clear()

  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")
  comments <- c("ddd", "ftf", "efr", "ssd")

  fasta_data <- data.frame(names = names, seqs = seqs, comments = comments)

  expect_error(data$add_sequences(fasta_data,
    sequence_names = "names",
    sequences = "sequences"
  ))
  data$add_sequences(fasta_data,
    sequence_names = "names",
    sequences = "seqs"
  )

  expect_equal(data$get_num_sequences(), 4)
  expect_equal(data$get_sequences(), seqs)
  data$clear()

  expect_error(data$add_sequences(fasta_data,
    sequences = "seqs",
    comments = "comments23"
  ))
  data$add_sequences(fasta_data,
    sequence_names = "names",
    sequences = "seqs", comments = "comments"
  )

  expect_equal(data$get_num_sequences(), 4)
  expect_equal(data$get_sequences(), seqs)

  data$clear()
  data$add_sequences(fasta_data,
    sequence_names = "names",
    comments = "comments"
  )
  expect_equal(data$get_num_sequences(), 4)
  expect_equal(data$get_sequences(), rep("", 4))
  data$clear()

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
  data$add_sequences(
    sequence_names = names, sequences = seqs,
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


  data$assign_sequence_abundance(
    data = NULL, ids, abundances,
    samples, treatments
  )

  # assign bins
  bins <- c("bin1", "bin2", "bin1", "bin2")
  data$assign_bins(bin_names = bins, sequence_names = names)

  expect_equal(data$get_num_bins(), 2)
  expect_equal(data$get_list()$otu_id, c("bin1", "bin1", "bin2", "bin2"))
  expect_equal(data$get_list()$seq_id, c("seq1", "seq3", "seq2", "seq4"))

  # get asv generated by dataset
  expect_equal(data$get_list("asv")$asv_id, c("ASV1", "ASV2", "ASV3", "ASV4"))
  expect_equal(data$get_list("asv")$seq_id, c("seq1", "seq2", "seq3", "seq4"))

  expect_equal(data$get_rabund()$otu_id, c("bin1", "bin2"))
  expect_equal(data$get_rabund()$abundance, c(1200, 120))

  expect_equal(data$get_bin_assignments()$bin_names, c(
    "bin1", "bin1", "bin1",
    "bin2", "bin2", "bin2"
  ))
  expect_equal(data$get_bin_assignments()$samples, c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  ))
  expect_equal(data$get_bin_assignments()$abundances, c(
    275, 425, 500,
    26, 40, 54
  ))

  expect_true(data$is_aligned())
  expect_equal(data$get_sequence_names(), names)
  expect_equal(data$get_sequences(), seqs)
  expect_equal(data$get_sequence_names("sample2"), names)
  expect_equal(data$get_sequence_names("sample3"), c("seq1", "seq2", "seq3"))
  expect_equal(data$get_sequence_names("sample4"), c("seq1", "seq2", "seq4"))
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

test_that("dataset - assign_sequence_abundance, remove_sequences", {
  data <- dataset$new("mydata")

  # missing data and names
  expect_error(data$assign_sequence_abundance())

  sequence_abundance <- readr::read_tsv(rdataset_example(
    "mothur2_count_table.tsv"
  ), show_col_types = FALSE)
  data$assign_sequence_abundance(sequence_abundance)

  expect_equal(data$get_num_sequences(TRUE), 2425)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_treatments(), 2)

  names(sequence_abundance) <- c("ids", "abunds", "groups", "time")

  expect_error(data$assign_sequence_abundance(sequence_abundance,
    sequence_names = "ids"
  ))
  expect_error(data$assign_sequence_abundance(sequence_abundance,
    sequence_names = "ids",
    abundances = "abunds",
    samples = "samples"
  ))
  expect_error(data$assign_sequence_abundance(sequence_abundance,
    sequence_names = "ids",
    abundances = "abunds",
    treatments = "treatments"
  ))

  data$assign_sequence_abundance(sequence_abundance,
    sequence_names = "ids",
    abundances = "abunds",
    samples = "groups", treatments = NULL
  )

  expect_equal(data$get_num_sequences(TRUE), 2425)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_treatments(), 0)


  data$assign_sequence_abundance(sequence_abundance,
    sequence_names = "ids",
    abundances = "abunds",
    samples = "groups", treatments = "time"
  )

  expect_equal(data$get_num_sequences(TRUE), 2425)
  expect_equal(data$get_num_samples(), 19)
  expect_equal(data$get_num_treatments(), 2)

  data$clear()

  names <- c("seq1", "seq2", "seq3", "seq4")
  abunds <- c(10, 20, 30)

  expect_error(data$assign_sequence_abundance(
    data = NULL,
    sequence_names = names,
    abundances = abunds
  ))

  expect_error(data$assign_sequence_abundance(
    data = NULL,
    sequence_names = NULL,
    abundances = abunds
  ))

  expect_error(data$assign_sequence_abundance(
    data = NULL,
    sequence_names = names,
    abundances = NULL
  ))
  data$clear()

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

  data$assign_sequence_abundance(sequence_names = names, abundances = rabunds)

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

  data$assign_sequence_abundance(data = NULL, ids, abundances, groups)

  expect_equal(data$get_num_sequences(), 1269)
  expect_equal(data$get_num_sequences(TRUE), 4)
  expect_equal(data$get_num_samples(), 3)
  expect_equal(data$get_num_treatments(), 0)

  data$assign_sequence_abundance(
    data = NULL, ids, abundances,
    groups, treatments
  )
  expect_equal(data$get_num_treatments(), 2)

  remove_sequences(data$data, seqs_to_remove, trash_codes)

  expect_equal(data$get_num_sequences(), 29)
  expect_equal(data$get_num_sequences(TRUE), 2)
  expect_equal(data$get_num_samples(), 2)
  expect_equal(data$get_num_treatments(), 1)
})

test_that("dataset - get_list get_rabund, get_bin_assignments", {
  dataset_t <- dataset$new("my_dataset")

  expect_error(dataset_t$assign_bins())

  seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  sequence_abundances <- c(10, 100, 1, 500, 25, 80)

  expect_error(dataset_t$assign_bins(data = NULL))
  expect_error(dataset_t$assign_bins(data = NULL, bin_names = bin_ids))
  expect_error(dataset_t$assign_bins(
    data = NULL, bin_names = bin_ids,
    sequence_names = c("not enough seqs")
  ))

  dataset_t$assign_bins(
    bin_names = bin_ids,
    abundances = sequence_abundances,
    sequence_names = seq_ids
  )
  # bins would look like:
  # label  bin1             bin2        bin3 ...
  # 0.03   seq1,seq2,seq4   seq3,seq6   seq5 ...
  # 0.03   110              525         80 ...

  list <- dataset_t$get_list()

  expect_equal(list$otu_id, bin_ids)
  expect_equal(list$seq_id, seq_ids)
  expect_equal(dataset_t$get_list("non_existance_bin_type"), data.frame())
  expect_equal(
    get_list_vector(dataset_t$data, "non_existance_bin_type"),
    character()
  )

  rabund <- dataset_t$get_rabund()

  abunds <- c(111, 525, 80)
  expect_equal(rabund$otu_id, unique(bin_ids))
  expect_equal(rabund$abundance, abunds)

  expect_equal(
    get_rabund_vector(dataset_t$data, "non_existance_bin_type"),
    numeric()
  )
  expect_equal(
    get_rabund_vector(dataset_t$data, "otu"),
    abunds
  )

  dataset_t <- dataset$new("my_dataset")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  samples <- c(
    "sample1", "sample2", "sample5",
    "sample1", "sample3", "sample1"
  )
  sample_abundances <- c(10, 100, 1, 500, 25, 80)
  dataset_t$assign_bins(
    bin_names = bin_ids, abundances = sample_abundances,
    samples = samples
  )

  shared <- dataset_t$get_bin_assignments()
  expect_equal(shared$bin_names, bin_ids)
  expect_equal(shared$abundances, sample_abundances)
  expect_equal(shared$samples, samples)
  expect_equal(dataset_t$get_num_bins(), 3)

  expect_equal(
    get_bin_assignments(dataset_t$data, "non_existance_bin_type"),
    data.frame()
  )
  expect_equal(get_shared_vector(dataset_t$data, "otu"), list(
    c(10, 100, 0, 1),
    c(500, 0, 25, 0),
    c(80, 0, 0, 0)
  ))
  expect_equal(length(get_shared_vector(dataset_t$data, "bad_type")), 0)

  expect_true(has_sample(dataset_t$data, "sample1"))
  expect_false(has_sample(dataset_t$data, "non_existant_sample"))

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

  dataset_t$assign_bins(
    data = NULL, bin_ids, sample_abundances, samples,
    seq_ids
  )

  expect_equal(dataset_t$get_num_bins(), 3)
  expect_equal(dataset_t$get_num_samples(), 6)
  expect_equal(
    dataset_t$get_sample_summary(TRUE)[[1]]$total,
    c(36, 25, 2, 20, 13, 4)
  )

  dataset_t$clear()

  bin_table <- readr::read_tsv(rdataset_example(
    "mothur2_bin_assignments_shared.tsv"
  ), show_col_types = FALSE)

  expect_error(dataset_t$assign_bins(
    data = bin_table,
    bin_names = "id", samples = "sample"
  ))

  dataset_t$assign_bins(bin_table)

  expect_equal(dataset_t$get_num_bins(), 531)
  expect_equal(dataset_t$get_num_samples(), 19)

  dataset_t$clear()

  bin_table <- readr::read_tsv(rdataset_example(
    "mothur2_bin_assignments_list.tsv"
  ), show_col_types = FALSE)

  dataset_t$assign_bins(bin_table,
    bin_names = "otu_id",
    sequence_names = "seq_id"
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

  url <- paste("https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset",
    "9_032012.pds.zip",
    collapse = ""
  )

  # assign taxonomy with reference
  dataset_t$assign_sequence_taxonomy(
    sequence_names = names, taxonomies = taxonomies,
    reference_name = "trainset9_032012.pds.zip",
    reference_note = "classification by mothur2 v1.0 using default options",
    reference_version = "9_032012", reference_url = url
  )

  references <- dataset_t$get_references()

  note <- "classification by mothur2 v1.0 using default options"
  expect_equal(nrow(references), 1)
  expect_equal(references[[1, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[1, 2]], "9_032012")
  expect_equal(references[[1, 3]], "sequence_classification")
  expect_equal(references[[1, 4]], note)
  expect_equal(references[[1, "url"]], url)

  report <- dataset_t$get_sequence_taxonomy_report()

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
  dataset_t$assign_sequence_taxonomy(data = NULL, names, taxonomies)
  report <- dataset_t$get_sequence_taxonomy_report()

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
  dataset_t$assign_sequence_taxonomy(data = NULL, names, taxonomies)
  report <- dataset_t$get_sequence_taxonomy_report()

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
  dataset_t$assign_sequence_taxonomy(data = NULL, names, taxonomies)
  report <- dataset_t$get_sequence_taxonomy_report()

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
  dataset_t$assign_bins(bin_names = bins, sequence_names = names)

  report <- dataset_t$get_bin_taxonomy_report()

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

  dataset_t$clear()
  bin_ids <- c("bin1", "bin2", "bin3", "bin4")
  taxonomies <- c(
    "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;"
  )

  expect_error(dataset_t$assign_bin_taxonomy(data = NULL, bin_ids, taxonomies))
  expect_error(dataset_t$assign_bin_taxonomy(data = "not_a_data.frame"))
  expect_error(dataset_t$assign_bins(bin_ids))
  expect_equal(dataset_t$get_contigs_assembly_report(), data.frame())
  expect_equal(dataset_t$get_sample_summary(), list())
  expect_equal(dataset_t$get_sequence_summary(), list())
  expect_false(dataset_t$has_sample("noSample"))

  abunds <- c(200, 40, 100, 5)
  dataset_t$assign_bins(bin_names = bin_ids, abundances = abunds)
  dataset_t$assign_bin_taxonomy(
    data = NULL, bin_ids, taxonomies,
    reference_name = "trainset9_032012.pds.zip",
    reference_note = "classification by mothur2 v1.0",
    reference_version = "9_032012",
    reference_url = url
  )

  references <- dataset_t$get_references()

  note <- "classification by mothur2 v1.0"
  expect_equal(nrow(references), 1)
  expect_equal(references[[1, "reference_name"]], "trainset9_032012.pds.zip")
  expect_equal(references[[1, "version"]], "9_032012")
  expect_equal(references[[1, "usage"]], "bin_classification")
  expect_equal(references[[1, "note"]], note)
  expect_equal(references[[1, "url"]], url)

  report <- dataset_t$get_bin_taxonomy_report()

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


  dataset_t$assign_bin_taxonomy(data = NULL, bin_ids, taxonomies)

  report <- dataset_t$get_bin_taxonomy_report()

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

  expect_equal(dataset_t$get_sequence_taxonomy_report(), data.frame())
  expect_equal(dataset_t$get_bin_assignments(), data.frame())

  dataset_t$clear()

  expect_error(dataset_t$assign_bin_taxonomy())

  table <- readr::read_tsv(rdataset_example("final.cons.taxonomy"),
    show_col_types = FALSE
  )
  bin_table <- readr::read_tsv(rdataset_example(
    "mothur2_bin_assignments_list.tsv"
  ), show_col_types = FALSE)

  dataset_t$assign_bins(bin_table,
    bin_names = "otu_id",
    sequence_names = "seq_id"
  )

  dataset_t$assign_bin_taxonomy(table,
    bin_names = "OTU",
    taxonomies = "Taxonomy"
  )

  table <- dataset_t$get_bin_taxonomy_report()

  expect_equal(table[2758, 1], "Otu460")
  expect_equal(table[2758, 3], "\"Bacteroidales\"")
  expect_equal(table[2758, 2], 4)

  expect_equal(table[2881, 1], "Otu481")
  expect_equal(table[2881, 3], "Bacteria")
  expect_equal(table[2881, 2], 1)
})

test_that("dataset - add_metadata, get_metadata", {
  dataset_t <- dataset$new("my_dataset")

  expect_equal(dataset_t$get_metadata(), data.frame())
  expect_error(dataset_t$add_metadata(c("bad_type")))

  metadata <- readr::read_tsv(rdataset_example("sample-metadata.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )

  dataset_t$add_metadata(metadata)
  metadata <- dataset_t$get_metadata()

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

  dataset_t$clear("metadata")
  metadata <- dataset_t$get_metadata()
  expect_equal(nrow(metadata), 0)
})

test_that("dataset - add_references, get_references", {
  dataset_t <- dataset$new("my_dataset")

  expect_equal(dataset_t$get_references(), data.frame())
  expect_error(dataset_t$add_references(reference = c("bad_type")))
  expect_error(dataset_t$add_references(data.frame()))
  expect_error(dataset_t$add_references())

  references <- dataset_t$get_references()
  expect_equal(nrow(references), 0)

  reference <- readr::read_csv(rdataset_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )

  mothur_url <- "https://github.com/mothur/mothur/releases/tag/v1.48.2"

  # add data.frame and single reference at the same time
  dataset_t$add_references(
    reference = reference,
    reference_name = "mothur software package",
    usage = "analysis of dataset",
    version = "1.48.2", url = mothur_url
  )

  references <- dataset_t$get_references()

  # random spot checks
  expect_equal(nrow(references), 3)
  expect_equal(references[[1, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[2, 1]], "silva.v4.fasta")
  expect_equal(references[[3, 1]], "mothur software package")
  expect_equal(references[[1, 4]], NA_character_)
  expect_equal(references[[2, 4]], "1.38.1")
  expect_equal(references[[3, 4]], "1.48.2")


  dataset_t <- dataset$new("my_dataset")

  # add single reference then dataframe
  dataset_t$add_references(
    reference_name = "mothur software package",
    usage = "analysis of dataset",
    version = "1.48.2", url = mothur_url,
    note = "This is my mothur note"
  )

  references <- dataset_t$get_references()
  expect_equal(nrow(references), 1)

  dataset_t$add_references(reference = reference)

  references <- dataset_t$get_references()
  expect_equal(nrow(references), 3)

  expect_equal(references[[2, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[3, 1]], "silva.v4.fasta")
  expect_equal(references[[1, 1]], "mothur software package")
  expect_equal(references[[2, 2]], NA_character_)
  expect_equal(references[[3, 2]], "1.38.1")
  expect_equal(references[[1, 4]], "This is my mothur note")

  dataset_t$clear("bad_type")
  expect_equal(nrow(dataset_t$get_references()), 3)

  dataset_t$clear("references")
  expect_equal(nrow(dataset_t$get_references()), 0)
})

test_that("dataset - add_alignment_report, get_alignment_report", {
  dataset_t <- dataset$new("my_dataset")

  expect_equal(dataset_t$get_alignment_report(), data.frame())
  expect_error(dataset_t$add_alignment_report(report = c("bad_type")))
  expect_error(dataset_t$add_alignment_report(data.frame()))
  expect_error(dataset_t$add_alignment_report())

  align_report <- readr::read_tsv(rdataset_example("alignment_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )

  expect_error(dataset_t$add_alignment_report(align_report, "badName"))
  dataset_t$add_alignment_report(align_report, "QueryName")

  align_report <- dataset_t$get_alignment_report()

  # random spot checks
  expect_equal(nrow(align_report), 5)
  expect_equal(align_report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))
  expect_equal(align_report[[2, 2]], 253)
  expect_equal(align_report[[3, 3]], "AF132257.1")
  expect_equal(align_report[[1, 4]], 293)
  expect_equal(align_report[[5, 8]], 1)
  expect_equal(align_report[[4, 4]], 292)

  dataset_t$clear("alignment_report")
  expect_equal(nrow(dataset_t$get_alignment_report()), 0)

  # no report added because of missing entries

  dataset_t$add_sequences(sequence_names = c("seq6", "seq7"))
  dataset_t$add_alignment_report(align_report, "QueryName")
  expect_equal(nrow(dataset_t$get_alignment_report()), 0)
})

test_that("dataset - add / get _contigs_assembly_report,", {
  dataset_t <- dataset$new("my_dataset")

  expect_equal(dataset_t$get_contigs_assembly_report(), data.frame())
  expect_error(dataset_t$add_contigs_assembly_report(report = c("bad_type")))
  expect_error(dataset_t$add_contigs_assembly_report(data.frame()))
  expect_error(dataset_t$add_contigs_assembly_report())

  report <- readr::read_tsv(rdataset_example("contigs_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  expect_error(dataset_t$add_contigs_assembly_report(report, "badName"))

  dataset_t$add_contigs_assembly_report(report, "Name")

  report <- dataset_t$get_contigs_assembly_report()

  # random spot checks
  expect_equal(nrow(report), 5)
  expect_equal(report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))
  expect_equal(report[[2, 2]], 252)
  expect_equal(report[[3, 3]], 249)
  expect_equal(report[[1, 4]], 2)
  expect_equal(round(report[[5, 8]], digits = 4), 0.0257)
  expect_equal(report[[4, 4]], 2)

  dataset_t$clear("contigs_assembly_report")
  expect_equal(nrow(dataset_t$get_contigs_assembly_report()), 0)
  expect_equal(dataset_t$get_num_sequences(), 5)
  dataset_t$add_contigs_assembly_report(report, "Name")

  report <- dataset_t$get_contigs_assembly_report()

  expect_equal(nrow(report), 5)
  expect_equal(report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))

  # no report added because of missing entries
  dataset_t$clear("contigs_assembly_report")
  dataset_t$add_sequences(sequence_names = c("seq6", "seq7"))
  dataset_t$add_contigs_assembly_report(report, "Name")
  expect_equal(nrow(dataset_t$get_contigs_assembly_report()), 0)
})

test_that("dataset - add / get _chimera_report,", {
  dataset_t <- dataset$new("my_dataset")

  expect_equal(dataset_t$get_chimera_report(), data.frame())
  expect_error(dataset_t$add_chimera_report(report = c("bad_type")))
  expect_error(dataset_t$add_chimera_report(data.frame()))

  report <- readr::read_tsv(rdataset_example("chimera_report.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  expect_error(dataset_t$add_chimera_report(report, "badName"))

  dataset_t$add_chimera_report(report, "Query")

  report <- dataset_t$get_chimera_report()

  # random spot checks
  expect_equal(nrow(report), 71)
  expect_equal(report[, 2], dataset_t$get_sequence_names())
  expect_equal(report[[8, 5]], 82.7)
  expect_equal(report[[8, 17]], "N")
  expect_equal(report[[67, 17]], "Y")

  dataset_t$clear("chimera_report")
  expect_equal(nrow(dataset_t$get_chimera_report()), 0)
  expect_equal(dataset_t$get_num_sequences(), 71)
  dataset_t$add_chimera_report(report, "Query")

  report <- dataset_t$get_chimera_report()

  expect_equal(nrow(report), 71)
  expect_equal(report[, 2], dataset_t$get_sequence_names())

  chimera_summary <- dataset_t$get_sequence_summary()

  expect_equal(ncol(chimera_summary$chimera_summary), 13)

  # no report added because of missing entries
  dataset_t$clear("chimera_report")
  dataset_t$add_sequences(sequence_names = c("seq6", "seq7"))
  dataset_t$add_chimera_report(report, "Query")
  expect_equal(nrow(dataset_t$get_chimera_report()), 0)
})

test_that("dataset - get_sequence_summary,", {
  dataset_t <- dataset$new("my_dataset")

  report <- readr::read_tsv(rdataset_example("contigs_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  dataset_t$add_contigs_assembly_report(report, "Name")

  report <- readr::read_tsv(rdataset_example("alignment_data.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  dataset_t$add_alignment_report(report, "QueryName")

  summary <- dataset_t$get_sequence_summary()

  expect_equal(summary$contigs_summary$MisMatches, c(0, 0, 1, 2, 7, 7, 7, 2))
  expect_equal(summary$contigs_summary$Overlap_End, rep(251, 8))
  expect_equal(summary$contigs_summary$Length[1], 252)
  expect_equal(summary$contigs_summary$Length[7], 253)

  expect_equal(summary$alignment_summary$QueryLength, c(
    252, 252, 253, 253,
    253, 253, 253, 252.6
  ))
  expect_equal(summary$alignment_summary$GapsInQuery, rep(0, 8))
  expect_equal(summary$alignment_summary$SearchScore[1], 57.55)
  expect_equal(summary$alignment_summary$SearchScore[7], 82.44)
})

test_that("dataset - add_sequence_tree / get_sequence_tree,", {
  # create tree from sequences
  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ACTGC", "ATTCC", "GTTGC", "ATGGC")
  dataset_t <- dataset$new()

  expect_error(dataset_t$add_sequence_tree(tree = c("bad_type")))

  dataset_t$add_sequences(sequence_names = names, sequences = seqs)
  expect_equal(dataset_t$get_sequence_tree(), NULL)

  # add full tree
  dataset_t <- dataset$new()
  dataset_t$add_sequences(sequence_names = names)

  l <- lapply(strsplit(seqs, split = ""), "[")
  names(l) <- names

  # should alert that the tree is missing reads, and not save it
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(dataset_t$get_sequence_names()), sort(tree$tip.label))
  expect_equal(tree$edge[, 1], c(5, 5, 5, 6, 6))
  expect_equal(tree$edge[, 2], c(4, 3, 6, 1, 2))
  expect_equal(
    round(tree$edge.length, digits = 2),
    c(0.26, 0.33, 0.07, 0.33, 0.26)
  )

  # remove seq and make sure tree is pruned as well
  remove_sequences(dataset_t$data, c("seq1"), c("bad"))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(dataset_t$get_sequence_names()), sort(tree$tip.label))
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
  dataset_t$add_sequences(sequence_names = names)
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))
  expect_equal(dataset_t$get_sequence_tree(), NULL)

  # add tree from file
  dataset_t <- dataset$new()
  dataset_t$add_sequence_tree(read.tree(rdataset_example("final.phylip.tre")))
  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(dataset_t$get_sequence_names()), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(2426, 2427, 2427, 2426, 2428))
  expect_equal(tree$edge[1:5, 2], c(2427, 1, 2, 2428, 2429))
  expect_equal(
    round(tree$edge.length[1:5], digits = 3),
    c(NaN, 0.004, 0.004, 0.000, 0.002)
  )

  # add tree with extra sequences, forcing prune
  dataset_t <- dataset$new()
  dataset_t$add_sequences(sequence_names = names)

  seqs <- c(seqs, "ACTGC")
  names <- c(names, "seq5")

  # add tree with extra sequences, forcing prune
  l <- lapply(strsplit(seqs, split = ""), "[")
  names(l) <- names
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(dataset_t$get_sequence_names()), sort(tree$tip.label))
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

  remove_samples(dataset_t$data, c("F3D1", "F3D141"))

  tree <- dataset_t$get_sample_tree()

  expect_equal(sort(dataset_t$get_samples()), sort(tree$tip.label))

  # add tree with all groups, prune tree on add
  dataset_t$add_sample_tree(sample_tree)

  tree <- dataset_t$get_sample_tree()

  # confirm pruning
  expect_equal(sort(dataset_t$get_samples()), sort(tree$tip.label))
})

test_that("dataset - assign_treatments,", {
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
  dataset_t$assign_treatments(design_table)

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
  dataset_t$assign_treatments(
    data = NULL, design_table[[1]],
    design_table[[2]]
  )

  expect_equal(dataset_t$get_num_samples(), 19)
  expect_equal(dataset_t$get_num_treatments(), 2)
})

test_that("dataset - assign_sequence_taxonomy,", {
  # create dataset without taxonomy assignments
  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  tax_table <- read_mothur_taxonomy(rdataset_example("final.taxonomy"))

  # no taxonomies yet
  expect_equal(dataset_t$get_sequence_taxonomy_report(), data.frame())

  expect_error(dataset_t$assign_sequence_taxonomy(tax_table,
    sequence_names = "not_valid"
  ))
  expect_error(dataset_t$assign_sequence_taxonomy(tax_table,
    treaments = "not_valid"
  ))
  expect_error(dataset_t$assign_sequence_taxonomy())
  expect_error(dataset_t$assign_sequence_taxonomy(
    data = NULL,
    sequence_names = "not_valid"
  ))
  expect_error(dataset_t$assign_sequence_taxonomy("not_a_data.frame"))

  # test with data.frame
  dataset_t$assign_sequence_taxonomy(tax_table)

  report <- dataset_t$get_sequence_taxonomy_report()

  expect_equal(nrow(dataset_t$get_sequence_taxonomy_report()), 14550)

  dataset_t <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  expect_equal(dataset_t$get_sequence_taxonomy_report(), data.frame())

  # test with samples and treatments
  dataset_t$assign_sequence_taxonomy(
    data = NULL, tax_table[[1]],
    tax_table[[2]]
  )

  report <- dataset_t$get_sequence_taxonomy_report()

  expect_equal(nrow(dataset_t$get_sequence_taxonomy_report()), 14550)
})

test_that("dataset - export,", {
  miseq <- miseq_sop_example()

  miseq_table <- miseq$export()

  table_names <- c(
    "sequence_data", "sequence_report",
    "sequence_abundance_table", "otu_bin_data",
    "otu_sequence_bin_assignments", "asv_bin_data",
    "asv_sequence_bin_assignments", "phylotype_bin_data",
    "phylotype_sequence_bin_assignments", "metadata",
    "references", "sequence_tree", "sample_tree"
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
    "reference_name", "usage", "url", "version", "note",
    "creation_date"
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
    miseq$get_sequence_names()
  )
  expect_equal(
    miseq_table$sequence_data$sequences,
    miseq$get_sequences()
  )
})
