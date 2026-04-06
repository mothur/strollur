# test "dataset"

test_that("dataset - intialize from read_mothur / print", {
  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    taxonomy = strollur_example("final.taxonomy.gz"),
    design = strollur_example("mouse.time.design"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    phylo_list = strollur_example("final.tx.list.gz"),
    asv_list = strollur_example("final.asv.list.gz"),
    dataset_name = "miseq_sop"
  )

  # add references, custom report and metadata
  contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
  add(
    data = dataset_t, table = contigs_report, type = "reports",
    report_type = "contigs_report", list(sequence_name = "Name")
  )

  metadata <- readRDS(strollur_example("miseq_metadata.rds"))

  add(data = dataset_t, table = metadata, type = "metadata")

  reference <- readr::read_csv(strollur_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )
  xdev_add_references(dataset_t, reference)

  actual <- dataset_t$report(type = "metadata")
  expect_equal(nrow(actual), 19)

  expect_equal(dataset_t$get_bin_types(), c("otu", "asv", "phylotype"))
  expect_equal(names(dataset_t, "dataset")[1], "miseq_sop")
  expect_equal(
    count(data = dataset_t, type = "sequences", distinct = TRUE),
    2425
  )
  expect_equal(count(data = dataset_t, type = "sequences"), 113963)
  expect_equal(count(data = dataset_t, type = "treatments"), 2)
  expect_equal(count(data = dataset_t, type = "samples"), 19)
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "otu"), 531)
  expect_equal(
    count(data = dataset_t, type = "bins", bin_type = "phylotype"),
    63
  )
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "asv"), 2425)

  seqs_summary <- summary(dataset_t, type = "sequences")

  expect_equal(seqs_summary$starts[1], 1)
  expect_equal(seqs_summary$ends[1], 375)
  expect_equal(seqs_summary$ambigs[2], 0)
  expect_equal(seqs_summary$numns[4], 0)
  expect_equal(seqs_summary$nbases[1], 249)
  expect_equal(seqs_summary$nbases[7], 256)
  expect_equal(seqs_summary$polymers[1], 3)
  expect_equal(seqs_summary$polymers[7], 6)

  seq_report <- report(dataset_t, type = "sequences")

  expect_equal(nrow(seq_report), 2425)
  expect_equal(ncol(seq_report), 7)

  # remove bin from "phylotype" list and confirm that it removes seqs from all
  # from all list types
  # "Phylo05"
  df <- abundance(dataset_t, type = "bins", bin_type = "phylotype")
  expect_equal(df[[2]][5], 5337)
  # "Phylo06"
  expect_equal(df[[2]][6], 715)

  phylo05 <- xdev_get_list_vector(dataset_t, "phylotype")[5]
  expect_equal(length(.split_at_char(phylo05)), 54)

  phylo06 <- xdev_get_list_vector(dataset_t, "phylotype")[6]
  expect_equal(length(.split_at_char(phylo06)), 47)

  xdev_remove_bins(
    dataset_t,
    c("Phylo05", "Phylo06"),
    c("test", "test"),
    "phylotype"
  )

  expect_equal(
    count(data = dataset_t, type = "bins", bin_type = "phylotype"),
    61
  )
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "otu"), 512)
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "asv"), 2324)
  expect_equal(count(data = dataset_t, type = "sequences"), 107911)
  expect_equal(
    count(data = dataset_t, type = "sequences", distinct = TRUE),
    2324
  )
  expect_equal(count(data = dataset_t, type = "treatments"), 2)
  expect_equal(count(data = dataset_t, type = "samples"), 19)
  expect_equal(count(
    data = dataset_t,
    type = "sequences",
    samples = c("F3D0")
  ), 5977)
  expect_equal(count(
    data = dataset_t,
    type = "sequences",
    samples = c("F3D1")
  ), 4467)

  # note that the number of seqs removed will be less that 297+266 because
  # some seqs are assigned to both samples and some seqs will be present in
  # other samples
  expect_equal(count(
    data = dataset_t, type = "sequences",
    distinct = TRUE, samples = "F3D0"
  ), 99)
  expect_equal(count(
    data = dataset_t, type = "sequences",
    distinct = TRUE, samples = "F3D1"
  ), 99)

  # remove samples
  xdev_remove_samples(dataset_t, c("F3D0", "F3D1"))
  expect_equal(count(data = dataset_t, type = "treatments"), 2)
  expect_equal(count(data = dataset_t, type = "samples"), 17)
  expect_equal(count(
    data = dataset_t,
    type = "bins",
    bin_type = "phylotype"
  ), 57)
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "otu"), 482)
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "asv"), 2124)
  expect_equal(count(data = dataset_t, type = "sequences"), 97467)
  expect_equal(count(
    data = dataset_t,
    type = "sequences",
    distinct = TRUE
  ), 2124)
  expect_equal(count(
    data = dataset_t, type = "sequences", distinct = TRUE,
    sample = "F3D0"
  ), 0)
  expect_equal(count(
    data = dataset_t, type = "sequences", distinct = TRUE,
    sample = "F3D1"
  ), 0)
  expect_equal(count(
    data = dataset_t, type = "sequences",
    sample = "F3D0"
  ), 0)
  expect_equal(count(
    data = dataset_t, type = "sequences", distinct = TRUE,
    sample = "F3D1"
  ), 0)

  # remove things just classified to bacteria, and things classified to
  # Bacteria;"Bacteroidetes"; with confidence less than 95
  xdev_remove_lineages(dataset_t, c(
    "Bacteria;Bacteria_unclassified;",
    "Bacteria(100);\"Bacteroidetes\"(95);"
  ))

  expect_equal(count(data = dataset_t, type = "treatments"), 2)
  expect_equal(count(data = dataset_t, type = "samples"), 17)
  expect_equal(count(
    data = dataset_t,
    type = "bins",
    bin_type = "phylotype"
  ), 57)
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "otu"), 475)
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "asv"), 2086)
  expect_equal(count(data = dataset_t, type = "sequences"), 97428)
  expect_equal(count(
    data = dataset_t,
    type = "sequences",
    distinct = TRUE
  ), 2086)

  expect_snapshot(
    waldo::compare(dataset_t$print(), dataset_t$print())
  )

  dataset_t$clear()
  expect_equal(count(data = dataset_t, type = "sequences"), 0)
})

test_that("dataset - intialize from dataset object", {
  temp <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    taxonomy = strollur_example("final.taxonomy.gz"),
    design = strollur_example("mouse.time.design"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    phylo_list = strollur_example("final.tx.list.gz"),
    asv_list = strollur_example("final.asv.list.gz"),
    dataset_name = "miseq_sop"
  )

  dataset_t <- strollur$new(
    name = "clone_of_miseq", dataset = temp,
    processors = 4
  )

  expect_equal(names(dataset_t, "dataset"), "clone_of_miseq")
  expect_equal(count(
    data = dataset_t,
    type = "sequences",
    distinct = TRUE
  ), 2425)
  expect_equal(count(data = dataset_t, type = "sequences"), 113963)
  expect_equal(count(data = dataset_t, type = "treatments"), 2)
  expect_equal(count(data = dataset_t, type = "samples"), 19)
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "otu"), 531)
  expect_equal(count(
    data = dataset_t,
    type = "bins",
    bin_type = "phylotype"
  ), 63)
  expect_equal(count(data = dataset_t, type = "bins", bin_type = "asv"), 2425)
})

test_that("dataset - addSeqs, assign samples", {
  data <- strollur$new("mydata")

  # missing data and names
  expect_error(xdev_add_sequences(data))

  fasta_data <- read_fasta(strollur_example("final.fasta.gz"))
  names(fasta_data) <- c("myNameTag", "mySeqTag")

  expect_error(xdev_add_sequences(data, fasta_data, NULL, "names", "mySeqTag"))
  expect_error(xdev_add_sequences(data, "not_a_data.frame"))

  xdev_add_sequences(data, fasta_data, NULL, "myNameTag", "mySeqTag")

  expect_equal(count(data = data, type = "sequences"), 2425)

  clear(data)

  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ATTGC", "ATTGC", "ATTGC", "ATTGC")
  comments <- c("ddd", "ftf", "efr", "ssd")

  fasta_data <- data.frame(names = names, seqs = seqs, comments = comments)

  expect_error(xdev_add_sequences(data, fasta_data,
    sequence_names = "names",
    sequences = "sequences"
  ))
  xdev_add_sequences(data, fasta_data, NULL, "names", "seqs")

  expect_equal(count(data = data, type = "sequences"), 4)
  expect_equal(xdev_get_sequences(data), seqs)
  clear(data)

  expect_error(xdev_add_sequences(data, fasta_data,
    sequences = "seqs",
    comments = "comments23"
  ))
  xdev_add_sequences(data, fasta_data, NULL, "names", "seqs", "comments")

  expect_equal(count(data = data, type = "sequences"), 4)
  expect_equal(xdev_get_sequences(data), seqs)

  clear(data)
  xdev_add_sequences(data, fasta_data, NULL, "names", "", "comments")

  expect_equal(count(data = data, type = "sequences"), 4)
  expect_equal(xdev_get_sequences(data), rep("", 4))
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
  xdev_add_sequences(
    data,
    data.frame(sequence_names = names, sequences = seqs),
    new_reference(
      "silva.bacteria.fasta", "1.38.1", "",
      "alignment by mothur2 v1.0", url
    )
  )

  references <- report(data, "references")

  expect_equal(nrow(references), 1)
  expect_equal(references[[1, "reference_names"]], "silva.bacteria.fasta")
  expect_equal(references[[1, "reference_versions"]], "1.38.1")
  expect_equal(references[[1, "reference_usages"]], "NA")
  expect_equal(references[[1, "reference_notes"]], "alignment by mothur2 v1.0")
  expect_equal(references[[1, "reference_urls"]], url)

  assign(
    data = data, table = data.frame(
      sequence_names = ids, abundances = abundances,
      samples = samples, treatments = treatments
    ), type = "sequence_abundance"
  )

  # assign bins
  bins <- c("bin1", "bin2", "bin1", "bin2")
  assign(
    data = data,
    table = data.frame(bin_names = bins, sequence_names = names),
    type = "bins"
  )

  list <- report(data = data, type = "sequence_bin_assignments")

  expect_equal(count(data = data, type = "bins"), 2)
  expect_equal(list$otu_id, c("bin1", "bin1", "bin2", "bin2"))
  expect_equal(list$seq_id, c("seq1", "seq3", "seq2", "seq4"))

  list <- report(
    data = data, type = "sequence_bin_assignments",
    bin_type = "asv"
  )

  # get asv generated by dataset
  expect_equal(list$asv_id, c("ASV1", "ASV2", "ASV3", "ASV4"))
  expect_equal(list$seq_id, c("seq1", "seq2", "seq3", "seq4"))

  expect_equal(abundance(data, type = "bins")[[1]], c("bin1", "bin2"))
  expect_equal(abundance(data, type = "bins")[[2]], c(1200, 120))

  df <- abundance(
    data = data, type = "bins",
    bin_type = "otu", by_sample = TRUE
  )
  expect_equal(df$bin_names, c(
    "bin1", "bin1", "bin1",
    "bin2", "bin2", "bin2"
  ))
  expect_equal(df$samples, c(
    "sample2", "sample3", "sample4",
    "sample2", "sample3", "sample4"
  ))
  expect_equal(df$abundances, c(
    275, 425, 500,
    26, 40, 54
  ))

  expect_true(is_aligned(data))
  expect_equal(names(data, "sequences"), names)
  expect_equal(xdev_get_sequences(data), seqs)
  expect_equal(
    names(data = data, type = "sequences", samples = c("sample2")),
    names
  )
  expect_equal(
    names(data, type = "sequences", samples = c("sample3")),
    c("seq1", "seq2", "seq3")
  )
  expect_equal(
    names(data, type = "sequences", samples = c("sample4")),
    c("seq1", "seq2", "seq4")
  )
  expect_equal(xdev_get_sequences(data, "sample2"), seqs)
  expect_equal(
    xdev_get_sequences(data, "sample3"),
    c("ATTGC", "ATTGC", "ATTGC")
  )
  expect_equal(
    xdev_get_sequences(data, "sample4"),
    c("ATTGC", "ATTGC", "ATTGC")
  )

  expect_equal(count(data = data, type = "samples"), 3)
  expect_equal(count(data = data, type = "treatments"), 2)
  expect_equal(names(data, "treatments"), c("early", "late"))
  expect_equal(names(data, "samples"), c("sample2", "sample3", "sample4"))

  sample_summary <- abundance(data, type = "samples")
  treatment_summary <- abundance(data, type = "treatments")

  expect_equal(treatment_summary$abundances, c(766, 554))
  expect_equal(sample_summary$abundances, c(301, 465, 554))

  # total
  expect_equal(count(data = data, type = "sequences"), 1320)
  # unique
  expect_equal(count(data = data, type = "sequences", distinct = TRUE), 4)
  # unique and sample
  expect_equal(count(
    data = data, type = "sequences",
    distinct = TRUE, samples = c("sample2")
  ), 0)
  expect_equal(count(
    data = data, type = "sequences",
    distinct = FALSE, samples = c("sample2")
  ), 301)
  expect_equal(count(
    data = data, type = "sequences",
    distinct = FALSE, samples = c("sample3")
  ), 465)

  expect_equal(count(
    data = data, type = "sequences",
    distinct = TRUE, samples = c("sample3")
  ), 0)
})

test_that("dataset - assign_sequence_abundance, remove_sequences", {
  data <- strollur$new("mydata")

  # missing data and names
  expect_error(xdev_assign_sequence_abundance(data))

  sequence_abundance <- readRDS(
    strollur_example("miseq_abundance_by_sample.rds")
  )

  assign(
    data = data, table = sequence_abundance, type = "sequence_abundance"
  )

  expect_equal(count(data = data, type = "sequences", distinct = TRUE), 2425)
  expect_equal(count(data = data, type = "samples"), 19)
  expect_equal(count(data = data, type = "treatments"), 2)

  names(sequence_abundance) <- c("ids", "abunds", "groups", "time")

  expect_error(xdev_assign_sequence_abundance(data, sequence_abundance, "ids"))

  # no treatments
  assign(
    data = data, table = sequence_abundance, type = "sequence_abundance",
    table_names = list(
      sequence_name = "ids",
      abundance = "abunds",
      sample = "groups"
    )
  )

  expect_equal(count(data = data, type = "sequences", distinct = TRUE), 2425)
  expect_equal(count(data = data, type = "samples"), 19)
  expect_equal(count(data = data, type = "treatments"), 0)

  assign(
    data = data, table = sequence_abundance, type = "sequence_abundance",
    table_names = list(
      sequence_name = "ids",
      abundance = "abunds",
      sample = "groups",
      treatment = "time"
    )
  )

  expect_equal(count(data = data, type = "sequences", distinct = TRUE), 2425)
  expect_equal(count(data = data, type = "samples"), 19)
  expect_equal(count(data = data, type = "treatments"), 2)

  clear(data)

  names <- c("seq1", "seq2", "seq3", "seq4")
  abunds <- c(10, 20, 30)

  expect_error(xdev_assign_sequence_abundance(
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

  xdev_assign_sequence_abundance(data, data.frame(
    sequence_names = names,
    abundances = rabunds
  ))

  expect_equal(count(data = data, type = "sequences"), 1269)
  expect_equal(count(data = data, type = "sequences", distinct = TRUE), 4)
  expect_equal(count(data = data, type = "samples"), 0)
  expect_equal(count(data = data, type = "treatments"), 0)

  missing_id <- c(
    "seq1", "seq1", "seq1",
    "seq2", "seq2",
    "seq3",
    "seq3"
  )

  xdev_assign_sequence_abundance(data, data.frame(
    sequence_names = ids,
    abundances = abundances,
    samples = groups
  ))

  expect_equal(count(data = data, type = "sequences"), 1269)
  expect_equal(count(data = data, type = "sequences", distinct = TRUE), 4)
  expect_equal(count(data = data, type = "samples"), 3)
  expect_equal(count(data = data, type = "treatments"), 0)

  xdev_assign_sequence_abundance(
    data, data.frame(
      sequence_names = ids, abundances = abundances,
      samples = groups, treatments = treatments
    )
  )
  expect_equal(count(data = data, type = "treatments"), 2)

  xdev_remove_sequences(data, seqs_to_remove, trash_codes)

  expect_equal(count(data = data, type = "sequences"), 29)
  expect_equal(count(data = data, type = "sequences", distinct = TRUE), 2)
  expect_equal(count(data = data, type = "samples"), 2)
  expect_equal(count(data = data, type = "treatments"), 1)
})

test_that("dataset - get_list get_rabund, get_bin_assignments", {
  dataset_t <- new_dataset("my_dataset")

  expect_error(assign(dataset_t))

  seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  sequence_abundances <- c(10, 100, 1, 500, 25, 80)

  assign(data = dataset_t, table = data.frame(
    bin_names = bin_ids,
    abundances = sequence_abundances,
    sequence_names = seq_ids
  ), type = "bins")
  # bins would look like:
  # label  bin1             bin2        bin3 ...
  # 0.03   seq1,seq2,seq4   seq3,seq6   seq5 ...
  # 0.03   110              525         80 ...

  list <- report(data = dataset_t, type = "sequence_bin_assignments")

  expect_equal(list$otu_id, bin_ids)
  expect_equal(list$seq_id, seq_ids)
  expect_equal(report(dataset_t, "non_existance_bin_type"), data.frame())
  expect_equal(
    xdev_get_list_vector(dataset_t, "non_existance_bin_type"),
    character()
  )

  rabund <- dataset_t$abundance(type = "bins", bin_type = "otu")

  abunds <- c(111, 525, 80)
  expect_equal(rabund$otu_id, unique(bin_ids))
  expect_equal(rabund$abundance, abunds)

  expect_equal(
    abundance(
      data = dataset_t, type = "bins",
      bin_type = "non_existance_bin_type"
    ),
    data.frame()
  )
  expect_equal(
    abundance(
      data = dataset_t, type = "bins",
      bin_type = "otu"
    )[[2]],
    abunds
  )

  dataset_t <- strollur$new("my_dataset")
  bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
  samples <- c(
    "sample1", "sample2", "sample5",
    "sample1", "sample3", "sample1"
  )
  sample_abundances <- c(10, 100, 1, 500, 25, 80)
  xdev_assign_bins(
    dataset_t,
    data.frame(
      bin_names = bin_ids, abundances = sample_abundances,
      samples = samples
    )
  )

  shared <- abundance(data = dataset_t, type = "bins", by_sample = TRUE)
  expect_equal(shared$bin_names, bin_ids)
  expect_equal(shared$abundances, sample_abundances)
  expect_equal(shared$samples, samples)
  expect_equal(count(data = dataset_t, type = "bins"), 3)

  expect_equal(
    abundance(
      data = dataset_t,
      type = "bins", bin_type = "non_existance_bin_type"
    ),
    data.frame()
  )

  expect_equal(shared[[2]], c(10, 100, 1, 500, 25, 80))
  expect_true(has_sample(dataset_t, "sample1"))
  expect_false(has_sample(dataset_t, "non_existant_sample"))

  dataset_t <- strollur$new("my_dataset")
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

  xdev_assign_bins(
    dataset_t,
    data.frame(
      bin_names = bin_ids, abundances = sample_abundances,
      samples = samples, sequence_names = seq_ids
    )
  )

  expect_equal(count(dataset_t, "bins"), 3)
  expect_equal(count(dataset_t, "samples"), 6)
  expect_equal(
    abundance(dataset_t, type = "samples")[["abundances"]],
    c(36, 25, 2, 20, 13, 4)
  )

  clear(dataset_t)

  bin_table <- readRDS(strollur_example("miseq_shared_otu.rds"))

  expect_error(assign(
    data = dataset_t, table = bin_table,
    type = "bins", bin_type = "otu",
    table_names = list(bin_name = "id")
  ))

  xdev_assign_bins(dataset_t, bin_table)

  expect_equal(count(dataset_t, "bins"), 531)
  expect_equal(count(dataset_t, "samples"), 19)

  clear(dataset_t)

  bin_table <- readRDS(strollur_example("miseq_list_otu.rds"))

  assign(
    data = dataset_t, table = bin_table, type = "bins",
    bin_type = "otu"
  )

  expect_equal(count(dataset_t, "bins"), 531)
  expect_equal(count(dataset_t, "sequences"), 2425)
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

  dataset_t <- strollur$new("my_dataset")

  url <- paste0(
    "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset",
    "9_032012.pds.zip"
  )

  # assign taxonomy with reference
  xdev_assign_sequence_taxonomy(
    dataset_t,
    data.frame(sequence_names = names, taxonomies = taxonomies),
    new_reference(
      "trainset9_032012.pds.zip", "9_032012",
      "classification by mothur2 v1.0 using default options", "",
      url
    )
  )

  references <- report(dataset_t, "references")

  note <- "classification by mothur2 v1.0 using default options"
  expect_equal(nrow(references), 1)
  expect_equal(references[[1, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[1, 2]], "9_032012")
  expect_equal(references[[1, 3]], note)
  expect_equal(references[[1, 4]], "NA")
  expect_equal(references[[1, 5]], url)

  report <- report(dataset_t, "sequence_taxonomy")

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
  #---------------------------------------------------------#
  dataset_t <- strollur$new("my_dataset")
  assign(data = dataset_t, table = data.frame(
    sequence_names = names,
    taxonomies = taxonomies
  ), type = "sequence_taxonomy")
  report <- report(dataset_t, "sequence_taxonomy")

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
  #---------------------------------------------------------#
  dataset_t <- strollur$new("my_dataset")
  xdev_assign_sequence_taxonomy(dataset_t, data.frame(
    sequence_names = names,
    taxonomies = taxonomies
  ))
  report <- report(dataset_t, "sequence_taxonomy")

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
  #---------------------------------------------------------#
  dataset_t <- strollur$new("my_dataset")
  xdev_assign_sequence_taxonomy(dataset_t, data.frame(
    sequence_names = names,
    taxonomies = taxonomies
  ))
  report <- report(dataset_t, "sequence_taxonomy")

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
  xdev_assign_bins(
    dataset_t,
    data.frame(bin_names = bins, sequence_names = names)
  )

  report <- report(dataset_t, "bin_taxonomy")

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
  #---------------------------------------------------------#
  clear(dataset_t)
  bin_ids <- c("bin1", "bin2", "bin3", "bin4")
  taxonomies <- c(
    "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
    "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
    "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
    "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;"
  )

  expect_error(assign(dataset_t, table = "not_a_data.frame"))
  expect_false(has_sample(dataset_t, "noSample"))

  abunds <- c(200, 40, 100, 5)
  xdev_assign_bins(
    dataset_t,
    data.frame(bin_names = bin_ids, abundances = abunds)
  )

  url <- paste0(
    "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset",
    "9_032012.pds.zip"
  )

  assign(
    data = dataset_t,
    table = data.frame(bin_names = bin_ids, taxonomies = taxonomies),
    type = "bin_taxonomy", bin_type = "otu",
    reference = new_reference(
      reference_name = "trainset9_032012.pds.zip",
      reference_version = "9_032012",
      reference_usage = "classification by mothur2 v1.0",
      reference_url = url
    )
  )

  references <- report(dataset_t, "references")

  expect_equal(nrow(references), 1)
  expect_equal(references[[1, "reference_names"]], "trainset9_032012.pds.zip")
  expect_equal(references[[1, "reference_versions"]], "9_032012")
  expect_equal(references[[1, "reference_notes"]], "NA")
  expect_equal(
    references[[1, "reference_usages"]],
    "classification by mothur2 v1.0"
  )
  expect_equal(references[[1, "reference_urls"]], url)

  report <- report(dataset_t, "bin_taxonomy")

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


  assign(data = dataset_t, table = data.frame(
    bin_names = bin_ids,
    taxonomies = taxonomies
  ), type = "bin_taxonomy")

  report <- report(dataset_t, "bin_taxonomy")

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

  expect_equal(report(dataset_t, "sequence_taxonomy"), data.frame())
  #---------------------------------------------------------#
  clear(dataset_t)

  table <- readr::read_tsv(strollur_example("final.cons.taxonomy"),
    show_col_types = FALSE
  )
  bin_table <- readRDS(strollur_example("miseq_list_otu.rds"))

  assign(
    data = dataset_t, table = bin_table, type = "bins",
    bin_type = "otu"
  )

  assign(
    data = dataset_t, table = table,
    type = "bin_taxonomy",
    bin_type = "otu", list(
      bin_name = "OTU",
      taxonomy = "Taxonomy"
    )
  )

  table <- report(dataset_t, "bin_taxonomy")

  expect_equal(table[2758, 1], "Otu460")
  expect_equal(table[2758, 3], "\"Bacteroidales\"")
  expect_equal(table[2758, 2], 4)

  expect_equal(table[2881, 1], "Otu481")
  expect_equal(table[2881, 3], "Bacteria")
  expect_equal(table[2881, 2], 1)
})

test_that("dataset - add_metadata, get_metadata", {
  dataset_t <- strollur$new("my_dataset")

  expect_equal(report(dataset_t, "metadata"), data.frame())

  metadata <- readr::read_tsv(strollur_example("sample-metadata.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )

  dataset_t$add(metadata, "metadata")
  metadata <- report(dataset_t, "metadata")

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

test_that("dataset - add_references, get_references", {
  dataset_t <- strollur$new("my_dataset")

  expect_equal(report(dataset_t, "references"), data.frame())
  expect_error(xdev_add_references(dataset_t, reference = c("bad_type")))
  expect_error(xdev_add_references(dataset_t, data.frame()))
  expect_error(xdev_add_references(dataset_t))

  references <- report(dataset_t, "references")
  expect_equal(nrow(references), 0)

  reference <- readr::read_csv(strollur_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )

  mothur_url <- "https://github.com/mothur/mothur/releases/tag/v1.48.2"

  # add data.frame and single reference at the same time
  xdev_add_references(dataset_t, reference)

  references <- report(dataset_t, "references")

  # random spot checks
  expect_equal(nrow(references), 2)
  expect_equal(references[[1, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[2, 1]], "silva.v4.fasta")
  expect_equal(references[[1, 4]], "NA")
  expect_equal(
    references[[2, 4]],
    "custom reference created by trimming silva.bacteria.fasta to the V4 region"
  )

  dataset_t <- strollur$new("my_dataset")

  # add single reference then dataframe
  ref <- data.frame(
    reference_name = "mothur software package",
    reference_version = "1.48.2",
    reference_usage = "analysis of dataset",
    reference_note = "This is my mothur note",
    reference_url = mothur_url
  )

  xdev_add_references(
    dataset_t, ref, "reference_name", "reference_version",
    "reference_usage", "reference_note", "reference_url"
  )

  references <- report(dataset_t, "references")
  expect_equal(nrow(references), 1)

  xdev_add_references(dataset_t, reference)

  references <- report(dataset_t, "references")
  expect_equal(nrow(references), 3)

  expect_equal(references[[2, 1]], "trainset9_032012.pds.zip")
  expect_equal(references[[3, 1]], "silva.v4.fasta")
  expect_equal(references[[1, 1]], "mothur software package")
  expect_equal(references[[2, 2]], "NA")
  expect_equal(references[[3, 2]], "1.38.1")
  expect_equal(references[[1, 4]], "This is my mothur note")
})

test_that("dataset - add_alignment_report, get_alignment_report", {
  dataset_t <- strollur$new("my_dataset")

  align_report <- readRDS(strollur_example("test_alignment_data.rds"))

  expect_error(xdev_add_report(
    dataset_t,
    align_report,
    "align_report",
    "badName"
  ))
  xdev_add_report(dataset_t, align_report, "align_report", "QueryName")

  align_report <- report(dataset_t, "align_report")

  # random spot checks
  expect_equal(nrow(align_report), 5)
  expect_equal(align_report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))
  expect_equal(align_report[[2, 2]], 253)
  expect_equal(align_report[[3, 3]], "AF132257.1")
  expect_equal(align_report[[1, 4]], 293)
  expect_equal(align_report[[5, 8]], 1)
  expect_equal(align_report[[4, 4]], 292)

  clear(dataset_t)
  expect_equal(report(dataset_t, "align_report"), data.frame())

  # no report added because of missing entries
  xdev_add_sequences(dataset_t, data.frame(sequence_names = c("seq6", "seq7")))
  xdev_add_report(dataset_t, align_report, "align_report", "QueryName")
  expect_equal(report(dataset_t, "align_report"), data.frame())
})

test_that("dataset - add / get _contigs_assembly_report,", {
  dataset_t <- strollur$new("my_dataset")

  report <- readRDS(strollur_example("test_contigs_data.rds"))
  expect_error(xdev_add_report(dataset_t, report, "contigs_report", "badName"))

  xdev_add_report(dataset_t, report, "contigs_report", "Name")

  report <- report(dataset_t, "contigs_report")

  # random spot checks
  expect_equal(nrow(report), 5)
  expect_equal(report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))
  expect_equal(report[[2, 2]], 252)
  expect_equal(report[[3, 3]], 249)
  expect_equal(report[[1, 4]], 2)
  expect_equal(round(report[[5, 8]], digits = 4), 0.0257)
  expect_equal(report[[4, 4]], 2)

  clear(dataset_t)
  xdev_add_report(dataset_t, report, "contigs_report", "Name")

  report <- report(dataset_t, "contigs_report")

  expect_equal(nrow(report), 5)
  expect_equal(report[, 1], c("seq1", "seq2", "seq3", "seq4", "seq5"))

  # no report added because of missing entries
  clear(dataset_t)
  xdev_add_sequences(dataset_t, data.frame(sequence_names = c("seq6", "seq7")))
  xdev_add_report(dataset_t, report, "contigs_report", "Name")
  expect_equal(length(report(dataset_t, "contigs_report")), 0)
})

test_that("dataset - add / get _chimera_report,", {
  dataset_t <- strollur$new("my_dataset")

  report <- readr::read_tsv(strollur_example("chimera_report.tsv"),
    col_names = TRUE, show_col_types = FALSE
  )
  expect_error(xdev_add_report(dataset_t, report, "chimera_report", "badName"))

  xdev_add_report(dataset_t, report, "chimera_report", "Query")

  report <- report(dataset_t, "chimera_report")

  # random spot checks
  expect_equal(nrow(report), 71)
  expect_equal(report[, 2], names(dataset_t, "sequences"))
  expect_equal(report[[8, 5]], 82.7)
  expect_equal(report[[8, 17]], "N")
  expect_equal(report[[67, 17]], "Y")

  clear(dataset_t)
  expect_equal(length(report(dataset_t, "chimera_report")), 0)
  xdev_add_report(dataset_t, report, "chimera_report", "Query")

  report <- report(dataset_t, "chimera_report")

  expect_equal(nrow(report), 71)
  expect_equal(report[, 2], names(dataset_t, "sequences"))
})

test_that("dataset - get_sequence_summary,", {
  dataset_t <- strollur$new("my_dataset")

  report <- readRDS(strollur_example("test_contigs_data.rds"))
  xdev_add_report(dataset_t, report, "contigs_report", "Name")

  report <- readRDS(strollur_example("test_alignment_data.rds"))
  xdev_add_report(dataset_t, report, "alignment_report", "QueryName")

  summary <- summary(dataset_t,
    type = "reports",
    report_type = "contigs_report"
  )

  expect_equal(summary$MisMatches, c(0, 0, 1, 2, 7, 7, 7, 2))
  expect_equal(summary$Overlap_End, rep(251, 8))
  expect_equal(summary$Length[1], 252)
  expect_equal(summary$Length[7], 253)

  summary <- summary(dataset_t,
    type = "reports",
    report_type = "alignment_report"
  )

  expect_equal(summary$QueryLength, c(
    252, 252, 253, 253,
    253, 253, 253, 252.6
  ))
  expect_equal(summary$GapsInQuery, rep(0, 8))
  expect_equal(summary$SearchScore[1], 57.55)
  expect_equal(summary$SearchScore[7], 82.44)
})

test_that("dataset - add_sequence_tree / get_sequence_tree,", {
  # create tree from sequences

  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ACTGC", "ATTCC", "GTTGC", "ATGGC")
  dataset_t <- strollur$new()

  expect_error(dataset_t$add_sequence_tree(tree = c("bad_type")))

  xdev_add_sequences(
    dataset_t,
    data.frame(sequence_names = names, sequences = seqs)
  )
  expect_equal(dataset_t$get_sequence_tree(), NULL)

  # add full tree
  dataset_t <- strollur$new()
  xdev_add_sequences(dataset_t, data.frame(sequence_names = names))

  l <- lapply(strsplit(seqs, split = ""), "[")
  names(l) <- names

  # should alert that the tree is missing reads, and not save it
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(names(dataset_t, "sequences")), sort(tree$tip.label))
  expect_equal(tree$edge[, 1], c(5, 5, 5, 6, 6))
  expect_equal(tree$edge[, 2], c(4, 3, 6, 1, 2))
  expect_equal(
    round(tree$edge.length, digits = 2),
    c(0.26, 0.33, 0.07, 0.33, 0.26)
  )

  # remove seq and make sure tree is pruned as well
  xdev_remove_sequences(dataset_t, c("seq1"), c("bad"))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(names(dataset_t, "sequences")), sort(tree$tip.label))
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
  dataset_t <- strollur$new()
  xdev_add_sequences(dataset_t, data.frame(sequence_names = names))
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))
  expect_equal(dataset_t$get_sequence_tree(), NULL)

  # add tree from file
  dataset_t <- strollur$new()
  dataset_t$add_sequence_tree(read.tree(
    strollur_example("final.phylip.tre.gz")
  ))
  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(names(dataset_t, "sequences")), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(2426, 2427, 2427, 2426, 2428))
  expect_equal(tree$edge[1:5, 2], c(2427, 1, 2, 2428, 2429))
  expect_equal(
    round(tree$edge.length[1:5], digits = 3),
    c(NaN, 0.004, 0.004, 0.000, 0.002)
  )

  # add tree with extra sequences, forcing prune
  dataset_t <- strollur$new()
  xdev_add_sequences(dataset_t, data.frame(sequence_names = names))

  seqs <- c(seqs, "ACTGC")
  names <- c(names, "seq5")

  # add tree with extra sequences, forcing prune
  l <- lapply(strsplit(seqs, split = ""), "[")
  names(l) <- names
  dataset_t$add_sequence_tree(nj(dist.dna(as.DNAbin(l))))

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(names(dataset_t, "sequences")), sort(tree$tip.label))
  expect_equal(tree$edge[, 1], c(5, 5, 5, 6, 6))
  expect_equal(tree$edge[, 2], c(4, 3, 6, 1, 2))
  expect_equal(
    round(tree$edge.length, digits = 2),
    c(0.26, 0.33, 0.07, 0.33, 0.26)
  )
})

test_that("dataset - add_sample_tree / get_sample_tree,", {
  dataset_t <- strollur$new()
  expect_null(dataset_t$get_sample_tree())

  sample_tree <- ape::read.tree(
    strollur_example("final.opti_mcc.jclass.ave.tre")
  )

  # should report no samples and not save
  dataset_t$add_sample_tree(sample_tree)
  expect_null(dataset_t$get_sample_tree())

  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    taxonomy = strollur_example("final.taxonomy.gz"),
    design = strollur_example("mouse.time.design"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    dataset_name = "miseq_sop"
  )

  expect_error(dataset_t$add_sample_tree(tree = c("bad_type")))

  sequence_tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))

  # should report missing samples since this is a sequence tree and not save
  dataset_t$add_sample_tree(sequence_tree)
  expect_null(dataset_t$get_sample_tree())

  dataset_t$add_sample_tree(sample_tree)

  tree <- dataset_t$get_sample_tree()

  expect_equal(sort(names(dataset_t, "samples")), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(20, 21, 22, 23, 24))

  xdev_remove_samples(dataset_t, c("F3D1", "F3D141"))

  tree <- dataset_t$get_sample_tree()

  expect_equal(sort(names(dataset_t, "samples")), sort(tree$tip.label))

  # add tree with all groups, prune tree on add
  dataset_t$add_sample_tree(sample_tree)

  tree <- dataset_t$get_sample_tree()

  # confirm pruning
  expect_equal(sort(names(dataset_t, "samples")), sort(tree$tip.label))
})

test_that("dataset - assign_treatments", {
  # create dataset without treatment assignments
  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    phylo_list = strollur_example("final.tx.list.gz"),
    asv_list = strollur_example("final.asv.list.gz"),
    dataset_name = "miseq_sop"
  )

  design_table <- readr::read_tsv(
    strollur_example(
      "mouse.time.design"
    ),
    show_col_types = FALSE
  )

  expect_equal(count(dataset_t, "samples"), 19)
  expect_equal(count(dataset_t, "treatments"), 0)

  expect_error(assign(
    data = dataset_t, table = design_table,
    type = "treatments",
    table_names = list(sample = "group")
  ))
  expect_error(xdev_assign_treatments(dataset_t))
  expect_error(xdev_assign_treatments(data = dataset_t, samples = "not_valid"))
  expect_error(xdev_assign_treatments(dataset_t, "not_a_data.frame"))

  # test with data.frame
  assign(data = dataset_t, table = design_table, type = "treatments")

  expect_equal(count(dataset_t, "samples"), 19)
  expect_equal(count(dataset_t, "treatments"), 2)

  report <- dataset_t$report(type = "sample_assignments")

  expect_equal(report$treatments, c(
    "Early", "Early", "Late", "Late", "Late",
    "Late", "Late", "Late", "Late", "Late",
    "Late", "Late", "Early", "Early", "Early",
    "Early", "Early", "Early", "Early"
  ))

  expect_equal(dataset_t$count(type = "samples"), 19)
  expect_equal(dataset_t$names(type = "treatments"), c("Early", "Late"))
  expect_equal(dataset_t$summary(type = "scrap"), data.frame())

  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    phylo_list = strollur_example("final.tx.list.gz"),
    asv_list = strollur_example("final.asv.list.gz"),
    dataset_name = "miseq_sop"
  )

  expect_equal(count(dataset_t, "samples"), 19)
  expect_equal(count(dataset_t, "treatments"), 0)

  # test with samples and treatments
  assign(data = dataset_t, table = design_table, type = "treatments")

  expect_equal(count(dataset_t, "samples"), 19)
  expect_equal(count(dataset_t, "treatments"), 2)
})

test_that("dataset - assign_sequence_taxonomy", {
  # create dataset without taxonomy assignments
  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    dataset_name = "miseq_sop"
  )

  tax_table <- read_mothur_taxonomy(strollur_example("final.taxonomy.gz"))

  # no taxonomies yet
  expect_equal(report(dataset_t, "sequence_taxonomy"), data.frame())

  expect_error(xdev_assign_sequence_taxonomy(
    data = dataset_t, table = tax_table,
    sequence_name = "not_valid"
  ))
  expect_error(xdev_assign_sequence_taxonomy(
    data = dataset_t, table = tax_table,
    sequence_name = "sequence_names", taxonomy = "not_valid"
  ))
  expect_error(xdev_assign_sequence_taxonomy(dataset_t))

  # test with data.frame
  dataset_t$assign(table = tax_table, type = "sequence_taxonomy")

  report <- report(dataset_t, "sequence_taxonomy")

  expect_equal(nrow(report(dataset_t, "sequence_taxonomy")), 14550)

  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    dataset_name = "miseq_sop"
  )

  expect_equal(report(dataset_t, "sequence_taxonomy"), data.frame())

  # test with samples and treatments
  xdev_assign_sequence_taxonomy(dataset_t, tax_table)

  report <- report(dataset_t, "sequence_taxonomy")

  expect_equal(nrow(report), 14550)
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
    "contigs_report", "sequence_tree", "sample_tree"
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
  expect_equal(
    names(miseq_table$sequence_abundance_table),
    sequence_at_names
  )

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
    names(miseq, "sequences")
  )
  expect_equal(
    miseq_table$sequence_data$sequences,
    xdev_get_sequences(miseq)
  )
})

test_that("dataset - assign_bin_representative_sequences", {
  # create dataset
  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    dataset_name = "miseq_sop"
  )

  # select first 531 seqs to be the representatives
  num_bins <- count(dataset_t, "bins")
  rep_names <- names(dataset_t, "sequences")[1:num_bins]
  bin_names <- names(dataset_t, "bins")

  assign(
    data = dataset_t,
    table = data.frame(
      bin_names = bin_names,
      sequence_names = rep_names
    ),
    type = "bin_representatives",
    bin_type = "otu"
  )

  df <- report(dataset_t, "bin_representatives")

  expect_equal(df[[1]], bin_names)
  expect_equal(df[[2]], rep_names)
  expect_equal(df[[3]], xdev_get_sequences(dataset_t)[1:num_bins])

  d <- strollur$new()
  expect_equal(report(d, "bin_representatives"), data.frame())
  expect_equal(report(d, "sample_assignments"), data.frame())

  dataset_t <- read_mothur(
    count = strollur_example("final.count_table.gz"),
    dataset_name = "miseq_sop"
  )
  expect_equal(report(dataset_t, "sample_assignments"), data.frame())
})
