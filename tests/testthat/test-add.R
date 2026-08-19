# tests add function

test_that("test add - piping", {
  data <- new_dataset()

  table <- data.frame(sequence_name = c("seq1", "seq2", "seq3"))

  abunds <- add(data, table) |> abundance()

  expect_equal(abunds$sequence_name, c("seq1", "seq2", "seq3"))
  expect_equal(abunds$abundance, rep(1, 3))
})

test_that("test add - errors", {
  data <- new_dataset()

  # not a strollur object
  x <- 10
  expect_error(add(x))

  fasta_data <- read_fasta(strollur_example("final.fasta.gz"))

  # test bad type
  expect_error(add(data, table = fasta_data, type = bad_type))

  # test bad report_type
  expect_error(
    add(data, table = fasta_data, type = report),
    "'report_type' is required when adding a report."
  )

  add(data, table = fasta_data)
  expect_equal(count(data, type = "sequence"), 2425)

  reference <- readr::read_csv(strollur_example("references.csv"),
    col_names = TRUE, show_col_types = FALSE
  )

  add(data, reference, "resource_reference")

  references <- report(data, "resource_reference")

  # random spot checks
  expect_equal(nrow(references), 2)
  expect_equal(references[[1, "name"]], "R phylotypr package")
  expect_equal(references[[2, "name"]], "silva.bacteria.fasta")
  expect_equal(references[[1, "note"]], "classification using Bayesian method")
  expect_equal(
    references[[2, "note"]],
    "alignment reference trimmed to V4 region"
  )
})

test_that("test add - sequence_tree, sample_tree", {
  names <- c("seq1", "seq2", "seq3", "seq4")
  seqs <- c("ACTGC", "ATTCC", "GTTGC", "ATGGC")

  dataset_t <- new_dataset()
  xdev_add_sequences(dataset_t, data.frame(sequence_name = names))

  l <- lapply(strsplit(seqs, split = ""), "[")
  names(l) <- names

  add(dataset_t, table = nj(dist.dna(as.DNAbin(l))), type = "sequence_tree")

  tree <- dataset_t$get_sequence_tree()

  expect_equal(sort(names(dataset_t, "sequence")), sort(tree$tip.label))
  expect_equal(tree$edge[, 1], c(5, 5, 5, 6, 6))
  expect_equal(tree$edge[, 2], c(4, 3, 6, 1, 2))
  expect_equal(
    round(tree$edge.length, digits = 2),
    c(0.26, 0.33, 0.07, 0.33, 0.26)
  )

  dataset_t <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    dataset_name = "miseq_sop"
  )

  sample_tree <- ape::read.tree(
    strollur_example("final.opti_mcc.jclass.ave.tre")
  )

  expect_error(add(dataset_t, table = c("bad_type"), type = "sample_tree"))
  add(dataset_t, table = tree, type = "sample_tree")
  expect_null(dataset_t$get_sample_tree())

  add(dataset_t, table = sample_tree, type = "sample_tree")

  tree <- dataset_t$get_sample_tree()

  expect_equal(sort(names(dataset_t, "sample")), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(20, 21, 22, 23, 24))

  xdev_remove_samples(dataset_t, c("F3D1", "F3D141"))

  tree <- dataset_t$get_sample_tree()

  expect_equal(sort(names(dataset_t, "sample")), sort(tree$tip.label))

  # add tree with all groups, prune tree on add
  add(dataset_t, table = sample_tree, type = "sample_tree")

  tree <- dataset_t$get_sample_tree()

  # confirm pruning
  expect_equal(sort(names(dataset_t, "sample")), sort(tree$tip.label))
})
