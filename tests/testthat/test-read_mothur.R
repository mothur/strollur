# test read_mothur_shared

test_that("test read_mothur - errors", {
  shared <- strollur_example("final.opti_mcc.shared")
  list <- strollur_example("final.opti_mcc.list.gz")

  # can't read otulist and otushared
  expect_error(read_mothur(
    otu_shared = shared,
    otu_list = list
  ))

  list <- strollur_example("final.asv.list.gz")

  # can't read otulist and otushared
  expect_error(read_mothur(
    asv_shared = shared,
    asv_list = list
  ))

  list <- strollur_example("final.asv.list.gz")

  # can't read otulist and otushared
  expect_error(read_mothur(
    phylo_shared = shared,
    phylo_list = list
  ))
})

test_that("test read_mothur", {
  # test complete dataset
  dataset <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    taxonomy = strollur_example("final.taxonomy.gz"),
    design = strollur_example("mouse.time.design"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    dataset_name = "miseq_sop"
  )

  expect_equal(names(dataset, "dataset"), "miseq_sop")
  expect_equal(count(dataset, "sequence", distinct = TRUE), 2425)
  expect_equal(count(dataset, "sequence"), 113963)
  expect_equal(count(dataset, "treatment"), 2)
  expect_equal(count(dataset, "sample"), 19)
  expect_equal(count(dataset, type = "bin", bin_type = "otu"), 531)

  # test count table only, no samples
  dataset <- read_mothur(count = strollur_example(
    "test_nogroups.count_table"
  ))

  expect_equal(count(dataset, "sequence", distinct = TRUE), 6)
  expect_equal(count(dataset, "sequence"), 295696)
  expect_equal(count(dataset, "treatment"), 0)
  expect_equal(count(dataset, "sample"), 0)
  expect_equal(count(dataset, type = "bin"), 0)

  # test shared file and cons_tax
  dataset <- read_mothur(
    otu_shared = strollur_example("final.opti_mcc.shared"),
    cons_taxonomy = strollur_example(
      "final.cons.taxonomy"
    )
  )

  expect_equal(count(dataset, "sequence"), 113963)
  expect_equal(count(dataset, "sample"), 19)
  expect_equal(count(dataset, type = "bin", bin_type = "otu"), 531)

  report <- report(dataset, "bin_taxonomy")

  expect_equal(report[2707, 1], "Otu452")
  expect_equal(report[2707, 2], 1)
  expect_equal(report[2707, 3], "Bacteria")
  expect_equal(report[2707, 4], 100)
  expect_equal(report[2708, 1], "Otu452")
  expect_equal(report[2708, 2], 2)
  expect_equal(report[2708, 3], "\"Bacteroidetes\"")
  expect_equal(report[2708, 4], 100)

  # error
  expect_error(read_mothur(rabund = "bad_parameter"))
  expect_error(read_mothur(otu_list = "non_existant_filename"))
})

test_that("test read_mothur shared", {
  # test complete dataset
  dataset <- read_mothur(
    asv_shared = strollur_example("final.opti_mcc.shared"),
    phylo_shared = strollur_example("final.opti_mcc.shared"),
    dataset_name = "miseq_sop"
  )

  expect_equal(count(dataset, type = "bin", bin_type = "asv"), 531)
  expect_equal(count(dataset, type = "bin", bin_type = "phylotype"), 531)
})

test_that("test read taxonomy files", {
  data <- read_mothur_cons_taxonomy(strollur_example(
    "final.cons.taxonomy"
  ))

  expect_equal(nrow(data), 531)
  expect_error(read_mothur_cons_taxonomy(taxonomy = "bad_parameter"))

  data <- read_mothur_taxonomy(strollur_example(
    "final.taxonomy.gz"
  ))

  expect_equal(nrow(data), 2425)
  expect_error(read_mothur_taxonomy(taxonomy = "bad_parameter"))
})

test_that("test read taxonomy files", {
  dataset <- read_mothur(
    fasta = strollur_example("final.fasta.gz"),
    count = strollur_example("final.count_table.gz"),
    taxonomy = strollur_example("final.taxonomy.gz"),
    design = strollur_example("mouse.time.design"),
    otu_list = strollur_example("final.opti_mcc.list.gz"),
    asv_list = strollur_example("final.asv.list.gz"),
    phylo_list = strollur_example("final.tx.list.gz"),
    sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
    sequence_tree = strollur_example("final.phylip.tre.gz"),
    dataset_name = "miseq_sop"
  )

  expect_equal(names(dataset, "dataset"), "miseq_sop")
  expect_equal(count(dataset, "sequence", distinct = TRUE), 2425)
  expect_equal(count(dataset, "sequence"), 113963)
  expect_equal(count(dataset, "treatment"), 2)
  expect_equal(count(dataset, "sample"), 19)
  expect_equal(count(dataset, type = "bin", bin_type = "otu"), 531)
  expect_equal(count(dataset, type = "bin", bin_type = "asv"), 2425)
  expect_equal(count(dataset, type = "bin", bin_type = "phylotype"), 63)

  tree <- dataset$get_sample_tree()

  expect_equal(sort(names(dataset, "sample")), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(20, 21, 22, 23, 24))

  tree <- dataset$get_sequence_tree()

  expect_equal(sort(names(dataset, "sequence")), sort(tree$tip.label))
})
