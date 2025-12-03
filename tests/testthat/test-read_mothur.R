# test read_mothur_shared

test_that("test read_mothur", {
  # test complete dataset
  dataset <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    dataset_name = "miseq_sop"
  )

  expect_equal(names(dataset, "dataset"), "miseq_sop")
  expect_equal(num(dataset, "sequences", distinct = TRUE), 2425)
  expect_equal(num(dataset, "sequences"), 113963)
  expect_equal(num(dataset, "treatments"), 2)
  expect_equal(num(dataset, "samples"), 19)
  expect_equal(num(dataset, type = "bins", bin_type = "otu"), 531)

  # test count table only, no samples
  dataset <- read_mothur(count = rdataset_example(
    "test_nogroups.count_table"
  ))

  expect_equal(num(dataset, "sequences", distinct = TRUE), 6)
  expect_equal(num(dataset, "sequences"), 295696)
  expect_equal(num(dataset, "treatments"), 0)
  expect_equal(num(dataset, "samples"), 0)
  expect_equal(num(dataset, type = "bins"), 0)

  # test shared file and cons_tax
  dataset <- read_mothur(
    otu_shared = rdataset_example("final.opti_mcc.shared"),
    cons_taxonomy = rdataset_example(
      "final.cons.taxonomy"
    )
  )

  expect_equal(num(dataset, "sequences"), 113963)
  expect_equal(num(dataset, "samples"), 19)
  expect_equal(num(dataset, type = "bins", bin_type = "otu"), 531)

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
    asv_shared = rdataset_example("final.opti_mcc.shared"),
    phylo_shared = rdataset_example("final.opti_mcc.shared"),
    dataset_name = "miseq_sop"
  )

  expect_equal(num(dataset, type = "bins", bin_type = "asv"), 531)
  expect_equal(num(dataset, type = "bins", bin_type = "phylotype"), 531)
})

test_that("test read taxonomy files", {
  data <- read_mothur_cons_taxonomy(rdataset_example(
    "final.cons.taxonomy"
  ))

  expect_equal(nrow(data), 531)
  expect_error(read_mothur_cons_taxonomy(taxonomy = "bad_parameter"))

  data <- read_mothur_taxonomy(rdataset_example(
    "final.taxonomy"
  ))

  expect_equal(nrow(data), 2425)
  expect_error(read_mothur_taxonomy(taxonomy = "bad_parameter"))
})

test_that("test read taxonomy files", {
  dataset <- read_mothur(
    fasta = rdataset_example("final.fasta"),
    count = rdataset_example("final.count_table"),
    taxonomy = rdataset_example("final.taxonomy"),
    design = rdataset_example("mouse.time.design"),
    otu_list = rdataset_example("final.opti_mcc.list"),
    asv_list = rdataset_example("final.asv.list"),
    phylo_list = rdataset_example("final.tx.list"),
    sample_tree = rdataset_example("final.opti_mcc.jclass.ave.tre"),
    sequence_tree = rdataset_example("final.phylip.tre"),
    dataset_name = "miseq_sop"
  )

  expect_equal(names(dataset, "dataset"), "miseq_sop")
  expect_equal(num(dataset, "sequences", distinct = TRUE), 2425)
  expect_equal(num(dataset, "sequences"), 113963)
  expect_equal(num(dataset, "treatments"), 2)
  expect_equal(num(dataset, "samples"), 19)
  expect_equal(num(dataset, type = "bins", bin_type = "otu"), 531)
  expect_equal(num(dataset, type = "bins", bin_type = "asv"), 2425)
  expect_equal(num(dataset, type = "bins", bin_type = "phylotype"), 63)

  tree <- dataset$get_sample_tree()

  expect_equal(sort(names(dataset, "samples")), sort(tree$tip.label))
  expect_equal(tree$edge[1:5, 1], c(20, 21, 22, 23, 24))

  tree <- dataset$get_sequence_tree()

  expect_equal(sort(names(dataset, "sequences")), sort(tree$tip.label))
})
