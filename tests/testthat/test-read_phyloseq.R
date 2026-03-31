test_that("test read_phyloseq", {
  list_file <- strollur_example("esophagus.fn.list")
  groups_file <- strollur_example("esophagus.good.groups")
  tree_file <- strollur_example("esophagus.tree")

  phyloseq_object <- phyloseq::import_mothur(
    list_file,
    groups_file,
    tree_file
  )
  data <- read_phyloseq(phyloseq_object)
  expect_true(count(data) == 675)

  list_file <- strollur_example("final.asv.list.gz")
  shared_file <- strollur_example("final.opti_mcc.shared")
  constaxonomy_file <- strollur_example("final.cons.taxonomy")

  phyloseq_object_w_tax <- phyloseq::import_mothur(
    mothur_list_file = list_file,
    mothur_shared_file =
      shared_file,
    mothur_constaxonomy_file =
      constaxonomy_file
  )

  data <- read_phyloseq(phyloseq_object_w_tax)
  row_data <- phyloseq::tax_table(phyloseq_object_w_tax)[9, ]
  rdata_set_data <- report(data, type = "sequence_taxonomy")
  expect_true(abundance(data)[1, ]$abundances == 12288)
  expect_true(abundance(data)[2, ]$abundances == 8892)
  expect_true(abundance(data)[3, ]$abundances == 7794)

  expect_true(nrow(rdata_set_data) == 3186)
  expect_true(rdata_set_data[49, 3] == row_data[[1]])
  expect_true(rdata_set_data[50, 3] == row_data[[2]])
  expect_true(rdata_set_data[51, 3] == row_data[[3]])
  expect_true(rdata_set_data[52, 3] == row_data[[4]])
  expect_true(rdata_set_data[53, 3] == row_data[[5]])
  expect_true(rdata_set_data[54, 3] == row_data[[6]])
})

test_that("read phyloseq fails if not given an phyloseq object", {
  empty <- c()
  expect_error(
    read_phyloseq(empty),
    paste0(
      "phyloseq_object has to an object created using ",
      "the phyloseq package."
    )
  )
})

test_that("Will error if phyloseq is not installed", {
  local_mocked_bindings(require_namespace = function(...) FALSE)
  expect_error(
    read_phyloseq(c()),
    "To use this functionality"
  )
})
