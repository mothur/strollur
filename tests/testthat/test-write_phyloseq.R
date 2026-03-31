test_that("write phyloseq creates phyloseq object", {
  list_file <- strollur_example("final.asv.list.gz")
  shared_file <- strollur_example("final.opti_mcc.shared")
  constaxonomy_file <- strollur_example("final.cons.taxonomy")

  phyloseq_object <- phyloseq::import_mothur(
    mothur_list_file = list_file,
    mothur_shared_file = shared_file,
    mothur_constaxonomy_file =
      constaxonomy_file
  )

  data <- read_phyloseq(phyloseq_object)

  recreated_phylo_object <- write_phyloseq(data)
  expect_true(inherits(recreated_phylo_object, "phyloseq"))
  expect_true(all(phyloseq::tax_table(phyloseq_object) ==
                    phyloseq::tax_table(recreated_phylo_object)))
  expect_true(all(phyloseq::otu_table(phyloseq_object) ==
                    phyloseq::otu_table(recreated_phylo_object)))
  expect_true(all(phyloseq::sample_data(data, FALSE) ==
                    phyloseq::sample_data(recreated_phylo_object, FALSE)))

  dat <- readRDS(strollur_example("GlobalPatterns.RDS"))
  phyloseq_object <- read_phyloseq(dat)
  recreated_phylo_object <- write_phyloseq(phyloseq_object)
  expect_true(all(phyloseq::tax_table(dat) ==
                    phyloseq::tax_table(recreated_phylo_object), na.rm = TRUE))
  expect_true(all(phyloseq::otu_table(dat) ==
                    phyloseq::otu_table(recreated_phylo_object)))
  expect_true(all(phyloseq::sample_data(dat) ==
                    phyloseq::sample_data(recreated_phylo_object)))
  expect_true(length(phyloseq::phy_tree(dat)$edge) ==
                length(phyloseq::phy_tree(recreated_phylo_object)$edge))


  miseq <- miseq_sop_example()
  recreated_phylo_object <- write_phyloseq(miseq)
  data <- read_phyloseq(recreated_phylo_object,
    treatment_column_name = "treatments"
  )
  expect_identical(
    report(miseq, "sample_assignments"),
    report(data, "sample_assignments")
  )
})


test_that("write phyloseq fails if not given a strollur object", {
  miseq <- c()
  expect_error(
    write_phyloseq(miseq),
    "The data parameter must be an object of type `strollur`."
  )
})


test_that("write phyloseq will fail if the dataset is empty", {
  empty_dataset <- strollur$new("")
  expect_error(write_phyloseq(empty_dataset))
})

test_that("Will error if phyloseq is not installed", {
  local_mocked_bindings(require_namespace = function(...) FALSE)
  expect_error(
    write_phyloseq(c()),
    "To use this functionality"
  )
})
