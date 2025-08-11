# tests clone of sequence_data object

test_that("clone - deep copy of sequence_data object", {
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

    dataset <- clone(temp)

    expect_equal(dataset$get_dataset_name(), "miseq_sop")
    expect_equal(dataset$get_num_sequences(TRUE), 2425)
    expect_equal(dataset$get_num_sequences(), 113963)
    expect_equal(dataset$get_num_treatments(), 2)
    expect_equal(dataset$get_num_samples(), 19)
    expect_equal(dataset$get_num_bins("otu"), 531)
    expect_equal(dataset$get_num_bins("phylotype"), 63)
    expect_equal(dataset$get_num_bins("asv"), 2425)

})
