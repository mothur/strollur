# test read_qiime2_feature_table

test_that("test read_qiime2_feature_table - errors", {
    expect_error(read_qiime2_feature_table("non_existant_filename"))
    expect_error(read_qiime2_feature_table(
        rdataset_example("taxonomy.qza"),
        file_root <- get_full_name("test-qiime2")
    ))
    unlink(get_full_name("test-qiime2"), recursive = TRUE)
})

test_that("test read_qiime2_feature_table", {
    bin_data <- read_qiime2_feature_table(rdataset_example("table.qza"))

    expect_equal(nrow(bin_data$data), 2290)
    expect_equal(ncol(bin_data$data), 3)

    expect_equal(bin_data$data[264, 2], 44)
    expect_equal(bin_data$data[264, 1], "b44621e5c80607cdfacbf7a81e1cbe41")
    expect_equal(bin_data$data[264, 3], "L2S357")

})
