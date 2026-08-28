test_that("write_quality - errors", {
    expect_error(write_quality("Bad_type"))

    # no file name with nameless dataset
    data <- new_dataset()
    expect_error(write_quality(data))
})

test_that("write_quality", {
    data <- strollur::new_dataset(dataset_name = "example_dataset")

    fastq_data <-
        strollur::read_fastq(fastq = strollur_example("tiny.fastq.gz"))

    strollur::add(data = data, table = fastq_data, type = "fastq")

    output <- write_quality(data, get_full_name("example_dataset"))

    df <- read_quality(output)

    remove_file(output)

    expect_equal(df[[1]], names(data, "sequence"))

    score_1_15 <- c(2, 29, 29, 32, 32, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)
    score_2_15 <- c(18, 32, 32, 30, 32, 33, 33, 35, 33, 37, 37, 33, 36, 38, 38)
    score_3_15 <- c(33, 32, 31, 33, 33, 33, 32, 33, 33, 37, 37, 37, 38, 38, 38)

    expect_equal(df$quality_score[[1]][1:15], score_1_15)
    expect_equal(df$quality_score[[2]][1:15], score_2_15)
    expect_equal(df$quality_score[[3]][1:15], score_3_15)

    data <- read_mothur(
        otu_list = strollur_example("final.opti_mcc.list.gz"),
        dataset_name = "data"
    )

    expect_equal(write_fastq(data), "no_fastq_data")
})
