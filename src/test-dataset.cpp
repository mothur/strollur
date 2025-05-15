#include <testthat.h>
#include "../inst/include/rdataset.h"

context("Dataset class C++ unit tests") {

    test_that("Tests constructor") {
        Dataset data("mydata", 1);

        expect_true(data.datasetName == "mydata");
        expect_false(data.isAligned);
        expect_true(data.numGroups == 0);
        expect_true(data.numUnique == 0);
    }

//     test_that("Tests clear") {
//         Dataset data("mydata", 1);
//
//         vector<string> names(3, "");
//         names[1] = "seq1";
//         names[2] = "seq2";
//         names[3] = "seq3";
//         vector<string> seqs(3, "");
//         seqs[1] = "ATTGGCTG";
//         seqs[2] = "ATCAGCTG";
//         seqs[3] = "TGTGGCTG";
//
//         data.addSeqs(names, seqs);
//
//         expect_true(data.numUnique == 3);
//
//         data.clear();
//
//         expect_true(data.numUnique == 0);
//     }
}
