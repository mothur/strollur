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

    test_that("Tests addSeqs, getFastaReport, clear") {
        Dataset data("mydata", 1);

        vector<string> names(3, "");
        names[0] = "seq1";
        names[1] = "seq2";
        names[2] = "seq3";
        vector<string> seqs(3, "");
        seqs[0] = "..ATGC-MGGT-AAA-TGC-NCT.";
        seqs[1] = "ATGCGGTAAATGCCT";
        seqs[2] = ".ATGCGGGGGTAAATGCCT.";

        vector<string> comments(3, "");
        comments[2] = "my very cool comment";

        data.addSeqs(names, seqs, comments);

        expect_true(data.numUnique == 3);

        vector<vector<int> > fastaReport = data.getFastaReport();

        // starts
        vector<int> expected(3, 0);
        expected[0] = 3;
        expected[1] = 1;
        expected[2] = 2;

        expect_true(expected == fastaReport[0]);

        // ends
        expected[0] = 23;
        expected[1] = 15;
        expected[2] = 19;

        expect_true(expected == fastaReport[1]);

        // lengths
        expected[0] = 17;
        expected[1] = 15;
        expected[2] = 18;

        expect_true(expected == fastaReport[2]);

        // ambigs
        expected[0] = 2;
        expected[1] = 0;
        expected[2] = 0;

        expect_true(expected == fastaReport[3]);

        // homopolymers
        expected[0] = 3;
        expected[1] = 3;
        expected[2] = 5;

        expect_true(expected == fastaReport[4]);

        // numNs
        expected[0] = 1;
        expected[1] = 0;
        expected[2] = 0;

        expect_true(expected == fastaReport[5]);

        data.clear();

        expect_true(data.numUnique == 0);
    }
}
