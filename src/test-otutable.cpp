#include <testthat.h>
#include "../inst/include/rdataset.h"

context("OtuTable class C++ unit tests") {

    test_that("Tests constructor") {
        OtuTable otuTable("0.03");

        expect_true(otuTable.label == "0.03");
        expect_true(otuTable.numOtus == 0);
        expect_true(otuTable.numUnique == -1);
    }

    test_that("Tests add") {
        OtuTable otuTable("0.03");

        vector<string> otuNames(10, "otu1");
        otuNames[1] = "otu2";
        otuNames[2] = "otu3";
        otuNames[3] = "otu4";
        otuNames[4] = "otu5";
        otuNames[5] = "otu6";
        otuNames[6] = "otu7";
        otuNames[7] = "otu8";
        otuNames[8] = "otu9";
        otuNames[9] = "otu10";
        vector<int> abundances(10, 10);

        // test adding otuNames and abundances (rabund)
        otuTable.add(otuNames, abundances);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.numOtus == 10);
        expect_true(otuTable.numUnique == -1);

        otuTable.clear();

        // test adding otuNames, seqNames, abundances (list)

        vector<string> seqNames(10, "");
        otuNames[0] = "otu1";   seqNames[0] = "seq1";
        otuNames[1] = "otu1";   seqNames[1] = "seq2";
        otuNames[2] = "otu1";   seqNames[2] = "seq3";
        otuNames[3] = "otu2";   seqNames[3] = "seq4";
        otuNames[4] = "otu2";   seqNames[4] = "seq5";
        otuNames[5] = "otu3";   seqNames[5] = "seq6";
        otuNames[6] = "otu4";   seqNames[6] = "seq7";
        otuNames[7] = "otu4";   seqNames[7] = "seq8";
        otuNames[8] = "otu4";   seqNames[8] = "seq9";
        otuNames[9] = "otu4";   seqNames[9] = "seq10";

        otuTable.add(otuNames, abundances, nullVector, seqNames);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.numOtus == 4);
        expect_true(otuTable.numUnique == 10);

        // test adding otuNames abundances, samples (shared)
        otuNames.resize(15, "otu1");
        abundances.resize(15, 10);
        vector<string> samples(15, "sample1");
        otuNames[0] = "otu1";   samples[0] = "sample1";  abundances[0] = 10;
        otuNames[1] = "otu1";   samples[1] = "sample2";  abundances[1] = 10;
        otuNames[2] = "otu1";   samples[2] = "sample4";  abundances[2] = 5;
        otuNames[3] = "otu1";   samples[3] = "sample5";  abundances[3] = 5;
        otuNames[4] = "otu2";   samples[4] = "sample1";  abundances[4] = 5;
        otuNames[5] = "otu2";   samples[5] = "sample2";  abundances[5] = 5;
        otuNames[6] = "otu2";   samples[6] = "sample4";  abundances[6] = 10;
        otuNames[7] = "otu3";   samples[7] = "sample1";  abundances[7] = 1;
        otuNames[8] = "otu3";   samples[8] = "sample3";  abundances[8] = 2;
        otuNames[9] = "otu3";   samples[9] = "sample5";  abundances[9] = 3;
        otuNames[10] = "otu3";   samples[10] = "sample6";  abundances[10] = 4;
        otuNames[11] = "otu4";   samples[11] = "sample1";  abundances[11] = 20;
        otuNames[12] = "otu4";   samples[12] = "sample2";  abundances[12] = 10;
        otuNames[13] = "otu4";   samples[13] = "sample4";  abundances[13] = 5;
        otuNames[14] = "otu4";   samples[14] = "sample5";  abundances[14] = 5;

        // add shared data
        otuTable.add(otuNames, abundances, samples);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.numOtus == 4);
        expect_true(otuTable.numUnique == 10);
        expect_true(otuTable.getNumSamples() == 6);

        vector<int> sampleTotals(6, 0);
        sampleTotals[0] = 36;
        sampleTotals[1] = 25;
        sampleTotals[2] = 2;
        sampleTotals[3] = 20;
        sampleTotals[4] = 13;
        sampleTotals[5] = 4;

        expect_true(otuTable.getSampleTotals() == sampleTotals);
    }
}
