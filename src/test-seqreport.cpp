#include <testthat.h>
#include "seqreport.h"

context("SeqReport class C++ unit tests") {

    test_that("Tests getReport") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCGGGGGTAAATGCCT.";

        vector<int> seqReport(6, 0);
        seqReport[0] = 3;
        seqReport[1] = 23;
        seqReport[2] = 17;
        seqReport[3] = 2;
        seqReport[4] = 3;
        seqReport[5] = 1;

        expect_true(report.getReport(seq) == seqReport);

        seqReport[0] = 1;
        seqReport[1] = 15;
        seqReport[2] = 15;
        seqReport[3] = 0;
        seqReport[4] = 3;
        seqReport[5] = 0;

        expect_true(report.getReport(seq2) == seqReport);

        seqReport[0] = 2;
        seqReport[1] = 19;
        seqReport[2] = 18;
        seqReport[3] = 0;
        seqReport[4] = 5;
        seqReport[5] = 0;

        expect_true(report.getReport(seq3) == seqReport);

    }

    test_that("Tests getStart") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCGGGGGTAAATGCCT.";

        expect_true(report.getStart(seq) == 3);
        expect_true(report.getStart(seq2) == 1);
        expect_true(report.getStart(seq3) == 2);
    }

    test_that("Tests getEnd") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCGGGGGTAAATGCCT.";

        expect_true(report.getEnd(seq) == 23);
        expect_true(report.getEnd(seq2) == 15);
        expect_true(report.getEnd(seq3) == 19);
    }

    test_that("Tests getNumbases") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCGGGGGTAAATGCCT.";

        expect_true(report.getNumbases(seq) == 17);
        expect_true(report.getNumbases(seq2) == 15);
        expect_true(report.getNumbases(seq3) == 18);
    }

    test_that("Tests getNumAmbigs") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCGGGGGTAAATGCCT.";

        expect_true(report.getNumAmbigs(seq) == 2);
        expect_true(report.getNumAmbigs(seq2) == 0);
        expect_true(report.getNumAmbigs(seq3) == 0);
    }

    test_that("Tests getLongestHomopolymer") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCAAATTGCCTGGGGG.";

        expect_true(report.getLongestHomopolymer(seq) == 3);
        expect_true(report.getLongestHomopolymer(seq2) == 3);
        expect_true(report.getLongestHomopolymer(seq3) == 5);
    }

    test_that("Tests getNumns") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCGGGGGTAAATGCCT.";

        expect_true(report.getNumns(seq) == 1);
        expect_true(report.getNumns(seq2) == 0);
        expect_true(report.getNumns(seq3) == 0);
    }


}
