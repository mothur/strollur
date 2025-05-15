#include <testthat.h>
#include "seqreport.h"

context("SeqReport class C++ unit tests") {

    test_that("Tests getReport - single seq") {
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

    test_that("Tests getReport - vector of seqs") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCGGGGGTAAATGCCT.";

        vector<string> seqs(3, "");
        seqs[0] = seq;
        seqs[1] = seq2;
        seqs[2] = seq3;

        // results[0] - starts, results[1] - ends, results[2] - numbases,
        // results[3] - ambigs, results[4] - polymers, results[5] - numns
        vector<vector<int>> results = report.getReport(seqs);

        vector<int> seqReport(3, 0);
        seqReport[0] = 3;
        seqReport[1] = 1;
        seqReport[2] = 2;

        expect_true(results[0] == seqReport);

        seqReport[0] = 23;
        seqReport[1] = 15;
        seqReport[2] = 19;

        expect_true(results[1] == seqReport);

        seqReport[0] = 17;
        seqReport[1] = 15;
        seqReport[2] = 18;

        expect_true(results[2] == seqReport);

        seqReport[0] = 2;
        seqReport[1] = 0;
        seqReport[2] = 0;

        expect_true(results[3] == seqReport);

        seqReport[0] = 3;
        seqReport[1] = 3;
        seqReport[2] = 5;

        expect_true(results[4] == seqReport);

        seqReport[0] = 1;
        seqReport[1] = 0;
        seqReport[2] = 0;

        expect_true(results[5] == seqReport);
    }

    test_that("Tests getReport - vector of seqs by reference") {
        SeqReport report;

        // includes gaps, ambiguous bases and N's
        string seq = "..ATGC-MGGT-AAA-TGC-NCT.";
        string seq2 = "ATGCGGTAAATGCCT";
        string seq3 = ".ATGCGGGGGTAAATGCCT.";

        vector<string> seqs(3, "");
        seqs[0] = seq;
        seqs[1] = seq2;
        seqs[2] = seq3;

        vector<int> starts, startsAfter;
        startsAfter.resize(3);
        startsAfter[0] = 3;
        startsAfter[1] = 1;
        startsAfter[2] = 2;

        vector<int> ends, endsAfter;
        endsAfter.resize(3);
        endsAfter[0] = 23;
        endsAfter[1] = 15;
        endsAfter[2] = 19;

        vector<int> numBases, numBasesAfter;
        numBasesAfter.resize(3);
        numBasesAfter[0] = 17;
        numBasesAfter[1] = 15;
        numBasesAfter[2] = 18;

        vector<int> ambigs, ambigsAfter;
        ambigsAfter.resize(3);
        ambigsAfter[0] = 2;
        ambigsAfter[1] = 0;
        ambigsAfter[2] = 0;

        vector<int> polymers, polymersAfter;
        polymersAfter.resize(3);
        polymersAfter[0] = 3;
        polymersAfter[1] = 3;
        polymersAfter[2] = 5;

        vector<int> numns, numnsAfter;
        numnsAfter.resize(3);
        numnsAfter[0] = 1;
        numnsAfter[1] = 0;
        numnsAfter[2] = 0;

        report.addReports(seqs, starts, ends, numBases,
                         ambigs, polymers, numns);

        expect_true(starts == startsAfter);
        expect_true(ends == endsAfter);
        expect_true(numBases == numBasesAfter);
        expect_true(ambigs == ambigsAfter);
        expect_true(polymers == polymersAfter);
        expect_true(numns == numnsAfter);
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
