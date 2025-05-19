// Unit tests for Summary class

#include <testthat.h>
#include "summary.h"
#include "utils.h"

context("Summary class C++ unit tests") {

    test_that("test summarizeFasta") {

        vector<int> dlengths = {250, 275, 233, 24, 240, 275};
        vector<int> dstarts = {1, 2, 1, 5, 50, 3};
        vector<int> dends = {213, 285, 243, 50, 300, 279};
        vector<int> dambigs = {1, 0, 0, 5, 2, 1};
        vector<int> dpolymers = {2, 1, 2, 2, 3, 7};
        vector<int> dnumns = {0, 0, 1, 0, 0, 2};

        vector<vector<int>> report;
        report.push_back(dstarts);
        report.push_back(dends);
        report.push_back(dlengths);
        report.push_back(dambigs);
        report.push_back(dpolymers);
        report.push_back(dnumns);

        vector<int> counts(6, 0);
        counts[0] = 1;
        counts[1] = 10;
        counts[2] = 15;
        counts[3] = 1;
        counts[4] = 20;
        counts[5] = 1;

        // single thread
        Summary* summary = new Summary(1);

        Rcpp::DataFrame results = summary->summarizeFasta(report, counts);
        Rcpp::NumericVector starts = results["starts"];
        Rcpp::NumericVector ends = results["ends"];
        Rcpp::NumericVector ambigs = results["ambigs"];
        Rcpp::NumericVector nbases = results["nbases"];
        Rcpp::NumericVector polymers = results["polymers"];
        Rcpp::StringVector numseqs = results["numseqs"];

        // minimum
        expect_true(starts[0] == 1);
        // maximum
        expect_true(starts[6] == 50);

        // minimum
        expect_true(ends[0] == 50);
        // maximum
        expect_true(ends[6] == 300);

        // minimum
        expect_true(ambigs[0] == 0);
        // maximum
        expect_true(ambigs[6] == 5);

        // 25%
        expect_true(nbases[2] == 233);
        // 75%
        expect_true(nbases[4] == 250);

        // 25%
        expect_true(polymers[2] == 2);
        // 75%
        expect_true(polymers[4] == 3);

        // Median
        expect_true(numseqs[3] == 25);
        // 97.5%
        expect_true(numseqs[5] == 47);
        // maximum
        expect_true(numseqs[6] == 48);

        delete summary;

        // multi thread
        summary = new Summary(10);

        results = summary->summarizeFasta(report, counts);
        starts = results["starts"];
        ends = results["ends"];
        ambigs = results["ambigs"];
        nbases = results["nbases"];
        polymers = results["polymers"];
        numseqs = results["numseqs"];

        // minimum
        expect_true(starts[0] == 1);
        // maximum
        expect_true(starts[6] == 50);

        // minimum
        expect_true(ends[0] == 50);
        // maximum
        expect_true(ends[6] == 300);

        // 25%
        expect_true(nbases[2] == 233);
        // 75%
        expect_true(nbases[4] == 250);

        // 25%
        expect_true(polymers[2] == 2);
        // 75%
        expect_true(polymers[4] == 3);

        // Median
        expect_true(numseqs[3] == 25);
        // 97.5%
        expect_true(numseqs[5] == 47);
        // maximum
        expect_true(numseqs[6] == 48);

        delete summary;
    }

    test_that("test summarizeContigs") {
        vector<int> expectedErrors(6, 0);

        vector<int> lengths(6, 250);
        lengths[1] = 223;
        lengths[3] = 248;
        lengths[5] = 253;

        vector<int> overlapLengths(6, 250);
        overlapLengths[1] = 223;
        overlapLengths[3] = 248;
        overlapLengths[5] = 3;

        vector<int> overlapStarts(6, 0);
        overlapStarts[1] = 3;
        overlapStarts[3] = 248;
        overlapStarts[5] = 1;

        vector<int> overlapEnds(6, 250);
        overlapEnds[1] = 223;
        overlapEnds[3] = 248;
        overlapEnds[5] = 235;

        vector<int> mismatch(6, 0);
        mismatch[1] = 3;
        mismatch[3] = 1;
        mismatch[5] = 108;

        vector<int> numNs(6, 0);
        numNs[1] = 3;
        numNs[3] = 1;
        numNs[5] = 108;

        vector<int> counts(6, 0);
        counts[0] = 1;
        counts[1] = 10;
        counts[2] = 15;
        counts[3] = 1;
        counts[4] = 20;
        counts[5] = 1;

        vector<vector<int>> report;
        report.push_back(lengths);
        report.push_back(overlapLengths);
        report.push_back(overlapStarts);
        report.push_back(overlapEnds);
        report.push_back(mismatch);
        report.push_back(numNs);

        // single thread
        Summary* summary = new Summary(1);

        Rcpp::DataFrame results = summary->summarizeContigs(report, counts);
        Rcpp::NumericVector ostarts = results["ostarts"];
        Rcpp::NumericVector oends = results["oends"];
        Rcpp::NumericVector olengths = results["olengths"];
        Rcpp::NumericVector nbases = results["lengths"];
        Rcpp::NumericVector mismatches = results["mismatches"];
        Rcpp::NumericVector numns = results["numns"];
        Rcpp::StringVector numseqs = results["numseqs"];

        // minimum
        expect_true(ostarts[0] == 0);
        // maximum
        expect_true(ostarts[6] == 248);

        // minimum
        expect_true(oends[0] == 223);
        // maximum
        expect_true(oends[6] == 250);

        // minimum
        expect_true(olengths[0] == 3);
        // maximum
        expect_true(olengths[6] == 250);

        // 25%
        expect_true(nbases[2] == 250);
        // 75%
        expect_true(nbases[4] == 250);

        // 25%
        expect_true(mismatches[2] == 0);
        // 75%
        expect_true(mismatches[4] == 1);

        // 25%
        expect_true(numns[2] == 0);
        // 75%
        expect_true(numns[4] == 1);

        // Median
        expect_true(numseqs[3] == 25);
        // 97.5%
        expect_true(numseqs[5] == 47);
        // maximum
        expect_true(numseqs[6] == 48);

        delete summary;

        // multi thread
        summary = new Summary(10);
        results = summary->summarizeContigs(report, counts);

        ostarts = results["ostarts"];
        oends = results["oends"];
        olengths = results["olengths"];
        nbases = results["lengths"];
        mismatches = results["mismatches"];
        numns = results["numns"];
        numseqs = results["numseqs"];

        // minimum
        expect_true(ostarts[0] == 0);
        // maximum
        expect_true(ostarts[6] == 248);

        // minimum
        expect_true(oends[0] == 223);
        // maximum
        expect_true(oends[6] == 250);

        // minimum
        expect_true(olengths[0] == 3);
        // maximum
        expect_true(olengths[6] == 250);

        // 25%
        expect_true(nbases[2] == 250);
        // 75%
        expect_true(nbases[4] == 250);

        // 25%
        expect_true(mismatches[2] == 0);
        // 75%
        expect_true(mismatches[4] == 1);

        // 25%
        expect_true(numns[2] == 0);
        // 75%
        expect_true(numns[4] == 1);

        // Median
        expect_true(numseqs[3] == 25);
        // 97.5%
        expect_true(numseqs[5] == 47);
        // maximum
        expect_true(numseqs[6] == 48);

        delete summary;
    }

    test_that("test summarizeAlign") {

        vector<float> searchScores(6, 100.0);
        searchScores[1] = 71.05;
        searchScores[3] = 85.7;
        searchScores[5] = 93.46;

        vector<float> simScores(6, 100.0);
        simScores[1] = 90.67;
        simScores[3] = 98.78;
        simScores[5] = 92.65;

        vector<int> longestInsert(6, 0);
        longestInsert[1] = 3;
        longestInsert[3] = 248;
        longestInsert[5] = 1;

        vector<int> counts(6, 0);
        counts[0] = 1;
        counts[1] = 10;
        counts[2] = 15;
        counts[3] = 1;
        counts[4] = 20;
        counts[5] = 1;

        vector<vector<float>> report;
        report.push_back(searchScores);
        report.push_back(simScores);

        // single thread
        Summary* summary = new Summary(1);

        Rcpp::DataFrame results = summary->summarizeAlign(report,
                                                          longestInsert,
                                                          counts);

        Rcpp::NumericVector sScores = results["search_scores"];
        Rcpp::NumericVector smScores = results["sim_scores"];
        Rcpp::NumericVector lInserts = results["longest_inserts"];
        Rcpp::StringVector numseqs = results["numseqs"];

        // minimum
        int sScores0 = sScores[0]*100;
        expect_true(sScores0 == 7105);
        // maximum
        int sScores6 = sScores[6]*100;
        expect_true(sScores6 == 10000);

        // minimum
        int smScores0 = smScores[0]*100;
        expect_true(smScores0 == 9066);
        // maximum
        int smScores6 = smScores[6]*100;
        expect_true(smScores6 == 10000);

        // minimum
        int lInserts0 = lInserts[0];
        expect_true(lInserts0 == 0);
        // maximum
        int lInserts6 = lInserts[6];
        expect_true(lInserts6 == 248);

        // Median
        expect_true(numseqs[3] == 25);
        // 97.5%
        expect_true(numseqs[5] == 47);
        // maximum
        expect_true(numseqs[6] == 48);

        delete summary;

        // multi thread
        summary = new Summary(10);
        results = summary->summarizeAlign(report, longestInsert, counts);

        sScores = results["search_scores"];
        smScores = results["sim_scores"];
        lInserts = results["longest_inserts"];
        numseqs = results["numseqs"];

        // minimum
        sScores0 = sScores[0]*100;
        expect_true(sScores0 == 7105);
        // maximum
        sScores6 = sScores[6]*100;
        expect_true(sScores6 == 10000);

        // minimum
        smScores0 = smScores[0]*100;
        expect_true(smScores0 == 9066);
        // maximum
        smScores6 = smScores[6]*100;
        expect_true(smScores6 == 10000);

        // minimum
        lInserts0 = lInserts[0];
        expect_true(lInserts0 == 0);
        // maximum
        lInserts6 = lInserts[6];
        expect_true(lInserts6 == 248);

        // Median
        expect_true(numseqs[3] == 25);
        // 97.5%
        expect_true(numseqs[5] == 47);
        // maximum
        expect_true(numseqs[6] == 48);

        delete summary;
    }
}
