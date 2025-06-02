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

    test_that("Tests addSeqs, getFastaReport, getFastaSummary, clear") {
        Dataset data("mydata", 4);

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

        Rcpp::DataFrame fastaReport = data.getFastaReport();

        // starts
        vector<int> expected(3, 0);
        expected[0] = 3;
        expected[1] = 1;
        expected[2] = 2;

        expect_true(names == Rcpp::as<vector<string>>(fastaReport[0]));

        expect_true(expected == Rcpp::as<vector<int>>(fastaReport[1]));

        // ends
        expected[0] = 23;
        expected[1] = 15;
        expected[2] = 19;

        expect_true(expected == Rcpp::as<vector<int>>(fastaReport[2]));

        // lengths
        expected[0] = 17;
        expected[1] = 15;
        expected[2] = 18;

        expect_true(expected == Rcpp::as<vector<int>>(fastaReport[3]));

        // ambigs
        expected[0] = 2;
        expected[1] = 0;
        expected[2] = 0;

        expect_true(expected == Rcpp::as<vector<int>>(fastaReport[4]));

        // homopolymers
        expected[0] = 3;
        expected[1] = 3;
        expected[2] = 5;

        expect_true(expected == Rcpp::as<vector<int>>(fastaReport[5]));

        // numNs
        expected[0] = 1;
        expected[1] = 0;
        expected[2] = 0;

        expect_true(expected == Rcpp::as<vector<int>>(fastaReport[6]));

        // report[0] = starts, report[1] = ends, report[2] = lengths,
        // report[3] = ambigs, report[4] = homopolymers, report[5] = num_ns
        Rcpp::DataFrame fastaSummary = data.getFastaSummary();

        // minimum / maximum start
        expect_true(1 == Rcpp::as<vector<int>>(fastaSummary[0])[0]);
        expect_true(3 == Rcpp::as<vector<int>>(fastaSummary[0])[6]);

        // minimum / maximum end
        expect_true(15 == Rcpp::as<vector<int>>(fastaSummary[1])[0]);
        expect_true(23 == Rcpp::as<vector<int>>(fastaSummary[1])[6]);

        // minimum / maximum lengths
        expect_true(15 == Rcpp::as<vector<int>>(fastaSummary[2])[0]);
        expect_true(18 == Rcpp::as<vector<int>>(fastaSummary[2])[6]);

        // minimum / maximum ambigs
        expect_true(0 == Rcpp::as<vector<int>>(fastaSummary[3])[0]);
        expect_true(2 == Rcpp::as<vector<int>>(fastaSummary[3])[6]);

        // minimum / maximum homopolymers
        expect_true(3 == Rcpp::as<vector<int>>(fastaSummary[4])[0]);
        expect_true(5 == Rcpp::as<vector<int>>(fastaSummary[4])[6]);

        // minimum / maximum numns
        expect_true(0 == Rcpp::as<vector<int>>(fastaSummary[5])[0]);
        expect_true(1 == Rcpp::as<vector<int>>(fastaSummary[5])[6]);

        data.clear();

        expect_true(data.numUnique == 0);
    }

    test_that("Tests addContigsReport, getContigsReport, getContigsSummary") {
        Dataset data("mydata", 4);

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

        // no contigs data yet
        Rcpp::DataFrame contigsReport = data.getContigsReport();

        expect_true(contigsReport.size() == 0);

        // add contigs data
        // oLengths, oStarts, oEnds, mismatches, expectedErrors
        vector<int> oLengths, oStarts, oEnds, mismatches;
        vector<double> ee;

        oLengths.resize(3, 5);
        oLengths[1] = 10;
        oLengths[2] = 15;

        oStarts.resize(3, 5);
        oStarts[1] = 2;
        oStarts[2] = 3;

        oEnds.resize(3, 10);
        oEnds[1] = 12;
        oEnds[2] = 18;

        mismatches.resize(3, 0);
        mismatches[1] = 2;
        mismatches[2] = 1;

        ee.resize(3, 0.0);
        ee[1] = 1.0;
        ee[2] = 2.0;

        data.addContigsReport(names, oLengths, oStarts, oEnds, mismatches, ee);

        contigsReport = data.getContigsReport();

        expect_true(contigsReport.size() == 8);

        expect_true(names == Rcpp::as<vector<string>>(contigsReport[0]));

        // lengths
        vector<int> expected(3, 0);
        expected[0] = 17;
        expected[1] = 15;
        expected[2] = 18;

        expect_true(expected == Rcpp::as<vector<int>>(contigsReport[1]));

        expect_true(oLengths == Rcpp::as<vector<int>>(contigsReport[2]));

        expect_true(oStarts == Rcpp::as<vector<int>>(contigsReport[3]));

        expect_true(oEnds == Rcpp::as<vector<int>>(contigsReport[4]));

        expect_true(mismatches == Rcpp::as<vector<int>>(contigsReport[5]));

        // numNs
        expected[0] = 1;
        expected[1] = 0;
        expected[2] = 0;

        expect_true(expected == Rcpp::as<vector<int>>(contigsReport[6]));

        expect_true(ee == Rcpp::as<vector<double>>(contigsReport[7]));

        Rcpp::DataFrame contigsSummary = data.getContigsSummary();

        // minimum / maximum ostart
        expect_true(2 == Rcpp::as<vector<int>>(contigsSummary[0])[0]);
        expect_true(5 == Rcpp::as<vector<int>>(contigsSummary[0])[6]);

        // minimum / maximum oend
        expect_true(10 == Rcpp::as<vector<int>>(contigsSummary[1])[0]);
        expect_true(18 == Rcpp::as<vector<int>>(contigsSummary[1])[6]);

        // minimum / maximum lengths
        expect_true(15 == Rcpp::as<vector<int>>(contigsSummary[2])[0]);
        expect_true(18 == Rcpp::as<vector<int>>(contigsSummary[2])[6]);

        // minimum / maximum olengths
        expect_true(5 == Rcpp::as<vector<int>>(contigsSummary[3])[0]);
        expect_true(15 == Rcpp::as<vector<int>>(contigsSummary[3])[6]);

        // minimum / maximum mismatches
        expect_true(0 == Rcpp::as<vector<int>>(contigsSummary[5])[0]);
        expect_true(2 == Rcpp::as<vector<int>>(contigsSummary[5])[6]);

        // minimum / maximum numns
        expect_true(0 == Rcpp::as<vector<int>>(contigsSummary[4])[0]);
        expect_true(1 == Rcpp::as<vector<int>>(contigsSummary[4])[6]);
    }

    test_that("Tests addAlignReport, getAlignReport, getAlignSummary") {
        Dataset data("mydata", 4);

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

        // no align data yet
        Rcpp::DataFrame alignReport = data.getAlignReport();

        expect_true(alignReport.size() == 0);

        // add align data
        // simScores, searchScores, longestInserts
        vector<int> longestInserts;
        vector<double> simScores, searchScores;

        longestInserts.resize(3, 1);
        longestInserts[1] = 2;
        longestInserts[2] = 2;

        searchScores.resize(3, 75.0);
        searchScores[1] = 86.0;
        searchScores[2] = 95.0;

        simScores.resize(3, 75.0);
        simScores[1] = 98.0;
        simScores[2] = 97.0;

        data.addAlignReport(names, searchScores, simScores, longestInserts);

        alignReport = data.getAlignReport();

        expect_true(alignReport.size() == 4);

        expect_true(names == Rcpp::as<vector<string>>(alignReport[0]));

        expect_true(searchScores == Rcpp::as<vector<double>>(alignReport[1]));

        expect_true(simScores == Rcpp::as<vector<double>>(alignReport[2]));

        expect_true(longestInserts == Rcpp::as<vector<int>>(alignReport[3]));

        Rcpp::DataFrame alignSummary = data.getAlignSummary();

        // minimum / maximum searchScores
        expect_true(75.0 == Rcpp::as<vector<double>>(alignSummary[0])[0]);
        expect_true(95.0 == Rcpp::as<vector<double>>(alignSummary[0])[6]);

        // minimum / maximum simScores
        expect_true(75.0 == Rcpp::as<vector<double>>(alignSummary[1])[0]);
        expect_true(98.0 == Rcpp::as<vector<double>>(alignSummary[1])[6]);

        // minimum / maximum longestInserts
        expect_true(1 == Rcpp::as<vector<int>>(alignSummary[2])[0]);
        expect_true(2 == Rcpp::as<vector<int>>(alignSummary[2])[6]);
    }

    test_that("Tests getNames and getSeqs") {
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

        expect_true(data.getNames() == names);
        expect_true(data.getSeqs() == seqs);

        // TODO - add tests for getNames and getSeqs with group
    }

    test_that("Tests addSeqs, assignSampleAbundance, getGroupTotals") {

        Dataset data("mydata", 1);

        vector<string> names(4, "");
        names[0] = "seq1";
        names[1] = "seq2";
        names[2] = "seq3";
        names[3] = "seq4";
        vector<string> seqs(4, "");
        seqs[0] = "..ATGC-MGGT-AAA-TGC-NCT.";
        seqs[1] = "ATGCGGTAAATGCCT";
        seqs[2] = ".ATGCGGGGGTAAATGCCT.";
        seqs[3] = ".ATGCGAMMTAAATGCCT.";

        vector<string> comments(4, "");
        comments[2] = "my very cool comment";

        data.addSeqs(names, seqs, comments);

        // mothur count file
        // Representative_Sequence     total   sample2	sample3	sample4
        // seq1	1150	250	400	500
        // seq2	115	25	40	50
        // seq3	50	25	25	0
        // seq4	5	1	0	4

        // as a sample table
        vector<string> ids(10, "");
        ids[0] = "seq1";
        ids[1] = "seq1";
        ids[2] = "seq1";
        ids[3] = "seq2";
        ids[4] = "seq2";
        ids[5] = "seq2";
        ids[6] = "seq3";
        ids[7] = "seq3";
        ids[8] = "seq4";
        ids[9] = "seq4";

        vector<string> groups(10, "");
        groups[0] = "sample2";
        groups[1] = "sample3";
        groups[2] = "sample4";
        groups[3] = "sample2";
        groups[4] = "sample3";
        groups[5] = "sample4";
        groups[6] = "sample2";
        groups[7] = "sample3";
        groups[8] = "sample2";
        groups[9] = "sample4";

        vector<int> abunds(10, 0);
        abunds[0] = 250;
        abunds[1] = 400;
        abunds[2] = 500;
        abunds[3] = 25;
        abunds[4] = 40;
        abunds[5] = 50;
        abunds[6] = 25;
        abunds[7] = 25;
        abunds[8] = 1;
        abunds[9] = 4;

        data.assignSampleAbundance(ids, abunds, groups);

        vector<int> groupTotals(3, 0);
        groupTotals[0] = 301;
        groupTotals[1] = 465;
        groupTotals[2] = 554;

        expect_true(data.getGroupTotals() == groupTotals);
        expect_true(data.getTotal() == 1320);

    }

}
