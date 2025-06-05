#include <testthat.h>
#include "../inst/include/rdataset.h"

context("Dataset class C++ unit tests") {

    test_that("Tests constructor") {
        Dataset data("mydata", 1);

        expect_true(data.datasetName == "mydata");
        expect_false(data.isAligned);
        expect_true(data.numSamples == 0);
        expect_true(data.numUnique == 0);
    }

    test_that("Tests addSeqs, getSequenceReport, getSequenceSummary, clear") {
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

        Rcpp::DataFrame seqReport = data.getSequenceReport();

        // starts
        vector<int> expected(3, 0);
        expected[0] = 3;
        expected[1] = 1;
        expected[2] = 2;

        expect_true(names == Rcpp::as<vector<string>>(seqReport[0]));

        expect_true(expected == Rcpp::as<vector<int>>(seqReport[1]));

        // ends
        expected[0] = 23;
        expected[1] = 15;
        expected[2] = 19;

        expect_true(expected == Rcpp::as<vector<int>>(seqReport[2]));

        // lengths
        expected[0] = 17;
        expected[1] = 15;
        expected[2] = 18;

        expect_true(expected == Rcpp::as<vector<int>>(seqReport[3]));

        // ambigs
        expected[0] = 2;
        expected[1] = 0;
        expected[2] = 0;

        expect_true(expected == Rcpp::as<vector<int>>(seqReport[4]));

        // homopolymers
        expected[0] = 3;
        expected[1] = 3;
        expected[2] = 5;

        expect_true(expected == Rcpp::as<vector<int>>(seqReport[5]));

        // numNs
        expected[0] = 1;
        expected[1] = 0;
        expected[2] = 0;

        expect_true(expected == Rcpp::as<vector<int>>(seqReport[6]));

        // report[0] = starts, report[1] = ends, report[2] = lengths,
        // report[3] = ambigs, report[4] = homopolymers, report[5] = num_ns
        Rcpp::List summary = data.getSequenceSummary();
        Rcpp::DataFrame df(summary["sequence_summary"]);

        // minimum / maximum start
        expect_true(1 == Rcpp::as<vector<int>>(df[0])[0]);
        expect_true(3 == Rcpp::as<vector<int>>(df[0])[6]);

        // minimum / maximum end
        expect_true(15 == Rcpp::as<vector<int>>(df[1])[0]);
        expect_true(23 == Rcpp::as<vector<int>>(df[1])[6]);

        // minimum / maximum lengths
        expect_true(15 == Rcpp::as<vector<int>>(df[2])[0]);
        expect_true(18 == Rcpp::as<vector<int>>(df[2])[6]);

        // minimum / maximum ambigs
        expect_true(0 == Rcpp::as<vector<int>>(df[3])[0]);
        expect_true(2 == Rcpp::as<vector<int>>(df[3])[6]);

        // minimum / maximum homopolymers
        expect_true(3 == Rcpp::as<vector<int>>(df[4])[0]);
        expect_true(5 == Rcpp::as<vector<int>>(df[4])[6]);

        // minimum / maximum numns
        expect_true(0 == Rcpp::as<vector<int>>(df[5])[0]);
        expect_true(1 == Rcpp::as<vector<int>>(df[5])[6]);

        // no treatment data
        Rcpp::DataFrame countTable = data.getSequenceAbundanceTable();

        expect_true(countTable.size() == 2);

        expect_true(names == Rcpp::as<vector<string>>(countTable[0]));
        vector<int> abunds(3, 1);
        expect_true(abunds == Rcpp::as<vector<int>>(countTable[1]));

        data.clear();

        expect_true(data.numUnique == 0);

        data.addSeqs(names, seqs, comments);

        expect_true(data.numUnique == 3);

        // as a sample table
        vector<string> ids(3, "");
        ids[0] = "seq1";
        ids[1] = "seq2";
        ids[2] = "seq3";

        abunds[0] = 250;
        abunds[1] = 400;
        abunds[2] = 500;

        data.assignSequenceAbundance(ids, abunds);
        expect_true(data.getTotal() == 1150);

        countTable = data.getSequenceAbundanceTable();

        expect_true(countTable.size() == 2);

        expect_true(ids == Rcpp::as<vector<string>>(countTable[0]));
        expect_true(abunds == Rcpp::as<vector<int>>(countTable[1]));
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

        Rcpp::List summary = data.getSequenceSummary();
        Rcpp::DataFrame df(summary["contigs_summary"]);

        // minimum / maximum ostart
        expect_true(2 == Rcpp::as<vector<int>>(df[0])[0]);
        expect_true(5 == Rcpp::as<vector<int>>(df[0])[6]);

        // minimum / maximum oend
        expect_true(10 == Rcpp::as<vector<int>>(df[1])[0]);
        expect_true(18 == Rcpp::as<vector<int>>(df[1])[6]);

        // minimum / maximum lengths
        expect_true(15 == Rcpp::as<vector<int>>(df[2])[0]);
        expect_true(18 == Rcpp::as<vector<int>>(df[2])[6]);

        // minimum / maximum olengths
        expect_true(5 == Rcpp::as<vector<int>>(df[3])[0]);
        expect_true(15 == Rcpp::as<vector<int>>(df[3])[6]);

        // minimum / maximum mismatches
        expect_true(0 == Rcpp::as<vector<int>>(df[5])[0]);
        expect_true(2 == Rcpp::as<vector<int>>(df[5])[6]);

        // minimum / maximum numns
        expect_true(0 == Rcpp::as<vector<int>>(df[4])[0]);
        expect_true(1 == Rcpp::as<vector<int>>(df[4])[6]);
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

        Rcpp::List summary = data.getSequenceSummary();
        Rcpp::DataFrame df(summary["align_summary"]);

        // minimum / maximum searchScores
        expect_true(75.0 == Rcpp::as<vector<double>>(df[0])[0]);
        expect_true(95.0 == Rcpp::as<vector<double>>(df[0])[6]);

        // minimum / maximum simScores
        expect_true(75.0 == Rcpp::as<vector<double>>(df[1])[0]);
        expect_true(98.0 == Rcpp::as<vector<double>>(df[1])[6]);

        // minimum / maximum longestInserts
        expect_true(1 == Rcpp::as<vector<int>>(df[2])[0]);
        expect_true(2 == Rcpp::as<vector<int>>(df[2])[6]);
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

        // TODO - add tests for getNames and getSeqs with samples
    }

    test_that("Tests addSeqs, assignSequenceAbundance, getSampleTotals, getTreatments, getTreatmentTotals, getSequenceAbundanceTable") {

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

        vector<string> samples(10, "");
        samples[0] = "sample2";
        samples[1] = "sample3";
        samples[2] = "sample4";
        samples[3] = "sample2";
        samples[4] = "sample3";
        samples[5] = "sample4";
        samples[6] = "sample2";
        samples[7] = "sample3";
        samples[8] = "sample2";
        samples[9] = "sample4";

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

        data.assignSequenceAbundance(ids, abunds, samples);

        vector<int> sampleTotals(3, 0);
        sampleTotals[0] = 301;
        sampleTotals[1] = 465;
        sampleTotals[2] = 554;

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getTotal() == 1320);

        // no treatment data
        Rcpp::DataFrame countTable = data.getSequenceAbundanceTable();

        expect_true(ids == Rcpp::as<vector<string>>(countTable[0]));
        expect_true(abunds == Rcpp::as<vector<int>>(countTable[1]));
        expect_true(samples == Rcpp::as<vector<string>>(countTable[2]));

        expect_true(countTable.size() == 3);

        // add with treatment assignments
        vector<string> treatments(10, "");
        treatments[0] = "early";
        treatments[1] = "early";
        treatments[2] = "late";
        treatments[3] = "early";
        treatments[4] = "early";
        treatments[5] = "late";
        treatments[6] = "early";
        treatments[7] = "early";
        treatments[8] = "early";
        treatments[9] = "late";

        data.assignSequenceAbundance(ids, abunds, samples, treatments);

        sampleTotals[0] = 301;
        sampleTotals[1] = 465;
        sampleTotals[2] = 554;

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getTotal() == 1320);

        vector<int> treatmentTotals(2, 0);
        treatmentTotals[0] = 766;
        treatmentTotals[1] = 554;

        expect_true(data.getTreatmentTotals() == treatmentTotals);

        vector<string> uniqueTreatments(2, "");
        uniqueTreatments[0] = "early";
        uniqueTreatments[1] = "late";

        expect_true(data.getTreatments() == uniqueTreatments);

        countTable = data.getSequenceAbundanceTable();

        expect_true(countTable.size() == 4);

        expect_true(ids == Rcpp::as<vector<string>>(countTable[0]));
        expect_true(abunds == Rcpp::as<vector<int>>(countTable[1]));
        expect_true(samples == Rcpp::as<vector<string>>(countTable[2]));
        expect_true(treatments == Rcpp::as<vector<string>>(countTable[3]));
    }

    test_that("Tests removeSeqs, getScrapReport, getScrapSummary") {

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
        // seq2	90	0	40	50
        // seq3	25	0	25	0
        // seq4	4	0	0	4

        // as a sample table
        vector<string> ids(7, "");
        ids[0] = "seq1";
        ids[1] = "seq1";
        ids[2] = "seq1";
        ids[3] = "seq2";
        ids[4] = "seq2";
        ids[5] = "seq3";
        ids[6] = "seq4";

        vector<string> samples(7, "");
        samples[0] = "sample2";
        samples[1] = "sample3";
        samples[2] = "sample4";
        samples[3] = "sample3";
        samples[4] = "sample4";
        samples[5] = "sample3";
        samples[6] = "sample4";

        vector<int> abunds(7, 0);
        abunds[0] = 250;
        abunds[1] = 400;
        abunds[2] = 500;
        abunds[3] = 40;
        abunds[4] = 50;
        abunds[5] = 25;
        abunds[6] = 4;

        // add with treatment assignments
        vector<string> treatments(7, "");
        treatments[0] = "early";
        treatments[1] = "early";
        treatments[2] = "early";
        treatments[3] = "late";
        treatments[4] = "late";
        treatments[5] = "late";
        treatments[6] = "late";

        data.assignSequenceAbundance(ids, abunds, samples, treatments);

        vector<int> sampleTotals(3, 0);
        sampleTotals[0] = 250;
        sampleTotals[1] = 465;
        sampleTotals[2] = 554;

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getTotal() == 1269);
        expect_true(data.numUnique == 4);
        expect_true(data.numSamples == 3);
        expect_true(data.numTreatments == 2);

        vector<int> treatmentTotals(2, 0);
        treatmentTotals[0] = 250;
        treatmentTotals[1] = 1019;

        expect_true(data.getTreatmentTotals() == treatmentTotals);

        vector<string> uniqueTreatments(2, "");
        uniqueTreatments[0] = "early";
        uniqueTreatments[1] = "late";

        expect_true(data.getTreatments() == uniqueTreatments);

        vector<string> seqsToRemove(2, "seq2");
        seqsToRemove[0] = "seq1";
        vector<string> trashCodes(2, "testRemoval");
        trashCodes[0] = "testRemoveSampleTreatment";

        // removes 2 seqs, 1 sample and 1 treatment
        data.removeSeqs(seqsToRemove, trashCodes);

        sampleTotals.resize(2, 0);
        sampleTotals[0] = 25;
        sampleTotals[1] = 4;

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getTotal() == 29);
        expect_true(data.numUnique == 2);
        expect_true(data.numSamples == 2);
        expect_true(data.numTreatments == 1);

        treatmentTotals.resize(1, 0);
        treatmentTotals[0] = 29;

        expect_true(data.getTreatmentTotals() == treatmentTotals);

        uniqueTreatments.resize(1, "");
        uniqueTreatments[0] = "late";

        expect_true(data.getTreatments() == uniqueTreatments);
        expect_true(data.getAbund("seq1") == 0);
        expect_true(data.getAbund("seq2") == 0);
        expect_true(data.getAbund("seq3") == 25);
        expect_true(data.getAbund("seq4") == 4);

        abunds.resize(2, 0);
        abunds[0] = 25;
        abunds[1] = 0;
        expect_true(data.getAbunds("seq3") == abunds);

        abunds[0] = 0;
        abunds[1] = 4;
        expect_true(data.getAbunds("seq4") == abunds);

        Rcpp::DataFrame scrapReport = data.getScrapReport();
        expect_true(scrapReport.size() == 2);

        expect_true(seqsToRemove == Rcpp::as<vector<string>>(scrapReport[0]));
        expect_true(trashCodes == Rcpp::as<vector<string>>(scrapReport[1]));

        Rcpp::DataFrame scrapSummary = data.getScrapSummary();
        expect_true(scrapSummary.size() == 3);

        sort(trashCodes.begin(), trashCodes.end());
        vector<int> uniqueCounts(2, 1);
        vector<int> totalCounts(2, 1);
        totalCounts[0] = 90;
        totalCounts[1] = 1150;

        vector<int> temp = Rcpp::as<vector<int>>(scrapSummary[2]);

        expect_true(trashCodes == Rcpp::as<vector<string>>(scrapSummary[0]));
        expect_true(uniqueCounts == Rcpp::as<vector<int>>(scrapSummary[1]));
        expect_true(totalCounts == Rcpp::as<vector<int>>(scrapSummary[2]));
    }

    test_that("Tests setAbundances") {

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

        expect_true(data.getTotal() == 4);
        expect_true(data.numUnique == 4);

        // add abundances with no samples
        vector<int> abunds(4, 0);
        abunds[0] = 250;
        abunds[1] = 400;
        abunds[2] = 500;
        abunds[3] = 10;

        data.assignSequenceAbundance(names, abunds);

        expect_true(data.getTotal() == 1160);
        expect_true(data.numUnique == 4);

        abunds[1] = 0;
        abunds[3] = 0;

        // set abundances in dataset with no samples or treatments
        data.setAbundance(names, abunds);

        expect_true(data.getTotal() == 750);
        expect_true(data.numUnique == 2);

        // mothur count file
        // Representative_Sequence     total   sample2	sample3	sample4
        // seq1	1150	250	400	500
        // seq2	90	0	40	50
        // seq3	25	0	25	0
        // seq4	4	0	0	4

        // as a sample table
        vector<string> ids(7, "");
        ids[0] = "seq1";
        ids[1] = "seq1";
        ids[2] = "seq1";
        ids[3] = "seq2";
        ids[4] = "seq2";
        ids[5] = "seq3";
        ids[6] = "seq4";

        vector<string> samples(7, "");
        samples[0] = "sample2";
        samples[1] = "sample3";
        samples[2] = "sample4";
        samples[3] = "sample3";
        samples[4] = "sample4";
        samples[5] = "sample3";
        samples[6] = "sample4";

        abunds.resize(7, 0);
        abunds[0] = 250;
        abunds[1] = 400;
        abunds[2] = 500;
        abunds[3] = 40;
        abunds[4] = 50;
        abunds[5] = 25;
        abunds[6] = 4;

        // add with treatment assignments
        vector<string> treatments(7, "");
        treatments[0] = "early";
        treatments[1] = "late";
        treatments[2] = "late";
        treatments[3] = "late";
        treatments[4] = "late";
        treatments[5] = "late";
        treatments[6] = "late";

        data.clear();

        data.addSeqs(names, seqs, comments);

        expect_true(data.getTotal() == 4);
        expect_true(data.numUnique == 4);

        data.assignSequenceAbundance(ids, abunds, samples, treatments);

        expect_true(data.getTotal() == 1269);
        expect_true(data.numUnique == 4);

        vector<vector<int>> parsedAbunds(4, vector<int>(3, 0));
        // removes sample2
        parsedAbunds[0][1] = 40; // sample3, late
        parsedAbunds[0][2] = 400; // sample4, late
        parsedAbunds[1][2] = 35; // sample4, late
        parsedAbunds[2][1] = 100; // sample3, late
        parsedAbunds[2][2] = 1; // sample4, late
        parsedAbunds[3][1] = 20; // sample3, late
        parsedAbunds[3][2] = 500; // sample4, late

        // set abundances in dataset with samples and treatments
        data.setAbundances(names, parsedAbunds);

        expect_true(data.getTotal() == 1096);
        expect_true(data.numUnique == 4);

        vector<int> sampleTotals(2, 0);
        sampleTotals[0] = 160;
        sampleTotals[1] = 936;

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.numSamples == 2);
        expect_true(data.numTreatments == 1);

        vector<int> treatmentTotals(1, 0);
        treatmentTotals[0] = 1096;

        expect_true(data.getTreatmentTotals() == treatmentTotals);

        vector<string> uniqueTreatments(1, "late");

        expect_true(data.getTreatments() == uniqueTreatments);
    }


}
