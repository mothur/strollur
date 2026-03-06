#include <testthat.h>
#include "../inst/include/strollur.h"
#include "dataset.h"

 context("Dataset class C++ unit tests") {

    test_that("Tests dataset.h utils") {
        vector<float> vector_with_dups(10, 1);

        expect_true(toSet(vector_with_dups).size() == 1);

        vector<bool> notAllTrue(2, true);

        expect_true(isTrue(notAllTrue));

        notAllTrue[0] = false;

        expect_false(isTrue(notAllTrue));
        expect_false(isFalse(notAllTrue));

        notAllTrue[1] = false;
        expect_true(isFalse(notAllTrue));

        int x;
        expect_error(convert("notANumber", x));

        vector<int> test1 = {1,2,3,4,5};
        vector<int> test2 = {3,2,1,5,4};

        expect_true(identical(test1, test2));

        test2[2] = 6;

        expect_false(identical(test1, test2));
    }

    test_that("Tests constructor") {
        Dataset data("mydata", 1);

        expect_true(data.datasetName == "mydata");
        expect_false(data.isAligned);
        expect_true(data.getNumSamples() == 0);
        expect_true(data.numUnique == 0);

        Dataset data2(data);
        data2.datasetName = "new_name";
        vector<string> names = {"new_seq"};

        data2.addSequences(names);

        expect_true(data2.datasetName == "new_name");
        expect_false(data2.isAligned);
        expect_true(data2.getNumSamples() == 0);
        expect_true(data2.numUnique == 1);
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

        data.addSequences(names, seqs, comments);

        expect_true(data.numUnique == 3);
        expect_false(data.isAligned);

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
        Rcpp::DataFrame df = data.getSummary();

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
        Rcpp::DataFrame countTable = data.getSequenceAbundances(true);

        expect_true(countTable.size() == 2);

        expect_true(names == Rcpp::as<vector<string>>(countTable[0]));
        vector<float> abunds(3, 1);
        expect_true(abunds == Rcpp::as<vector<float>>(countTable[1]));

        data.clear();

        expect_true(data.numUnique == 0);

        data.addSequences(names, seqs, comments);

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

        countTable = data.getSequenceAbundances(true);

        expect_true(countTable.size() == 2);

        expect_true(ids == Rcpp::as<vector<string>>(countTable[0]));
        expect_true(abunds == Rcpp::as<vector<float>>(countTable[1]));
    }

    test_that("Tests addSeqs, assignSequenceAbundance, getSampleTotals, getTreatments, getTreatmentTotals, getSequenceAbundanceTable") {

        Dataset data("mydata", 1);

        vector<string> names(4, "");
        names[0] = "seq1"; names[1] = "seq2"; names[2] = "seq3"; names[3] = "seq4";

        vector<string> seqs(4, "");
        seqs[0] = "..ATGC-MGGT-AAA-TGC-NCT.";
        seqs[1] = "ATGCGGTAAATGCCT";
        seqs[2] = ".ATGCGGGGGTAAATGCCT.";
        seqs[3] = ".ATGCGAMMTAAATGCCT.";

        vector<string> comments(4, "");
        comments[2] = "my very cool comment";

        data.addSequences(names, seqs, comments);

        Dataset data2(data);

        expect_true(data2.getTotal() == 4);

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

        vector<float> abunds(10, 0);
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

        vector<string> badNames;
        expect_error(data.assignSequenceAbundance(badNames, abunds, samples));

        data.assignSequenceAbundance(ids, abunds, samples);

        vector<double> sampleTotals(3, 0);
        sampleTotals[0] = 301;
        sampleTotals[1] = 465;
        sampleTotals[2] = 554;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("samples")[1]) == sampleTotals);
        expect_true(data.getTotal() == 1320);

        vector<vector<string>> parsed = data.getSequenceNamesBySample(unique(samples));

        vector<string> uniqueIds = unique(ids);
        expect_true(parsed[0] == uniqueIds);
        uniqueIds.pop_back();
        expect_true(parsed[1] == uniqueIds);
        uniqueIds[2] = "seq4";
        expect_true(parsed[2] == uniqueIds);

        parsed = data.getSequencesBySample(unique(samples));

        expect_true(parsed[0] == seqs);
        seqs.pop_back();
        expect_true(parsed[1] == seqs);
        seqs[2] = ".ATGCGAMMTAAATGCCT.";
        expect_true(parsed[2] == seqs);

        // no treatment data
        Rcpp::DataFrame countTable = data.getSequenceAbundances(true);

        expect_true(ids == Rcpp::as<vector<string>>(countTable[0]));
        expect_true(abunds == Rcpp::as<vector<float>>(countTable[1]));
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

        expect_true(Rcpp::as<vector<double>>(data.getTotals("samples")[1]) == sampleTotals);
        expect_true(data.getTotal() == 1320);

        vector<double> treatmentTotals(2, 0);
        treatmentTotals[0] = 766;
        treatmentTotals[1] = 554;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("treatments")[1]) == treatmentTotals);

        vector<string> uniqueTreatments(2, "");
        uniqueTreatments[0] = "early";
        uniqueTreatments[1] = "late";

        expect_true(data.getTreatments() == uniqueTreatments);

        // assign treatments differently
        vector<string> newTreatments(3, "early");
        newTreatments[2] = "late";

        vector<string> newSamples(3, "");
        newSamples[0] = "sample2";
        newSamples[1] = "sample3";
        newSamples[2] = "sample4";

        data.assignTreatments(newSamples, newTreatments);

        expect_true(Rcpp::as<vector<double>>(data.getTotals("samples")[1]) == sampleTotals);
        expect_true(data.getTotal() == 1320);
        expect_true(Rcpp::as<vector<double>>(data.getTotals("treatments")[1]) == treatmentTotals);
        expect_true(data.getTreatments() == uniqueTreatments);

        countTable = data.getSequenceAbundances(true);

        expect_true(countTable.size() == 4);

        expect_true(ids == Rcpp::as<vector<string>>(countTable[0]));
        expect_true(abunds == Rcpp::as<vector<float>>(countTable[1]));
        expect_true(samples == Rcpp::as<vector<string>>(countTable[2]));
        expect_true(treatments == Rcpp::as<vector<string>>(countTable[3]));
    }

    test_that("Tests mergeSeqs, getScrapReport, getScrapSummary") {

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

        data.addSequences(names, seqs, comments);

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

        vector<float> abunds(7, 0);
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

        data.assignSequenceAbundance(ids, abunds, samples, treatments);

        expect_true(data.numUnique == 4);
        expect_true(data.getTotal() == 1269);

        vector<float> seqTotals(4, 0);
        seqTotals[0] = 1150;
        seqTotals[1] = 90;
        seqTotals[2] = 25;
        seqTotals[3] = 4;

        Rcpp::DataFrame df = data.getSequenceAbundances(false);
        expect_true(Rcpp::as<vector<float>>(df[1]) == seqTotals);

        vector<string> seqsToMerge(2, "seq2");
        seqsToMerge[0] = "seq1";

        // merges seq1 and seq2
        data.mergeSequences(seqsToMerge, "testMerge");

        seqTotals.resize(3);
        seqTotals[0] = 1240;
        seqTotals[1] = 25;
        seqTotals[2] = 4;

        df = data.getSequenceAbundances(false);
        expect_true(Rcpp::as<vector<float>>(df[1]) == seqTotals);
        expect_true(data.numUnique == 3);
        expect_true(data.getTotal() == 1269);

        // assign otus and merge seqs in different otus, should error
        vector<string> otus(3, "otu1");
        otus[1] = "otu2";
        otus[2] = "otu3";

        names.pop_back();
        names[0] = "seq1";
        names[1] = "seq3";
        names[2] = "seq4";

        seqsToMerge[1] = "seq3";

        data.assignBins(otus, nullFloatVector, nullVector, names);
        expect_error(data.mergeSequences(seqsToMerge));
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

        data.addSequences(names, seqs, comments);

        // mothur count file
        // Representative_Sequence     total   sample2	sample3	sample4
        // seq1	1150	250	400	500
        // seq2	90	    0	40	50
        // seq3	25	    0	25	0
        // seq4	4	    0	0	4

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

        vector<float> abunds(7, 0);
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

        data.assignSequenceAbundance(ids, abunds, samples, treatments);

        vector<double> sampleTotals(3, 0);
        sampleTotals[0] = 250;
        sampleTotals[1] = 465;
        sampleTotals[2] = 554;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("samples")[1]) == sampleTotals);
        expect_true(data.getTotal() == 1269);
        expect_true(data.numUnique == 4);
        expect_true(data.getNumSamples() == 3);
        expect_true(data.getNumTreatments() == 2);
        expect_true(data.hasSample("sample2"));
        expect_true(data.hasSample("sample3"));
        expect_true(data.hasSample("sample4"));
        expect_false(data.hasSample("badSample"));

        vector<double> treatmentTotals(2, 0);
        treatmentTotals[0] = 250;
        treatmentTotals[1] = 1019;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("treatments")[1]) == treatmentTotals);

        vector<string> uniqueTreatments(2, "");
        uniqueTreatments[0] = "early";
        uniqueTreatments[1] = "late";

        expect_true(data.getTreatments() == uniqueTreatments);

        vector<string> seqsToRemove(2, "seq2");
        seqsToRemove[0] = "seq1";
        vector<string> trashCodes(2, "testRemoval");
        trashCodes[0] = "testRemoveSampleTreatment";

        vector<string> badNames;
        expect_error(data.removeSequences(badNames, trashCodes));

        // removes 2 seqs, 1 sample and 1 treatment
        data.removeSequences(seqsToRemove, trashCodes);

        sampleTotals.resize(2, 0);
        sampleTotals[0] = 25;
        sampleTotals[1] = 4;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("samples")[1]) == sampleTotals);
        expect_true(data.getTotal() == 29);
        expect_true(data.numUnique == 2);
        expect_true(data.getNumSamples() == 2);
        expect_true(data.getNumTreatments() == 1);

        treatmentTotals.resize(1, 0);
        treatmentTotals[0] = 29;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("treatments")[1]) == treatmentTotals);

        uniqueTreatments.resize(1, "");
        uniqueTreatments[0] = "late";

        expect_true(data.getTreatments() == uniqueTreatments);
        vector<float> actual = Rcpp::as<vector<float>>(data.getSequenceAbundances()[1]);
        vector<float> expected = {25,4};
        expect_true(actual == expected);

        abunds.resize(1, 0);
        abunds[0] = 25;

        actual = Rcpp::as<vector<float>>(data.getSequenceAbundances(true)[1]);
        expect_true(actual == expected);

        Rcpp::DataFrame scrapReport = data.getScrapReport();
        expect_true(scrapReport.size() == 2);

        expect_true(seqsToRemove == Rcpp::as<vector<string>>(scrapReport[0]));
        expect_true(trashCodes == Rcpp::as<vector<string>>(scrapReport[1]));

        Rcpp::DataFrame scrapSummary = data.getScrapSummary();
        expect_true(scrapSummary.size() == 4);

        sort(trashCodes.begin(), trashCodes.end());
        vector<double> uniqueCounts(2, 1);
        vector<double> totalCounts(2, 1);
        totalCounts[0] = 90;
        totalCounts[1] = 1150;

        expect_true(trashCodes == Rcpp::as<vector<string>>(scrapSummary[1]));
        expect_true(uniqueCounts == Rcpp::as<vector<double>>(scrapSummary[2]));
        expect_true(totalCounts == Rcpp::as<vector<double>>(scrapSummary[3]));

        Rcpp::DataFrame df = data.getSummary("scrap");

        expect_true(trashCodes == Rcpp::as<vector<string>>(df[1]));
        expect_true(uniqueCounts == Rcpp::as<vector<double>>(df[2]));
        expect_true(totalCounts == Rcpp::as<vector<double>>(df[3]));
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

        data.addSequences(names, seqs, comments);

        expect_true(data.getTotal() == 4);
        expect_true(data.numUnique == 4);

        // add abundances with no samples
        vector<float> abunds(4, 0);
        abunds[0] = 250;
        abunds[1] = 400;
        abunds[2] = 500;
        abunds[3] = 10;

        data.assignSequenceAbundance(names, abunds);

        expect_true(data.getTotal() == 1160);
        expect_true(data.numUnique == 4);

        abunds[1] = 0;
        abunds[3] = 0;

        vector<string> badNames;
        expect_error(data.setAbundance(badNames, abunds));

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

        data.addSequences(names, seqs, comments);

        expect_true(data.getTotal() == 4);
        expect_true(data.numUnique == 4);

        data.assignSequenceAbundance(ids, abunds, samples, treatments);

        expect_true(data.getTotal() == 1269);
        expect_true(data.numUnique == 4);

        vector<vector<float>> parsedAbunds(4, vector<float>(3, 0));
        // removes sample2
        parsedAbunds[0][1] = 40; // sample3, late
        parsedAbunds[0][2] = 400; // sample4, late
        parsedAbunds[1][2] = 35; // sample4, late
        parsedAbunds[2][1] = 100; // sample3, late
        parsedAbunds[2][2] = 1; // sample4, late
        parsedAbunds[3][1] = 20; // sample3, late
        parsedAbunds[3][2] = 500; // sample4, late

        expect_error(data.setAbundances(badNames, parsedAbunds));

        // set abundances in dataset with samples and treatments
        data.setAbundances(names, parsedAbunds);

        expect_true(data.getTotal() == 1096);
        expect_true(data.numUnique == 4);

        vector<double> sampleTotals(2, 0);
        sampleTotals[0] = 160;
        sampleTotals[1] = 936;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("samples")[1]) == sampleTotals);
        expect_true(data.getNumSamples() == 2);
        expect_true(data.getNumTreatments() == 1);

        vector<double> treatmentTotals(1, 0);
        treatmentTotals[0] = 1096;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("treatments")[1]) == treatmentTotals);

        vector<string> uniqueTreatments(1, "late");

        expect_true(data.getTreatments() == uniqueTreatments);
    }

    test_that("Tests setSequences") {

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

        data.addSequences(names, seqs, comments);

        expect_true(data.getSequences() == seqs);
        expect_false(data.isAligned);

        seqs.clear();
        seqs.resize(4, ".AAATTT-CC-G.");
        comments[2] = "my_additional_comment";
        comments[3] = "newComment";

        vector<string> badNames;
        expect_error(data.setSequences(badNames, seqs, comments));

        // test with all seqs
        data.setSequences(names, seqs, comments);

        expect_true(data.getSequences() == seqs);
        expect_true(data.isAligned);
        expect_true(data.getTotal() == 4);

        seqs.pop_back(); names.pop_back();
        seqs.pop_back(); names.pop_back();
        seqs[0] = "AATT";

        // test with 2 seqs, no comments
        data.setSequences(names, seqs);

        seqs.push_back(".AAATTT-CC-G.");
        seqs.push_back(".AAATTT-CC-G.");

        expect_true(data.getSequences() == seqs);
        expect_false(data.isAligned);
        expect_true(data.getTotal() == 4);
    }

        // Bin tests
     test_that("Tests assignBins, getBinIds, getList, getRabund, getShared") {

        Dataset data("mydata", 1);

        expect_true(data.datasetName == "mydata");
        expect_false(data.isAligned);
        expect_true(data.getNumSamples() == 0);
        expect_true(data.numUnique == 0);

        expect_true(data.getBinAbundances().size() == 0);

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
        vector<float> abundances(10, 10);

        // test adding otuNames and abundances (rabund)
        data.assignBins(otuNames, abundances);

        expect_true(data.getTotal() == 100);
        expect_true(data.getNumBins() == 10);
        expect_true(data.getBinIds() == otuNames);

        expect_true(data.getList().size() == 0);
        Rcpp::DataFrame rabund = data.getBinAbundances();
        expect_true(otuNames == Rcpp::as<vector<string>>(rabund[0]));
        expect_true(abundances == Rcpp::as<vector<float>>(rabund[1]));
     }

    test_that("Tests mergeBins - abundance only") {

        Dataset data("mydata", 1);

        expect_true(data.datasetName == "mydata");
        expect_false(data.isAligned);
        expect_true(data.getNumSamples() == 0);
        expect_true(data.numUnique == 0);

        vector<string> otuNames = {"otu1", "otu2", "otu3", "otu4", "otu5",
                                   "otu6", "otu7", "otu8", "otu9", "otu10"};
        vector<float> abundances(10, 10);

        // test adding otuNames and abundances (rabund)
        data.assignBins(otuNames, abundances);

        expect_true(data.getTotal() == 100);
        expect_true(data.getNumBins() == 10);
        expect_true(data.getBinIds() == otuNames);

        // no sequence data was given
        expect_true(data.getList().size() == 0);
        vector<float> actual = Rcpp::as<vector<float>>(data.getBinAbundances()[1]);
        expect_true(actual == abundances);

        // merge otus without seqs or samples
        vector<string> otusToMerge(4, "otu1");
        otusToMerge[1] = "otu2";
        otusToMerge[2] = "otu4";
        otusToMerge[3] = "otu6";

        data.mergeBins(otusToMerge);

        expect_true(data.getTotal() == 100);
        expect_true(data.getNumBins() == 7);
        vector<float> binAbundances = Rcpp::as<vector<float>>(data.getBinAbundances()[1]);
        expect_true(binAbundances[0] == 40); // otu1
        expect_true(binAbundances[1] == 10); // otu3
        expect_true(binAbundances.size() == 7);
    }

    test_that("Tests mergeBins - with sequences") {

        Dataset data("mydata", 1);

        // test adding otuNames, seqNames, abundances (list)
        vector<string> otuNames = {"otu1", "otu1", "otu1", "otu2", "otu2",
                                   "otu3", "otu4", "otu4", "otu4", "otu4"};

        vector<string> seqNames = {"seq1", "seq2", "seq3", "seq4", "seq5",
                                   "seq6", "seq7", "seq8", "seq9", "seq10"};
        vector<float> abundances(10, 10);

        data.assignBins(otuNames, abundances, nullVector, seqNames);

        expect_true(data.getTotal() == 100);
        expect_true(data.getNumBins() == 4);
        expect_true(data.numUnique == 10);

        vector<float> expected = {30,20,10,40};
        vector<float> actual = Rcpp::as<vector<float>>(data.getBinAbundances()[1]);
        vector<string> otuIds = {"otu1", "otu2", "otu3", "otu4"};

        expect_true(actual == expected);
        expect_true(data.getBinIds() == otuIds);

        // test merge with seqids
        vector<string> otusToMerge = {"otu2", "otu4"};
        data.mergeBins(otusToMerge);
        otuIds.pop_back();

        expect_true(data.getTotal() == 100);
        expect_true(data.getNumBins() == 3);
        expect_true(data.numUnique == 10);

        vector<string> listVector = data.getListVector();
        expect_true(listVector[0] == "seq1,seq2,seq3"); // otu1
        expect_true(listVector[1] == "seq10,seq4,seq5,seq7,seq8,seq9"); // otu2
        expect_true(listVector[2] == "seq6"); // otu2
        expect_true(listVector.size() == 3);

        expected = {30,60,10};
        actual = Rcpp::as<vector<float>>(data.getBinAbundances()[1]);
        expect_true(actual == expected);
        expect_true(data.getTotal() == 100);
    }

    test_that("Tests removeBins, getScrapReport, getScrapSummary") {

        Dataset data("mydata", 1);

        // test adding otuNames, seqNames, abundances (list)
        vector<string> otuNames = {"otu1", "otu1", "otu1", "otu2", "otu2",
                                   "otu3", "otu4", "otu4", "otu4", "otu4"};

        vector<string> seqNames = {"seq1", "seq2", "seq3", "seq4", "seq5",
                                   "seq6", "seq7", "seq8", "seq9", "seq10"};
        vector<float> abundances(10, 10);

        data.assignBins(otuNames, abundances, nullVector, seqNames);

        vector<string> otusToMerge = {"otu2", "otu4"};
        data.mergeBins(otusToMerge);

        // remove sequence that will remove otu
        vector<string> seqToRemove(1, "seq6");
        vector<string> seqReason(1, "removeBin");

        data.removeSequences(seqToRemove, seqReason);

        expect_true(data.getTotal() == 90);
        expect_true(data.getNumBins() == 2);
        expect_true(data.numUnique == 9);

        vector<float> actual = Rcpp::as<vector<float>>(data.getBinAbundances()[1]);
        vector<float> expected = {30,60};
        expect_true(actual == expected);
        expect_true(actual.size() == 2);

        // remove otu by setting abundance to 0
        vector<string> testRemove(1, "otu2");
        vector<string> testTrash(1, "zeroedOTU");

        data.removeBins(testRemove, testTrash);

        expect_true(data.getTotal() == 30);
        expect_true(data.getNumBins() == 1);
        expect_true(data.numUnique == 3);

        actual = Rcpp::as<vector<float>>(data.getBinAbundances()[1]);
        expected = {30};
        expect_true(actual == expected);
        expect_true(actual.size() == 1);
    }

    test_that("Tests removeSamples") {

        Dataset data;

        // test adding otuNames, seqNames, abundances (list)
        vector<string> seqNames(16, "");
        vector<string> otuNames(16, "");
        vector<float> abundances(16, 10);
        vector<string> samples(16, "");
        otuNames[0] = "otu1"; seqNames[0] = "seq1";  samples[0] = "sample1";  abundances[0] = 10;
        otuNames[1] = "otu1"; seqNames[1] = "seq2";  samples[1] = "sample2";  abundances[1] = 10;
        otuNames[2] = "otu1"; seqNames[2] = "seq3";  samples[2] = "sample4";  abundances[2] = 5;
        otuNames[3] = "otu1"; seqNames[3] = "seq3";  samples[3] = "sample5";  abundances[3] = 5;
        otuNames[4] = "otu2"; seqNames[4] = "seq4";  samples[4] = "sample1";  abundances[4] = 5;
        otuNames[5] = "otu2"; seqNames[5] = "seq4"; samples[5] = "sample2";  abundances[5] = 5;
        otuNames[6] = "otu2"; seqNames[6] = "seq5";  samples[6] = "sample1";  abundances[6] = 10;
        otuNames[7] = "otu2"; seqNames[7] = "seq6";  samples[7] = "sample1";  abundances[7] = 10;
        otuNames[8] = "otu2"; seqNames[8] = "seq7";  samples[8] = "sample2";  abundances[8] = 10;
        otuNames[9] = "otu2"; seqNames[9] = "seq8";  samples[9] = "sample4";  abundances[9] = 10;
        otuNames[10] = "otu2"; seqNames[10] = "seq9";  samples[10] = "sample4";  abundances[10] = 5;
        otuNames[11] = "otu2"; seqNames[11] = "seq9"; samples[11] = "sample5";  abundances[11] = 5;
        otuNames[12] = "otu3"; seqNames[12] = "seq10";  samples[12] = "sample1";  abundances[12] = 1;
        otuNames[13] = "otu3"; seqNames[13] = "seq10";  samples[13] = "sample3";  abundances[13] = 2;
        otuNames[14] = "otu3"; seqNames[14] = "seq10";  samples[14] = "sample5";  abundances[14] = 3;
        otuNames[15] = "otu3"; seqNames[15] = "seq10";  samples[15] = "sample6";  abundances[15] = 4;

        data.assignBins(otuNames, abundances, samples, seqNames);

        vector<string> uniqueSamples(6, "");
        uniqueSamples[0] = "sample1";
        uniqueSamples[1] = "sample2";
        uniqueSamples[2] = "sample3";
        uniqueSamples[3] = "sample4";
        uniqueSamples[4] = "sample5";
        uniqueSamples[5] = "sample6";

        // assign treatments differently
        vector<string> treatments(6, "early");
        treatments[3] = "late";
        treatments[4] = "late";
        treatments[5] = "late";

        data.assignTreatments(uniqueSamples, treatments);

        Rcpp::DataFrame binAbunds = data.getBinAbundances();
        vector<float> actual = Rcpp::as<vector<float>>(binAbunds[1]);
        vector<float> expected = {30,60,10};
        expect_true(actual == expected);

        vector<double> sampleTotals(6, 0);
        sampleTotals[0] = 36;
        sampleTotals[1] = 25;
        sampleTotals[2] = 2;
        sampleTotals[3] = 20;
        sampleTotals[4] = 13;
        sampleTotals[5] = 4;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("samples")[1]) == sampleTotals);
        expect_true(data.getSamples() == uniqueSamples);

        vector<double> treatmentTotals(2, 0);
        treatmentTotals[0] = 63;
        treatmentTotals[1] = 37;

        expect_true(Rcpp::as<vector<double>>(data.getTotals("treatments")[1]) == treatmentTotals);
        expect_true(data.getTreatments() == unique(treatments));

        expect_true(data.getTotal() == 100);
        expect_true(data.getNumBins() == 3);
        expect_true(data.numUnique == 10);
        expect_true(data.getNumSamples() == 6);

        vector<string> samplesToRemove(2, "sample5");
        samplesToRemove[1] = "sample6";
        data.removeSamples(samplesToRemove);

        expect_true(data.getTotal() == 83);
        expect_true(data.getNumBins() == 3);
        expect_true(data.numUnique == 10);
        expect_true(data.getNumSamples() == 4);
    }

    test_that("Tests assignBins and assignSequenceAbundance together") {

        Dataset data;

        // test adding bins
        vector<string> seqNames(10, "");
        vector<string> otuNames(10, "");
        otuNames[0] = "otu1";   seqNames[0] = "seq1";
        otuNames[1] = "otu1";   seqNames[1] = "seq2";
        otuNames[2] = "otu1";   seqNames[2] = "seq3";
        otuNames[3] = "otu2";   seqNames[3] = "seq4";
        otuNames[4] = "otu2";   seqNames[4] = "seq5";
        otuNames[5] = "otu3";   seqNames[5] = "seq6";
        otuNames[6] = "otu4";   seqNames[6] = "seq10";
        otuNames[7] = "otu4";   seqNames[7] = "seq7";
        otuNames[8] = "otu4";   seqNames[8] = "seq8";
        otuNames[9] = "otu4";   seqNames[9] = "seq9";

        data.assignBins(otuNames, nullFloatVector, nullVector, seqNames);

        expect_true(data.getTotal() == 10);
        expect_true(data.getNumBins() == 4);
        expect_true(data.numUnique == 10);

        Rcpp::DataFrame binAbunds = data.getBinAbundances();
        vector<float> actual = Rcpp::as<vector<float>>(binAbunds[1]);
        vector<float> expected = {3,2,1,4};
        vector<string> list = data.getListVector();

        expect_true(actual == expected);
        expect_true(list[0] == "seq1,seq2,seq3");
        expect_true(list[1] == "seq4,seq5");
        expect_true(list[2] == "seq6");
        expect_true(list[3] == "seq10,seq7,seq8,seq9");

        // then assigning sequence abundance
        seqNames.resize(15);
        vector<string> samples(15, "");
        vector<float> abundances(15, 0);
        seqNames[0] = "seq1";   samples[0] = "sample1";  abundances[0] = 10;
        seqNames[1] = "seq1";   samples[1] = "sample2";  abundances[1] = 10;
        seqNames[2] = "seq2";   samples[2] = "sample1";  abundances[2] = 5;
        seqNames[3] = "seq2";   samples[3] = "sample2";  abundances[3] = 5;
        seqNames[4] = "seq3";   samples[4] = "sample1";  abundances[4] = 5;
        seqNames[5] = "seq3";   samples[5] = "sample2";  abundances[5] = 5;
        seqNames[6] = "seq4";   samples[6] = "sample2";  abundances[6] = 10;
        seqNames[7] = "seq5";   samples[7] = "sample1";  abundances[7] = 1;
        seqNames[8] = "seq6";   samples[8] = "sample1";  abundances[8] = 2;
        seqNames[9] = "seq6";   samples[9] = "sample2";  abundances[9] = 3;
        seqNames[10] = "seq7";   samples[10] = "sample1";  abundances[10] = 4;
        seqNames[11] = "seq8";   samples[11] = "sample1";  abundances[11] = 20;
        seqNames[12] = "seq9";   samples[12] = "sample1";  abundances[12] = 10;
        seqNames[13] = "seq9";   samples[13] = "sample2";  abundances[13] = 5;
        seqNames[14] = "seq10";   samples[14] = "sample2";  abundances[14] = 5;

        data.assignSequenceAbundance(seqNames, abundances, samples);

        expect_true(data.getTotal() == 100);
        expect_true(data.getNumBins() == 4);
        expect_true(data.numUnique == 10);

        binAbunds = data.getBinAbundances();
        actual = Rcpp::as<vector<float>>(binAbunds[1]);
        expected = {40,11,5,44};
        expect_true(actual == expected);
        expect_true(actual.size() == 4);
    }

}
