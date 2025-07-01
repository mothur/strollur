#include <testthat.h>
#include "../inst/include/rdataset.h"
#include "dataset.h"

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

        data.addSequences(names, seqs, comments);

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

        vector<string> badNames;
        badNames = names;
        badNames[0] = "badName";

        expect_error(data.addContigsReport(badNames, oLengths, oStarts,
                                           oEnds, mismatches, ee));

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

        data.clear();
        data.addSequences(names, seqs, comments);
        ee.pop_back();
        // warning given for size mismatches and no data added
        data.addContigsReport(names, oLengths, oStarts,
                              oEnds, mismatches, ee);

        // no contigs data because of size mismatch
        contigsReport = data.getContigsReport();

        expect_true(contigsReport.size() == 0);

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

        data.addSequences(names, seqs, comments);

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

        vector<string> badNames;
        badNames = names;
        badNames[0] = "badName";

        expect_error(data.addAlignReport(badNames, searchScores,
                                          simScores, longestInserts));

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

        data.clear();
        data.addSequences(names, seqs, comments);
        searchScores.pop_back();
        // warning given for size mismatches and no data added
        data.addAlignReport(names, searchScores, simScores, longestInserts);

        // no contigs data because of size mismatch
        alignReport = data.getAlignReport();

        expect_true(alignReport.size() == 0);
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

        data.addSequences(names, seqs, comments);

        expect_true(data.getNames() == names);
        expect_true(data.getSequences() == seqs);

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

        data.addSequences(names, seqs, comments);

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

        vector<string> badNames;
        expect_error(data.assignSequenceAbundance(badNames, abunds, samples));

        data.assignSequenceAbundance(ids, abunds, samples);

        vector<int> sampleTotals(3, 0);
        sampleTotals[0] = 301;
        sampleTotals[1] = 465;
        sampleTotals[2] = 554;

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getTotal() == 1320);

        vector<vector<string>> parsed = data.getNamesBySample(unique(samples));

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

        // assign treatments differently
        vector<string> newTreatments(3, "early");
        newTreatments[2] = "late";

        vector<string> newSamples(3, "");
        newSamples[0] = "sample2";
        newSamples[1] = "sample3";
        newSamples[2] = "sample4";

        data.assignTreatments(newSamples, newTreatments);

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getTotal() == 1320);
        expect_true(data.getTreatmentTotals() == treatmentTotals);
        expect_true(data.getTreatments() == uniqueTreatments);

        countTable = data.getSequenceAbundanceTable();

        expect_true(countTable.size() == 4);

        expect_true(ids == Rcpp::as<vector<string>>(countTable[0]));
        expect_true(abunds == Rcpp::as<vector<int>>(countTable[1]));
        expect_true(samples == Rcpp::as<vector<string>>(countTable[2]));
        expect_true(treatments == Rcpp::as<vector<string>>(countTable[3]));
        expect_true(data.getOtuIds() == nullVector);
        expect_true(data.getOtuAbundances("otu1") == nullIntVector);
        expect_true(data.getOtuAbundance("otu1") == 0);
        expect_true(data.getOtu("otu1") == "");

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
        treatments[1] = "late";
        treatments[2] = "late";
        treatments[3] = "late";
        treatments[4] = "late";
        treatments[5] = "late";
        treatments[6] = "late";

        data.assignSequenceAbundance(ids, abunds, samples, treatments);

        expect_true(data.numUnique == 4);
        expect_true(data.getTotal() == 1269);

        vector<int> seqTotals(4, 0);
        seqTotals[0] = 1150;
        seqTotals[1] = 90;
        seqTotals[2] = 25;
        seqTotals[3] = 4;

        expect_true(data.getSequenceAbundances() == seqTotals);

        vector<string> seqsToMerge(2, "seq2");
        seqsToMerge[0] = "seq1";

        // merges seq1 and seq2
        data.mergeSequences(seqsToMerge, "testMerge");

        seqTotals.resize(3);
        seqTotals[0] = 1240;
        seqTotals[1] = 25;
        seqTotals[2] = 4;

        expect_true(data.getSequenceAbundances() == seqTotals);
        expect_true(data.numUnique == 3);
        expect_true(data.getTotal() == 1269);
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
        treatments[1] = "late";
        treatments[2] = "late";
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
        expect_true(data.hasSample("sample2"));
        expect_true(data.hasSample("sample3"));
        expect_true(data.hasSample("sample4"));
        expect_false(data.hasSample("badSample"));

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

        vector<string> badNames;
        expect_error(data.removeSequences(badNames, trashCodes));

        // removes 2 seqs, 1 sample and 1 treatment
        data.removeSequences(seqsToRemove, trashCodes);

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
        expect_true(data.getAbundance("seq1") == 0);
        expect_true(data.getAbundance("seq2") == 0);
        expect_true(data.getAbundance("seq3") == 25);
        expect_true(data.getAbundance("seq4") == 4);

        abunds.resize(2, 0);
        abunds[0] = 25;
        abunds[1] = 0;
        expect_true(data.getAbundances("seq3") == abunds);

        abunds[0] = 0;
        abunds[1] = 4;
        expect_true(data.getAbundances("seq4") == abunds);

        Rcpp::DataFrame scrapReport = data.getScrapReport();
        expect_true(scrapReport.size() == 2);

        expect_true(seqsToRemove == Rcpp::as<vector<string>>(scrapReport[0]));
        expect_true(trashCodes == Rcpp::as<vector<string>>(scrapReport[1]));

        Rcpp::List list = data.getScrapSummary();
        Rcpp::DataFrame scrapSummary(list["sequence_scrap_summary"]);
        expect_true(scrapSummary.size() == 3);

        sort(trashCodes.begin(), trashCodes.end());
        vector<int> uniqueCounts(2, 1);
        vector<int> totalCounts(2, 1);
        totalCounts[0] = 90;
        totalCounts[1] = 1150;

        expect_true(trashCodes == Rcpp::as<vector<string>>(scrapSummary[0]));
        expect_true(uniqueCounts == Rcpp::as<vector<int>>(scrapSummary[1]));
        expect_true(totalCounts == Rcpp::as<vector<int>>(scrapSummary[2]));

        Rcpp::List summary = data.getSequenceSummary();
        Rcpp::DataFrame df(summary["scrap_summary"]);

        expect_true(trashCodes == Rcpp::as<vector<string>>(df[0]));
        expect_true(uniqueCounts == Rcpp::as<vector<int>>(df[1]));
        expect_true(totalCounts == Rcpp::as<vector<int>>(df[2]));

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

        vector<vector<int>> parsedAbunds(4, vector<int>(3, 0));
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

    // Otu tests
    test_that("Tests assignOtus, getOtuIds, getList, getRabund, getShared") {

        Dataset data("mydata", 1);

        expect_true(data.datasetName == "mydata");
        expect_false(data.isAligned);
        expect_true(data.numSamples == 0);
        expect_true(data.numUnique == 0);

        expect_true(data.getRAbund().size() == 0);

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
        data.assignOtus(otuNames, abundances);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 10);
        expect_true(data.getOtuIds() == otuNames);

        expect_true(data.getList().size() == 0);
        Rcpp::DataFrame rabund = data.getRAbund();
        expect_true(otuNames == Rcpp::as<vector<string>>(rabund[0]));
        expect_true(abundances == Rcpp::as<vector<int>>(rabund[1]));

        // no sequence data was given
        expect_true(data.getOtu("otu1") == "");
        expect_true(data.getOtuAbundance("otu1") == 10);
        vector<int> temp(1, 10);
        expect_true(data.getOtuAbundances("otu1") == temp);

        data.clear();

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

        data.assignOtus(otuNames, abundances, nullVector, seqNames);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 4);
        expect_true(data.numUnique == 10);

        expect_true(data.getOtu("otu1") == "seq1,seq2,seq3");
        expect_true(data.getOtu("otu2") == "seq4,seq5");
        expect_true(data.getOtu("otu3") == "seq6");
        expect_true(data.getOtu("otu4") == "seq7,seq8,seq9,seq10");
        expect_true(data.getListVector()[0] == "seq1,seq2,seq3");
        expect_true(data.getListVector()[1] == "seq4,seq5");
        expect_true(data.getListVector()[2] == "seq6");
        expect_true(data.getListVector()[3] == "seq7,seq8,seq9,seq10");
        expect_true(data.getOtuAbundance("otu1") == 30);
        expect_true(data.getOtuAbundance("otu2") == 20);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 40);
        temp[0] = 20;
        expect_true(data.getOtuAbundances("otu2") == temp);
        expect_true(data.getShared().size() == 0);

        vector<string> otuIds(4, "otu1");
        otuIds[1] = "otu2";
        otuIds[2] = "otu3";
        otuIds[3] = "otu4";

        expect_true(data.getOtuIds() == otuIds);

        Rcpp::DataFrame list = data.getList();
        expect_true(otuNames == Rcpp::as<vector<string>>(list[0]));
        expect_true(seqNames == Rcpp::as<vector<string>>(list[1]));

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
        data.assignOtus(otuNames, abundances, samples);

        vector<string> uniqueSamples(6, "");
        uniqueSamples[0] = "sample1";
        uniqueSamples[1] = "sample2";
        uniqueSamples[2] = "sample3";
        uniqueSamples[3] = "sample4";
        uniqueSamples[4] = "sample5";
        uniqueSamples[5] = "sample6";

        // assign treatments differently
        vector<string> treatments(6, "early");

        expect_error(data.assignTreatments(uniqueSamples, unique(treatments)));

        data.assignTreatments(uniqueSamples, treatments);

        vector<int> treatmentTotals(1, 0);
        treatmentTotals[0] = 100;

        expect_true(data.getTreatmentTotals() == treatmentTotals);
        expect_true(data.getTreatments() == unique(treatments));

        treatments[3] = "late";
        treatments[4] = "late";
        treatments[5] = "late";

        data.assignTreatments(uniqueSamples, treatments);

        treatmentTotals.resize(2, 0);
        treatmentTotals[0] = 63;
        treatmentTotals[1] = 37;

        expect_true(data.getTreatmentTotals() == treatmentTotals);
        expect_true(data.getTreatments() == unique(treatments));

        Rcpp::DataFrame shared = data.getShared();
        expect_true(otuNames == Rcpp::as<vector<string>>(shared[0]));
        expect_true(abundances == Rcpp::as<vector<int>>(shared[1]));
        expect_true(samples == Rcpp::as<vector<string>>(shared[2]));

        uniqueSamples.push_back("SampleNotInDataset");
        treatments.push_back("badEntry");

        data.assignTreatments(uniqueSamples, treatments);

        // ignored bad entry
        expect_true(data.getTreatmentTotals() == treatmentTotals);
        expect_true(data.numTreatments == 2);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 4);
        expect_true(data.numUnique == 10);
        expect_true(data.numSamples == 6);

        vector<int> sampleTotals(6, 0);
        sampleTotals[0] = 36;
        sampleTotals[1] = 25;
        sampleTotals[2] = 2;
        sampleTotals[3] = 20;
        sampleTotals[4] = 13;
        sampleTotals[5] = 4;

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getOtuIds() == otuIds);
        expect_true(data.getSamples() == unique(samples));

        expect_true(data.getOtu("otu4") == "seq7,seq8,seq9,seq10");
        expect_true(data.getOtuAbundance("otu1") == 30);
        expect_true(data.getOtuAbundance("otu2") == 20);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 40);
        expect_true(data.getOtuAbundances("badotu") == nullIntVector);
        temp.resize(6, 0);
        temp[0] = 5;
        temp[1] = 5;
        temp[3] = 10;
        expect_true(data.getOtuAbundances("otu2") == temp);

        otuNames.clear();
        expect_error(data.assignOtus(otuNames, abundances));
    }

   test_that("Tests mergeOtus, removeOtus, getScrapReport, getScrapSummary") {

        Dataset data("mydata", 1);

        expect_true(data.datasetName == "mydata");
        expect_false(data.isAligned);
        expect_true(data.numSamples == 0);
        expect_true(data.numUnique == 0);

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
        data.assignOtus(otuNames, abundances);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 10);
        expect_true(data.getOtuIds() == otuNames);

        // no sequence data was given
        expect_true(data.getOtu("otu1") == "");
        expect_true(data.getOtuAbundance("otu1") == 10);
        vector<int> temp(1, 10);
        expect_true(data.getOtuAbundances("otu1") == temp);

        // merge otus without seqs or samples
        vector<string> otusToMerge(4, "otu1");
        otusToMerge[1] = "otu2";
        otusToMerge[2] = "otu4";
        otusToMerge[3] = "otu6";
        data.mergeOtus(otusToMerge);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 7);
        expect_true(data.getOtuAbundance("otu1") == 40);
        expect_true(data.getOtuAbundance("otu2") == 0);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 0);
        expect_true(data.getOtuAbundance("otu6") == 0);
        data.clear();

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

        data.assignOtus(otuNames, abundances, nullVector, seqNames);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 4);
        expect_true(data.numUnique == 10);

        expect_true(data.getOtu("otu1") == "seq1,seq2,seq3");
        expect_true(data.getOtu("otu2") == "seq4,seq5");
        expect_true(data.getOtu("otu3") == "seq6");
        expect_true(data.getOtu("otu4") == "seq7,seq8,seq9,seq10");
        expect_true(data.getOtuAbundance("otu1") == 30);
        expect_true(data.getOtuAbundance("otu2") == 20);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 40);
        temp[0] = 20;
        expect_true(data.getOtuAbundances("otu2") == temp);

        vector<string> otuIds(4, "otu1");
        otuIds[1] = "otu2";
        otuIds[2] = "otu3";
        otuIds[3] = "otu4";

        expect_true(data.getOtuIds() == otuIds);

        // test merge with seqids
        otusToMerge.resize(2);
        otusToMerge[0] = "otu2";
        otusToMerge[1] = "otu4";

        data.mergeOtus(otusToMerge);
        otuIds.pop_back();

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 3);
        expect_true(data.numUnique == 10);

        expect_true(data.getOtu("otu1") == "seq1,seq2,seq3");
        expect_true(data.getOtu("otu2") == "seq4,seq5,seq7,seq8,seq9,seq10");
        expect_true(data.getOtu("otu3") == "seq6");
        expect_true(data.getOtu("otu4") == "");
        expect_true(data.getOtuAbundance("otu1") == 30);
        expect_true(data.getOtuAbundance("otu2") == 60);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 0);
        temp[0] = 60;
        expect_true(data.getOtuAbundances("otu2") == temp);
        expect_true(data.getTotal() == 100);

        // remove sequence that will remove otu
        vector<string> seqToRemove(1, "seq6");
        vector<string> seqReason(1, "removeOtu");

        data.removeSequences(seqToRemove, seqReason);

        expect_true(data.getTotal() == 90);
        expect_true(data.numOtus == 2);
        expect_true(data.numUnique == 9);
        expect_true(data.getOtuAbundance("otu3") == 0);
        expect_true(data.getOtu("otu3") == "");

        // remove otu by setting abundance to 0
        vector<string> testRemove(1, "otu2");
        vector<int> r(6, 0);
        vector<vector<int>> abundsRemove; abundsRemove.push_back(r);

        data.setOtuAbundances(testRemove, abundsRemove, "zeroedOTU");

        expect_true(data.getTotal() == 30);
        expect_true(data.numOtus == 1);
        expect_true(data.numUnique == 3);
        expect_true(data.getOtuAbundance("otu2") == 0);
        expect_true(data.getOtu("otu2") == "");

        // test adding otuNames abundances, samples and seqids
        data.clear();
        otuNames.resize(16, "otu1");
        abundances.resize(16, 10);
        seqNames.resize(16, "");
        vector<string> samples(16, "sample1");
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

        data.assignOtus(otuNames, abundances, samples, seqNames);

        vector<string> uniqueSamples(6, "");
        uniqueSamples[0] = "sample1";
        uniqueSamples[1] = "sample2";
        uniqueSamples[2] = "sample3";
        uniqueSamples[3] = "sample4";
        uniqueSamples[4] = "sample5";
        uniqueSamples[5] = "sample6";

        // assign treatments differently
        vector<string> treatments(6, "early");

        expect_error(data.assignTreatments(uniqueSamples, unique(treatments)));

        data.assignTreatments(uniqueSamples, treatments);

        expect_true(data.getOtuAbundance("otu1") == 30);
        expect_true(data.getOtuAbundance("otu2") == 60);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 0);

        vector<int> treatmentTotals(1, 0);
        treatmentTotals[0] = 100;

        expect_true(data.getTreatmentTotals() == treatmentTotals);
        expect_true(data.getTreatments() == unique(treatments));

        vector<int> sampleTotals(6, 0);
        sampleTotals[0] = 36;
        sampleTotals[1] = 25;
        sampleTotals[2] = 2;
        sampleTotals[3] = 20;
        sampleTotals[4] = 13;
        sampleTotals[5] = 4;

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getSamples() == uniqueSamples);

        treatments[3] = "late";
        treatments[4] = "late";
        treatments[5] = "late";

        data.assignTreatments(uniqueSamples, treatments);

        treatmentTotals.resize(2, 0);
        treatmentTotals[0] = 63;
        treatmentTotals[1] = 37;

        expect_true(data.getTreatmentTotals() == treatmentTotals);
        expect_true(data.getTreatments() == unique(treatments));

        uniqueSamples.push_back("SampleNotInDataset");
        treatments.push_back("badEntry");

        data.assignTreatments(uniqueSamples, treatments);

        // ignored bad entry
        expect_true(data.getTreatmentTotals() == treatmentTotals);
        expect_true(data.numTreatments == 2);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 3);
        expect_true(data.numUnique == 10);
        expect_true(data.numSamples == 6);

        expect_true(data.getSampleTotals() == sampleTotals);
        expect_true(data.getOtuIds() == otuIds);
        expect_true(data.getSamples() == unique(samples));

        expect_true(data.getOtu("otu4") == "");
        expect_true(data.getOtu("otu2") == "seq4,seq5,seq6,seq7,seq8,seq9");
        expect_true(data.getOtu("otu1") == "seq1,seq2,seq3");
        expect_true(data.getOtu("otu3") == "seq10");
        expect_true(data.getOtuAbundance("otu1") == 30);
        expect_true(data.getOtuAbundance("otu2") == 60);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 0);
        temp.resize(6, 0);
        temp[0] = 25;
        temp[1] = 15;
        temp[3] = 15;
        temp[4] = 5;
        expect_true(data.getOtuAbundances("otu2") == temp);
        expect_true(data.getOtuAbundances("badotu") == nullIntVector);

        auto countMatrix = data.getSeqsAbundsBySample();
        auto goodNames = data.getNames();
        vector<int> abunds(6, 0);
        //seq1
        abunds[0] = 10;
        expect_true(countMatrix[0] == abunds);
        // seq10
        abunds[0] = 1; abunds[2] = 2; abunds[4] = 3; abunds[5] = 4;
        expect_true(countMatrix[1] == abunds);
        // seq2
        abunds[0] = 0; abunds[1] = 10;
        abunds[2] = 0; abunds[4] = 0; abunds[5] = 0;
        expect_true(countMatrix[2] == abunds);

        // remove otu
        vector<string> otusToRemove(2, "otu1");
        otusToRemove[1] = "non_existant_otu";
        vector<string> reasonsToRemove(2, "badOtu");
        data.removeOtus(otusToRemove, reasonsToRemove);

        expect_true(data.getOtu("otu1") == "");
        expect_true(data.getOtuAbundance("otu1") == 0);
        expect_true(data.getOtuAbundance("otu2") == 60);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 0);

        expect_true(data.getTotal() == 70);
        expect_true(data.numOtus == 2);
        expect_true(data.numUnique == 7);
        expect_true(data.numSamples == 6);

        // getScrapReport
        Rcpp::DataFrame scrapReport = data.getScrapReport("otu");
        expect_true(scrapReport.size() == 2);

        otusToRemove.clear();
        otusToRemove.push_back("otu1");
        reasonsToRemove.clear();
        reasonsToRemove.push_back("badOtu");

        expect_true(otusToRemove == Rcpp::as<vector<string>>(scrapReport[0]));
        expect_true(reasonsToRemove == Rcpp::as<vector<string>>(scrapReport[1]));

        Rcpp::List list = data.getScrapSummary();
        Rcpp::DataFrame scrapSummary(list["otu_scrap_summary"]);
        expect_true(scrapSummary.size() == 3);

        vector<int> uniqueCounts(1, 1);
        vector<int> totalCounts(1, 30);

        expect_true(reasonsToRemove == Rcpp::as<vector<string>>(scrapSummary[0]));
        expect_true(uniqueCounts == Rcpp::as<vector<int>>(scrapSummary[1]));
        expect_true(totalCounts == Rcpp::as<vector<int>>(scrapSummary[2]));

        otuNames.clear();
        expect_error(data.assignOtus(otuNames, abundances));
    }

    test_that("Tests setOtuAbundance, setOtuAbundances") {

        Dataset data("mydata", 1);

        expect_true(data.datasetName == "mydata");
        expect_false(data.isAligned);
        expect_true(data.numSamples == 0);
        expect_true(data.numUnique == 0);

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
        data.assignOtus(otuNames, abundances);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 10);
        expect_true(data.numUnique == 0);
        expect_true(data.getOtuIds() == otuNames);

        // no sequence data was given
        expect_true(data.getOtu("otu1") == "");
        expect_true(data.getOtuAbundance("otu1") == 10);
        vector<int> temp(1, 10);
        expect_true(data.getOtuAbundances("otu1") == temp);

        // set abundance of otus
        vector<string> otusToChange(4, "otu1");
        otusToChange[1] = "otu2";
        otusToChange[2] = "otu4";
        otusToChange[3] = "otu6";
        vector<int> abundsToChange(4, 0);
        abundsToChange[1] = 20;
        abundsToChange[2] = 30;

        data.setOtuAbundance(otusToChange, abundsToChange, "zeroAbundance");

        expect_true(data.getTotal() == 110);
        expect_true(data.numOtus == 8);
        expect_true(data.numUnique == 0);
        expect_true(data.getOtuAbundance("otu1") == 0);
        expect_true(data.getOtuAbundance("otu2") == 20);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 30);
        expect_true(data.getOtuAbundance("otu6") == 0);
        data.clear();

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

        data.assignOtus(otuNames, abundances, nullVector, seqNames);

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 4);
        expect_true(data.numUnique == 10);

        expect_true(data.getOtu("otu1") == "seq1,seq2,seq3");
        expect_true(data.getOtu("otu2") == "seq4,seq5");
        expect_true(data.getOtu("otu3") == "seq6");
        expect_true(data.getOtu("otu4") == "seq7,seq8,seq9,seq10");
        expect_true(data.getOtuAbundance("otu1") == 30);
        expect_true(data.getOtuAbundance("otu2") == 20);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 40);
        temp[0] = 20;
        expect_true(data.getOtuAbundances("otu2") == temp);

        otusToChange.resize(2);
        otusToChange[0] = "otu1";
        otusToChange[1] = "otu4";
        abundsToChange.resize(2);
        abundsToChange[0] = 70;
        abundsToChange[1] = 0;

        data.setOtuAbundance(otusToChange, abundsToChange, "zeroAbundance");

        expect_true(data.getTotal() == 100);
        expect_true(data.numOtus == 3);
        expect_true(data.numUnique == 6);
        expect_true(data.getOtuAbundance("otu1") == 70);
        expect_true(data.getOtuAbundance("otu2") == 20);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 0);

        data.clear();

        // test adding otuNames abundances, samples (shared)
        otuNames.resize(12, "otu1");
        abundances.resize(12, 10);
        vector<string> samples(12, "sample1");
        otuNames[0] = "otu1";   samples[0] = "sample1";  abundances[0] = 10;
        otuNames[1] = "otu1";   samples[1] = "sample2";  abundances[1] = 10;
        otuNames[2] = "otu1";   samples[2] = "sample4";  abundances[2] = 5;
        otuNames[3] = "otu1";   samples[3] = "sample5";  abundances[3] = 5;
        otuNames[4] = "otu2";   samples[4] = "sample1";  abundances[4] = 25;
        otuNames[5] = "otu2";   samples[5] = "sample2";  abundances[5] = 15;
        otuNames[6] = "otu2";   samples[6] = "sample4";  abundances[6] = 15;
        otuNames[7] = "otu2";  samples[7] = "sample5";  abundances[7] = 5;
        otuNames[8] = "otu3";   samples[8] = "sample1";  abundances[8] = 1;
        otuNames[9] = "otu3";   samples[9] = "sample3";  abundances[9] = 2;
        otuNames[10] = "otu3";   samples[10] = "sample5";  abundances[10] = 3;
        otuNames[11] = "otu3";   samples[11] = "sample6";  abundances[11] = 4;

        // add shared data
        data.assignOtus(otuNames, abundances, samples);

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

        expect_true(data.getOtuAbundance("otu1") == 30);
        expect_true(data.getOtuAbundance("otu2") == 60);
        expect_true(data.getOtuAbundance("otu3") == 10);
        expect_true(data.getOtuAbundance("otu4") == 0);

        vector<int> treatmentTotals(2, 0);
        treatmentTotals[0] = 63;
        treatmentTotals[1] = 37;

        expect_true(data.getTreatmentTotals() == treatmentTotals);
        expect_true(data.getTreatments() == unique(treatments));

        otusToChange.resize(2);
        otusToChange[0] = "otu1";
        otusToChange[1] = "otu3";
        // remove otu1, modify otu3 (removing sample 6)
        vector<vector<int>> abunds(2, vector<int>(6,0));
        abunds[1][0] = 1;
        abunds[1][2] = 2;
        abunds[1][4] = 3;

        data.setOtuAbundances(otusToChange, abunds);

        expect_true(data.getOtuAbundance("otu1") == 0);
        expect_true(data.getOtuAbundance("otu2") == 60);
        expect_true(data.getOtuAbundance("otu3") == 6);
        expect_true(data.getOtuAbundance("otu4") == 0);
        expect_true(data.getTotal() == 66);
        expect_true(data.numOtus == 2);
        expect_true(data.numSamples == 5);
        expect_true(data.getTotal("sample6") == 0);
    }
}
