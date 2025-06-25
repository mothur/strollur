#include <testthat.h>
#include "../inst/include/rdataset.h"

context("OtuTable class C++ unit tests") {

    test_that("Tests constructor") {
        OtuTable otuTable("0.03");

        expect_true(otuTable.label == "0.03");
        expect_true(otuTable.numOtus == 0);
        expect_true(otuTable.numUnique == -1);
    }

    test_that("Tests add, get, getAbundance, getAbundances") {
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
        // no sequence data was given
        expect_true(otuTable.get("otu1") == "");
        expect_true(otuTable.getAbundance("otu1") == 10);
        vector<int> temp(1, 10);
        expect_true(otuTable.getAbundances("otu1") == temp);

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

        expect_true(otuTable.get("otu1") == "seq1,seq2,seq3");
        expect_true(otuTable.get("otu2") == "seq4,seq5");
        expect_true(otuTable.get("otu3") == "seq6");
        expect_true(otuTable.get("otu4") == "seq7,seq8,seq9,seq10");
        expect_true(otuTable.getAbundance("otu1") == 30);
        expect_true(otuTable.getAbundance("otu2") == 20);
        expect_true(otuTable.getAbundance("otu3") == 10);
        expect_true(otuTable.getAbundance("otu4") == 40);
        temp[0] = 20;
        expect_true(otuTable.getAbundances("otu2") == temp);

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

        expect_true(otuTable.get("otu4") == "seq7,seq8,seq9,seq10");
        expect_true(otuTable.getAbundance("otu1") == 30);
        expect_true(otuTable.getAbundance("otu2") == 20);
        expect_true(otuTable.getAbundance("otu3") == 10);
        expect_true(otuTable.getAbundance("otu4") == 40);
        temp.resize(6, 0);
        temp[0] = 5;
        temp[1] = 5;
        temp[3] = 10;
        expect_true(otuTable.getAbundances("otu2") == temp);

        vector<int> sampleTotals(6, 0);
        sampleTotals[0] = 36;
        sampleTotals[1] = 25;
        sampleTotals[2] = 2;
        sampleTotals[3] = 20;
        sampleTotals[4] = 13;
        sampleTotals[5] = 4;

        expect_true(otuTable.getSampleTotals() == sampleTotals);
    }

    test_that("Tests getScrapReport, getRAbundVector, getRAbund, getListVector, getList, getSharedVector, hasSample") {
        OtuTable otuTable("0.03");

        Rcpp::DataFrame rabund = otuTable.getRAbund();
        expect_true(rabund.size() == 0);

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
        // no sequence data was given
        expect_true(otuTable.get("otu1") == "");
        expect_true(otuTable.getAbundance("otu1") == 10);
        vector<int> temp(1, 10);
        expect_true(otuTable.getAbundances("otu1") == temp);

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

        vector<string> uniqueOtuNames(4);
        uniqueOtuNames[0] = "otu1";
        uniqueOtuNames[1] = "otu2";
        uniqueOtuNames[2] = "otu3";
        uniqueOtuNames[3] = "otu4";

        otuTable.add(otuNames, abundances, nullVector, seqNames);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.numOtus == 4);
        expect_true(otuTable.numUnique == 10);

        vector<string> otuListVector(4, "");
        otuListVector[0] = "seq1,seq2,seq3";
        otuListVector[1] = "seq4,seq5";
        otuListVector[2] = "seq6";
        otuListVector[3] = "seq7,seq8,seq9,seq10";

        vector<int> otuRabundVector(4, 0);
        otuRabundVector[0] = 30;
        otuRabundVector[1] = 20;
        otuRabundVector[2] = 10;
        otuRabundVector[3] = 40;

        expect_true(otuTable.get("otu1") == otuListVector[0]);
        expect_true(otuTable.get("otu2") == otuListVector[1]);
        expect_true(otuTable.get("otu3") == otuListVector[2]);
        expect_true(otuTable.get("otu4") == otuListVector[3]);
        expect_true(otuTable.getListVector() == otuListVector);
        expect_true(otuTable.getAbundance("otu1") == otuRabundVector[0]);
        expect_true(otuTable.getAbundance("otu2") == otuRabundVector[1]);
        expect_true(otuTable.getAbundance("otu3") == otuRabundVector[2]);
        expect_true(otuTable.getAbundance("otu4") == otuRabundVector[3]);
        expect_true(otuTable.getRAbundVector() == otuRabundVector);
        temp[0] = 20;
        expect_true(otuTable.getAbundances("otu2") == temp);

        Rcpp::DataFrame list = otuTable.getList();
        expect_true(otuNames == Rcpp::as<vector<string>>(list[0]));
        expect_true(seqNames == Rcpp::as<vector<string>>(list[1]));

        Rcpp::DataFrame shared = otuTable.getShared();
        expect_true(shared.size() == 0);

        rabund = otuTable.getRAbund();
        expect_true(uniqueOtuNames == Rcpp::as<vector<string>>(rabund[0]));
        expect_true(otuRabundVector == Rcpp::as<vector<int>>(rabund[1]));

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

        vector<vector<int> > sharedVector(4, vector<int>(6, 0));
        // otu1
        sharedVector[0][0] = 10;
        sharedVector[0][1] = 10;
        sharedVector[0][3] = 5;
        sharedVector[0][4] = 5;
        // otu2
        sharedVector[1][0] = 5;
        sharedVector[1][1] = 5;
        sharedVector[1][3] = 10;
        // otu3
        sharedVector[2][0] = 1;
        sharedVector[2][2] = 2;
        sharedVector[2][4] = 3;
        sharedVector[2][5] = 4;
        // otu4
        sharedVector[3][0] = 20;
        sharedVector[3][1] = 10;
        sharedVector[3][3] = 5;
        sharedVector[3][4] = 5;

        // add shared data
        otuTable.add(otuNames, abundances, samples);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.numOtus == 4);
        expect_true(otuTable.numUnique == 10);
        expect_true(otuTable.getNumSamples() == 6);

        expect_true(otuTable.get("otu4") == otuListVector[3]);
        expect_true(otuTable.getAbundance("otu1") == otuRabundVector[0]);
        expect_true(otuTable.getAbundance("otu2") == otuRabundVector[1]);
        expect_true(otuTable.getAbundance("otu3") == otuRabundVector[2]);
        expect_true(otuTable.getAbundance("otu4") == otuRabundVector[3]);
        expect_true(otuTable.getSharedVector() == sharedVector);
        expect_true(otuTable.hasSample("sample4"));
        expect_false(otuTable.hasSample("badSample"));
        temp.resize(6, 0);
        temp[0] = 5;
        temp[1] = 5;
        temp[3] = 10;
        expect_true(otuTable.getAbundances("otu2") == temp);

        rabund = otuTable.getRAbund();
        expect_true(uniqueOtuNames == Rcpp::as<vector<string>>(rabund[0]));
        expect_true(otuRabundVector == Rcpp::as<vector<int>>(rabund[1]));

        shared = otuTable.getShared();
        expect_true(otuNames == Rcpp::as<vector<string>>(shared[0]));
        expect_true(abundances == Rcpp::as<vector<int>>(shared[1]));
        expect_true(samples == Rcpp::as<vector<string>>(shared[2]));

        vector<int> sampleTotals(6, 0);
        sampleTotals[0] = 36;
        sampleTotals[1] = 25;
        sampleTotals[2] = 2;
        sampleTotals[3] = 20;
        sampleTotals[4] = 13;
        sampleTotals[5] = 4;

        expect_true(otuTable.getSampleTotals() == sampleTotals);

        expect_true(otuTable.getScrapReport().size() == 0);
        expect_true(otuTable.getScrapSummary().size() == 0);

        expect_error(otuTable.setAbundance(nullVector, sampleTotals));
        expect_error(otuTable.setAbundances(nullVector, sharedVector));
    }


}
