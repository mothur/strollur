#include <testthat.h>
#include "../inst/include/rdataset.h"

context("BinTable class C++ unit tests") {

    test_that("Tests constructor") {
        BinTable otuTable("0.03");

        expect_true(otuTable.label == "0.03");
        expect_true(otuTable.getNumBins() == 0);
    }

    test_that("Tests add, get, getAbundance, getAbundances") {
        BinTable otuTable("0.03");

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
        AbundTable ab;
        otuTable.assignAbundance(otuNames, abundances,
                                 nullVector, nullIntVector, ab);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.getNumBins() == 10);

        // no sequence data was given
        expect_true(otuTable.getAbundance("otu1") == 10);
        vector<float> temp(1, 10);
        expect_true(otuTable.getAbundances("otu1") == temp);

        otuTable.clear();

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
        otuTable.assignAbundance(otuNames, abundances, samples,
                                 nullIntVector, ab);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.getNumBins() == 4);
        expect_true(otuTable.getNumSamples() == 6);

        expect_true(otuTable.getAbundance("otu1") == 30);
        expect_true(otuTable.getAbundance("otu2") == 20);
        expect_true(otuTable.getAbundance("otu3") == 10);
        expect_true(otuTable.getAbundance("otu4") == 40);
        temp.resize(6, 0);
        temp[0] = 5;
        temp[1] = 5;
        temp[3] = 10;
        expect_true(otuTable.getAbundances("otu2") == temp);

        vector<double> sampleTotals(6, 0);
        sampleTotals[0] = 36;
        sampleTotals[1] = 25;
        sampleTotals[2] = 2;
        sampleTotals[3] = 20;
        sampleTotals[4] = 13;
        sampleTotals[5] = 4;

        expect_true(otuTable.getSampleTotals() == sampleTotals);
    }

    test_that("Tests getScrapReport, getRAbundVector, getRAbund, getSharedVector, hasSample") {
        BinTable otuTable("0.03");

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
        vector<float> abundances(10, 10);

        // test adding otuNames and abundances (rabund)
        AbundTable ab;
        otuTable.assignAbundance(otuNames, abundances,
                                 nullVector, nullIntVector, ab);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.getNumBins() == 10);
        expect_true(otuTable.getAbundance("otu1") == 10);
        vector<float> temp(1, 10);
        expect_true(otuTable.getAbundances("otu1") == temp);

        otuTable.clear();

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

        vector<vector<float> > sharedVector(4, vector<float>(6, 0));
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
        otuTable.assignAbundance(otuNames, abundances, samples,
                                 nullIntVector, ab);

        expect_true(otuTable.getTotal() == 100);
        expect_true(otuTable.getNumBins() == 4);
        expect_true(otuTable.getNumSamples() == 6);

        vector<float> otuRabundVector(4, 0);
        otuRabundVector[0] = 30;
        otuRabundVector[1] = 20;
        otuRabundVector[2] = 10;
        otuRabundVector[3] = 40;

        vector<string> uniqueBinNames(4);
        uniqueBinNames[0] = "otu1";
        uniqueBinNames[1] = "otu2";
        uniqueBinNames[2] = "otu3";
        uniqueBinNames[3] = "otu4";

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
        expect_true(uniqueBinNames == Rcpp::as<vector<string>>(rabund[0]));
        expect_true(otuRabundVector == Rcpp::as<vector<float>>(rabund[1]));

        Rcpp::DataFrame shared = otuTable.getShared();
        expect_true(otuNames == Rcpp::as<vector<string>>(shared[0]));
        expect_true(abundances == Rcpp::as<vector<float>>(shared[1]));
        expect_true(samples == Rcpp::as<vector<string>>(shared[2]));

        vector<double> sampleTotals(6, 0);
        sampleTotals[0] = 36;
        sampleTotals[1] = 25;
        sampleTotals[2] = 2;
        sampleTotals[3] = 20;
        sampleTotals[4] = 13;
        sampleTotals[5] = 4;

        expect_true(otuTable.getSampleTotals() == sampleTotals);

        expect_true(otuTable.getScrapReport().size() == 0);
        expect_true(otuTable.getScrapSummary().size() == 0);

        expect_error(otuTable.setAbundance(nullVector, otuRabundVector));
        expect_error(otuTable.setAbundances(nullVector, sharedVector));
    }
}
