#include <testthat.h>
#include "../inst/include/strollur.h"
#include "dataset.h"

 context("Report class C++ unit tests") {

    test_that("test Report") {

        Report report;

        expect_false(report.hasReport);

        vector<string> names = {"seq1", "seq2", "seq3", "seq4", "seq5"};
        vector<int> ints = {1,2,3,4,5};
        vector<double> doubs = {10.5,20.5,30.5,40.5,50.5};
        vector<bool> logs = {true, true, false, true, false};

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("sequence_names") = names,
            Rcpp::_["intCol"] = ints,
            Rcpp::_["numCol"] = doubs,
            Rcpp::_["logCol"] = logs);

        df.attr("sequence_name") = "sequence_names";

        report.addReport(df);

        expect_true(report.hasReport);

        Rcpp::DataFrame df2 = report.getReport(toSet(names));

        vector<string> resultStr = Rcpp::as<vector<string>>(df2[0]);

        expect_true(resultStr.size() == 5);
        expect_true(resultStr == names);
        expect_true(Rcpp::as<vector<int>>(df2[1]) == ints);
        expect_true(Rcpp::as<vector<double>>(df2[2]) == doubs);
        expect_true(Rcpp::as<vector<bool>>(df2[3]) == logs);

        vector<float> flo = {1,2,3,4,5};
        Rcpp::DataFrame summary = report.summarizeReport(toSet(names), 1, flo);

        vector<int> summaryResults = {1,2,3,4,5,5,5,3};

        vector<float> summaryFlo = Rcpp::as<vector<float>>(summary[1]);
        flo = {10.5,20.5,30.5,40.5,50.5,50.5,50.5,37.1667};

        // * 100 rounded
        for (size_t i = 0; i < flo.size(); i++) {
            summaryFlo[i] = ceil(summaryFlo[i]*100);
            flo[i] = ceil(flo[i]*100);
        }

        expect_true(Rcpp::as<vector<int>>(summary[0]) == summaryResults);
        expect_true(summaryFlo == flo);

        // remove seq5
        names.pop_back();
        ints.pop_back();
        doubs.pop_back();
        logs.pop_back();

        df2 = report.getReport(toSet(names));
        expect_true(Rcpp::as<vector<string>>(df2[0]).size() == 4);
        expect_true(Rcpp::as<vector<string>>(df2[0]) == names);
        expect_true(Rcpp::as<vector<int>>(df2[1]) == ints);
        expect_true(Rcpp::as<vector<double>>(df2[2]) == doubs);
        expect_true(Rcpp::as<vector<bool>>(df2[3]) == logs);

        // get report asking for more names than present
        names.push_back("missingName");
        df2 = report.getReport(toSet(names));

        expect_false(report.hasReport);
    }
 }
