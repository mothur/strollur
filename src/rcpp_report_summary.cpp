#include <Rcpp.h>
#include "summary.h"

//' @name summarize_reports
//' @title summarize_reports
//' @rdname summarize_reports
//' @param report DataFrame, containing all numeric columns
//' @param count NumericalVector, containing abundances of sequences
//' @param processors Integer, number of cores to use. Default = all available
//' @description Summarizes a report
//' @return DataFrame with summary values
//[[Rcpp::export]]
Rcpp::DataFrame  summarize_reports(Rcpp::DataFrame& report,
                                    Rcpp::IntegerVector& count,
                                    int processors) {
     // create Summary object
     Summary* summary = new Summary(processors);

     Rcpp::CharacterVector colNames = report.attr("names");
     vector<string> names = Rcpp::as<vector<string>>(colNames);
     vector<int> c = Rcpp::as<vector<int>>(count);

     Rcpp::DataFrame results = summary->summarize(report, c, names);

     delete summary;

     return results;
 }
