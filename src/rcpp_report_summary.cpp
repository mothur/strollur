#include <Rcpp.h>
#include "summary.h"

//' @name summarize_reports
//' @title summarize_reports
//' @rdname summarize_reports
//' @param report DataFrame, containing all numeric columns
//' @param count NumericalVector, containing abundances of sequences
//' @param processors Integer, number of cores to use. Default = all available
//' @description Summarizes a report
//' @examples
//'
//' alignment_report <- readr::read_tsv(rdataset_example("alignment_data.tsv"),
//'    col_names = TRUE, show_col_types = FALSE)
//'
//' # for this simple example we assume all the sequences in the report are
//' # unique, therefore we set their abundance to 1.
//' counts <- rep(1, 5)
//'
//' # extract numeric columns to summarize
//' numeric_report <- alignment_report[sapply(alignment_report, is.numeric)]
//'
//' report_summary <- summarize_reports(numeric_report, counts, 10)
//'
//' @return DataFrame with summary values
//[[Rcpp::export]]
Rcpp::DataFrame  summarize_reports(Rcpp::DataFrame& report,
                                    Rcpp::NumericVector& count,
                                    int processors) {
     // create Summary object
     Summary* summary = new Summary(processors);

     Rcpp::CharacterVector colNames = report.attr("names");
     vector<string> names = Rcpp::as<vector<string>>(colNames);
     vector<float> c = Rcpp::as<vector<float>>(count);

     Rcpp::DataFrame results = summary->summarize(report, c, names);

     delete summary;

     return results;
 }
