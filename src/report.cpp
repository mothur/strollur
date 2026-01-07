#include "../inst/include/rdataset.h"
#include "dataset.h"
#include "summary.h"

/******************************************************************************/
Report::Report() {
    sequence_name = "";
    sequence_name_col = -1;
    hasStr = false;
    hasInt = false;
    hasNum = false;
    hasLog = false;
    hasColumnNames = false;
    hasReport = false;
    numRows = 0;
}
/******************************************************************************/
void Report::clear() {
    sequence_name = "";
    sequence_name_col = -1;
    hasStr = false;
    hasInt = false;
    hasNum = false;
    hasLog = false;
    hasColumnNames = false;
    hasReport = false;
    numRows = 0;

    columnNames.clear();
    strColumns.clear();
    intColumns.clear();
    numColumns.clear();
    logColumns.clear();
}
/******************************************************************************/
void Report::addReport(Rcpp::DataFrame& report) {
    clear();

    vector<string> cnames;

    // save columnNames
    if (report.hasAttribute("names")) {
        cnames = Rcpp::as<vector<string>>(report.attr("names"));
        hasColumnNames = true;
    }

    // save sequence_name
    if (report.hasAttribute("sequence_name")) {
        sequence_name = Rcpp::as<string>(report.attr("sequence_name"));
    }

    numRows = report.nrow();

    // save columns
    for (int i = 0; i < report.size(); i++) {

        bool addColName = true;

        // vector<int>
        if (TYPEOF(report[i]) == INTSXP) {
            intColumns[i] = Rcpp::as<vector<int> >(report[i]);
            hasInt = true;
        }
        // vector<double>
        else if (TYPEOF(report[i]) == REALSXP) {
            numColumns[i] = Rcpp::as<vector<double> >(report[i]);
            hasNum = true;
        }
        // vector<string>
        else if (TYPEOF(report[i]) == STRSXP) {
            strColumns[i] = Rcpp::as<vector<string> >(report[i]);
            if (hasColumnNames) {
                if (sequence_name == cnames[i]) {
                    sequence_name_col = i;
                }
            }
            hasStr = true;
        }
        // vector<bool>
        else if (TYPEOF(report[i]) == LGLSXP) {
            logColumns[i] = Rcpp::as<vector<bool> >(report[i]);
            hasLog = true;
        }else{
            addColName = false;
        }

        // only save column names of columns we are saving
        if (addColName) {
            hasReport = true;
            if (hasColumnNames) {
                columnNames[i] = cnames[i];
            }else{
                columnNames[i] = toString(i);
            }
        }
    }
}
/******************************************************************************/
Rcpp::DataFrame Report::getReport(set<string> datasetNames) {

    Rcpp::DataFrame data = Rcpp::DataFrame::create();

    pruneReport(datasetNames);

    if (hasReport) {
        // assemble dataframe in order
        for (auto itCol = columnNames.begin(); itCol != columnNames.end(); itCol++) {
            bool foundNext = false;

            if (hasStr) {
                foundNext = addNextColumn(data, strColumns, itCol->first);
            }

            if (hasInt && !foundNext) {
                foundNext = addNextColumn(data, intColumns, itCol->first);
            }

            if (hasNum && !foundNext) {
                foundNext = addNextColumn(data, numColumns, itCol->first);
            }

            if (hasLog && !foundNext) {
                foundNext = addNextColumn(data, logColumns, itCol->first);
            }
        }

        if (hasColumnNames) {
            data.attr("names") = getValues(columnNames);
        }

        if (sequence_name != "") {
            data.attr("sequence_name") = sequence_name;
        }
    }

    return data;
}
/******************************************************************************/
Rcpp::DataFrame Report::summarizeReport(set<string> datasetNames,
                                        int processors, vector<float> counts) {

    Rcpp::DataFrame results = Rcpp::DataFrame::create();

    pruneReport(datasetNames);

    if (hasReport) {
        // create Summary object
        Summary* summary = new Summary(processors);

        Rcpp::DataFrame report = Rcpp::DataFrame::create();
        vector<string> reportNames;

        for (auto it = intColumns.begin(); it != intColumns.end(); it++) {
            report.push_back(it->second);
            reportNames.push_back(columnNames[it->first]);
        }

        for (auto it = numColumns.begin(); it != numColumns.end(); it++) {
            report.push_back(it->second);
            reportNames.push_back(columnNames[it->first]);
        }

        results = summary->summarize(report, counts, reportNames);

        delete summary;
    }

    return results;
}
/******************************************************************************/
void Report::pruneReport(set<string> datasetNames) {

    // metadata case
    if (datasetNames.size() == 0) { return; }

    if (datasetNames.size() > strColumns[sequence_name_col].size()) {
        // error, report is missing names. This can happen if you add sequences
        // after adding a report
        string message = "Your dataset contains more sequences than your report";
        message += ", removing report.";
        RcppThread::Rcout << endl << message << endl;
        clear();
    }else {

        // create a boolean filter, assume guilty
        vector<bool> filter(numRows, false);

        vector<string> reportSequenceNames = strColumns[sequence_name_col];

        for (int i = 0; i < reportSequenceNames.size(); i++) {
            // if the name is in the dataset, set filter pos to true
            if (datasetNames.count(reportSequenceNames[i]) != 0) {
                filter[i] = true;
            }
        }

        // sanity check
        int trueCount = count(filter.begin(), filter.end(), true);

        if (trueCount < datasetNames.size()) {
            // error, report is missing names. This can happen if you add sequences
            // after adding a report
            string message = "Your dataset contains more sequences than your report";
            message += ", removing report.";
            RcppThread::Rcout << endl << message << endl;
            clear();
            return;
        }

        // filter all columns
        numRows = trueCount;
        if (hasInt) {
            for (auto it = intColumns.begin(); it != intColumns.end(); it++) {
                it->second = select(it->second, filter);
            }
        }
        if (hasStr) {
            for (auto it = strColumns.begin(); it != strColumns.end(); it++) {
                it->second = select(it->second, filter);
            }
        }
        if (hasNum) {
            for (auto it = numColumns.begin(); it != numColumns.end(); it++) {
                it->second = select(it->second, filter);
            }
        }
        if (hasLog) {
            for (auto it = logColumns.begin(); it != logColumns.end(); it++) {
                it->second = select(it->second, filter);
            }
        }
    }
}
/******************************************************************************/
