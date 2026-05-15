#ifndef SRC_RCPP_DATASET
#define SRC_RCPP_DATASET

#include <Rcpp.h>
#include "../inst/include/strollur.h"
#include "rcpp_xint_xdev_functions.h"
#include "dataset.h"

/******************************************************************************/
//' @title clear
//' @description
//' Clear data from a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' clear(data)
//'
//' @return an updated \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @export
//[[Rcpp::export]]
Rcpp::Environment clear(Rcpp::Environment& data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->clear();

    data["raw"] = R_NilValue;
    data["sequence_tree"] = R_NilValue;
    data["sample_tree"] = R_NilValue;

    return data;
}
/******************************************************************************/
//' @title export_dataset
//' @description
//' Export all data from a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 2)
//' export_dataset(dataset)
//'
//' @return Rcpp::List, containing the data in the 'Dataset
//' @export
//[[Rcpp::export]]
Rcpp::List export_dataset(Rcpp::Environment data) {

    Rcpp::XPtr<Dataset> d = data["data"];
    Rcpp::List results = d.get()->exportDataset();
    vector<string> resultNames = results.names();

    Rcpp::Function get = data["get_sequence_tree"];
    Rcpp::List sequence_tree = get();

    if (sequence_tree.size() != 0) {
        results.push_back(sequence_tree);
        resultNames.push_back("sequence_tree");
    }

    get = data["get_sample_tree"];
    Rcpp::List sample_tree = get();

    if (sample_tree.size() != 0) {
        results.push_back(sample_tree);
        resultNames.push_back("sample_tree");
    }

    results.attr("names") = resultNames;
    results.attr("strollur_version") = "1.0.0";
    results.attr("dataset_name") = d.get()->datasetName;

    return results;
}
/******************************************************************************/
//' @title get_bin_types
//' @description
//' Get bin table types of a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @examples
//'
//' data <- miseq_sop_example()
//' get_bin_types(data)
//'
//' @return vector of strings
//' @export
//[[Rcpp::export]]
vector<string> get_bin_types(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getBinTypes();
}
/******************************************************************************/
//' @title has_sample
//' @description
//' Determine if a given sample is in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//' @param sample a string containing the name of a sample.
//' @examples
//'
//' data <- miseq_sop_example()
//' has_sample(data, "F3D0")
//' has_sample(data, "not a valid sample")
//'
//' @return boolean indicating whether the dataset has a given sample
//' @export
//[[Rcpp::export]]
bool has_sample(Rcpp::Environment data, string sample) {
     Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->hasSample(sample);
}
/******************************************************************************/
//' @title has_sequence_strings
//' @description
//' Determine if a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object contains sequence nucleotide strings.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//' @examples
//'
//' data <- miseq_sop_example()
//' has_sequence_strings(data)
//'
//' @return boolean indicating whether the dataset has sequence nucleotide
//' strings.
//' @export
//[[Rcpp::export]]
bool has_sequence_strings(Rcpp::Environment data) {
     Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->hasSeqs();
}
/******************************************************************************/
//' @title is_aligned
//' @description
//' Determine if a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object contains aligned sequences.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @examples
//'
//' dataset <- miseq_sop_example()
//' is_aligned(dataset)
//'
//' @return Boolean
//' @export
//[[Rcpp::export]]
bool is_aligned(Rcpp::Environment data) {
     Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->isAligned;
}
/******************************************************************************/
#endif
