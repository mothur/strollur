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
