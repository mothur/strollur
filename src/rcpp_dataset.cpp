#ifndef SRC_RCPP_DATASET
#define SRC_RCPP_DATASET

#include <Rcpp.h>
#include "../inst/include/strollur.h"
#include "rcpp_xint_xdev_functions.h"
#include "dataset.h"

/******************************************************************************/
//' @title new_dataset
//' @description
//' Create a new \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param dataset_name string, a string containing the dataset name.
//' Default = ""
//' @param processors integer, number of cores to use during summary functions.
//' Default = all available
//' @examples
//'
//' data <- new_dataset()
//'
//' # to create a new dataset named "soil" and allow for all available
//' # processors during summary functions, run the following:
//'
//' data <- new_dataset(dataset_name = "soil")
//'
//' # to create a new dataset named "soil" and allow for 2
//' # processors during summary functions, run the following:
//'
//' data <- new_dataset(dataset_name = "soil", processors = 2)
//'
//' @returns a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @export
//' @seealso The 'new' method in the \href{https://mothur.org/strollur/reference/strollur.html}{strollur} class
//[[Rcpp::export]]
Rcpp::Environment new_dataset(string dataset_name = "",
                              Rcpp::Nullable<int> processors = R_NilValue) {

    // strollur$new()
    Rcpp::Environment strollur_env("package:strollur");
    Rcpp::Environment dataset_class_env = strollur_env["strollur"];
    Rcpp::Function constructor = dataset_class_env["new"];

    int num_proc = 1;
    if (processors.isNotNull()) {
        num_proc = Rcpp::as<int>(processors);
    }else{
        Rcpp::Function getProc = dataset_class_env["get_available_processors"];
        num_proc = Rcpp::as<int>(getProc());
    }

    Rcpp::Environment data = constructor(dataset_name, num_proc, R_NilValue);
    return data;
}
/******************************************************************************/
//' @title new_reference
//' @description
//' Create a reference you can add to your dataset
//'
//' @param reference_name, a string containing the name of the reference used
//' in the preparing of the sequences. For example: 'silva.bacteria.fasta'.
//'
//' @param reference_version, a string containing the version of the reference
//' used in the preparing of the sequences. For example: '1.38.1'. Default = "".
//'
//' @param reference_usage, a string containing the usage of the reference in
//' your analysis. For example: 'alignment using mothur2' Default = NULL.
//'
//' @param reference_note, a string containing the any additional notes about
//' the reference. Default = "".
//'
//' @param reference_url, a string containing a web address where the
//' reference may be downloaded. Default = "".
//' @examples
//'
//' reference <- new_reference("silva.bacteria.fasta",
//'                            "1.38.1",
//'                            "alignment by mothur2 v1.0 using default options",
//'                            "",
//'                            "https://mothur.org/wiki/silva_reference_files/")
//'
//' @returns a list
//' @export
//' @seealso [add()]
//[[Rcpp::export]]
Rcpp::List new_reference(string reference_name,
                         string reference_version = "",
                         string reference_usage = "",
                         string reference_note = "",
                         string reference_url = "") {

     return Rcpp::List::create(
         Rcpp::Named("reference_name") = reference_name,
         Rcpp::Named("reference_version") = reference_version,
         Rcpp::Named("reference_usage") = reference_usage,
         Rcpp::Named("reference_note") = reference_note,
         Rcpp::Named("reference_url") = reference_url
     );
 }
/******************************************************************************/
//' @title copy_dataset
//' @description
//' Create a new \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object from an existing dataset.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @examples
//'
//' miseq <- miseq_sop_example()
//'
//' # to create a new dataset that is a copy of miseq
//'
//' data <- copy_dataset(miseq)
//'
//' @returns a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @seealso The 'new' method in the \href{https://mothur.org/strollur/reference/strollur.html}{strollur} class
//' @export
//[[Rcpp::export]]
Rcpp::Environment copy_dataset(Rcpp::Environment data) {

    // strollur$new()
    Rcpp::Environment strollur_env("package:strollur");
    Rcpp::Environment dataset_class_env = strollur_env["strollur"];
    Rcpp::Function constructor = dataset_class_env["new"];
    Rcpp::XPtr<Dataset> d = data["data"];

    Rcpp::Environment copy = constructor(d.get()->datasetName,
                                         d.get()->processors,
                                         data);

    return copy;
}
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
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void clear(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->clear();

    data["raw"] = R_NilValue;
    data["sequence_tree"] = R_NilValue;
    data["sample_tree"] = R_NilValue;
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
