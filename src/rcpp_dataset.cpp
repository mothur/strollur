#ifndef SRC_RCPP_DATASET
#define SRC_RCPP_DATASET

#include <Rcpp.h>
#include "../inst/include/rdataset.h"
#include "rcpp_xint_xdev_functions.h"
#include "dataset.h"

/******************************************************************************/
//' @title get_available_processors
//' @name get_available_processors
//' @description
//' Get the number of available cores
// [[Rcpp::export]]
int get_available_processors() {
     // Use Rcpp::Environment and Rcpp::Function to call R code from C++.
     Rcpp::Environment parallelly_env = Rcpp::Environment::namespace_env("parallelly");
     Rcpp::Function availableCores = parallelly_env["availableCores"];

     // Call the R function and return the result.
     return Rcpp::as<int>(availableCores());
}
/******************************************************************************/
//' @title new_dataset
//' @description
//' Create a new \link{dataset} object
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
//' @returns a \link{dataset} object
//' @seealso The 'new' method in the \link{dataset} class
//[[Rcpp::export]]
Rcpp::Environment new_dataset(string dataset_name = "",
                              Rcpp::Nullable<int> processors = R_NilValue) {

    // dataset$new()
    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Environment dataset_class_env = rdataset_env["dataset"];
    Rcpp::Function constructor = dataset_class_env["new"];

    int num_proc = 1;
    if (processors.isNotNull()) {
        num_proc = Rcpp::as<int>(processors);
    }else{
        num_proc = get_available_processors();
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
//' Create a new \link{dataset} object from an existing dataset.
//'
//' @param data, a \link{dataset} object
//' @examples
//'
//' miseq <- miseq_sop_example()
//'
//' # to create a new dataset that is a copy of miseq
//'
//' data <- copy_dataset(miseq)
//'
//' @returns a \link{dataset} object
//' @seealso The 'new' method in the \link{dataset} class
//[[Rcpp::export]]
Rcpp::Environment copy_dataset(Rcpp::Environment data) {

    // dataset$new()
    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Environment dataset_class_env = rdataset_env["dataset"];
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
//' Clear data from a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param tags a vector of strings containing the items you wish to clear.
//' Options are 'sequence_data', 'bin_data', 'metadata',
//' 'references', 'sequence_tree', 'sample_tree' and 'reports'. By default,
//' everything is cleared.
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' clear(data)
//'
//[[Rcpp::export]]
void clear(Rcpp::Environment data,
           Rcpp::CharacterVector tags = Rcpp::CharacterVector::create()) {
    Rcpp::XPtr<Dataset> d = data["data"];
    data["raw"] = R_NilValue;

    vector<string> t = Rcpp::as<vector<string>>(tags);
    bool hasTags = false;
    if (t.size() > 0) { hasTags = true; }

    d.get()->clear(t);

    Rcpp::Function clear = data["clear"];
    if (!hasTags) {
        clear();
    }else {
        // remove tags for c++ back end
        vector<string> goodTags;
        for (string tag : t) {
            if ((tag != "references") && (tag != "sequence_data") &&
                (tag != "bin_data")) {
                goodTags.push_back(tag);
            }
        }
        if (goodTags.size() != 0) { clear(goodTags); }
    }
}
/******************************************************************************/
//' @title export_dataset
//' @description
//' Export all data from an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param tags a vector of strings containing the items you wish to export.
//' Options are 'sequence_data' and 'bin_data', 'metadata',
//' 'references', 'sequence_tree', 'sample_tree', and 'reports'.
//' By default, everything is exported.
//'
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 2)
//' export_dataset(dataset)
//'
//' @return Rcpp::List, containing the data in the 'Dataset
//[[Rcpp::export]]
Rcpp::List export_dataset(Rcpp::Environment data,
                          Rcpp::CharacterVector tags = Rcpp::CharacterVector::create()) {

    vector<string> t = Rcpp::as<vector<string>>(tags);
    bool hasTags = false;
    if (t.size() > 0) { hasTags = true; }

    Rcpp::XPtr<Dataset> d = data["data"];
    Rcpp::List results = d.get()->exportDataset(t);
    vector<string> resultNames = results.names();

    if ((!hasTags) || (vectorContains(t, "sequence_tree"))) {

        Rcpp::Function get = data["get_sequence_tree"];
        Rcpp::List sequence_tree = get();

        if (sequence_tree.size() != 0) {
            results.push_back(sequence_tree);
            resultNames.push_back("sequence_tree");
        }
    }

    if ((!hasTags) || (vectorContains(t, "sample_tree"))) {

        Rcpp::Function get = data["get_sample_tree"];
        Rcpp::List sample_tree = get();

        if (sample_tree.size() != 0) {
            results.push_back(sample_tree);
            resultNames.push_back("sample_tree");
        }
    }

    results.attr("names") = resultNames;
    results.attr("rdataset_version") = "1.0.0";
    results.attr("dataset_name") = d.get()->datasetName;

    return results;
}
/******************************************************************************/
//' @title get_bin_representative_sequences
//' @description
//' Get the representative sequences of the bins in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param bin_type, string indicating the type of clusters. Default = "otu".
//' @examples
//'
//'   miseq <- miseq_sop_example()
//'
//'   num_bins <- count(data = miseq, type = "bins", bin_type = "otu")
//'
//'   # For examples sake, select first 531 sequences to be the representatives
//'   table <- data.frame(bin_names = names(data = miseq,
//'                                        type = "bins",
//'                                        bin_type = "otu"),
//'                       sequence_names = names(data = miseq,
//'                                             type = "sequences")[1:num_bins]
//'                       )
//'
//'   assign(data = miseq, table = table,
//'          type = "bin_representatives", bin_type = "otu")
//'
//'
//'   get_bin_representative_sequences(data = miseq, bin_type= "otu")
//'
//' @return data.frame containing 2 columns representative_names and
//'  representative_sequences
//[[Rcpp::export]]
Rcpp::DataFrame get_bin_representative_sequences(Rcpp::Environment data,
                                                 string bin_type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getBinRepresentativeSequences(bin_type);
}
/******************************************************************************/
//' @title get_bin_types
//' @description
//' Get bin table types of a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @examples
//'
//' data <- miseq_sop_example()
//' get_bin_types(data)
//'
//' @return vector of strings
//[[Rcpp::export]]
vector<string> get_bin_types(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getBinTypes();
}
/******************************************************************************/
//' @title get_list
//' @description
//' Get data frame containing sequence bin assignments of a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//' data <- miseq_sop_example()
//' get_list(data)
//'
//' @return 2 column data.frame sequences assigned to bins
//[[Rcpp::export]]
Rcpp::DataFrame get_list(Rcpp::Environment data, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getList(type);
}
/******************************************************************************/
//' @title get_list_vector
//' @description
//' Get vector of strings containing the sequences bin data
//'
//' @param data, a \link{dataset} object
//'
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//' data <- miseq_sop_example()
//' get_list_vector(data)
//'
//' @return vector of strings containing the names of the sequences in each bin
//' separated by commas
//[[Rcpp::export]]
vector<string> get_list_vector(Rcpp::Environment data, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getListVector(type);
}
/******************************************************************************/
//' @title report
//' @description
//' Get a data.frame containing the given report in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param type, string containing the type of report you would like. Options
//' include: "sequence_data", "sequence_taxonomy", "bin_taxonomy",
//' "sequence_scrap", "bin_scrap", "metadata", "references". If you have added
//' custom reports for alignment, contigs_assembly or chimeras, you can get those
//' as well. Default = "sequence_data".
//'
//' @param bin_type, string containing the bin type you would like a bin_taxonomy
//' report for. Default = "otu".
//'
//' @examples
//'
//' # First let's create a dataset from the \href{https://mothur.org/wiki/miseq_sop/}{MiSeq_SOP}
//'
//' miseq <- miseq_sop_example()
//'
//' # To get a report about the FASTA data
//'
//' sequence_report <- report(data = miseq, type = "sequence_data")
//' head(sequence_report, n = 10)
//'
//' # To get a report about sequence classifications
//'
//' sequence_taxonomy_report <- report(data = miseq,
//'                                        type = "sequence_taxonomy")
//' head(sequence_taxonomy_report, n = 10)
//'
//' # To get a report about bin classifications for 'otu' data
//'
//' otu_taxonomy_report <- report(data = miseq,
//'                                        type = "bin_taxonomy",
//'                                        bin_type = "otu")
//' head(otu_taxonomy_report, n = 10)
//'
//' # To get a report about bin classifications for 'asv' data
//'
//' asv_taxonomy_report <- report(data = miseq, type = "bin_taxonomy",
//'                               bin_type = "asv")
//' head(asv_taxonomy_report, n = 10)
//'
//' # To get a report about bin classifications for 'phylotype' data
//'
//' phylotype_taxonomy_report <- report(data = miseq, type = "bin_taxonomy",
//'                                     bin_type = "phylotype")
//' head(phylotype_taxonomy_report, n = 10)
//'
//' # To get a report about the sequences removed during your analysis:
//'
//' scrapped_sequence_report <- report(data = miseq, type = "sequence_scrap")
//'
//' # To get a report about the "otu" bins removed during your analysis:
//'
//' scrapped_otu_report <- report(data = miseq, type = "bin_scrap",
//'                               bin_type = "otu")
//'
//' # To get a report about the "phylotype" bins removed during your analysis:
//'
//' scrapped_phylotype_report <- report(data = miseq, type = "bin_scrap",
//'                                     bin_type = "phylotype")
//'
//' # To get the metadata associated with your data:
//'
//' metadata <- report(data = miseq, type = "metadata")
//'
//' # To get the resource references associated with your data:
//'
//' references <- report(data = miseq, type = "references")
//'
//' # To get our custom report containing the contigs assembly data:
//'
//' contigs_report <- report(data = miseq, type = "contigs_report")
//' head(contigs_report, n = 10)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame report(Rcpp::Environment data, string type = "sequence_data",
                       string bin_type = "otu") {

     Rcpp::XPtr<Dataset> d = data["data"];

     // sequence_data reports contain the starts, ends, ambigs,...
     if (type == "sequence_data") {
        return d.get()->getSequenceReport();
     }
     // sequence classification report
     else if (type == "sequence_taxonomy") {
        return d.get()->getSequenceTaxonomyReport();
     }
     // bin classification report
     else if (type == "bin_taxonomy") {
        return d.get()->getBinTaxonomyReport(bin_type);
     }
     // sequence_scrap report
     else if (type == "sequence_scrap") {
         return d.get()->getScrapReport("sequence");
     }
     // bin_scrap report
     else if (type == "bin_scrap") {
         return d.get()->getScrapReport(bin_type);
     }
     // metadata
     else if (type == "metadata") {
         return d.get()->getMetadata();
     }
     // references
     else if (type == "references") {
         return d.get()->getReferences();
     }
     else {
         // custom reports like alignreport, contigs report and chimera reports
        return d.get()->getReports(type);
     }

     // empty report
     return Rcpp::DataFrame::create();
}
/******************************************************************************/
//' @title get_sample_treatment_assignments
//' @description
//' Get treatment assignments for samples in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' get_sample_treatment_assignments(data)
//'
//' @return 2 column data.frame containing the sample treatment assignments in a
//'  \link{dataset} object
//[[Rcpp::export]]
Rcpp::DataFrame get_sample_treatment_assignments(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getSampleTreatmentAssignments();
}
/******************************************************************************/
//' @title get_sequences
//' @description
//' Get the nucleotide strings for each sequence in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param sample, a string containing the name of the sample you
//' would like sequence names for. For all samples in dataset, sample = "".
//' @examples
//'
//'  data <- miseq_sop_example()
//'  get_sequences(data)
//'
//' @return vector of string containing nucleotide strings of the sequences in
//' a \link{dataset} object
//[[Rcpp::export]]
vector<string> get_sequences(Rcpp::Environment data, string sample = "") {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getSequences(sample);
}
/******************************************************************************/
//' @title get_sequence_summary
//' @description
//' Get a summary of the sequence report data, as well as reports of containing
//' scrapped data.
//'
//' @param data, a \link{dataset} object
//' @examples
//'
//'  data <- miseq_sop_example()
//'
//'  # Sequence summary, nothing has been scrapped
//'
//'  get_sequence_summary(data)
//'
//'  # Sequence summary, after removing sample 'F3D0'
//'
//'  xdev_remove_samples(data, c("F3D0"))
//'  get_sequence_summary(data)
//'
//' @return list of data.frames containing the 'sequence_summary' table and
//' 'scrap_summary' table if sequences have been removed
//[[Rcpp::export]]
Rcpp::List get_sequence_summary(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequenceSummary();
}
/******************************************************************************/
//' @title has_sample
//' @description
//' Determine if a given sample is in a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//' @param sample a string containing the name of a sample.
//' @examples
//'
//' data <- miseq_sop_example()
//' has_sample(data, "F3D0")
//' has_sample(data, "not a valid sample")
//'
//' @return boolean indicating whether the dataset has a given sample
//[[Rcpp::export]]
bool has_sample(Rcpp::Environment data, string sample) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->hasSample(sample);
}
/******************************************************************************/
//' @title has_sequence_strings
//' @description
//' Determine if a \link{dataset} object contains sequence nucleotide strings.
//'
//' @param data, a \link{dataset} object.
//' @examples
//'
//' data <- miseq_sop_example()
//' has_sequence_strings(data)
//'
//' @return boolean indicating whether the dataset has sequence nucleotide
//' strings.
//[[Rcpp::export]]
bool has_sequence_strings(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
   return d.get()->hasSeqs();
}
/******************************************************************************/
//' @title is_aligned
//' @description
//' Determine if a \link{dataset} object contains aligned sequences.
//'
//' @param data, a \link{dataset} object
//' @examples
//'
//' dataset <- miseq_sop_example()
//' is_aligned(dataset)
//'
//' @return Boolean
//[[Rcpp::export]]
bool is_aligned(Rcpp::Environment data) {
     Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->isAligned;
}
/******************************************************************************/
#endif
