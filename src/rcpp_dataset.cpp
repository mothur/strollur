#include <Rcpp.h>
#include "../inst/include/rdataset.h"
#include "dataset.h"



/******************************************************************************/
//' @title new_dataset
//' @description
//' Create a pointer to an instance of the 'Dataset' c++ class.
//'
//' The 'Dataset' class is the c++ implementation of the R6 'sequence_data'
//' object. This class allows package developers access to additional
//' functionality. 'Dataset' stores nucleotide sequences, abundance, sample and
//' treatment assignments, taxonomic classifications, asv / otu clusters. It
//' creates various reports and summaries. It is designed to facilitate data
//' transfer and access across multiple R packages.
//'
//' @param dataset_name string, a string containing the dataset name
//' @param processors integer, number of cores to use during summary functions.
//' @examples
//'
//' # to create a new dataset and allow for 4 processors during summary,
//' # run the following:
//'
//' dataset <- new_dataset("soil", 4)
//'
//' @returns Rcpp::XPtr<Dataset> pointer to an instance of the 'Dataset' c++
//'  class.
//' @seealso [sequence_data$new()]
//[[Rcpp::export]]
Rcpp::XPtr<Dataset> new_dataset(string dataset_name, int processors) {
    Dataset* data = new Dataset(dataset_name, processors);
    return Rcpp::XPtr<Dataset>(data);
}
/******************************************************************************/
//' @title copy_dataset
//' @description
//' Create a pointer to an instance of the 'Dataset' c++ class.
//'
//' The 'Dataset' class is the c++ implementation of the R6 'sequence_data'
//' object. This class allows package developers access to additional
//' functionality. 'Dataset' stores nucleotide sequences, abundance, sample and
//' treatment assignments, taxonomic classifications, asv / otu clusters. It
//' creates various reports and summaries. It is designed to facilitate data
//' transfer and access across multiple R packages.
//'
//' @param dataset an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' # to create a new dataset and allow for 4 processors during summary,
//' # run the following:
//'
//' dataset1 <- new_dataset("soil", 4)
//'
//' # to create a new dataset that is a copy of dataset1
//'
//' dataset2 <- copy_dataset(dataset1)
//'
//' @returns Rcpp::XPtr<Dataset> pointer to an instance of the 'Dataset' c++
//'  class.
//' @seealso [sequence_data$new()]
//[[Rcpp::export]]
Rcpp::XPtr<Dataset> copy_dataset(Rcpp::XPtr<Dataset> dataset) {
    Dataset* data = new Dataset(*(dataset.get()));
    return Rcpp::XPtr<Dataset>(data);
}
/******************************************************************************/
//' @title add_sequences
//' @description
//' Add sequence data to an instance of the 'Dataset' class
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sequence_names a vector of strings containing sequence names
//' @param sequences a vector of strings containing sequence nucleotide strings
//' @param comments a vector of strings containing sequence comments
//' @examples
//'
//'  dataset <- new_dataset("miseq_sop", 4)
//'  sequences <- read_fasta(rdataset_example("final.fasta"))
//'  add_sequences(dataset, sequences$sequence_names, sequences$sequences, "")
//'
//' @seealso [sequence_data$add_sequences()]
//' @return double containing the number of sequences added
//[[Rcpp::export]]
double add_sequences(Rcpp::XPtr<Dataset> data,
                   const vector<string> sequence_names,
                   vector<string> sequences,
                   vector<string> comments) {

    double numSeqsAdded = 0;

    // all three parameters are given
    if ((sequence_names.size() == sequences.size()) &&
        (sequence_names.size() == comments.size())) {
        data.get()->addSequences(sequence_names, sequences, comments);
    }else  {
        if (sequence_names.size() != 1) {
            // no comments
            if (comments.size() == 1) {
                comments = nullVector;
            }
            // no sequences
            if (sequences.size() == 1) {
                sequences = nullVector;
            }
        }
        numSeqsAdded = data.get()->addSequences(sequence_names,
                                sequences, comments);
    }

    return numSeqsAdded;
}
/******************************************************************************/
//' @title assign_bins
//' @description
//' Add bin assignments to an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_names a vector strings containing of bin labels
//' @param abundances a vector of integers containing abundances. Note: You must
//'  provide either abundances or seq_ids.
//' @param samples a vector of strings containing sample assignments
//' @param sequence_names a vector of strings containing sequence names.
//' Note: You must provide either abundances or seq_ids.
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   # To assign sequences to bins:
//'
//'   dataset <- new_dataset("miseq_sop", 4)
//'
//'   otu_data <- read_mothur_list(rdataset_example(
//'                             "final.opti_mcc.list"))
//'   assign_bins(dataset, otu_data$bin_names, 0, "", otu_data$sequence_names)
//'
//'   # To add abundance only bin assignments:
//'
//'   dataset <- new_dataset("miseq_sop", 4)
//'
//'   otu_data <- read_mothur_rabund(rdataset_example(
//'                             "final.opti_mcc.rabund"))
//'   assign_bins(dataset, otu_data$bin_names, otu_data$abundances, "", "")
//'
//'   # To add abundance bin assignments parsed by sample:
//'
//'   dataset <- new_dataset("miseq_sop", 4)
//'   bin_table <- readr::read_tsv(rdataset_example(
//'                                "mothur2_bin_assignments_shared.tsv"))
//'
//'   assign_bins(dataset, bin_table$bin_names, bin_table$abundances,
//'                  bin_table$samples, "")
//'
//'   # To assign sequences to bins with their abundances parsed by sample:
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin1", "bin1", "bin1",
//'                "bin2", "bin2", "bin2",
//'                "bin3", "bin3")
//'   seq_ids <- c("seq1", "seq1", "seq1", "seq2", "seq4", "seq4",
//'                "seq3", "seq3", "seq6",
//'                "seq5", "seq5")
//'   samples <- c("sample1", "sample2", "sample5",
//'                "sample1", "sample3", "sample4",
//'                "sample2", "sample3", "sample1",
//'                "sample1", "sample6")
//'   abundances <- c(10, 100, 1, 500, 25, 80, 20, 5, 60, 15, 50)
//'   assign_bins(dataset, bin_ids, abundances, samples, seq_ids, "otu")
//'
//' @seealso [sequence_data$assign_bins()]
//' @return double containing the number of bins assigned
//[[Rcpp::export]]
double assign_bins(Rcpp::XPtr<Dataset> data, const vector<string> bin_names,
                   vector<int> abundances, vector<string> samples,
                   vector<string> sequence_names, string type = "otu") {

    double numBinsAssigned = 0;

    set<int> lengths;
    lengths.insert(bin_names.size());

    // empty values
    if (bin_names.size() != 1) {
        // no seqIDs
        if (sequence_names.size() == 1) {
            sequence_names = nullVector;
        }
        // no samples
        if (samples.size() == 1) {
            samples = nullVector;
        }
        // no abundances
        if (abundances.size() == 1) {
            abundances = nullIntVector;
        }
    }

    // make sure vector lengths all match
    lengths.insert(sequence_names.size());
    lengths.insert(abundances.size());
    lengths.insert(samples.size());
    lengths.erase(0);

    if (lengths.size() != 1) {
        string message = "[ERROR]: assign_bins expect lengths of inputs";
        message += " to match.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    if ((abundances.size() == 0) && (sequence_names.size() == 0)) {
        string message = "[ERROR]: You must provide either abundances or ";
        message += "sequence_names to assign bins.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    // if you have list assignments, don't allow setting bin abundances
    if (data.get()->hasListAssignments(type) && ((abundances != nullIntVector) ||
        (samples != nullVector))) {
        string message = "[ERROR]: You cannot assign abundance and sample data";
        message += " to bins that have sequence assignments. This could cause ";
        message += "inconsistencies.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }else{
        numBinsAssigned = data.get()->assignBins(bin_names, abundances, samples,
                 sequence_names, type);
    }
    return numBinsAssigned;
}
/******************************************************************************/
//' @title assign_bin_taxonomy
//' @description
//' Assign bin classifications to an instance of the 'Dataset' class.
//'
//' Note, if you assign sequence taxonomies and assign bins, 'Dataset' will find
//'  the concensus taxonomy for each bin for you.
//'
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_names a vector of strings containing bin names
//' @param taxonomies a vector of strings containing bin classifications
//' @param type a string indicating the type of clusters. Default = "otu"
//' @examples
//'
//' otu_data <- read_mothur_cons_taxonomy(rdataset_example(
//'                         "final.cons.taxonomy"))
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_bins(dataset, otu_data$bin_names, otu_data$abundances, "", "")
//' assign_bin_taxonomy(dataset, otu_data$bin_names, otu_data$taxonomies)
//'
//' @seealso [sequence_data$assign_bin_taxonomy()]
//' @return double containing the number of bins assigned
//[[Rcpp::export]]
double assign_bin_taxonomy(Rcpp::XPtr<Dataset> data, vector<string>& bin_names,
                         vector<string>& taxonomies, string type = "otu") {

    if (data.get()->getNumBins(type) == 0) {
        string message = "[ERROR]: No bin data for type " + type + ", please ";
        message += " assign bins using the 'assign_bins' function then try ";
        message += "again.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    return data.get()->assignBinTaxonomy(bin_names, taxonomies, type);
}
/******************************************************************************/
//' @title assign_sequence_abundance
//' @description
//' Set sequence abundance and optionally assign sample and treatment data to
//'  an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sequence_names a vector of strings containing sequence names
//' @param abundances a vector of integers containing sequence abundances
//' @param samples a vector of strings containing sample assignments
//' @param treatments a vector of strings containing treatment assignments
//' @examples
//'
//' sequence_names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
//'              "seq2", "seq3", "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4",
//'             "sample2", "sample3", "sample4",
//'            "sample2", "sample3",
//'            "sample4")
//' abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, sequence_names, abundances, samples, "")
//'
//' @seealso [sequence_data$assign_sequence_abundance()]
//' @return double containing the number of sequences assigned
//[[Rcpp::export]]
double assign_sequence_abundance(Rcpp::XPtr<Dataset> data,
                               vector<string>& sequence_names,
                               vector<int>& abundances,
                               vector<string>& samples,
                               vector<string>& treatments) {

    if (sequence_names.size() != abundances.size()) {
        string message = "[ERROR]: The names and abundances must be the same";
        message += " length.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    vector<string> unique_names = unique(sequence_names);
    vector<string> dataset_names = data.get()->getSequenceNames();

    // add seqs if needed
    if (dataset_names.size() == 0) {
        data.get()->addSequences(unique_names);
    }else {
        // sanity check, make sure names are present in dataset
        if (!identical(unique_names, dataset_names)) {
            string message = "[ERROR]: You must provide assignments for all";
            message += " sequences in your dataset.";
            RcppThread::Rcerr << endl << message << endl;
            throw Rcpp::exception(message.c_str());
        }
    }

    if (sequence_names.size() != 1) {
        // no samples
        if (samples.size() == 1) {
            samples = nullVector;
        }

        // no treatments
        if (treatments.size() == 1) {
            treatments = nullVector;
        }
    }

    return data.get()->assignSequenceAbundance(sequence_names, abundances,
             samples, treatments);
}
/******************************************************************************/
//' @title assign_sequence_taxonomy
//' @description
//' Assign sequence classifications to an instance of the 'Dataset' class.
//'
//' Note, if you assign sequence taxonomies and assign bins, 'Dataset' will find
//'  the concensus taxonomy for each bin for you.
//'
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sequence_names a vector of strings containing sequence names
//' @param taxonomies a vector of strings containing sequence classifications
//' @examples
//'
//' sequence_names <- c("seq1", "seq2", "seq3", "seq4")
//' taxonomies <- c("Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
//'               "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
//'                "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
//'            "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;")
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_taxonomy(dataset, sequence_names, taxonomies)
//'
//' @seealso [sequence_data$assign_sequence_taxonomy()]
//' @return double containing the number of sequences assigned
//[[Rcpp::export]]
double assign_sequence_taxonomy(Rcpp::XPtr<Dataset> data,
                              vector<string>& sequence_names,
                              vector<string>& taxonomies) {

    // make sure names is same size as taxonomies
    if (sequence_names.size() != taxonomies.size()) {
        string message = "[ERROR]: Size mismatch. names and taxonomies must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    return data.get()->assignSequenceTaxonomy(sequence_names, taxonomies);
}
/******************************************************************************/
//' @title assign_treatments
//' @description
//' Assign samples to treatments in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param samples a vector of strings containing sample names
//' @param treatments a vector of strings containing treatment names
//' @examples
//'
//' names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
//'              "seq2", "seq3", "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4",
//'             "sample2", "sample3", "sample4",
//'            "sample2", "sample3",
//'            "sample4")
//' abundances <- c(250, 400, 500,
//'                25, 40, 50,
//'                25, 25,
//'                4)
//' treatments <- c("early", "early", "late")
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, names, abundances, samples, "")
//' assign_treatments(dataset, unique(samples), treatments)
//'
//' @seealso [sequence_data$assign_treatments()]
//' @return double containing the number of bins assigned
//[[Rcpp::export]]
double assign_treatments(Rcpp::XPtr<Dataset> data, vector<string>& samples,
                       vector<string>& treatments) {

    // check to make sure samples and treatments are same length
    if (samples.size() != treatments.size()) {
        string message = "[ERROR]: The samples and treatments must be the same";
        message += " length.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    if (data.get()->getNumSamples() == 0) {
        string message = "[ERROR]: You cannot assign treatments, your dataset";
        message += " does not include sample data.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    // make sure every sample in dataset is assigned a treatment
    if (!identical(data.get()->getSamples(), unique(samples))) {
        string message = "[ERROR]: You must provide treatment assignments for";
        message += " all samples in your dataset.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    return data.get()->assignTreatments(samples, treatments);
}
/******************************************************************************/
//' @title clear
//' @description
//' Clear all data from an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param tags a vector of strings containing the items you wish to clear.
//' Options are 'sequence_data' and 'bin_data'. By default, everything is
//' cleared.
//'
//' @examples
//' dataset <- miseq_sop_example()
//' clear(dataset$data, "")
//'
//[[Rcpp::export]]
void clear(Rcpp::XPtr<Dataset> data, vector<string> tags) {

    if (tags.size() == 1) {
        if (tags[0] == "") {
            tags = nullVector;
        }
    }

    return data.get()->clear(tags);
}
/******************************************************************************/
//' @title export_dataset
//' @description
//' Export all data from an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param tags a vector of strings containing the items you wish to export.
//' Options are 'sequence_data' and 'bin_data'. By default, everything is
//'  exported.
//'
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 2)
//' export_dataset(dataset, "")
//'
//' @return Rcpp::List, containing the data in the 'Dataset
//[[Rcpp::export]]
Rcpp::List export_dataset(Rcpp::XPtr<Dataset> data, vector<string> tags) {

    if (tags.size() == 1) {
        if (tags[0] == "") {
            tags = nullVector;
        }
    }

     return data.get()->exportDataset(tags);
}
/******************************************************************************/
//' @title get_bin
//' @description
//' Get the names of the sequences in a given bin in an instance of the
//' 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_name, string containing the bin name
//' @param type a string indicating the type of clusters. Default = "otu".
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
//' dataset <- new_dataset("my_dataset", 4)
//' assign_bins(dataset, bin_ids, 0, "", seq_ids)
//' get_bin(dataset, "bin1")
//'
//' @return String, containing names of the sequences in a given bin
//[[Rcpp::export]]
string get_bin(Rcpp::XPtr<Dataset> data, string bin_name, string type = "otu") {
    return data.get()->getBin(bin_name, type);
}
/******************************************************************************/
//' @title get_bin_abundance
//' @description
//' Get the abundance of a given bin in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_name, string containing the bin name
//' @param type a string indicating the type of clusters. Default = "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin1", "bin1", "bin1",
//'                "bin2", "bin2", "bin2",
//'                "bin3", "bin3")
//'   seq_ids <- c("seq1", "seq1", "seq1", "seq2", "seq4", "seq4",
//'                "seq3", "seq3", "seq6",
//'                "seq5", "seq5")
//'   samples <- c("sample1", "sample2", "sample5",
//'                "sample1", "sample3", "sample4",
//'                "sample2", "sample3", "sample1",
//'                "sample1", "sample6")
//'   abundances <- c(10, 100, 1, 500, 25, 80, 20, 5, 60, 15, 50)
//'   assign_bins(dataset, bin_ids, abundances, samples, seq_ids)
//'   get_bin_abundance(dataset, "bin1")
//'
//' @return Integer, containing the abundance of a given bin
//[[Rcpp::export]]
int get_bin_abundance(Rcpp::XPtr<Dataset> data,
                       string bin_name, string type = "otu") {
     return data.get()->getBinAbundance(bin_name, type);
}
/******************************************************************************/
//' @title get_bin_abundances
//' @description
//' Get the abundance of a given bin parsed by sample in an instance of the
//' 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_name, string containing the bin name
//' @param type, string indicating the type of clusters. Default = "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin1", "bin1", "bin1",
//'                "bin2", "bin2", "bin2",
//'                "bin3", "bin3")
//'   seq_ids <- c("seq1", "seq1", "seq1", "seq2", "seq4", "seq4",
//'                "seq3", "seq3", "seq6",
//'                "seq5", "seq5")
//'   samples <- c("sample1", "sample2", "sample5",
//'                "sample1", "sample3", "sample4",
//'                "sample2", "sample3", "sample1",
//'                "sample1", "sample6")
//'   abundances <- c(10, 100, 1, 500, 25, 80, 20, 5, 60, 15, 50)
//'
//'   assign_bins(dataset, bin_ids, abundances, samples, seq_ids)
//'   get_bin_abundances(dataset, "bin1")
//'
//' @return vector of integers, containing the abundance of a given bin parsed
//' by sample
//[[Rcpp::export]]
vector<int> get_bin_abundances(Rcpp::XPtr<Dataset> data,
                       string bin_name, string type = "otu") {
    return data.get()->getBinAbundances(bin_name, type);
}
/******************************************************************************/
//' @title get_bin_taxonomy_report
//' @description
//' Get the bin classifications of an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type, string indicating the type of bin clusters. Default = "otu".
//' @examples
//'
//' bin_ids <- c("bin1", "bin2", "bin3", "bin4")
//' abunds <- c(200, 40, 100, 5)
//' taxonomies <- c("Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
//'               "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
//'                "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
//'            "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;")
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_bins(dataset, bin_ids, abunds, "", "")
//' assign_bin_taxonomy(dataset, bin_ids, taxonomies)
//' get_bin_taxonomy_report(dataset)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_bin_taxonomy_report(Rcpp::XPtr<Dataset> data,
                                        string type = "otu") {
    return data.get()->getBinTaxonomyReport(type);
}
/******************************************************************************/
//' @title get_dataset_name
//' @description
//' Get the name of an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//' dataset <- new_dataset("my_dataset", 4)
//' get_dataset_name(dataset)
//'
//' @return String, containing the name of the dataset
//[[Rcpp::export]]
string get_dataset_name(Rcpp::XPtr<Dataset> data) {
    return data.get()->datasetName;
}
/******************************************************************************/
//' @title get_list
//' @description
//' Get data frame containing sequence bin assignments
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
//'   sequence_abundances <- c(10, 100, 1, 500, 25, 80)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   assign_bins(dataset, bin_ids, sequence_abundances, "", seq_ids)
//'
//'   # (list) bins would look like:
//'   # bin1             bin2        bin3
//'   # seq1,seq2,seq4   seq3,seq6   seq5
//'
//'   get_list(dataset)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_list(Rcpp::XPtr<Dataset> data, string type = "otu") {
     return data.get()->getList(type);
}
/******************************************************************************/
//' @title get_list_vector
//' @description
//' Get vector of strings containing the sequences bin data
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
//'   sequence_abundances <- c(10, 100, 1, 500, 25, 80)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   assign_bins(dataset, bin_ids, sequence_abundances, "", seq_ids)
//'
//'   # (list) bins would look like:
//'   # bin1             bin2        bin3
//'   # seq1,seq2,seq4   seq3,seq6   seq5
//'
//'   get_list_vector(dataset)
//'
//' @return vector of strings containing the sequences in each bin separated
//' by commas
//[[Rcpp::export]]
vector<string> get_list_vector(Rcpp::XPtr<Dataset> data, string type = "otu") {
     return data.get()->getListVector(type);
}
/******************************************************************************/
//' @title get_num_processors
//' @description
//' Get the number of processors used to summarize an instance of the
//'  'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//' dataset <- new_dataset("my_dataset", 4)
//' get_num_processors(dataset)
//'
//' @return Integer, containing number of processors
//[[Rcpp::export]]
int get_num_processors(Rcpp::XPtr<Dataset> data) {
     return data.get()->processors;
}
/******************************************************************************/
//' @title get_num_bins
//' @description
//' Get the number of bins of a specific type in an instance of the 'Dataset'
//'  class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type a string indicating the type of clusters. Default = "otu".
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin2", "bin3")
//' abundances <- c(110, 525, 80)
//' assign_bins(dataset, bin_ids, abundances, "", "")
//' get_num_bins(dataset)
//'
//' @return Integer, the number of bins of a specific type in an instance of
//' the 'Dataset' class.
//[[Rcpp::export]]
int get_num_bins(Rcpp::XPtr<Dataset> data, string type = "otu") {
    return data.get()->getNumBins(type);
}
/******************************************************************************/
//' @title get_num_samples
//' @description
//' Get the number of samples in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' samples <- c("sample1", "sample2", "sample5", "sample1", "sample3", "sample1")
//' sample_abundances <- c(10, 100, 1, 500, 25, 80)
//' assign_bins(dataset, bin_ids, sample_abundances, samples, "")
//' get_num_samples(dataset)
//'
//' @return Integer, the number of samples in an instance of the 'Dataset' class.
//[[Rcpp::export]]
int get_num_samples(Rcpp::XPtr<Dataset> data) {
     return data.get()->getNumSamples();
}
/******************************************************************************/
//' @title get_num_sequences
//' @description
//' Get the number of sequences in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param distinct Boolean. When distinct is TRUE the number of unique
//' sequence is returned.
//' @param sample, string containing the name of the sample you want number of
//'  sequences for.
//' @return An integer
//[[Rcpp::export]]
long long get_num_sequences(Rcpp::XPtr<Dataset> data, bool distinct = false,
                            string sample = "") {

    if (distinct) {
        return data.get()->getUniqueTotal(sample);
    }

    return data.get()->getTotal(sample);
}
/******************************************************************************/
//' @title get_num_treatments
//' @description
//' Get the number of treatments in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' samples <- c("sample1", "sample2", "sample5", "sample1", "sample3", "sample1")
//' treatments <- c("early", "early", "late", "early", "late", "early")
//' sample_abundances <- c(10, 100, 1, 500, 25, 80)
//' assign_bins(dataset, bin_ids, sample_abundances, samples, treatments)
//' get_num_treatments(dataset)
//'
//' @return Integer, the number of treatments in an instance of the 'Dataset' class.
//[[Rcpp::export]]
int get_num_treatments(Rcpp::XPtr<Dataset> data) {
     return data.get()->getNumTreatments();
}
/******************************************************************************/
//' @title get_rabund
//' @description
//' Get data.frame containing bin abundance data
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
//'   sequence_abundances <- c(10, 100, 1, 500, 25, 80)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   assign_bins(dataset, bin_ids, sequence_abundances, "", seq_ids)
//'
//'   # (rabund) bins would look like:
//'   # bin1  bin2  bin3
//'   # 111   525   80
//'
//'   get_rabund(dataset)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_rabund(Rcpp::XPtr<Dataset> data, string type = "otu") {
    return data.get()->getRAbund(type);
}
/******************************************************************************/
//' @title get_rabund_vector
//' @description
//' Get vector of integers containing bin abundance data
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   seq_ids <- c("seq1", "seq2", "seq4", "seq3", "seq6", "seq5")
//'   sequence_abundances <- c(10, 100, 1, 500, 25, 80)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   assign_bins(dataset, bin_ids, sequence_abundances, "", seq_ids)
//'
//'   # (rabund) bins would look like:
//'   # bin1  bin2  bin3
//'   # 111   525   80
//'
//'   get_rabund_vector(dataset)
//'
//' @return vector of integers containing each bins abundance
//[[Rcpp::export]]
vector<int> get_rabund_vector(Rcpp::XPtr<Dataset> data, string type = "otu") {
    return data.get()->getRAbundVector(type);
}
/******************************************************************************/
//' @title get_samples
//' @description
//' Get the samples in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' samples <- c("sample1", "sample2", "sample5", "sample1", "sample3", "sample1")
//' sample_abundances <- c(10, 100, 1, 500, 25, 80)
//' assign_bins(dataset, bin_ids, sample_abundances, samples, "")
//' get_samples(dataset)
//'
//' @return vector of strings containing the names of the samples in an instance
//'  of the 'Dataset' class.
//[[Rcpp::export]]
vector<string> get_samples(Rcpp::XPtr<Dataset> data) {
     return data.get()->getSamples();
}
/******************************************************************************/
//' @title get_sample_totals
//' @description
//' Get the number of sequences in each sample in an instance of the 'Dataset'
//'  class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' samples <- c("sample1", "sample2", "sample5", "sample1", "sample3", "sample1")
//' sample_abundances <- c(10, 100, 1, 500, 25, 80)
//' assign_bins(dataset, bin_ids, sample_abundances, samples, "")
//' get_sample_totals(dataset)
//'
//' @return vector of integers containing the number of sequences in each
//' sample in an instance of the 'Dataset' class.
//[[Rcpp::export]]
vector<int> get_sample_totals(Rcpp::XPtr<Dataset> data) {
   return data.get()->getSampleTotals();
}
/******************************************************************************/
//' @title get_scrap_report
//' @description
//' Get a scrap report containing sequences and bins eliminated from an instance
//'  of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type a string indicating the type of scrap report you would like.
//'  Default = 'sequence'.
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin1", "bin1", "bin1",
//'                "bin2", "bin2", "bin2", "bin3", "bin3")
//'   seq_ids <- c("seq1", "seq1", "seq1", "seq2", "seq4", "seq4",
//'                "seq3", "seq3", "seq6", "seq5", "seq5")
//'   samples <- c("sample1", "sample2", "sample5", "sample1", "sample3",
//'                "sample4", "sample2", "sample3", "sample1", "sample1",
//'                "sample6")
//'   abundances <- c(10, 100, 1, 500, 25, 80, 20, 5, 60, 15, 50)
//'   assign_bins(dataset, bin_ids, abundances, samples, seq_ids, "otu")
//'
//'   remove_bins(dataset, c("bin1"), c("bad_bin"))
//'
//'   sequence_scrap_report <- get_scrap_report(dataset, "sequence")
//'   otu_scrap_report <- get_scrap_report(dataset, "otu")
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_scrap_report(Rcpp::XPtr<Dataset> data,
                               string type = "sequence") {

    return data.get()->getScrapReport(type);
}
/******************************************************************************/
//' @title get_sequence_abundances
//' @description
//' Get the total abundance for each sequence in an instance of the 'Dataset'
//'  class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
//'              "seq2", "seq3", "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4",
//'             "sample2", "sample3", "sample4",
//'            "sample2", "sample3",
//'            "sample4")
//' abundances <- c(250, 400, 500,
//'                25, 40, 50,
//'                25, 25,
//'                4)
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, names, abundances, samples, "")
//' get_sequence_abundances(dataset)
//'
//' @return vector of integers containing the total abundance for each sequence
//'  in the 'Dataset' class.
//[[Rcpp::export]]
vector<int> get_sequence_abundances(Rcpp::XPtr<Dataset> data) {
     return data.get()->getSequenceAbundances();
}
/******************************************************************************/
//' @title get_sequence_abundances_by_sample
//' @description
//' Get the abundances of each sequence in an instance of the 'Dataset' class
//' parsed by sample.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
//'              "seq2", "seq3", "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4",
//'             "sample2", "sample3", "sample4",
//'            "sample2", "sample3",
//'            "sample4")
//' abundances <- c(250, 400, 500,
//'                25, 40, 50,
//'                25, 25,
//'                4)
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, names, abundances, samples, "")
//' get_sequence_abundances_by_sample(dataset)
//'
//' @return 2D vector of integers ([num_seqs][num_samples]) containing the
//' abundances of each sequence in an instance of the 'Dataset' class parsed by
//' sample.
//[[Rcpp::export]]
vector<vector<int> > get_sequence_abundances_by_sample(Rcpp::XPtr<Dataset> data) {
     return data.get()->getSeqsAbundsBySample();
}
/******************************************************************************/
//' @title get_sequence_abundance_table
//' @description
//' Get the abundances of each sequence in an instance of the 'Dataset' class
//' parsed by sample.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' samples <- c("sample1", "sample2", "sample5", "sample1", "sample3", "sample1")
//' treatments <- c("early", "early", "late", "early", "late", "early")
//' sample_abundances <- c(10, 100, 1, 500, 25, 80)
//' assign_bins(dataset, bin_ids, sample_abundances, samples, treatments)
//'
//' get_sequence_abundance_table(dataset)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_sequence_abundance_table(Rcpp::XPtr<Dataset> data) {
    return data.get()->getSequenceAbundanceTable();
}
/******************************************************************************/
//' @title get_sequence_names
//' @description
//' Get the names of the sequences in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sample a string containing the name of the sample you
//' would like sequence names for. For all samples in dataset, sample = "".
//' @examples
//'
//' sequence_names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
//'              "seq2", "seq3", "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4",
//'             "sample2", "sample3", "sample4",
//'            "sample2", "sample3",
//'            "sample4")
//' abundances <- c(250, 400, 500,
//'                25, 40, 50,
//'                25, 25,
//'                4)
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, sequence_names, abundances, samples, "")
//' get_sequence_names(dataset)
//'
//' @return vector of string containing the names of the sequences in the
//'  'Dataset' class.
//[[Rcpp::export]]
vector<string> get_sequence_names(Rcpp::XPtr<Dataset> data, string sample = "") {
     return data.get()->getSequenceNames(sample);
}
/******************************************************************************/
//' @title get_sequence_names_by_sample
//' @description
//' Get the names of the sequences in an instance of the 'Dataset' class parsed
//' by sample.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param samples a vector of strings containing the names of the samples you
//' would like sequence names for. For all samples in dataset, samples = "".
//' @examples
//'
//' names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
//'              "seq2", "seq3", "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4",
//'             "sample2", "sample3", "sample4",
//'            "sample2", "sample3",
//'            "sample4")
//' abundances <- c(250, 400, 500,
//'                25, 40, 50,
//'                25, 25,
//'                4)
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, names, abundances, samples, "")
//' get_sequence_names_by_sample(dataset, "")
//'
//' @return 2D vector of strings ([num_seqs][num_samples]) containing the
//' names of each sequence in an instance of the 'Dataset' class parsed by
//' sample.
//[[Rcpp::export]]
vector<vector<string> > get_sequence_names_by_sample(
        Rcpp::XPtr<Dataset> data, vector<string> samples) {

    // no samples given
    if (samples.size() == 1) {
        if (samples[0] == "") {
            samples = nullVector;
        }
    }
    return data.get()->getSequenceNamesBySample(samples);
}
/******************************************************************************/
//' @title get_sequences
//' @description
//' Get the nucleotide strings for each sequence in an instance of the 'Dataset'
//'  class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sample a string containing the name of the sample you
//' would like sequence names for. For all samples in dataset, sample = "".
//' @examples
//'
//'  dataset <- new_dataset("miseq_sop", 4)
//'  sequences <- read_fasta(rdataset_example("final.fasta"))
//'  add_sequences(dataset, sequences$sequence_names, sequences$sequences, "")
//'  get_sequences(dataset)
//'
//' @return vector of string containing nucleotide strings of the sequences in
//' the 'Dataset' class.
//[[Rcpp::export]]
vector<string> get_sequences(Rcpp::XPtr<Dataset> data, string sample = "") {
     return data.get()->getSequences(sample);
}
/******************************************************************************/
//' @title get_sequences_by_sample
//' @description
//' Get the nucleotide strings for each sequence in an instance of the 'Dataset'
//' class parsed by sample.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param samples a vector of strings containing the names of the samples you
//' would like sequence names for. For all samples in dataset, samples = "".
//' @examples
//'
//' seq_names <- c("seq1", "seq2", "seq3", "seq4")
//' sequences <- c("ATTGC", "ACTGC", "AGTGC", "TTTGC")
//'
//' names <- c("seq1", "seq1", "seq2", "seq2",
//'             "seq3", "seq4", "seq4")
//' abundances <- c(200, 50, 400, 10, 25, 100,  425)
//' samples <- c("sample1", "sample2", "sample1", "sample2",
//'              "sample1", "sample1", "sample2")
//'
//' dataset <- new_dataset("my_dataset", 4)
//' add_sequences(dataset, seq_names, sequences, "")
//' assign_sequence_abundance(dataset, names, abundances, samples, "")
//' get_sequences_by_sample(dataset, "")
//'
//' @return 2D vector of strings ([num_seqs][num_samples]) containing the
//' nucleotide strings for each sequence in an instance of the 'Dataset' class
//' parsed by sample.
//[[Rcpp::export]]
vector<vector<string> > get_sequences_by_sample(
         Rcpp::XPtr<Dataset> data, vector<string> samples) {

     // no samples given
     if (samples.size() == 1) {
         if (samples[0] == "") {
             samples = nullVector;
         }
     }
     return data.get()->getSequencesBySample(samples);
}
/******************************************************************************/
//' @title get_sequence_report
//' @description
//' Get sequence report data: starts, ends, lengths, ambigs, longest
//' homopolymers and numns.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//'  dataset <- new_dataset("miseq_sop", 4)
//'  sequences <- read_fasta(rdataset_example("final.fasta"))
//'  add_sequences(dataset, sequences$sequence_names, sequences$sequences, "")
//'  get_sequence_report(dataset)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_sequence_report(Rcpp::XPtr<Dataset> data) {
    return data.get()->getSequenceReport();
}
/******************************************************************************/
//' @title get_sequence_summary
//' @description
//' Get a summary of the sequence report data, as well as reports of containing
//' scrapped data.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//'  dataset <- new_dataset("miseq_sop", 4)
//'  sequences <- read_fasta(rdataset_example("final.fasta"))
//'  add_sequences(dataset, sequences$sequence_names, sequences$sequences, "")
//'  get_sequence_summary(dataset)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::List get_sequence_summary(Rcpp::XPtr<Dataset> data) {
    return data.get()->getSequenceSummary();
}
/******************************************************************************/
//' @title get_sequence_taxonomy_report
//' @description
//' Get the sequence classifications of an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' seq_ids <- c("seq1", "seq2", "seq3", "seq4")
//' taxonomies <- c("Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;",
//'               "Bacteria;Proteobacteria;Betaproteobacteria;Neisseriales;",
//'                "Bacteria;Firmicutes;Bacilli;Lactobacillales;",
//'            "Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;")
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_taxonomy(dataset, seq_ids, taxonomies)
//' get_sequence_taxonomy_report(dataset)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_sequence_taxonomy_report(Rcpp::XPtr<Dataset> data) {
    return data.get()->getSequenceTaxonomyReport();
}
/******************************************************************************/
//' @title get_bin_assignments
//' @description
//' Get data.frame containing bin abundance data
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   samples <- c("sample1", "sample2", "sample5",
//'    "sample1", "sample3", "sample1")
//'   sample_abundances <- c(10, 100, 1, 500, 25, 80)
//'   assign_bins(dataset, bin_ids, sample_abundances, samples, "")
//'
//'   # (shared) bins would look like:
//'   # sample   bin1   bin2   bin3
//'   # sample1  10     500    80
//'   # sample2  100    0      0
//'   # sample3  0      25     0
//'   # sample5  1      0      0
//'
//'   get_bin_assignments(dataset)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_bin_assignments(Rcpp::XPtr<Dataset> data, string type = "otu") {
     return data.get()->getShared(type);
}
/******************************************************************************/
//' @title get_shared_vector
//' @description
//' Get 2D vector of integers containing bin abundance data by sample
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   samples <- c("sample1", "sample2", "sample5",
//'    "sample1", "sample3", "sample1")
//'   sample_abundances <- c(10, 100, 1, 500, 25, 80)
//'   assign_bins(dataset, bin_ids, sample_abundances, samples, "")
//'
//'   # (shared) bins would look like:
//'   # sample   bin1   bin2   bin3
//'   # sample1  10     500    80
//'   # sample2  100    0      0
//'   # sample3  0      25     0
//'   # sample5  1      0      0
//'
//'   get_shared_vector(dataset)
//'
//' @return 2D vector of integers ([num_bins][num_samples]) containing
//' the abundances of each bin parsed by sample.
//[[Rcpp::export]]
vector<vector<int> > get_shared_vector(Rcpp::XPtr<Dataset> data, string type = "otu") {
     return data.get()->getSharedVector(type);
}
/******************************************************************************/
//' @title get_treatments
//' @description
//' Get the treatments in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' samples <- c("sample1", "sample2", "sample5", "sample1", "sample3", "sample1")
//' treatments <- c("early", "early", "late", "early", "late", "early")
//' sample_abundances <- c(10, 100, 1, 500, 25, 80)
//' assign_bins(dataset, bin_ids, sample_abundances, samples, treatments)
//' get_treatments(dataset)
//'
//' @return vector of strings containing the names of the treatments in an
//' instance of the 'Dataset' class.
//[[Rcpp::export]]
vector<string> get_treatments(Rcpp::XPtr<Dataset> data) {
     return data.get()->getTreatments();
}
/******************************************************************************/
//' @title get_treatment_totals
//' @description
//' Get the number of sequences in each treatment in an instance of the 'Dataset'
//'  class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' samples <- c("sample1", "sample2", "sample5", "sample1", "sample3", "sample1")
//' treatments <- c("early", "early", "late", "early", "late", "early")
//' sample_abundances <- c(10, 100, 1, 500, 25, 80)
//' assign_bins(dataset, bin_ids, sample_abundances, samples, treatments)
//' get_treatment_totals(dataset)
//'
//' @return vector of integers containing the number of sequences in each
//' treatment in an instance of the 'Dataset' class.
//[[Rcpp::export]]
vector<int> get_treatment_totals(Rcpp::XPtr<Dataset> data) {
     return data.get()->getTreatmentTotals();
}
/******************************************************************************/
//' @title has_sample
//' @description
//' Determine if a given sample is in an instance of the 'Dataset'
//'  class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sample a string containing the name of a sample.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//' samples <- c("sample1", "sample2", "sample5", "sample1", "sample3", "sample1")
//' sample_abundances <- c(10, 100, 1, 500, 25, 80)
//' treatments <- c("early", "early", "late", "early", "late", "early")
//' assign_bins(dataset, bin_ids, sample_abundances, samples, treatments)
//' has_sample(dataset, "sample2")
//'
//' @return boolean indicating whether the dataset has a given sample
//[[Rcpp::export]]
bool has_sample(Rcpp::XPtr<Dataset> data, string sample) {
     return data.get()->hasSample(sample);
}
/******************************************************************************/
//' @title has_sequence_strings
//' @description
//' Determine if an instance of the 'Dataset' class contains sequence
//' nucleotide strings.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' seq_names <- c("seq1", "seq2", "seq3", "seq4")
//' sequences <- c("ATTGC", "ACTGC", "AGTGC", "TTTGC")
//'
//' dataset <- new_dataset("my_dataset", 4)
//' add_sequences(dataset, seq_names, sequences, "")
//' has_sequence_strings(dataset)
//'
//' @return boolean indicating whether the dataset has sequence nucleotide
//' strings.
//[[Rcpp::export]]
bool has_sequence_strings(Rcpp::XPtr<Dataset> data) {
   return data.get()->hasSeqs();
}
/******************************************************************************/
//' @title merge_bins
//' @description
//' Merge bins in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_names a vector of strings containing the names of the bins you
//' would like merge. The resulting merged bin will be stored in the first
//' bin_id in the vector.
//' @param reason a string indicating why you are merging bins
//' @param type a string indicating the type of bin clusters. Default = "otu"
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin2", "bin3")
//'   abundances <- c(110, 525, 80)
//'   assign_bins(dataset, bin_ids, abundances, "", "", "otu")
//'   bins_to_merge <- c("bin1", "bin3")
//'   merge_bins(dataset, bins_to_merge)
//'
//[[Rcpp::export]]
void merge_bins(Rcpp::XPtr<Dataset> data, vector<string> bin_names,
                string reason = "merged", string type = "otu") {
    data.get()->mergeBins(bin_names, reason, type);
}
/******************************************************************************/
//' @title merge_sequences
//' @description
//' Merge sequences in an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sequence_names a vector of strings containing the names of the sequences you
//' sequences you would like merge. The resulting merged sequence will be stored
//' in the first sequence name in the vector.
//' @param reason a string indicating why you are merging sequences.
//' Default = "merged"
//' @examples
//'
//' sequence_names <- c("seq1", "seq1", "seq1", "seq2", "seq2",
//'              "seq2", "seq3", "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4",
//'             "sample2", "sample3", "sample4",
//'            "sample2", "sample3", "sample4")
//' abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, sequence_names, abundances, samples, "")
//'
//' seqs_to_merge <- c("seq1", "seq4")
//' merge_sequences(dataset, seqs_to_merge, "identical")
//'
//[[Rcpp::export]]
void merge_sequences(Rcpp::XPtr<Dataset> data, vector<string> sequence_names,
                 string reason = "merged") {
     data.get()->mergeSequences(sequence_names, reason);
}
/******************************************************************************/
//' @title remove_bins
//' @description
//' Remove bins from an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_names a vector of strings containing the names of the bins you
//' would like removed.
//' @param trash_tags a vector of strings containing the reasons you are
//' removing each bin
//' @param type a string indicating the type of clusters.
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin2", "bin3")
//'   abundances <- c(110, 525, 80)
//'   assign_bins(dataset, bin_ids, abundances, "", "")
//'   bins_to_remove <- c("bin1")
//'   trash_tag <- c("bad_bin")
//'   remove_bins(dataset, bins_to_remove, trash_tag)
//'
//[[Rcpp::export]]
void remove_bins(Rcpp::XPtr<Dataset> data, vector<string> bin_names,
                 vector<string> trash_tags, string type = "otu") {
    data.get()->removeBins(bin_names, trash_tags, type);
}
/******************************************************************************/
//' @title remove_lineages
//' @description
//' Remove contaminants from an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param contaminants vector of strings containing the taxonomies you would
//' like to remove
//' @param trash_tag a string containing reason you are removing the lineages.
//' Default = "contaminant".
//' @examples
//' dataset <- read_mothur(fasta = rdataset_example("final.fasta"),
//'                        count = rdataset_example("final.count_table"),
//'                        taxonomy = rdataset_example("final.taxonomy"),
//'                        design = rdataset_example("mouse.time.design"),
//'                        otu_list = rdataset_example("final.opti_mcc.list"),
//'                        dataset_name = "miseq_sop")
//'
//' contaminants <- c("Chloroplast", "Mitochondria", "unknown", "Archaea",
//'  "Eukaryota")
//'
//' remove_lineages(dataset$data, contaminants)
//'
//[[Rcpp::export]]
void remove_lineages(Rcpp::XPtr<Dataset> data, vector<string> contaminants,
                     string trash_tag = "contaminant") {
    data.get()->removeLineages(contaminants, trash_tag);
}
/******************************************************************************/
//' @title remove_samples
//' @description
//' Remove samples from an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param samples vector of strings containing the names of the samples to
//' @examples
//' dataset <- read_mothur(fasta = rdataset_example("final.fasta"),
//'                        count = rdataset_example("final.count_table"),
//'                        taxonomy = rdataset_example("final.taxonomy"),
//'                        design = rdataset_example("mouse.time.design"),
//'                        otu_list = rdataset_example("final.opti_mcc.list"),
//'                        dataset_name = "miseq_sop")
//'
//' get_num_samples(dataset$data)
//'
//' # To remove samples 'F3D0' and 'F3D1'
//'
//' remove_samples(dataset$data, c("F3D0", "F3D1"))
//'
//' get_num_samples(dataset$data)
//'
//[[Rcpp::export]]
void remove_samples(Rcpp::XPtr<Dataset> data, vector<string> samples) {
    data.get()->removeSamples(samples);
}
/******************************************************************************/
//' @title remove_sequences
//' @description
//' Remove sequences from an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param names vector of strings containing the names of the sequences to remove
//' @param trash_tags vector of strings containing the reasons for the sequences
//'  removals
//' @examples
//'
//' seq_names <- c("seq1", "seq2", "seq3", "seq4")
//'
//' dataset <- new_dataset("my_dataset", 4)
//' add_sequences(dataset, seq_names, "", "")
//'
//' names_to_remove <- c("seq1", "seq3")
//' trash_tags <- c("screen_seqs-too_short", "screen_seqs-contains_ns")
//'
//' remove_sequences(dataset, names_to_remove, trash_tags)
//'
//[[Rcpp::export]]
void remove_sequences(Rcpp::XPtr<Dataset> data,
                      vector<string> names, vector<string> trash_tags) {
    data.get()->removeSequences(names, trash_tags);
}
/******************************************************************************/
//' @title set_abundance
//' @description
//' Set sequence abundances for an instance of the 'Dataset' class without
//' sample data
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sequence_names a vector of strings containing sequence names
//' @param sequence_abundances vector of integers containing the abundances of
//' each sequence.
//' @param reason a string containing the trash tag to be applied to any sequences
//'  set to 0 abundance. Default = "update".
//' @examples
//'
//' sequence_names <- c("seq1", "seq2", "seq3",  "seq4")
//' abundances <- c(1250, 65, 50, 4)
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, sequence_names, abundances, "", "")
//' get_sequence_abundances(dataset)
//'
//' seqs_to_update <- c("seq1", "seq3")
//' new_abunds <- c(1000, 100)
//'
//' set_abundance(dataset, seqs_to_update, new_abunds)
//' get_sequence_abundances(dataset)
//'
//[[Rcpp::export]]
void set_abundance(Rcpp::XPtr<Dataset> data,
                   vector<string> sequence_names, vector<int> sequence_abundances,
                   string reason = "update") {

    if (get_num_samples(data) == 0) {
        data.get()->setAbundance(sequence_names, sequence_abundances, reason);
    }else{
        string message = "[ERROR]: You cannot set the total sequence abundance";
        message += " for sequences whose abundances are parsed by sample. ";
        message += "Try 'set_abundances' instead of 'set_abundance'.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

}
/******************************************************************************/
//' @title set_abundances
//' @description
//' Set sequence abundances for an instance of the 'Dataset' class with
//' sample data
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sequence_names a vector of strings containing sequence names
//' @param abundances 2D vector of integers ([num_seqs][num_samples]) containing
//' the abundances of each sequence parsed by sample.
//' @param reason a string containing the trash tag to be applied to any sequences
//'  set to 0 abundance. Default = "update".
//' @examples
//'
//' sequence_names <- c("seq1", "seq1", "seq1", "seq2", "seq2", "seq2", "seq3",
//'                     "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4", "sample2", "sample3",
//'              "sample4", "sample2", "sample3", "sample4")
//' abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)
//'
//' dataset <- new_dataset("my_dataset", 4)
//' assign_sequence_abundance(dataset, sequence_names, abundances, samples, "")
//'
//' seqs_to_update <- c("seq4")
//' new_abunds <- list(c(20, 10, 4))
//'
//' set_abundances(dataset, seqs_to_update, new_abunds)
//'
//[[Rcpp::export]]
void set_abundances(Rcpp::XPtr<Dataset> data,
                   vector<string> sequence_names, vector<vector<int>> abundances,
                   string reason = "update") {
    if (get_num_samples(data) != 0) {
        data.get()->setAbundances(sequence_names, abundances, reason);
    }else {
        string message = "[ERROR]: You cannot set parsed sequence abundances ";
        message += "when your dataset does not include sample data. ";
        message += "Try 'set_abundance' instead of 'set_abundances'.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }
}
/******************************************************************************/
//' @title set_bin_abundance
//' @description
//' Set bin abundances for an instance of the 'Dataset' class without sample data
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_names a vector strings containing of bin names to set the
//' abundances for.
//' @param abunds vector of integers containing the abundances of each bin.
//' @param reason a string containing the trash tag to be applied to any bins
//'  set to 0 abundance. Default = "update".
//' @param type a string indicating the type of clusters. Default = "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin2", "bin3")
//'   abundances <- c(110, 525, 80)
//'   assign_bins(dataset, bin_ids, abundances, "", "")
//'
//'   bins <- c("bin1", "bin2")
//'   new_abunds <- c(300, 250)
//'
//'   set_bin_abundance(dataset, bins, new_abunds)
//'   get_bin_abundance(dataset, "bin1")
//'
//[[Rcpp::export]]
void set_bin_abundance(Rcpp::XPtr<Dataset> data,
                    vector<string> bin_names, vector<int> abunds,
                    string reason = "update", string type = "otu") {

    if (data.get()->hasListAssignments(type)) {
        string message = "[ERROR]: You cannot set the bin abundance for bin ";
        message += "clusters with sequence assignments.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());

    }else if (get_num_samples(data) != 0) {
        string message = "[ERROR]: You cannot set the total bin abundance";
        message += " for bins whose abundances are parsed by sample. ";
        message += "Try 'set_bin_abundances' instead of 'set_bin_abundance'.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }else{
        data.get()->setBinAbundance(bin_names, abunds, reason, type);
    }
}
/******************************************************************************/
//' @title set_bin_abundances
//' @description
//' Set bin abundances for an instance of the 'Dataset' class with sample data
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param bin_names a vector strings containing of bin names to set the
//' abundances for.
//' @param abunds 2D vector of integers ([num_seqs][num_samples]) containing the
//' abundances of each bin parsed by sample.
//' @param reason a string containing the trash tag to be applied to any bins
//'  set to 0 abundance. Default = "update".
//' @param type a string indicating the type of clusters. Default = "otu".
//' @examples
//'
//'   dataset <- new_dataset("my_dataset", 4)
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   samples <- c("sample1", "sample2", "sample5",
//'    "sample1", "sample3", "sample1")
//'   sample_abundances <- c(10, 100, 1, 500, 25, 80)
//'   assign_bins(dataset, bin_ids, sample_abundances, samples, "")
//'
//'   # bin1's abundances parsed by sample: c(10,100,0,1)
//'   old_bin1_abunds <- get_bin_abundances(dataset, "bin1")
//'
//'   new_bin1_abunds <- list(c(10,50,0,0))
//'   bins <- c("bin1")
//'
//'   set_bin_abundances(dataset, bins, new_bin1_abunds)
//'   get_bin_abundances(dataset, "bin1")
//'
//[[Rcpp::export]]
void set_bin_abundances(Rcpp::XPtr<Dataset> data,
                     vector<string> bin_names, vector<vector<int>> abunds,
                     string reason = "update", string type = "otu") {

    if (data.get()->hasListAssignments(type)) {
        string message = "[ERROR]: You cannot set the bin abundance for bin ";
        message += "clusters with sequence assignments.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }else if (get_num_samples(data) == 0) {
        string message = "[ERROR]: You cannot set parsed bin abundances ";
        message += "when your dataset does not include sample data. ";
        message += "Try 'set_bin_abundance' instead of 'set_bin_abundances'.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }else{
        data.get()->setBinAbundances(bin_names, abunds, reason, type);
    }
}
/******************************************************************************/
//' @title set_sequences
//' @description
//' Set sequence nucleotide strings for an instance of the 'Dataset' class
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param sequence_names a vector of strings containing sequence names
//' @param sequences a vector of strings containing sequence nucleotide strings
//' @param comments a vector of strings containing sequence comments
//' @examples
//'
//' sequence_names <- c("seq1", "seq2", "seq3", "seq4")
//' sequences <- c("ATTGC", "ACTGC", "AGTGC", "TTTGC")
//'
//' dataset <- new_dataset("my_dataset", 4)
//'
//' add_sequences(dataset, sequence_names, "", "")
//' set_sequences(dataset, sequence_names, sequences, "")
//'
//[[Rcpp::export]]
void set_sequences(Rcpp::XPtr<Dataset> data,
                   vector<string> sequence_names, vector<string> sequences,
                   vector<string> comments) {

    if (sequence_names.size() != 1) {
        // no comments
        if (comments.size() == 1) {
            comments = nullVector;
        }
    }

    data.get()->setSequences(sequence_names, sequences, comments);
}
/******************************************************************************/
//' @title set_dataset_name
//' @description
//' Set the name of an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param dataset_name, a string containing the desired name
//' @examples
//' dataset <- new_dataset("my_dataset", 4)
//' set_dataset_name(dataset, "new_dataset_name")
//'
//[[Rcpp::export]]
void set_dataset_name(Rcpp::XPtr<Dataset> data, string dataset_name) {
    data.get()->datasetName = dataset_name;
}
/******************************************************************************/
//' @title set_num_processors
//' @description
//' Set the number of processors used to summarize an instance of the
//'  'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @param processors, a integer containing the desired number of processors
//' @examples
//' dataset <- new_dataset("my_dataset", 4)
//' set_num_processors(dataset, 8)
//'
//[[Rcpp::export]]
void set_num_processors(Rcpp::XPtr<Dataset> data, int processors) {
    data.get()->processors = processors;
}

/******************************************************************************/
//' @title is_aligned
//' @description
//' Determine if the instance of the 'Dataset' class contains aligned sequences.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//' dataset <- new_dataset("my_dataset", 4)
//' is_aligned(dataset)
//'
//' @return Boolean
//[[Rcpp::export]]
bool is_aligned(Rcpp::XPtr<Dataset> data) {
    return data.get()->isAligned;
}
/******************************************************************************/
//' @title load_dataset
//' @description
//' Load an instance of the 'Dataset' class from serialized raw data.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//' 'Dataset' c++ class.
//' @param raw a Rcpp::RawVector containing the serialize 'Dataset' class data.
//' @examples
//' dataset <- new_dataset("my_dataset", 4)
//' raw <- serialize(dataset)
//' load_dataset(dataset, raw)
//'
//[[Rcpp::export]]
void load_dataset(Rcpp::XPtr<Dataset> data, Rcpp::RawVector raw) {
     data.get()->loadFromSerialized(raw);
}
/******************************************************************************/
//' @title serialize
//' @description
//' Serialize an instance of the 'Dataset' class.
//' @param data an Rcpp::XPtr<Dataset> pointer to an instance of the
//'  'Dataset' c++ class.
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 4)
//' raw <- serialize(dataset)
//' load_dataset(dataset, raw)
//'
//' @return Rcpp::RawVector
//[[Rcpp::export]]
Rcpp::RawVector serialize(Rcpp::XPtr<Dataset> data) {
     return data.get()->serializeDataset();
}
/******************************************************************************/


