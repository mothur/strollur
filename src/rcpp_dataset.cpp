#include <Rcpp.h>
#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
SEXP fill_required_parameters(const Rcpp::DataFrame df,
                              const string& given_column_name,
                              string type = "string") {

    if (!df.containsElementNamed(given_column_name.c_str())) {
        string message = "Expected a data.frame column named " +
            given_column_name +" to be provided.";

        throw Rcpp::exception(message.c_str());
    }

    return df[given_column_name];
}
/******************************************************************************/
SEXP fill_optional_parameters(const Rcpp::DataFrame df,
                              const string& default_column_name,
                              const string& given_column_name,
                              string type = "string") {


        // Check if the column exists in the data frame
        if (df.containsElementNamed(given_column_name.c_str())) {
            return df[given_column_name];
        }else if (df.containsElementNamed(default_column_name.c_str())) {
            return df[default_column_name];
        }else if ((given_column_name == "") ||
                        (given_column_name == default_column_name)) {
            //fall through
        }else {
            string message = "Expected a data.frame column named " +
                given_column_name +" to be provided.";

            throw Rcpp::exception(message.c_str());
        }

    if (type == "string") {
        return Rcpp::CharacterVector(0);
    }else if ((type == "float") || (type == "double")) {
        return Rcpp::NumericVector(0);
    }else{
        string message = "Unsupported column type for conversion.";
        throw Rcpp::exception(message.c_str());
    }

    return Rcpp::CharacterVector(0);
}
/******************************************************************************/
//' @title new_dataset
//' @description
//' Create a new \link{dataset} object
//'
//' @param dataset_name string, a string containing the dataset name.
//' Default = ""
//' @param processors integer, number of cores to use during summary functions.
//' Default = 2
//' @examples
//'
//' data <- new_dataset()
//'
//' # to create a new dataset named "soil" and allow for 2 processors during
//' # summary functions, run the following:
//'
//' data <- new_dataset("soil", 2)
//'
//' # to create a new dataset named "soil" and allow for all available
//' # processors during summary functions, run the following:
//'
//' data <- new_dataset("soil", get_available_processors())
//'
//' @returns a \link{dataset} object
//' @seealso The 'new' method in the \link{dataset} class
//[[Rcpp::export]]
Rcpp::Environment new_dataset(string dataset_name = "", int processors = 2) {

    // dataset$new()
    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Environment dataset_class_env = rdataset_env["dataset"];
    Rcpp::Function constructor = dataset_class_env["new"];

    Rcpp::Environment data = constructor(dataset_name, processors, R_NilValue);
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
//' @seealso [add_references()]
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
//' @title add_references
//' @description
//' Add resource references to a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing reference_names, reference_versions
//' (optional), reference_usages (optional), reference_notes (optional), and
//' reference_urls (optional).
//'
//' @param reference_name, a string containing the name of the column in 'table'
//' that contains the reference names. Default column name is 'reference_names'.
//' @param reference_version, a string containing the name of the column in
//' 'table' that contains the reference versions. Default column name is
//' 'reference_versions'.
//' @param reference_usage, a string containing the name of the column in
//' 'table' that contains the reference usages. Default column name is
//'  reference_usages'.
//' @param reference_note, a string containing the name of the column in
//' 'table' that contains the reference notes. Default column name is
//'  reference_notes'.
//' @param reference_url, a string containing the name of the column in
//' 'table' that contains the reference urls. Default column name is
//'  reference_urls'.
//' @examples
//'
//' data <- new_dataset("just for fun", 2)
//' reference_table <- readr::read_csv(rdataset_example("references.csv"),
//'                              col_names = TRUE, show_col_types = FALSE)
//' add_references(data, reference_table)
//'
//' @return double containing the number of references added
//[[Rcpp::export]]
double add_references(Rcpp::Environment data,
                     Rcpp::DataFrame table,
                     string reference_name = "reference_names",
                     string reference_version = "reference_versions",
                     string reference_usage = "reference_usages",
                     string reference_note = "reference_notes",
                     string reference_url = "reference_urls") {

    vector<string> reference_names, reference_versions, reference_usages;
    vector<string> reference_notes, reference_urls;

    reference_names = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                        reference_name));
    reference_versions = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                  "reference_versions",
                                                                  reference_version));
    reference_usages = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                           "reference_usages",
                                                                           reference_usage));
    reference_notes = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                           "reference_notes",
                                                                           reference_note));
    reference_urls = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                       "reference_urls",
                                                                       reference_url));


    bool hasVersions = false, hasUsages = false, hasNotes = false, hasUrls = false;
    if (reference_versions.size() == reference_names.size()) {
        hasVersions = true;
    }
    if (reference_usages.size() == reference_names.size()) {
        hasUsages = true;
    }
    if (reference_notes.size() == reference_names.size()) {
        hasNotes = true;
    }
    if (reference_urls.size() == reference_names.size()) {
        hasUrls = true;
    }

    vector<Reference> refs;
    for (int i = 0; i < reference_names.size(); i++) {

        Reference ref(reference_names[i]);
        if (hasVersions) {
            ref.version = reference_versions[i];
        }
        if (hasUsages) {
            ref.usage = reference_usages[i];
        }
        if (hasNotes) {
            ref.note = reference_notes[i];
        }
        if (hasUrls) {
            ref.url = reference_urls[i];
        }
        refs.push_back(ref);
    }

    Rcpp::XPtr<Dataset> d = data["data"];

    double numAdded = d.get()->addReferences(refs);

    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Function message = rdataset_env["added_message"];
    message(numAdded, "resource references");

    return numAdded;
}
/******************************************************************************/
//' @title add_sequences
//' @description
//' Add sequence data to a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing names, sequences(optional) and
//' comments(optional).
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//' @param sequence_name, a string containing the name of the column in 'table'
//' that contains the sequence names. Default column name is 'sequence_names'.
//' @param sequence, a string containing the name of the column in 'table' that
//' contains the sequence nucleotide strings. Default column name is
//' 'sequences'.
//' @param comment, a string containing the name of the column in
//' 'table' that contains the sequence comments. Default column name is
//' 'comments'.
//'
//' @examples
//'
//'  data <- new_dataset("miseq_sop", 2)
//'  fasta_data <- read_fasta(rdataset_example("final.fasta"))
//'  add_sequences(data, fasta_data)
//'
//' # With the additional parameters to add information about the reference
//'
//'  data <- new_dataset("miseq_sop", 2)
//'  fasta_data <- read_fasta(rdataset_example("final.fasta"))
//'
//'  add_sequences(data, fasta_data,
//'                new_reference("silva.bacteria.fasta",
//'                "1.38.1",
//'                "alignment by mothur2 v1.0 using default options",
//'                "https://mothur.org/wiki/silva_reference_files/"))
//'
//' # You can also add references using the 'add_references' function.
//'
//' @return double containing the number of sequences added
//[[Rcpp::export]]
double add_sequences(Rcpp::Environment data,
                     Rcpp::DataFrame table,
                     Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                     string sequence_name = "sequence_names",
                     string sequence = "sequences",
                     string comment = "comments") {

    vector<string> sequence_names, sequences, comments;
    sequence_names = Rcpp::as<vector<string>>(
        fill_required_parameters(table, sequence_name));
    sequences = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                  "sequences",
                                                                  sequence));
    comments = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                  "comments",
                                                                  comment));

    Rcpp::XPtr<Dataset> d = data["data"];

    Reference ref;
    if (reference.isNotNull()) {

        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        ref.name = Rcpp::as<string>(ref_list["reference_name"]);
        ref.version = Rcpp::as<string>(ref_list["reference_version"]);
        ref.usage = Rcpp::as<string>(ref_list["reference_usage"]);
        ref.note = Rcpp::as<string>(ref_list["reference_note"]);
        ref.url = Rcpp::as<string>(ref_list["reference_url"]);
    }

    double numAdded = d.get()->addSequences(sequence_names,
                            sequences, comments, ref);

    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Function message = rdataset_env["added_message"];
    message(numAdded);

    return numAdded;
}
/******************************************************************************/
//' @title assign_bins
//' @description
//' Add bin assignments to a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing bin_data assignments
//' @param type a string indicating the type of bin assignments. Default "otu".
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//'
//' @param bin_name, a string containing the name of the column in 'table' that
//' contains the bin names. Default column name is 'bin_names'.
//' @param abundance, a string containing the name of the column in 'table'
//' that contains the bin abundances. Default column name is 'abundances'. Note:
//'  You must provide either abundances or sequence_names in the table.
//' @param sample, a string containing the name of the column in 'table' that
//' contains the sample names for datasets where the abundances are broken down
//' by sample. Default column name is 'samples'.
//' @param sequence_name, a string containing the name of the column in 'table'
//' that contains the sequence names. Default column name is 'sequence_names'.
//' Note: You must provide either abundances or sequence_names in the table.
//'
//' @examples
//'
//'   # To assign sequences to bins:
//'
//'   data <- new_dataset("miseq_sop", 2)
//'   otu_data <- read_mothur_list(rdataset_example("final.opti_mcc.list"))
//'
//'   assign_bins(data, otu_data)
//'
//'   # To add abundance only bin assignments:
//'
//'   data <- new_dataset("miseq_sop", 2)
//'   otu_data <- read_mothur_rabund(rdataset_example("final.opti_mcc.rabund"))
//'
//'   assign_bins(data, otu_data)
//'
//'   # To add abundance bin assignments parsed by sample:
//'
//'   data <- new_dataset("miseq_sop", 2)
//'   otu_data <- readr::read_tsv(rdataset_example(
//'                                 "mothur2_bin_assignments_shared.tsv"))
//'
//'   assign_bins(data, otu_data)
//'
//' @return double containing the number of bins assigned
//[[Rcpp::export]]
double assign_bins(Rcpp::Environment data,
                   Rcpp::DataFrame table,
                   string type = "otu",
                   Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                   string bin_name = "bin_names",
                   string abundance = "abundances",
                   string sample = "samples",
                   string sequence_name = "sequence_names") {

    vector<string> bin_names, samples, sequence_names;
    vector<float> abundances;

    // fill vectors with columns from table
    bin_names = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                  bin_name));
    samples = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                "samples",
                                                                sample));
    sequence_names = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                 "sequence_names",
                                                                 sequence_name));
    abundances = Rcpp::as<vector<float>>(fill_optional_parameters(table,
                                                                  "abundances",
                                                                  abundance,
                                                                  "float"));

    if ((abundances.size() == 0) && (sequence_names.size() == 0)) {
        string message = "[ERROR]: You must provide either abundances or ";
        message += "sequence_names to assign bins.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    double numBinsAssigned = 0;
    Rcpp::XPtr<Dataset> d = data["data"];

    set<int> lengths;
    lengths.insert(bin_names.size());

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

    // if you have list assignments, don't allow setting bin abundances
    if (d.get()->hasListAssignments(type) && ((abundances.size() != 0) ||
        (samples.size() != 0))) {
        string message = "[ERROR]: You cannot assign abundance and sample data";
        message += " to bins that have sequence assignments. This could cause ";
        message += "inconsistencies.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }else{
        numBinsAssigned = d.get()->assignBins(bin_names, abundances, samples,
                                   sequence_names, type);
    }

    string tag = " " + type +" bins.";
    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Function message = rdataset_env["assigned_message"];
    message(numBinsAssigned, tag);

    Reference ref;
    if (reference.isNotNull()) {

        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        ref.name = Rcpp::as<string>(ref_list["reference_name"]);
        ref.version = Rcpp::as<string>(ref_list["reference_version"]);
        ref.usage = Rcpp::as<string>(ref_list["reference_usage"]);
        ref.note = Rcpp::as<string>(ref_list["reference_note"]);
        ref.url = Rcpp::as<string>(ref_list["reference_url"]);
        vector<Reference> refs;
        refs.push_back(ref);

        d.get()->addReferences(refs);
    }

    return numBinsAssigned;
}
/******************************************************************************/
//' @title assign_bin_representative_sequences
//' @description
//' Assign representative sequences to bins.
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing bin representative assignments
//' @param type a string indicating the type of bin assignments. Default "otu".
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//'
//' @param bin_name, a string containing the name of the column in 'table' that
//' contains the bin names. Default column name is 'bin_names'.
//' @param sequence_name a string containing the name of the column in 'table' that
//' contains the bin names. Default column name is 'sequence_names'.
//'
//' @examples
//'
//'   miseq <- miseq_sop_example()
//'
//'   num_bins <- get_num_bins(miseq, "otu")
//'
//'   # For examples sake, select first 531 sequences to be the representatives
//'   table <- data.frame(bin_names = get_bin_names(miseq, "otu"),
//'                       sequence_names = get_sequence_names(miseq)[1:num_bins]
//'                       )
//'
//'   assign_bin_representative_sequences(miseq, table, "otu")
//'
//' @return double containing the number of representative sequences assigned
//[[Rcpp::export]]
double assign_bin_representative_sequences(Rcpp::Environment data,
                                           Rcpp::DataFrame table,
                                           string type = "otu",
                                           Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                                           string bin_name = "bin_names",
                                           string sequence_name = "sequence_names") {

    vector<string> bin_names,  sequence_names;
    bin_names = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                  bin_name));
    sequence_names = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                       sequence_name));

    Rcpp::XPtr<Dataset> d = data["data"];

    double numAssigned = d.get()->assignBinRepresentativeSequences(bin_names,
             sequence_names, type);

    string tag = " " + type +" bin representative sequences.";
    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Function message = rdataset_env["assigned_message"];
    message(numAssigned, tag);

    Reference ref;
    if (reference.isNotNull()) {

        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        ref.name = Rcpp::as<string>(ref_list["reference_name"]);
        ref.version = Rcpp::as<string>(ref_list["reference_version"]);
        ref.usage = Rcpp::as<string>(ref_list["reference_usage"]);
        ref.note = Rcpp::as<string>(ref_list["reference_note"]);
        ref.url = Rcpp::as<string>(ref_list["reference_url"]);
        vector<Reference> refs;
        refs.push_back(ref);

        d.get()->addReferences(refs);
    }

    return numAssigned;
}
/******************************************************************************/
//' @title assign_bin_taxonomy
//' @description
//' Assign bin classifications to a \link{dataset} object
//'
//' Note, if you assign sequence taxonomies and assign bins, 'Dataset' will find
//'  the concensus taxonomy for each bin for you.
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing bin taxonomy assignments
//' @param type a string indicating the type of bin assignments. Default "otu".
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//'
//' @param bin_name, a string containing the name of the column in 'table' that
//' contains the bin names. Default column name is 'bin_names'.
//' @param taxonomy, a string containing the name of the column in 'table' that
//' contains the bin taxonomies. Default column name is 'taxonomies'.
//'
//' @examples
//'
//' otu_data <- read_mothur_cons_taxonomy(rdataset_example(
//'                         "final.cons.taxonomy"))
//'
//' data <- new_dataset("my_dataset", 2)
//' assign_bins(data, otu_data)
//' assign_bin_taxonomy(data, otu_data)
//'
//' @return double containing the number of bins assigned
//[[Rcpp::export]]
double assign_bin_taxonomy(Rcpp::Environment data,
                           Rcpp::DataFrame table,
                           string type = "otu",
                           Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                           string bin_name = "bin_names",
                           string taxonomy = "taxonomies") {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (d.get()->getNumBins(type) == 0) {
        string message = "[ERROR]: No bin data for type " + type + ", please ";
        message += " assign bins using the 'assign_bins' function then try ";
        message += "again.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    vector<string> bin_names, taxonomies;
    bin_names = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                  bin_name));
    taxonomies = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                   taxonomy));

    double numAssigned = d.get()->assignBinTaxonomy(bin_names,
                                                    taxonomies, type);

    string tag = " " + type +" bin taxonomies.";
    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Function message = rdataset_env["assigned_message"];
    message(numAssigned, tag);

    Reference ref;
    if (reference.isNotNull()) {

        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        ref.name = Rcpp::as<string>(ref_list["reference_name"]);
        ref.version = Rcpp::as<string>(ref_list["reference_version"]);
        ref.usage = Rcpp::as<string>(ref_list["reference_usage"]);
        ref.note = Rcpp::as<string>(ref_list["reference_note"]);
        ref.url = Rcpp::as<string>(ref_list["reference_url"]);
        vector<Reference> refs;
        refs.push_back(ref);

        d.get()->addReferences(refs);
    }

    return numAssigned;
}
/******************************************************************************/
//' @title assign_sequence_abundance
//' @description
//' Set sequence abundance and optionally assign sample and treatment data to a
//'  \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing sequence abundance assignments
//'
//' @param sequence_name, a string containing the name of the column in 'table'
//'  that contains the sequence names. Default column name is 'sequence_names'.
//' @param abundance, a string containing the name of the column in 'table'
//'  that contains the sequence abundances. Default column name is 'abundances'.
//' @param sample, a string containing the name of the column in 'table'
//' that contains the sequence samples. Default column name is 'samples'.
//' (Optional)
//' @param treatment, a string containing the name of the column in 'table'
//' that contains the sequence treatments. Default column name is 'treatments'.
//'
//' @examples
//'
//' data <- new_dataset("my_dataset", 2)
//' sequence_abundance <- readr::read_tsv(rdataset_example(
//'                                       "mothur2_count_table.tsv"),
//'                                       show_col_types = FALSE)
//'
//' assign_sequence_abundance(data, sequence_abundance, "names")
//'
//' @return double containing the number of sequences assigned
//[[Rcpp::export]]
double assign_sequence_abundance(Rcpp::Environment data,
                                 Rcpp::DataFrame table,
                               string sequence_name = "sequence_names",
                               string abundance = "abundances",
                               string sample = "samples",
                               string treatment = "treatments") {

    vector<string> sequence_names, samples, treatments;
    vector<float> abundances;

    sequence_names = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                       sequence_name));
    samples = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                "samples",
                                                                sample));
    treatments = Rcpp::as<vector<string>>(fill_optional_parameters(table,
                                                                "treatments",
                                                                treatment));
    abundances = Rcpp::as<vector<float>>(fill_optional_parameters(table,
                                                                   "abundances",
                                                                   abundance,
                                                                   "float"));
    if (sequence_names.size() != abundances.size()) {
        string message = "[ERROR]: The names and abundances must be the same";
        message += " length.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    Rcpp::XPtr<Dataset> d = data["data"];

    vector<string> unique_names = unique(sequence_names);
    vector<string> dataset_names = d.get()->getSequenceNames();

    // add seqs if needed
    if (dataset_names.size() == 0) {
        d.get()->addSequences(unique_names);
    }else {
        // sanity check, make sure names are present in dataset
        if (!identical(unique_names, dataset_names)) {
            string message = "[ERROR]: You must provide assignments for all";
            message += " sequences in your dataset.";
            RcppThread::Rcerr << endl << message << endl;
            throw Rcpp::exception(message.c_str());
        }
    }

    double numAssigned = d.get()->assignSequenceAbundance(sequence_names,
                                                          abundances,
                                                          samples, treatments);

    string tag = " sequence abundances.";
    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Function message = rdataset_env["assigned_message"];
    message(numAssigned, tag);

    return numAssigned;
}
/******************************************************************************/
//' @title assign_sequence_taxonomy
//' @description
//' Assign sequence classifications to a \link{dataset} object
//'
//' Note, if you assign sequence taxonomies and assign bins, 'Dataset' will find
//'  the concensus taxonomy for each bin for you.
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing sequence taxonomy assignments
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//'
//' @param sequence_name, a string containing the name of the column in 'table'
//' that contains the sequence names. Default column name is 'sequence_names'.
//' @param taxonomy, a string containing the name of the column in 'table' that
//' contains the sequence taxonomies. Default column name is 'taxonomies'.
//'
//' @examples
//'
//' sequence_classifications <- read_mothur_taxonomy(rdataset_example(
//'                         "final.taxonomy"))
//'
//' data <- new_dataset("my_dataset", 2)
//'
//' assign_sequence_taxonomy(data, sequence_classifications)
//'
//' # With the reference parameter you can add information about the reference
//' # you used to classify your sequences. You can also add references using the
//' # 'add_references' function.
//'
//' reference <- new_reference("trainset9_032012.pds.zip", "9_032012",
//'               "classification by mothur2 v1.0 using default options", "",
//' "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip")
//'
//' assign_sequence_taxonomy(data, sequence_classifications, reference)
//'
//' @return double containing the number of sequence assigned
//[[Rcpp::export]]
double assign_sequence_taxonomy(Rcpp::Environment data,
                                Rcpp::DataFrame table,
                                Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                                string sequence_name = "sequence_names",
                                string taxonomy = "taxonomies") {

    vector<string> sequence_names, taxonomies;
    sequence_names = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                       sequence_name));
    taxonomies = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                    taxonomy));
    // make sure names is same size as taxonomies
    if (sequence_names.size() != taxonomies.size()) {
        string message = "Size mismatch. names and taxonomies must be";
        message += " the same size.";
        throw Rcpp::exception(message.c_str());
    }

    Rcpp::XPtr<Dataset> d = data["data"];

    double numAssigned = d.get()->assignSequenceTaxonomy(sequence_names,
                               taxonomies);

    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Function message = rdataset_env["assigned_message"];
    message(numAssigned, " sequence taxonomies.");

    Reference ref;
    if (reference.isNotNull()) {

        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        ref.name = Rcpp::as<string>(ref_list["reference_name"]);
        ref.version = Rcpp::as<string>(ref_list["reference_version"]);
        ref.usage = Rcpp::as<string>(ref_list["reference_usage"]);
        ref.note = Rcpp::as<string>(ref_list["reference_note"]);
        ref.url = Rcpp::as<string>(ref_list["reference_url"]);
        vector<Reference> refs;
        refs.push_back(ref);

        d.get()->addReferences(refs);
    }

    return numAssigned;
}
/******************************************************************************/
//' @title assign_treatments
//' @description
//' Assign samples to treatments in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing sample treatment assignments
//'
//' @param sample, a string containing the name of the column in 'table'
//' that contains the samples. Default column name is 'samples'.
//' @param treatment, a string containing the name of the column in 'table'
//' that contains the treatments. Default column name is 'treatments'.
//'
//' @examples
//'
//' data <- new_dataset("my_dataset", 2)
//' sequence_abundance <- readr::read_tsv(rdataset_example(
//'                                       "mothur2_count_table.tsv"),
//'                                       show_col_types = FALSE)
//'
//' assign_sequence_abundance(data, sequence_abundance, "names")
//'
//' sample_assignments <- readr::read_table(file = rdataset_example("mouse.time.design"),
//'                              col_names = TRUE, show_col_types = FALSE)
//'
//' assign_treatments(data, sample_assignments)
//'
//' @return double containing the number of samples assigned to treatments
//[[Rcpp::export]]
double assign_treatments(Rcpp::Environment data,
                         Rcpp::DataFrame table,
                         string sample = "samples",
                         string treatment = "treatments") {

    vector<string> samples, treatments;
    samples = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                sample));
    treatments = Rcpp::as<vector<string>>(fill_required_parameters(table,
                                                                   treatment));

    // check to make sure samples and treatments are same length
    if (samples.size() != treatments.size()) {
        string message = "[ERROR]: The samples and treatments must be the same";
        message += " length.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    Rcpp::XPtr<Dataset> d = data["data"];

    if (d.get()->getNumSamples() == 0) {
        string message = "[ERROR]: You cannot assign treatments, your dataset";
        message += " does not include sample data.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    // make sure every sample in dataset is assigned a treatment
    if (!identical(d.get()->getSamples(), unique(samples))) {
        string message = "[ERROR]: You must provide treatment assignments for";
        message += " all samples in your dataset.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    double numAssigned = d.get()->assignTreatments(samples, treatments);

    Rcpp::Environment rdataset_env("package:rdataset");
    Rcpp::Function message = rdataset_env["assigned_message"];
    message(numAssigned, " samples to treatments.");

    return numAssigned;
}
/******************************************************************************/
//' @title clear
//' @description
//' Clear data from a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param tags a vector of strings containing the items you wish to clear.
//' Options are 'sequence_data', 'bin_data', 'metadata',
//' 'references', 'sequence_tree', 'sample_tree', 'alignment_report',
//' 'contigs_assembly_report' and 'chimera_report'. By default, everything
//'  is cleared.
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
//' Options are 'sequence_data' and 'bin_data'. By default, everything is
//'  exported.
//'
//' @examples
//'
//' dataset <- new_dataset("my_dataset", 2)
//' export_dataset(dataset, c(""))
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

    if ((!hasTags) || (vectorContains(t, "metadata"))) {
        Rcpp::Function get = data["get_metadata"];

        Rcpp::DataFrame metadata = get();
        if (metadata.size() != 0) {
            results.push_back(metadata);
            resultNames.push_back("metadata");
        }
    }

    if ((!hasTags) || (vectorContains(t, "alignment_report"))) {

        Rcpp::Function get = data["get_alignment_report"];
        Rcpp::DataFrame alignment_report = get();

        if (alignment_report.size() != 0) {
            results.push_back(alignment_report);
            resultNames.push_back("alignment_report");
        }
    }

    if ((!hasTags) || (vectorContains(t, "chimera_report"))) {

        Rcpp::Function get = data["get_chimera_report"];
        Rcpp::DataFrame chimera_report = get();

        if (chimera_report.size() != 0) {
            results.push_back(chimera_report);
            resultNames.push_back("chimera_report");
        }
    }

    if ((!hasTags) || (vectorContains(t, "contigs_assembly_report"))) {

        Rcpp::Function get = data["get_contigs_assembly_report"];
        Rcpp::DataFrame contigs_report = get();

        if (contigs_report.size() != 0) {
            results.push_back(contigs_report);
            resultNames.push_back("contigs_assembly_report");
        }
    }

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
//' @title get_bin
//' @description
//' Get the names of the sequences in a given bin in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param bin_name, string containing the bin name
//' @param type a string indicating the type of clusters. Default = "otu".
//' @examples
//'
//' data <- miseq_sop_example()
//' get_bin(data, "Otu005", "otu")
//'
//' @return String, containing names of the sequences in a given bin
//[[Rcpp::export]]
string get_bin(Rcpp::Environment data,
               string bin_name,
               string type = "otu") {

    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getBin(bin_name, type);
}
/******************************************************************************/
//' @title get_bin_abundance
//' @description
//' Get the abundance of a given bin in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param bin_name, string containing the bin name
//' @param type a string indicating the type of clusters. Default = "otu".
//' @examples
//'
//' data <- miseq_sop_example()
//' get_bin_abundance(data, "Otu005", "otu")
//'
//' @return double, containing the abundance of a given bin
//[[Rcpp::export]]
double get_bin_abundance(Rcpp::Environment data,
                       string bin_name, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getBinAbundance(bin_name, type);
}
/******************************************************************************/
//' @title get_bin_abundances
//' @description
//' Get the abundance of a given bin parsed by sample in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param bin_name, string containing the bin name
//' @param type, string indicating the type of clusters. Default = "otu".
//' @examples
//'
//' data <- miseq_sop_example()
//' get_bin_abundances(data, "Otu005", "otu")
//'
//' @return vector containing the abundance of a given bin parsed by sample
//[[Rcpp::export]]
vector<float> get_bin_abundances(Rcpp::Environment data,
                       string bin_name, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getBinAbundances(bin_name, type);
}
/******************************************************************************/
//' @title get_bin_names
//' @description
//' Get the names of the bins in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param type, string indicating the type of clusters. Default = "otu".
//' @examples
//'
//'   miseq <- miseq_sop_example()
//'   get_bin_names(miseq, "phylotype")
//'
//' @return vector containing the names of bins
//[[Rcpp::export]]
vector<string> get_bin_names(Rcpp::Environment data, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getBinIds(type);
}
/******************************************************************************/
//' @title get_bin_representative_sequences
//' @description
//' Get the representative sequences of the bins in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param type, string indicating the type of clusters. Default = "otu".
//' @examples
//'
//'   miseq <- miseq_sop_example()
//'
//'   num_bins <- get_num_bins(miseq, "otu")
//'
//'   # For examples sake, select first 531 sequences to be the representatives
//'   table <- data.frame(bin_names = get_bin_names(miseq, "otu"),
//'                       sequence_names = get_sequence_names(miseq)[1:num_bins]
//'                       )
//'
//'   assign_bin_representative_sequences(miseq, table, "otu")
//'
//'
//'   get_bin_representative_sequences(miseq, "otu")
//'
//' @return data.frame containing 2 columns representative_names and
//'  representative_sequences
//[[Rcpp::export]]
Rcpp::DataFrame get_bin_representative_sequences(Rcpp::Environment data,
                                                 string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getBinRepresentativeSequences(type);
}
/******************************************************************************/
//' @title get_bin_taxonomy_report
//' @description
//' Get the bin classifications of a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param type, string indicating the type of bin clusters. Default = "otu".
//' @examples
//'
//' data <- miseq_sop_example()
//' get_bin_taxonomy_report(data, "otu")
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_bin_taxonomy_report(Rcpp::Environment data,
                                        string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getBinTaxonomyReport(type);
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
//' @title get_dataset_name
//' @description
//' Get the name of a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @examples
//' data <- new_dataset("my_dataset", 2)
//' get_dataset_name(data)
//'
//' @return String, containing the name of the dataset
//[[Rcpp::export]]
string get_dataset_name(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->datasetName;
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
//' @title get_num_processors
//' @description
//' Get the number of processors used to summarize a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @examples
//' data <- new_dataset("my_dataset", 2)
//' get_num_processors(data)
//'
//' @return Integer, containing number of processors
//[[Rcpp::export]]
int get_num_processors(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->processors;
}
/******************************************************************************/
//' @title get_num_bins
//' @description
//' Get the number of bins of a specific type in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param type a string indicating the type of clusters. Default = "otu".
//' @examples
//'
//' data <- miseq_sop_example()
//' get_num_bins(data, "otu")
//' get_num_bins(data, "asv")
//' get_num_bins(data, "phylotype")
//'
//' @return Integer, the number of bins of a specific type in a \link{dataset}
//'  object
//[[Rcpp::export]]
int get_num_bins(Rcpp::Environment data, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getNumBins(type);
}
/******************************************************************************/
//' @title get_num_samples
//' @description
//' Get the number of samples in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @examples
//'
//' data <- new_dataset()
//' bin_table <- readr::read_tsv(rdataset_example(
//'                               "mothur2_bin_assignments_shared.tsv"),
//'                               show_col_types = FALSE)
//' assign_bins(data, bin_table)
//' get_num_samples(data)
//'
//' @return Integer, the number of samples in a \link{dataset} object
//[[Rcpp::export]]
int get_num_samples(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getNumSamples();
}
/******************************************************************************/
//' @title get_num_sequences
//' @description
//' Get the number of sequences in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param distinct Boolean. When distinct is TRUE the number of unique
//' sequence is returned. Default = FALSE.
//' @param sample, string containing the name of the sample you want number of
//'  sequences for.
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' get_num_sequences(data)
//' get_num_sequences(data, TRUE)
//' get_num_sequences(data, FALSE, "F3D0")
//' get_num_sequences(data, TRUE, "F3D0")
//'
//' @return double, the number of sequences in a \link{dataset} object
//[[Rcpp::export]]
double get_num_sequences(Rcpp::Environment data, bool distinct = false,
                            string sample = "") {
    Rcpp::XPtr<Dataset> d = data["data"];

    if (distinct) {
        return d.get()->getUniqueTotal(sample);
    }

    return d.get()->getTotal(sample);
}
/******************************************************************************/
//' @title get_num_treatments
//' @description
//' Get the number of treatments in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' get_num_treatments(data)
//'
//' @return Integer, the number of treatments in a \link{dataset} object
//[[Rcpp::export]]
int get_num_treatments(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getNumTreatments();
}
/******************************************************************************/
//' @title get_rabund
//' @description
//' Get data.frame containing bin abundance data in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   data <- miseq_sop_example()
//'   get_rabund(data)
//'
//' @return a 2 column data.frame containing bin names and bin abundances
//[[Rcpp::export]]
Rcpp::DataFrame get_rabund(Rcpp::Environment data, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getRAbund(type);
}
/******************************************************************************/
//' @title get_rabund_vector
//' @description
//' Get vector containing bin abundance data in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   data <- miseq_sop_example()
//'   get_rabund_vector(data)
//'
//' @return vector containing each bins abundance
//[[Rcpp::export]]
vector<float> get_rabund_vector(Rcpp::Environment data, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getRAbundVector(type);
}
/******************************************************************************/
//' @title get_references
//' @description
//' Get a table containing resource references in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_references(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getReferences();
}
/******************************************************************************/
//' @title get_samples
//' @description
//' Get the samples in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' get_samples(data)
//'
//' @return vector of strings containing the names of the samples in a
//' \link{dataset} object
//[[Rcpp::export]]
vector<string> get_samples(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getSamples();
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
//' @title get_sample_totals
//' @description
//' Get the number of sequences in each sample in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' get_sample_totals(data)
//'
//' @return vector containing the number of sequences in each sample in a
//' \link{dataset} object
//[[Rcpp::export]]
vector<double> get_sample_totals(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
   return d.get()->getSampleTotals();
}
/******************************************************************************/
//' @title get_scrap_report
//' @description
//' Get a scrap report containing sequences and bins eliminated from a
//' \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param type a string indicating the type of scrap report you would like.
//'  Default = 'sequence'.
//' @examples
//'
//'   data <- miseq_sop_example()
//'   remove_bins(data, c("Otu005"), c("bad_bin"))
//'
//'   sequence_scrap_report <- get_scrap_report(data, "sequence")
//'   otu_scrap_report <- get_scrap_report(data, "otu")
//'
//' @return data.frame containing sequences or otus removed from the
//' \link{dataset} object during analysis
//[[Rcpp::export]]
Rcpp::DataFrame get_scrap_report(Rcpp::Environment data,
                               string type = "sequence") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getScrapReport(type);
}
/******************************************************************************/
//' @title get_sequence_abundances
//' @description
//' Get the total abundance for each sequence in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' get_sequence_abundances(data)
//'
//' @return vector containing the total abundance for each sequence a
//' \link{dataset} object
//[[Rcpp::export]]
vector<float> get_sequence_abundances(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getSequenceAbundances();
}
/******************************************************************************/
//' @title get_sequence_abundances_by_sample
//' @description
//' Get the abundances of each sequence in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' get_sequence_abundances_by_sample(data)
//'
//' @return 2D vector ([num_seqs][num_samples]) containing the
//' abundances of each sequence in an instance of the 'Dataset' class parsed by
//' sample.
//[[Rcpp::export]]
vector<vector<float> > get_sequence_abundances_by_sample(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getSeqsAbundsBySample();
}
/******************************************************************************/
//' @title get_sequence_abundance_table
//' @description
//' Get the abundances of each sequence in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @examples
//'
//' data <- miseq_sop_example()
//' get_sequence_abundance_table(data)
//'
//' @return data.frame containing sequence abundance data in a \link{dataset}
//' object
//[[Rcpp::export]]
Rcpp::DataFrame get_sequence_abundance_table(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequenceAbundanceTable();
}
/******************************************************************************/
//' @title get_sequence_names
//' @description
//' Get the names of the sequences in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param sample a string containing the name of the sample you
//' would like sequence names for. For all samples in dataset, sample = "".
//'
//' @examples
//'
//'
//' data <- miseq_sop_example()
//'
//' # to get the names of all the sequences in the dataset
//'
//' get_sequence_names(data)
//'
//' # to get the names of the sequences in sample 'F3D0' from the dataset
//'
//' get_sequence_names(data, "F3D0")
//'
//' @return vector of string containing the names of the sequences a
//' \link{dataset} object
//[[Rcpp::export]]
vector<string> get_sequence_names(Rcpp::Environment data, string sample = "") {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getSequenceNames(sample);
}
/******************************************************************************/
//' @title get_sequence_names_by_sample
//' @description
//' Get the names of the sequences in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param samples a vector of strings containing the names of the samples you
//' would like sequence names for. By default all samples are included.
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' get_sequence_names_by_sample(data)
//'
//' @return 2D vector of strings ([num_seqs][num_samples]) containing the
//' names of each sequence a \link{dataset} object parsed by sample.
//[[Rcpp::export]]
vector<vector<string> > get_sequence_names_by_sample(Rcpp::Environment data,
                    Rcpp::CharacterVector samples = Rcpp::CharacterVector::create()) {

    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequenceNamesBySample(Rcpp::as<vector<string>>(samples));
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
//' @title get_sequences_by_sample
//' @description
//' Get the nucleotide strings for each sequence in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param samples a vector of strings containing the names of the samples you
//' would like sequence names for. By default all samples are included.
//' @examples
//'
//' data <- miseq_sop_example()
//' get_sequences_by_sample(data)
//'
//' @return 2D vector of strings ([num_seqs][num_samples]) containing the
//' nucleotide strings for each sequence a \link{dataset} object parsed by sample.
//[[Rcpp::export]]
vector<vector<string> > get_sequences_by_sample(Rcpp::Environment data,
                 Rcpp::CharacterVector samples = Rcpp::CharacterVector::create()) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequencesBySample(Rcpp::as<vector<string>>(samples));
}
/******************************************************************************/
//' @title get_sequence_report
//' @description
//' Get sequence report data: starts, ends, lengths, ambigs, longest
//' homopolymers and numns.
//'
//' @param data, a \link{dataset} object
//' @examples
//'
//'  data <- miseq_sop_example()
//'  get_sequence_report(data)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_sequence_report(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequenceReport();
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
//'  remove_samples(data, c("F3D0"))
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
//' @title get_sequence_taxonomy_report
//' @description
//' Get the sequence classifications of a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @examples
//'
//' data <- miseq_sop_example()
//' get_sequence_taxonomy_report(data)
//'
//' @return data.frame
//[[Rcpp::export]]
Rcpp::DataFrame get_sequence_taxonomy_report(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequenceTaxonomyReport();
}
/******************************************************************************/
//' @title get_bin_assignments
//' @description
//' Get data.frame containing bin abundance data of a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'   data <- miseq_sop_example()
//'   get_bin_assignments(data, "otu")
//'
//' @return data.frame containing 2, 3 or 4 columns: bin_names, abundances,
//' samples (if assigned), and treatments (if assigned)
//[[Rcpp::export]]
Rcpp::DataFrame get_bin_assignments(Rcpp::Environment data, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getShared(type);
}
/******************************************************************************/
//' @title get_shared_vector
//' @description
//' Get 2D vector containing bin abundance data by sample in a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//'  data <- miseq_sop_example()
//'  get_shared_vector(data)
//'
//' @return 2D vector ([num_bins][num_samples]) containing the abundances of
//' each bin parsed by sample.
//[[Rcpp::export]]
vector<vector<float> > get_shared_vector(Rcpp::Environment data,
                                         string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getSharedVector(type);
}
/******************************************************************************/
//' @title get_treatments
//' @description
//' Get the treatments in a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//' @examples
//'
//'  data <- miseq_sop_example()
//'  get_treatments(data)
//'
//' @return vector of strings containing the names of the treatments in a
//' \link{dataset} object
//[[Rcpp::export]]
vector<string> get_treatments(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getTreatments();
}
/******************************************************************************/
//' @title get_treatment_totals
//' @description
//' Get the number of sequences in each treatment in a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//' @examples
//'
//' data <- miseq_sop_example()
//' get_treatment_totals(data)
//'
//' @return vector containing the number of sequences in each treatment in a
//' \link{dataset} object
//[[Rcpp::export]]
vector<double> get_treatment_totals(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
     return d.get()->getTreatmentTotals();
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
//' @title merge_bins
//' @description
//' Designed with package integration in mind, the merge bins function allows
//' you to merge bins in a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//'
//' @param bin_names, a vector of strings containing the names of the bins you
//' would like merge. The resulting merged bin will be stored in the first
//' bin_id in the vector.
//' @param reason, a string indicating why you are merging bins. Default =
//' "merged".
//' @param type, a string indicating the type of bin clusters. Default = "otu"
//'
//' @examples
//'
//'  data <- miseq_sop_example()
//'
//'  # to merge otu5 and otu6
//'
//'  bins_to_merge <- c("Otu005", "Otu006")
//'
//'  merge_bins(data, bins_to_merge)
//'
//'  # If you look at the scrap report, you will see Otu006 with the trash code
//'  # set to "merged".
//'
//'  get_scrap_report(data, "bin")
//'
//[[Rcpp::export]]
void merge_bins(Rcpp::Environment data, vector<string> bin_names,
                string reason = "merged", string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->mergeBins(bin_names, reason, type);
}
/******************************************************************************/
//' @title merge_sequences
//' @description
//' Designed with package integration in mind, the merge sequences function
//' allows you to merge sequences in a \link{dataset} object.
//'
//' @param data, a \link{dataset} object.
//'
//' @param sequence_names, a vector of strings containing the names of the
//' sequences you would like merge. The resulting merged sequence will be stored
//' in the first sequence name in the vector.
//' @param reason a string indicating why you are merging sequences.
//' Default = "merged"
//'
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' get_num_sequences(data)
//'
//' # For the sake of example let's merge the first 3 sequences from
//' # miseq_sop_example:
//'
//' seqs_to_merge <- c("M00967_43_000000000-A3JHG_1_2101_16474_12783",
//'                    "M00967_43_000000000-A3JHG_1_1113_12711_3318",
//'                    "M00967_43_000000000-A3JHG_1_2108_14707_9807")
//'
//' merge_sequences(data, seqs_to_merge)
//'
//' # If you look at the scrap report, you will see the second two sequence
//' # names, listed with the trash code set to "merged".
//'
//' get_scrap_report(data)
//'
//' # You can see from the get_num_sequences function that the merged sequence's
//' # abundances are added to the first sequence.
//'
//' get_num_sequences(data)
//'
//[[Rcpp::export]]
void merge_sequences(Rcpp::Environment data, vector<string> sequence_names,
                 string reason = "merged") {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->mergeSequences(sequence_names, reason);
}
/******************************************************************************/
//' @title remove_bins
//' @description
//' Designed with package integration in mind, the remove bins function allows
//' you to remove bins from a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//'
//' @param bin_names, a vector of strings containing the names of the bins you
//' would like removed.
//' @param trash_tags, a vector of strings containing the reasons you are
//' removing each bin
//' @param type a string indicating the type of clusters.
//' @examples
//'
//'   data <- new_dataset("my_dataset", 2)
//'
//'   bin_names <- c("bin1", "bin2", "bin3")
//'   abundances <- c(110, 525, 80)
//'
//'   assign_bins(data, data.frame(bin_names = bin_names,
//'                                abundances = abundances))
//'
//'   get_num_bins(data)
//'
//'   bins_to_remove <- c("bin1")
//'   trash_tag <- c("bad_bin")
//'
//'   remove_bins(data, bins_to_remove, trash_tag)
//'
//'   get_num_bins(data)
//'
//[[Rcpp::export]]
void remove_bins(Rcpp::Environment data, vector<string> bin_names,
                 vector<string> trash_tags, string type = "otu") {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->removeBins(bin_names, trash_tags, type);
}
/******************************************************************************/
//' @title remove_lineages
//' @description
//' Designed with package integration in mind, the remove lineages function
//' allows you to remove contaminents from a \link{dataset}
//'
//' @param data, a \link{dataset} object.
//'
//' @param contaminants, vector of strings containing the taxonomies you would
//' like to remove
//' @param trash_tag, a string containing reason you are removing the lineages.
//' Default = "contaminant".
//'
//' @examples
//' data <- read_mothur(fasta = rdataset_example("final.fasta"),
//'                        count = rdataset_example("final.count_table"),
//'                        taxonomy = rdataset_example("final.taxonomy"),
//'                        design = rdataset_example("mouse.time.design"),
//'                        otu_list = rdataset_example("final.opti_mcc.list"),
//'                        dataset_name = "miseq_sop")
//'
//' contaminants <- c("Chloroplast", "Mitochondria", "unknown", "Archaea",
//'  "Eukaryota")
//'
//' remove_lineages(data, contaminants)
//'
//[[Rcpp::export]]
void remove_lineages(Rcpp::Environment data, vector<string> contaminants,
                     string trash_tag = "contaminant") {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->removeLineages(contaminants, trash_tag);
}
/******************************************************************************/
//' @title remove_samples
//' @description
//' Remove samples from a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//'
//' @param samples, vector of strings containing the names of the samples to
//' remove.
//'
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' get_num_samples(data)
//'
//' # To remove samples 'F3D0' and 'F3D1'
//'
//' remove_samples(data, c("F3D0", "F3D1"))
//'
//' get_num_samples(data)
//'
//[[Rcpp::export]]
void remove_samples(Rcpp::Environment data, vector<string> samples) {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->removeSamples(samples);
}
/******************************************************************************/
//' @title remove_sequences
//' @description
//' Designed with package integration in mind, the remove sequences function
//' allows you to remove sequences from a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//'
//' @param sequence_names, vector of strings containing the names of the
//' sequences to remove
//' @param trash_tags vector of strings containing the reasons for the sequences
//' removals
//'
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' get_num_sequences(data)
//'
//' # For the sake of example let's remove the first 3 sequences from
//' # miseq_sop_example:
//'
//' seqs_to_remove <- c("M00967_43_000000000-A3JHG_1_2101_16474_12783",
//'                    "M00967_43_000000000-A3JHG_1_1113_12711_3318",
//'                    "M00967_43_000000000-A3JHG_1_2108_14707_9807")
//' trash_codes <- c("example", "removing", "sequences")
//'
//' remove_sequences(data, seqs_to_remove, trash_codes)
//'
//' # If you look at the scrap report, you the sequences names, listed with the
//' # trash codes set to "example", "removing", "sequences".
//'
//' get_scrap_report(data)
//'
//' # You can see from the get_num_sequences function that the removed
//' # sequence's abundances are removed from the dataset.
//'
//' get_num_sequences(data)
//'
//[[Rcpp::export]]
void remove_sequences(Rcpp::Environment data,
                      vector<string> sequence_names,
                      vector<string> trash_tags) {

    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->removeSequences(sequence_names, trash_tags);
}
/******************************************************************************/
//' @title set_abundance
//' @description
//' Designed with package integration in mind, the set abundance function
//' allows you to change the abundances of sequences in a \link{dataset} object
//' without samples.
//'
//' @param data, a \link{dataset} object
//'
//' @param sequence_names, a vector of strings containing sequence names
//' @param sequence_abundances, vector containing the abundances of each
//' sequence.
//' @param reason, a string containing the trash tag to be applied to any
//' sequences set to 0 abundance. Default = "update".
//'
//' @examples
//'
//' names <- c("seq1", "seq2", "seq3",  "seq4")
//' abunds <- c(1250, 65, 50, 4)
//'
//' data <- new_dataset("my_dataset", 2)
//'
//' assign_sequence_abundance(data, data.frame(sequence_names = names,
//'                                            abundances = abunds))
//' get_sequence_abundances(data)
//'
//' seqs_to_update <- c("seq1", "seq3")
//' new_abunds <- c(1000, 100)
//'
//' set_abundance(data, seqs_to_update, new_abunds)
//'
//' get_sequence_abundances(data)
//'
//[[Rcpp::export]]
void set_abundance(Rcpp::Environment data,
                   vector<string> sequence_names,
                   vector<float> sequence_abundances,
                   string reason = "update") {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (get_num_samples(data) == 0) {
        d.get()->setAbundance(sequence_names, sequence_abundances, reason);
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
//' Designed with package integration in mind, the set abundances function
//' allows you to change the abundances of sequences in a \link{dataset} object
//' with samples.
//'
//' @param data, a \link{dataset} object
//'
//' @param sequence_names, a vector of strings containing sequence names
//' @param abundances, 2D vector ([num_seqs][num_samples]) containing
//' the abundances of each sequence parsed by sample.
//' @param reason, a string containing the trash tag to be applied to any
//' sequences set to 0 abundance. Default = "update".
//' @examples
//'
//' data <- new_dataset("my_dataset", 2)
//'
//' sequence_names <- c("seq1", "seq1", "seq1", "seq2", "seq2", "seq2", "seq3",
//'                     "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4", "sample2", "sample3",
//'              "sample4", "sample2", "sample3", "sample4")
//' abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)
//'
//' assign_sequence_abundance(data, data.frame(sequence_names = sequence_names,
//'                                            abundances = abundances,
//'                                            samples = samples))
//'
//' seqs_to_update <- c("seq4")
//' new_abunds <- list(c(20, 10, 4))
//'
//' set_abundances(data, seqs_to_update, new_abunds)
//'
//[[Rcpp::export]]
void set_abundances(Rcpp::Environment data,
                   vector<string> sequence_names,
                   vector<vector<float>> abundances,
                   string reason = "update") {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (get_num_samples(data) != 0) {
        d.get()->setAbundances(sequence_names, abundances, reason);
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
//' Designed with package integration in mind, the set bin abundance function
//' allows you to change the abundances of bins in a \link{dataset} object
//' without sample data.
//'
//' @param data, a \link{dataset} object
//'
//' @param bin_names, a vector strings containing of bin names to set the
//' abundances for.
//' @param abundances, vector containing the abundances of each bin.
//' @param type, a string indicating the type of clusters. Default = "otu".
//' @param reason, a string containing the trash tag to be applied to any bins
//'  set to 0 abundance. Default = "update".
//'
//' @examples
//'   # For example sake, let's create a dataset with 3 bins:
//'
//'   data <- new_dataset("my_dataset", 2)
//'
//'   bin_ids <- c("bin1", "bin2", "bin3")
//'   abundances <- c(110, 525, 80)
//'
//'   assign_bins(data, data.frame(bin_names = bin_ids,
//'                                abundances = abundances))
//'
//'   get_bin_abundance(data, "bin1")
//'   get_bin_abundance(data, "bin2")
//'
//'   # Now we can use set_bin_abundance to change the abundances of bin1 and
//'   # bin2
//'
//'   bins <- c("bin1", "bin2")
//'   new_abunds <- c(300, 250)
//'
//'   set_bin_abundance(data, bins, new_abunds)
//'
//'   get_bin_abundance(data, "bin1")
//'   get_bin_abundance(data, "bin2")
//'
//[[Rcpp::export]]
void set_bin_abundance(Rcpp::Environment data,
                       vector<string> bin_names,
                       vector<float> abundances,
                       string type = "otu",
                       string reason = "update") {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (d.get()->hasListAssignments(type)) {
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
        d.get()->setBinAbundance(bin_names, abundances, reason, type);
    }
}
/******************************************************************************/
//' @title set_bin_abundances
//' @description
//' Designed with package integration in mind, the set bin abundances function
//' allows you to change the abundances of bins in a \link{dataset} object
//' with sample data.
//'
//' @param data, a \link{dataset} object
//'
//' @param bin_names, a vector strings containing of bin names to set the
//' abundances for.
//' @param abundances, 2D vector ([num_seqs][num_samples]) containing the
//' abundances of each bin parsed by sample.
//' @param type a string indicating the type of clusters. Default = "otu".
//' @param reason, a string containing the trash tag to be applied to any bins
//'  set to 0 abundance. Default = "update".
//'
//' @examples
//'
//'   # For example sake, let's create a dataset with 3 bins:
//'
//'   data <- new_dataset("my_dataset", 2)
//'
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   samples <- c("sample1", "sample2", "sample5", "sample1", "sample3",
//'                "sample1")
//'   sample_abundances <- c(10, 100, 1, 500, 25, 80)
//'
//'   assign_bins(data, data.frame(bin_names = bin_ids,
//'                                abundances = sample_abundances,
//'                                samples = samples))
//'
//'   # You can see bin1's abundances parsed by sample using get_bin_abundances:
//'   get_bin_abundances(data, "bin1")
//'
//'   # You can change bin1's abundances as follows:
//'
//'   new_bin1_abunds <- list(c(10,50,0,0))
//'   bins <- c("bin1")
//'
//'   set_bin_abundances(data, bins, new_bin1_abunds)
//'
//'   get_bin_abundances(data, "bin1")
//'
//[[Rcpp::export]]
void set_bin_abundances(Rcpp::Environment data,
                        vector<string> bin_names,
                        vector<vector<float>> abundances,
                        string type = "otu", string reason = "update") {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (d.get()->hasListAssignments(type)) {
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
        d.get()->setBinAbundances(bin_names, abundances, reason, type);
    }
}
/******************************************************************************/
//' @title set_sequences
//' @description
//' Designed with package integration in mind, the set sequences function
//' allows you to change the nucleotide strings of sequences in a \link{dataset}
//' object. For example, set_sequences may be used after alignment to overwrite
//' the unaligned sequences with aligned sequences.
//'
//' @param data, a \link{dataset} object
//' @param sequence_names, a vector of strings containing sequence names
//' @param sequences, a vector of strings containing sequence nucleotide strings
//' @param comments, a vector of strings containing sequence comments.
//' (Optional)
//'
//' @examples
//'
//' data <- new_dataset("my_dataset", 2)
//'
//' add_sequences(data, data.frame(sequence_names = c("seq1", "seq2",
//'                                                   "seq3", "seq4")))
//'
//' set_sequences(data, c("seq1", "seq2","seq3", "seq4"),
//'                     c("ATTGC", "ACTGC", "AGTGC", "TTTGC"))
//'
//[[Rcpp::export]]
void set_sequences(Rcpp::Environment data,
                   vector<string> sequence_names,
                   vector<string> sequences,
                   Rcpp::CharacterVector comments = Rcpp::CharacterVector::create()) {

    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->setSequences(sequence_names, sequences,
                          Rcpp::as<vector<string>>(comments));
}
/******************************************************************************/
//' @title set_dataset_name
//' @description
//' Set the name of a \link{dataset} object.
//'
//' @param data, a \link{dataset} object
//' @param dataset_name, a string containing the desired name
//'
//' @examples
//'
//' data <- new_dataset("my_dataset", 2)
//' set_dataset_name(data, "new_dataset_name")
//'
//[[Rcpp::export]]
void set_dataset_name(Rcpp::Environment data, string dataset_name) {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->datasetName = dataset_name;
}
/******************************************************************************/
//' @title set_num_processors
//' @description
//' Set the number of processors used to summarize a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param processors, a integer containing the desired number of processors
//' @examples
//'
//' data <- new_dataset("my_dataset", 2)
//' set_num_processors(data, 1)
//'
//[[Rcpp::export]]
void set_num_processors(Rcpp::Environment data, int processors) {
    Rcpp::XPtr<Dataset> d = data["data"];
    d.get()->processors = processors;
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
//' @title new_pointer
//' @name new_pointer
//' @description
//' For internal use only, create an instance of the C++ 'Dataset' class.
//' @param dataset_name, string containing dataset name
//' @param processors, number of processors to use
//' @return pointer to an instance of the C++ 'Dataset' class.
//' @keywords internal
//[[Rcpp::export]]
Rcpp::XPtr<Dataset> new_pointer(string dataset_name = "", int processors = 1) {
     Dataset* d = new Dataset(dataset_name, processors);
     return Rcpp::XPtr<Dataset>(d);
}
/******************************************************************************/
//' @title copy_pointer
//' @name copy_pointer
//' @description
//' For internal use only, copy an instance of the C++ 'Dataset' class.
//' @param data, a \link{dataset} object
//' @return pointer to an instance of the C++ 'Dataset' class.
//' @keywords internal
//[[Rcpp::export]]
Rcpp::XPtr<Dataset> copy_pointer(Rcpp::Environment data) {
    Rcpp::XPtr<Dataset> d = data["data"];
    Dataset* copy = new Dataset(*d.get());
    return Rcpp::XPtr<Dataset>(copy);
}
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
//' @title deserialize_dobject
//' @name deserialize_dobject
//' @description
//' For internal use only, deserialize_dobject an instance of the C++ 'Dataset'
//'  class.
//' @param data, a \link{dataset} object
//[[Rcpp::export]]
void deserialize_dobject(Rcpp::Environment data) {
     data["data"] = new_pointer();
     Rcpp::XPtr<Dataset> d = data["data"];
     const Rcpp::RawVector raw = data["raw"];
     d.get()->loadFromSerialized(raw);
}
/******************************************************************************/
//' @title serialize_dobject
//' @name serialize_dobject
//' @description
//' For internal use only, serialize_dobject an instance of the C++ 'Dataset'
//' class.
//' @param data, a \link{dataset} object
//[[Rcpp::export]]
void serialize_dobject(Rcpp::Environment data) {
     Rcpp::XPtr<Dataset> d = data["data"];
     data["raw"] = d.get()->serializeDataset();
 }
/******************************************************************************/
