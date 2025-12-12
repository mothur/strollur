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
//' @title add_report
//' @description
//' Add a report to a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//'
//' @param table, a data.frame containing your report.
//'
//' @param type, a string containing the type of report. Options include:
//' "metadata" and custom report tags. Default = "metadata".
//'
//' @param sequence_name, a string containing the name of the column in 'table'
//' that contains the sequence names. This is used for custom reports, metadata
//' does not require a sequence_name column. Default column name is 'sequence_names'.
//' @examples
//'
//' # To add a custom report including your contigs assembly data
//'
//' data <- new_dataset("just for fun", 2)
//' contigs_report <- readr::read_tsv(rdataset_example("final.contigs_report"),
//'    col_names = TRUE, show_col_types = FALSE)
//'
//' add_report(data, contigs_report, "contigs_report", "Name")
//'
//' # To add metadata related to your study
//'
//' metadata <- readr::read_tsv(rdataset_example("mouse.dpw.metadata"),
//'                             col_names = TRUE, show_col_types = FALSE)
//'
//' add_report(data, metadata, "metadata")
//'
//[[Rcpp::export]]
void add_report(Rcpp::Environment data,
                       Rcpp::DataFrame table,
                       string type = "metadata",
                       string sequence_name = "sequence_names") {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (type == "metadata") {
        d.get()->addMetadata(table);
        Rcpp::Environment rdataset_env("package:rdataset");
        Rcpp::Function message = rdataset_env["added_message"];
        message(R_NilValue, type);
    }
    else{

        vector<string> dfNames = Rcpp::as<vector<string>>(table.attr("names"));

        // do we have a sequence_names column
        if (!vectorContains(dfNames, sequence_name)) {
            string message = "The report must include a column containing sequence";
            message += " names. " + sequence_name + " is not a named column in your report.";
            throw Rcpp::exception(message.c_str());
        }else {

            // sequence names in table
            vector<string> sequenceNames = Rcpp::as<vector<string>>(
                xint_fill_required_parameters(table, sequence_name));

            // if we don't have any sequences yet, add them
            if (d.get()->getTotal() == 0) {
                vector<string> sequences, comments;
                Reference ref;

                d.get()->addSequences(sequenceNames,
                    sequences, comments, ref);
            }else {
                // we have sequences already, make sure there is a report row for
                // each sequence in dataset
                vector<string> datasetSeqNames = d.get()->getSequenceNames();

                // find seqs in dataset and not in report
                vector<string> missingSeqs = setDiff(datasetSeqNames, sequenceNames);

                if (missingSeqs.size() == 0) {
                    // preserve order of dataset
                    Rcpp::Environment rdataset_env("package:rdataset");
                    Rcpp::Function sort = rdataset_env["sort_dataframe"];

                    table = sort(table, datasetSeqNames, sequence_name);
                }else {
                    string message = "Your report does not contain an entry for ";
                    message += "every sequence in your dataset, ignoring report. ",
                    RcppThread::Rcout << endl << message << endl;
                    return;
                }
            }

            // save name column
            table.attr("sequence_name") = sequence_name;

            // add report to dataset
            d.get()->addReport(table, type);

            Rcpp::Environment rdataset_env("package:rdataset");
            Rcpp::Function message = rdataset_env["added_message"];
            message(R_NilValue, "a " + type);
        }
    }
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

    reference_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                        reference_name));
    reference_versions = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                  "reference_versions",
                                                                  reference_version));
    reference_usages = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                           "reference_usages",
                                                                           reference_usage));
    reference_notes = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                           "reference_notes",
                                                                           reference_note));
    reference_urls = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
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
        xint_fill_required_parameters(table, sequence_name));
    sequences = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                  "sequences",
                                                                  sequence));
    comments = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
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
//' @param bin_type a string indicating the type of bin assignments. Default "otu".
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
//'   assign_bins(data = data, table = otu_data, bin_type = "otu")
//'
//'   # To add abundance only bin assignments:
//'
//'   data <- new_dataset("miseq_sop", 2)
//'   otu_data <- read_mothur_rabund(rdataset_example("final.opti_mcc.rabund"))
//'
//'   assign_bins(data = data, table = otu_data, bin_type = "otu")
//'
//'   # To add abundance bin assignments parsed by sample:
//'
//'   data <- new_dataset("miseq_sop", 2)
//'   otu_data <- readr::read_tsv(rdataset_example(
//'                                 "mothur2_bin_assignments_shared.tsv"))
//'
//'   assign_bins(data = data, table = otu_data, bin_type = "otu")
//'
//' @return double containing the number of bins assigned
//[[Rcpp::export]]
double assign_bins(Rcpp::Environment data,
                   Rcpp::DataFrame table,
                   string bin_type = "otu",
                   Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                   string bin_name = "bin_names",
                   string abundance = "abundances",
                   string sample = "samples",
                   string sequence_name = "sequence_names") {

    vector<string> bin_names, samples, sequence_names;
    vector<float> abundances;

    // fill vectors with columns from table
    bin_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                  bin_name));
    samples = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                "samples",
                                                                sample));
    sequence_names = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                 "sequence_names",
                                                                 sequence_name));
    abundances = Rcpp::as<vector<float>>(xint_fill_optional_parameters(table,
                                                                  "abundances",
                                                                  abundance,
                                                                  "float"));

    if ((abundances.size() == 0) && (sequence_names.size() == 0)) {
        string message = "You must provide either abundances or ";
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
    if (d.get()->hasListAssignments(bin_type) && ((abundances.size() != 0) ||
        (samples.size() != 0))) {
        string message = "[ERROR]: You cannot assign abundance and sample data";
        message += " to bins that have sequence assignments. This could cause ";
        message += "inconsistencies.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }else{
        numBinsAssigned = d.get()->assignBins(bin_names, abundances, samples,
                                   sequence_names, bin_type);
    }

    string tag = " " + bin_type +" bins.";
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
//' @param bin_type a string indicating the type of bin assignments. Default "otu".
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
//'   assign_bin_representative_sequences(data = miseq,
//'                                       table = table,
//'                                       bin_type = "otu")
//'
//' @return double containing the number of representative sequences assigned
//[[Rcpp::export]]
double assign_bin_representative_sequences(Rcpp::Environment data,
                                           Rcpp::DataFrame table,
                                           string bin_type = "otu",
                                           Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                                           string bin_name = "bin_names",
                                           string sequence_name = "sequence_names") {

    vector<string> bin_names,  sequence_names;
    bin_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                  bin_name));
    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                       sequence_name));

    Rcpp::XPtr<Dataset> d = data["data"];

    double numAssigned = d.get()->assignBinRepresentativeSequences(bin_names,
             sequence_names, bin_type);

    string tag = " " + bin_type +" bin representative sequences.";
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
    bin_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                  bin_name));
    taxonomies = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
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

    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                       sequence_name));
    samples = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                "samples",
                                                                sample));
    treatments = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                "treatments",
                                                                treatment));
    abundances = Rcpp::as<vector<float>>(xint_fill_optional_parameters(table,
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
    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                       sequence_name));
    taxonomies = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
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
    samples = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                sample));
    treatments = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
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
//'   assign_bin_representative_sequences(data = miseq,
//'                                       table = table,
//'                                       bin_type = "otu")
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

