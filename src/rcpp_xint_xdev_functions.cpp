

#include "rcpp_xint_xdev_functions.h"

/******************************************************************************/
SEXP xint_fill_required_parameters(const Rcpp::DataFrame df,
                                   const string& given_column_name,
                                   string type) {

    if (!df.containsElementNamed(given_column_name.c_str())) {
        string message = "Expected a data.frame column named " +
            given_column_name +" to be provided.";

        throw Rcpp::exception(message.c_str());
    }

    return df[given_column_name];
}
/******************************************************************************/
SEXP xint_fill_optional_parameters(const Rcpp::DataFrame df,
                                   const string& default_column_name,
                                   const string& given_column_name,
                                   string type) {


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
    }

    return Rcpp::CharacterVector(0);
}
/******************************************************************************/
Rcpp::DataFrame xdev_abundance(Rcpp::Environment data,
                           string type,
                           string bin_type,
                           bool by_sample) {

     Rcpp::XPtr<Dataset> d = data["data"];

    if (type == "sequences") {
        // data.frame containing 2, 3 or 4 columns: sequence_names, abundances,
        //' samples (if assigned), and treatments (if assigned)
        return d.get()->getSequenceAbundances(by_sample);
    }
    else if (type == "bins") {
        // data.frame containing 2, 3 or 4 columns: bin_names, abundances,
        //' samples (if assigned), and treatments (if assigned)
        return d.get()->getBinAbundances(bin_type, by_sample);
    }
    else if ((type == "samples") || (type == "treatments")){
        return d.get()->getTotals(type);
    }
    else {
        string message = type + " is not a valid type for the abundance function";
        message += ". Types include: 'bins', 'sequences', samples' and 'treatments'.";
        RcppThread::Rcout << endl << message << endl;
    }

    // empty
    return Rcpp::DataFrame::create();
}
/******************************************************************************/
double xdev_add_references(Rcpp::Environment data,
                           Rcpp::DataFrame table,
                           string reference_name,
                           string reference_version,
                           string reference_usage,
                           string reference_note,
                           string reference_url,
                           bool verbose) {

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

    if (verbose) {
        Rcpp::Environment strollur_env("package:strollur");
        Rcpp::Function message = strollur_env["added_message"];
        message(numAdded, "resource references");
    }

    return numAdded;
}
/******************************************************************************/
void xdev_add_report(Rcpp::Environment data,
                     Rcpp::DataFrame table,
                     string type,
                     string sequence_name,
                     bool verbose) {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (type == "metadata") {
        d.get()->addMetadata(table);
        if (verbose) {
            Rcpp::Environment strollur_env("package:strollur");
            Rcpp::Function message = strollur_env["added_message"];
            message(R_NilValue, type);
        }
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
                    Rcpp::Environment strollur_env("package:strollur");
                    Rcpp::Function sort = strollur_env["sort_dataframe"];

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

            if (verbose) {
                Rcpp::Environment strollur_env("package:strollur");
                Rcpp::Function message = strollur_env["added_message"];
                message(R_NilValue, "a " + type);
            }
        }
    }
}
/******************************************************************************/
double xdev_add_sequences(Rcpp::Environment data,
                     Rcpp::DataFrame table,
                     Rcpp::Nullable<Rcpp::List> reference,
                     string sequence_name,
                     string sequence,
                     string comment,
                     bool verbose) {

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

    if (verbose) {
        Rcpp::Environment strollur_env("package:strollur");
        Rcpp::Function message = strollur_env["added_message"];
        message(numAdded);
    }

    return numAdded;
}
/******************************************************************************/
double xdev_assign_bins(Rcpp::Environment data,
                        Rcpp::DataFrame table,
                        string bin_type,
                        Rcpp::Nullable<Rcpp::List> reference,
                        string bin_name,
                        string abundance,
                        string sample,
                        string sequence_name,
                        bool verbose) {

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

    // if you have list assignments, don't allow setting bin abundances
    if (d.get()->hasListAssignments() && ((abundances.size() != 0) || (samples.size() != 0))) {
        string message = "[ERROR]: You cannot assign abundance and sample data";
        message += " to bins that have sequence assignments. This could cause ";
        message += "inconsistencies.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }else{
        numBinsAssigned = d.get()->assignBins(bin_names, abundances, samples,
                                sequence_names, bin_type);
    }

    if (verbose) {
        string tag = " " + bin_type +" bins.";
        Rcpp::Environment strollur_env("package:strollur");
        Rcpp::Function message = strollur_env["assigned_message"];
        message(numBinsAssigned, tag);
    }

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
double xdev_assign_bin_representative_sequences(Rcpp::Environment data,
                                                Rcpp::DataFrame table,
                                                string bin_type,
                                                Rcpp::Nullable<Rcpp::List> reference,
                                                string bin_name,
                                                string sequence_name,
                                                bool verbose) {

    vector<string> bin_names,  sequence_names;
    bin_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                       bin_name));
    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                            sequence_name));

    Rcpp::XPtr<Dataset> d = data["data"];

    double numAssigned = d.get()->assignBinRepresentativeSequences(bin_names,
                               sequence_names, bin_type);

    if (verbose) {
        string tag = " " + bin_type +" bin representative sequences.";
        Rcpp::Environment strollur_env("package:strollur");
        Rcpp::Function message = strollur_env["assigned_message"];
        message(numAssigned, tag);
    }

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
double xdev_assign_bin_taxonomy(Rcpp::Environment data,
                                Rcpp::DataFrame table,
                                string bin_type,
                                Rcpp::Nullable<Rcpp::List> reference,
                                string bin_name,
                                string taxonomy,
                                bool verbose) {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (d.get()->getNumBins(bin_type) == 0) {
        string message = "[ERROR]: No bin data for type " + bin_type + ", please ";
        message += " assign bins using the 'assign' function then try ";
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
                               taxonomies, bin_type);

    if (verbose) {
        string tag = " " + bin_type +" bin taxonomies.";
        Rcpp::Environment strollur_env("package:strollur");
        Rcpp::Function message = strollur_env["assigned_message"];
        message(numAssigned, tag);
    }

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
double xdev_assign_sequence_abundance(Rcpp::Environment data,
                                      Rcpp::DataFrame table,
                                      string sequence_name,
                                      string abundance,
                                      string sample,
                                      string treatment,
                                      bool verbose) {

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

    if (verbose) {
        string tag = " sequence abundances.";
        Rcpp::Environment strollur_env("package:strollur");
        Rcpp::Function message = strollur_env["assigned_message"];
        message(numAssigned, tag);
    }

    return numAssigned;
}
/******************************************************************************/
double xdev_assign_sequence_taxonomy(Rcpp::Environment data,
                                     Rcpp::DataFrame table,
                                     Rcpp::Nullable<Rcpp::List> reference,
                                     string sequence_name,
                                     string taxonomy,
                                     bool verbose) {

    vector<string> sequence_names, taxonomies;
    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                            sequence_name));
    taxonomies = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                        taxonomy));
    Rcpp::XPtr<Dataset> d = data["data"];

    double numAssigned = d.get()->assignSequenceTaxonomy(sequence_names,
                               taxonomies);

    if (verbose) {
        Rcpp::Environment strollur_env("package:strollur");
        Rcpp::Function message = strollur_env["assigned_message"];
        message(numAssigned, " sequence taxonomies.");
    }

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
double xdev_assign_treatments(Rcpp::Environment data,
                              Rcpp::DataFrame table,
                              string sample,
                              string treatment,
                              bool verbose) {

    vector<string> samples, treatments;
    samples = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                     sample));
    treatments = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                        treatment));

    Rcpp::XPtr<Dataset> d = data["data"];

    if (d.get()->getNumSamples() == 0) {
        string message = "[ERROR]: You cannot assign treatments, your dataset";
        message += " does not include sample data.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    // make sure every sample in dataset is assigned a treatment
    if (!identical(d.get()->getSamples(), unique(samples))) {
        string message = "You must provide treatment assignments for";
        message += " all samples in your dataset.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    double numAssigned = d.get()->assignTreatments(samples, treatments);

    if (verbose) {
        Rcpp::Environment strollur_env("package:strollur");
        Rcpp::Function message = strollur_env["assigned_message"];
        message(numAssigned, " samples to treatments.");
    }

    return numAssigned;
}
/******************************************************************************/
double xdev_count(Rcpp::Environment data,
                  string type,
                  string bin_type,
                  Rcpp::Nullable<Rcpp::List> samples,
                  bool distinct) {

    Rcpp::XPtr<Dataset> d = data["data"];

    vector<string> s;
    if (samples.isNotNull()) {
        s = Rcpp::as<vector<string>>(samples);
    }

    // types include "sequences", "samples", "treatments", "bins"

    if (type == "sequences") {
        // no sequence data and asked for samples. we can't find the number of
        // sequences in the given samples without sequence data because there is
        // not a correlation between otu sample abundance counts and # of seqs
        if (!d.get()->hasSequenceData && !s.empty()) {
            string message = "Unable to find the number of sequences in the ";
            message += "samples requested without sequence data.";
            throw Rcpp::exception(message.c_str());
        }

        if (distinct) {
            return d.get()->getUniqueTotal(s);
        }

        return d.get()->getTotal(s);
    }
    else if (type == "samples") {
        return d.get()->getNumSamples();
    }
    else if (type == "treatments") {
        return d.get()->getNumTreatments();
    }
    else if (type == "bins") {
        if (!s.empty()) {
            if (d.get()->hasSamples(s)) {
                return d.get()->getNumBins(bin_type,
                             s, distinct);
            }else {
                string message = "Your dataset does not include all the ";
                message += "samples requested, ignoring.";
                RcppThread::Rcout << endl << message << endl;
            }
        }else {
            return d.get()->getNumBins(bin_type);
        }
    }
    else if (type == "references") {
        return d.get()->getNumResourceReferences();
    }else{
        string message = "Invalid type. Types include: 'sequences', 'samples'";
        message += ", 'treatments' and 'bins'";
        throw Rcpp::exception(message.c_str());
    }

    return 0;
}
/******************************************************************************/
vector<vector<float> > xdev_get_abundances_by_sample(Rcpp::Environment data,
                                                     Rcpp::CharacterVector samples) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequenceAbundanceBySample(Rcpp::as<vector<string>>(samples));
}

/******************************************************************************/
vector<string> xdev_get_list_vector(Rcpp::Environment data,
                                    string type) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getListVector(type);
}
/******************************************************************************/
vector<vector<string> > xdev_get_by_sample(Rcpp::Environment data,
                                           string type,
                                           Rcpp::CharacterVector samples) {

    Rcpp::XPtr<Dataset> d = data["data"];

    if (type == "sequence_names") {
        return d.get()->getSequenceNamesBySample(Rcpp::as<vector<string>>(samples));
    }
    else if (type == "sequences") {
        return d.get()->getSequencesBySample(Rcpp::as<vector<string>>(samples));
    }
    else {
        string message = "Invalid type. Types include: 'sequence_names' and ";
        message += "'sequences'";
        throw Rcpp::exception(message.c_str());
    }
    return null2DVector;
}
/******************************************************************************/
vector<string> xdev_get_sequences(Rcpp::Environment data, string sample) {
    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequences(sample);
}
/******************************************************************************/
void xdev_merge_bins(Rcpp::Environment data, vector<string> bin_names,
                      string reason, string bin_type) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->mergeBins(bin_names, reason, bin_type);
}
/******************************************************************************/
void xdev_merge_sequences(Rcpp::Environment data, vector<string> sequence_names,
                           string reason) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->mergeSequences(sequence_names, reason);
}
/******************************************************************************/
const vector<string> xdev_names(Rcpp::Environment data,
                          string type,
                          string bin_type,
                          Rcpp::Nullable<Rcpp::List> samples,
                          bool distinct) {

    Rcpp::XPtr<Dataset> d = data["data"];

    vector<string> s;
    if (samples.isNotNull()) {
        s = Rcpp::as<vector<string>>(samples);
    }

    // types -> "dataset", "sequences", "bins", "samples",
    //            "treatments", "reports"

    vector<string> names;
    if (type == "sequences") {
        return d.get()->getSequenceNames(s, distinct);
    }
    else if (type == "samples") {
        return d.get()->getSamples();
    }
    else if (type == "treatments") {
        return d.get()->getTreatments();
    }
    else if (type == "bins") {
        if (!s.empty()) {
            if (d.get()->hasSamples(s)) {
                return d.get()->getBinIds(bin_type,
                             s, distinct);
            }else {
                string message = "Your dataset does not include all the ";
                message += "samples requested, ignoring.";
                RcppThread::Rcout << endl << message << endl;
            }
        }else {
            return d.get()->getBinIds(bin_type,
                         nullVector, distinct);
        }
    }
    else if (type == "reports") {
        return d.get()->getReportTypes();
    }
    else if (type == "dataset") {
        names.push_back(d.get()->datasetName);
    }else{
        string message = "Invalid type. Types include: 'dataset', 'sequences'";
        message += ", 'bins', 'samples', 'treatments' and 'reports'";
        throw Rcpp::exception(message.c_str());
    }

    return names;
}
/******************************************************************************/
void xdev_remove_bins(Rcpp::Environment data, vector<string> bin_names,
                       vector<string> trash_tags, string bin_type) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeBins(bin_names, trash_tags, bin_type);
}
/******************************************************************************/
void xdev_remove_lineages(Rcpp::Environment data, vector<string> contaminants,
                           string reason) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeLineages(contaminants, reason);
}
/******************************************************************************/
void xdev_remove_samples(Rcpp::Environment data, vector<string> samples,
                         string reason) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeSamples(samples, reason);
}
/******************************************************************************/
void xdev_remove_sequences(Rcpp::Environment data,
                            vector<string> sequence_names,
                            vector<string> trash_tags) {

     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeSequences(sequence_names, trash_tags);
}
/******************************************************************************/
Rcpp::DataFrame xdev_report(Rcpp::Environment data, string type,
                            string bin_type) {

    Rcpp::XPtr<Dataset> d = data["data"];

    // sequence_data reports contain the starts, ends, ambigs,...
    if (type == "sequences") {
        return d.get()->getSequenceReport();
    }
    // sequence fasta data
    else if (type == "fasta") {
        return d.get()->getFastaReport();
    }
    // sequence bin assignments report
    else if (type == "sequence_bin_assignments") {
        return d.get()->getList(bin_type);
    }
    // sample treatment assignments report
    else if (type == "sample_assignments") {
        return d.get()->getSampleTreatmentAssignments();
    }
    // representative sequences assignments report
    else if (type == "bin_representatives") {
        return d.get()->getBinRepresentativeSequences(bin_type);
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
void xdev_set_abundance(Rcpp::Environment data,
                         vector<string> sequence_names,
                         vector<float> sequence_abundances,
                         string reason) {

     Rcpp::XPtr<Dataset> d = data["data"];

     if (d.get()->getNumSamples() == 0) {
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
void xdev_set_abundances(Rcpp::Environment data,
                          vector<string> sequence_names,
                          vector<vector<float>> abundances,
                          string reason) {

     Rcpp::XPtr<Dataset> d = data["data"];

     if (d.get()->getNumSamples() != 0) {
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
void xdev_set_sequences(Rcpp::Environment data,
                         vector<string> sequence_names,
                         vector<string> sequences,
                         Rcpp::CharacterVector comments) {

     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->setSequences(sequence_names, sequences,
           Rcpp::as<vector<string>>(comments));
}
/******************************************************************************/
void xdev_set_dataset_name(Rcpp::Environment data, string dataset_name = "") {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->datasetName = dataset_name;
}
/******************************************************************************/
void xdev_set_num_processors(Rcpp::Environment data, int processors = 1) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->processors = processors;
}
/******************************************************************************/
Rcpp::DataFrame xdev_summarize(Rcpp::Environment data, string type,
                               Rcpp::Nullable<Rcpp::CharacterVector> report_type) {

    Rcpp::XPtr<Dataset> d = data["data"];

    string rType = "";
    if (report_type.isNotNull()) {
        rType = Rcpp::as<string>(report_type);

        // make sure its a valid report type
        vector<string> reportOptions = d.get()->getReportTypes();
        if (!vectorContains(reportOptions, rType)) {
            string message = rType + " is not a valid report_type option. ";
            if (!reportOptions.empty()) {
                message += "Options include: " + toString(reportOptions, ',') + ".";
            }
            RcppThread::Rcerr << endl << message << endl;
            throw Rcpp::exception(message.c_str());
        }
    }

    vector<string> typeOptions = {"sequences", "reports", "scrap"};
    if (!vectorContains(typeOptions, type)) {
        string message = type + " is not a valid type option. Options include:";
        message += " 'sequences', 'reports' and 'scrap'.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }


    return d.get()->getSummary(type, rType);
}
/******************************************************************************/
Rcpp::XPtr<Dataset> xint_copy_pointer(Rcpp::Environment data) {
     Rcpp::XPtr<Dataset> d = data["data"];
     Dataset* copy = new Dataset(*d.get());
     return Rcpp::XPtr<Dataset>(copy);
}
/******************************************************************************/
Rcpp::XPtr<Dataset> xint_new_pointer(string dataset_name = "", int processors = 1) {
     Dataset* d = new Dataset(dataset_name, processors);
     return Rcpp::XPtr<Dataset>(d);
}
/******************************************************************************/
void xint_deserialize_dobject(Rcpp::Environment data) {
     data["data"] = xint_new_pointer();
     Rcpp::XPtr<Dataset> d = data["data"];
     const Rcpp::RawVector raw = data["raw"];
     d.get()->loadFromSerialized(raw);
}
/******************************************************************************/
void xint_serialize_dobject(Rcpp::Environment data) {
     Rcpp::XPtr<Dataset> d = data["data"];
     data["raw"] = d.get()->serializeDataset();
}
/******************************************************************************/


