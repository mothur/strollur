

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
    }else{
        string message = "Unsupported column type for conversion.";
        throw Rcpp::exception(message.c_str());
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
    }else{
        string message = "Invalid type. Types include: 'sequences', 'samples'";
        message += ", 'treatments' and 'bins'";
        throw Rcpp::exception(message.c_str());
    }

    return 0;
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
void xdev_merge_bins(Rcpp::Environment data, vector<string> bin_names,
                      string reason, string type) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->mergeBins(bin_names, reason, type);
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
                       vector<string> trash_tags, string type) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeBins(bin_names, trash_tags, type);
}
/******************************************************************************/
void xdev_remove_lineages(Rcpp::Environment data, vector<string> contaminants,
                           string trash_tag) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeLineages(contaminants, trash_tag);
}
/******************************************************************************/
void xdev_remove_samples(Rcpp::Environment data, vector<string> samples) {
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeSamples(samples);
}
/******************************************************************************/
void xdev_remove_sequences(Rcpp::Environment data,
                            vector<string> sequence_names,
                            vector<string> trash_tags) {

     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeSequences(sequence_names, trash_tags);
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
void xdev_set_bin_abundance(Rcpp::Environment data,
                             vector<string> bin_names,
                             vector<float> abundances,
                             string type, string reason) {

     Rcpp::XPtr<Dataset> d = data["data"];

     if (d.get()->hasListAssignments(type)) {
         string message = "[ERROR]: You cannot set the bin abundance for bin ";
         message += "clusters with sequence assignments.";
         RcppThread::Rcerr << endl << message << endl;
         throw Rcpp::exception(message.c_str());

     }else if (d.get()->getNumSamples() != 0) {
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
void xdev_set_bin_abundances(Rcpp::Environment data,
                              vector<string> bin_names,
                              vector<vector<float>> abundances,
                              string type, string reason) {

     Rcpp::XPtr<Dataset> d = data["data"];

     if (d.get()->hasListAssignments(type)) {
         string message = "[ERROR]: You cannot set the bin abundance for bin ";
         message += "clusters with sequence assignments.";
         RcppThread::Rcerr << endl << message << endl;
         throw Rcpp::exception(message.c_str());
     }else if (d.get()->getNumSamples() == 0) {
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


