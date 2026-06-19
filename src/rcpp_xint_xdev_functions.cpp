

#include "rcpp_xint_xdev_functions.h"
/******************************************************************************/
void xint_added_message(double num, string tag) {
    string message = "Added ";
    if (num != -1) {
        message += toString(num) + " " + tag + ".";
    } else {
        message += tag + ".";
    }

    Rcpp::Rcout << message << endl;
}
/******************************************************************************/
void xint_assigned_message(double num, string tag) {
    string message = "Assigned " + toString(num) + tag;
    Rcpp::Rcout << message << endl;
}
/******************************************************************************/
void xint_updated_message(double num, string tag) {
    string message = "Updated ";
    if (num != -1) {
        message += toString(num) + " " + tag + ".";
    } else {
        message += tag + ".";
    }

    Rcpp::Rcout << message << endl;
}
/******************************************************************************/
void fillReference(Reference& ref, Rcpp::List ref_list) {
    ref.name = Rcpp::as<string>(ref_list["name"]);
    ref.vendor = Rcpp::as<string>(ref_list["vendor"]);
    ref.version = Rcpp::as<string>(ref_list["version"]);
    ref.usage = Rcpp::as<string>(ref_list["usage"]);
    ref.note = Rcpp::as<string>(ref_list["note"]);
    ref.method_url = Rcpp::as<string>(ref_list["method_url"]);
    ref.documentation_url = Rcpp::as<string>(ref_list["documentation_url"]);
    ref.parameter = Rcpp::as<string>(ref_list["parameter"]);
    ref.citation = Rcpp::as<string>(ref_list["citation"]);
    ref.checks();
}
/******************************************************************************/
SEXP xint_fill_required_parameters(const Rcpp::DataFrame& df,
                                   const string& given_column_name,
                                   const string& type) {

    if (!df.containsElementNamed(given_column_name.c_str())) {
        const string message = "Expected a data.frame column named " +
            given_column_name +" to be provided.";

        throw Rcpp::exception(message.c_str());
    }

    return df[given_column_name];
}
/******************************************************************************/
SEXP xint_fill_optional_parameters(const Rcpp::DataFrame& df,
                                   const string& default_column_name,
                                   const string& given_column_name,
                                   const string& type) {


    // Check if the column exists in the data frame
    if (df.containsElementNamed(given_column_name.c_str())) {
        return df[given_column_name];
    }else if (df.containsElementNamed(default_column_name.c_str())) {
        return df[default_column_name];
    }else if ((given_column_name == "") ||
        (given_column_name == default_column_name)) {
        //fall through
    }else {
        const string message = "Expected a data.frame column named " +
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
Rcpp::DataFrame xdev_abundance(const Rcpp::Environment& data,
                           const string& type,
                           const string& bin_type,
                           const bool by_sample) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];

    if (type == "sequence") {
        // data.frame containing 2, 3 or 4 columns: sequence_names, abundances,
        //' samples (if assigned), and treatments (if assigned)
        return d.get()->getSequenceAbundances(by_sample);
    }
    else if (type == "bin") {
        // data.frame containing 2, 3 or 4 columns: bin_names, abundances,
        //' samples (if assigned), and treatments (if assigned)
        return d.get()->getBinAbundances(bin_type, by_sample);
    }
    else if ((type == "sample") || (type == "treatment")){
        return d.get()->getTotals(type);
    }
    else {
        string message = type + " is not a valid type for the abundance function";
        message += ". Types include: 'bin', 'sequence', sample' and 'treatment'.";
        Rcpp::Rcout << endl << message << endl;
    }

    // empty
    return Rcpp::DataFrame::create();
}
/******************************************************************************/
Rcpp::Environment xdev_add_references(const Rcpp::Environment& data,
                           const Rcpp::DataFrame& table,
                           const string& name,
                           const string& vendor,
                           const string& version,
                           const string& usage,
                           const string& note,
                           const string& method_url,
                           const string& documentation_url,
                           const string& parameter,
                           const string& citation,
                           bool verbose) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    vector<string> reference_vendors, reference_names, reference_versions;
    vector<string> reference_usages, reference_notes, reference_method_urls, reference_urls;
    vector<string> reference_parameters, reference_citations;

    reference_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                             name));
    reference_vendors = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                                "vendor",
                                                                                vendor));
    reference_versions = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                                "version",
                                                                                version));
    reference_usages = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                              "usage",
                                                                              usage));
    reference_notes = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                              "note",
                                                                              note));
    reference_method_urls = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                             "method_url",
                                                                             method_url));
    reference_urls = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                            "documentation_url",
                                                                            documentation_url));

    reference_parameters = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                    "parameter",
                                                    parameter));

    reference_citations = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table, "citation",
                                                                                citation));
    bool hasParameters = false, hasVersions = false, hasUsages = false, hasNotes = false;
    bool hasMethods = false, hasUrls = false, hasVendors = false, hasCitations = false;

    if (reference_vendors.size() == reference_names.size()) {
        hasVendors = true;
    }
    if (reference_versions.size() == reference_names.size()) {
        hasVersions = true;
    }
    if (reference_method_urls.size() == reference_names.size()) {
        hasMethods = true;
    }
    if (reference_usages.size() == reference_names.size()) {
        hasUsages = true;
    }
    if (reference_notes.size() == reference_names.size()) {
        hasNotes = true;
    }
    if (reference_method_urls.size() == reference_names.size()) {
        hasMethods = true;
    }
    if (reference_urls.size() == reference_names.size()) {
        hasUrls = true;
    }
    if (reference_parameters.size() == reference_names.size()) {
        hasParameters = true;
    }
    if (reference_citations.size() == reference_names.size()) {
        hasCitations = true;
    }

    vector<Reference> refs;
    for (size_t i = 0; i < reference_names.size(); i++) {

        Reference ref(reference_names[i]);
        if (hasVendors) {
            ref.vendor = reference_vendors[i];
        }
        if (hasVersions) {
            ref.version = reference_versions[i];
        }
        if (hasUsages) {
            ref.usage = reference_usages[i];
        }
        if (hasNotes) {
            ref.note = reference_notes[i];
        }
        if (hasMethods) {
            ref.method_url = reference_method_urls[i];
        }
        if (hasUrls) {
            ref.documentation_url = reference_urls[i];
        }
        if (hasParameters) {
            ref.parameter = reference_parameters[i];
        }
        if (hasCitations) {
            ref.citation = reference_citations[i];
        }
        ref.checks();
        refs.push_back(ref);
    }

    Rcpp::XPtr<Dataset> d = data["data"];

    double numAdded = d.get()->addReferences(refs);

    if (verbose) {
        if (numAdded == 0) {
            xint_updated_message(-1, "resource references");
        }else {
            xint_added_message(numAdded, "resource references");
        }
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_add_report(const Rcpp::Environment& data,
                     Rcpp::DataFrame table,
                     const string& type,
                     const string& sequence_name,
                     bool verbose) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    Rcpp::XPtr<Dataset> d = data["data"];

    if (sequence_name != "none") {
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
                d.get()->addSequences(sequenceNames,
                      sequences, comments);
            }else {
                // we have sequences already, make sure there is a report row for
                // each sequence in dataset
                vector<string> datasetSeqNames = d.get()->getSequenceNames();

                // find seqs in dataset and not in report
                vector<string> missingSeqs = setDiff(datasetSeqNames, sequenceNames);

                if (missingSeqs.size() == 0) {
                    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("strollur");
                    Rcpp::Function sort = pkg["sort_dataframe"];
                    table = sort(table, datasetSeqNames, sequence_name);
                }else {
                    string message = "Your report does not contain an entry for ";
                    message += "every sequence in your dataset, ignoring report. ",
                        Rcpp::Rcout << endl << message << endl;
                    return data;
                }
            }

            // save name column
            table.attr("sequence_name") = sequence_name;

            // add report to dataset
            d.get()->addReport(table, type);
        }
    }else {
        //add generic non sequence report
        d.get()->addReport(table, type);
    }

    if (verbose) {
        xint_added_message(-1, "a " + type + " report");
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_add_sequences(const Rcpp::Environment& data,
                     const Rcpp::DataFrame& table,
                     Rcpp::Nullable<Rcpp::List> reference,
                     const string& sequence_name,
                     const string& sequence,
                     const string& comment,
                     const bool verbose) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    vector<string> sequence_names, sequences, comments;
    sequence_names = Rcpp::as<vector<string>>(
        xint_fill_required_parameters(table, sequence_name));
    sequences = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                       "sequence",
                                                                       sequence));
    comments = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                      "comment",
                                                                      comment));

    Rcpp::XPtr<Dataset> d = data["data"];

    const double numAdded = d.get()->addSequences(sequence_names,
                                  sequences, comments);

    if (verbose) {
        xint_added_message(numAdded);
    }

    Reference ref;
    if (reference.isNotNull()) {

        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);
        fillReference(ref, ref_list);

        vector<Reference> refs;
        refs.push_back(ref);

        double numAdded = d.get()->addReferences(refs);

        if (verbose) {
            if (numAdded == 0) {
                xint_updated_message(1, "resource references");
            }else {
                xint_added_message(1, "resource references");
            }
        }
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_assign_bins(const Rcpp::Environment& data,
                        const Rcpp::DataFrame& table,
                        const string& bin_type,
                        Rcpp::Nullable<Rcpp::List> reference,
                        const string& bin_name,
                        const string& abundance,
                        const string& sample,
                        const string& sequence_name,
                        bool verbose) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    vector<string> bin_names, samples, sequence_names;
    vector<float> abundances;

    // fill vectors with columns from table
    bin_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                       bin_name));
    samples = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                     "sample",
                                                                     sample));
    sequence_names = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                            "sequence_name",
                                                                            sequence_name));
    abundances = Rcpp::as<vector<float>>(xint_fill_optional_parameters(table,
                                                                       "abundance",
                                                                       abundance,
                                                                       "float"));

    if ((abundances.empty()) && (sequence_names.empty())) {
        string message = "You must provide either abundances or ";
        message += "sequence_names to assign bins.";
        throw Rcpp::exception(message.c_str());
    }

    double numBinsAssigned = 0;
    Rcpp::XPtr<Dataset> d = data["data"];

    // if you have list assignments, don't allow setting bin abundances
    if (d.get()->hasListAssignments() && ((abundances.size() != 0) || (samples.size() != 0))) {
        string message = "[ERROR]: You cannot assign abundance and sample data";
        message += " to bins that have sequence assignments. This could cause ";
        message += "inconsistencies.";
        throw Rcpp::exception(message.c_str());
    }else{
        numBinsAssigned = d.get()->assignBins(bin_names, abundances, samples,
                                sequence_names, bin_type);
    }

    if (verbose) {
        string tag = " " + bin_type +" bins.";
        xint_assigned_message(numBinsAssigned, tag);
    }

    if (reference.isNotNull()) {
        Reference ref;
        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        fillReference(ref, ref_list);

        vector<Reference> refs;
        refs.push_back(ref);

        double numAdded = d.get()->addReferences(refs);

        if (verbose) {
            if (numAdded == 0) {
                xint_updated_message(1, "resource references");
            }else {
                xint_added_message(1, "resource references");
            }
        }
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_assign_bin_representative_sequences(const Rcpp::Environment& data,
                                                const Rcpp::DataFrame& table,
                                                const string& bin_type,
                                                Rcpp::Nullable<Rcpp::List> reference,
                                                const string& bin_name,
                                                const string& sequence_name,
                                                bool verbose) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    vector<string> bin_names,  sequence_names;
    bin_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                       bin_name));
    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                            sequence_name));

    const Rcpp::XPtr<Dataset> d = data["data"];

    const double numAssigned = d.get()->assignBinRepresentativeSequences(bin_names,
                               sequence_names, bin_type);

    if (verbose) {
        const string tag = " " + bin_type +" bin representative sequences.";
        xint_assigned_message(numAssigned, tag);
    }


    if (reference.isNotNull()) {
        Reference ref;
        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        fillReference(ref, ref_list);

        vector<Reference> refs;
        refs.push_back(ref);

        d.get()->addReferences(refs);

        if (verbose) {
            xint_added_message(1, "resource references");
        }
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_assign_bin_taxonomy(const Rcpp::Environment& data,
                                const Rcpp::DataFrame& table,
                                const string& bin_type,
                                Rcpp::Nullable<Rcpp::List> reference,
                                const string& bin_name,
                                const string& taxonomy,
                                bool verbose) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    Rcpp::XPtr<Dataset> d = data["data"];

    if (d.get()->getNumBins(bin_type) == 0) {
        string message = "[ERROR]: No bin data for type " + bin_type + ", please ";
        message += " assign bins using the 'assign' function then try ";
        message += "again.";
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
        xint_assigned_message(numAssigned, tag);
    }


    if (reference.isNotNull()) {
        Reference ref;
        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        fillReference(ref, ref_list);

        vector<Reference> refs;
        refs.push_back(ref);

        d.get()->addReferences(refs);

        if (verbose) {
            xint_added_message(1, "resource references");
        }
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_assign_sequence_abundance(const Rcpp::Environment& data,
                                      const Rcpp::DataFrame& table,
                                      const string& sequence_name,
                                      const string& abundance,
                                      const string& sample,
                                      const string& treatment,
                                      bool verbose) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    vector<string> sequence_names, samples, treatments;
    vector<float> abundances;

    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                            sequence_name));
    samples = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                     "sample",
                                                                     sample));
    treatments = Rcpp::as<vector<string>>(xint_fill_optional_parameters(table,
                                                                        "treatment",
                                                                        treatment));
    abundances = Rcpp::as<vector<float>>(xint_fill_optional_parameters(table,
                                                                       "abundance",
                                                                       abundance,
                                                                       "float"));
    if (sequence_names.size() != abundances.size()) {
        string message = "[ERROR]: The names and abundances must be the same";
        message += " length.";
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
            throw Rcpp::exception(message.c_str());
        }
    }

    double numAssigned = d.get()->assignSequenceAbundance(sequence_names,
                               abundances,
                               samples, treatments);

    if (verbose) {
        string tag = " sequence abundances.";
        xint_assigned_message(numAssigned, tag);
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_assign_sequence_taxonomy(const Rcpp::Environment& data,
                                     const Rcpp::DataFrame& table,
                                     Rcpp::Nullable<Rcpp::List> reference,
                                     const string& sequence_name,
                                     const string& taxonomy,
                                     const bool verbose) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    vector<string> sequence_names, taxonomies;
    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                            sequence_name));
    taxonomies = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                        taxonomy));
    const Rcpp::XPtr<Dataset> d = data["data"];

    const double numAssigned = d.get()->assignSequenceTaxonomy(sequence_names,
                               taxonomies);

    if (verbose) {
        xint_assigned_message(numAssigned, " sequence taxonomies.");
    }


    if (reference.isNotNull()) {

        Reference ref;
        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        fillReference(ref, ref_list);

        vector<Reference> refs;
        refs.push_back(ref);

        d.get()->addReferences(refs);

        if (verbose) {
            xint_added_message(1, "resource references");
        }
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_assign_sequence_taxonomy_tidy(const Rcpp::Environment& data,
                                          const Rcpp::DataFrame& table,
                                          Rcpp::Nullable<Rcpp::List> reference,
                                          const string& sequence_name,
                                          const string& level,
                                          const string& taxonomy,
                                          const string& confidence,
                                          const bool verbose){
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    vector<string> sequence_names, taxonomies;
    vector<float> confidences;
    vector<int> levels;

    sequence_names = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                            sequence_name));
    taxonomies = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                        taxonomy));
    levels = Rcpp::as<vector<int>>(xint_fill_required_parameters(table,
                                                                 level));
    confidences = Rcpp::as<vector<float>>(xint_fill_required_parameters(table,
                                                                        confidence));

    const Rcpp::XPtr<Dataset> d = data["data"];

    const double numAssigned = d.get()->assignSequenceTaxonomyTidy(sequence_names,
                                     levels, taxonomies, confidences);

    if (verbose) {
        xint_assigned_message(numAssigned, " sequence taxonomies.");
    }


    if (reference.isNotNull()) {

        Reference ref;
        Rcpp::List ref_list = Rcpp::as<Rcpp::List>(reference);

        fillReference(ref, ref_list);

        vector<Reference> refs;
        refs.push_back(ref);

        d.get()->addReferences(refs);

        if (verbose) {
            xint_added_message(1, "resource references");
        }
    }

    return data;
}
/******************************************************************************/
Rcpp::Environment xdev_assign_treatments(const Rcpp::Environment& data,
                              const Rcpp::DataFrame& table,
                              const string& sample,
                              const string& treatment,
                              bool verbose) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    vector<string> samples, treatments;
    samples = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                     sample));
    treatments = Rcpp::as<vector<string>>(xint_fill_required_parameters(table,
                                                                        treatment));

    const Rcpp::XPtr<Dataset> d = data["data"];

    if (d.get()->getNumSamples() == 0) {
        string message = "[ERROR]: You cannot assign treatments, your dataset";
        message += " does not include sample data.";
        throw Rcpp::exception(message.c_str());
    }

    // make sure every sample in dataset is assigned a treatment
    if (!identical(d.get()->getSamples(), unique(samples))) {
        string message = "You must provide treatment assignments for";
        message += " all samples in your dataset.";
        throw Rcpp::exception(message.c_str());
    }

    const double numAssigned = d.get()->assignTreatments(samples, treatments);

    if (verbose) {
        xint_assigned_message(numAssigned, " samples to treatments.");
    }

    return data;
}
/******************************************************************************/
double xdev_count(const Rcpp::Environment& data,
                  const string& type,
                  const string& bin_type,
                  Rcpp::Nullable<Rcpp::List> samples,
                  bool distinct) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];

    vector<string> s;
    if (samples.isNotNull()) {
        s = Rcpp::as<vector<string>>(samples);
    }

    // types include "sequence", "sample", "treatment", "bin"

    if (type == "sequence") {
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
    else if (type == "sample") {
        return d.get()->getNumSamples();
    }
    else if (type == "treatment") {
        return d.get()->getNumTreatments();
    }
    else if (type == "bin") {
        if (!s.empty()) {
            if (d.get()->hasSamples(s)) {
                return d.get()->getNumBins(bin_type,
                             s, distinct);
            }else {
                string message = "Your dataset does not include all the ";
                message += "samples requested, ignoring.";
                Rcpp::Rcout << endl << message << endl;
            }
        }else {
            return d.get()->getNumBins(bin_type);
        }
    }
    else if (type == "resource_reference") {
        return d.get()->getNumResourceReferences();
    }else{
        string message = "Invalid type. Types include: 'sequence', 'sample'";
        message += ", 'treatment', 'bin' and 'resource_reference'.";
        throw Rcpp::exception(message.c_str());
    }

    return 0;
}
/******************************************************************************/
Rcpp::List xdev_export_dataset(const Rcpp::Environment& data) {

    Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->exportDataset();
}
/******************************************************************************/
vector<vector<float> > xdev_get_abundances_by_sample(const Rcpp::Environment& data,
                                                     const Rcpp::CharacterVector& samples) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];

    return d.get()->getSequenceAbundanceBySample(Rcpp::as<vector<string>>(samples));
}

/******************************************************************************/
vector<string> xdev_get_list_vector(const Rcpp::Environment& data,
                                    const string& type) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getListVector(type);
}
/******************************************************************************/
vector<vector<string> > xdev_get_by_sample(const Rcpp::Environment& data,
                                           const string& type,
                                           const Rcpp::CharacterVector& samples,
                                           bool degap) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];

    if (type == "sequence_name") {
        return d.get()->getSequenceNamesBySample(Rcpp::as<vector<string>>(samples));
    }
    else if (type == "sequence") {
        return d.get()->getSequencesBySample(Rcpp::as<vector<string>>(samples),
                                             degap);
    }
    else {
        string message = "Invalid type. Types include: 'sequence_name' and ";
        message += "'sequence'";
        throw Rcpp::exception(message.c_str());
    }
    return null2DVector;
}
/******************************************************************************/
vector<string> xdev_get_sequences(const Rcpp::Environment& data,
                                  const string& sample,
                                  bool degap) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getSequences(sample, degap);
}
/******************************************************************************/
bool xdev_has_sequence_taxonomy(const Rcpp::Environment& data) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->hasSequenceTaxonomy;
}
/******************************************************************************/
Rcpp::Environment xdev_merge_bins(const Rcpp::Environment& data, const vector<string>& bin_names,
                     const string& reason, const string& bin_type) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->mergeBins(bin_names, reason, bin_type);

     return data;
}
/******************************************************************************/
Rcpp::Environment xdev_merge_sequences(const Rcpp::Environment& data,
                          const vector<string>& sequence_names,
                          const string& reason) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

     const Rcpp::XPtr<Dataset> d = data["data"];

     d.get()->mergeSequences(sequence_names, reason);

     return data;
}
/******************************************************************************/
vector<string> xdev_names(const Rcpp::Environment& data,
                          const string& type,
                          const string& bin_type,
                          Rcpp::Nullable<Rcpp::List> samples,
                          const bool distinct) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];

    vector<string> s;
    if (samples.isNotNull()) {
        s = Rcpp::as<vector<string>>(samples);
    }

    // types -> "dataset", "sequence", "bin", "sample",
    //            "treatment", "report"

    vector<string> names;
    if (type == "sequence") {
        return d.get()->getSequenceNames(s, distinct);
    }
    else if (type == "sample") {
        return d.get()->getSamples();
    }
    else if (type == "treatment") {
        return d.get()->getTreatments();
    }
    else if (type == "bin") {
        if (!s.empty()) {
            if (d.get()->hasSamples(s)) {
                return d.get()->getBinIds(bin_type,
                             s, distinct);
            }else {
                string message = "Your dataset does not include all the ";
                message += "samples requested, ignoring.";
                Rcpp::Rcout << endl << message << endl;
            }
        }else {
            return d.get()->getBinIds(bin_type,
                         nullVector, distinct);
        }
    }
    else if (type == "report") {
        return d.get()->getReportTypes();
    }
    else if (type == "dataset") {
        names.push_back(d.get()->datasetName);
    }else{
        string message = "Invalid type. Types include: 'dataset', 'sequence'";
        message += ", 'bin', 'sample', 'treatment' and 'report'";
        throw Rcpp::exception(message.c_str());
    }

    return names;
}
/******************************************************************************/
Rcpp::Environment xdev_remove_bins(const Rcpp::Environment& data, const vector<string>& bin_names,
                      const vector<string>& trash_tags, const string& bin_type) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const  Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeBins(bin_names, trash_tags, bin_type);
     return data;
}
/******************************************************************************/
Rcpp::Environment xdev_remove_lineages(const Rcpp::Environment& data,
                          const vector<string>& contaminants,
                          const string& reason) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeLineages(contaminants, reason);
     return data;
}
/******************************************************************************/
Rcpp::Environment xdev_remove_samples(const Rcpp::Environment& data,
                         const vector<string>& samples,
                         const string& reason) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeSamples(samples, reason);
     return data;
}
/******************************************************************************/
Rcpp::Environment xdev_remove_sequences(const Rcpp::Environment& data,
                           const vector<string>& sequence_names,
                           const vector<string>& trash_tags) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

     const Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->removeSequences(sequence_names, trash_tags);
     return data;
}
/******************************************************************************/
Rcpp::DataFrame xdev_report(const Rcpp::Environment& data, const string& type,
                            const string& bin_type) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    const Rcpp::XPtr<Dataset> d = data["data"];

    // sequence_data reports contain the starts, ends, ambigs,...
    if (type == "sequence") {
        return d.get()->getSequenceReport();
    }
    // sequence fasta data
    else if (type == "fasta") {
        return d.get()->getFastaReport();
    }
    // sequence bin assignments report
    else if (type == "sequence_bin_assignment") {
        return d.get()->getList(bin_type);
    }
    // sample treatment assignments report
    else if (type == "sample_assignment") {
        return d.get()->getSampleTreatmentAssignments();
    }
    // representative sequences assignments report
    else if (type == "bin_representative") {
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
    // references
    else if (type == "resource_reference") {
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
Rcpp::Environment xdev_set_abundance(const Rcpp::Environment& data,
                        const vector<string>& sequence_names,
                        const vector<float>& sequence_abundances,
                        const string& reason) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

     const Rcpp::XPtr<Dataset> d = data["data"];

     if (d.get()->getNumSamples() == 0) {
         d.get()->setAbundance(sequence_names, sequence_abundances, reason);
     }else{
         string message = "[ERROR]: You cannot set the total sequence abundance";
         message += " for sequences whose abundances are parsed by sample. ";
         message += "Try 'set_abundances' instead of 'set_abundance'.";
         throw Rcpp::exception(message.c_str());
     }
     return data;
}
/******************************************************************************/
Rcpp::Environment xdev_set_abundances(const Rcpp::Environment& data,
                          const vector<string>& sequence_names,
                          const vector<vector<float>>& abundances,
                          const string& reason) {

    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

     const Rcpp::XPtr<Dataset> d = data["data"];

     if (d.get()->getNumSamples() != 0) {
         d.get()->setAbundances(sequence_names, abundances, reason);
     }else {
         string message = "[ERROR]: You cannot set parsed sequence abundances ";
         message += "when your dataset does not include sample data. ";
         message += "Try 'set_abundance' instead of 'set_abundances'.";
         throw Rcpp::exception(message.c_str());
     }
     return data;
}
/******************************************************************************/
Rcpp::Environment xdev_set_sequences(const Rcpp::Environment& data,
                         const vector<string>& sequence_names,
                         const vector<string>& sequences,
                         const Rcpp::CharacterVector& comments) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

     const Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->setSequences(sequence_names, sequences,
           Rcpp::as<vector<string>>(comments));
     return data;
}
/******************************************************************************/
void xdev_set_dataset_name(const Rcpp::Environment& data, const string& dataset_name = "") {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }
     Rcpp::XPtr<Dataset> d = data["data"];
     d.get()->datasetName = dataset_name;
}
/******************************************************************************/
Rcpp::DataFrame xint_get_scrap_summary(const Rcpp::Environment& data) {
    if (!data.inherits("strollur")) {
        const string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }
    const Rcpp::XPtr<Dataset> d = data["data"];
    return d.get()->getScrapSummary();
}
/******************************************************************************/
Rcpp::XPtr<Dataset> xint_copy_pointer(const Rcpp::Environment& data) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

     const Rcpp::XPtr<Dataset> d = data["data"];
     Dataset* copy = new Dataset(*d.get());
     return Rcpp::XPtr<Dataset>(copy);
}
/******************************************************************************/
Rcpp::XPtr<Dataset> xint_new_pointer(const string& dataset_name = "") {
     Dataset* d = new Dataset(dataset_name);
     return Rcpp::XPtr<Dataset>(d);
}
/******************************************************************************/
bool xint_is_equal(Rcpp::Environment data, Rcpp::Environment data2) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }
    if (!data2.inherits("strollur")) {
        string message = "data2 must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }

    Rcpp::XPtr<Dataset> d = data["data"];
    Rcpp::XPtr<Dataset> d2 = data2["data"];

    return d.get()->isEqual(*d2.get());
}
/******************************************************************************/
void xint_deserialize_dobject(Rcpp::Environment data) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }
     data["data"] = xint_new_pointer();
     const Rcpp::XPtr<Dataset> d = data["data"];
     const Rcpp::RawVector raw = data["raw"];
     d.get()->loadFromSerialized(raw);
}
/******************************************************************************/
Rcpp::RawVector xint_serialize_dobject(Rcpp::Environment data) {
    if (!data.inherits("strollur")) {
        string message = "data must be a strollur object.";
        throw Rcpp::exception(message.c_str());
    }
     Rcpp::XPtr<Dataset> d = data["data"];
     Rcpp::RawVector raw = d.get()->serializeDataset();
     data["raw"] = raw;
     return raw;
}
/******************************************************************************/


