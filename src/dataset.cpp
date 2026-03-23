
#include "../inst/include/strollur.h"
#include "seqreport.h"
#include "summary.h"
#include "dataset.h"
#include "utils.h"

/******************************************************************************/
Dataset::Dataset() {
    datasetName = "";
    isAligned = false;
    hasSequenceData = false;
    hasSequenceTaxonomy = false;
    hasList = false;
    numUnique = 0;
    uniqueBad = 0;
    alignmentLength = 0;
    processors = 1;
}
/******************************************************************************/
Dataset::Dataset(string n, int proc) : datasetName(n) {
    isAligned = false;
    hasSequenceData = false;
    hasSequenceTaxonomy = false;
    hasList = false;
    numUnique = 0;
    uniqueBad = 0;
    alignmentLength = 0;
    processors = proc;
}
/******************************************************************************/
Dataset::Dataset(const Dataset& dataset) {

    // public
    datasetName = dataset.datasetName;
    isAligned = dataset.isAligned;
    alignmentLength = dataset.alignmentLength;
    hasSequenceData = dataset.hasSequenceData;
    hasSequenceTaxonomy = dataset.hasSequenceTaxonomy;
    numUnique = dataset.numUnique;
    processors = dataset.processors;

    // private
    names = dataset.names;
    seqs = dataset.seqs;
    comments = dataset.comments;
    trashCodes = dataset.trashCodes;

    starts = dataset.starts;
    ends = dataset.ends;
    lengths = dataset.lengths;
    ambigs = dataset.ambigs;
    polymers = dataset.polymers;
    numns = dataset.numns;

    taxonomies = dataset.taxonomies;
    seqIndex = dataset.seqIndex;
    tableSeqs = dataset.tableSeqs;

    badAccnos = dataset.badAccnos;
    uniqueBad = dataset.uniqueBad;
    hasList = dataset.hasList;

    count.clone(dataset.count);

    for (int i = 0; i < dataset.binTables.size(); i++) {
        BinTable table;
        table.clone(dataset.binTables[i]);
        binTables.push_back(table);
    }

    metadata = dataset.metadata;
    reports = dataset.reports;
}
/******************************************************************************/
void Dataset::loadFromSerialized(Rcpp::RawVector serializedDataset) {
    std::string serialized_data(reinterpret_cast<const char*>(serializedDataset.begin()),
                                serializedDataset.size());
    std::stringstream ss(serialized_data);
    {
        cereal::BinaryInputArchive iarchive(ss);
        iarchive(*this);
    }
}
/******************************************************************************/
Dataset::~Dataset() {}
/******************************************************************************/
void Dataset::clear() {
        isAligned = false;
        hasSequenceData = false;
        hasSequenceTaxonomy = false;
        hasList = false;
        numUnique = 0;
        uniqueBad = 0;
        alignmentLength = 0;

        // sequence data
        names.clear();
        seqs.clear();
        comments.clear();
        trashCodes.clear();

        // sequence summary data
        starts.clear();
        ends.clear();
        lengths.clear();
        ambigs.clear();
        polymers.clear();
        numns.clear();

        // sequence taxonomy assignments
        taxonomies.clear();

        // maps sequence name to index in vectors
        seqIndex.clear();
        tableSeqs.clear();

        badAccnos.clear();
        count.clear();
        binTables.clear();
        references.clear();
        reports.clear();
        metadata.clear();
}
/******************************************************************************/
Rcpp::List Dataset::exportDataset(vector<string> tags){

    Rcpp::List results = Rcpp::List::create();
    vector<string> resultsLabels;

    set<string> t;
    bool hasTags = false;
    if (tags.size() != 0) {
        t = toSet(tags);
        hasTags = true;

        // check for tags for data you dont have and warn
        if (!hasSequenceData) {

            if (setContains(t, "sequence_data")) {
                string message = "[WARNING]: The dataset does not include ";
                message += " sequence data, ignoring 'sequence_data' tag.";
                RcppThread::Rcout << endl << message << endl;
            }
        }

        if (binTables.empty()) {

            if (setContains(t, "bin_data")) {
                string message = "[WARNING]: The dataset does not include ";
                message += " bin data, ignoring 'bin_data' tag.";
                RcppThread::Rcout << endl << message << endl;
            }
        }

        if (!metadata.hasReport) {

            if (setContains(t, "metadata")) {
                string message = "[WARNING]: The dataset does not include ";
                message += " metadata, ignoring 'metadata' tag.";
                RcppThread::Rcout << endl << message << endl;
            }
        }

        if (reports.empty()) {
            if (setContains(t, "reports")) {
                string message = "[WARNING]: The dataset does not include ";
                message += " reports, ignoring 'reports' tag.";
                RcppThread::Rcout << endl << message << endl;
            }
        }
    }

    if (hasSequenceData) {

        if (!hasTags || setContains(t, "sequence_data")) {

            // sequence data.frame
            // ids, names, seqs, comments(optional),
            // trashCodes, taxonomies(optional), tableSeqs
            Rcpp::DataFrame sequenceData = Rcpp::DataFrame::create();
            vector<string> sequenceDataLabels;

            sequenceData.push_back(getIndexes(names));
            sequenceDataLabels.push_back("sequence_ids");
            sequenceData.push_back(names);
            sequenceDataLabels.push_back("sequence_names");

            if (!allBlank(seqs)) {
                sequenceData.push_back(seqs);
                sequenceDataLabels.push_back("sequences");
            }
            if (!allBlank(comments)) {
                sequenceData.push_back(comments);
                sequenceDataLabels.push_back("comments");
            }
            if (!allBlank(taxonomies)) {
                sequenceData.push_back(taxonomies);
                sequenceDataLabels.push_back("taxonomies");
            }
            if (!allBlank(trashCodes)) {
                sequenceData.push_back(trashCodes);
                sequenceDataLabels.push_back("trash_codes");
            }
            sequenceData.push_back(tableSeqs);
            sequenceDataLabels.push_back("include_sequence");

            sequenceData.attr("names") = sequenceDataLabels;

            results.push_back(sequenceData);
            resultsLabels.push_back("sequence_data");

            // only create sequence report if its not blank
            if (!allBlank(seqs)) {

                // sequence report data.frame
                // starts, ends, lengths, ambigs, polymers, numns
                Rcpp::DataFrame sequenceReport = Rcpp::DataFrame::create(
                    Rcpp::Named("sequence_ids") = getIndexes(names),
                    Rcpp::_["starts"] = starts,
                    Rcpp::_["ends"] = ends,
                    Rcpp::_["lengths"] = lengths,
                    Rcpp::_["ambigs"] = ambigs,
                    Rcpp::_["longest_homopolymers"] = polymers,
                    Rcpp::_["num_ns"] = numns);

                results.push_back(sequenceReport);
                resultsLabels.push_back("sequence_report");
            }

            // count_data(id, abundance, sample, treatment)
            results.push_back(count.getAbundanceTable(names,
                                                      getIndexes(names),
                                                      "sequence", false));
            resultsLabels.push_back("sequence_abundance_table");
        }
    }

    if (!hasTags || setContains(t, "bin_data")) {
        // sequence bin table
        for (int i = 0; i < binTables.size(); i++) {
            // create bin taxonomy if needed
            fillTaxReport(binTables[i].label);

            Rcpp::List binList = binTables[i].exportBinTable(count);
            vector<string> binListNames = binList.attr("names");

            for (int j = 0; j < binList.size(); j++) {
                results.push_back(binList[j]);
                resultsLabels.push_back(binTables[i].label+"_"+binListNames[j]);
            }
        }
    }

    if (!hasTags || setContains(t, "references")) {
        if (!references.empty()) {
            results.push_back(getReferences());
            resultsLabels.push_back("references");
        }
    }

    if (!hasTags || setContains(t, "metadata")) {
        if (metadata.hasReport) {
            results.push_back(getMetadata());
            resultsLabels.push_back("metadata");
        }
    }

    if (!hasTags || setContains(t, "reports")) {
        if (!reports.empty()) {
            for (auto it = reports.begin(); it != reports.end(); it++) {
                results.push_back(getReports(it->first));
                resultsLabels.push_back(it->first);
            }
        }
    }

    results.attr("names") = resultsLabels;
    return results;
}
/******************************************************************************/
double Dataset::addReferences(const vector<Reference>& refs) {
    references.insert(references.end(), refs.begin(), refs.end());
    return refs.size();
}
/******************************************************************************/
void Dataset::addReport(Rcpp::DataFrame& report, string type) {
    auto it = reports.find(type);

    if (it == reports.end()) {
        // this is a new report
        reports[type].addReport(report);
    }else{
        it->second.addReport(report);
    }
}
/******************************************************************************/
void Dataset::addMetadata(Rcpp::DataFrame& data){
    metadata.addReport(data);
}
/******************************************************************************/
double Dataset::addSequences(const vector<string>& n,
                             vector<string> s,
                             vector<string> c,
                             Reference reference) {

    // must provide the same number of names and seqs
    if (s.size() == 0) {
        s.resize(n.size(), "");
    }

    // add to seqIndex
    int numSeqs = names.size();
    vector<int> countNames;
    for (int i = 0; i < n.size(); i++) {
        countNames.push_back(numSeqs);
        seqIndex[n[i]] = numSeqs;
        numSeqs++;
    }

    // add to count
    count.add(countNames);
    countNames.clear();

    // add to names
    names.insert(names.end(), n.begin(), n.end());

    // add to seqs
    seqs.insert(seqs.end(), s.begin(), s.end());

    // innocent
    vector<bool> ts(n.size(), true);
    tableSeqs.insert(tableSeqs.end(), ts.begin(), ts.end());

    if (c.size() == 0) {
        c.resize(n.size(), "");
    }
    // add to comments
    comments.insert(comments.end(), c.begin(), c.end());

    // blank trash code because we assume seqs are "good"
    vector<string> t(n.size(), "");
    trashCodes.insert(trashCodes.end(), t.begin(), t.end());

   // add calcs for starts, ends, lengths, ambigs, polymers, numns
    SeqReport report;
    report.addReports(s, starts, ends, lengths, ambigs, polymers, numns);

    // set isAligned and aligned length
    getAlignedLength();

    // add to "good" sequence count - giving preference to binTable
    numUnique += n.size();

    hasSequenceData = true;

    if (reference.name != "") {
        references.push_back(reference);
    }

    return names.size();
}
/******************************************************************************/
double Dataset::assignBins(const vector<string>& binIds,
                            vector<float> abunds,
                            vector<string> samples,
                            vector<string> seqIds, const string type) {


    double numBinsAdded = 0;
    int binTableIndex = getBinTableIndex(type);

    // new table type
    if (binTableIndex == -1) {
        BinTable binTable;
        binTable.label = type;
        binTableIndex = binTables.size();
        binTables.push_back(binTable);
    }

    if (seqIds.empty()) {
        seqIds = binIds;

        // if you have already added sequences,
        // and are now adding shared or rabund data
        // make sure the binIds are in the dataset.
        if (hasSequenceData) {

            for (int i = 0; i < binIds.size(); i++) {

                auto it = seqIndex.find(binIds[i]);
                if (it == seqIndex.end()) {
                    string message = "The dataset does not contain a ";
                    message += "sequence named " + binIds[i] + ". If you have ";
                    message += "already added sequence data, you cannot add abundance";
                    message += " only bin assignments.\n";
                    RcppThread::Rcerr << endl << message << endl;
                    throw Rcpp::exception(message.c_str());
                }
            }
        }
    }else{
        hasList = true;
    }

    // if you don't have sequence data, but are adding seq_names for bins,
    // add seqs, add abundances if provided
    if (!hasSequenceData) {

        vector<string> uniqueSeqIds = unique(seqIds);
        addSequences(uniqueSeqIds);

        // if the abundances are provided then use them
        if (!abunds.empty()) {
            assignSequenceAbundance(seqIds, abunds, samples);
        }
    }
    abunds.clear();
    samples.clear();

    numBinsAdded = binTables[binTableIndex].assignBins(count, binIds,
                                                       getIndexes(seqIds));

    return numBinsAdded;
}
/******************************************************************************/
double Dataset::assignBinRepresentativeSequences(const vector<string>& binNames,
                                                 const vector<string>& repNames,
                                                 const string type){
    double repAssigned = 0;

    if (hasBinTable(type)) {

        if (!hasSequenceData) {
            addSequences(repNames);
        }

        repAssigned = binTables[getBinTableIndex(type)].assignRepresentativeSequences(binNames, getIndexes(repNames));
    }

    return repAssigned;
}
/******************************************************************************/
double Dataset::assignBinTaxonomy(const vector<string>& binIds,
                                  const vector<string>& taxs,
                                  const string type) {

    double numBinTaxonomiesAdded = 0;

    if (hasBinTable(type)) {
        numBinTaxonomiesAdded = binTables[getBinTableIndex(type)].assignTaxonomy(binIds, taxs);
    }

    return numBinTaxonomiesAdded;
}
/******************************************************************************/
double Dataset::assignTreatments(const vector<string>& samples,
                                 const vector<string>& treatments) {
    return count.assignTreatments(samples, treatments);
}
/******************************************************************************/
// names, abundances, samples(optional), treatments(optional)
// assumes same size
double Dataset::assignSequenceAbundance(vector<string>& ids,
                                        const vector<float>& abunds,
                                        const vector<string> samples,
                                        const vector<string> treatments) {

    double numSeqsAssigned = 0;

    vector<string> uniqueNames = unique(ids);

    // are there assignments for all seqs in the dataset
    if (uniqueNames.size() != numUnique){
        string message = "The dataset contains ";
        message += toString(numUnique) + " sequences, but you assigned ";
        message += toString(uniqueNames.size()) + " sequences. All sequences ";
        message += "in the dataset must be assigned abundances.\n\n";
        throw Rcpp::exception(message.c_str());
    }

    vector<int> idIndexes = getIndexes(ids);

    numSeqsAssigned = count.assignAbundance(idIndexes, abunds,
                                            samples, treatments);

    return numSeqsAssigned;
}
/******************************************************************************/
double Dataset::assignSequenceTaxonomy(const vector<string>& n,
                                       const vector<string>& t){

    double numSeqsTaxonomyAssigned = 0;

    if (!hasSequenceData) {
        addSequences(n);
    }

    // allocate space
    if (taxonomies.size() != names.size()) {
        taxonomies.resize(names.size(), "");
    }

    for (int i = 0; i < n.size(); i++) {

        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {

            int index = it->second;
            numSeqsTaxonomyAssigned++;

            // update taxonomy
            taxonomies[index] = t[i];
        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    hasSequenceTaxonomy = true;

    for (int i = 0; i < binTables.size(); i++) {
        binTables[i].runClassify = true;
    }
    return numSeqsTaxonomyAssigned;
}
/******************************************************************************/
/*
 id         level  taxon                 confidence
 seq1     1       Bacteria             100.0
 seq1     2      "Acidobacteria"  99.8
 seq1     3      Holophagae        99.8
 seq1     4      Holophagales      95.0
 seq1     5      Holophagaceae   90.0
 seq1     6      Holophaga           87.0
 seq2 ...
 */
Rcpp::DataFrame Dataset::fillTaxReport(string mode) {

    vector<string> ids, taxes;

    if (mode != "sequence") {
        if (hasBinTable(mode)) {
            int binTableIndex = getBinTableIndex(mode);

            // no sequence or bin taxonomy
            if ((allBlank(taxonomies)) &&
                (!binTables[binTableIndex].hasBinTaxonomy)) {
                return Rcpp::DataFrame::create();
            }else {
                ids = binTables[binTableIndex].getIds(count);
                taxes = binTables[binTableIndex].getTaxonomies(taxonomies, count);
            }
        }else{
            return Rcpp::DataFrame::create();
        }
    }else {
        ids = getSequenceNames();
        taxes = select(taxonomies, tableSeqs);
    }

    Utils util;
    int maxLevel = 1;

    vector<vector<string> > taxons(taxes.size());
    vector<vector<int> > confidences(taxes.size());
    bool hasConfidences = true;

    // split taxonomy and confidences by level
    for (int i = 0; i < taxes.size(); i++) {
        int numLevels = split(taxes[i], ';',
                              back_inserter(taxons[i]));

        if (numLevels > maxLevel) { maxLevel = numLevels; }

        confidences[i] = util.removeConfidences(taxons[i]);

        if (hasConfidences) {
            if (sum(confidences[i]) == 0) {
                hasConfidences = false;
            }
        }
    }

    taxes.clear();
    vector<string> dfIds(taxons.size()*maxLevel);
    vector<string> dfTaxs(taxons.size()*maxLevel);
    vector<int> levels(taxons.size()*maxLevel);
    vector<int> dfConfidences(taxons.size()*maxLevel);

    for (int i = 0; i < taxons.size(); i++) {
        // extend taxonomies to the same level by adding unclassifieds
        util.addUnclassifieds(taxons[i], confidences[i], maxLevel);

        for (int j = 0; j < maxLevel; j++) {
            dfIds[i*maxLevel+j] = ids[i];
            dfTaxs[i*maxLevel+j] = taxons[i][j];
            levels[i*maxLevel+j] = j+1;
            dfConfidences[i*maxLevel+j] = confidences[i][j];
        }
    }

    taxons.clear();
    confidences.clear();

    if (!hasConfidences) {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("id") = dfIds,
            Rcpp::_["level"] = levels,
            Rcpp::_["taxon"] = dfTaxs);
        return df;
    }

    Rcpp::DataFrame df = Rcpp::DataFrame::create(
        Rcpp::Named("id") = dfIds,
        Rcpp::_["level"] = levels,
        Rcpp::_["taxon"] = dfTaxs,
        Rcpp::_["confidence"] = dfConfidences);

    return df;
}
/******************************************************************************/
int Dataset::getAlignedLength() {
    set<int> seqLengths;
    for (int i = 0; i < seqs.size(); i++) {
        // is this a "good" seq
        if (tableSeqs[i]) {
            seqLengths.insert(seqs[i].length());
        }
    }

    if (seqLengths.size() == 1) {
        isAligned = true;
        alignmentLength = *(seqLengths.begin());

        if (alignmentLength == 0) {
            isAligned = false;
            alignmentLength = -1;
        }
    }else{
        isAligned = false;
        alignmentLength = -1;
    }
    return alignmentLength;
}
/******************************************************************************/
// fasta data.frame 2 or 3 columns, sequence_names, sequences, comments
const Rcpp::DataFrame Dataset::getFastaReport() {
    vector<string> n(numUnique, "");
    vector<string> s(numUnique, "");
    vector<string> c(numUnique, "");

    int index = 0;
    for (int i = 0; i < tableSeqs.size(); i++) {
        if (tableSeqs[i]) {
            n[index] = names[i];
            s[index] = seqs[i];
            c[index] = comments[i];
            index++;
        }
    }

    if (allBlank(c)) {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("sequence_names") = n,
            Rcpp::_["sequences"] = s);

        return df;
    }else{
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("sequence_names") = n,
            Rcpp::_["sequences"] = s,
            Rcpp::_["comments"] = c);

        return df;
    }
    return Rcpp::DataFrame::create();
}
/******************************************************************************/
// returns indexes of "good" seqs in table
const vector<int> Dataset::getIncludedNamesIndexes() {
    vector<int> included;
    for (int i = 0; i < tableSeqs.size(); i++) {
        if (tableSeqs[i]) {
            included.push_back(i);
        }
    }
    return included;
}
/******************************************************************************/
const vector<int> Dataset::getIndexes(const vector<string>& ids) {
    vector<int> indexes(ids.size(), -1);

    for (int i = 0; i < ids.size(); i++) {

        auto it = seqIndex.find(ids[i]);
        if (it != seqIndex.end()) {
            indexes[i] = it->second;
        }else{
            string message = "The dataset does not contain a ";
            message += "sequence named " + ids[i] + ".\n";
            RcppThread::Rcerr << endl << message << endl;
            throw Rcpp::exception(message.c_str());
        }
    }
    return indexes;
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getList(string type) {
    if (hasBinTable(type) && hasList) {
        return binTables[getBinTableIndex(type)].getList(names);
    }else if ((type == "asv") && (hasSequenceData)) {
        vector<string> seqIds = getSequenceNames();
        vector<string> asvIds(seqIds.size(), "");
        for (int i = 0; i < asvIds.size(); i++) {
            asvIds[i] = "ASV" + toString(i+1);
        }
        assignBins(asvIds, nullFloatVector, nullVector, seqIds, "asv");
        return binTables[getBinTableIndex(type)].getList(names);
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getMetadata() {
    return metadata.getReport(nullSet);
}
/******************************************************************************/
const vector<string> Dataset::getListVector(string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getListVector(names);
    }

    return nullVector;
}
/******************************************************************************/
const vector<string> Dataset::getSequenceNames(vector<string> samples,
                                               bool distinct){
    vector<string> included;

    // get all "good" names in dataset
    if (samples.empty())  {
        for (int i = 0; i < tableSeqs.size(); i++) {
            if (tableSeqs[i]) {
                included.push_back(names[i]);
            }
        }
    // get all "good" names in specific set of samples
    }else {

        // all seqs
        for (int i = 0; i < tableSeqs.size(); i++) {
            // if "good" seq
            if (tableSeqs[i]) {

                if (!distinct) {
                    // include all the requested samples, but may have
                    // additional samples present
                    if (count.hasSamples(samples, i)) {
                        included.push_back(names[i]);
                    }
                }else {
                    // sequences must have ONLY the samples requested
                    if (identical(count.getSamples(i), samples)) {
                        included.push_back(names[i]);
                    }
                }
            }
        }
    }

    return included;
}
/******************************************************************************/
const vector<vector<string> > Dataset::getSequenceNamesBySample(vector<string> samples){
    vector<vector<string> > result;

    // return all samples if none specified
    if (samples.size() == 0) {
        samples = getSamples();
    }

    for (int i = 0; i < samples.size(); i++) {
        vector<string> s;
        s.push_back(samples[i]);
        result.push_back(getSequenceNames(s));
    }

    return result;
}
/******************************************************************************/
const int Dataset::getNumSamples() {
    return count.getNumSamples();
}
/******************************************************************************/
const int Dataset::getNumTreatments() {
    return count.getNumTreatments();
}
/******************************************************************************/
const int Dataset::getNumResourceReferences() {
    return references.size();
}
/******************************************************************************/
// TODO document in module_exports.R
const int Dataset::getNumBins(string type, vector<string> samples,
                              bool distinct) {
    int numBins = 0;

    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getNumBins(count, samples, distinct);
    }

    return numBins;
}
/******************************************************************************/
const vector<string> Dataset::getBinIds(string type,
                                        vector<string> samples, bool distinct) {

    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getIds(count, samples, distinct);
    }
    return nullVector;
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getBinRepresentativeSequences(string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getRepresentativeSequences(names, seqs);
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getBinTaxonomyReport(string type) {
    return (fillTaxReport(type));
}
/******************************************************************************/
const vector<string> Dataset::getBinTypes() {
    vector<string> types;

    for (size_t i = 0; i < binTables.size(); i++) {
        types.push_back(binTables[i].label);
    }
    return types;
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getBinAbundances(string bin_type, bool bySample) {

    if (hasBinTable(bin_type)) {
        if (!bySample) {
            return binTables[getBinTableIndex(bin_type)].getRAbund(count);
        }else {
            return binTables[getBinTableIndex(bin_type)].getShared(count);
        }
    }

    return Rcpp::DataFrame::create();
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getReports(string type) {
    Rcpp::DataFrame reportResults = Rcpp::DataFrame::create();

    if (!reports.empty()) {
        auto it = reports.find(type);

        if (it != reports.end()) {
            return it->second.getReport(toSet(getSequenceNames()));
        }else{
            string message = "Your dataset does not include a report named ";
            message += type + ", ignoring request.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    return reportResults;
}
/******************************************************************************/
const vector<string> Dataset::getReportTypes() {

    //custom report types
    vector<string> reportTypes = getKeys(reports);

    if (hasSeqs()) {
        //reportTypes.push_back("sequence_data");

        if (badAccnos.size() != 0) {
            reportTypes.push_back("sequence_scrap");
        }
    }

    for (int i = 0; i < binTables.size(); i++) {
        if (binTables[i].getScrapReport().size() != 0) {
            reportTypes.push_back("bin_scrap");
            break;
        }
    }

    return reportTypes;
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getReferences() {

    if (!references.empty()) {

        vector<string> refNames(references.size(), "NA");
        vector<string> refNotes(references.size(), "NA");
        vector<string> refUrls(references.size(), "NA");
        vector<string> refUsages(references.size(), "NA");
        vector<string> refVersions(references.size(), "NA");

        for (size_t i = 0; i < references.size(); i++) {
            if (references[i].name != "") {
                refNames[i] = references[i].name;
            }
            if (references[i].version != "") {
                refVersions[i] = references[i].version;
            }
            if (references[i].note != "") {
                refNotes[i] = references[i].note;
            }
            if (references[i].url != "") {
                refUrls[i] = references[i].url;
            }
            if (references[i].usage != "") {
                refUsages[i] = references[i].usage;
            }
        }

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("reference_names") = refNames,
            Rcpp::_["reference_versions"] = refVersions,
            Rcpp::_["reference_usages"] = refUsages,
            Rcpp::_["reference_notes"] = refNotes,
            Rcpp::_["reference_urls"] = refUrls);

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
const vector<string> Dataset::getSamples(){
    return count.getSamples();
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getSampleTreatmentAssignments() {
    map<string, string> sampleToTreatment = count.getSampleTreatmentAssignments();

    if (sampleToTreatment.size() != 0) {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("samples") = getKeys(sampleToTreatment),
            Rcpp::_["treatments"] = getValues(sampleToTreatment));
        return df;
    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// id, trashCode
const Rcpp::DataFrame Dataset::getScrapReport(string mode) {

    if (mode == "sequence") {
        if (badAccnos.size() != 0) {
            vector<string> badNames(uniqueBad, "");
            vector<string> badCodes(uniqueBad, "");

            int next = 0;
            for (int i = 0; i < trashCodes.size(); i++) {
                if (trashCodes[i] != "") {
                    badNames[next] = names[i];
                    // remove last comma
                    trashCodes[i].pop_back();
                    badCodes[next] = trashCodes[i];
                    next++;
                }
            }

            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named("id") = badNames,
                Rcpp::_["trash_code"] = badCodes);
            return df;
        }
    }else {
        if (hasBinTable(mode)) {
            return binTables[getBinTableIndex(mode)].getScrapReport();
        }
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// type, trashCode, uniqueCount, totalCount
const Rcpp::DataFrame Dataset::getScrapSummary() {

    vector<string> types, codes;
    vector<float> uniqueCounts, totalCounts;

    if (badAccnos.size() != 0) {
        types.resize(badAccnos.size());
        codes.resize(badAccnos.size());
        uniqueCounts.resize(badAccnos.size());
        totalCounts.resize(badAccnos.size());

        int index = 0;
        for (auto it = badAccnos.begin(); it != badAccnos.end(); it++) {
            codes[index] = it->first;
            uniqueCounts[index] = it->second[0];
            totalCounts[index] = it->second[1];
            types[index] = "sequence";
            index++;
        }
    }

    if (binTables.size() != 0) {
        for (int i = 0; i < binTables.size(); i++) {
            Rcpp::DataFrame df = binTables[i].getScrapSummary();
            if (df.size() == 4) {
                vector<string> b = Rcpp::as<vector<string>>(df[0]);
                types.insert(types.end(), b.begin(), b.end());

                b = Rcpp::as<vector<string>>(df[1]);
                codes.insert(codes.end(), b.begin(), b.end());

                vector<float> f = Rcpp::as<vector<float>>(df[2]);
                uniqueCounts.insert(uniqueCounts.end(), f.begin(), f.end());

                f = Rcpp::as<vector<float>>(df[3]);
                totalCounts.insert(totalCounts.end(), f.begin(), f.end());
            }
        }
    }

    if (!types.empty()) {
        return Rcpp::DataFrame::create(
            Rcpp::Named("type") = types,
            Rcpp::_["trash_code"] = codes,
            Rcpp::_["unique"] = uniqueCounts,
            Rcpp::_["total"] = totalCounts);
    }

    return Rcpp::DataFrame::create();
}
/******************************************************************************/
// total abundance for each sequence
const Rcpp::DataFrame Dataset::getSequenceAbundances(bool bySample){

    if (!bySample) {
        vector<float> abunds = count.getTotalAbundances(getIncludedNamesIndexes());

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("sequence_names") = select(names, tableSeqs),
            Rcpp::_["abundances"] = abunds);

        return df;
    }

    // ids, abundances, sample(optional), treatment(optional)
    return count.getAbundanceTable(select(names, tableSeqs),
                                    getIncludedNamesIndexes());
}
/******************************************************************************/
const vector<vector<float> > Dataset::getSequenceAbundanceBySample(vector<string> samples) {
    vector<int> ids = getIncludedNamesIndexes();
    return count.getAbundanceBySample(ids, samples);
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getSequenceTable(string sample) {

    if (hasSequenceData) {

        vector<string> samples;

        if(sample != "") {
            samples.push_back(sample);
        }

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("sequence_names") = getSequenceNames(samples),
            Rcpp::_["sequences"] = getSequences(sample)
        );

        return df;
    }

    return Rcpp::DataFrame::create();
}
/******************************************************************************/
const vector<string> Dataset::getSequences(string sample){
    vector<string> included;

    if (sample == "") {
        for (int i = 0; i < tableSeqs.size(); i++) {
            if (tableSeqs[i]) {
                included.push_back(seqs[i]);
            }
        }
    }else{
        if (count.hasSample(sample)) {
            // all seqs
            for (int i = 0; i < tableSeqs.size(); i++) {
                // if "good" seq
                if (tableSeqs[i]) {
                    // if this sequence is in sample
                    if (count.hasSample(sample, i)) {
                        included.push_back(seqs[i]);
                    }
                }
            }
        }
    }

    return included;
}
/******************************************************************************/
const vector<vector<string> > Dataset::getSequencesBySample(vector<string> samples){
    vector<vector<string> > result;

    // return all samples if none specified
    if (samples.size() == 0) {
        samples = getSamples();
    }

    for (int i = 0; i < samples.size(); i++) {
        result.push_back(getSequences(samples[i]));
    }

    return result;
}
/******************************************************************************/
// fasta summary data: starts, ends, lengths, ambigs, polymers, numns
const Rcpp::DataFrame Dataset::getSequenceReport(){

    Rcpp::DataFrame df = Rcpp::DataFrame::create(
        Rcpp::Named("id") = select(names, tableSeqs),
        Rcpp::_["start"] = select(starts, tableSeqs),
        Rcpp::_["end"] = select(ends, tableSeqs),
        Rcpp::_["length"] = select(lengths, tableSeqs),
        Rcpp::_["ambig"] = select(ambigs, tableSeqs),
        Rcpp::_["longest_homopolymer"] = select(polymers, tableSeqs),
        Rcpp::_["num_n"] = select(numns, tableSeqs));

    return df;
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getSummary(string type, string reportType) {

    Rcpp::DataFrame result = Rcpp::DataFrame::create();

    if (type == "sequences") {
        if (hasSeqs()) {
            Summary* summary = new Summary(processors);

            vector<vector<int> > report;

            report.push_back(select(starts, tableSeqs));
            report.push_back(select(ends, tableSeqs));
            report.push_back(select(lengths, tableSeqs));
            report.push_back(select(ambigs, tableSeqs));
            report.push_back(select(polymers, tableSeqs));
            report.push_back(select(numns, tableSeqs));

            result = summary->summarizeFasta(
                report, count.getTotalAbundances(getIncludedNamesIndexes()));

            delete summary;
        }
    }else if (type == "reports") {
        auto it = reports.find(reportType);

        // do we have this report type
        if (it != reports.end()) {
            Rcpp::DataFrame df = getSequenceAbundances();
            result = it->second.summarizeReport(
                    toSet(getSequenceNames()), processors, Rcpp::as<vector<float>>(df[1]));
        }
    }else if (type == "scrap") {
        result = getScrapSummary();
    }

    return result;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getSequenceTaxonomyReport() {
    if (hasSequenceTaxonomy){
        return (fillTaxReport("sequence"));
    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
const vector<string> Dataset::getTreatments(){
    return count.getTreatments();
}
/******************************************************************************/
const double Dataset::getTotal(vector<string> samples){

    if (samples.empty()) {
        return count.getTotal();
    }else {
       if (hasSequenceData) {

           double seqTotal = 0;

           // all seqs
           for (int i = 0; i < tableSeqs.size(); i++) {

               // if "good" seq
               if (tableSeqs[i]) {

                   // include all the requested samples, but may have
                   // additional samples present
                   if (count.hasSamples(samples, i)) {
                       seqTotal += count.getAbundance(i, samples);
                   }
               }
           }

           return seqTotal;
       }
    }

    return 0;
}
/******************************************************************************/
const Rcpp::DataFrame Dataset::getTotals(string type){

    vector<double> totals;
    vector<string> names;

    if (type == "samples") {
        totals = count.getSampleTotals();
        names = getSamples();
    }else if (type == "treatments") {
        totals = count.getTreatmentTotals();
        names = getTreatments();
    }

    if (!totals.empty()) {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(type) = names,
            Rcpp::_["abundances"] = totals);
        return df;
    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
const double Dataset::getUniqueTotal(vector<string> samples){
    if (samples.empty()) {
        return numUnique;
    }
    return getSequenceNames(samples, true).size();
}
/******************************************************************************/
const bool Dataset::hasBinTable(string type) {
    for (int i = 0; i < binTables.size(); i++) {
        if (binTables[i].label == type) {
            return true;
        }
    }
    return false;
}
/******************************************************************************/
const int Dataset::getBinTableIndex(string type) {
    for (int i = 0; i < binTables.size(); i++) {
        if (binTables[i].label == type) {
            return i;
        }
    }
    return -1;
}
/******************************************************************************/
const bool Dataset::hasSample(string sample){
    return count.hasSample(sample);
}
/******************************************************************************/
const bool Dataset::hasSamples(vector<string> samples) {

    int numFound = 0;
    for (string sample : samples) {
        if (hasSample(sample)) {
            numFound++;
        }
    }

    if (numFound == samples.size()) {
        return true;
    }
    return false;
}
/******************************************************************************/
const bool Dataset::hasSeqs() {
    if (seqs.size() == 0) {
        return false;
    }

    string id;
    if (allIdentical(seqs, id)) {
        if (id == "") { return false; }
    }
    return true;
}
/******************************************************************************/
void Dataset::mergeSequences(const vector<string>& ids, string reason){
    if (ids.size() != 1) {

        vector<int> indexes = getIndexes(ids);

        // sanity check: if you have assigned bins, make sure the sequences
        // are in the same bin
        if (binTables.size() != 0) {
            // innocent until proven guilty
            vector<bool> okToMerge(binTables.size(), true);

            int tableIndex = 0;
            string message = "Can not merge sequences assigned";
            message += " to different bins.";
            for (int i = 0; i < binTables.size(); i++) {

                // if seqs were assigned to otus
                if (!binTables[i].okToMerge(indexes)) {
                    okToMerge[tableIndex] = false;
                    message += " The '" + binTables[i].label + "' bin table has";
                    message += " sequences assigned to different bins.";
                }
                tableIndex++;
            }

            if (!isTrue(okToMerge)) {
                RcppThread::Rcerr << endl << message << endl;
                throw Rcpp::exception(message.c_str());
            }
        }

        count.merge(indexes);

        for (int i = 1; i < indexes.size(); i++) {
            // no need to update the sample and treatment counts
            removeSequence(indexes[i], reason, false, false);
        }
    }
}
/******************************************************************************/
void Dataset::mergeBins(const vector<string>& ids, string reason, string type){
    if (ids.size() != 1) {
        if (hasBinTable(type)) {
            binTables[getBinTableIndex(type)].merge(ids, reason);
        }
    }
}
/******************************************************************************/
void Dataset::removeLineages(const vector<string>& contaminants, string trashTag) {
    // if no sequence taxonomy or clusters
    if (!hasSequenceTaxonomy && (binTables.size() == 0)) { return; }

    Utils util;
    vector<vector<string> > conTax(contaminants.size());
    vector<vector<int> > conConfidenceThreshold(contaminants.size());
    vector<bool> conHasConfidences(contaminants.size(), false);

    // split contaminant taxonomies and confidences by level
    for (int i = 0; i < contaminants.size(); i++) {

        split(contaminants[i], ';', back_inserter(conTax[i]));
        conConfidenceThreshold[i] = util.removeConfidences(conTax[i]);

        if (sum(conConfidenceThreshold[i]) != 0) {
            conHasConfidences[i] = true;
        }
    }

    if (hasSequenceTaxonomy) {
        // you have sequence taxonomies, remove contaminants
        // and reclassify bins (below)

        // remove seqs assigned to this taxonomy
        for (int i = 0; i < taxonomies.size(); i++) {

            if (tableSeqs[i]) {
                vector<string> userTax;
                split(taxonomies[i], ';', back_inserter(userTax));
                vector<int> userConfidences = util.removeConfidences(userTax);

                // if this seq is a contaminant, remove it
                if (util.searchTax(userTax, userConfidences,
                                   conHasConfidences,
                                   conTax, conConfidenceThreshold)) {
                    removeSequence(i, trashTag, true, true);
                }
            }
        }
    }

    for (int i = 0; i < binTables.size(); i++) {

        vector<string> binIds = binTables[i].getIds(count);
        vector<string> binTaxes;

        // if you have only have bin classifications, and no seq classifications
        if (!hasSequenceTaxonomy) {
            if (binTables[i].hasBinTaxonomy) {
                binTaxes = binTables[i].getTaxonomies(binTaxes, count);
            }
        }else{
            binTaxes = binTables[i].getTaxonomies(taxonomies, count);
        }

        // remove bins assigned to this taxonomy
        for (int j = 0; j < binTaxes.size(); j++) {

            vector<string> userTax;
            split(binTaxes[j], ';', back_inserter(userTax));
            vector<int> userConfidences = util.removeConfidences(userTax);

            // if this seq is a contaminant, remove it
            if (util.searchTax(userTax, userConfidences,
                               conHasConfidences,
                               conTax, conConfidenceThreshold)) {
                vector<string> binToRemove(1, binIds[j]);
                vector<string> trashTags(1, trashTag);

                removeBins(binToRemove, trashTags, binTables[i].label);
            }
        }
    }

    count.updateTotals();
}
/******************************************************************************/
void Dataset::removeBins(const vector<string>& namesToRemove,
                         const vector<string>& trashTags, string type){

    if (binTables.size() == 0) { return; }

    if (hasBinTable(type)) {

        if (namesToRemove.size() != trashTags.size()) {
            string message = "Size mismatch. You must provide a trash";
            message += " code for each bin.";
            RcppThread::Rcerr << endl << message << endl;
            throw Rcpp::exception(message.c_str());
        }

        for (int i = 0; i < namesToRemove.size(); i++) {
            vector<int> seqsToRemove = binTables[getBinTableIndex(type)].remove(namesToRemove[i], trashTags[i]);

            // remove any sequences from removed bin
            for (int seq : seqsToRemove) {
                removeSequence(seq, trashTags[i], true, false);
            }

            // remove from other lists
            for (int j = 0; j < binTables.size(); j++) {
                if (binTables[j].label != type) {
                    for (int seq : seqsToRemove) {
                        binTables[j].removeSeq(seq, trashTags[i]);
                    }
                }
            }
        }

        count.updateTotals();
    }
}
/******************************************************************************/
void Dataset::removeSamples(const vector<string>& samples, string reason) {
    AbundTable dupCount(count);

    dupCount.removeSamples(samples);

    // // remove samples from count
    // remove any seqs only assigned to these samples
    for (int i = 0; i < names.size(); i++) {
        // included seq
        if (tableSeqs[i]) {
            if (sum(dupCount.getAbundances(i)) == 0) {
                removeSequence(i, reason, true, false);

                // remove from list otus
                for (int j = 0; j < binTables.size(); j++) {
                    binTables[j].removeSeq(i, reason);
                }
            }
        }
    }

    count.removeSamples(samples);
}
/******************************************************************************/
void Dataset::removeSequence(const int index, const string reasons,
                             bool update, bool removeFromBin) {
    // remove from tableSeqs and add trashCode
    tableSeqs[index] = false;
    trashCodes[index] += reasons + ",";
    numUnique--;

    // remove from counts
    int abund = 1;
    if (update) {
        if ((binTables.size() != 0) && removeFromBin) {

            // remove from list otus
            for (int i = 0; i < binTables.size(); i++) {
               binTables[i].removeSeq(index, reasons);
            }
        }

        abund = count.remove(index);
    }else{
        abund = count.getAbundance(index);
    }

    // add to badAccnos
    vector<string> theseReasons;
    split(reasons, ',', back_inserter(theseReasons));

    for (int j = 0; j < theseReasons.size(); j++) {
        auto itBad = badAccnos.find(theseReasons[j]);

        if (itBad != badAccnos.end()) {
            // update counts of trashCode
            itBad->second[0]++;
            itBad->second[1] += abund;
        }else{
            // add new trashCode
            vector<double> badAbunds(2, 1);
            badAbunds[1] = abund;
            badAccnos[theseReasons[j]] = badAbunds;
        }
    }

    // update uniqueBad
    uniqueBad++;
}
/******************************************************************************/
void Dataset::removeSequences(const vector<string>& namesToRemove,
                              const vector<string>& trashTags){

    if (namesToRemove.size() != trashTags.size()) {
        string message = "Size mismatch. You must provide a trash";
        message += " code for each sequence.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    for (int i = 0; i < namesToRemove.size(); i++) {
        auto it = seqIndex.find(namesToRemove[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            removeSequence(index, trashTags[i], true, true);
        }else{
            string message = "[WARNING]: " + namesToRemove[i] + " is not in ";
            message += "your dataset, ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count.updateTotals();
}
/******************************************************************************/
// for datasets without samples
void Dataset::setAbundance(const vector<string>& n, const vector<float>& abunds,
                            string reason){

    if (n.size() != abunds.size()) {
        string message = "Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    float diff = 0;

    for (int i = 0; i < n.size(); i++) {
        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            if (abunds[i] == 0) {
                removeSequence(index, reason, true, true);
            }else{
                diff += count.setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count.updateTotals();

    if (diff > 0) {
        // may not be removing the whole sequence perhaps just reducing
        auto itBad = badAccnos.find(reason);

        if (itBad != badAccnos.end()) {
            // update counts of trashCode
            itBad->second[1] += diff;
        }else{
            // add new trashCode
            vector<double> badAbunds(2, 0);
            badAbunds[1] = diff;
            badAccnos[reason] = badAbunds;
        }
    }
}
/******************************************************************************/
// for datasets with samples
void Dataset::setAbundances(const vector<string>& n,
                            const vector<vector<float>>& abunds,
                            string reason){

    if (n.size() != abunds.size()) {
        string message = "Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    float diff = 0;

    for (int i = 0; i < n.size(); i++) {
        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            if (sum(abunds[i]) == 0) {
                removeSequence(index, reason, true, true);
            }else{
                diff += count.setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count.updateTotals();

    if (diff > 0) {
        // may not be removing the whole sequence perhaps just reducing
        auto itBad = badAccnos.find(reason);

        if (itBad != badAccnos.end()) {
            // update counts of trashCode
            itBad->second[1] += diff;
        }else{
            // add new trashCode
            vector<double> badAbunds(2, 0);
            badAbunds[1] = diff;
            badAccnos[reason] = badAbunds;
        }
    }

}
/******************************************************************************/
void Dataset::setSequences(const vector<string>& n, const vector<string>& s,
                           const vector<string> c){
    if (n.size() != s.size()) {
        string message = "Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    bool hasComments = false;
    if (c.size() != 0) {
        if (c.size() != n.size()) {
            string message = "Size mismatch. When providing comments,";
            message += " ids and comments must be the same size.";
            throw Rcpp::exception(message.c_str());
        }
        hasComments = true;
    }

    SeqReport report;

    for (int i = 0; i < n.size(); i++) {

        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {

            int index = it->second;

            // update sequence
            seqs[index] = s[i];

            // update start, end, numbase, ambig, polymer, numn
            vector<int> reportResults = report.getReport(s[i]);

            starts[index] = reportResults[0];
            ends[index] = reportResults[1];
            lengths[index] = reportResults[2];
            ambigs[index] = reportResults[3];
            polymers[index] = reportResults[4];
            numns[index] = reportResults[5];

            // update comment
            if (hasComments) {
                comments[index] += " " + c[i];
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    // set isAligned and aligned length
    getAlignedLength();
}
/******************************************************************************/
const Rcpp::RawVector Dataset::serializeDataset() {
    std::stringstream ss; // Use stringstream for in-memory serialization
    {
        cereal::BinaryOutputArchive oarchive(ss);
        oarchive(*this);
    }

    // Convert stringstream content to Rcpp::RawVector
    std::string serialized_data = ss.str();
    Rcpp::RawVector retval(serialized_data.size());
    std::copy(serialized_data.begin(), serialized_data.end(), retval.begin());
    return retval;
}
/******************************************************************************/


