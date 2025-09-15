
#include "../inst/include/rdataset.h"
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

    count.clone(dataset.count);

    for (int i = 0; i < dataset.binTables.size(); i++) {
        BinTable table;
        table.clone(dataset.binTables[i]);
        binTables.push_back(table);
    }
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
void Dataset::clear(vector<string> tags) {

    // clear all
    if (tags.empty()) {
        isAligned = false;
        hasSequenceData = false;
        hasSequenceTaxonomy = false;
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
    }else {

        for (string tag : tags) {
            if (tag == "sequence_taxonomy") {
                hasSequenceTaxonomy = false;
                taxonomies.clear();
            }else if (tag == "bin_taxonomy") {
                for (int i = 0; i < binTables.size(); i++) {
                    binTables[i].clear("taxonomy");
                }
            }else if (tag == "bin_assignment") {
                binTables.clear();
            }
        }
    }
}
/******************************************************************************/
Rcpp::List Dataset::exportDataset(){
    Rcpp::List result;

    // sequence data.frame
    // name, seqs, comments

    // sequence report data.frame
    // starts, ends, lengths, ambigs, polymers, numns

    // sequence trashCodes
    // name, trashCode

    // sequence abundance table
    // id, abundance, sample, treatment

    // sequence bin table
    // id, abundance, sample, seq_id

    // sequence taxonomy table

    // bin taxonomy table

    // bin trashCodes

    // align report

    // contigs report

    // metadata




    return result;
}
/******************************************************************************/
void Dataset::addSequences(vector<string> n, vector<string> s, vector<string> c) {

    // must provide the same number of names and seqs
    if (s.size() == 0) {
        s.resize(names.size(), "");
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

    if (s.size() == 0) {
        s.resize(names.size(), "");
    }

    // add to seqs
    seqs.insert(seqs.end(), s.begin(), s.end());

    // innocent
    vector<bool> ts(n.size(), true);
    tableSeqs.insert(tableSeqs.end(), ts.begin(), ts.end());

    if (c.size() == 0) {
        c.resize(names.size(), "");
    }
    // add to comments
    comments.insert(comments.end(), c.begin(), c.end());

    // blank trash code because we assume seqs are "good"
    vector<string> t(names.size(), "");
    trashCodes.insert(trashCodes.end(), t.begin(), t.end());

   // add calcs for starts, ends, lengths, ambigs, polymers, numns
    SeqReport report;
    report.addReports(s, starts, ends, lengths, ambigs, polymers, numns);

    // set isAligned and aligned length
    getAlignedLength();

    // add to "good" sequence count - giving preference to binTable
    numUnique += names.size();

    hasSequenceData = true;
}
/******************************************************************************/
void Dataset::assignBins(vector<string> binIds,
                        vector<int> abunds,
                        vector<string> samples,
                        vector<string> seqIds, string type) {



    bool useSeqIds = false;
    if (!seqIds.empty()) { useSeqIds = true;  }

    if (!useSeqIds && (binIds.size() != abunds.size())) {
        string message = "[ERROR]: Size mismatch. bin_ids and abunds must be";
        message += " the same size.";
        throw Rcpp::exception(message.c_str());
    }

    int binTableIndex = getBinTableIndex(type);

    // new table type
    if (binTableIndex == -1) {
        BinTable binTable;
        binTable.label = type;
        binTableIndex = binTables.size();
        binTables.push_back(binTable);
    }

    if (useSeqIds) {

        // if you don't have sequence data, but are adding seq_names for bins,
        // add seqs, add abundances if provided
        if (!hasSequenceData) {

            vector<string> uniqueSeqIds = unique(seqIds);
            addSequences(uniqueSeqIds);

            // if the abundances are provided then use them
            if (!abunds.empty()) {
                assignSequenceAbundance(seqIds, abunds, samples);
            }
        }else {
            if (!abunds.empty() || !samples.empty()) {
                string message = "[WARNING]: Ignoring sample abundances, using";
                message += " sequence sample abundances already assigned.";
                RcppThread::Rcout << endl << message << endl;
            }
        }
        abunds.clear();
        samples.clear();

        binTables[binTableIndex].assignAbundance(binIds, abunds, samples,
                                  getIndexes(seqIds), count);
    }else{
        // add bins - shared
        binTables[binTableIndex].assignAbundance(binIds, abunds, samples,
                                  nullIntVector, count);
    }
}
/******************************************************************************/
void Dataset::assignBinTaxonomy(vector<string> binIds, vector<string> taxs,
                                string type) {

    if (hasBinTable(type)) {
        if (binIds.size() != taxs.size()) {
            string message = "[ERROR]: Size mismatch. bin_ids and taxonomies ";
            message += "must be the same size.";
            throw Rcpp::exception(message.c_str());
        }

        binTables[getBinTableIndex(type)].assignTaxonomy(binIds, taxs);
    }
}
/******************************************************************************/
void Dataset::assignTreatments(vector<string> samples, vector<string> treatments) {
    count.assignTreatments(samples, treatments);

    // for each type of binTable, pass new sample assignments
    for (int i = 0; i < binTables.size(); i++) {
        binTables[i].assignTreatments(samples, treatments);
    }
}
/******************************************************************************/
// names, abundances, samples(optional), treatments(optional)
// assumes same size
void Dataset::assignSequenceAbundance(vector<string> ids,
                           vector<int> abunds,
                           vector<string> samples,
                           vector<string> treatments) {

    vector<string> uniqueNames = unique(ids);

    // are there assignments for all seqs in the dataset
    if (uniqueNames.size() != numUnique){
        string message = "[ERROR]: The dataset contains ";
        message += toString(numUnique) + " sequences, but you assigned ";
        message += toString(uniqueNames.size()) + " sequences. All sequences ";
        message += "in the dataset must be assigned abundances.\n\n";
        throw Rcpp::exception(message.c_str());
    }

    vector<int> idIndexes = getIndexes(ids);

    count.assignAbundance(idIndexes, abunds, samples, treatments);

    // update list assignments if they include sequence ids
    for (int i = 0; i < binTables.size(); i++) {
        if (binTables[i].hasListAssignments) {
            binTables[i].assignAbundance(nullVector, nullIntVector,
                                         nullVector, nullIntVector,
                                         count, true);
        }
    }
}
/******************************************************************************/
void Dataset::assignSequenceTaxonomy(vector<string> n, vector<string> t){

    if (n.size() != t.size()) {
        string message = "[ERROR]: Size mismatch. ids and taxonomies must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

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

            // update taxonomy
            taxonomies[index] = t[i];
        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    hasSequenceTaxonomy = true;
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
            ids = binTables[binTableIndex].getIds();
            taxes = binTables[binTableIndex].getTaxonomies(taxonomies, count);

            // taxes is empty if no bin classifications have been added
            if (taxes.empty()) {
                Rcpp::DataFrame empty = Rcpp::DataFrame::create();
                return empty;
            }
        }else{
            Rcpp::DataFrame empty = Rcpp::DataFrame::create();
            return empty;
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
int Dataset::getAbundance(string name){
    int abund = 0;
    auto it = seqIndex.find(name);
    if (it != seqIndex.end()) {
        if (tableSeqs[it->second]) {
            abund = count.getAbundance(it->second, "");
        }
    }
    return abund;
}
/******************************************************************************/
vector<int> Dataset::getAbundances(string name){
    vector<int> abunds;
    auto it = seqIndex.find(name);
    if (it != seqIndex.end()) {
        if (tableSeqs[it->second]) {
            abunds = count.getAbundances(it->second);
        }
    }
    return abunds;
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
// returns indexes of "good" seqs in table
vector<int> Dataset::getIncludedNamesIndexes() {
    vector<int> included;
    for (int i = 0; i < tableSeqs.size(); i++) {
        if (tableSeqs[i]) {
            included.push_back(i);
        }
    }
    return included;
}
/******************************************************************************/
vector<int> Dataset::getIndexes(vector<string>& ids) {
    vector<int> indexes(ids.size(), -1);
    map<string, int>::iterator it;

    for (int i = 0; i < ids.size(); i++) {

        it = seqIndex.find(ids[i]);
        if (it != seqIndex.end()) {
            indexes[i] = it->second;
        }else{
            string message = "[ERROR]: The dataset does not contain a ";
            message += "sequence named " + ids[i] + ".\n";
            RcppThread::Rcerr << endl << message << endl;
            throw Rcpp::exception(message.c_str());
        }
    }
    return indexes;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getList(string type) {
    if (hasBinTable(type)) {

        return binTables[getBinTableIndex(type)].getList(names);

    }else if ((type == "asv") && (hasSequenceData)) {
        vector<string> seqIds = getSequenceNames();
        vector<string> asvIds(seqIds.size(), "");
        for (int i = 0; i < asvIds.size(); i++) {
            asvIds[i] = "ASV" + toString(i+1);
        }
        assignBins(asvIds, nullIntVector, nullVector, seqIds, "asv");
        return binTables[getBinTableIndex(type)].getList(names);
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
vector<string> Dataset::getListVector(string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getListVector(names);
    }

    return nullVector;
}
/******************************************************************************/
vector<string> Dataset::getSequenceNames(string sample){
    vector<string> included;

    // get all "good" names in dataset
    if (sample == "")  {
        for (int i = 0; i < tableSeqs.size(); i++) {
            if (tableSeqs[i]) {
                included.push_back(names[i]);
            }
        }
    // get all "good" names in specific sample
    }else {
        if (count.hasSample(sample)) {
            // all seqs
            for (int i = 0; i < tableSeqs.size(); i++) {
                // if "good" seq
                if (tableSeqs[i]) {
                    // if this sequence is in sample
                    if (count.hasSample(sample, i)) {
                        included.push_back(names[i]);
                    }
                }
            }
        }
    }

    return included;
}
/******************************************************************************/
vector<vector<string> > Dataset::getSequenceNamesBySample(vector<string> samples){
    vector<vector<string> > result;

    // return all samples if none specified
    if (samples.size() == 0) {
        samples = getSamples();
    }

    for (int i = 0; i < samples.size(); i++) {
        result.push_back(getSequenceNames(samples[i]));
    }

    return result;
}
/******************************************************************************/
int Dataset::getNumSamples() {
    if (binTables.size() != 0) {
        return binTables[0].getNumSamples();
    }

    return count.getNumSamples();
}
/******************************************************************************/
int Dataset::getNumTreatments() {
    if (binTables.size() != 0) {
        return binTables[0].getNumTreatments();
    }

    return count.getNumTreatments();
}
/******************************************************************************/
// TODO document in module_exports.R
int Dataset::getNumBins(string type) {
    int numBins = 0;

    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].numBins;
    }

    return numBins;
}
/******************************************************************************/
// total abundance for a given binID
int Dataset::getBinAbundance(string binId, string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getAbundance(binId, "");
    }
    return 0;
}
/******************************************************************************/
// abundances for given binID broken down by sample
vector<int> Dataset::getBinAbundances(string binId, string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getAbundances(binId);
    }
    return nullIntVector;
}
/******************************************************************************/
// string containing sequence names for given binID
string Dataset::getBin(string binId, string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].get(binId, names);
    }
    return "";
}
/******************************************************************************/
vector<string> Dataset::getBinIds(string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getIds();
    }
    return nullVector;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getBinTaxonomyReport(string type) {
    return (fillTaxReport(type));
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getRAbund(string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getRAbund();
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
vector<int> Dataset::getRAbundVector(string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getRAbundVector();
    }
    return nullIntVector;
}
/******************************************************************************/
// sample functions
vector<string> Dataset::getSamples(){
    // maybe we just provided an binTable with no seqIds
    if (binTables.size() != 0) {
        return binTables[0].getSamples();
    }
    return count.getSamples();
}
/******************************************************************************/
vector<int> Dataset::getSampleTotals(){
    if (binTables.size() != 0) {
        return binTables[0].getSampleTotals();
    }
    return count.getSampleTotals();
}
/******************************************************************************/
// id, trashCode
Rcpp::DataFrame Dataset::getScrapReport(string mode) {

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
// trashCode, uniqueCount, totalCount
Rcpp::List Dataset::getScrapSummary() {

    Rcpp::List list = Rcpp::List::create();
    vector<string> listNames;

    if (badAccnos.size() != 0) {
        vector<string> codes(badAccnos.size(), "");
        vector<int> uniqueCounts(badAccnos.size(), 0);
        vector<int> totalCounts(badAccnos.size(), 0);

        int index = 0;
        for (auto it = badAccnos.begin(); it != badAccnos.end(); it++) {
            codes[index] = it->first;
            uniqueCounts[index] = it->second[0];
            totalCounts[index] = it->second[1];
            index++;
        }

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("trash_code") = codes,
            Rcpp::_["unique_count"] = uniqueCounts,
            Rcpp::_["total_count"] = totalCounts);

        list.push_back(df);
        listNames.push_back("sequence_scrap_summary");
    }

    if (binTables.size() != 0) {
        for (int i = 0; i < binTables.size(); i++) {
            list.push_back(binTables[i].getScrapSummary());
            string tag = binTables[i].label + "_scrap_summary";
            listNames.push_back(tag);
        }
    }

    list.attr("names") = listNames;

    return list;
}
/******************************************************************************/
// ids, abundances, sample(optional), treatment(optional)
// This table represents mothur's count and design files.
Rcpp::DataFrame Dataset::getSequenceAbundanceTable() {
    return count.getAbundanceTable(select(names, tableSeqs),
                                            getIncludedNamesIndexes());
}
/******************************************************************************/
// total abundance for each sequence
vector<int> Dataset::getSequenceAbundances(){
    return count.getTotalAbundances(getIncludedNamesIndexes());
}
/******************************************************************************/
vector<vector<int>> Dataset::getSeqsAbundsBySample(){
    return count.getAbundances(getIncludedNamesIndexes());
}
/******************************************************************************/
vector<string> Dataset::getSequences(string sample){
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
vector<vector<string> > Dataset::getSequencesBySample(vector<string> samples){
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
Rcpp::DataFrame Dataset::getSequenceReport(){

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
Rcpp::List Dataset::getSequenceSummary() {

    Rcpp::List result = Rcpp::List::create();

    if (!hasSeqs()) {
        return result;
    }

    vector<string> result_names;

    Summary* summary = new Summary(processors);

    vector<vector<int> > report;

    report.push_back(select(starts, tableSeqs));
    report.push_back(select(ends, tableSeqs));
    report.push_back(select(lengths, tableSeqs));
    report.push_back(select(ambigs, tableSeqs));
    report.push_back(select(polymers, tableSeqs));
    report.push_back(select(numns, tableSeqs));

    Rcpp::DataFrame seqResults = summary->summarizeFasta(
        report, count.getTotalAbundances(getIncludedNamesIndexes()));

    result.push_back(seqResults);
    result_names.push_back("sequence_summary");

    delete summary;

    Rcpp::DataFrame scrap = getScrapSummary();
    if (scrap.size() != 0) {
        result.push_back(scrap);
        result_names.push_back("scrap_summary");
    }

    result.attr("names") = result_names;
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
Rcpp::DataFrame Dataset::getShared(string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getShared();
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
vector<vector<int> > Dataset::getSharedVector(string type) {
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].getSharedVector();
    }
    return null2DIntVector;
}
/******************************************************************************/
vector<string> Dataset::getTreatments(){
    // maybe we just provided an binTable with no seqIds
    if (binTables.size() != 0) {
        return binTables[0].getTreatments();
    }
    return count.getTreatments();
}
/******************************************************************************/
vector<int> Dataset::getTreatmentTotals(){
    // maybe we just provided an binTable with no seqIds
    if (binTables.size() != 0) {
        return binTables[0].getTreatmentTotals();
    }
    return count.getTreatmentTotals();
}
/******************************************************************************/
long long Dataset::getTotal(string sample){

    if (binTables.size() != 0) {
        return binTables[0].getTotal(sample);
    }

    return count.getTotal(sample);
}
/******************************************************************************/
long long Dataset::getUniqueTotal(string sample){
    if (sample == "") {
        return numUnique;
    }
    return getSequenceNames(sample).size();
}
/******************************************************************************/
bool Dataset::hasBinTable(string type) {
    for (int i = 0; i < binTables.size(); i++) {
        if (binTables[i].label == type) {
            return true;
        }
    }
    return false;
}
/******************************************************************************/
int Dataset::getBinTableIndex(string type) {
    for (int i = 0; i < binTables.size(); i++) {
        if (binTables[i].label == type) {
            return i;
        }
    }
    return -1;
}
/******************************************************************************/
bool Dataset::hasSample(string sample){
    if (binTables.size() != 0) {
        return binTables[0].hasSample(sample);
    }
    return count.hasSample(sample);
}
/******************************************************************************/
bool Dataset::hasListAssignments(string type){
    if (hasBinTable(type)) {
        return binTables[getBinTableIndex(type)].hasListAssignments;
    }
    return false;
}
/******************************************************************************/
bool Dataset::hasSeqs() {
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
void Dataset::mergeSequences(vector<string> ids, string reason){
    if (ids.size() != 1) {

        vector<int> indexes = getIndexes(ids);

        // sanity check: if you have assigned bins, make sure the sequences
        // are in the same bin
        if (binTables.size() != 0) {
            // innocent until proven guilty
            vector<bool> okToMerge(binTables.size(), true);

            int tableIndex = 0;
            string message = "[ERROR]: can not merge sequences assigned";
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
void Dataset::mergeBins(vector<string> ids, string reason, string type){
    if (ids.size() != 1) {
        if (hasBinTable(type)) {
            binTables[getBinTableIndex(type)].merge(ids);
        }
    }
}
/******************************************************************************/
void Dataset::removeLineages(vector<string> contaminants, string trashTag) {
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

        vector<string> binIds = binTables[i].getIds();
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
        binTables[i].updateTotals();
    }

    count.updateTotals();
}
/******************************************************************************/
void Dataset::removeBins(vector<string> namesToRemove,
                              vector<string> trashTags, string type){

    if (binTables.size() == 0) { return; }

    if (hasBinTable(type)) {

        if (namesToRemove.size() != trashTags.size()) {
            string message = "[ERROR]: Size mismatch. You must provide a trash";
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
                    if (binTables[j].hasListAssignments) {
                        for (int seq : seqsToRemove) {
                            binTables[j].remove(seq, count, trashTags[i], true);
                        }
                    }else {
                        binTables[j].remove(namesToRemove[i], trashTags[i]);
                    }
                }
                binTables[j].updateTotals();
            }
        }

        count.updateTotals();
    }
}
/******************************************************************************/
void Dataset::removeSamples(vector<string> samples) {
    count.removeSamples(samples);

    // // remove samples from count
    // remove any seqs only assigned to these samples
    for (int i = 0; i < names.size(); i++) {
        // included seq
        if (tableSeqs[i]) {
            if (sum(count.getAbundances(i)) == 0) {
                removeSequence(i, "removedSamples", true, false);

                // remove from list otus
                for (int j = 0; j < binTables.size(); j++) {
                    if (binTables[j].hasListAssignments) {
                        binTables[j].remove(i, count, "removedSamples", false);
                    }
                }
            }
        }
    }

    // // remove samples from binTables
    for (int i = 0; i < binTables.size(); i++) {
        binTables[i].removeSamples(samples);
        binTables[i].updateTotals();
    }
    count.updateTotals();
}
/******************************************************************************/
void Dataset::removeSequence(int index, string reasons,
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
                if (binTables[i].hasListAssignments) {
                    binTables[i].remove(index, count, reasons, update);
                }
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
            vector<int> badAbunds(2, 1);
            badAbunds[1] = abund;
            badAccnos[theseReasons[j]] = badAbunds;
        }
    }

    // update uniqueBad
    uniqueBad++;
}
/******************************************************************************/
void Dataset::removeSequences(vector<string> namesToRemove,
                         vector<string> trashTags){

    if (namesToRemove.size() != trashTags.size()) {
        string message = "[ERROR]: Size mismatch. You must provide a trash";
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

    if (binTables.size() != 0) {
        for (int i = 0; i < binTables.size(); i++) {
            if (binTables[i].hasListAssignments) {
                binTables[i].updateTotals();
            }
        }
    }
}
/******************************************************************************/
// for datasets without samples
void Dataset::setAbundance(vector<string> n, vector<int> abunds,
                            string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";

    for (int i = 0; i < n.size(); i++) {
        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            if (abunds[i] == 0) {
                removeSequence(index, reason, true, true);
            }else{
                count.setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count.updateTotals();

    // update list assignments if they include sequence ids
    for (int i = 0; i < binTables.size(); i++) {
        if (binTables[i].hasListAssignments) {
            binTables[i].assignAbundance(nullVector, nullIntVector,
                                         nullVector, nullIntVector,
                                         count, true);
        }
    }
}
/******************************************************************************/
// for datasets with samples
void Dataset::setAbundances(vector<string> n, vector<vector<int>> abunds,
                       string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";

    for (int i = 0; i < n.size(); i++) {
        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            if (sum(abunds[i]) == 0) {
                removeSequence(index, reason, true, true);
            }else{
                count.setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count.updateTotals();

    // update list assignments if they include sequence ids
    for (int i = 0; i < binTables.size(); i++) {
        if (binTables[i].hasListAssignments) {
            binTables[i].assignAbundance(nullVector, nullIntVector,
                                         nullVector, nullIntVector,
                                         count, true);
        }
    }
}
/******************************************************************************/
// for datasets without samples
void Dataset::setBinAbundance(vector<string> binIDS, vector<int> abunds,
                     string reason, string type) {
    if (hasBinTable(type)) {

        if (binTables[getBinTableIndex(type)].hasListAssignments) {
            string message = "[WARNING]: cannot assign bin abundance for bins ";
            message += " with sequences assigned to them. Doing so could cause ";
            message += "inconsistencies, ignoring bin assignments.";
            RcppThread::Rcout << endl << message << endl;
            return;
        }

        vector<int> seqsToRemove = binTables[getBinTableIndex(type)].setAbundance(binIDS, abunds,
                                                             reason);
        // remove any sequences from removed bin from dataset
        for (int seq : seqsToRemove) {
            removeSequence(seq, reason, true, false);
        }

        // remove from other lists
        for (int i = 0; i < binTables.size(); i++) {
            if (binTables[i].hasListAssignments && (binTables[i].label != type)) {
                for (int seq : seqsToRemove) {
                    binTables[i].remove(seq, count, reason, true);
                }
            }
            binTables[i].updateTotals();
        }
    }
}
/******************************************************************************/
// for datasets with samples
void Dataset::setBinAbundances(vector<string> binIDS,
                               vector<vector<int> > abunds,
                               string reason, string type) {
    if (hasBinTable(type)) {

        if (binTables[getBinTableIndex(type)].hasListAssignments) {
            string message = "[WARNING]: cannot assign bin abundance for bins ";
            message += " with sequences assigned to them. Doing so could cause ";
            message += "inconsistencies, ignoring bin assignments.";
            RcppThread::Rcout << endl << message << endl;
            return;
        }

        vector<int> seqsToRemove = binTables[getBinTableIndex(type)].
        setAbundances(binIDS, abunds, reason);

        // remove any sequences from removed bin
        for (int seq : seqsToRemove) {
            removeSequence(seq, reason, true, false);
        }

        // remove from other lists
        for (int i = 0; i < binTables.size(); i++) {
            if (binTables[i].hasListAssignments && (binTables[i].label != type)) {
                for (int seq : seqsToRemove) {
                    binTables[i].remove(seq, count, reason, true);
                }
            }
            binTables[i].updateTotals();
        }
    }
}
/******************************************************************************/
void Dataset::setSequences(vector<string> n, vector<string> s,
                      vector<string> c){
    if (n.size() != s.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    bool hasComments = false;
    if (c.size() != 0) {
        if (c.size() != n.size()) {
            string message = "[ERROR]: Size mismatch. When providing comments,";
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
Rcpp::RawVector Dataset::serializeDataset() {
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


