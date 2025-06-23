
#include "../inst/include/rdataset.h"
#include "seqreport.h"
#include "summary.h"

/******************************************************************************/
Dataset::Dataset(string n, int proc) : datasetName(n) {
    isAligned = false;
    hasContigsData = false;
    hasAlignData = false;
    hasOtuData = false;
    hasSequenceData = false;
    numSamples = 0;
    numTreatments = 0;
    numOtus = 0;
    label = "";
    numUnique = 0;
    uniqueBad = 0;
    alignmentLength = 0;
    processors = proc;
    count = new AbundTable();
    otuTable = nullptr;
}
/******************************************************************************/
Dataset::~Dataset() {
    delete count;
    if (otuTable != nullptr) { delete otuTable; }
}
/******************************************************************************/
SEXP Dataset::getPointer() {
    Rcpp::XPtr<Dataset> ptr(this);
    return ptr;
}
/******************************************************************************/
void Dataset::clear() {
    isAligned = false;
    hasContigsData = false;
    hasAlignData = false;
    hasSequenceData = false;
    numSamples = 0;
    numTreatments = 0;
    numOtus = 0;
    label = "";
    numUnique = 0;
    uniqueBad = 0;
    alignmentLength = 0;

    // sequence data
    names.clear();
    seqs.clear();
    comments.clear();
    trashCodes.clear();

    // contigs report
    olengths.clear();
    ostarts.clear();
    oends.clear();
    mismatches.clear();
    ee.clear();

    // sequence summary data
    starts.clear();
    ends.clear();
    lengths.clear();
    ambigs.clear();
    polymers.clear();
    numns.clear();

    // alignment report
    searchScore.clear();
    simScore.clear();
    longestInsert.clear();

    // sequence taxonomy assignments
    taxonomies.clear();

    // maps sequence name to index in vectors
    seqIndex.clear();
    tableSeqs.clear();

    badAccnos.clear();

    count->clear();

    // if you have an otuTable then otuTable.clear();
    if (hasOtuData) { otuTable->clear(); hasOtuData = false; }
}
/******************************************************************************/
Rcpp::List Dataset::exportDataset(){
    Rcpp::List result;
    // TODO
    return result;
}
/******************************************************************************/
void Dataset::addSequences(vector<string> n, vector<string> s, vector<string> c) {

    // must provide the same number of names and seqs
    if (n.size() != s.size()) { return; }

    // add to seqIndex
    int numSeqs = names.size();
    vector<int> countNames;
    for (int i = 0; i < n.size(); i++) {
        countNames.push_back(numSeqs);
        seqIndex[n[i]] = numSeqs;
        numSeqs++;
    }

    // add to count
    count->add(countNames);
    countNames.clear();

    // add to names
    names.insert(names.end(), n.begin(), n.end());
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

    // add to "good" sequence count - giving preference to otuTable
    if (!hasOtuData) {
        numUnique += names.size();
    }

    hasSequenceData = true;
}
/******************************************************************************/
void Dataset::assignOtuAbundance(string label, vector<string> otuIds,
                        vector<int> abunds, vector<string> samples,
                        vector<string> seqIds) {

    if (otuIds.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. otu_ids and abunds must be";
        message += " the same size.";
        throw Rcpp::exception(message.c_str());
    }

    hasOtuData = true;

    // new table
    if (otuTable == nullptr) {
        otuTable = new OtuTable(label);
    }else{
        // assume they had "list" data and are adding "shared" data
        otuTable->setLabel(label);
    }

    // sanity checks - R6 object passes blank strings if seqIds == NULL
    if (!seqIds.empty()) {
        string id = "";
        if (allIdentical(seqIds, id)) {
            if (id == "") { seqIds.clear(); }
        }
    }

    // sanity checks - R6 object passes blank strings if samples == NULL
    if (!samples.empty()) {
        string id = "";
        if (allIdentical(samples, id)) {
            if (id == "") { samples.clear(); }
        }
    }

    // add otus
    otuTable->add(otuIds, abunds, samples, seqIds);

    // sanity check unique totals
    if (numUnique != 0) {
        // numUnique == -1 in otuTable when no seqIds are given. This is because
        // without seqIds we assume the abundances are not unique.
        if ((otuTable->numUnique != numUnique) && (otuTable->numUnique != -1)){
            string message = "[ERROR]: The number of sequences in your dataset";
            message += " is " + toString(numUnique) + ", but the number of ";
            message += " sequences in your otu table is ";
            message += toString(otuTable->numUnique) + ". Removing otu data.";
            RcppThread::Rcout << endl << message << endl;
            otuTable->clear();
        }
    }else {
        numUnique = otuTable->numUnique;
    }

    numOtus = otuTable->numOtus;
    label = otuTable->label;
    numSamples = otuTable->getNumSamples();
    numTreatments = otuTable->getNumTreatments();
}
/******************************************************************************/
void Dataset::assignTreatments(vector<string> samples, vector<string> treatments) {
    if (hasOtuData) {
        otuTable->assignTreatments(samples, treatments);
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
        return;
    }
    count->assignTreatments(samples, treatments);
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
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

    count->assignAbundance(idIndexes, abunds, samples, treatments);

    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
}
/******************************************************************************/
// align_seqs will create searchScores, simScores and longestInserts
void Dataset::addAlignReport(vector<string>& n, vector<double>& ss,
                    vector<double>& sims,
                    vector<int>& li) {

    set<int> sizes;
    sizes.insert(n.size());
    sizes.insert(ss.size());
    sizes.insert(sims.size());
    sizes.insert(li.size());
    sizes.insert(numUnique);

    if (sizes.size() > 1) {
        string message = "[WARNING]: You must provide align report info for ";
        message += "each sequence in the dataset, ignoring align report data.";
        RcppThread::Rcout << endl << message << endl;
    }else {
        // create space
        searchScore.resize(names.size(), 0);
        simScore.resize(names.size(), 0);
        longestInsert.resize(names.size(), 0);
        hasAlignData = true;

        for (int i = 0; i < n.size(); i++) {

            auto it = seqIndex.find(n[i]);
            if (it != seqIndex.end()) {
                int index = it->second;
                searchScore[index] = ss[index];
                simScore[index] = sims[index];
                longestInsert[index] = li[index];
            }else{
                string message = "[ERROR]: The dataset does not contain a ";
                message += "sequence named " + n[i] + ", quitting.\n";
                throw Rcpp::exception(message.c_str());
            }
        }
    }

}
/******************************************************************************/
// make_contigs will create overlapLengths, overlapStarts, overlapEnds,
// mismatches, and expectedErrors
void Dataset::addContigsReport(vector<string>& n, vector<int>& ol,
                      vector<int>& os,
                      vector<int>& oe,
                      vector<int>& m,
                      vector<double>& e) {

    set<int> sizes;
    sizes.insert(n.size());
    sizes.insert(ol.size());
    sizes.insert(os.size());
    sizes.insert(oe.size());
    sizes.insert(m.size());
    sizes.insert(e.size());
    sizes.insert(numUnique);

    if (sizes.size() > 1) {
        string message = "[WARNING]: You must provide contigs report info for ";
        message += "each sequence in the dataset, ignoring contigs report ";
        message += "data.";
        RcppThread::Rcout << endl << message << endl;
    }else {
        hasContigsData = true;
        // create space
        olengths.resize(names.size(), 0);
        ostarts.resize(names.size(), 0);
        oends.resize(names.size(), 0);
        mismatches.resize(names.size(), 0);
        ee.resize(names.size(), 0);

        for (int i = 0; i < n.size(); i++) {

            auto it = seqIndex.find(n[i]);
            if (it != seqIndex.end()) {
                int index = it->second;
                olengths[index] = ol[index];
                ostarts[index] = os[index];
                oends[index] = oe[index];
                mismatches[index] = m[index];
                ee[index] = e[index];
            }else{
                string message = "[ERROR]: The dataset does not contain a ";
                message += "sequence named " + n[i] + ", quitting.\n";
                throw Rcpp::exception(message.c_str());
            }
        }
    }
}
/******************************************************************************/
int Dataset::getAbundance(string name, string sample){
    int abund = 0;
    auto it = seqIndex.find(name);
    if (it != seqIndex.end()) {
        if (tableSeqs[it->second]) {
            abund = count->getAbundance(it->second, sample);
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
            abunds = count->getAbundances(it->second);
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
    }else{
        isAligned = false;
        alignmentLength = -1;
    }
    return alignmentLength;
}
/******************************************************************************/
// search_score, sim_score, longest_insert
Rcpp::DataFrame Dataset::getAlignReport(){

    if (hasAlignData) {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("id") = select(names, tableSeqs),
            Rcpp::_["search_score"] = select(searchScore, tableSeqs),
            Rcpp::_["sim_score"] = select(simScore, tableSeqs),
            Rcpp::_["longest_insert"] = select(longestInsert, tableSeqs));

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// contigs sumary data: olengths, ostarts, oends, mismatches, ee
// report[0] = length, report[1] = overlap_length, report[2] = overlap_start,
// report[3] = overlap_end, report[4] = mismatches, report[5] = num_ns,
// report[6] = ee
Rcpp::DataFrame Dataset::getContigsReport(){

    if (hasContigsData){
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("id") = select(names, tableSeqs),
            Rcpp::_["length"] = select(lengths, tableSeqs),
            Rcpp::_["overlap_length"] = select(olengths, tableSeqs),
            Rcpp::_["overlap_start"] = select(ostarts, tableSeqs),
            Rcpp::_["overlap_end"] = select(oends, tableSeqs),
            Rcpp::_["mismatches"] = select(mismatches, tableSeqs),
            Rcpp::_["num_ns"] = select(numns, tableSeqs),
            Rcpp::_["ee"] = select(ee, tableSeqs));

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
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
            throw Rcpp::exception(message.c_str());
        }
    }
    return indexes;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getList() {
    if (hasOtuData) {
        return otuTable->getList();
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
vector<string> Dataset::getNames(string sample){
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
        if (count->hasSample(sample)) {
            // all seqs
            for (int i = 0; i < tableSeqs.size(); i++) {
                // if "good" seq
                if (tableSeqs[i]) {
                    // if this sequence is in sample
                    if (count->hasSample(sample, i)) {
                        included.push_back(names[i]);
                    }
                }
            }
        }
    }

    return included;
}
/******************************************************************************/
vector<vector<string> > Dataset::getNamesBySample(vector<string> samples){
    vector<vector<string> > result;
    // TODO
    return result;
}
/******************************************************************************/
// total abundance for a given outID, optional sample
int Dataset::getOtuAbundance(string otuId, string sample) {
    if (hasOtuData) {
        return otuTable->getAbundance(otuId, sample);
    }
    return 0;
}
/******************************************************************************/
// abundances for given otuID broken down by sample
vector<int> Dataset::getOtuAbundances(string otuId) {
    if (hasOtuData) {
        return otuTable->getAbundances(otuId);
    }
    return nullIntVector;
}
/******************************************************************************/
// string containing sequence names for given otuID
string Dataset::getOtu(string otuId) {
    if (hasOtuData) {
        return otuTable->get(otuId);
    }
    return "";
}
/******************************************************************************/
vector<string> Dataset::getOtuIds() {
    if (hasOtuData) {
        return otuTable->getOtuIds();
    }
    return nullVector;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getRAbund() {
    if (hasOtuData) {
        return otuTable->getRAbund();
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// sample functions
vector<string> Dataset::getSamples(){
    if (hasOtuData) {
       return otuTable->getSamples();
    }
    return count->getSamples();
}
/******************************************************************************/
vector<int> Dataset::getSampleTotals(){
    if (hasOtuData) {
        return otuTable->getSampleTotals();
    }
    return count->getSampleTotals();
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
    }else if (mode == "otu") {
        if (hasOtuData) {
            return otuTable->getScrapReport();
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

    if (hasOtuData) {
        list.push_back(otuTable->getScrapSummary());
        listNames.push_back("otu_scrap_summary");
    }

    list.attr("names") = listNames;

    return list;
}
/******************************************************************************/
// ids, abundances, sample(optional), treatment(optional)
// This table represents mothur's count and design files.
Rcpp::DataFrame Dataset::getSequenceAbundanceTable() {
    return count->getAbundanceTable(select(names, tableSeqs),
                                            getIncludedNamesIndexes());
}
/******************************************************************************/
// total abundance for each sequence
vector<int> Dataset::getSequenceAbundances(){
    return count->getTotalAbundances(getIncludedNamesIndexes());
}
/******************************************************************************/
vector<vector<int>> Dataset::getSeqsAbundsBySample(){
    vector<vector<int> > result;
    // TODO
    return result;
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
        if (count->hasSample(sample)) {
            // all seqs
            for (int i = 0; i < tableSeqs.size(); i++) {
                // if "good" seq
                if (tableSeqs[i]) {
                    // if this sequence is in sample
                    if (count->hasSample(sample, i)) {
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
    // TODO
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

    if (!hasSequenceData) {
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
        report, count->getTotalAbundances(getIncludedNamesIndexes()));

    result.push_back(seqResults);
    result_names.push_back("sequence_summary");

    if (hasContigsData) {

        vector<vector<int> > report;
        report.push_back(select(lengths, tableSeqs));
        report.push_back(select(olengths, tableSeqs));
        report.push_back(select(ostarts, tableSeqs));
        report.push_back(select(oends, tableSeqs));
        report.push_back(select(mismatches, tableSeqs));
        report.push_back(select(numns, tableSeqs));

        Rcpp::DataFrame contigsResults = summary->summarizeContigs(
            report, count->getTotalAbundances(getIncludedNamesIndexes()));

        result.push_back(contigsResults);
        result_names.push_back("contigs_summary");
    }

    if (hasAlignData) {
        vector<vector<float> > report;
        report.push_back(select(searchScore, tableSeqs));
        report.push_back(select(simScore, tableSeqs));

        Rcpp::DataFrame alignResults = summary->summarizeAlign(
            report, select(longestInsert, tableSeqs),
            count->getTotalAbundances(getIncludedNamesIndexes()));

        result.push_back(alignResults);
        result_names.push_back("align_summary");
    }

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
vector<string> Dataset::getTreatments(){
    if (hasOtuData) {
        return otuTable->getTreatments();
    }
    return count->getTreatments();
}
/******************************************************************************/
vector<int> Dataset::getTreatmentTotals(){
    if (hasOtuData) {
        return otuTable->getTreatmentTotals();
    }
    return count->getTreatmentTotals();
}
/******************************************************************************/
long long Dataset::getTotal(string sample){
    if (hasOtuData) {
        return otuTable->getTotal(sample);
    }
    return count->getTotal(sample);
}
/******************************************************************************/
long long Dataset::getUniqueTotal(string sample){
    if (sample == "") {
        return numUnique;
    }
    return getNames(sample).size();
}
/******************************************************************************/
bool Dataset::hasSample(string sample){
    return count->hasSample(sample);
}
/******************************************************************************/
void Dataset::mergeSequences(vector<string> ids, string reason){
    if (ids.size() != 1) {

        vector<int> indexes = getIndexes(ids);
        count->merge(indexes);

        for (int i = 1; i < indexes.size(); i++) {
            // no need to update the sample and treatment counts
            removeSequence(indexes[i], reason, false);
        }
    }
}
/******************************************************************************/
void Dataset::mergeOtus(vector<string> ids, string reason){
    if (ids.size() != 1) {
        if (hasOtuData) {
            otuTable->merge(ids);
            numOtus = otuTable->numOtus;
        }
    }
}
/******************************************************************************/
void Dataset::removeOtus(vector<string> namesToRemove,
                              vector<string> trashTags){

    if (!hasOtuData) { return; }

    if (namesToRemove.size() != trashTags.size()) {
        string message = "[ERROR]: Size mismatch. You must provide a trash";
        message += " code for each otu.";
        throw Rcpp::exception(message.c_str());
    }

    for (int i = 0; i < namesToRemove.size(); i++) {
        otuTable->remove(namesToRemove[i], trashTags[i]);
    }

    otuTable->updateTotals();
    numSamples = otuTable->getNumSamples();
    numTreatments = otuTable->getNumTreatments();
    numOtus = otuTable->numOtus;

    if (otuTable->numUnique != -1) {
        numUnique = otuTable->numUnique;
    }
}

/******************************************************************************/
void Dataset::removeSequence(int index, string reasons, bool update) {
    // remove from tableSeqs and add trashCode
    tableSeqs[index] = false;
    trashCodes[index] += reasons + ",";
    numUnique--;

    // remove from counts
    int abund = 1;
    if (update) {
        abund = count->remove(index);
    }else{
        abund = count->getAbundance(index);
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
        throw Rcpp::exception(message.c_str());
    }

    for (int i = 0; i < namesToRemove.size(); i++) {
        auto it = seqIndex.find(namesToRemove[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            removeSequence(index, trashTags[i]);
        }else{
            string message = "[WARNING]: " + namesToRemove[i] + " is not in ";
            message += "your dataset, ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count->updateTotals();
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();

}
/******************************************************************************/
// for datasets without samples
void Dataset::setAbundance(vector<string> n, vector<int> abunds,
                            string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";

    for (int i = 0; i < n.size(); i++) {
        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            if (abunds[i] == 0) {
                removeSequence(index, reason);
            }else{
                count->setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count->updateTotals();
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
}
/******************************************************************************/
// for datasets with samples
void Dataset::setAbundances(vector<string> n, vector<vector<int>> abunds,
                       string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";

    for (int i = 0; i < n.size(); i++) {
        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            if (sum(abunds[i]) == 0) {
                removeSequence(index, reason);
            }else{
                count->setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count->updateTotals();
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
}
/******************************************************************************/
// for datasets without samples
void Dataset::setOtuAbundance(vector<string> otuIDS, vector<int> abunds,
                     string reason) {
    if (hasOtuData) {
        otuTable->setAbundance(otuIDS, abunds, reason);

        otuTable->updateTotals();
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
        numOtus = otuTable->numOtus;

        if (otuTable->numUnique != -1) {
            numUnique = otuTable->numUnique;
        }
    }
}
/******************************************************************************/
// for datasets with samples
void Dataset::setOtuAbundances(vector<string> otuIDS,
                               vector<vector<int> > abunds,
                               string reason) {
    if (hasOtuData) {
        otuTable->setAbundances(otuIDS, abunds, reason);

        otuTable->updateTotals();
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
        numOtus = otuTable->numOtus;

        if (otuTable->numUnique != -1) {
            numUnique = otuTable->numUnique;
        }
    }
}
/******************************************************************************/
void Dataset::setSequences(vector<string> n, vector<string> s,
                      vector<string> c){
    if (n.size() != s.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
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


