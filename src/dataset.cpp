
#include "../inst/include/rdataset.h"
#include "seqreport.h"
#include "summary.h"

/******************************************************************************/
Dataset::Dataset(string n, int proc) : datasetName(n) {
    isAligned = false;
    hasContigsData = false;
    hasAlignData = false;
    numSamples = 0;
    numTreatments = 0;
    numUnique = 0;
    totalBad = 0;
    uniqueBad = 0;
    alignmentLength = 0;
    processors = proc;
    count = new SeqAbundTable();
}
/******************************************************************************/
Dataset::~Dataset() {
    delete count;
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
    numSamples = 0;
    numTreatments = 0;
    numUnique = 0;
    totalBad = 0;
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
}
/******************************************************************************/
Rcpp::List Dataset::exportDataset(){
    Rcpp::List result;
    // TODO
    return result;
}
/******************************************************************************/
void Dataset::addSeqs(vector<string> n, vector<string> s, vector<string> c) {

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
    count->addSeqs(countNames);
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

    // add to "good" sequence count
    numUnique += names.size();

   // add calcs for starts, ends, lengths, ambigs, polymers, numns
    SeqReport report;
    report.addReports(s, starts, ends, lengths, ambigs, polymers, numns);
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
        RcppThread::Rcout << endl << message;
        return;
    }

    vector<int> idIndexes = getIndexes(ids);

    count->assignSequenceAbundance(idIndexes, abunds, samples, treatments);

    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
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
            RcppThread::Rcout << endl << message;
        }
    }
    return indexes;
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
                RcppThread::Rcout << endl << message;
                return;
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
                RcppThread::Rcout << endl << message;
                return;
            }
        }
    }
}
/******************************************************************************/
// ids, abundances, sample(optional), treatment(optional)
// This table represents mothur's count and design files.
Rcpp::DataFrame Dataset::getSequenceAbundanceTable() {
    return count->getSequenceAbundanceTable(select(names, tableSeqs),
                                            getIncludedNamesIndexes());
}
/******************************************************************************/
// sample functions
vector<string> Dataset::getSamples(){
    return count->getSamples();
}
/******************************************************************************/
vector<string> Dataset::getTreatments(){
    return count->getTreatments();
}
/******************************************************************************/
vector<int> Dataset::getSampleTotals(){
    return count->getSampleTotals();
}
/******************************************************************************/
vector<int> Dataset::getTreatmentTotals(){
    return count->getTreatmentTotals();
}
/******************************************************************************/
long long Dataset::getTotal(string sample){
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
vector<string> Dataset::getNames(string sample){
    vector<string> included;

    if (sample == "") {
        for (int i = 0; i < tableSeqs.size(); i++) {
            if (tableSeqs[i]) {
                included.push_back(names[i]);
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
vector<string> Dataset::getSeqs(string sample){
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
vector<vector<string> > Dataset::getSeqsBySample(vector<string> samples){
    vector<vector<string> > result;
    // TODO
    return result;
}
/******************************************************************************/
void Dataset::removeSeqs(vector<string> names, vector<string> trashTags){
    // TODO
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
        report, count->getSeqsAbunds(getIncludedNamesIndexes()));

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
            report, count->getSeqsAbunds(getIncludedNamesIndexes()));

        result.push_back(contigsResults);
        result_names.push_back("contigs_summary");
    }

    if (hasAlignData) {
        vector<vector<float> > report;
        report.push_back(select(searchScore, tableSeqs));
        report.push_back(select(simScore, tableSeqs));

        Rcpp::DataFrame alignResults = summary->summarizeAlign(
            report, select(longestInsert, tableSeqs),
            count->getSeqsAbunds(getIncludedNamesIndexes()));

        result.push_back(alignResults);
        result_names.push_back("align_summary");
    }

    delete summary;

    result.attr("names") = result_names;
    return result;
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
void Dataset::mergeSeqs(vector<string> names, string reason, string sample){
    // TODO
}
/******************************************************************************/
int Dataset::getAbund(string name, string sample){
    int abund = 0;
    auto it = seqIndex.find(name);
    if (it != seqIndex.end()) {
        abund = count->getAbund(it->second, sample);
    }
    return abund;
}
/******************************************************************************/
vector<int> Dataset::getAbunds(string name){
    vector<int> abunds;
    auto it = seqIndex.find(name);
    if (it != seqIndex.end()) {
       abunds = count->getAbunds(it->second);
    }
    return abunds;
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
// total abundance for each sequence
vector<int> Dataset::getSeqsAbunds(){
    return count->getSeqsAbunds(getIncludedNamesIndexes());
}
/******************************************************************************/
vector<vector<int>> Dataset::getSeqsAbundsBySample(){
    vector<vector<int> > result;
    // TODO
    return result;
}
/******************************************************************************/
void Dataset::setSeqs(vector<string> n, vector<string> s,
             vector<string> c){
    names = n;
    seqs = s;
    if (c.size() == 0) {
        c.resize(names.size(), "");
    }
    comments = c;
}
/******************************************************************************/

void Dataset::setAbundances(vector<string> names, vector<int> abunds,
                       string reason){
    // TODO
}
/******************************************************************************/
int Dataset::getAlignedLength(vector<string> seqs) {
    int length = -1;


    return length;
}
/******************************************************************************/


