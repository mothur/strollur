
#include "../inst/include/rdataset.h"
#include "seqreport.h"
#include "summary.h"

/******************************************************************************/
Dataset::Dataset(string n, int proc) : datasetName(n) {
    isAligned = false;
    numGroups = 0;
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
    numGroups = 0;
    numUnique = 0;
    totalBad = 0;
    uniqueBad = 0;

    // fasta data
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

    // fasta summary data
    starts.clear();
    ends.clear();
    lengths.clear();
    ambigs.clear();
    polymers.clear();
    numns.clear();

    // alignment report
    search_score.clear();
    sim_score.clear();
    longest_insert.clear();

    // sequence taxonomy assignments
    taxonomies.clear();

    // maps sequence name to index in vectors
    seqIndex.clear();

    bad_accnos.clear();

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
string Dataset::print(){
    string result = "";
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
// names, abundances, groups(optional)
// assumes same size, R side will check user inputs for us
void Dataset::assignSampleAbundance(vector<string> ids,
                           vector<int> abunds,
                           vector<string> samples) {

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

    count->assignSampleAbundance(idIndexes, abunds, samples);

    numGroups = count->getNumGroups();
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
void Dataset::addAlignReport(vector<string> n, vector<double> ss,
                    vector<double> sims,
                    vector<double> li) {}
/******************************************************************************/
// make_contigs will create overlapLengths, overlapStarts, overlapEnds,
// mismatches, and expectedErrors
void Dataset::addContigsReport(vector<string> n, vector<int> ol,
                      vector<int> os,
                      vector<int> oe,
                      vector<int> m,
                      vector<double> e) {}
/******************************************************************************/
// abundance functions
Rcpp::DataFrame Dataset::getSequenceAbundanceTable() {
    Rcpp::DataFrame results;
    // TODO
    return results;
}
/******************************************************************************/
// group functions
vector<string> Dataset::getGroups(string name){
    int index = -1;
    if (name != "") {
        index = seqIndex[name];
    }
    return count->getGroups(index);
}
/******************************************************************************/
vector<int> Dataset::getGroupTotals(string name){
    return count->getGroupTotals();
}
/******************************************************************************/
long long Dataset::getTotal(string group){
    return count->getTotal(group);
}
/******************************************************************************/
bool Dataset::hasGroup(string group){
    return count->hasGroup(group);
}
/******************************************************************************/
vector<string> Dataset::getNames(string group){
    vector<string> included;

    if (group == "") {
        for (int i = 0; i < tableSeqs.size(); i++) {
            if (tableSeqs[i]) {
                included.push_back(names[i]);
            }
        }
    }else{
        if (count->hasGroup(group)) {
            // all seqs
            for (int i = 0; i < tableSeqs.size(); i++) {
                // if "good" seq
                if (tableSeqs[i]) {
                    // if this sequence is in group
                    if (count->hasGroup(group, i)) {
                        included.push_back(names[i]);
                    }
                }
            }
        }
    }

    return included;
}
/******************************************************************************/
vector<vector<string> > Dataset::getNamesBySample(vector<string> group){
    vector<vector<string> > result;
    // TODO
    return result;
}
/******************************************************************************/
vector<string> Dataset::getSeqs(string group){
    vector<string> included;

    if (group == "") {
        for (int i = 0; i < tableSeqs.size(); i++) {
            if (tableSeqs[i]) {
                included.push_back(seqs[i]);
            }
        }
    }else{
        if (count->hasGroup(group)) {
            // all seqs
            for (int i = 0; i < tableSeqs.size(); i++) {
                // if "good" seq
                if (tableSeqs[i]) {
                    // if this sequence is in group
                    if (count->hasGroup(group, i)) {
                        included.push_back(seqs[i]);
                    }
                }
            }
        }
    }

    return included;
}
/******************************************************************************/
vector<vector<string> > Dataset::getSeqsBySample(vector<string> group){
    vector<vector<string> > result;
    // TODO
    return result;
}
/******************************************************************************/
void Dataset::reinstateSeqs(vector<string> trashTags){
    // TODO
}
/******************************************************************************/
void Dataset::removeSeqs(vector<string> names, vector<string> trashTags){
    // TODO
}
/******************************************************************************/
// fasta summary data: starts, ends, lengths, ambigs, polymers, numns
vector<vector<int>> Dataset::getFastaReport(){
    vector<vector<int> > results;

    results.push_back(select(starts, tableSeqs));
    results.push_back(select(ends, tableSeqs));
    results.push_back(select(lengths, tableSeqs));
    results.push_back(select(ambigs, tableSeqs));
    results.push_back(select(polymers, tableSeqs));
    results.push_back(select(numns, tableSeqs));

    return results;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getFastaSummary() {
    // create Summary object
    Summary* summary = new Summary(processors);

    Rcpp::DataFrame results = summary->summarizeFasta(
        getFastaReport(), count->getSeqsAbunds(getIncludedNamesIndexes()));

    delete summary;

    return results;
}
/******************************************************************************/
vector<vector<double>> Dataset::getContigsReport(){
    vector<vector<double> > result;
    // TODO
    return result;
}
/******************************************************************************/
vector<vector<double>> Dataset::getAlignReport(){
    vector<vector<double> > result;
    // TODO
    return result;
}
/******************************************************************************/
void Dataset::mergeSeqs(vector<string> names, string reason, string group){
    // TODO
}
/******************************************************************************/
int Dataset::getAbund(string name, string group){
    int abund = 0;
    auto it = seqIndex.find(name);
    if (it != seqIndex.end()) {
        abund = count->getAbund(it->second, group);
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


