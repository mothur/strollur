
#include "../inst/include/rdataset.h"

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

   // TODO add calcs for starts, ends, lengths, ambigs, polymers, numns



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
    vector<string> result;
    // TODO
    return result;
}
/******************************************************************************/
vector<int> Dataset::getGroupTotals(string name){
    vector<int> result;
    // TODO
    return result;
}
/******************************************************************************/
long long Dataset::getTotal(string group){
    long long result = 0;
    // TODO
    return result;
}
/******************************************************************************/
bool Dataset::hasGroup(string group){
    bool result = false;
    // TODO
    return result;
}
/******************************************************************************/
vector<string> Dataset::getNames(string group){
    //vector<string> result;
    // TODO
   // return result;
   return names;
}
/******************************************************************************/
vector<vector<string> > Dataset::getNamesBySample(vector<string> group){
    vector<vector<string> > result;
    // TODO
    return result;
}
/******************************************************************************/
vector<string> Dataset::getSeqs(string group){
    // vector<string> result;
    // // TODO
    // return result;
    return seqs;
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
vector<vector<int>> Dataset::getFastaReport(){
    vector<vector<int> > result;
    // TODO
    return result;
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

void Dataset::assignSampleAbundance(vector<string>,
                                    vector<string>,
                                    vector<int>){
    // TODO
}
/******************************************************************************/
int Dataset::getAbund(string name, string group){
    int result = 0;
    // TODO
    return result;
}
/******************************************************************************/
// abundances for seq broken down by sample
vector<int> Dataset::getAbunds(string name){
    vector<int>  result;
    // TODO
    return result;
}
/******************************************************************************/
// total abundance for each sequence
vector<int> Dataset::getSeqsAbunds(){
    vector<int>  result;
    // TODO
    return result;
}
/******************************************************************************/
vector<vector<int>> Dataset::getSeqsAbundsBySample(){
    vector<vector<int> > result;
    // TODO
    return result;
}

/******************************************************************************/
//vector<double> searchScores,
//vector<double> simScores,
//vector<double> longestInserts
// align_seqs will set searchScores, simScores and longestInserts
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


