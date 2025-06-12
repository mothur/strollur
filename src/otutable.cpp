#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
OtuTable::OtuTable(string l) {
    numOtus = 0;
    numUnique = -1;
    uniqueBad = 0;
    label = l;
    count = new AbundTable();
}
/******************************************************************************/
OtuTable::~OtuTable() { delete count; }
/******************************************************************************/
void OtuTable::clear() {
    numOtus = 0;
    numUnique = -1;
    uniqueBad = 0;
    label = "";

    otuIndex.clear();
    seqOtuIndex.clear();
    tableOtus.clear();
    sequenceOtus.clear();
    otuNames.clear();
    badAccnos.clear();

    count->clear();
}
/******************************************************************************/
void OtuTable::add(vector<string> otuIds, vector<int> abundance,
         vector<string> samples, vector<string> seqIds){

    bool hasSeqIds = false;
    if (seqIds.size() != 0) {
        hasSeqIds = true;
    }

    // when we have seq ids, we need to merge the otus and abundances
    vector<string> uniqueOtuIds;
    map<string, int> mergedAbunds;

    for (int i = 0; i < otuIds.size(); i++) {

        string otuName = otuIds[i];

        auto it = otuIndex.find(otuName);

        // new otu
        if (it == otuIndex.end()) {
            otuIndex[otuName] = numOtus;
            numOtus++;
            // all otus are initially assumed "good"
            tableOtus.push_back(true);
            otuNames.push_back(otuName);

            if (hasSeqIds) {
                uniqueOtuIds.push_back(otuName);
                mergedAbunds[otuName] = abundance[i];
            }
        }else if (hasSeqIds) {
            // existing otu
            mergedAbunds[otuName] += abundance[i];
        }
    }

    // map seqid to otuid
    if (hasSeqIds) {
        seqOtuIndex.clear();
        numUnique = 0;

        for (int i = 0; i < otuIds.size(); i++) {

            string otuName = otuIds[i];

            auto itSeq = seqOtuIndex.find(seqIds[i]);

            if (itSeq == seqOtuIndex.end()) {
                seqOtuIndex[seqIds[i]] = otuIndex[otuName];
                numUnique++;
            }else{
                if (otuNames[otuIndex[otuName]] != otuName) {
                    string message = "[WARNING]: You previously assigned ";
                    message += seqIds[i] + " to " + otuNames[otuIndex[otuName]];
                    message += ", ignoring assignment to " + otuName + ".";
                    RcppThread::Rcout << endl << message << endl;
                }
            }
        }
        sort(uniqueOtuIds.begin(), uniqueOtuIds.end());
        abundance = getValues(mergedAbunds);

        count->assignAbundance(getIndexes(uniqueOtuIds), abundance);
    }else{
        count->assignAbundance(getIndexes(otuIds), abundance, samples);
    }
}
/******************************************************************************/
void OtuTable::assignTreatments(vector<string> samples,
                      vector<string> treatments) {
    //TODO
}
/******************************************************************************/
vector<int> OtuTable::getIndexes(vector<string>& ids) {
    vector<int> indexes(ids.size(), -1);
    map<string, int>::iterator it;

    for (int i = 0; i < ids.size(); i++) {

        it = otuIndex.find(ids[i]);
        if (it != otuIndex.end()) {
            indexes[i] = it->second;
        }else{
            string message = "[ERROR]: The dataset does not contain an ";
            message += "otu named " + ids[i] + ".\n";
            throw Rcpp::exception(message.c_str());
        }
    }
    return indexes;
}
/******************************************************************************/
// names of OTUs
vector<string> OtuTable::getOtuIds(){
    return select(otuNames, tableOtus);
}
/******************************************************************************/
// vector string containing sequence names for each otu
vector<string> OtuTable::getList(){
    return select(sequenceOtus, tableOtus);
}
/******************************************************************************/
// vector of total abundances for each outID
vector<int> OtuTable::getRAbund(){
    vector<int> otus;
    //TODO
    return otus;
}
/******************************************************************************/
// abundances for each OTU broken down by sample
vector<vector<int> > OtuTable::getShared(){
    vector<vector<int> > otus;
    //TODO
    return otus;
}
/******************************************************************************/
// string containing sequence names for given otuID
string OtuTable::get(string otuID, string sample){
    string otu = "";
    //TODO
    return otu;
}
/******************************************************************************/
// total abundance for a given outID, optional sample
int OtuTable::getAbundance(string otuID, string sample){
    int abund = -1;
    //TODO
    return abund;
}
/******************************************************************************/
// abundances for given otuID broken down by sample
vector<int> OtuTable::getAbundances(string otuID){
    vector<int> abunds;
    //TODO
    return abunds;
}
/******************************************************************************/
// sample functions
int OtuTable::getNumSamples(){
    return count->getNumSamples();
}
/******************************************************************************/
int OtuTable::getNumTreatments(){
    return count->getNumTreatments();
}
/******************************************************************************/
vector<string> OtuTable::getSamples(){
    return count->getSamples();
}
/******************************************************************************/
vector<int> OtuTable::getSampleTotals(){
    return count->getSampleTotals();
}
/******************************************************************************/
bool OtuTable::hasSample(string sample){
    return count->hasSample(sample);
}
/******************************************************************************/
vector<string> OtuTable::getTreatments(){
    return count->getTreatments();
}
/******************************************************************************/
// vector containing total abundance for each treatment
vector<int> OtuTable::getTreatmentTotals(){
    return count->getTreatmentTotals();
}
/******************************************************************************/
// total number of sequences
int OtuTable::getTotal(string sample){
    return count->getTotal();
}
/******************************************************************************/
void OtuTable::merge(vector<string> otuIDS, string reason){
    //TODO
}
/******************************************************************************/
// remove given outID
void OtuTable::remove(string otuID){
    //TODO
}
/******************************************************************************/
// for datasets without samples
void OtuTable::setAbundance(vector<string> otuIDS, vector<int> abunds,
                  string reason){
    //TODO
}
/******************************************************************************/
// for datasets with samples
void OtuTable::setAbundances(vector<string> otuIDS, vector<vector<int>> abunds,
                   string reason){
    //TODO
}
/******************************************************************************/
void OtuTable::setLabel(string l){
    label = l;
}
/******************************************************************************/

