#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
OtuTable::OtuTable(string l) {
    hasSeqIds = false;
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
    hasSeqIds = false;

    otuIndex.clear();
    tableOtus.clear();
    sequenceOtus.clear();
    otuNames.clear();
    badAccnos.clear();

    count->clear();
}
/******************************************************************************/
void OtuTable::add(vector<string> otuIds, vector<int> abundance,
         vector<string> samples, vector<string> seqIds){

    bool useSeqIds = false;
    if (seqIds.size() != 0) {
        useSeqIds = true;
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

            if (useSeqIds) {
                uniqueOtuIds.push_back(otuName);
                mergedAbunds[otuName] = abundance[i];
            }
        }else if (useSeqIds) {
            // existing otu
            mergedAbunds[otuName] += abundance[i];
        }
    }

    // map seqid to otuid
    if (useSeqIds) {
        // maps seqName to otu index
        map<string, int> seqOtuIndex;
        numUnique = 0;

        sequenceOtus.resize(sequenceOtus.size()+uniqueOtuIds.size(), "");

        for (int i = 0; i < seqIds.size(); i++) {

            string otuName = otuIds[i];
            int index = otuIndex[otuName];

            auto itSeq = seqOtuIndex.find(seqIds[i]);

            if (itSeq == seqOtuIndex.end()) {
                seqOtuIndex[seqIds[i]] = index;
                numUnique++;
            }else{
                // sanity check, each sequence name must be unique
                if (otuNames[index] != otuName) {
                    string message = "[WARNING]: You previously assigned ";
                    message += seqIds[i] + " to " + otuNames[index];
                    message += ", ignoring assignment to " + otuName + ".";
                    RcppThread::Rcout << endl << message << endl;
                }
            }

            // add to sequenceOtus
            if (sequenceOtus[index] != "") {
                sequenceOtus[index] += "," + seqIds[i];
            }else{
                sequenceOtus[index] = seqIds[i];
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
int OtuTable::getIndex(string& id) {
    auto it = otuIndex.find(id);
    if (it != otuIndex.end()) {
        return it->second;
    }else{
        string message = "[ERROR]: The dataset does not contain an ";
        message += "otu named " + id + ".\n";
        throw Rcpp::exception(message.c_str());
    }
}
/******************************************************************************/
vector<int> OtuTable::getIndexes(vector<string>& ids) {
    vector<int> indexes(ids.size(), -1);

    for (int i = 0; i < ids.size(); i++) {
        indexes[i] = getIndex(ids[i]);
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
string OtuTable::get(string otuId){

    if (hasSeqIds) {
        int otuIndex = getIndex(otuId);

        // if this is a "good" otu
        if (tableOtus[otuIndex]) {
            return sequenceOtus[otuIndex];
        }
    }
    return "";
}
/******************************************************************************/
// total abundance for a given outID, optional sample
int OtuTable::getAbundance(string otuId, string sample){
    int otuIndex = getIndex(otuId);

    // if this is a "good" otu
    if (tableOtus[otuIndex]) {
        return count->getAbundance(otuIndex, sample);
    }
    return 0;
}
/******************************************************************************/
// abundances for given otuID broken down by sample
vector<int> OtuTable::getAbundances(string otuId){
    int otuIndex = getIndex(otuId);

    // if this is a "good" otu
    if (tableOtus[otuIndex]) {
        return count->getAbundances(otuIndex);
    }
    return nullIntVector;
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

