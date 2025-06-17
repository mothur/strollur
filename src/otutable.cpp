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
    trashCodes.clear();

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
            trashCodes.push_back("");

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
    count->assignTreatments(samples, treatments);
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
vector<int> OtuTable::getIndexes(vector<string> ids) {

    if (ids.empty()) {
        ids = getOtuIds();
    }
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
    return count->getTotalAbundances(getIndexes());
}
/******************************************************************************/
// abundances for each OTU broken down by sample
vector<vector<int> > OtuTable::getShared(){

    vector<int> goodOtus = getIndexes();

    vector<vector<int> > otus(goodOtus.size());

    for (int i = 0; i < goodOtus.size(); i++) {
        otus[i] = count->getAbundances(goodOtus[i]);
    }

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
int OtuTable::getNumNames(string& names){
    if(names == ""){ return 0; }

    int count = 1;
    for_each(names.begin(), names.end(),[&count](char n){
        if(n == ','){ count++; }
    });

    return count;
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
// id, trashCode
Rcpp::DataFrame OtuTable::getScrapReport() {

    if (badAccnos.size() != 0) {
        vector<string> badNames(uniqueBad, "");
        vector<string> badCodes(uniqueBad, "");

        int next = 0;
        for (int i = 0; i < trashCodes.size(); i++) {
            if (trashCodes[i] != "") {
                badNames[next] = otuNames[i];
                // remove last comma
                trashCodes[i].pop_back();
                badCodes[next] = trashCodes[i];
                next++;
            }
        }

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("otu_id") = badNames,
            Rcpp::_["trash_code"] = badCodes);
        return df;
    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// trashCode, otuCount, abundanceCount
Rcpp::DataFrame OtuTable::getScrapSummary() {

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
            Rcpp::_["otu_count"] = uniqueCounts,
            Rcpp::_["total_abundance"] = totalCounts);

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
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
    if (otuIDS.size() != 1) {

        vector<int> indexes = getIndexes(otuIDS);
        count->merge(indexes);

        for (int i = 1; i < otuIDS.size(); i++) {
            // no need to update the sample and treatment counts
            remove(otuIDS[i], reason, false);
        }
    }
}
/******************************************************************************/
// remove given outID
void OtuTable::remove(string otuID, string reason, bool update){

    auto it = otuIndex.find(otuID);

    // otu is not in dataset
    if (it == otuIndex.end()) {
        string message = "[WARNING]: " + otuID + " is not in ";
        message += "your dataset, ignoring.";
        RcppThread::Rcout << endl << message << endl;
    }else{

        int index = it->second;

        // remove from tableOtus and add trashCode
        tableOtus[index] = false;
        trashCodes[index] += reason;

        numOtus--;
        if (hasSeqIds) {
            string seqsInOtu = sequenceOtus[index];
            numUnique -= getNumNames(seqsInOtu);
        }

        // remove from counts
        int abund = 1;
        if (update) {
            abund = count->remove(index);
        }else{
            abund = count->getAbundance(index);
        }

        auto itBad = badAccnos.find(reason);

        if (itBad != badAccnos.end()) {
            // update counts of trashCode
            itBad->second[0]++;
        }else{
            // add new trashCode
            vector<int> badAbunds(2, 1);
            badAbunds[1] = abund;
            badAccnos[reason] = badAbunds;
        }

        // update uniqueBad
        uniqueBad++;
    }
}
/******************************************************************************/
// for datasets without samples
void OtuTable::setAbundance(vector<string> n,
                            vector<int> abunds, string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. When setting an otu ";
        message += "abundance, you must provide an abundance for each otu name.";
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";

    for (int i = 0; i < n.size(); i++) {
        auto it = otuIndex.find(n[i]);

        if (it != otuIndex.end()) {
            int index = it->second;

            if (abunds[i] == 0) {
                remove(n[i], reason);
            }else{
                count->setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }
}
/******************************************************************************/
// for datasets with samples
void OtuTable::setAbundances(vector<string> n, vector<vector<int>> abunds,
                   string reason){
    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch.  When setting an otu ";
        message += "abundances, you must provide an abundance vector ";
        message += "for each otu name.";
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";

    for (int i = 0; i < n.size(); i++) {
        auto it = otuIndex.find(n[i]);

        if (it != otuIndex.end()) {
            int index = it->second;

            if (sum(abunds[i]) == 0) {
                remove(n[i], reason);
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
}
/******************************************************************************/
void OtuTable::setLabel(string l){
    label = l;
}
/******************************************************************************/
void OtuTable::updateTotals() {
    count->updateTotals();
}
/******************************************************************************/

