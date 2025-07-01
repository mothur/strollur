#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
OtuTable::OtuTable(string l) {
    numOtus = 0;
    uniqueBad = 0;
    label = l;
    otuCount = new AbundTable();
}
/******************************************************************************/
OtuTable::~OtuTable() { delete otuCount; }
/******************************************************************************/
void OtuTable::clear() {
    numOtus = 0;
    uniqueBad = 0;
    label = "";

    otuIndex.clear();
    tableOtus.clear();
    otuNames.clear();
    badAccnos.clear();
    trashCodes.clear();

    otuCount->clear();
}
/******************************************************************************/
void OtuTable::assignAbundance(vector<string> otuIds, vector<int> abundance,
         vector<string> samples){

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
        }
    }

    otuCount->assignAbundance(getIndexes(otuIds), abundance, samples);
}
/******************************************************************************/
void OtuTable::assignTreatments(vector<string> samples,
                      vector<string> treatments) {
    otuCount->assignTreatments(samples, treatments);
}

/******************************************************************************/
vector<int> OtuTable::getGoodIndexes() {

    vector<string> ids = getOtuIds();

    vector<int> indexes(ids.size(), -1);

    for (int i = 0; i < ids.size(); i++) {
        indexes[i] = otuIndex[ids[i]];
    }
    return indexes;
}
/******************************************************************************/
// names of OTUs
vector<string> OtuTable::getOtuIds(){
    return select(otuNames, tableOtus);
}
/******************************************************************************/
Rcpp::DataFrame OtuTable::getRAbund() {
    if (numOtus != 0) {
        string tag = label + "_id";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(tag) = select(otuNames, tableOtus),
            Rcpp::_["abundance"] = otuCount->getTotalAbundances(getGoodIndexes())
        );

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// vector of total abundances for each outID
vector<int> OtuTable::getRAbundVector(){
    return otuCount->getTotalAbundances(getGoodIndexes());
}
/******************************************************************************/
Rcpp::DataFrame OtuTable::getShared() {
    if ((numOtus != 0) && (getNumSamples() != 0)) {

        vector<string> ids = getOtuIds();
        vector<int> indexes(ids.size(), -1);

        for (int i = 0; i < ids.size(); i++) {
            indexes[i] = otuIndex[ids[i]];
        }
        Rcpp::DataFrame df = otuCount->getAbundanceTable(ids, indexes, false);
        vector<string> names(3);
        names[0] = label + "_id";
        names[1] = "abundance";
        names[2] = "sample";
        df.attr("names") = names;

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// abundances for each OTU broken down by sample
vector<vector<int> > OtuTable::getSharedVector(){

    vector<int> goodOtus = getGoodIndexes();

    vector<vector<int> > otus(goodOtus.size());

    for (int i = 0; i < goodOtus.size(); i++) {
        otus[i] = otuCount->getAbundances(goodOtus[i]);
    }

    return otus;
}
/******************************************************************************/
// total abundance for a given outID, optional sample
int OtuTable::getAbundance(string otuId, string sample){

    auto it = otuIndex.find(otuId);
    if (it != otuIndex.end()) {
        // if this is a "good" otu
        if (tableOtus[it->second]) {
            return otuCount->getAbundance(it->second, sample);
        }
    }

    return 0;
}
/******************************************************************************/
// abundances for given otuID broken down by sample
vector<int> OtuTable::getAbundances(string otuId){

    auto it = otuIndex.find(otuId);
    if (it != otuIndex.end()) {
        // if this is a "good" otu
        if (tableOtus[it->second]) {
            return otuCount->getAbundances(it->second);
        }
    }

    return nullIntVector;
}
/******************************************************************************/
vector<int> OtuTable::getIndexes(vector<string>& otuIDS) {

    bool done = false;
    int firstGoodIndex = 0;
    map<string, int>::iterator it;
    vector<string> validOtuIds;

    while (!done) {
        it = otuIndex.find(otuIDS[firstGoodIndex]);
        firstGoodIndex++;

        // done if we find otu
        if (it != otuIndex.end()) {
            done = true;
        }else if (firstGoodIndex >= otuIDS.size()){
            // done if out of otus to check
            done = true;
        }
    }

    vector<int> indexes;

    // you have a valid otuId
    if (it != otuIndex.end()) {
        indexes.push_back(it->second);
        validOtuIds.push_back(otuIDS[firstGoodIndex]);

        for (int i = firstGoodIndex; i < otuIDS.size(); i++) {
            it = otuIndex.find(otuIDS[i]);

            if (it != otuIndex.end()) {
                indexes.push_back(it->second);
                validOtuIds.push_back(otuIDS[i]);
            }
        }
    }

    otuIDS = validOtuIds;
    return indexes;
}
/******************************************************************************/
// sample functions
int OtuTable::getNumSamples(){
    return otuCount->getNumSamples();
}
/******************************************************************************/
int OtuTable::getNumTreatments(){
    return otuCount->getNumTreatments();
}
/******************************************************************************/
vector<string> OtuTable::getSamples(){
    return otuCount->getSamples();
}
/******************************************************************************/
vector<int> OtuTable::getSampleTotals(){
    return otuCount->getSampleTotals();
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
                if (trashCodes[i][trashCodes[i].length()-1] == ',') {
                    // remove last comma
                    trashCodes[i].pop_back();
                }
                badCodes[next] = trashCodes[i];
                next++;
            }
        }
        string tag = label + "_id";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(tag) = badNames,
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
        string tag = label + "_count";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("trash_code") = codes,
            Rcpp::_[tag] = uniqueCounts,
            Rcpp::_["total_abundance"] = totalCounts);

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
bool OtuTable::hasSample(string sample){
    return otuCount->hasSample(sample);
}
/******************************************************************************/
bool OtuTable::hasId(string otuId) {

    auto it = otuIndex.find(otuId);

    // otu is in dataset
    if (it != otuIndex.end()) {
        // otu is "good"
        if (tableOtus[it->second]) {
            return true;
        }
    }
    return false;
}
/******************************************************************************/
vector<string> OtuTable::getTreatments(){
    return otuCount->getTreatments();
}
/******************************************************************************/
// vector containing total abundance for each treatment
vector<int> OtuTable::getTreatmentTotals(){
    return otuCount->getTreatmentTotals();
}
/******************************************************************************/
// total number of sequences
int OtuTable::getTotal(string sample){
    return otuCount->getTotal(sample);
}
/******************************************************************************/
void OtuTable::merge(vector<string> otuIDS, string reason){
    if (otuIDS.size() != 1) {

        vector<int> indexes = getIndexes(otuIDS);

        otuCount->merge(indexes);

        for (int i = 1; i < indexes.size(); i++) {

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

        // remove from counts
        int abund = 1;
        if (update) {
            abund = otuCount->remove(index);
        }else{
            abund = otuCount->getAbundance(index);
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
vector<string> OtuTable::setAbundance(vector<string> n,
                            vector<int> abunds, string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. When setting an otu ";
        message += "abundance, you must provide an abundance for each otu name.";
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";
    vector<string> otusRemoved;

    for (int i = 0; i < n.size(); i++) {
        auto it = otuIndex.find(n[i]);

        if (it != otuIndex.end()) {
            int index = it->second;

            if (abunds[i] == 0) {
                remove(n[i], reason);
                otusRemoved.push_back(n[i]);
            }else{
                otuCount->setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }
    return otusRemoved;
}
/******************************************************************************/
// for datasets with samples
vector<string> OtuTable::setAbundances(vector<string> n, vector<vector<int>> abunds,
                   string reason){
    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch.  When setting an otu ";
        message += "abundances, you must provide an abundance vector ";
        message += "for each otu name.";
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";
    vector<string> otusRemoved;

    for (int i = 0; i < n.size(); i++) {
        auto it = otuIndex.find(n[i]);

        if (it != otuIndex.end()) {
            int index = it->second;

            if (sum(abunds[i]) == 0) {
                remove(n[i], reason);
                otusRemoved.push_back(n[i]);
            }else{
                otuCount->setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    otuCount->updateTotals();

    return otusRemoved;
}
/******************************************************************************/
void OtuTable::setLabel(string l){
    label = l;
}
/******************************************************************************/
void OtuTable::updateTotals() {
    otuCount->updateTotals();
}
/******************************************************************************/

