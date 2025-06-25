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
Rcpp::DataFrame OtuTable::getList() {
    if (hasSeqIds) {
        vector<string> ids, seqids;
        int total = 0;

        for (int i = 0; i < tableOtus.size(); i++) {
            if (tableOtus[i]) {
                //parse string by delim and store in vector
                total += split(sequenceOtus[i], ',',
                                    back_inserter(seqids));
                ids.resize(total, otuNames[i]);
            }
        }

        string tag = label + "_id";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(tag) = ids,
            Rcpp::_["seq_id"] = seqids);

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// vector string containing sequence names for each otu
vector<string> OtuTable::getListVector(){
    return select(sequenceOtus, tableOtus);
}
/******************************************************************************/
Rcpp::DataFrame OtuTable::getRAbund() {
    if (numOtus != 0) {
        string tag = label + "_id";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(tag) = select(otuNames, tableOtus),
            Rcpp::_["abundance"] = count->getTotalAbundances(getGoodIndexes())
        );

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// vector of total abundances for each outID
vector<int> OtuTable::getRAbundVector(){
    return count->getTotalAbundances(getGoodIndexes());
}
/******************************************************************************/
Rcpp::DataFrame OtuTable::getShared() {
    if ((numOtus != 0) && (getNumSamples() != 0)) {

        vector<string> ids = getOtuIds();
        vector<int> indexes(ids.size(), -1);

        for (int i = 0; i < ids.size(); i++) {
            indexes[i] = otuIndex[ids[i]];
        }
        Rcpp::DataFrame df = count->getAbundanceTable(ids, indexes, false);
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
        otus[i] = count->getAbundances(goodOtus[i]);
    }

    return otus;
}
/******************************************************************************/
// string containing sequence names for given otuID
string OtuTable::get(string otuId){

    if (hasSeqIds) {
        auto it = otuIndex.find(otuId);
        if (it != otuIndex.end()) {
            if (tableOtus[it->second]) {
                return sequenceOtus[it->second];
            }
        }
    }
    return "";
}
/******************************************************************************/
// total abundance for a given outID, optional sample
int OtuTable::getAbundance(string otuId, string sample){

    auto it = otuIndex.find(otuId);
    if (it != otuIndex.end()) {
        // if this is a "good" otu
        if (tableOtus[it->second]) {
            return count->getAbundance(it->second, sample);
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
            return count->getAbundances(it->second);
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
    return count->getTotal(sample);
}
/******************************************************************************/
void OtuTable::merge(vector<string> otuIDS, string reason){
    if (otuIDS.size() != 1) {

        vector<int> indexes = getIndexes(otuIDS);

        count->merge(indexes);

        for (int i = 1; i < indexes.size(); i++) {
            if (hasSeqIds) {
                sequenceOtus[indexes[0]] += "," + sequenceOtus[indexes[i]];
                sequenceOtus[indexes[i]] = "";
            }

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

