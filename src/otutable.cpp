#include "../inst/include/rdataset.h"
#include "dataset.h"


/******************************************************************************/
OtuTable::OtuTable(AbundTable* c, map<string, int>* sI, string l) {
    hasSeqIds = false;
    numOtus = 0;
    numUnique = -1;
    uniqueBad = 0;
    label = l;
    otuCount = new AbundTable();
    seqCount = c;
    seqIndex = sI;
}
/******************************************************************************/
OtuTable::OtuTable(string l) {
    hasSeqIds = false;
    numOtus = 0;
    numUnique = -1;
    uniqueBad = 0;
    label = l;
    otuCount = new AbundTable();
    seqCount = nullptr;
}
/******************************************************************************/
OtuTable::~OtuTable() { delete otuCount; }
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

    otuCount->clear();
}
/******************************************************************************/
void OtuTable::add(vector<string> otuIds, vector<int> abundance,
         vector<string> samples, vector<string> seqIds){

    if (abundance.empty() && seqIds.empty()) {
        // dataset should catch this
        return;
    }

    bool useSeqIds = false;
    if (seqIds.size() != 0) {
        useSeqIds = true;
        hasSeqIds = true;
    }

    vector<string> uniqueOtuIds;
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
            uniqueOtuIds.push_back(otuName);
        }
    }

    // if the user provides seqIds, we need to calculate the otu abunds by sample
    // This can be done from the abundance and sample provided as parameters
    // or if they are not provided, abundance and sample data can be retrieved
    // from the seqCount variable from 'dataset'.
    if (useSeqIds) {

        bool hasAbundance = false;
        if (abundance.size() != 0) {
            hasAbundance = true;
        }

        bool hasSamples = false;
        if (samples.size() != 0) {
            hasSamples = true;
        }

        // otuIndex -> (sampleName -> abundance)
        map<int, map<string, int>> otuAbunds;

        vector<string> allSamples;
        bool seqCountSamples = false;
        if (!hasAbundance) {
            allSamples = seqCount->getSamples();
            if (!allSamples.empty()) {
                seqCountSamples = true;
            }
        }

        // maps seqName to otu index
        map<string, int> seqOtuIndex;
        numUnique = 0;

        sequenceOtus.resize(sequenceOtus.size()+uniqueOtuIds.size(), "");

        for (int i = 0; i < seqIds.size(); i++) {

            string otuName = otuIds[i];
            int index = otuIndex[otuName];

            auto itSeq = seqOtuIndex.find(seqIds[i]);

            // first time we are seeing this seqName, add to otu
            if (itSeq == seqOtuIndex.end()) {
                seqOtuIndex[seqIds[i]] = index;
                numUnique++;

                // add to sequenceOtus - (list file)
                if (sequenceOtus[index] != "") {
                    sequenceOtus[index] += "," + seqIds[i];
                }else{
                    sequenceOtus[index] = seqIds[i];
                }
            }

            auto it = otuAbunds.find(index);

            // first time seeing this otu, initialize sample counts
            if (it == otuAbunds.end()) {
                map<string, int> thisSeqsSampleAbunds;
                // inputs from parameters
                if (hasAbundance) {
                   if (hasSamples) {
                       thisSeqsSampleAbunds[samples[i]] = abundance[i];
                   }else{
                       thisSeqsSampleAbunds["total"] = abundance[i];
                   }
                }else{
                    // inputs from dataset
                    int thisIndex = (*seqIndex)[seqIds[i]];
                    vector<int> abunds = seqCount->getAbundances(thisIndex);

                    // seqCount has samples
                    if (seqCountSamples) {
                        for (int j = 0; j < abunds.size(); j++) {
                            if (abunds[j] != 0) {
                                thisSeqsSampleAbunds[allSamples[j]] = abunds[j];
                            }
                        }
                    }else{
                        thisSeqsSampleAbunds["total"] = abunds[0];
                    }
                }
                otuAbunds[index] = thisSeqsSampleAbunds;
            }else{
                // inputs from parameters
                if (hasAbundance) {
                    if (hasSamples) {
                        // have we seen this sample before?
                        auto itSample = it->second.find(samples[i]);
                        if (itSample == it->second.end()) {
                            it->second[samples[i]] = abundance[i];
                        }else{
                            itSample->second += abundance[i];
                        }
                    }else{
                        it->second["total"] += abundance[i];
                    }
                }else{
                    // inputs from dataset
                    int thisIndex = (*seqIndex)[seqIds[i]];
                    vector<int> abunds = seqCount->getAbundances(thisIndex);

                    // seqCount has samples
                    if (seqCountSamples) {
                        for (int j = 0; j < abunds.size(); j++) {
                            if (abunds[j] != 0) {
                                // have we seen this sample before?
                                auto itSample = it->second.find(allSamples[j]);
                                if (itSample == it->second.end()) {
                                    it->second[allSamples[j]] = abunds[j];
                                }else{
                                    itSample->second += abunds[j];
                                }
                            }
                        }
                    }else{
                        it->second["total"] += abunds[0];
                    }
                }
            }
        }

        // abundances are parsed by samples
        if (hasSamples || seqCountSamples) {
            vector<int> otuIndexes;
            abundance.clear();
            samples.clear();

            for (auto it = otuAbunds.begin(); it != otuAbunds.end(); it++) {

                map<string, int> abunds = it->second;
                for (auto itSamples = abunds.begin(); itSamples != abunds.end();
                        itSamples++) {
                    otuIndexes.push_back(it->first);
                    samples.push_back(itSamples->first);
                    abundance.push_back(itSamples->second);
                }
            }

            otuCount->assignAbundance(otuIndexes, abundance, samples);

            // inputs from dataset
            if (seqCountSamples) {
                // if there are treatments add them to otu
                if (seqCount->getNumTreatments() != 0) {
                    map<string, string> sampleTreatments =
                        seqCount->getSampleTreatmentAssignments();

                    vector<string> samples(sampleTreatments.size());
                    vector<string> treatments(sampleTreatments.size());

                    int index = 0;
                    for (auto itTreatments = sampleTreatments.begin();
                         itTreatments != sampleTreatments.end(); itTreatments++) {
                        samples[index] = itTreatments->first;
                        treatments[index] = itTreatments->second;
                        index++;
                    }

                    otuCount->assignTreatments(samples, treatments);
                }
            }

        }else{
            vector<int> totalAbunds(otuAbunds.size(), 0);
            vector<int> otuIndexes(otuAbunds.size(), 0);

            int index = 0;
            for (auto it = otuAbunds.begin(); it != otuAbunds.end(); it++) {
                otuIndexes[index] = it->first;
                totalAbunds[index] = it->second["total"];
                index++;
            }
            otuCount->assignAbundance(otuIndexes, totalAbunds);
        }
    }else{
        otuCount->assignAbundance(getIndexes(otuIds), abundance, samples);
    }
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
// TODO - complete removeSequences function to allow for removing seqs from dataset and otutables

void OtuTable::removeSequences(vector<string> ids, vector<string> reasons) {
    if (hasSeqIds) {

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
                otuCount->setAbundance(index, abunds[i]);
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
                otuCount->setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    otuCount->updateTotals();
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

