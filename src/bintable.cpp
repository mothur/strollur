#include "../inst/include/rdataset.h"
#include "dataset.h"
#include "phylotree.h"

/******************************************************************************/
BinTable::BinTable() {
    uniqueBad = 0;
    label = "";
    hasListAssignments = false;
    hasBinTaxonomy = false;
    runClassify = false;
}
/******************************************************************************/
BinTable::BinTable(string l) {
    uniqueBad = 0;
    label = l;
    hasListAssignments = false;
    hasBinTaxonomy = false;
    runClassify = false;
}
/******************************************************************************/
BinTable::BinTable(const BinTable& binTable) {
    label = binTable.label;
    hasListAssignments = binTable.hasListAssignments;
    hasBinTaxonomy = binTable.hasBinTaxonomy;

    binIndex = binTable.binIndex;
    tableBins = binTable.tableBins;
    binNames = binTable.binNames;
    trashCodes = binTable.trashCodes;
    taxonomies = binTable.taxonomies;
    runClassify = binTable.runClassify;

    binList = binTable.binList;
    seqBins = binTable.seqBins;

    badAccnos = binTable.badAccnos;
    uniqueBad = binTable.uniqueBad;

    binCount.clone(binTable.binCount);
}
/******************************************************************************/
BinTable::~BinTable() {}
/******************************************************************************/
void BinTable::clear(string tag) {

    // clear all
    if (tag == "") {
        uniqueBad = 0;
        label = "";
        hasListAssignments = false;
        hasBinTaxonomy = false;
        runClassify = false;

        binIndex.clear();
        tableBins.clear();
        binNames.clear();
        badAccnos.clear();
        trashCodes.clear();
        seqBins.clear();
        binList.clear();
        taxonomies.clear();
        binCount.clear();

    }else if (tag == "taxonomy") {
        hasBinTaxonomy = false;
        runClassify = false;
        taxonomies.clear();
    }

}
/******************************************************************************/
void BinTable::clone(const BinTable& binTable) {
    label = binTable.label;
    hasListAssignments = binTable.hasListAssignments;
    hasBinTaxonomy = binTable.hasBinTaxonomy;

    binIndex = binTable.binIndex;
    tableBins = binTable.tableBins;
    binNames = binTable.binNames;
    trashCodes = binTable.trashCodes;
    taxonomies = binTable.taxonomies;
    runClassify = binTable.runClassify;

    binList = binTable.binList;
    seqBins = binTable.seqBins;

    badAccnos = binTable.badAccnos;
    uniqueBad = binTable.uniqueBad;

    binCount.clone(binTable.binCount);
}
/******************************************************************************/
double BinTable::assignAbundance(vector<string> binIds,
                                 const vector<float>& abundance,
                                 const vector<string>& samples,
                                 const vector<int>& seqIds,
                                 AbundTable& count, bool update){

    double otusAssigned = 0;

    if (update) {
        // seqs have been assigned to bins, we are using count to update
        // the bin abundance. This is a case where (mothur1) list file was
        // added, then sequence abundances are assigned. We want to update the
        // bin abundance with the new seq abundances.

        // we need to calculate the bin abunds by sample
        vector<string> allSamples = count.getSamples();

        // does the sequence abundance include samples
        bool hasSamples = false;
        if (!allSamples.empty()) {
            hasSamples = true;
        }

        // binName -> (sampleName -> abundance)
        map<int, map<string, float>> binAbunds;

        // seq2 -> bin1
        for (auto it = seqBins.begin(); it != seqBins.end(); it++) {
            updateBinAbunds(binAbunds, count, it->second, it->first, allSamples);
        }

        otusAssigned = updateBins(binAbunds, count, hasSamples);
        runClassify = true;
   }else{

        // just abundances
        if (seqIds.empty()) {

            int numBins = tableBins.size()-1;
            if (numBins == -1) { numBins = 0; }

            for (int i = 0; i < binIds.size(); i++) {

                string binName = binIds[i];

                auto it = binIndex.find(binName);

                // new bin
                if (it == binIndex.end()) {
                    binIndex[binName] = numBins;
                    numBins++;
                    // all bins are initially assumed "good"
                    tableBins.push_back(true);
                    binNames.push_back(binName);
                    trashCodes.push_back("");
                }
            }

            vector<int> binIndexes = getIndexes(binIds);
            otusAssigned = binCount.assignAbundance(binIndexes,
                                                    abundance, samples);

        }else {

            // we are assigning seqIds to bins, so remove all old bin data
            // to avoid inconsistencies
            string oldLabel = label;
            clear();
            int numBins = 0;
            label = oldLabel;

            // we need to calculate the bin abunds by sample
            vector<string> allSamples = count.getSamples();

            // does the sequence abundance include samples
            bool hasSamples = false;
            if (!allSamples.empty()) {
                hasSamples = true;
            }

            // binName -> (sampleName -> abundance)
            map<int, map<string, float>> binAbunds;

            for (int i = 0; i < seqIds.size(); i++) {

                auto itIndex = binIndex.find(binIds[i]);
                int index = numBins;
                int seqIndex = seqIds[i];

                // new bin
                if (itIndex == binIndex.end()) {
                    binIndex[binIds[i]] = numBins;
                    numBins++;
                    // all bins are initially assumed "good"
                    tableBins.push_back(true);
                    binNames.push_back(binIds[i]);
                    trashCodes.push_back("");
                    set<int> temp;
                    binList.push_back(temp);
                }else{
                    index = itIndex->second;
                }

                auto itSeq = seqBins.find(seqIndex);
                bool firstTimeSeq = false;

                // first time we are seeing this seqName, add seq to bin
                if (itSeq == seqBins.end()) {
                    firstTimeSeq = true;
                    seqBins[seqIndex] = index;
                    binList[index].insert(seqIndex);
                }

                updateBinAbunds(binAbunds, count, index, seqIndex, allSamples,
                                firstTimeSeq);
            }

            otusAssigned = updateBins(binAbunds, count, hasSamples);

            hasListAssignments = true;
            runClassify = true;
        }
    }

    return otusAssigned;
}
/******************************************************************************/
double BinTable::assignTaxonomy(const vector<string>& binIds,
                                const vector<string>& taxs){

    double numBinsAssigned = 0;

    // make space for assignments
    taxonomies.resize(binNames.size(), "");

    for (int i = 0; i < binIds.size(); i++) {

        string binName = binIds[i];

        auto it = binIndex.find(binName);

        // valid bin
        if (it != binIndex.end()) {
            numBinsAssigned++;
            taxonomies[it->second] = taxs[i];
        }
    }
    runClassify = false;
    hasBinTaxonomy = true;
    return numBinsAssigned;
}
/******************************************************************************/
double BinTable::assignTreatments(const vector<string>& samples,
                                  const vector<string>& treatments) {
    return binCount.assignTreatments(samples, treatments);
}
/******************************************************************************/
void BinTable::classify(const vector<string>& taxs,
                        AbundTable& count) {
    taxonomies.resize(binNames.size(), "");
    for (int i = 0; i < binNames.size(); i++) {
        if (tableBins[i]) {
            taxonomies[i] = classifyBin(i, taxs, count);
        }
    }
    hasBinTaxonomy = true;
    runClassify = false;
}
/******************************************************************************/
string BinTable::classifyBin(const int binId,
                             const vector<string>& tax,
                             AbundTable& count){
    string binTax = "";

    if (hasListAssignments && (tax.size() != 0)) {

        set<int> binSeqs = binList[binId];

        PhyloTree phylo;
        int size = 0;

        // add all "good" seqs in bin to phylotree
        for (int index : binSeqs) {
            int seqAbund = count.getAbundance(index);
            size += seqAbund;
            phylo.addSeqToTree(tax[index], seqAbund);
        }

        TaxNode currentNode = phylo.getRoot();

        while (!isZero(currentNode.total)) {

            TaxNode bestChild;
            int bestChildSize = 0;

            //go through children
            for (auto itChild = currentNode.children.begin();
                 itChild != currentNode.children.end(); itChild++) {

                TaxNode temp = phylo.get(itChild->second);

                // select taxonomy with most seqs assigned to it
                if (temp.total > bestChildSize) {
                    bestChild = phylo.get(itChild->second);
                    bestChildSize = temp.total;
                }
            }

            //phylotree adds an extra unknown so we want to remove that
            if (bestChild.name == "unknown") { bestChildSize--; }

            int consensusConfidence = ceil((bestChildSize /
                                           (float) size) * 100);

            if (bestChild.name != "") {
                binTax += bestChild.name + "(" + toString(consensusConfidence) + ");";
            }
            //move down a level
            currentNode = bestChild;
        }
    }

    if (binTax == "") {  binTax = "unknown;";  }

    return binTax;
}
/******************************************************************************/
Rcpp::List BinTable::exportBinTable() {
    Rcpp::List results = Rcpp::List::create();
    vector<string> resultLabels;

    Rcpp::DataFrame binData = Rcpp::DataFrame::create();
    vector<string> binDataLabels;

    vector<int> binIndexes = getIndexes(binNames);

    // ids, bin_names, bin_abundance, bin_taxonomy, trash_codes, binTable
    binData.push_back(binIndexes);
    binDataLabels.push_back("bin_ids");
    binData.push_back(binNames);
    binDataLabels.push_back("bin_names");
    binData.push_back(binCount.getTotalAbundances(binIndexes));
    binDataLabels.push_back("abundances");

    if (!allBlank(taxonomies)) {
        binData.push_back(taxonomies);
        binDataLabels.push_back("taxonomies");
    }

    if (!allBlank(trashCodes)) {
        binData.push_back(trashCodes);
        binDataLabels.push_back("trash_codes");
    }

    binData.push_back(tableBins);
    binDataLabels.push_back("include_bin");
    binData.attr("names") = binDataLabels;

    results.push_back(binData);
    resultLabels.push_back("bin_data");

    // bin_id, seq_id
    if (hasListAssignments) {
        Rcpp::DataFrame binSeqAssignments = Rcpp::DataFrame::create(
            Rcpp::Named("bin_ids") = getValues(seqBins),
            Rcpp::_["sequence_ids"] = getKeys(seqBins));
        results.push_back(binSeqAssignments);
        resultLabels.push_back("sequence_bin_assignments");
    }else{
        // bin_id, abund, sample, treatment
        results.push_back(binCount.getAbundanceTable(binNames,
                                                     binIndexes,
                                                  "bin", false));
        resultLabels.push_back("bin_abundance_table");
    }

    results.attr("names") = resultLabels;

    return results;
}
/******************************************************************************/
const vector<int> BinTable::getGoodIndexes() {

    vector<string> ids = getIds();

    vector<int> indexes(ids.size(), -1);

    for (int i = 0; i < ids.size(); i++) {
        indexes[i] = binIndex[ids[i]];
    }
    return indexes;
}
/******************************************************************************/
// string containing seqs in bin, comma separated
const string BinTable::get(const string binName,
                           const vector<string>& seqNames) {
    string bin = "";

    // you have assigned sequences to bins
    if (hasListAssignments) {

        auto it = binIndex.find(binName);
        if (it != binIndex.end()) {
            // "good" bin
            if (tableBins[it->second]) {
                set<int> seqs = binList[it->second];

                for (int seq : seqs) {
                    if (bin != "") {
                        bin += "," + seqNames[seq];
                    }else{
                        bin = seqNames[seq];
                    }
                }
            }
        }
    }
    return bin;
}
/******************************************************************************/
const vector<string> BinTable::getListVector(const vector<string>& seqNames) {
    vector<string> listVector;

    // you have assigned sequences to bins
    if (hasListAssignments) {
        for (int i = 0; i < binNames.size(); i++) {
            if (tableBins[i]) {
                listVector.push_back(get(binNames[i], seqNames));
            }
        }
    }
    return listVector;
}
/******************************************************************************/
// 2 column dataframe - bin_id, seq_id
const Rcpp::DataFrame BinTable::getList(const vector<string>& seqNames){

    if (hasListAssignments) {

        vector<string> ids, seqids;
        for (int i = 0; i < binNames.size(); i++) {

            // if this is a "good" bin
            if (tableBins[i]) {
                set<int> seqs = binList[i];

                for (int seq : seqs) {
                    ids.push_back(binNames[i]);
                    seqids.push_back(seqNames[seq]);
                }
            }
        }

        string tag = label + "_id";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(tag.c_str()) = ids,
            Rcpp::_["seq_id"] = seqids);
        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// names of OTUs
const vector<string> BinTable::getIds(){
    return select(binNames, tableBins);
}
/******************************************************************************/
vector<string> BinTable::getTaxonomies(const vector<string>& tax,
                                       AbundTable& count) {
    if (runClassify) {
        classify(tax, count);
    }

    if (taxonomies.size() != 0) {
        return select(taxonomies, tableBins);
    }

    return nullVector;
}
/******************************************************************************/
const Rcpp::DataFrame BinTable::getRAbund() {
    if (getNumBins()!= 0) {
        string tag = label + "_id";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(tag) = select(binNames, tableBins),
            Rcpp::_["abundance"] = binCount.getTotalAbundances(getGoodIndexes())
        );

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// vector of total abundances for each outID
const vector<float> BinTable::getRAbundVector(){
    return binCount.getTotalAbundances(getGoodIndexes());
}
/******************************************************************************/
const Rcpp::DataFrame BinTable::getShared() {
    if ((getNumBins() != 0) && (getNumSamples() != 0)) {

        vector<string> ids = getIds();
        vector<int> indexes(ids.size(), -1);

        for (int i = 0; i < ids.size(); i++) {
            indexes[i] = binIndex[ids[i]];
        }
        Rcpp::DataFrame df = binCount.getAbundanceTable(ids, indexes, "bin");
        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// abundances for each OTU broken down by sample
const vector<vector<float> > BinTable::getSharedVector(){

    vector<int> goodBins = getGoodIndexes();

    vector<vector<float> > bins(goodBins.size());

    for (int i = 0; i < goodBins.size(); i++) {
        bins[i] = binCount.getAbundances(goodBins[i]);
    }

    return bins;
}
/******************************************************************************/
// total abundance for a given binID, optional sample
const float BinTable::getAbundance(const string binId, const string sample){

    auto it = binIndex.find(binId);
    if (it != binIndex.end()) {
        // if this is a "good" bin
        if (tableBins[it->second]) {
            return binCount.getAbundance(it->second, sample);
        }
    }

    return 0;
}
/******************************************************************************/
// abundances for given binID broken down by sample
const vector<float> BinTable::getAbundances(const string binId){

    auto it = binIndex.find(binId);
    if (it != binIndex.end()) {
        // if this is a "good" bin
        if (tableBins[it->second]) {
            return binCount.getAbundances(it->second);
        }
    }

    return nullFloatVector;
}
/******************************************************************************/
const int BinTable::getNumBins() {
    return accumulate(tableBins.begin(), tableBins.end(), 0);
}
/******************************************************************************/
const vector<int> BinTable::getIndexes(vector<string>& binIDS) {

    bool done = false;
    int firstGoodIndex = 0;
    map<string, int>::iterator it;
    vector<string> validOtuIds;

    while (!done) {
        it = binIndex.find(binIDS[firstGoodIndex]);

        // done if we find bin
        if (it != binIndex.end()) {
            done = true;
        }else if (firstGoodIndex >= binIDS.size()){
            // done if out of bins to check
            done = true;
        }else{
            firstGoodIndex++;
        }
    }

    vector<int> indexes;

    // you have a valid binId
    if (it != binIndex.end()) {

        for (int i = firstGoodIndex; i < binIDS.size(); i++) {
            it = binIndex.find(binIDS[i]);

            if (it != binIndex.end()) {
                indexes.push_back(it->second);
                validOtuIds.push_back(binIDS[i]);
            }
        }
    }

    binIDS = validOtuIds;
    return indexes;
}
/******************************************************************************/
// sample functions
const int BinTable::getNumSamples(){
    return binCount.getNumSamples();
}
/******************************************************************************/
const int BinTable::getNumTreatments(){
    return binCount.getNumTreatments();
}
/******************************************************************************/
const vector<string> BinTable::getSamples(){
    return binCount.getSamples();
}
/******************************************************************************/
const vector<double> BinTable::getSampleTotals(){
    return binCount.getSampleTotals();
}
/******************************************************************************/
const map<string, string> BinTable::getSampleTreatmentAssignments() {
    return binCount.getSampleTreatmentAssignments();
}
/******************************************************************************/
// id, trashCode
const Rcpp::DataFrame BinTable::getScrapReport() {

    if (badAccnos.size() != 0) {
        vector<string> badNames(uniqueBad, "");
        vector<string> badCodes(uniqueBad, "");

        int next = 0;
        for (int i = 0; i < trashCodes.size(); i++) {
            if (trashCodes[i] != "") {
                badNames[next] = binNames[i];
                if (trashCodes[i][trashCodes[i].length()-1] == ',') {
                    // remove last comma
                    trashCodes[i].pop_back();
                }
                badCodes[next] = trashCodes[i];
                next++;
            }
        }

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("id") = badNames,
            Rcpp::_["trash_code"] = badCodes);
        return df;
    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// trashCode, binCount, abundanceCount
const Rcpp::DataFrame BinTable::getScrapSummary() {

    if (badAccnos.size() != 0) {
        vector<string> codes(badAccnos.size(), "");
        vector<float> uniqueCounts(badAccnos.size(), 0);
        vector<float> totalCounts(badAccnos.size(), 0);

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
const bool BinTable::hasSample(const string sample){
    return binCount.hasSample(sample);
}
/******************************************************************************/
const vector<string> BinTable::getTreatments(){
    return binCount.getTreatments();
}
/******************************************************************************/
// vector containing total abundance for each treatment
const vector<double> BinTable::getTreatmentTotals(){
    return binCount.getTreatmentTotals();
}
/******************************************************************************/
// total number of sequences
const double BinTable::getTotal(const string sample){
    return binCount.getTotal(sample);
}
/******************************************************************************/
const bool BinTable::okToMerge(const vector<int>& seqIds) {

    // if seqs were assigned to otus
    if (hasListAssignments) {

        int binNumber = -1;
        auto it = seqBins.find(seqIds[0]);
        if (it != seqBins.end()) {
            binNumber = it->second;
        }

        for (int i = 1; i < seqIds.size(); i++) {
            it = seqBins.find(seqIds[i]);

            if (it != seqBins.end()) {
                if (it->second != binNumber) {
                    return false;
                }
            }
        }
    }
    return true;
}
/******************************************************************************/
void BinTable::merge(vector<string> binIDS, const string reason){
    if (binIDS.size() != 1) {

        vector<int> indexes = getIndexes(binIDS);
        binCount.merge(indexes);

        for (int i = 1; i < indexes.size(); i++) {

            // no need to update the sample and treatment counts
            vector<int> seqIds = remove(binIDS[i], reason, false);

            // move seqs to merged bin
            if (hasListAssignments) {

                // set each seqs bin assignment
                for (int seq : seqIds) {
                    seqBins[seq] = indexes[0];
                    binList[indexes[0]].insert(seq);
                }

                binList[indexes[i]].clear();
            }
        }
    }
}
/******************************************************************************/
void BinTable::remove(const int seqId, AbundTable& count,
                      const string reason, bool update) {

    // if seqs were assigned to otus
    if (hasListAssignments) {

        runClassify = true;

        // remove from otu
        auto it = seqBins.find(seqId);

        if (it != seqBins.end()) {

            int index = it->second;

            // does removing the seq remove a bin
            binList[index].erase(seqId);

            // remove bin, if empty
            if (binList[index].size() == 0) {
                // remove from tableBins and add trashCode
                tableBins[index] = false;
                trashCodes[index] += reason;

                // remove from counts
                float abund = 1;
                if (update) {
                    abund = binCount.remove(index);
                }else{
                    abund = binCount.getAbundance(index);
                }

                auto itBad = badAccnos.find(reason);

                if (itBad != badAccnos.end()) {
                    // update counts of trashCode
                    itBad->second[0]++;
                }else{
                    // add new trashCode
                    vector<double> badAbunds(2, 1);
                    badAbunds[1] = abund;
                    badAccnos[reason] = badAbunds;
                }

                // update uniqueBad
                uniqueBad++;

            }else{
                // update bin counts using count table
                if (update) {
                    if (binCount.hasSamples()) {

                        vector<float> seqAbunds = count.getAbundances(seqId);
                        vector<string> seqSamples = count.getSamples();

                        for (int i = 0; i < seqAbunds.size(); i++) {

                            if (!isZero(seqAbunds[i]) &&
                                binCount.hasSample(seqSamples[i])) {

                                float orig = binCount.getAbundance(index, seqSamples[i]);
                                binCount.setAbundance(index,
                                                       (orig-seqAbunds[i]),
                                                       seqSamples[i]);
                            }
                        }
                    }else{
                        float seqAbund = count.getAbundance(seqId);
                        float binAbund = binCount.getAbundance(index);
                        binCount.setAbundance(index, (binAbund-seqAbund));
                    }
                }
            }

            seqBins.erase(seqId);
        }
    }
}
/******************************************************************************/
void BinTable::removeSamples(const vector<string>& samples) {
   binCount.removeSamples(samples);

   // remove any bins only assigned to these samples
   for (int i = 0; i < binNames.size(); i++) {
       // included bin
       if (tableBins[i]) {
           if (isZero(sum(binCount.getAbundances(i)))) {
               remove(binNames[i], "removedSamples", true);
           }
       }
   }
}
/******************************************************************************/
// remove given outID, returns seqs removed
vector<int> BinTable::remove(string binID, string reason, bool update){

    auto it = binIndex.find(binID);

    // bin is not in dataset
    if (it == binIndex.end()) {
        string message = "[WARNING]: " + binID + " is not in ";
        message += "your dataset, ignoring.";
        RcppThread::Rcout << endl << message << endl;
        return nullIntVector;
    }

    int index = it->second;

    // remove from tableBins and add trashCode
    tableBins[index] = false;
    trashCodes[index] += reason;

    // remove from counts
    float abund = 1;
    if (update) {
        abund = binCount.remove(index);
    }else{
        abund = binCount.getAbundance(index);
    }

    auto itBad = badAccnos.find(reason);

    if (itBad != badAccnos.end()) {
        // update counts of trashCode
        itBad->second[0]++;
    }else{
        // add new trashCode
        vector<double> badAbunds(2, 1);
        badAbunds[1] = abund;
        badAccnos[reason] = badAbunds;
    }

    // update uniqueBad
    uniqueBad++;

    if (hasListAssignments) {
        return toVector(binList[index]);
    }

    return nullIntVector;
}
/******************************************************************************/
// for datasets without samples
vector<int> BinTable::setAbundance(const vector<string>& n,
                                   const vector<float>& abunds,
                                   string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. When setting an bin ";
        message += "abundance, you must provide an abundance for each bin name.";
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";
    vector<int> seqsRemoved;

    for (int i = 0; i < n.size(); i++) {
        auto it = binIndex.find(n[i]);

        if (it != binIndex.end()) {
            int index = it->second;

            if (isZero(abunds[i])) {
                vector<int> thisBinsSeqsToRemove = remove(n[i], reason);
                seqsRemoved.insert(seqsRemoved.end(),
                                   thisBinsSeqsToRemove.begin(),
                                   thisBinsSeqsToRemove.end());
            }else{
                binCount.setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }
    return seqsRemoved;
}

/******************************************************************************/
// for datasets with samples
vector<int> BinTable::setAbundances(const vector<string>& n,
                                    const vector<vector<float>>& abunds,
                                    string reason){
    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch.  When setting an bin ";
        message += "abundances, you must provide an abundance vector ";
        message += "for each bin name.";
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";
    vector<int> seqsRemoved;

    for (int i = 0; i < n.size(); i++) {
        auto it = binIndex.find(n[i]);

        if (it != binIndex.end()) {
            int index = it->second;

            if (isZero(sum(abunds[i]))) {
                vector<int> thisBinsSeqsToRemove = remove(n[i], reason);
                seqsRemoved.insert(seqsRemoved.end(),
                                   thisBinsSeqsToRemove.begin(),
                                   thisBinsSeqsToRemove.end());
            }else{
                binCount.setAbundance(index, abunds[i]);
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    binCount.updateTotals();

    return seqsRemoved;
}
/******************************************************************************/
void BinTable::updateBinAbunds(map<int, map<string, float>>& binAbunds,
               AbundTable& count, const int bIndex, const int seqIndex,
               const vector<string>& allSamples, bool firstTimeSeq) {

    bool hasSamples = !(allSamples.empty());

    auto it = binAbunds.find(bIndex);

    // first time seeing this bin, initialize sample counts
    if (it == binAbunds.end()) {
        map<string, float> thisSeqsSampleAbunds;

        // inputs from dataset
        vector<float> abunds = count.getAbundances(seqIndex);

        // count has samples
        if (hasSamples) {
            for (int j = 0; j < abunds.size(); j++) {
                if (!isZero(abunds[j])) {
                    thisSeqsSampleAbunds[allSamples[j]] = abunds[j];
                }
            }
        }else{
            thisSeqsSampleAbunds["total"] = abunds[0];
        }

        binAbunds[bIndex] = thisSeqsSampleAbunds;
    }else{
        // inputs from dataset
        if (firstTimeSeq) {
            vector<float> abunds = count.getAbundances(seqIndex);

            // count has samples
            if (hasSamples) {
                for (int j = 0; j < abunds.size(); j++) {
                    if (!isZero(abunds[j])) {
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
/******************************************************************************/
double BinTable::updateBins(map<int, map<string, float>>& binAbunds,
                          AbundTable& count, bool hasSamples) {

    double numBinsUpdated = 0;

    if (hasSamples) {

        vector<int> ids;
        vector<float> abundance;
        vector<string> samples;
        for (auto it = binAbunds.begin(); it != binAbunds.end(); it++) {

            map<string, float> abundances = it->second;
            for (auto itSamples = abundances.begin();
                 itSamples != abundances.end(); itSamples++) {
                ids.push_back(it->first);
                samples.push_back(itSamples->first);
                abundance.push_back(itSamples->second);
            }
        }

        numBinsUpdated = binCount.assignAbundance(ids, abundance, samples);

        // if there are treatments add them to bin
        if (count.getNumTreatments() != 0) {
            map<string, string> sampleTreatments =
                count.getSampleTreatmentAssignments();

            vector<string> samples(sampleTreatments.size());
            vector<string> treatments(sampleTreatments.size());

            int index = 0;
            for (auto itTreatments = sampleTreatments.begin();
                 itTreatments != sampleTreatments.end(); itTreatments++) {
                samples[index] = itTreatments->first;
                treatments[index] = itTreatments->second;
                index++;
            }
            assignTreatments(samples, treatments);
        }

    }else{
        vector<int> ids;
        vector<float> abundance;
        for (auto it = binAbunds.begin(); it != binAbunds.end(); it++) {
            ids.push_back(it->first);
            abundance.push_back(it->second["total"]);
        }

        numBinsUpdated = binCount.assignAbundance(ids, abundance);
    }
    return numBinsUpdated;
}
/******************************************************************************/
void BinTable::updateTotals() {
    binCount.updateTotals();
}
/******************************************************************************/

