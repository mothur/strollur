#include "../inst/include/strollur.h"
#include "dataset.h"
#include "phylotree.h"

/******************************************************************************/
BinTable::BinTable() {
    uniqueBad = 0;
    label = "";
    hasBinTaxonomy = false;
    hasBinReps = false;
    runClassify = false;
}
/******************************************************************************/
BinTable::BinTable(const BinTable& binTable) {
    label = binTable.label;
    hasBinTaxonomy = binTable.hasBinTaxonomy;
    hasBinReps = binTable.hasBinReps;

    binIndex = binTable.binIndex;
    tableBins = binTable.tableBins;
    binNames = binTable.binNames;
    trashCodes = binTable.trashCodes;
    taxonomies = binTable.taxonomies;
    repSequences = binTable.repSequences;
    runClassify = binTable.runClassify;
    originalBinAbunds = binTable.originalBinAbunds;

    binList = binTable.binList;
    seqBins = binTable.seqBins;

    badAccnos = binTable.badAccnos;
    uniqueBad = binTable.uniqueBad;
}
/******************************************************************************/
BinTable::~BinTable() {}
/******************************************************************************/
void BinTable::clear() {
        uniqueBad = 0;
        label = "";
        hasBinTaxonomy = false;
        hasBinReps = false;
        runClassify = false;

        binIndex.clear();
        tableBins.clear();
        binNames.clear();
        badAccnos.clear();
        trashCodes.clear();
        seqBins.clear();
        binList.clear();
        taxonomies.clear();
        repSequences.clear();
        originalBinAbunds.clear();
}
/******************************************************************************/
void BinTable::clone(const BinTable& binTable) {
    label = binTable.label;
    hasBinTaxonomy = binTable.hasBinTaxonomy;
    hasBinReps = binTable.hasBinReps;

    binIndex = binTable.binIndex;
    tableBins = binTable.tableBins;
    binNames = binTable.binNames;
    trashCodes = binTable.trashCodes;
    taxonomies = binTable.taxonomies;
    repSequences = binTable.repSequences;
    runClassify = binTable.runClassify;
    originalBinAbunds = binTable.originalBinAbunds;

    binList = binTable.binList;
    seqBins = binTable.seqBins;

    badAccnos = binTable.badAccnos;
    uniqueBad = binTable.uniqueBad;
}
/******************************************************************************/
double BinTable::assignBins(AbundTable& count, vector<string> binIds,
                                 const vector<int>& seqIds){
    // we are assigning seqIds to bins, so remove all old bin data
    // to avoid inconsistencies
    string oldLabel = label;
    clear();
    int numBins = 0;
    label = oldLabel;
    runClassify = true;

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

        // first time we are seeing this seqName, add seq to bin
        if (itSeq == seqBins.end()) {
            seqBins[seqIndex] = index;
            binList[index].insert(seqIndex);
        }
    }

    originalBinAbunds = getAbundances(count);

    return numBins;
}
/******************************************************************************/
double BinTable::assignRepresentativeSequences(const vector<string>& binNames,
                                                 const vector<int>& repNames){
    double repAssigned = 0;

    // set as unassigned
    repSequences.resize(tableBins.size(), -1);

    for (int i = 0; i < binNames.size(); i++) {

        string binName = binNames[i];

        auto it = binIndex.find(binName);

        // valid bin
        if (it != binIndex.end()) {
            repAssigned++;
            repSequences[it->second] = repNames[i];
        }
    }
    hasBinReps = true;

    return repAssigned;
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

    if (tax.size() != 0) {

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
Rcpp::List BinTable::exportBinTable(AbundTable& count) {
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
    binData.push_back(getAbundances(count, false));
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
    Rcpp::DataFrame binSeqAssignments = Rcpp::DataFrame::create(
        Rcpp::Named("bin_ids") = getValues(seqBins),
        Rcpp::_["sequence_ids"] = getKeys(seqBins));
    results.push_back(binSeqAssignments);
    resultLabels.push_back("sequence_bin_assignments");


    if (hasBinReps) {
        Rcpp::DataFrame binSeqReps = Rcpp::DataFrame::create(
            Rcpp::Named("bin_ids") = binIndexes,
            Rcpp::_["sequence_ids"] = repSequences);
        results.push_back(binSeqReps);
        resultLabels.push_back("bin_representative_sequences");
    }
    results.attr("names") = resultLabels;

    return results;
}
/******************************************************************************/
const vector<int> BinTable::getGoodIndexes(AbundTable& count) {

    vector<string> ids = getIds(count);

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

    return bin;
}
/******************************************************************************/
const vector<string> BinTable::getListVector(const vector<string>& seqNames) {
    vector<string> listVector;

    for (int i = 0; i < binNames.size(); i++) {
        if (tableBins[i]) {
            listVector.push_back(get(binNames[i], seqNames));
        }
    }

    return listVector;
}
/******************************************************************************/
// 2 column dataframe - bin_id, seq_id
const Rcpp::DataFrame BinTable::getList(const vector<string>& seqNames){

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
/******************************************************************************/
// names of OTUs
const vector<string> BinTable::getIds(AbundTable& count,
                                      vector<string> samples, bool distinct){

    vector<string> results;

    // no sample given, return all "good" bin names
    if (samples.empty()) {
        return select(binNames, tableBins);
    }else {
        // index of sample we are looking for
        vector<int> sampleIndexes(samples.size(), -1);
        if (distinct) {
            vector<string> sampleNames = count.getSamples();
            int next = 0;
            for (int i = 0; i < sampleNames.size(); i++) {
                // is this sample in the samples requested
                if (vectorContains(samples, sampleNames[i])) {
                    sampleIndexes[next] = i;
                    next++;
                }
            }
        }

        // for every bin
        for (int i = 0; i < binNames.size(); i++) {

            auto it = binIndex.find(binNames[i]);

            if (it != binIndex.end()) {

                // if this is a "good" bin
                if (tableBins[it->second]) {

                    int numSamplesFound = 0;
                    for (string sample : samples) {
                        vector<string> s; s.push_back(sample);
                        if (getAbundance(count, it->second, s) != 0) {
                            numSamplesFound++;
                        }
                    }

                    // includes ONLY sequences from these samples
                    if (distinct) {
                        const vector<float> sampleAbunds = getAbundances(count, it->second);

                        float theseSamples = 0;
                        for (int index : sampleIndexes) {
                            theseSamples += sampleAbunds[index];
                        }
                        // if all the sequences come from this sample, save name
                        if (isEqual(sum(sampleAbunds), theseSamples) &&
                                    (numSamplesFound == samples.size())) {
                            results.push_back(binNames[i]);
                        }
                    }else {
                        // all samples requested are present
                        if (numSamplesFound == samples.size()) {
                            results.push_back(binNames[i]);
                        }
                    }
                }
            }
        }

    }
    return results;
}
/******************************************************************************/
// 3 column dataframe - bin_id, abundance, sample
const Rcpp::DataFrame BinTable::getRepresentativeSequences(const vector<string>& seqNames,
                                                 const vector<string>& seqs) {

    if (hasBinReps) {

        vector<string> ids, seqids, seqDNA;
        for (int i = 0; i < binNames.size(); i++) {

            // if this is a "good" bin
            if (tableBins[i]) {
                int repSeq = repSequences[i];

                ids.push_back(binNames[i]);
                if (repSeq != -1) {
                    seqids.push_back(seqNames[repSeq]);
                    seqDNA.push_back(seqs[repSeq]);
                }else{
                    seqids.push_back("NA");
                    seqDNA.push_back("NA");
                }
            }
        }

        string tag = label + "_names";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(tag.c_str()) = ids,
            Rcpp::_["representative_names"] = seqids,
            Rcpp::_["representative_sequences"] = seqDNA);
        return df;
    }

    return Rcpp::DataFrame::create();
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
const Rcpp::DataFrame BinTable::getRAbund(AbundTable& count) {
    if (getNumBins(count)!= 0) {
        string tag = label + "_id";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named(tag) = select(binNames, tableBins),
            Rcpp::_["abundance"] = getAbundances(count)
        );

        return df;
    }
    return Rcpp::DataFrame::create();
}
/******************************************************************************/
// vector of total abundances for each outID
const vector<float> BinTable::getRAbundVector(AbundTable& count){
    return getAbundances(count);
}
/******************************************************************************/
const Rcpp::DataFrame BinTable::getShared(AbundTable& count) {
    if ((getNumBins(count) != 0) && (count.getNumSamples() != 0)) {

        vector<string> ids = getIds(count);
        vector<int> indexes(ids.size(), -1);

        for (int i = 0; i < ids.size(); i++) {
            indexes[i] = binIndex[ids[i]];
        }
        return getAbundanceTable(count);
    }

    return Rcpp::DataFrame::create();
}
/******************************************************************************/
// abundances for each OTU broken down by sample
const vector<vector<float> > BinTable::getSharedVector(AbundTable& count){

    vector<int> goodBins = getGoodIndexes(count);

    vector<vector<float> > bins(goodBins.size());

    for (int i = 0; i < goodBins.size(); i++) {
        bins[i] = getAbundances(count, goodBins[i]);
    }

    return bins;
}
/******************************************************************************/
// total abundance for a given binID, optional sample
const float BinTable::getAbundance(AbundTable& count,
                                   const string binId,
                                   vector<string> samples){

    auto it = binIndex.find(binId);

    if (it != binIndex.end()) {
        return getAbundance(count, it->second, samples);
    }

    return 0;
}
/******************************************************************************/
// total abundance for a given binID, optional sample
const float BinTable::getAbundance(AbundTable& count,
                                   const int binId,
                                   vector<string> samples){
    // if this is a "good" bin
    if (tableBins[binId]) {

        set<int> seqsInBin = binList[binId];
        float binAbund = 0;

        // for each seq in bin
        for (auto itSeq = seqsInBin.begin();
             itSeq != seqsInBin.end(); itSeq++) {

            binAbund += count.getAbundance(*itSeq, samples);
        }

        return binAbund;
    }

    return 0;
}
/******************************************************************************/
// total abundance for each bin
const vector<float> BinTable::getAbundances(AbundTable& count, bool onlyGood) {

    vector<string> namesToInclude;

    if (onlyGood) {
        namesToInclude = getIds(count);
    }else{
        namesToInclude = binNames;
    }

    vector<float> abunds(namesToInclude.size(), 0);

    // for each bin
    for (int i = 0; i < namesToInclude.size(); i++) {

        int index = binIndex[namesToInclude[i]];
        set<int> seqsInBin = binList[index];

        // for each seq in bin
        for (int seq : seqsInBin) {
            // add seq abunds to abunds[i]
            abunds[i] += count.getAbundance(seq);
        }
    }
    return abunds;
}
/******************************************************************************/
// abundances for given binID broken down by sample
const vector<float> BinTable::getAbundances(AbundTable& count,
                                            const string binName){
    auto it = binIndex.find(binName);
    if (it != binIndex.end()) {
        return getAbundances(count, it->second);
    }
    return nullFloatVector;
}
/******************************************************************************/
// abundances for given binID broken down by sample
const vector<float> BinTable::getAbundances(AbundTable& count,
                                            const int binId){
    // if this is a "good" bin
    if (tableBins[binId]) {

        vector<float> abundances(count.getNumSamples(), 0);

        set<int> seqsInBin = binList[binId];

        // for each seq in bin
        for (int seq : seqsInBin) {

            // get seq abunds
            vector<float> seqAbunds = count.getAbundances(seq);

            // add seq abunds to binAbunds
            sum(abundances, seqAbunds);
        }

        return abundances;
    }

    return nullFloatVector;
}
/******************************************************************************/
const Rcpp::DataFrame BinTable::getAbundanceTable(AbundTable& count,
                                                  const bool useNames) {

    if (count.hasSamplesData()) {

        bool hasTreatments = (count.getNumTreatments() != 0);
        vector<string> allSamples = count.getSamples();
        map<string, string> treatmentAssignments = count.getSampleTreatmentAssignments();

        vector<string> ids;
        vector<int> ids_index;
        vector<float> abunds;
        vector<string> samples;
        vector<string> treaments;

        vector<int> binIds = getGoodIndexes(count);

        // for each bin index
        for (int index : binIds) {
            string name = binNames[index];

            // binAbundances parsed by sample
            vector<float> binAbunds = getAbundances(count, index);

            for (int i = 0; i < binAbunds.size(); i++) {
                // nonzero sample
                if (binAbunds[i] != 0) {
                    if (useNames) {
                        ids.push_back(name);
                    }else{
                        ids_index.push_back(index);
                    }
                    abunds.push_back(binAbunds[i]);
                    samples.push_back(allSamples[i]);
                    if (hasTreatments) {
                        treaments.push_back(treatmentAssignments[allSamples[i]]);
                    }
                }
            }
        }

        if (useNames) {
            if (hasTreatments) {
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::Named("bin_names") = ids,
                    Rcpp::_["abundances"] = abunds,
                    Rcpp::_["samples"] = samples,
                    Rcpp::_["treatments"] = treaments);
                return df;
            }else{
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::Named("bin_names") = ids,
                    Rcpp::_["abundances"] = abunds,
                    Rcpp::_["samples"] = samples);
                return df;
            }
        }else {
            if (hasTreatments) {
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::Named("bin_ids") = ids_index,
                    Rcpp::_["abundances"] = abunds,
                    Rcpp::_["samples"] = samples,
                    Rcpp::_["treatments"] = treaments);
                return df;
            }else{
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::Named("bin_ids") = ids_index,
                    Rcpp::_["abundances"] = abunds,
                    Rcpp::_["samples"] = samples);
                return df;
            }
        }
    }else{

        if (useNames) {
            // no sample information
            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named("bin_names") = getIds(count),
                Rcpp::_["abundances"] = getAbundances(count));
            return df;
        }else{
            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named("bin_ids") = getGoodIndexes(count),
                Rcpp::_["abundances"] = getAbundances(count));
            return df;
        }

    }

    return Rcpp::DataFrame::create();
}
/******************************************************************************/
const int BinTable::getNumBins(AbundTable& count, vector<string> samples,
                               bool distinct) {
    if (!samples.empty()) {
        return getIds(count, samples, distinct).size();
    }

    // all "good" bins
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

    return Rcpp::DataFrame::create();
}
/******************************************************************************/
// type, trashCode, binCount, abundanceCount
const Rcpp::DataFrame BinTable::getScrapSummary() {

    vector<string> types, codes;
    vector<float> uniqueCounts, totalCounts;

    if (badAccnos.size() != 0) {
        types.resize(badAccnos.size());
        codes.resize(badAccnos.size());
        uniqueCounts.resize(badAccnos.size());
        totalCounts.resize(badAccnos.size());

        int index = 0;
        for (auto it = badAccnos.begin(); it != badAccnos.end(); it++) {
            codes[index] = it->first;
            uniqueCounts[index] = it->second[0];
            totalCounts[index] = it->second[1];
            types[index] = label;
            index++;
        }
    }

    if (!types.empty()) {
        return Rcpp::DataFrame::create(
            Rcpp::Named("type") = types,
            Rcpp::_["trash_code"] = codes,
            Rcpp::_["unique"] = uniqueCounts,
            Rcpp::_["total"] = totalCounts);
    }

    return Rcpp::DataFrame::create();
}
/******************************************************************************/
const bool BinTable::okToMerge(const vector<int>& seqIds) {

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

    return true;
}
/******************************************************************************/
void BinTable::merge(vector<string> binIDS, const string reason){
    if (binIDS.size() != 1) {

        vector<int> indexes = getIndexes(binIDS);

        for (int i = 1; i < indexes.size(); i++) {

            // no need to update the sample and treatment counts
            vector<int> seqIds = remove(binIDS[i], reason);

            // set each seqs bin assignment
            for (int seq : seqIds) {
                seqBins[seq] = indexes[0];
                binList[indexes[0]].insert(seq);
            }

            binList[indexes[i]].clear();
        }
    }
}
/******************************************************************************/
void BinTable::removeSeq(const int seqId, const string reason) {

    // if seqs were assigned to otus
    runClassify = true;

    // remove from otu
    auto it = seqBins.find(seqId);

    if (it != seqBins.end()) {

        int index = it->second;

        // does removing the seq remove a bin
        binList[index].erase(seqId);

        // remove bin, if empty
        if (binList[index].size() == 0) {
            remove(index, reason);
        }
    }

    if (hasBinReps) {
        auto it = std::find(repSequences.begin(), repSequences.end(), seqId);
        if (it != repSequences.end()) {
            *it = -1;
        }
    }
}
/******************************************************************************/
// remove given outID, returns seqs removed
vector<int> BinTable::remove(string binID, string reason){

    auto it = binIndex.find(binID);

    // bin is not in dataset
    if (it == binIndex.end()) {
        string message = "[WARNING]: " + binID + " is not in ";
        message += "your dataset, ignoring.";
        RcppThread::Rcout << endl << message << endl;
        return nullIntVector;
    }

    return (remove(it->second, reason));
}
/******************************************************************************/
// remove given outID, returns seqs removed
vector<int> BinTable::remove(int index, string reason){

    // remove from tableBins and add trashCode
    tableBins[index] = false;
    trashCodes[index] += reason;

    auto itBad = badAccnos.find(reason);

    if (itBad != badAccnos.end()) {
        // update counts of trashCode
        itBad->second[0]++;
        itBad->second[1] += originalBinAbunds[index];
    }else{
        // add new trashCode
        vector<double> badAbunds(2, 1);
        badAbunds[1] = originalBinAbunds[index];
        badAccnos[reason] = badAbunds;
    }

    // update uniqueBad
    uniqueBad++;

    if (hasBinReps) {
        repSequences[index] = -1;
    }

    vector<int> seqsToRemoveFromOtherBinTables = toVector(binList[index]);
    binList[index].clear();

    return seqsToRemoveFromOtherBinTables;
}
/******************************************************************************/

