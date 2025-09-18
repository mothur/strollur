#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
AbundTable::AbundTable() {
    hasSampleData = false;
    hasTreatments = false;
    total = 0;
    numSamples = 0;
    numTreatments = 0;
}
/******************************************************************************/
AbundTable::~AbundTable() {}
/******************************************************************************/
void AbundTable::clear() {
    hasSampleData = false;
    hasTreatments = false;
    total = 0;
    numSamples = 0;
    numTreatments = 0;

    counts.clear();
    sampleIndex.clear();
    sampleNames.clear();
    sampleTotals.clear();
    tableSamples.clear();

    treatmentIndex.clear();
    treatmentTotals.clear();
    tableTreatments.clear();
    sampleTreatment.clear();
}
/******************************************************************************/
void AbundTable::clone(const AbundTable& abundTable) {
    counts = abundTable.counts;

    total = abundTable.total;
    numSamples = abundTable.numSamples;
    numTreatments = abundTable.numTreatments;

    sampleIndex = abundTable.sampleIndex;
    sampleNames = abundTable.sampleNames;
    sampleTotals = abundTable.sampleTotals;
    tableSamples = abundTable.tableSamples;

    treatmentIndex = abundTable.treatmentIndex;
    treatmentTotals = abundTable.treatmentTotals;
    tableTreatments = abundTable.tableTreatments;
    sampleTreatment = abundTable.sampleTreatment;

    hasSampleData = abundTable.hasSampleData;
    hasTreatments = abundTable.hasTreatments;
}
/******************************************************************************/
// the names are the indexes in dataset
double AbundTable::add(vector<int>& names) {

    counts.resize(counts.size()+names.size());

    for (int i = 0; i < names.size(); i++) {
        sampleAbunds thisSeq(0, 1);
        counts[names[i]] = thisSeq;
    }

    // make space for total abundance
    tableSamples.clear();
    tableSamples.push_back(true);
    numSamples = 1;

    tableTreatments.clear();
    tableTreatments.push_back(true);
    numTreatments = 1;

    total += names.size();
    return names.size();
}
/******************************************************************************/
// private function
void AbundTable::addSamples(vector<string> samples) {
    // store samples alphabetically
    sort(samples.begin(), samples.end());
    sampleNames = samples;

    // all samples start off as "good"
    tableSamples.resize(samples.size(), true);
    sampleTotals.resize(samples.size(), 0);

    for (int i = 0; i < samples.size(); i++) {
        sampleIndex[samples[i]] = i;
    }

    hasSampleData = true;
}
/******************************************************************************/
// private function
void AbundTable::addTreatments(vector<string> treatments) {
    // store samples alphabetically
    sort(treatments.begin(), treatments.end());

    // all samples start off as "good"
    tableTreatments.resize(treatments.size(), true);
    treatmentTotals.resize(treatments.size(), 0);

    for (int i = 0; i < treatments.size(); i++) {
        treatmentIndex[treatments[i]] = i;
    }

    hasTreatments = true;
}
/******************************************************************************/
// names, abundances, samples (optional), treatment (optional)
double AbundTable::assignAbundance(vector<int> names,
                           vector<int> abunds,
                           vector<string> samples,
                           vector<string> treatments) {

    clear();

    if ((samples.size() == 0) && (treatments.size() == 0)) {
        hasSampleData = false;
        tableSamples.push_back(true);
        tableTreatments.push_back(true);
        numSamples = 1; // placeholder for sparse abundances
        numTreatments = 1;
    }else{
        vector<string> uniqueSamples = unique(samples);
        vector<string> uniqueTreatments = unique(treatments);

        if (uniqueSamples.size() >= 1) {
            addSamples(uniqueSamples);
            numSamples = uniqueSamples.size();

            if (uniqueTreatments.size() >= 1) {
                addTreatments(uniqueTreatments);
                numTreatments = uniqueTreatments.size();
            }else {
                tableTreatments.push_back(true);
                numTreatments = 1;
            }
        }
    }

    // update counts, sampleTotals and total
    if (hasSampleData) {
        vector<int> uniqueNames = unique(names);

        // fill in long format
        counts.resize(uniqueNames.size());
        for (int i = 0; i < samples.size(); i++) {
            int index = sampleIndex[samples[i]];
            counts[names[i]].sampleIndex.push_back(index);
            counts[names[i]].abunds.push_back(abunds[i]);
            sampleTotals[index] += abunds[i];
        }

        // convert to sparse
        total = sum(sampleTotals);
    }else{

        counts.resize(names.size());
        for (int i = 0; i < names.size(); i++) {
            sampleAbunds thisSeq(0, abunds[i]);
            counts[names[i]] = thisSeq;
        }
        total = sum(abunds);
    }

    if (hasTreatments) {
        // create sample to treatment map
        for (int i = 0; i < samples.size(); i++) {
            sampleTreatment[sampleIndex[samples[i]]] = treatments[i];
        }

        // fill treatment totals
        for (auto it = sampleIndex.begin(); it != sampleIndex.end(); it++) {
            string treatment = sampleTreatment[it->second];
            int tindex = treatmentIndex[treatment];

            treatmentTotals[tindex] += sampleTotals[it->second];
        }
    }

    return counts.size();
}
/******************************************************************************/
double AbundTable::assignTreatments(vector<string> s,
                                  vector<string> t) {

    double numTreatmentsAssigned = 0;

   if (hasSampleData) {
       if (hasTreatments) {
           treatmentIndex.clear();
           treatmentTotals.clear();
           tableTreatments.clear();
           sampleTreatment.clear();
       }

       vector<string> uniqueTreatments = unique(t);

       if (uniqueTreatments.size() > 1) {
           addTreatments(uniqueTreatments);
           numTreatments = uniqueTreatments.size();
       }else if ((uniqueTreatments.size() == 1) &&
           (uniqueTreatments[0] != "")) {
           addTreatments(uniqueTreatments);
           numTreatments = 1;
       }

        // create sample to treatment map
       for (int i = 0; i < s.size(); i++) {
           auto itSample = sampleIndex.find(s[i]);


           // valid sample
           if (itSample != sampleIndex.end()) {
               sampleTreatment[itSample->second] = t[i];
               numTreatmentsAssigned++;
           }else {
               string message = "[WARNING]: The dataset does not contain sample ,";
               message += s[i] + "'. Ignoring '" + s[i] + "'.";
               RcppThread::Rcout << endl << message << endl;
           }
       }

       // fill treatment totals
       for (auto it = sampleIndex.begin(); it != sampleIndex.end(); it++) {
           string treatment = sampleTreatment[it->second];
           int tindex = treatmentIndex[treatment];

           treatmentTotals[tindex] += sampleTotals[it->second];
       }

       // check for treatments with 0 abundance, this can happen if
       // you assign and sample not in dataset to a treatment
       for (int i = 0; i < treatmentTotals.size(); i++) {
           if (treatmentTotals[i] == 0) {
               tableTreatments[i] = false;
               numTreatments--;
           }
       }
   }

   return numTreatmentsAssigned;
}
/******************************************************************************/
// total abundance for sequence.
// If sample provided, then abundance for sequence in sample
int AbundTable::getAbundance(int name, string sample) {

    int abund = 0;

    if (sample == "") {
        if (hasSampleData) {
            vector<int> abunds = getAbundances(name);
            abund = sum(abunds);
        }else{
            abund = counts[name].abunds[0];
        }
    }else if (hasSample(sample)) {
        int gIndex = sampleIndex[sample];

        // if the sequence does not have this sample, -1 returned
        int thisSamplesIndex = getSparseIndex(name, gIndex);

        if (thisSamplesIndex != -1) {
            abund = counts[name].abunds[thisSamplesIndex];
        }
    }

    return abund;
}
/******************************************************************************/
// abundances by sample, in the same order as the samples
vector<int> AbundTable::getAbundances(int name) {
    vector<int> abunds;

    if (hasSampleData) {
        // all samples not just "good" ones
        vector<int> allAbunds(tableSamples.size(), 0);

        sampleAbunds data = counts[name];

        // data -> sampleIndex(2,5), abunds(100, 50)
        // becomes abunds(0,0,100,0,0,50)
        for (int i = 0; i < data.sampleIndex.size(); i++) {
            allAbunds[data.sampleIndex[i]] = data.abunds[i];
        }

        // only include "good" samples
        abunds = select(allAbunds, tableSamples);
    }else{
        abunds.push_back(counts[name].abunds[0]);
    }

    return abunds;
}
/******************************************************************************/
vector<vector<int>> AbundTable::getAbundances(vector<int> ids) {
    vector<vector<int>> results(ids.size());

    for (int i = 0; i < ids.size(); i++) {
        results[i] = getAbundances(ids[i]);
    }

    return results;
}
/******************************************************************************/
vector<string> AbundTable::getSamples() {
    vector<string> samples;

    if (hasSampleData) {
        // want all "good" samples
        return select(sampleNames, tableSamples);
    }

    return samples;
}
/******************************************************************************/
vector<int> AbundTable::getSampleIndexes() {
    vector<int> sampleIndexes;

    if (hasSampleData) {

        vector<string> s = getSamples();
        sampleIndexes.resize(s.size(), 0);
        for (int i = 0; i < s.size(); i++) {
            sampleIndexes[i] = sampleIndex[s[i]];
        }

    }else {
        sampleIndexes.resize(1, 0);
    }

    return sampleIndexes;
}
/******************************************************************************/
vector<int> AbundTable::getSampleTotals() {
    vector<int> totals;
    if (hasSampleData) {
        totals = select(sampleTotals, tableSamples);
    }
    return totals;
}
/******************************************************************************/
map<string, string> AbundTable::getSampleTreatmentAssignments() {
    map<string, string> results;

    if (hasTreatments) {
        for (int i = 0; i < sampleNames.size(); i++) {
            string treatmentName = sampleTreatment[sampleIndex[sampleNames[i]]];
            results[sampleNames[i]] = treatmentName;
        }
    }

    return results;
}
/******************************************************************************/
int AbundTable::getNumSamples() {
    if (hasSampleData) {
        return numSamples;
    }
    return 0;
}
/******************************************************************************/
int AbundTable::getNumTreatments() {
    if (hasTreatments) {
        return numTreatments;
    }
    return 0;
}
/******************************************************************************/
// vector containing total abundance for each sequence
vector<int> AbundTable::getTotalAbundances(vector<int> names) {

    vector<int> abunds(names.size(), 0);

    for (int i = 0; i < names.size(); i++) {
        abunds[i] = getAbundance(names[i]);
    }
    return abunds;
}
/******************************************************************************/
Rcpp::DataFrame AbundTable::getAbundanceTable(const vector<string>& outputNames,
                                          const vector<int>& names,
                                          const string tag, const bool useNames) {


    if (hasSampleData) {

        vector<string> ids;
        vector<int> ids_index;
        vector<int> abunds;
        vector<string> samples;
        vector<string> treaments;

        for (int i = 0; i < names.size(); i++) {
            string name = outputNames[i];

            sampleAbunds data = counts[names[i]];

            // data -> sampleIndex(2,5), abunds(100, 50)
            // becomes abunds(0,0,100,0,0,50)
            for (int j = 0; j < data.sampleIndex.size(); j++) {
                //allAbunds[data.sampleIndex[i]] = data.abunds[i];
                // if this sample "good"
                if (tableSamples[data.sampleIndex[j]]) {
                    if (useNames) {
                        ids.push_back(name);
                    }else{
                        ids_index.push_back(names[i]);
                    }
                    abunds.push_back(data.abunds[j]);
                    samples.push_back(sampleNames[data.sampleIndex[j]]);

                    if (hasTreatments) {
                        treaments.push_back(sampleTreatment[data.sampleIndex[j]]);
                    }
                }else if (!useNames) { // want eliminated groups too for export
                    ids_index.push_back(names[i]);
                    abunds.push_back(data.abunds[j]);
                    samples.push_back(sampleNames[data.sampleIndex[j]]);

                    if (hasTreatments) {
                        treaments.push_back(sampleTreatment[data.sampleIndex[j]]);
                    }
                }
            }
        }

        if (useNames) {
            if (hasTreatments) {
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::Named(tag+"_names") = ids,
                    Rcpp::_["abundances"] = abunds,
                    Rcpp::_["samples"] = samples,
                    Rcpp::_["treatments"] = treaments);
                return df;
            }else{
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::Named(tag+"_names") = ids,
                    Rcpp::_["abundances"] = abunds,
                    Rcpp::_["samples"] = samples);
                return df;
            }
        }else {
            if (hasTreatments) {
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::Named(tag+"_ids") = ids_index,
                    Rcpp::_["abundances"] = abunds,
                    Rcpp::_["samples"] = samples,
                    Rcpp::_["treatments"] = treaments);
                return df;
            }else{
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::Named(tag+"_ids") = ids_index,
                    Rcpp::_["abundances"] = abunds,
                    Rcpp::_["samples"] = samples);
                return df;
            }
        }
    }else{

        if (useNames) {
            // no sample information
            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named(tag+"_names") = outputNames,
                Rcpp::_["abundances"] = getTotalAbundances(names));
            return df;
        }else{
            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named(tag+"_ids") = names,
                Rcpp::_["abundances"] = getTotalAbundances(names));
            return df;
        }

    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
int AbundTable::getSparseIndex(int name, int sample) {
    int index = -1;

    for (int i = 0; i < counts[name].sampleIndex.size(); i++) {
        if (counts[name].sampleIndex[i] == sample) { return i; }
    }

    return index;
}
/******************************************************************************/
int AbundTable::getTotal(string sample) {
    if (sample != "") {
        auto it = sampleIndex.find(sample);

        // is this a sample in the table
        if (it != sampleIndex.end()) {
            // is it "good"
            if (tableSamples[it->second]) {
                return sampleTotals[it->second];
            }
        }
        // dataset does not include the sample
        return 0;
    }
    return total;
}
/******************************************************************************/
vector<string> AbundTable::getTreatments() {
    vector<string> treatments;

    if (hasTreatments) {
        for (auto it = treatmentIndex.begin(); it != treatmentIndex.end(); it++) {
            if (tableTreatments[it->second]) {
                treatments.push_back(it->first);
            }
        }
    }

    return treatments;
}
/******************************************************************************/
vector<int> AbundTable::getTreatmentTotals() {
    vector<int> totals;
    if (hasTreatments) {
        totals = select(treatmentTotals, tableTreatments);
    }
    return totals;
}
/******************************************************************************/
bool AbundTable::hasSample(string sample, int name) {
    auto it = sampleIndex.find(sample);

    // is this a sample in the table
    if (it != sampleIndex.end()) {

        // no name provided
        if (name == -1) {
            // is it "good"
            return tableSamples[it->second];
        }else{
            // if the sequence does not have this sample, -1 returned
            int thisSamplesIndex = getSparseIndex(name, it->second);

            if (thisSamplesIndex != -1) {
                return true;
            }
        }
    }

    return false;
}
/******************************************************************************/
// adds sequences counts of idsToMerge[1-n] into idsToMerge[0]
void AbundTable::merge(vector<int> idsToMerge) {
    if (idsToMerge.size() > 1) {

        int keeperSeq = idsToMerge[0];
        vector<int> kAbunds = getAbundances(keeperSeq);

        for (int i = 1; i < idsToMerge.size(); i++) {
            int dupSeq = idsToMerge[i];
            vector<int> dupAbunds = getAbundances(dupSeq);
            sum(kAbunds, dupAbunds);
        }

        sampleAbunds newKeeper(getSampleIndexes(), kAbunds, "full");
        counts[keeperSeq] = newKeeper;
    }
}
/******************************************************************************/
int AbundTable::remove(int name) {

    int abund = 0;

    if (hasSampleData) {
        vector<int> abunds = getAbundances(name);
        abund = sum(abunds);

        // remove seq from sample / treatment totals
        if (abund != 0) {
            int index = 0;
            // sampleTotals may be larger than diffAbunds
            // if samples have been removed
            for (int i = 0; i < sampleTotals.size(); i++) {
                if (tableSamples[i]) {
                    sampleTotals[i] -= abunds[index];
                    index++;
                }
            }
        }

        // remove seq from total
        total -= abund;
    }else{
        // remove seq from total
        abund = getAbundance(name);
        total -= abund;
    }

    return abund;
}
/******************************************************************************/
void AbundTable::removeSamples(vector<string> samples) {
    if (hasSampleData) {

        for (int i = 0; i < samples.size(); i++) {

            auto it = sampleIndex.find(samples[i]);

            if (it != sampleIndex.end()) {
                tableSamples[it->second] = false;
                numSamples--;
                sampleTotals[it->second] = 0;
            }
        }

        // if no samples left
        if (isFalse(tableSamples)) {
            hasSampleData = false;
        }

        updateTotals();
    }
}
/******************************************************************************/
// set abundance parsed by sample - for datasets with samples
void AbundTable::setAbundance(int name, vector<int> abunds) {

    vector<int> origAbunds = getAbundances(name);

    sampleAbunds newAbunds(getSampleIndexes(), abunds, "full");
    counts[name] = newAbunds;

    if (hasSampleData) {
        int index = 0;
        // sampleTotals may be larger than abunds
        // if samples have been removed
        int diffAbund = 0;
        for (int i = 0; i < sampleTotals.size(); i++) {
            if (tableSamples[i]) {
                int diff = origAbunds[index] - abunds[index];
                sampleTotals[i] -= diff;
                diffAbund += diff;
                index++;
            }
        }

        total -= diffAbund;
    }else{
        total -= (origAbunds[0] - abunds[0]);
    }
}
/******************************************************************************/
// set abundance - for datasets without samples
void AbundTable::setAbundance(int name, int abund, string sample) {

    if (!hasSampleData && (sample == "")) {
        int origAbund = getAbundance(name);

        sampleAbunds thisCount(0, abund);
        counts[name] = thisCount;

        total -= (origAbund - abund);
    }else if (hasSampleData && (sample != "")) {
        // valid sample
        auto it = sampleIndex.find(sample);

        if (it != sampleIndex.end()) {
            int index = it->second;

            // "good" sample
            if (tableSamples[index]) {
                // if the sequence does not have this sample, -1 returned
                int thisSamplesIndex = getSparseIndex(name, index);

                int orig = 0;
                if (thisSamplesIndex != -1) {
                    // update existing abunds
                    orig = counts[name].abunds[thisSamplesIndex];
                    counts[name].abunds[thisSamplesIndex] = abund;
                }else{
                    // add new abund
                    counts[name].abunds.push_back(abund);
                    counts[name].sampleIndex.push_back(index);
                }

                int diff = orig - abund;
                sampleTotals[index] -= diff;
                total -= diff;
            }
        }
    }
}
/******************************************************************************/
void AbundTable::updateTotals() {

    if (hasSampleData) {

    // does removing this sequence remove a sample
    for (int i = 0; i < sampleTotals.size(); i++) {

        if (tableSamples[i] && (sampleTotals[i] == 0)) {
            tableSamples[i] = false;
            numSamples--;
        }
    }

    total = sum(sampleTotals);

    if (hasTreatments) {

        // reset to 0
        fill(treatmentTotals.begin(), treatmentTotals.end(), 0);

        // fill treatment totals
        for (auto it = sampleIndex.begin(); it != sampleIndex.end(); it++) {

            // if this sample is "good", add to treatment totals
            if (tableSamples[it->second]) {
                string treatment = sampleTreatment[it->second];
                int tindex = treatmentIndex[treatment];

                treatmentTotals[tindex] += sampleTotals[it->second];
            }
        }

        // does removing this sequence remove a treatment
        for (int i = 0; i < treatmentTotals.size(); i++) {
            if (tableTreatments[i] && (treatmentTotals[i] == 0)) {
                tableTreatments[i] = false;
                numTreatments--;
            }
        }
    }
    }
}
/******************************************************************************/
