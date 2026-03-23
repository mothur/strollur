#include "../inst/include/strollur.h"
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
double AbundTable::add(const vector<int>& names) {

    counts.resize(counts.size()+names.size());

    for (const int& name : names) {
        const sampleAbunds thisSeq(0, 1);
        counts[name] = thisSeq;
    }

    // make space for total abundance
    tableSamples.clear();
    tableSamples.push_back(true);
    numSamples = 1;

    tableTreatments.clear();
    tableTreatments.push_back(true);
    numTreatments = 1;

    const auto nameSize = static_cast<double>(names.size());
    total += nameSize;
    return nameSize;
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

    for (size_t i = 0; i < samples.size(); i++) {
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

    for (int i = 0; i < static_cast<int>(treatments.size()); i++) {
        treatmentIndex[treatments[i]] = i;
    }

    hasTreatments = true;
}
/******************************************************************************/
// names, abundances, samples (optional), treatment (optional)
double AbundTable::assignAbundance(vector<int>& names,
                                   const vector<float>& abunds,
                                   const vector<string>& samples,
                                   const vector<string>& treatments) {

    clear();

    if (samples.empty() && treatments.empty()) {
        hasSampleData = false;
        tableSamples.push_back(true);
        tableTreatments.push_back(true);
        numSamples = 1; // placeholder for sparse abundances
        numTreatments = 1;
    }else{
        const vector<string> uniqueSamples = unique(samples);
        const vector<string> uniqueTreatments = unique(treatments);

        if (!uniqueSamples.empty()) {
            addSamples(uniqueSamples);
            numSamples = uniqueSamples.size();

            if (!uniqueTreatments.empty()) {
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
        const vector<int> uniqueNames = unique(names);

        // fill in long format
        counts.resize(uniqueNames.size());
        for (size_t i = 0; i < samples.size(); i++) {
            int index = sampleIndex[samples[i]];
            counts[names[i]].sampleIndex.push_back(index);
            counts[names[i]].abunds.push_back(abunds[i]);
            sampleTotals[index] += abunds[i];
        }

        // convert to sparse
        total = sum(sampleTotals);
    }else{

        counts.resize(names.size());
        for (size_t i = 0; i < names.size(); i++) {
            const sampleAbunds thisSeq(0, abunds[i]);
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

    return static_cast<double>(counts.size());
}
/******************************************************************************/
double AbundTable::assignTreatments(const vector<string>& samples,
                                    const vector<string>& treatments) {

    double numTreatmentsAssigned = 0;

   if (hasSampleData) {
       if (hasTreatments) {
           treatmentIndex.clear();
           treatmentTotals.clear();
           tableTreatments.clear();
           sampleTreatment.clear();
       }

       const vector<string> uniqueTreatments = unique(treatments);

       if (uniqueTreatments.size() > 1) {
           addTreatments(uniqueTreatments);
           numTreatments = uniqueTreatments.size();
       }else if ((uniqueTreatments.size() == 1) &&
           (uniqueTreatments[0] != "")) {
           addTreatments(uniqueTreatments);
           numTreatments = 1;
       }

        // create sample to treatment map
       for (size_t i = 0; i < samples.size(); i++) {
           auto itSample = sampleIndex.find(samples[i]);


           // valid sample
           if (itSample != sampleIndex.end()) {
               sampleTreatment[itSample->second] = treatments[i];
               numTreatmentsAssigned++;
           }else {
               string message = "[WARNING]: The dataset does not contain sample, '";
               message += samples[i] + "'. Ignoring '" + samples[i] + "'.";
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
       for (size_t i = 0; i < treatmentTotals.size(); i++) {
           if (isZero(treatmentTotals[i])) {
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
const float AbundTable::getAbundance(const int name, const vector<string>& samples) {

    float abund = 0;

    if (samples.empty()) {
        if (hasSampleData) {
            const vector<float> abunds = getAbundances(name);
            abund = sum(abunds);
        }else{
            abund = counts[name].abunds[0];
        }
    }else if (hasSamples(samples)) {

        for (const string& sample : samples) {
            const int gIndex = sampleIndex[sample];

            // if the sequence does not have this sample, -1 returned
            const int thisSamplesIndex = getSparseIndex(name, gIndex);

            if (thisSamplesIndex != -1) {
                abund += counts[name].abunds[thisSamplesIndex];
            }
        }
    }

    return abund;
}
/******************************************************************************/
// abundances by sample, in the same order as the samples
const vector<float> AbundTable::getAbundances(const int id) const {
    vector<float> abunds;

    if (hasSampleData) {
        // all samples not just "good" ones
        vector<float> allAbunds(tableSamples.size(), 0);

        const sampleAbunds data = counts[id];

        // data -> sampleIndex(2,5), abunds(100, 50)
        // becomes abunds(0,0,100,0,0,50)
        for (size_t i = 0; i < data.sampleIndex.size(); i++) {
            allAbunds[data.sampleIndex[i]] = data.abunds[i];
        }

        // only include "good" samples
        abunds = select(allAbunds, tableSamples);
    }else{
        abunds.push_back(counts[id].abunds[0]);
    }

    return abunds;
}
/******************************************************************************/
const vector<vector<float>> AbundTable::getAbundances(const vector<int>& ids) const {
    vector<vector<float>> results(ids.size());

    for (size_t i = 0; i < ids.size(); i++) {
        results[i] = getAbundances(ids[i]);
    }

    return results;
}
/******************************************************************************/
const vector<vector<float>> AbundTable::getAbundanceBySample(const vector<int>& ids,
                                                 vector<string> samplesToSelect) {

    if (samplesToSelect.empty()) {
        samplesToSelect = getSamples();
    }

    vector<vector<float>> results(samplesToSelect.size());

    if (!hasSampleData) { return results; }

    for (const int& id : ids) {
        sampleAbunds data = counts[id];

        // data -> sampleIndex(2,5), abunds(100, 50)
        // becomes abunds(0,0,100,0,0,50)
        for (size_t i = 0; i < data.sampleIndex.size(); i++) {

            // this is a "good" sample
            if (tableSamples[data.sampleIndex[i]]) {
                // samplesIndexInResults -> sampleNames[data.sampleIndex[i]]]
                results[sampleIndex[sampleNames[data.sampleIndex[i]]]].push_back(data.abunds[i]);
            }
        }
    }

    return results;

}
/******************************************************************************/
const vector<string> AbundTable::getSamples(const int name) {
    vector<string> samples;

    if (hasSampleData) {

        // want all "good" samples
        if (name == -1) {
            return select(sampleNames, tableSamples);
        }else {
            const vector<int> sampleIndexes = counts[name].sampleIndex;
            for (const int& index : sampleIndexes) {
                // only count "good" samples
                if (tableSamples[index]) {
                    samples.push_back(sampleNames[index]);
                }
            }
        }
    }

    return samples;
}
/******************************************************************************/
const vector<int> AbundTable::getSampleIndexes() {
    vector<int> sampleIndexes;

    if (hasSampleData) {

        const vector<string> samples = getSamples();
        sampleIndexes.resize(samples.size(), 0);
        for (size_t i = 0; i < samples.size(); i++) {
            sampleIndexes[i] = sampleIndex[samples[i]];
        }

    }else {
        sampleIndexes.resize(1, 0);
    }

    return sampleIndexes;
}
/******************************************************************************/
const vector<double> AbundTable::getSampleTotals() const{
    vector<double> totals;
    if (hasSampleData) {
        totals = select(sampleTotals, tableSamples);
    }
    return totals;
}
/******************************************************************************/
const map<string, string> AbundTable::getSampleTreatmentAssignments() const {
    map<string, string> results;

    if (hasTreatments) {
        for (const auto & sampleName : sampleNames) {
            int index = sampleIndex.at(sampleName);
            if (tableSamples[index]) {
                const string treatmentName = sampleTreatment.at(index);
                results[sampleName] = treatmentName;
            }
        }
    }

    return results;
}
/******************************************************************************/
const int AbundTable::getNumSamples(const int name) {
    if (hasSampleData) {
        if (name == -1) {
           return numSamples;
        }else {
           vector<int> sampleIndexes = counts[name].sampleIndex;
           int num = 0;
           for (const int& index : sampleIndexes) {
               // only count "good" samples
               if (tableSamples[index]) {
                   num++;
               }
           }
           return num;
        }
    }
    return 0;
}
/******************************************************************************/
const int AbundTable::getNumTreatments() const {
    if (hasTreatments) {
        return numTreatments;
    }
    return 0;
}
/******************************************************************************/
// vector containing total abundance for each sequence
const vector<float> AbundTable::getTotalAbundances(const vector<int>& names) {

    vector<float> abunds(names.size(), 0);

    for (size_t i = 0; i < names.size(); i++) {
        abunds[i] = getAbundance(names[i]);
    }
    return abunds;
}
/******************************************************************************/
const Rcpp::DataFrame AbundTable::getAbundanceTable(const vector<string>& outputNames,
                                          const vector<int>& names,
                                          const string& tag, const bool useNames) {


    if (hasSampleData) {

        vector<string> ids;
        vector<int> ids_index;
        vector<float> abunds;
        vector<string> samples;
        vector<string> treaments;

        for (size_t i = 0; i < names.size(); i++) {
            const string& name = outputNames[i];

            const sampleAbunds data = counts[names[i]];

            // data -> sampleIndex(2,5), abunds(100, 50)
            // becomes abunds(0,0,100,0,0,50)
            for (size_t j = 0; j < data.sampleIndex.size(); j++) {
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
const int AbundTable::getSparseIndex(const int name, const int sample) {
    constexpr int index = -1;

    // if this is a removed sample
    if (!tableSamples[sample]) {
        return index;
    }

    for (size_t i = 0; i < counts[name].sampleIndex.size(); i++) {
        if (counts[name].sampleIndex[i] == sample) { return i; }
    }

    return index;
}
/******************************************************************************/
const double AbundTable::getTotal(const string& sample) const {
    if (!sample.empty()) {
        const auto it = sampleIndex.find(sample);

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
const vector<string> AbundTable::getTreatments() const {
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
const vector<double> AbundTable::getTreatmentTotals() {
    vector<double> totals;
    if (hasTreatments) {
        totals = select(treatmentTotals, tableTreatments);
    }
    return totals;
}
/******************************************************************************/
const bool AbundTable::hasSample(const string sample, const int name) {
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
const bool AbundTable::hasSamples(const vector<string> samples, const int name) {

    // must have all samples asked for
    for (string sample : samples) {
        if (!hasSample(sample, name)) {
            return false;
        }
    }
    return true;
}
/******************************************************************************/
// adds sequences counts of idsToMerge[1-n] into idsToMerge[0]
void AbundTable::merge(const vector<int>& idsToMerge) {
    if (idsToMerge.size() > 1) {

        int keeperSeq = idsToMerge[0];
        vector<float> kAbunds = getAbundances(keeperSeq);

        for (int i = 1; i < idsToMerge.size(); i++) {
            int dupSeq = idsToMerge[i];
            vector<float> dupAbunds = getAbundances(dupSeq);
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
        vector<float> abunds = getAbundances(name);
        abund = sum(abunds);

        // remove seq from sample / treatment totals
        if (!isZero(abund)) {
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
void AbundTable::removeSamples(const vector<string>& samples) {
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
// set abundance
void AbundTable::setAbundance(const int name, const vector<float>& abunds) {

    vector<float> origAbunds = getAbundances(name);

    sampleAbunds newAbunds(getSampleIndexes(), abunds, "full");
    counts[name] = newAbunds;

    if (hasSampleData) {
        int index = 0;
        // sampleTotals may be larger than abunds
        // if samples have been removed
        float diffAbund = 0;
        for (int i = 0; i < sampleTotals.size(); i++) {
            if (tableSamples[i]) {
                float diff = origAbunds[index] - abunds[index];
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
// set abundance - for single sample
void AbundTable::setAbundance(int name, float abund, string sample) {

    if (!hasSampleData && (sample == "")) {
        float origAbund = getAbundance(name);

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

                float orig = 0;
                if (thisSamplesIndex != -1) {
                    // update existing abunds
                    orig = counts[name].abunds[thisSamplesIndex];
                    counts[name].abunds[thisSamplesIndex] = abund;
                }else{
                    // add new abund
                    counts[name].abunds.push_back(abund);
                    counts[name].sampleIndex.push_back(index);
                }

                float diff = orig - abund;
                sampleTotals[index] -= diff;
                total -= diff;
            }
        }
    }
}
/******************************************************************************/
void AbundTable::setAbundances(const map<int, map<string, float>>& binAbunds) {

    // just abunds no parse by sample
    if (hasSampleData) {
        for (auto it = binAbunds.begin(); it != binAbunds.end(); it ++) {

            vector<float> abunds(numSamples, 0);
            for (auto itSample = it->second.begin();
                 itSample != it->second.end(); itSample++) {

                 // valid sample
                 auto itIndex = sampleIndex.find(itSample->first);

                if (itIndex != sampleIndex.end()) {
                    int index = itIndex->second;

                    abunds[index] = itSample->second;
                }
            }

            setAbundance(it->first, abunds);
        }
    }else{
        for (auto it = binAbunds.begin(); it != binAbunds.end(); it ++) {
            setAbundance(it->first, it->second.begin()->second);
        }
    }

    updateTotals();
}
/******************************************************************************/
void AbundTable::updateTotals() {

    if (hasSampleData) {

    // does removing this sequence remove a sample
    for (int i = 0; i < sampleTotals.size(); i++) {

        if (tableSamples[i] && (isZero(sampleTotals[i]))) {
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
            if (tableTreatments[i] && (isZero(treatmentTotals[i]))) {
                tableTreatments[i] = false;
                numTreatments--;
            }
        }
    }
    }
}
/******************************************************************************/
