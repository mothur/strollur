#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
SeqAbundTable::SeqAbundTable() {
    hasSampleData = false;
    hasTreatments = false;
    total = 0;
    numSamples = 0;
    numTreatments = 0;
}
/******************************************************************************/
SeqAbundTable::~SeqAbundTable() {}
/******************************************************************************/
void SeqAbundTable::clear() {
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
// the names are the indexes in dataset
void SeqAbundTable::addSeqs(vector<int>& names) {

    if (hasSampleData) {
        string message = "[ERROR]: The dataset contains sample information, ";
        message += "you must provide by sample data.\n\n";
        RcppThread::Rcout << endl << message;
    }else{

        counts.resize(counts.size()+names.size());

        for (int i = 0; i < names.size(); i++) {
            seqCount thisSeq(0, 1);
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
    }
}
/******************************************************************************/
// private function
void SeqAbundTable::addSamples(vector<string> samples) {
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
void SeqAbundTable::addTreatments(vector<string> treatments) {
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
void SeqAbundTable::assignSequenceAbundance(vector<int> names,
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

        if (uniqueSamples.size() > 1) {
            addSamples(uniqueSamples);
            numSamples = uniqueSamples.size();

            if (uniqueTreatments.size() > 1) {
                addTreatments(uniqueTreatments);
                numTreatments = uniqueTreatments.size();
            }else if ((uniqueTreatments.size() == 1) &&
                (uniqueTreatments[0] != "")) {
                addTreatments(uniqueTreatments);
                numTreatments = 1;
            }else {
                tableTreatments.push_back(true);
                numTreatments = 1;
            }
        }else if ((uniqueSamples.size() == 1) && (uniqueSamples[0] != "")) {
            addSamples(uniqueSamples);
            numSamples = 1;
            numTreatments = 1;
        }else{
            // empty strings passed in for samples, ignore
            hasSampleData = false;
            hasTreatments = false;
            tableSamples.push_back(true);
            tableTreatments.push_back(true);
            numSamples = 1; // placeholder for sparse abundances
            numTreatments = 1;
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
        total = accumulate(sampleTotals.begin(), sampleTotals.end(), 0);
    }else{

        counts.resize(names.size());
        for (int i = 0; i < names.size(); i++) {
            seqCount thisSeq(0, abunds[i]);
            counts[names[i]] = thisSeq;
        }
        total = accumulate(abunds.begin(), abunds.end(), 0);
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
}
/******************************************************************************/
// total abundance for sequence.
// If sample provided, then abundance for sequence in sample
int SeqAbundTable::getAbund(int name, string sample) {

    int abund = 0;

    if (sample == "") {
        if (hasSampleData) {
            vector<int> abunds = getAbunds(name);
            abund = accumulate(abunds.begin(), abunds.end(), 0);
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
vector<int> SeqAbundTable::getAbunds(int name) {
    vector<int> abunds;

    if (hasSampleData) {
        // all samples not just "good" ones
        vector<int> allAbunds(tableSamples.size(), 0);

        seqCount data = counts[name];

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
vector<string> SeqAbundTable::getSamples() {
    vector<string> samples;

    if (hasSampleData) {
        // want all "good" samples
        return select(sampleNames, tableSamples);
    }

    return samples;
}
/******************************************************************************/
vector<int> SeqAbundTable::getSampleTotals() {
    vector<int> totals;
    if (hasSampleData) {
        totals = select(sampleTotals, tableSamples);
    }
    return totals;
}
/******************************************************************************/
int SeqAbundTable::getNumSamples() {
    if (hasSampleData) {
        return numSamples;
    }
    return 0;
}
/******************************************************************************/
int SeqAbundTable::getNumTreatments() {
    if (hasTreatments) {
        return numTreatments;
    }
    return 0;
}
/******************************************************************************/
// vector containing total abundance for each sequence
vector<int> SeqAbundTable::getSeqsAbunds(vector<int> names) {

    vector<int> abunds(names.size(), 0);

    for (int i = 0; i < names.size(); i++) {
        abunds[i] = getAbund(names[i]);
    }
    return abunds;
}
/******************************************************************************/
Rcpp::DataFrame SeqAbundTable::getSequenceAbundanceTable(vector<string> outputNames,
                                          vector<int> names) {


    if (hasSampleData) {
        vector<string> ids;
        vector<int> abunds;
        vector<string> samples;
        vector<string> treaments;

        for (int i = 0; i < names.size(); i++) {
            string name = outputNames[i];

            seqCount data = counts[names[i]];

            // data -> sampleIndex(2,5), abunds(100, 50)
            // becomes abunds(0,0,100,0,0,50)
            for (int j = 0; j < data.sampleIndex.size(); j++) {
                //allAbunds[data.sampleIndex[i]] = data.abunds[i];
                // if this sample "good"
                if (tableSamples[data.sampleIndex[j]]) {
                    ids.push_back(name);
                    abunds.push_back(data.abunds[j]);
                    samples.push_back(sampleNames[data.sampleIndex[j]]);

                    if (hasTreatments) {
                        treaments.push_back(sampleTreatment[data.sampleIndex[j]]);
                    }
                }
            }
        }

        if (hasTreatments) {
            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named("id") = ids,
                Rcpp::_["abundance"] = abunds,
                Rcpp::_["sample"] = samples,
                Rcpp::_["treatment"] = treaments);
            return df;
        }else{
            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named("id") = ids,
                Rcpp::_["abundance"] = abunds,
                Rcpp::_["sample"] = samples);
            return df;
        }
    }else{
        // no sample information
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("id") = outputNames,
            Rcpp::_["abundance"] = getSeqsAbunds(names));
        return df;
    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
int SeqAbundTable::getSparseIndex(int name, int sample) {
    int index = -1;

    for (int i = 0; i < counts[name].sampleIndex.size(); i++) {
        if (counts[name].sampleIndex[i] == sample) { return i; }
    }

    return index;
}
/******************************************************************************/
int SeqAbundTable::getTotal(string sample) {
    if (sample != "") {
        auto it = sampleIndex.find(sample);

        // is this a sample in the table
        if (it != sampleIndex.end()) {
            // is it "good"
            if (tableSamples[it->second]) {
                return sampleTotals[it->second];
            }
        }
        RcppThread::Rcout << endl <<
            "[ERROR]: The dataset does not include sample: " +
            sample + ".\n\n";
        return 0;
    }
    return total;
}
/******************************************************************************/
vector<string> SeqAbundTable::getTreatments() {
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
vector<int> SeqAbundTable::getTreatmentTotals() {
    vector<int> totals;
    if (hasTreatments) {
        totals = select(treatmentTotals, tableTreatments);
    }
    return totals;
}
/******************************************************************************/
bool SeqAbundTable::hasSample(string sample, int name) {
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
int SeqAbundTable::removeSeq(int name) {

    int abund = 0;

    if (hasSampleData) {
        vector<int> abunds = getAbunds(name);
        abund = accumulate(abunds.begin(), abunds.end(), 0);

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
        abund = getAbund(name);
        total -= abund;
    }

    return abund;
}
/******************************************************************************/
void SeqAbundTable::updateTotals() {

    if (hasSampleData) {

    // does removing this sequence remove a sample
    for (int i = 0; i < sampleTotals.size(); i++) {
        if (sampleTotals[i] == 0) {
           tableSamples[i] = false;
            numSamples--;
        }
    }

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
            if (treatmentTotals[i] == 0) {
                tableTreatments[i] = false;
                numTreatments--;
            }
        }
    }
    }
}
/******************************************************************************/
