#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
SeqAbundTable::SeqAbundTable() {
    hasGroupData = false;
    total = 0;
    numGroups = 0;
}
/******************************************************************************/
SeqAbundTable::~SeqAbundTable() {}
/******************************************************************************/
void SeqAbundTable::clear() {
    hasGroupData = false;
    total = 0;
    numGroups = 0;

    counts.clear();
    groupIndex.clear();
    groupTotals.clear();
    tableGroups.clear();
}
/******************************************************************************/
// the names are the indexes in dataset
void SeqAbundTable::addSeqs(vector<int>& names) {

    if (hasGroupData) {
        RcppThread::Rcout << endl <<
            "[ERROR]: The dataset contains group information, you must provide by sample data.\n\n";
    }else{

        counts.resize(counts.size()+names.size());

        for (int i = 0; i < names.size(); i++) {
            seqCount thisSeq(0, 1);
            counts[names[i]] = thisSeq;
        }

        // make space for total abundance
        tableGroups.clear();
        tableGroups.push_back(true);
        numGroups = 1;

        total += names.size();
    }

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
// total abundance for sequence.
// If group provided, then abundance for sequence in sample
int SeqAbundTable::getAbund(int name, string group) {

    int abund = 0;

    if (group == "") {
        if (hasGroupData) {
            vector<int> abunds = getAbunds(name);
            abund = accumulate(abunds.begin(), abunds.end(), 0);
        }else{
            abund = counts[name].abunds[0];
        }
    }else if (hasGroup(group)) {
        int gIndex = groupIndex[group];

        // if the sequence does not have this group, -1 returned
        int thisGroupsIndex = getSparseIndex(name, gIndex);

        if (thisGroupsIndex != -1) {
            abund = counts[name].abunds[thisGroupsIndex];
        }
    }

    return abund;
}
/******************************************************************************/
int SeqAbundTable::getSparseIndex(int name, int group) {
    int index = -1;

    for (int i = 0; i < counts[name].sampleIndex.size(); i++) {
        if (counts[name].sampleIndex[i] == group) { return i; }
    }

    return index;
}
/******************************************************************************/
// abundances by sample, in the same order as the groups
vector<int> SeqAbundTable::getAbunds(int name) {
    vector<int> abunds;

    if (hasGroupData) {
        // all groups not just "good" ones
        vector<int> allAbunds(tableGroups.size(), 0);

        seqCount data = counts[name];

        // data -> sampleIndex(2,5), abunds(100, 50)
        // becomes abunds(0,0,100,0,0,50)
        for (int i = 0; i < data.sampleIndex.size(); i++) {
            allAbunds[data.sampleIndex[i]] = data.abunds[i];
        }

        // only include "good" groups
        abunds = select(allAbunds, tableGroups);
    }else{
        abunds.push_back(counts[name].abunds[0]);
    }

    return abunds;
}
/******************************************************************************/
bool SeqAbundTable::hasGroup(string group) {
    auto it = groupIndex.find(group);

    // is this a group in the table
    if (it != groupIndex.end()) {
        // is it "good"
        return tableGroups[it->second];
    }

    return false;
}
/******************************************************************************/
int SeqAbundTable::getTotal(string group) {
    if (group != "") {
        auto it = groupIndex.find(group);

        // is this a group in the table
        if (it != groupIndex.end()) {
            // is it "good"
            if (tableGroups[it->second]) {
                return groupTotals[it->second];
            }
        }
        RcppThread::Rcout << endl <<
            "[ERROR]: The dataset does not include group: " +
             group + ".\n\n";
        return 0;
    }
    return total;
}
/******************************************************************************/
