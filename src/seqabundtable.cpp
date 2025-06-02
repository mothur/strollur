#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
SeqAbundTable::SeqAbundTable() {
    hasGroupData = false;
    hasTreatments = false;
    total = 0;
    numGroups = 0;
    numTreatments = 0;
}
/******************************************************************************/
SeqAbundTable::~SeqAbundTable() {}
/******************************************************************************/
void SeqAbundTable::clear() {
    hasGroupData = false;
    hasTreatments = false;
    total = 0;
    numGroups = 0;
    numTreatments = 0;

    counts.clear();
    groupIndex.clear();
    groupTotals.clear();
    tableGroups.clear();

    treatmentIndex.clear();
    treatmentTotals.clear();
    tableTreatments.clear();
    groupTreatment.clear();
}
/******************************************************************************/
// the names are the indexes in dataset
void SeqAbundTable::addSeqs(vector<int>& names) {

    if (hasGroupData) {
        string message = "[ERROR]: The dataset contains group information, ";
        message += "you must provide by sample data.\n\n";
        RcppThread::Rcout << endl << message;
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

        tableTreatments.clear();
        tableTreatments.push_back(true);
        numTreatments = 1;

        total += names.size();
    }
}
/******************************************************************************/
// private function
void SeqAbundTable::addGroups(vector<string> groups) {
    // store groups alphabetically
    sort(groups.begin(), groups.end());

    // all groups start off as "good"
    tableGroups.resize(groups.size(), true);
    groupTotals.resize(groups.size(), 0);

    for (int i = 0; i < groups.size(); i++) {
        groupIndex[groups[i]] = i;
    }

    hasGroupData = true;
}
/******************************************************************************/
// private function
void SeqAbundTable::addTreatments(vector<string> treatments) {
    // store groups alphabetically
    sort(treatments.begin(), treatments.end());

    // all groups start off as "good"
    tableTreatments.resize(treatments.size(), true);
    treatmentTotals.resize(treatments.size(), 0);

    for (int i = 0; i < treatments.size(); i++) {
        treatmentIndex[treatments[i]] = i;
    }

    hasTreatments = true;
}
/******************************************************************************/
// names, abundances, groups (optional), treatment (optional)
void SeqAbundTable::assignSampleAbundance(vector<int> names,
                           vector<int> abunds,
                           vector<string> groups,
                           vector<string> treatments) {

    clear();

    if ((groups.size() == 0) && (treatments.size() == 0)) {
        hasGroupData = false;
        tableGroups.push_back(true);
        tableTreatments.push_back(true);
        numGroups = 1; // placeholder for sparse abundances
        numTreatments = 1;
    }else{
        vector<string> uniqueGroups = unique(groups);
        vector<string> uniqueTreatments = unique(treatments);

        if (uniqueGroups.size() > 1) {
            addGroups(uniqueGroups);
            numGroups = uniqueGroups.size();

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
        }else if ((uniqueGroups.size() == 1) && (uniqueGroups[0] != "")) {
            addGroups(uniqueGroups);
            numGroups = 1;
            numTreatments = 1;
        }else{
            // empty strings passed in for groups, ignore
            hasGroupData = false;
            hasTreatments = false;
            tableGroups.push_back(true);
            tableTreatments.push_back(true);
            numGroups = 1; // placeholder for sparse abundances
            numTreatments = 1;
        }
    }


    // update counts, groupTotals and total
    if (hasGroupData) {
        vector<int> uniqueNames = unique(names);

        // fill in long format
        counts.resize(uniqueNames.size());
        for (int i = 0; i < groups.size(); i++) {
            int index = groupIndex[groups[i]];
            counts[names[i]].sampleIndex.push_back(index);
            counts[names[i]].abunds.push_back(abunds[i]);
            groupTotals[index] += abunds[i];
        }

        // convert to sparse
        total = accumulate(groupTotals.begin(), groupTotals.end(), 0);
    }else{

        counts.resize(names.size());
        for (int i = 0; i < names.size(); i++) {
            seqCount thisSeq(0, abunds[i]);
            counts[names[i]] = thisSeq;
        }
        total = accumulate(abunds.begin(), abunds.end(), 0);
    }

    if (hasTreatments) {
        // create group to treatment map
        for (int i = 0; i < groups.size(); i++) {
            groupTreatment[groups[i]] = treatments[i];
        }

        // fill treatment totals
        for (auto it = groupIndex.begin(); it != groupIndex.end(); it++) {
            string treatment = groupTreatment[it->first];
            int tindex = treatmentIndex[treatment];

            treatmentTotals[tindex] += groupTotals[it->second];
        }
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
bool SeqAbundTable::hasGroup(string group, int name) {
    auto it = groupIndex.find(group);

    // is this a group in the table
    if (it != groupIndex.end()) {

        // no name provided
        if (name == -1) {
            // is it "good"
            return tableGroups[it->second];
        }else{
            // if the sequence does not have this group, -1 returned
            int thisGroupsIndex = getSparseIndex(name, it->second);

            if (thisGroupsIndex != -1) {
                return true;
            }
        }
    }

    return false;
}
/******************************************************************************/
vector<int> SeqAbundTable::getGroupTotals() {
    vector<int> totals;
    if (hasGroupData) {
        totals = select(groupTotals, tableGroups);
    }
    return totals;
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
int SeqAbundTable::getNumGroups() {
    if (hasGroupData) {
        return numGroups;
    }
    return 0;
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
vector<string> SeqAbundTable::getGroups(int name) {
    vector<string> groups;

    if (hasGroupData) {
        // want all "good" groups
        if (name == -1) {
            for (auto it = groupIndex.begin(); it != groupIndex.end(); it++) {
                if (tableGroups[it->second]) {
                    groups.push_back(it->first);
                }
            }
        }else{
            // want all "good" groups for sequence
            vector<int> thisSeqsGroupIndexes = counts[name].sampleIndex;

            for (auto it = groupIndex.begin(); it != groupIndex.end(); it++) {
                // if "good" group
                if (tableGroups[it->second]) {

                    auto itFind = find(thisSeqsGroupIndexes.begin(),
                                       thisSeqsGroupIndexes.end(), it->second);
                    if (itFind != thisSeqsGroupIndexes.end()) {
                        groups.push_back(it->first);
                    }
                }
            }

        }
    }

    return groups;
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
