#include "../inst/include/rdataset.h"

/******************************************************************************/
SeqAbundTable::SeqAbundTable() {
    hasGroupData = false;
    total = 0;
}
/******************************************************************************/
SeqAbundTable::~SeqAbundTable() {}
/******************************************************************************/
void SeqAbundTable::clear() {
    hasGroupData = false;
    total = 0;

    counts.clear();
    groupIndex.clear();
    groupTotals.clear();
    tableGroups.clear();
}
/******************************************************************************/
void SeqAbundTable::addSeqs(vector<int>& names) {

    if (hasGroupData) {

    }else{

    }
    //TODO
}
/******************************************************************************/
