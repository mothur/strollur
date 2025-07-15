
#include "../inst/include/rdataset.h"
#include "seqreport.h"
#include "summary.h"
#include "dataset.h"
#include "utils.h"
#include "phylotree.h"

/******************************************************************************/
Dataset::Dataset(string n, int proc) : datasetName(n) {
    isAligned = false;
    hasContigsData = false;
    hasAlignData = false;
    hasOtuData = false;
    hasSequenceData = false;
    hasSequenceTaxonomy = false;
    hasOtuTaxonomy = false;
    runClassifyOtu = true;
    numSamples = 0;
    numTreatments = 0;
    numOtus = 0;
    label = "";
    numUnique = 0;
    uniqueBad = 0;
    alignmentLength = 0;
    processors = proc;
    count = new AbundTable();
    otuTable = nullptr;
}
/******************************************************************************/
Dataset::~Dataset() {
    delete count;
    if (otuTable != nullptr) { delete otuTable; }
}
/******************************************************************************/
void Dataset::clear() {
    isAligned = false;
    hasContigsData = false;
    hasAlignData = false;
    hasSequenceData = false;
    hasSequenceTaxonomy = false;
    hasOtuTaxonomy = false;
    runClassifyOtu = true;
    numSamples = 0;
    numTreatments = 0;
    numOtus = 0;
    label = "";
    numUnique = 0;
    uniqueBad = 0;
    alignmentLength = 0;

    // sequence data
    names.clear();
    seqs.clear();
    comments.clear();
    trashCodes.clear();

    // contigs report
    olengths.clear();
    ostarts.clear();
    oends.clear();
    mismatches.clear();
    ee.clear();

    // sequence summary data
    starts.clear();
    ends.clear();
    lengths.clear();
    ambigs.clear();
    polymers.clear();
    numns.clear();

    // alignment report
    searchScore.clear();
    simScore.clear();
    longestInsert.clear();

    // sequence taxonomy assignments
    taxonomies.clear();

    // maps sequence name to index in vectors
    seqIndex.clear();
    tableSeqs.clear();

    badAccnos.clear();

    count->clear();

    // if you have an otuTable then otuTable.clear();
    if (hasOtuData) {
        otuTable->clear();
        otuTable->setLabel(label);
        seqOtus.clear();
        list.clear();
        hasOtuData = false;
    }
}
/******************************************************************************/
Rcpp::List Dataset::exportDataset(){
    Rcpp::List result;

    // sequence data.frame
    // name, seqs, comments

    // sequence report data.frame
    // starts, ends, lengths, ambigs, polymers, numns

    // sequence trashCodes
    // name, trashCode

    // sequence abundance table
    // id, abundance, sample, treatment

    // sequence otu table
    // id, abundance, sample, seq_id

    // sequence taxonomy table

    // otu taxonomy table

    // otu trashCodes

    // align report

    // contigs report

    // metadata




    return result;
}
/******************************************************************************/
void Dataset::addSequences(vector<string> n, vector<string> s, vector<string> c) {

    // must provide the same number of names and seqs
    if (s.size() == 0) {
        s.resize(names.size(), "");
    }

    // add to seqIndex
    int numSeqs = names.size();
    vector<int> countNames;
    for (int i = 0; i < n.size(); i++) {
        countNames.push_back(numSeqs);
        seqIndex[n[i]] = numSeqs;
        numSeqs++;
    }

    // add to count
    count->add(countNames);
    countNames.clear();

    // add to names
    names.insert(names.end(), n.begin(), n.end());
    // add to seqs
    seqs.insert(seqs.end(), s.begin(), s.end());

    // innocent
    vector<bool> ts(n.size(), true);
    tableSeqs.insert(tableSeqs.end(), ts.begin(), ts.end());

    if (c.size() == 0) {
        c.resize(names.size(), "");
    }
    // add to comments
    comments.insert(comments.end(), c.begin(), c.end());

    // blank trash code because we assume seqs are "good"
    vector<string> t(names.size(), "");
    trashCodes.insert(trashCodes.end(), t.begin(), t.end());

   // add calcs for starts, ends, lengths, ambigs, polymers, numns
    SeqReport report;
    report.addReports(s, starts, ends, lengths, ambigs, polymers, numns);

    // set isAligned and aligned length
    getAlignedLength();

    // add to "good" sequence count - giving preference to otuTable
    numUnique += names.size();

    hasSequenceData = true;
}
/******************************************************************************/
void Dataset::assignOtus(vector<string> otuIds,
                        vector<int> abunds,
                        vector<string> samples,
                        vector<string> seqIds, string type) {


    // sanity checks - R6 object passes blank strings if seqIds == NULL
    bool useSeqIds = false;
    if (!seqIds.empty()) {
        useSeqIds = true;
        string id = "";
        if (allIdentical(seqIds, id)) {
            if (id == "") { seqIds.clear(); useSeqIds = false; }
        }
    }

    // sanity checks - R6 object passes blank strings if samples == NULL
    if (!samples.empty()) {
        string id = "";
        if (allIdentical(samples, id)) {
            if (id == "") { samples.clear(); }
        }
    }

    // sanity checks - R6 object passes all 0's if abundance == NULL
    if (!abunds.empty()) {
        int id = 0;
        if (allIdentical(abunds, id)) {
            if (id == 0) { abunds.clear();  }
        }
    }


    if (!useSeqIds && (otuIds.size() != abunds.size())) {
        string message = "[ERROR]: Size mismatch. otu_ids and abunds must be";
        message += " the same size.";

        throw Rcpp::exception(message.c_str());
    }

    // new table
    if (otuTable == nullptr) {
        otuTable = new OtuTable("otu");
    }else{
        otuTable->clear();
        otuTable->setLabel("otu");
    }

    if (useSeqIds) {

        // if you don't have sequence data, but are adding seq_names for otus,
        // add seqs, add abundances if provided
        if (!hasSequenceData) {

            vector<string> uniqueSeqIds = unique(seqIds);
            addSequences(uniqueSeqIds);

            // if the abundances are provided then use them
            if (!abunds.empty()) {
                assignSequenceAbundance(seqIds, abunds, samples);
            }
        }else {
            if (!abunds.empty() || !samples.empty()) {
                string message = "[WARNING]: Ignoring sample abundances, using";
                message += " sequence sample abundances already assigned.";
                RcppThread::Rcout << endl << message << endl;
            }
        }
        abunds.clear();
        samples.clear();

        // if the user provides seqIds,
        // we need to calculate the otu abunds by sample
        vector<string> allSamples = count->getSamples();

        // does the sequence abundance include samples
        bool hasSamples = false;
        if (!allSamples.empty()) {
            hasSamples = true;
        }

        // otuName -> (sampleName -> abundance)
        map<string, map<string, int>> otuAbunds;

        // maps seqName to otuName
        map<string, string> seqOtu;
        seqOtus.resize(names.size(), "");

        for (int i = 0; i < seqIds.size(); i++) {

            string otuName = otuIds[i];

            auto itSeq = seqOtu.find(seqIds[i]);
            bool firstTimeSeq = false;

            // first time we are seeing this seqName, add seq to otu
            if (itSeq == seqOtu.end()) {
                firstTimeSeq = true;
                seqOtu[seqIds[i]] = otuName;

                auto seqI = seqIndex.find(seqIds[i]);

                // this seq is in the dataset
                if (seqI != seqIndex.end()) {

                    seqOtus[seqI->second] = otuName;
                    auto itList = list.find(otuName);

                    // add to list - (list file)
                    if (itList != list.end()) {

                        (itList->second).push_back(seqI->second);
                    }else{
                        vector<int> temp; temp.push_back(seqI->second);
                        list[otuName] = temp;
                    }
                }else{
                    continue;
                }
            }

            auto it = otuAbunds.find(otuName);

            // first time seeing this otu, initialize sample counts
            if (it == otuAbunds.end()) {
                map<string, int> thisSeqsSampleAbunds;

                // inputs from dataset
                int thisIndex = seqIndex[seqIds[i]];
                vector<int> abunds = count->getAbundances(thisIndex);

                // count has samples
                if (hasSamples) {
                    for (int j = 0; j < abunds.size(); j++) {
                        if (abunds[j] != 0) {
                            thisSeqsSampleAbunds[allSamples[j]] = abunds[j];
                        }
                    }
                }else{
                    thisSeqsSampleAbunds["total"] = abunds[0];
                }

                otuAbunds[otuName] = thisSeqsSampleAbunds;
            }else{
                // inputs from dataset
                if (firstTimeSeq) {

                    int thisIndex = seqIndex[seqIds[i]];
                    vector<int> abunds = count->getAbundances(thisIndex);

                    // count has samples
                    if (hasSamples) {
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
        otuIds.clear();
        if (hasSamples) {

            for (auto it = otuAbunds.begin(); it != otuAbunds.end(); it++) {

                map<string, int> abundances = it->second;
                for (auto itSamples = abundances.begin();
                     itSamples != abundances.end(); itSamples++) {
                    otuIds.push_back(it->first);
                    samples.push_back(itSamples->first);
                    abunds.push_back(itSamples->second);
                }
            }

            otuTable->assignAbundance(otuIds, abunds, samples);

            // inputs from dataset
            if (hasSamples) {
                // if there are treatments add them to otu
                if (count->getNumTreatments() != 0) {
                    map<string, string> sampleTreatments =
                        count->getSampleTreatmentAssignments();

                    vector<string> samples(sampleTreatments.size());
                    vector<string> treatments(sampleTreatments.size());

                    int index = 0;
                    for (auto itTreatments = sampleTreatments.begin();
                         itTreatments != sampleTreatments.end(); itTreatments++) {
                        samples[index] = itTreatments->first;
                        treatments[index] = itTreatments->second;
                        index++;
                    }

                    otuTable->assignTreatments(samples, treatments);
                }
            }

        }else{

            for (auto it = otuAbunds.begin(); it != otuAbunds.end(); it++) {
                otuIds.push_back(it->first);
                abunds.push_back(it->second["total"]);
            }

            otuTable->assignAbundance(otuIds, abunds);
        }

        hasOtuData = true;

        // if the sequences have been assigned taxonomy,
        // assign the otu taxonomy
        if (hasSequenceTaxonomy) { classifyOtus(); }
    }else{
        // add otus - shared
        otuTable->assignAbundance(otuIds, abunds, samples);
    }

    hasOtuData = true;
    numOtus = otuTable->numOtus;
    label = otuTable->label;
    numSamples = otuTable->getNumSamples();
    numTreatments = otuTable->getNumTreatments();
}
/******************************************************************************/
void Dataset::assignOtuTaxonomy(vector<string> otuIds, vector<string> taxs) {

    if (hasOtuData) {
        if (otuIds.size() != taxs.size()) {
            string message = "[ERROR]: Size mismatch. otu_ids and taxonomies ";
            message += "must be the same size.";
            throw Rcpp::exception(message.c_str());
        }

        otuTable->assignTaxonomy(otuIds, taxs);
        hasOtuTaxonomy = true;
        runClassifyOtu = false;
    }
}
/******************************************************************************/
void Dataset::assignTreatments(vector<string> samples, vector<string> treatments) {
    count->assignTreatments(samples, treatments);
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
    if (hasOtuData) {
        otuTable->assignTreatments(samples, treatments);
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
    }
}
/******************************************************************************/
// names, abundances, samples(optional), treatments(optional)
// assumes same size
void Dataset::assignSequenceAbundance(vector<string> ids,
                           vector<int> abunds,
                           vector<string> samples,
                           vector<string> treatments) {

    vector<string> uniqueNames = unique(ids);

    // are there assignments for all seqs in the dataset
    if (uniqueNames.size() != numUnique){
        string message = "[ERROR]: The dataset contains ";
        message += toString(numUnique) + " sequences, but you assigned ";
        message += toString(uniqueNames.size()) + " sequences. All sequences ";
        message += "in the dataset must be assigned abundances.\n\n";
        throw Rcpp::exception(message.c_str());
    }

    vector<int> idIndexes = getIndexes(ids);

    count->assignAbundance(idIndexes, abunds, samples, treatments);

    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
}
/******************************************************************************/
void Dataset::assignSequenceTaxonomy(vector<string> n, vector<string> t){

    if (!hasSequenceData) {
        addSequences(n);
    }

    // allocate space
    if (taxonomies.size() != names.size()) {
        taxonomies.resize(names.size(), "");
    }

    if (n.size() != t.size()) {
        string message = "[ERROR]: Size mismatch. ids and taxonomies must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    for (int i = 0; i < n.size(); i++) {

        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {

            int index = it->second;

            // update taxonomy
            taxonomies[index] = t[i];
        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    hasSequenceTaxonomy = true;

    if (hasOtuData) { classifyOtus(); }
}
/******************************************************************************/
// align_seqs will create searchScores, simScores and longestInserts
void Dataset::addAlignReport(vector<string>& n, vector<double>& ss,
                    vector<double>& sims,
                    vector<int>& li) {

    if (names.size() == 0) {
        addSequences(n);
    }

    set<int> sizes;
    sizes.insert(n.size());
    sizes.insert(ss.size());
    sizes.insert(sims.size());
    sizes.insert(li.size());
    sizes.insert(numUnique);

    if (sizes.size() > 1) {
        string message = "[WARNING]: You must provide align report info for ";
        message += "each sequence in the dataset, ignoring align report data.";
        RcppThread::Rcout << endl << message << endl;
    }else {
        // create space
        searchScore.resize(names.size(), 0);
        simScore.resize(names.size(), 0);
        longestInsert.resize(names.size(), 0);
        hasAlignData = true;

        for (int i = 0; i < n.size(); i++) {

            auto it = seqIndex.find(n[i]);
            if (it != seqIndex.end()) {
                int index = it->second;
                searchScore[index] = ss[index];
                simScore[index] = sims[index];
                longestInsert[index] = li[index];
            }else{
                string message = "[ERROR]: The dataset does not contain a ";
                message += "sequence named " + n[i] + ", quitting.\n";
                RcppThread::Rcerr << endl << message << endl;
                throw Rcpp::exception(message.c_str());
            }
        }
    }

}
/******************************************************************************/
// make_contigs will create overlapLengths, overlapStarts, overlapEnds,
// mismatches, and expectedErrors
void Dataset::addContigsReport(vector<string>& n, vector<int>& ol,
                      vector<int>& os,
                      vector<int>& oe,
                      vector<int>& m,
                      vector<double>& e) {

    if (names.size() == 0) {
        addSequences(n);
    }

    set<int> sizes;
    sizes.insert(n.size());
    sizes.insert(ol.size());
    sizes.insert(os.size());
    sizes.insert(oe.size());
    sizes.insert(m.size());
    sizes.insert(e.size());
    sizes.insert(numUnique);

    if (sizes.size() > 1) {
        string message = "[WARNING]: You must provide contigs report info for ";
        message += "each sequence in the dataset, ignoring contigs report ";
        message += "data.";
        RcppThread::Rcout << endl << message << endl;
    }else {
        hasContigsData = true;
        // create space
        olengths.resize(names.size(), 0);
        ostarts.resize(names.size(), 0);
        oends.resize(names.size(), 0);
        mismatches.resize(names.size(), 0);
        ee.resize(names.size(), 0);

        for (int i = 0; i < n.size(); i++) {

            auto it = seqIndex.find(n[i]);
            if (it != seqIndex.end()) {
                int index = it->second;
                olengths[index] = ol[index];
                ostarts[index] = os[index];
                oends[index] = oe[index];
                mismatches[index] = m[index];
                ee[index] = e[index];
            }else{
                string message = "[ERROR]: The dataset does not contain a ";
                message += "sequence named " + n[i] + ", quitting.\n";
                RcppThread::Rcerr << endl << message << endl;
                throw Rcpp::exception(message.c_str());
            }
        }
    }
}
/******************************************************************************/
string Dataset::classifyOtu(string otuName){
    string otuTax = "";

    auto it = list.find(otuName);

    // valid otu
    if (it != list.end()) {

        vector<int> otuSeqs = it->second;

        PhyloTree phylo;
        int size = 0;

        // add all "good" seqs in otu to phylotree
        for (int index : otuSeqs) {

            // is a "good" seq
            if (tableSeqs[index]) {
                int seqAbund = count->getAbundance(index);
                size += seqAbund;
                phylo.addSeqToTree(taxonomies[index], seqAbund);
            }
        }

        TaxNode currentNode = phylo.getRoot();

        while (currentNode.total != 0) {

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
                otuTax += bestChild.name + "(" + toString(consensusConfidence) + ");";
            }
            //move down a level
            currentNode = bestChild;
        }
    }

    if (otuTax == "") {  otuTax = "unknown;";  }

    return otuTax;
}
/******************************************************************************/
void Dataset::classifyOtus() {
    if (hasOtuData) {
        // get names of "good" otus
        vector<string> otuNames = otuTable->getOtuIds();
        vector<string> otuClassifications(otuNames.size(), "");
        for (int i = 0; i < otuNames.size(); i++) {
            otuClassifications[i] = classifyOtu(otuNames[i]);
        }
        otuTable->assignTaxonomy(otuNames, otuClassifications);

        hasOtuTaxonomy = true;
        numOtus = otuTable->numOtus;
        label = otuTable->label;
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
        runClassifyOtu = false;
    }
}
/******************************************************************************/
/*
 id         level  taxon                 confidence
 seq1     1       Bacteria             100.0
 seq1     2      "Acidobacteria"  99.8
 seq1     3      Holophagae        99.8
 seq1     4      Holophagales      95.0
 seq1     5      Holophagaceae   90.0
 seq1     6      Holophaga           87.0
 seq2 ...
 */
Rcpp::DataFrame Dataset::fillTaxReport(string mode) {

    vector<string> ids, taxes;

    if (mode == "otu") {
        ids = otuTable->getOtuIds();
        taxes = otuTable->getTaxonomies();
    }else {
        ids = getNames();
        taxes = select(taxonomies, tableSeqs);
    }

    Utils util;
    int maxLevel = 1;

    vector<vector<string> > taxons(taxes.size());
    vector<vector<int> > confidences(taxes.size());
    bool hasConfidences = true;

    // split taxonomy and confidences by level
    for (int i = 0; i < taxes.size(); i++) {
        int numLevels = split(taxes[i], ';',
                              back_inserter(taxons[i]));

        if (numLevels > maxLevel) { maxLevel = numLevels; }

        confidences[i] = util.removeConfidences(taxons[i]);

        if (hasConfidences) {
            if (sum(confidences[i]) == 0) {
                hasConfidences = false;
            }
        }
    }

    taxes.clear();
    vector<string> dfIds(taxons.size()*maxLevel);
    vector<string> dfTaxs(taxons.size()*maxLevel);
    vector<int> levels(taxons.size()*maxLevel);
    vector<int> dfConfidences(taxons.size()*maxLevel);

    for (int i = 0; i < taxons.size(); i++) {
        // extend taxonomies to the same level by adding unclassifieds
        util.addUnclassifieds(taxons[i], confidences[i], maxLevel);

        for (int j = 0; j < maxLevel; j++) {
            dfIds[i*maxLevel+j] = ids[i];
            dfTaxs[i*maxLevel+j] = taxons[i][j];
            levels[i*maxLevel+j] = j+1;
            dfConfidences[i*maxLevel+j] = confidences[i][j];
        }
    }

    taxons.clear();
    confidences.clear();


    if (!hasConfidences) {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("id") = dfIds,
            Rcpp::_["level"] = levels,
            Rcpp::_["taxon"] = dfTaxs);
        return df;
    }

    Rcpp::DataFrame df = Rcpp::DataFrame::create(
        Rcpp::Named("id") = dfIds,
        Rcpp::_["level"] = levels,
        Rcpp::_["taxon"] = dfTaxs,
        Rcpp::_["confidence"] = dfConfidences);

    return df;
}
/******************************************************************************/
int Dataset::getAbundance(string name){
    int abund = 0;
    auto it = seqIndex.find(name);
    if (it != seqIndex.end()) {
        if (tableSeqs[it->second]) {
            abund = count->getAbundance(it->second, "");
        }
    }
    return abund;
}
/******************************************************************************/
vector<int> Dataset::getAbundances(string name){
    vector<int> abunds;
    auto it = seqIndex.find(name);
    if (it != seqIndex.end()) {
        if (tableSeqs[it->second]) {
            abunds = count->getAbundances(it->second);
        }
    }
    return abunds;
}
/******************************************************************************/
int Dataset::getAlignedLength() {
    set<int> seqLengths;
    for (int i = 0; i < seqs.size(); i++) {
        // is this a "good" seq
        if (tableSeqs[i]) {
            seqLengths.insert(seqs[i].length());
        }
    }

    if (seqLengths.size() == 1) {
        isAligned = true;
        alignmentLength = *(seqLengths.begin());
    }else{
        isAligned = false;
        alignmentLength = -1;
    }
    return alignmentLength;
}
/******************************************************************************/
// search_score, sim_score, longest_insert
Rcpp::DataFrame Dataset::getAlignReport(){

    if (hasAlignData) {
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("id") = select(names, tableSeqs),
            Rcpp::_["search_score"] = select(searchScore, tableSeqs),
            Rcpp::_["sim_score"] = select(simScore, tableSeqs),
            Rcpp::_["longest_insert"] = select(longestInsert, tableSeqs));

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// contigs sumary data: olengths, ostarts, oends, mismatches, ee
// report[0] = length, report[1] = overlap_length, report[2] = overlap_start,
// report[3] = overlap_end, report[4] = mismatches, report[5] = num_ns,
// report[6] = ee
Rcpp::DataFrame Dataset::getContigsReport(){

    if (hasContigsData){
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("id") = select(names, tableSeqs),
            Rcpp::_["length"] = select(lengths, tableSeqs),
            Rcpp::_["overlap_length"] = select(olengths, tableSeqs),
            Rcpp::_["overlap_start"] = select(ostarts, tableSeqs),
            Rcpp::_["overlap_end"] = select(oends, tableSeqs),
            Rcpp::_["mismatches"] = select(mismatches, tableSeqs),
            Rcpp::_["num_n"] = select(numns, tableSeqs),
            Rcpp::_["ee"] = select(ee, tableSeqs));

        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// returns indexes of "good" seqs in table
vector<int> Dataset::getIncludedNamesIndexes() {
    vector<int> included;
    for (int i = 0; i < tableSeqs.size(); i++) {
        if (tableSeqs[i]) {
            included.push_back(i);
        }
    }
    return included;
}
/******************************************************************************/
vector<int> Dataset::getIndexes(vector<string>& ids) {
    vector<int> indexes(ids.size(), -1);
    map<string, int>::iterator it;

    for (int i = 0; i < ids.size(); i++) {

        it = seqIndex.find(ids[i]);
        if (it != seqIndex.end()) {
            indexes[i] = it->second;
        }else{
            string message = "[ERROR]: The dataset does not contain a ";
            message += "sequence named " + ids[i] + ".\n";
            RcppThread::Rcerr << endl << message << endl;
            throw Rcpp::exception(message.c_str());
        }
    }
    return indexes;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getList() {
    if (hasOtuData && list.size() != 0) {
        vector<string> otuNames = otuTable->getOtuIds();

        vector<string> ids, seqids;
        for (int i = 0; i < otuNames.size(); i++) {
            auto it = list.find(otuNames[i]);

            if (it != list.end()) {
                // add each sequence name
                for (int j = 0; j < it->second.size(); j++) {
                    int index = it->second[j];
                    if (tableSeqs[index]) {
                        seqids.push_back(names[index]);
                        ids.push_back(otuNames[i]);
                    }

                }
            }
        }

        string tag = otuTable->label + "_id";
        Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named(tag.c_str()) = ids,
                Rcpp::_["seq_id"] = seqids);
        return df;
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
vector<string> Dataset::getListVector() {
    if (hasOtuData && list.size() != 0) {
        vector<string> otuNames = otuTable->getOtuIds();

        vector<string> otus;
        for (int i = 0; i < otuNames.size(); i++) {
            string otu = getOtu(otuNames[i]);
            // is this a "good" otu
            if (otu != "") {
                otus.push_back(otu);
            }
        }
        return otus;
    }

    return nullVector;
}
/******************************************************************************/
vector<string> Dataset::getNames(string sample){
    vector<string> included;

    // get all "good" names in dataset
    if (sample == "")  {
        for (int i = 0; i < tableSeqs.size(); i++) {
            if (tableSeqs[i]) {
                included.push_back(names[i]);
            }
        }
    // get all "good" names in specific sample
    }else {
        if (count->hasSample(sample)) {
            // all seqs
            for (int i = 0; i < tableSeqs.size(); i++) {
                // if "good" seq
                if (tableSeqs[i]) {
                    // if this sequence is in sample
                    if (count->hasSample(sample, i)) {
                        included.push_back(names[i]);
                    }
                }
            }
        }
    }

    return included;
}
/******************************************************************************/
vector<vector<string> > Dataset::getNamesBySample(vector<string> samples){
    vector<vector<string> > result;

    for (int i = 0; i < samples.size(); i++) {
        result.push_back(getNames(samples[i]));
    }

    return result;
}
/******************************************************************************/
// total abundance for a given outID
int Dataset::getOtuAbundance(string otuId) {
    if (hasOtuData) {
        return otuTable->getAbundance(otuId, "");
    }
    return 0;
}
/******************************************************************************/
// abundances for given otuID broken down by sample
vector<int> Dataset::getOtuAbundances(string otuId) {
    if (hasOtuData) {
        return otuTable->getAbundances(otuId);
    }
    return nullIntVector;
}
/******************************************************************************/
// string containing sequence names for given otuID
string Dataset::getOtu(string otuId) {
    if (hasOtuData) {
        // provided sequence otu data
        if (list.size() != 0) {
            // is this a good otu
            if (otuTable->hasId(otuId)) {

                auto it = list.find(otuId);

                vector<string> otuSeqNames;
                vector<int> otuSeqIndexes = it->second;
                for (int j = 0; j < otuSeqIndexes.size(); j++) {
                    if (tableSeqs[otuSeqIndexes[j]]) {
                        otuSeqNames.push_back(names[otuSeqIndexes[j]]);
                    }
                }
                return toString(otuSeqNames, ',');
            }
        }
    }
    return "";
}
/******************************************************************************/
vector<string> Dataset::getOtuIds() {
    if (hasOtuData) {
        return otuTable->getOtuIds();
    }
    return nullVector;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getOtuTaxonomyReport() {
    if (hasOtuTaxonomy){
        // this is set if you have removed sequences after classifying
        if (runClassifyOtu) {
            classifyOtus();
        }
        return (fillTaxReport("otu"));
    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getRAbund() {
    if (hasOtuData) {
        return otuTable->getRAbund();
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// sample functions
vector<string> Dataset::getSamples(){
    if (hasOtuData) {
       return otuTable->getSamples();
    }
    return count->getSamples();
}
/******************************************************************************/
vector<int> Dataset::getSampleTotals(){
    if (hasOtuData) {
        return otuTable->getSampleTotals();
    }
    return count->getSampleTotals();
}
/******************************************************************************/
// id, trashCode
Rcpp::DataFrame Dataset::getScrapReport(string mode) {

    if (mode == "sequence") {
        if (badAccnos.size() != 0) {
            vector<string> badNames(uniqueBad, "");
            vector<string> badCodes(uniqueBad, "");

            int next = 0;
            for (int i = 0; i < trashCodes.size(); i++) {
                if (trashCodes[i] != "") {
                    badNames[next] = names[i];
                    // remove last comma
                    trashCodes[i].pop_back();
                    badCodes[next] = trashCodes[i];
                    next++;
                }
            }

            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::Named("id") = badNames,
                Rcpp::_["trash_code"] = badCodes);
            return df;
        }
    }else if (mode == "otu") {
        if (hasOtuData) {
            return otuTable->getScrapReport();
        }
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
// trashCode, uniqueCount, totalCount
Rcpp::List Dataset::getScrapSummary() {

    Rcpp::List list = Rcpp::List::create();
    vector<string> listNames;

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

        Rcpp::DataFrame df = Rcpp::DataFrame::create(
            Rcpp::Named("trash_code") = codes,
            Rcpp::_["unique_count"] = uniqueCounts,
            Rcpp::_["total_count"] = totalCounts);

        list.push_back(df);
        listNames.push_back("sequence_scrap_summary");
    }

    if (hasOtuData) {
        list.push_back(otuTable->getScrapSummary());
        listNames.push_back("otu_scrap_summary");
    }

    list.attr("names") = listNames;

    return list;
}
/******************************************************************************/
// ids, abundances, sample(optional), treatment(optional)
// This table represents mothur's count and design files.
Rcpp::DataFrame Dataset::getSequenceAbundanceTable() {
    return count->getAbundanceTable(select(names, tableSeqs),
                                            getIncludedNamesIndexes());
}
/******************************************************************************/
// total abundance for each sequence
vector<int> Dataset::getSequenceAbundances(){
    return count->getTotalAbundances(getIncludedNamesIndexes());
}
/******************************************************************************/
vector<vector<int>> Dataset::getSeqsAbundsBySample(){
    return count->getAbundances(getIncludedNamesIndexes());
}
/******************************************************************************/
vector<string> Dataset::getSequences(string sample){
    vector<string> included;

    if (sample == "") {
        for (int i = 0; i < tableSeqs.size(); i++) {
            if (tableSeqs[i]) {
                included.push_back(seqs[i]);
            }
        }
    }else{
        if (count->hasSample(sample)) {
            // all seqs
            for (int i = 0; i < tableSeqs.size(); i++) {
                // if "good" seq
                if (tableSeqs[i]) {
                    // if this sequence is in sample
                    if (count->hasSample(sample, i)) {
                        included.push_back(seqs[i]);
                    }
                }
            }
        }
    }

    return included;
}
/******************************************************************************/
vector<vector<string> > Dataset::getSequencesBySample(vector<string> samples){
    vector<vector<string> > result;

    for (int i = 0; i < samples.size(); i++) {
        result.push_back(getSequences(samples[i]));
    }

    return result;
}
/******************************************************************************/
// fasta summary data: starts, ends, lengths, ambigs, polymers, numns
Rcpp::DataFrame Dataset::getSequenceReport(){

    Rcpp::DataFrame df = Rcpp::DataFrame::create(
        Rcpp::Named("id") = select(names, tableSeqs),
        Rcpp::_["start"] = select(starts, tableSeqs),
        Rcpp::_["end"] = select(ends, tableSeqs),
        Rcpp::_["length"] = select(lengths, tableSeqs),
        Rcpp::_["ambig"] = select(ambigs, tableSeqs),
        Rcpp::_["longest_homopolymer"] = select(polymers, tableSeqs),
        Rcpp::_["num_n"] = select(numns, tableSeqs));

    return df;
}
/******************************************************************************/
Rcpp::List Dataset::getSequenceSummary() {

    Rcpp::List result = Rcpp::List::create();

    if (!hasSeqs()) {
        return result;
    }

    vector<string> result_names;

    Summary* summary = new Summary(processors);

    vector<vector<int> > report;

    report.push_back(select(starts, tableSeqs));
    report.push_back(select(ends, tableSeqs));
    report.push_back(select(lengths, tableSeqs));
    report.push_back(select(ambigs, tableSeqs));
    report.push_back(select(polymers, tableSeqs));
    report.push_back(select(numns, tableSeqs));

    Rcpp::DataFrame seqResults = summary->summarizeFasta(
        report, count->getTotalAbundances(getIncludedNamesIndexes()));

    result.push_back(seqResults);
    result_names.push_back("sequence_summary");

    if (hasContigsData) {

        vector<vector<int> > report;
        report.push_back(select(lengths, tableSeqs));
        report.push_back(select(olengths, tableSeqs));
        report.push_back(select(ostarts, tableSeqs));
        report.push_back(select(oends, tableSeqs));
        report.push_back(select(mismatches, tableSeqs));
        report.push_back(select(numns, tableSeqs));

        Rcpp::DataFrame contigsResults = summary->summarizeContigs(
            report, count->getTotalAbundances(getIncludedNamesIndexes()));

        result.push_back(contigsResults);
        result_names.push_back("contigs_summary");
    }

    if (hasAlignData) {
        vector<vector<float> > report;
        report.push_back(select(searchScore, tableSeqs));
        report.push_back(select(simScore, tableSeqs));

        Rcpp::DataFrame alignResults = summary->summarizeAlign(
            report, select(longestInsert, tableSeqs),
            count->getTotalAbundances(getIncludedNamesIndexes()));

        result.push_back(alignResults);
        result_names.push_back("align_summary");
    }

    delete summary;

    Rcpp::DataFrame scrap = getScrapSummary();
    if (scrap.size() != 0) {
        result.push_back(scrap);
        result_names.push_back("scrap_summary");
    }

    result.attr("names") = result_names;
    return result;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getSequenceTaxonomyReport() {
    if (hasSequenceTaxonomy){
        return (fillTaxReport("sequence"));
    }

    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
Rcpp::DataFrame Dataset::getShared() {
    if (hasOtuData) {
        return otuTable->getShared();
    }
    Rcpp::DataFrame empty = Rcpp::DataFrame::create();
    return empty;
}
/******************************************************************************/
vector<string> Dataset::getTreatments(){
    if (hasOtuData) {
        return otuTable->getTreatments();
    }
    return count->getTreatments();
}
/******************************************************************************/
vector<int> Dataset::getTreatmentTotals(){
    if (hasOtuData) {
        return otuTable->getTreatmentTotals();
    }
    return count->getTreatmentTotals();
}
/******************************************************************************/
long long Dataset::getTotal(string sample){
    if (hasOtuData) {
        return otuTable->getTotal(sample);
    }
    return count->getTotal(sample);
}
/******************************************************************************/
long long Dataset::getUniqueTotal(string sample){
    if (sample == "") {
        return numUnique;
    }
    return getNames(sample).size();
}
/******************************************************************************/
bool Dataset::hasSample(string sample){
    return count->hasSample(sample);
}
/******************************************************************************/
bool Dataset::hasSeqs() {
    if (seqs.size() == 0) {
        return false;
    }

    string id;
    if (allIdentical(seqs, id)) {
        if (id == "") { return false; }
    }
    return true;
}
/******************************************************************************/
void Dataset::mergeSequences(vector<string> ids, string reason){
    if (ids.size() != 1) {

        vector<int> indexes = getIndexes(ids);

        // sanity check: if you have assigned otus, make sure the sequences
        // are in the same otu
        if (hasOtuData) {
            string otu = seqOtus[indexes[0]];
            for (int i = 1; i < indexes.size(); i++) {
                if (seqOtus[indexes[i]] != otu) {
                    string message = "[ERROR]: can not merge sequences assigned";
                    message += " to different otus.";
                    RcppThread::Rcerr << endl << message << endl;
                    throw Rcpp::exception(message.c_str());
                }
            }
        }

        count->merge(indexes);

        for (int i = 1; i < indexes.size(); i++) {
            // no need to update the sample and treatment counts
            removeSequence(indexes[i], reason, false);
        }
    }
}
/******************************************************************************/
void Dataset::mergeOtus(vector<string> ids, string reason){
    if (ids.size() != 1) {
        if (hasOtuData) {

            // if we have assigned sequences to otus
            if ((list.size() != 0) && (ids.size() != 1)) {

                string mergeOtu = "";
                vector<int> mergedSeqs;
                for (int i = 0; i < ids.size(); i++) {
                    // is this a good otu
                    if (otuTable->hasId(ids[i])) {
                        if (mergeOtu == "") {
                            mergeOtu = ids[i];
                            mergedSeqs = list[ids[i]];
                        }else {
                            vector<int> thisOtusSeqs = list[ids[i]];
                            mergedSeqs.insert(mergedSeqs.end(),
                                              thisOtusSeqs.begin(),
                                              thisOtusSeqs.end());
                            for (int j = 0; j < thisOtusSeqs.size(); j++) {
                                seqOtus[thisOtusSeqs[j]] = mergeOtu;
                            }
                        }
                    }
                }

                if (mergeOtu != "") {
                    list[mergeOtu] = mergedSeqs;
                }
            }

            otuTable->merge(ids);
            numOtus = otuTable->numOtus;
        }
    }
}
/******************************************************************************/
void Dataset::removeLineages(vector<string> contaminants, string trashTag) {

    if (!hasSequenceTaxonomy && !hasOtuTaxonomy ) { return; }

    Utils util;
    vector<vector<string> > conTax(contaminants.size());
    vector<vector<int> > conConfidenceThreshold(contaminants.size());
    vector<bool> conHasConfidences(contaminants.size(), false);

    // split contaminant taxonomies and confidences by level
    for (int i = 0; i < contaminants.size(); i++) {

        split(contaminants[i], ';', back_inserter(conTax[i]));
        conConfidenceThreshold[i] = util.removeConfidences(conTax[i]);

        if (sum(conConfidenceThreshold[i]) != 0) {
            conHasConfidences[i] = true;
        }
    }

    bool foundContaminants = false;
    if (hasSequenceTaxonomy) {
        // you have sequence taxonomies, remove contaminants
        // and reclassify otus (below)

        // remove seqs assigned to this taxonomy
        for (int i = 0; i < taxonomies.size(); i++) {

            if (tableSeqs[i]) {
                vector<string> userTax;
                split(taxonomies[i], ';', back_inserter(userTax));
                vector<int> userConfidences = util.removeConfidences(userTax);

                // if this seq is a contaminant, remove it
                if (util.searchTax(userTax, userConfidences,
                                   conHasConfidences,
                                   conTax, conConfidenceThreshold)) {
                    removeSequence(i, trashTag);
                    foundContaminants = true;
                }
            }
        }
    }else if (!hasSequenceTaxonomy && hasOtuTaxonomy) {
        // if you have only have otu classifications, remove otus

        vector<string> otuIds = otuTable->getOtuIds();
        vector<string> otuTaxes = otuTable->getTaxonomies();

        // remove otus assigned to this taxonomy
        for (int i = 0; i < otuTaxes.size(); i++) {
            vector<string> userTax;
            split(otuTaxes[i], ';', back_inserter(userTax));
            vector<int> userConfidences = util.removeConfidences(userTax);

            // if this seq is a contaminant, remove it
            if (util.searchTax(userTax, userConfidences,
                               conHasConfidences,
                               conTax, conConfidenceThreshold)) {
                otuTable->remove(otuIds[i], trashTag);
            }
        }
    }

    count->updateTotals();
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();

    // if you have otu classifications and sequence classifications
    // removing contaminant should force a reclassification of the otus
    // using the new data
    if (hasSequenceTaxonomy && hasOtuTaxonomy && foundContaminants) {
        classifyOtus();
    }

    if (hasOtuTaxonomy) {
        otuTable->updateTotals();
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
        numOtus = otuTable->numOtus;
    }
}
/******************************************************************************/
void Dataset::removeOtus(vector<string> namesToRemove,
                              vector<string> trashTags){

    if (!hasOtuData) { return; }

    if (namesToRemove.size() != trashTags.size()) {
        string message = "[ERROR]: Size mismatch. You must provide a trash";
        message += " code for each otu.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    for (int i = 0; i < namesToRemove.size(); i++) {
        otuTable->remove(namesToRemove[i], trashTags[i]);

        if (list.size() != 0) {
            auto it = list.find(namesToRemove[i]);

            if (it != list.end()) {
                // remove any sequences from removed otu
                for (int j = 0; j < it->second.size(); j++) {
                    removeSequence(it->second[j], trashTags[i], true);
                }
            }
        }
    }

    count->updateTotals();
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();

    otuTable->updateTotals();
    numSamples = otuTable->getNumSamples();
    numTreatments = otuTable->getNumTreatments();
    numOtus = otuTable->numOtus;
}
/******************************************************************************/
void Dataset::removeSequence(int index, string reasons, bool update) {
    // remove from tableSeqs and add trashCode
    tableSeqs[index] = false;
    trashCodes[index] += reasons + ",";
    numUnique--;
    runClassifyOtu = true;

    // remove from counts
    int abund = 1;
    if (update) {
        abund = count->remove(index);
        if (hasOtuData) {
            // seqIds assigned to otus
            if (list.size() != 0) {
                string otuName = seqOtus[index];

                if (otuTable->hasId(otuName)) {

                    // sequence abundance parsed by sample
                    vector<int> abunds = count->getAbundances(index);
                    // otuAbundances parsed by sample
                    vector<int> otuAbunds = otuTable->getAbundances(otuName);
                    // subtract this seqs abunds from otuAbunds
                    subtract(otuAbunds, abunds);
                    // wrap
                    vector<string> otuNames(1, otuName);
                    vector<vector<int>> otuAbundances(1);
                    otuAbundances[0] = otuAbunds;
                    // set otus new abundances
                    otuTable->setAbundances(otuNames, otuAbundances, reasons);
                }
            }
        }
    }else{
        abund = count->getAbundance(index);
    }

    // add to badAccnos
    vector<string> theseReasons;
    split(reasons, ',', back_inserter(theseReasons));

    for (int j = 0; j < theseReasons.size(); j++) {
        auto itBad = badAccnos.find(theseReasons[j]);

        if (itBad != badAccnos.end()) {
            // update counts of trashCode
            itBad->second[0]++;
            itBad->second[1] += abund;
        }else{
            // add new trashCode
            vector<int> badAbunds(2, 1);
            badAbunds[1] = abund;
            badAccnos[theseReasons[j]] = badAbunds;
        }
    }

    // update uniqueBad
    uniqueBad++;
}
/******************************************************************************/
void Dataset::removeSequences(vector<string> namesToRemove,
                         vector<string> trashTags){

    if (namesToRemove.size() != trashTags.size()) {
        string message = "[ERROR]: Size mismatch. You must provide a trash";
        message += " code for each sequence.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    for (int i = 0; i < namesToRemove.size(); i++) {
        auto it = seqIndex.find(namesToRemove[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            removeSequence(index, trashTags[i]);
        }else{
            string message = "[WARNING]: " + namesToRemove[i] + " is not in ";
            message += "your dataset, ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    count->updateTotals();
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();

    if (hasOtuData && (list.size() != 0)) {
        otuTable->updateTotals();
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
        numOtus = otuTable->numOtus;
    }
}
/******************************************************************************/
// for datasets without samples
void Dataset::setAbundance(vector<string> n, vector<int> abunds,
                            string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";

    for (int i = 0; i < n.size(); i++) {
        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            if (abunds[i] == 0) {
                removeSequence(index, reason);
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
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
}
/******************************************************************************/
// for datasets with samples
void Dataset::setAbundances(vector<string> n, vector<vector<int>> abunds,
                       string reason){

    if (n.size() != abunds.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    reason += ",";

    for (int i = 0; i < n.size(); i++) {
        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {
            int index = it->second;

            if (sum(abunds[i]) == 0) {
                removeSequence(index, reason);
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
    numSamples = count->getNumSamples();
    numTreatments = count->getNumTreatments();
}
/******************************************************************************/
// for datasets without samples
void Dataset::setOtuAbundance(vector<string> otuIDS, vector<int> abunds,
                     string reason) {
    if (hasOtuData) {
        vector<string> otusToRemove = otuTable->setAbundance(otuIDS, abunds,
                                                             reason);
        if (list.size() != 0) {
            for (int i = 0; i < otusToRemove.size(); i++) {

                auto it = list.find(otusToRemove[i]);

                if (it != list.end()) {
                    // remove any sequences from removed otu
                    for (int j = 0; j < it->second.size(); j++) {
                        removeSequence(it->second[j], reason, true);
                    }
                }
            }
        }

        otuTable->updateTotals();
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
        numOtus = otuTable->numOtus;
    }
}
/******************************************************************************/
// for datasets with samples
void Dataset::setOtuAbundances(vector<string> otuIDS,
                               vector<vector<int> > abunds,
                               string reason) {
    if (hasOtuData) {
        vector<string> otusToRemove = otuTable->setAbundances(otuIDS, abunds,
                                                              reason);

        if (list.size() != 0) {
            for (int i = 0; i < otusToRemove.size(); i++) {

                auto it = list.find(otusToRemove[i]);

                if (it != list.end()) {
                    // remove any sequences from removed otu
                    for (int j = 0; j < it->second.size(); j++) {
                        removeSequence(it->second[j], reason, true);
                    }
                }
            }
        }

        otuTable->updateTotals();
        numSamples = otuTable->getNumSamples();
        numTreatments = otuTable->getNumTreatments();
        numOtus = otuTable->numOtus;
    }
}
/******************************************************************************/
void Dataset::setSequences(vector<string> n, vector<string> s,
                      vector<string> c){
    if (n.size() != s.size()) {
        string message = "[ERROR]: Size mismatch. ids and sequences must be";
        message += " the same size.";
        RcppThread::Rcerr << endl << message << endl;
        throw Rcpp::exception(message.c_str());
    }

    bool hasComments = false;
    if (c.size() != 0) {
        if (c.size() != n.size()) {
            string message = "[ERROR]: Size mismatch. When providing comments,";
            message += " ids and comments must be the same size.";
            throw Rcpp::exception(message.c_str());
        }
        hasComments = true;
    }

    SeqReport report;

    for (int i = 0; i < n.size(); i++) {

        auto it = seqIndex.find(n[i]);

        if (it != seqIndex.end()) {

            int index = it->second;

            // update sequence
            seqs[index] = s[i];

            // update start, end, numbase, ambig, polymer, numn
            vector<int> reportResults = report.getReport(s[i]);

            starts[index] = reportResults[0];
            ends[index] = reportResults[1];
            lengths[index] = reportResults[2];
            ambigs[index] = reportResults[3];
            polymers[index] = reportResults[4];
            numns[index] = reportResults[5];

            // update comment
            if (hasComments) {
                comments[index] += " " + c[i];
            }

        }else{
            string message = "[WARNING]: " + n[i] + " is not in your dataset,";
            message += " ignoring.";
            RcppThread::Rcout << endl << message << endl;
        }
    }

    // set isAligned and aligned length
    getAlignedLength();
}
/******************************************************************************/


