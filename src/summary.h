//
//  summary.hpp
//  Mothur
//
//  Modified by Sarah Westcott on 5/14/25.
//  Copyright © 2025 Schloss Lab. All rights reserved.
//

#ifndef summary_hpp
#define summary_hpp

#include "utils.h"

class Summary {

public:

    Summary(int p);
    ~Summary() = default;

    // report[0] = lengths, report[1] = starts, report[2] = ends,
    // report[3] = ambigs, report[4] = homopolymers, report[5] = num_ns
    Rcpp::DataFrame summarizeFasta(vector<vector<int>>& report,
                                   vector<int> counts = nullIntVector);

    // report[0] = length, report[1] = overlap_length, report[2] = overlap_start,
    // report[3] = overlap_end, report[4] = mismatches, report[5] = num_ns
    Rcpp::DataFrame summarizeContigs(vector<vector<int>>& report,
                                     vector<int> counts = nullIntVector);

    // align[0] = search_scores, align[1] = sim_scores, longest_inserts, abunds
    Rcpp::DataFrame summarizeAlign(vector<vector<float>>& align,
                                   vector<int>& inserts,
                                   vector<int> counts = nullIntVector);

    // fasta - returns value at criteria level
    long long getStart(int value) {
        return (getValue(startPosition, value)); }
    long long getEnd(int value) {
        return (getValue(endPosition, value)); }
    long long getAmbig(int value) {
        return (getValue(ambigBases, value)); }
    long long getHomop(int value) {
        return (getValue(longHomoPolymer, value)); }
    long long getLength(int value) {
        return (getValue(seqLength, value)); }
    long long getNumNs(int value) {
        return (getValue(numNs, value)); }

    long long getOStart(int value) {
        return (getValue(ostartPosition, value)); }
    long long getOEnd(int value) {
        return (getValue(oendPosition, value)); }
    long long getMisMatches(int value) {
        return (getValue(misMatches, value)); }
    long long getOLength(int value) {
        return (getValue(oseqLength, value)); }

private:
    int processors;
    long long total, numUniques;
    bool hasCount;
    map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;
    map<int, long long> ostartPosition;
    map<int, long long> oendPosition;
    map<int, long long> oseqLength;
    map<int, long long> misMatches;
    map<int, long long> numNs;
    map<float, long long> sims;
    map<float, long long> scores;
    map<int, long long> inserts;

    void createThreadsFasta(vector<vector<int>>& fasta,
                                          vector<int>& counts);
    void createThreadsContigs(vector<vector<int>>& contigs,
                                       vector<int>& counts);
    void createThreadsAlign(vector<vector<float>>& align,
                            vector<int>& inserts,
                            vector<int>& counts);

    vector<int> getValues(map<int, long long>& positions);
    vector<float> getValues(map<float, long long>& positions);
    long long getValue(map<int, long long>& positions, double);
    vector<long long> getDefaults();

    //fasta
    //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    vector<int> getStart() { return (getValues(startPosition)); }
    vector<int> getEnd() { return (getValues(endPosition)); }
    vector<int> getAmbig() { return (getValues(ambigBases)); } //return
    vector<int> getLength() { return (getValues(seqLength)); } //returns
    vector<int> getHomop() { return (getValues(longHomoPolymer)); }

    // contigs
    // contigs overlap start - returns vector of 8 locations.
    // (min, 2.5, 25, 50, 75, 97.5, max, mean)
    vector<int> getOStart() { return (getValues(ostartPosition)); }
    vector<int> getOEnd() { return (getValues(oendPosition)); }
    vector<int> getOLength() { return (getValues(oseqLength)); }
    vector<int> getMisMatches() { return (getValues(misMatches)); }
    vector<int> getNumNs() { return (getValues(numNs)); }

    // align
    // align overlap start - returns vector of 8 locations.
    // (min, 2.5, 25, 50, 75, 97.5, max, mean)
    vector<float> getSims() { return (getValues(sims)); }
    vector<float> getScores() { return (getValues(scores)); }
    vector<int> getNumInserts() { return (getValues(inserts)); }
};
/******************************************************************************/

#endif /* summary_hpp */
