#ifndef RCPP_rdataset_H_GEN_
#define RCPP_rdataset_H_GEN_

// io libraries
#include <string.h>
#include <iostream>

// containers
#include <vector>
#include <map>
//#include <set>
#include <string>

// Rcpp
#include <Rcpp.h>

using namespace std;

const vector<string> nullVector;  // used to pass blank vector
const vector< vector<string> > null2DVector;  // used to pass blank vector
const vector<int> nullIntVector;  // used to pass blank ints
const vector<bool> nullBoolVector;  // used to pass blank ints
const vector<char> nullCharVector;  // used to pass blank char
const vector<double> nullDoubleVector;  // used to pass blank double
const map<int, int> nullIntMap;
const pair<string, string> nullStringPair("", "");

//// [[Rcpp::plugins(cpp11)]]
/******************************************************************************/

class Dataset {

public:

    Dataset(string n);
    ~Dataset() = default;

    // public fields exposed through RCPP_MODULE
    string datasetName;
    bool isAligned;
    int numGroups;
    long long numUnique;

    // ********** public functions exposed through RCPP_MODULE ********** //

    SEXP getPointer();
    void clear();
    Rcpp::List exportDataset();
    string print();

    // add seqs
    void addSeqs(vector<string> n, vector<string> s,
                 vector<string> c = nullVector);

    // align_seqs will create searchScores, simScores and longestInserts
    void addAlignReport(vector<string> n, vector<double> ss,
                          vector<double> sims, vector<double> li);

    // make_contigs will create overlapLengths, overlapStarts, overlapEnds,
    // mismatches, and expectedErrors
    void addContigsReport(vector<string> n, vector<int> ol,
                          vector<int> os, vector<int> oe,
                          vector<int> m, vector<double> e);

    // names, groups, abundances
    void assignSampleAbundance(vector<string>, vector<string>, vector<int>);

    // **** functions for summarizing dataset **** //
    // fasta summary data: starts, ends, lengths, ambigs, polymers, numns
    vector<vector<int>> getFastaReport();
    // contigs sumary data: olengths, ostarts, oends, mismatches, ee
    vector<vector<double>> getContigsReport();
    // align summary data: search_score, sim_score, longest_insert
    vector<vector<double>> getAlignReport();
    // 3 columns: id, group, abundance
    Rcpp::DataFrame getSequenceAbundanceTable();

    int getAbund(string name, string group = "");
    // abundances for seq broken down by sample
    vector<int> getAbunds(string name);
    // total abundance for each sequence
    vector<int> getSeqsAbunds();
    // vector[5][1] contains the abundance of seq5 in sample1
    vector<vector<int>> getSeqsAbundsBySample();

    // group functions
    vector<string> getGroups(string name = "");
    vector<int> getGroupTotals(string name = "");
    long long getTotal(string group = "");
    bool hasGroup(string group);

    // fasta sequence data
    vector<string> getNames(string group = "");
    vector<vector<string> > getNamesBySample(vector<string> group);
    vector<string> getSeqs(string group = "");
    vector<vector<string> > getSeqsBySample(vector<string> group);

    // modifiers
    void reinstateSeqs(vector<string> trashTags);
    void removeSeqs(vector<string> names, vector<string> trashTags);
    void mergeSeqs(vector<string>, string reason = "merged", string group = "");

    // set sequence string and optionally comments
    void setSeqs(vector<string> names, vector<string> sequences,
                 vector<string> comments = nullVector);

    // set abundances
    void setAbundances(vector<string> names, vector<int> abunds,
                       string reason = "merged");

private:
    // fasta data
    vector<string> names, seqs, comments, trashCodes;

    // contigs report
    vector<double> olengths, ostarts, oends, mismatches, ee;

    // fasta summary data
    vector<int> starts, ends, lengths, ambigs, polymers, numns;

    // alignment report
    vector<double> search_score, sim_score, longest_insert;

    // sequence taxonomy assignments
    vector<string> taxonomies;

    // maps sequence name to index in vectors
    map<string, int> seqIndex;

    // map reason for deletion to vector containing unique and total counts
    // example: "pre_cluster" ->  c(10,  230) means precluster removed 10 unique
    // sequences that represented 230 total sequences.
    map<string, vector<int> > bad_accnos;
    long long total_bad, unique_bad, total;

    // count table data
    // SeqAbundTable* count;

    // otu table, asv / list / shared / constaxonomy
    // OTUTable* otuTable;

};

/******************************************************************************/

#endif

