#ifndef RCPP_rdataset_H_GEN_
#define RCPP_rdataset_H_GEN_

// io libraries
#include <string.h>
#include <iostream>

// containers
#include <vector>
#include <map>
#include <string>

// Rcpp
#include <Rcpp.h>
#include <RcppThread.h>
#include <cli/progress.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

using namespace std;

const vector<string> nullVector;  // used to pass blank vector
const vector< vector<string> > null2DVector;  // used to pass blank vector
const vector<int> nullIntVector;  // used to pass blank ints
const vector<bool> nullBoolVector;  // used to pass blank ints
const vector<double> nullDoubleVector;  // used to pass blank double

/******************************************************************************/
/*
 * The 'seqCount' struct will store abundance data for sequences in sparse form.
 *
 * sampleIndex contains the group indexes
 * abunds contains the abundances
 *
 * Let's assume 5 samples in the dataset, and the following are abundances by
 *  sample for a seq1
 *
 *  seq1 <- c(20, 0, 35, 0, 5) meaning seq1 has abundance of 20 in sample1
 *                                     seq1 has abundance of 0 in sample2
 *                                     seq1 has abundance of 35 in sample3
 *                                     seq1 has abundance of 0 in sample4
 *                                     seq1 has abundance of 5 in sample5
 *
 * would be stored in a sparse format
 *
 * sampleIndex <- c(0,2,4) - zero indexing of c++
 * abunds <- c(20, 35, 5)
 *
 */
struct seqCount {
    vector<int> sampleIndex;
    vector<int> abunds;

    seqCount() {}
    ~seqCount() {}

    // no groups constructor
    seqCount(int i, int a) {
        sampleIndex.push_back(i);
        abunds.push_back(a);
    }

    // sparse format constructor
    seqCount(vector<int> i, vector<int> a) : sampleIndex(i), abunds(a) {}

    // full format constructor
    seqCount(vector<int> fullAbunds) {
        for (int i = 0; i < fullAbunds.size(); i++){
            if (fullAbunds[i] != 0) {
                sampleIndex.push_back(i);
                abunds.push_back(fullAbunds[i]);
            }
        }
    }
};

/******************************************************************************/
/*
 * The 'SeqAbundTable' class will store abundance data for sequences.
 *
 * The "names" of the sequence are stored as indexes to save space.
 */

class SeqAbundTable {

public:

    SeqAbundTable();
    ~SeqAbundTable();

    void clear();
    string print();
    // 2 or 3 columns: id, abundance, group (optional -
    //                                   added when table includes group data)
    // used to export SeqAbundTable
    Rcpp::DataFrame getSequenceAbundanceTable(vector<string> outputNames,
                                              vector<int> names = nullIntVector);

    // names, sets abundance to 1
    void addSeqs(vector<int>& names);

    // names, abundances, groups (optional)
    void assignSampleAbundance(vector<int> names,
                               vector<int> abunds,
                               vector<string> groups = nullVector);

    // set abundance parsed by sample - for datasets WITH samples
    void setAbundance(int name, vector<int> abunds);
    // set abundance - for datasets WITHOUT samples
    void setAbundance(int name, int abund);

    // removes sequence from total and group totals
    void removeSeq(int name);
    // adds sequence count back into total and group totals
    void reinstateSeq(int name);
    // adds sequences counts of seqsToMerge[1-n] into seqsToMerge[0], optional
    // group will only merge counts for that sample
    void mergeSeqs(vector<int> seqsToMerge, string group = "");

    // vector containing total abundance for each sequence
    vector<int> getSeqsAbunds(vector<int> names);
    // total abundance for sequence, if group provided then abundance for that
    // sequence in that sample
    int getAbund(int name, string group = "");
    // abundances by sample, in the same order as the groups
    vector<int> getAbunds(int name);
    // total number of sequences
    int getTotal(string group = "");

    int getNumGroups();
    // vector containing total abundance for each sample
    vector<int> getGroupTotals();
    // vector containing names of samples
    // if name is provided then the names of samples where the seq is present
    vector<string> getGroups(int name = -1);
    // does the table contain a group
    // if name provided, does the sequence have this group
    bool hasGroup(string group, int name = -1);
    // does the table have group information
    bool hasGroups();

private:

    vector<seqCount> counts;
    // numGroups is 1 for datasets without groups
    // numGroups equals the number of "good" groups in dataset
    int total, numGroups;

    // sample name to index.
    map<string, int> groupIndex;
    // total abundance for each sample
    vector<int> groupTotals;
    // are samples "present" in table
    vector<bool> tableGroups;

    bool hasGroupData;

    int getSparseIndex(int, int);
    void addGroups(vector<string> groups);

};

/******************************************************************************/
/*
 * The 'Dataset' class will store data for DNA analysis.
 */

class Dataset {

public:

    Dataset(string n, int proc);
    ~Dataset();

    // public fields exposed through RCPP_MODULE
    string datasetName;
    bool isAligned;
    int numGroups;
    long long numUnique;
    int processors;

    // ********** public functions exposed through RCPP_MODULE ********** //

    SEXP getPointer();
    void clear();
    Rcpp::List exportDataset();
    string print();

    // add seqs
    void addSeqs(vector<string> n, vector<string> s,
                 vector<string> c);

    // align_seqs will create searchScores, simScores and longestInserts
    void addAlignReport(vector<string> n, vector<double> ss,
                          vector<double> sims, vector<double> li);

    // make_contigs will create overlapLengths, overlapStarts, overlapEnds,
    // mismatches, and expectedErrors
    void addContigsReport(vector<string> n, vector<int> ol,
                          vector<int> os, vector<int> oe,
                          vector<int> m, vector<double> e);

    // names, abundances, groups(optional)
    void assignSampleAbundance(vector<string> names,
                               vector<int> abunds,
                               vector<string> groups);

    // **** functions for summarizing dataset **** //
    // fasta summary data: starts, ends, lengths, ambigs, polymers, numns
    vector<vector<int>> getFastaReport();
    // contigs sumary data: olengths, ostarts, oends, mismatches, ee
    Rcpp::DataFrame getFastaSummary();
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
                 vector<string> comments);

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

    // maps sequence name to index in summary vectors
    map<string, int> seqIndex;
    vector<bool> tableSeqs;

    // map reason for deletion to vector containing unique and total counts
    // example: "pre_cluster" ->  c(10,  230) means precluster removed 10 unique
    // sequences that represented 230 total sequences.
    map<string, vector<int> > bad_accnos;
    int totalBad, uniqueBad;
    int alignmentLength;

    // count table data
    SeqAbundTable* count;

    // otu table, asv / list / shared / constaxonomy
    // OTUTable* otuTable;

    // if unaligned, returns -1
    int getAlignedLength(vector<string>);
    vector<int> getIncludedNamesIndexes();
    vector<int> getIndexes(vector<string>&);

};

/******************************************************************************/

#endif

