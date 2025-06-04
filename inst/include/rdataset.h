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
 * sampleIndex contains the sample indexes
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

    // no samples constructor
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

    // 2 or 3 columns: id, abundance, sample (optional -
    //                                   added when table includes sample data)
    // used to export SeqAbundTable
    Rcpp::DataFrame getSequenceAbundanceTable(vector<string> outputNames,
                                              vector<int> names);

    // names, sets abundance to 1
    void addSeqs(vector<int>& names);

    // names, abundances, samples (optional), treatments (optional)
    void assignSequenceAbundance(vector<int> names,
                               vector<int> abunds,
                               vector<string> samples = nullVector,
                               vector<string> treatments = nullVector);

    // set abundance parsed by sample - for datasets WITH samples
    void setAbundance(int name, vector<int> abunds);
    // set abundance - for datasets WITHOUT samples
    void setAbundance(int name, int abund);

    // removes sequence, returns abund. Be sure to run updateTotals after.
    // totals are not updated in function for time savings when removing multiple
    // sequences. Only calc totals once rather than after each removal.
    int removeSeq(int name);
    void updateTotals();
    // adds sequences counts of seqsToMerge[1-n] into seqsToMerge[0], optional
    // sample will only merge counts for that sample
    void mergeSeqs(vector<int> seqsToMerge, string sample = "");

    // vector containing total abundance for each sequence
    vector<int> getSeqsAbunds(vector<int> names);
    // total abundance for sequence, if sample is provided then abundance for
    // that sequence in that sample
    int getAbund(int name, string sample = "");
    // abundances by sample, in the same order as the samples
    vector<int> getAbunds(int name);
    // total number of sequences
    int getTotal(string sample = "");

    int getNumSamples();
    int getNumTreatments();
    // vector containing total abundance for each sample
    vector<int> getSampleTotals();
    // vector containing total abundance for each treatment
    vector<int> getTreatmentTotals();
    // vector containing names of samples
    vector<string> getSamples();
    vector<string> getTreatments();
    // does the table contain a sample
    // if name provided, does the sequence have this sample
    bool hasSample(string sample, int name = -1);
    // does the table have sample information
    bool hasSamples();

private:

    vector<seqCount> counts;
    // numSamples is 1 for datasets without samples
    // numSamples equals the number of "good" samples in dataset
    int total, numSamples, numTreatments;

    // sample name to index.
    map<string, int> sampleIndex;
    vector<string> sampleNames;
    // total abundance for each sample
    vector<int> sampleTotals;
    // are samples "present" in table
    vector<bool> tableSamples;

    // sample name to index.
    map<string, int> treatmentIndex;
    // total abundance for each sample
    vector<int> treatmentTotals;
    // are samples "present" in table
    vector<bool> tableTreatments;
    // sample index to treatment
    map<int, string> sampleTreatment;

    bool hasSampleData, hasTreatments;

    int getSparseIndex(int, int);
    void addSamples(vector<string> samples);
    void addTreatments(vector<string> treatments);
    void updateSampleTotals(vector<int> diffAbunds);

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
    bool hasContigsData, hasAlignData;

    int numSamples, numTreatments;
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
    void addAlignReport(vector<string>& n, vector<double>& ss,
                          vector<double>& sims, vector<int>& li);

    // make_contigs will create overlapLengths, overlapStarts, overlapEnds,
    // mismatches, and expectedErrors
    void addContigsReport(vector<string>& n, vector<int>& ol,
                          vector<int>& os, vector<int>& oe,
                          vector<int>& m, vector<double>& e);

    // names, abundances, samples(optional), treatments(optional)
    void assignSequenceAbundance(vector<string> names,
                               vector<int> abunds,
                               vector<string> samples = nullVector,
                               vector<string> treatments = nullVector);

    // **** functions for summarizing dataset **** //
    // sequence report: starts, ends, lengths, ambigs, polymers, numns
    Rcpp::DataFrame getSequenceReport();
    // sequence summary summarizes sequence, contigs and align reports
    Rcpp::List getSequenceSummary();
    // trashCode, uniqueCount, totalCount
    Rcpp::DataFrame getScrapSummary();
    Rcpp::DataFrame getScrapReport();
    // contigs report data: lengths, olengths, ostarts, oends, mismatches,
    //                      numns, ee
    Rcpp::DataFrame getContigsReport();
    // align summary data: search_score, sim_score, longest_insert
    Rcpp::DataFrame getAlignReport();
    // 3 columns: id, sample, abundance
    Rcpp::DataFrame getSequenceAbundanceTable();

    int getAbund(string name, string sample = "");
    // abundances for seq broken down by sample
    vector<int> getAbunds(string name);
    // total abundance for each sequence
    vector<int> getSeqsAbunds();
    // vector[5][1] contains the abundance of seq5 in sample1
    vector<vector<int>> getSeqsAbundsBySample();

    // sample functions
    vector<string> getSamples();
    vector<string> getTreatments();
    vector<int> getSampleTotals();
    vector<int> getTreatmentTotals();
    long long getTotal(string sample = "");
    long long getUniqueTotal(string sample = "");
    bool hasSample(string sample);

    // fasta sequence data
    vector<string> getNames(string sample = "");
    vector<vector<string> > getNamesBySample(vector<string> samples);
    vector<string> getSeqs(string sample = "");
    vector<vector<string> > getSeqsBySample(vector<string> samples);

    // modifiers
    void removeSeqs(vector<string> names, vector<string> trashTags);
    void mergeSeqs(vector<string>, string reason = "merged",
                   string sample = "");

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
    vector<int> olengths, ostarts, oends, mismatches;
    vector<double> ee;

    // fasta summary data
    vector<int> starts, ends, lengths, ambigs, polymers, numns;

    // alignment report
    vector<int> longestInsert;
    vector<float> searchScore, simScore;

    // sequence taxonomy assignments
    vector<string> taxonomies;

    // maps sequence name to index in summary vectors
    map<string, int> seqIndex;
    vector<bool> tableSeqs;

    // map reason for deletion to vector containing unique and total counts
    // example: "pre_cluster" ->  c(10,  230) means precluster removed 10 unique
    // sequences that represented 230 total sequences.
    map<string, vector<int> > badAccnos;
    int uniqueBad;
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

