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
 * The 'sampleAbunds' struct will store abundance data for samples
 *    in sparse form.
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
 * This structure is also used to store OTU abundance data for samples
 *
 * otu1 <- c(10, 0, 0, 250, 1) meaning otu1 has abundance of 10 in sample1
 *                                     otu1 has abundance of 0 in sample2
 *                                     otu1 has abundance of 0 in sample3
 *                                     otu1 has abundance of 250 in sample4
 *                                     otu1 has abundance of 1 in sample5
 *
 */
struct sampleAbunds {
    vector<int> sampleIndex;
    vector<int> abunds;

    sampleAbunds() {}
    ~sampleAbunds() {}

    // no samples constructor
    sampleAbunds(int i, int a) {
        sampleIndex.push_back(i);
        abunds.push_back(a);
    }

    // sparse format constructor
    sampleAbunds(vector<int> i, vector<int> a) : sampleIndex(i), abunds(a) {}

    // full format constructor
    sampleAbunds(vector<int> fullAbunds) {
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
 * The 'AbundTable' class will store abundance data for samples.
 *
 * It is used by 'Dataset' to store sequence abundances by sample / treatment.
 * It is used by 'OTUTable' to store OTU abundances by sample.
 *
 * The "names" of the sequences or OTUs are stored as indexes to save space.
 */

class AbundTable {

public:

    AbundTable();
    ~AbundTable();

    void clear();

    // 2, 3 or 4 columns: id, abundance, sample (optional -
    //                                   added when table includes sample data),
    //                                   treatment (optional -
    //                                   added when table includes treatment data),
    // used to export AbundTable
    Rcpp::DataFrame getAbundanceTable(vector<string> outputNames,
                                              vector<int> names,
                                              bool includeTreatments = true);

    // names, sets abundance to 1
    void add(vector<int>& names);

    void assignTreatments(vector<string> samples,
                          vector<string> treatments);

    // names, abundances, samples (optional), treatments (optional)
    void assignAbundance(vector<int> names,
                         vector<int> abunds,
                         vector<string> samples = nullVector,
                         vector<string> treatments = nullVector);

    // set abundance parsed by sample - for datasets WITH samples
    void setAbundance(int name, vector<int> abunds);
    // set abundance - for datasets WITHOUT samples
    void setAbundance(int name, int abund);

    // removes id, returns abund. Be sure to run updateTotals after.
    // totals are not updated in function for time savings when removing multiple
    // ids. Only calc totals once rather than after each removal.
    int remove(int name);
    void updateTotals();
    // adds counts of idsToMerge[1-n] into idsToMerge[0]
    void merge(vector<int> idsToMerge);

    // vector containing total abundance for each id
    vector<int> getTotalAbundances(vector<int> names);
    // total abundance for sequence, if sample is provided then abundance for
    // that sequence in that sample
    int getAbundance(int name, string sample = "");
    // abundances by sample for id, (in the same order as the samples)
    vector<int> getAbundances(int id);
    // abundances by sample for ids
    vector<vector<int>> getAbundances(vector<int> ids);
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
    // maps sampleName to treatmentName
    map<string, string> getSampleTreatmentAssignments();

    // does the table contain a sample
    // if name provided, does the sequence have this sample
    bool hasSample(string sample, int name = -1);
    // does the table have sample information
    bool hasSamples();

private:

    vector<sampleAbunds> counts;
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
 * The 'OtuTable' class will store OTU data.
 */

class OtuTable {

public:

    OtuTable(string label);
    ~OtuTable();

    int numOtus;
    string label;

    // OTUIDs, abundances, samples(optional)
    void assignAbundance(vector<string> otuIDs, vector<int> abundance,
             vector<string> samples = nullVector);

    void assignTreatments(vector<string> samples,
                          vector<string> treatments);

    void clear();

    // names of OTUs
    vector<string> getOtuIds();
    // 2 column dataframe - otu_id, abundance
    Rcpp::DataFrame getRAbund();
    // vector of total abundances for each outID
    vector<int> getRAbundVector();
    // 3 column dataframe - otu_id, abundance, sample
    Rcpp::DataFrame getShared();
    // abundances for each OTU broken down by sample
    vector<vector<int> > getSharedVector();

    // total abundance for a given outID, optional sample
    int getAbundance(string otuID, string sample = "");
    // abundances for given otuID broken down by sample
    vector<int> getAbundances(string otuID);

    // sample functions
    int getNumSamples();
    int getNumTreatments();
    vector<string> getSamples();
    vector<int> getSampleTotals();
    Rcpp::DataFrame getScrapReport();
    // trashCode, otuCount, abundanceCount
    Rcpp::DataFrame getScrapSummary();
    vector<string> getTreatments();
    // vector containing total abundance for each treatment
    vector<int> getTreatmentTotals();
    // total number of sequences
    int getTotal(string sample = "");
    bool hasSample(string sample);
    bool hasId(string otuId);

    void merge(vector<string> otuIDS, string reason = "merged");
    // remove given outID
    void remove(string otuID, string reason, bool update = true);

    // for datasets without samples
    vector<string> setAbundance(vector<string> otuIDS, vector<int> abunds,
                      string reason = "merged");
    // for datasets with samples
    vector<string> setAbundances(vector<string> otuIDS,
                                 vector<vector<int>> abunds,
                       string reason = "merged");

    void setLabel(string label);

    void updateTotals();

private:

    // maps otu name to index
    map<string, int> otuIndex;
    // filter for "good" otus
    vector<bool> tableOtus;
    vector<string> otuNames, trashCodes;

    // map reason for deletion to vector containing the number of otus removed
    // for that reason, and the total number of sequences removed by those otus
    // example: "cons_taxonomy" ->  c(10,  230) means cons_taxonomy merged
    // 10 OTUs that represented 230 sequences.
    map<string, vector<int> > badAccnos;
    int uniqueBad;

    // count table data
    // otu abundances, parsed by sample
    AbundTable* otuCount;

    vector<int> getGoodIndexes();
    vector<int> getIndexes(vector<string>&);

};
/******************************************************************************/
/*
 * The 'Dataset' class will store sequence data for DNA analysis.
 */

class Dataset {

public:

    Dataset(string n, int proc);
    ~Dataset();

    // public fields exposed through RCPP_MODULE
    string datasetName, label;

    bool isAligned;
    // -1 if unaligned
    int alignmentLength;
    bool hasContigsData, hasAlignData, hasOtuData, hasSequenceData;

    int numSamples, numTreatments, numOtus;
    long long numUnique;
    int processors;

    // ********** public functions exposed through RCPP_MODULE ********** //

    SEXP getPointer();
    void clear();
    Rcpp::List exportDataset();
    string print();

    // add seqs
    void addSequences(vector<string> n, vector<string> s,
                 vector<string> c = nullVector);

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
    // otuIDS, abundances, samples(optional), seqIDs(optional), label (optional)
    void assignOtus(vector<string> otuIDS,
                                 vector<int> abunds = nullIntVector,
                                 vector<string> samples = nullVector,
                                 vector<string> seqIDs = nullVector);

    void assignTreatments(vector<string> samples,
                          vector<string> treatments);

    // **** functions for summarizing dataset **** //
    // sequence report: starts, ends, lengths, ambigs, polymers, numns
    Rcpp::DataFrame getSequenceReport();
    // sequence summary summarizes sequence, contigs and align reports
    Rcpp::List getSequenceSummary();
    // trashCode, uniqueCount, totalCount
    Rcpp::List getScrapSummary();
    Rcpp::DataFrame getScrapReport(string mode = "sequence");
    // contigs report data: lengths, olengths, ostarts, oends, mismatches,
    //                      numns, ee
    Rcpp::DataFrame getContigsReport();
    // align summary data: search_score, sim_score, longest_insert
    Rcpp::DataFrame getAlignReport();
    // 3 columns: id, sample, abundance
    Rcpp::DataFrame getSequenceAbundanceTable();

    int getAbundance(string name, string sample = "");
    // abundances for seq broken down by sample
    vector<int> getAbundances(string name);
    // total abundance for a given outID, optional sample
    int getOtuAbundance(string otuID, string sample = "");
    // abundances for given otuID broken down by sample
    vector<int> getOtuAbundances(string otuID);
    // string containing sequence names for given otuID
    string getOtu(string otuID);
    // total abundance for each sequence
    vector<int> getSequenceAbundances();
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
    vector<string> getSequences(string sample = "");
    vector<vector<string> > getSequencesBySample(vector<string> samples);

    // modifiers
    void removeOtus(vector<string> otuIDs, vector<string> trashTags);
    void removeSequences(vector<string> names, vector<string> trashTags);
    void mergeOtus(vector<string> otuIDS, string reason = "merged");
    void mergeSequences(vector<string>, string reason = "merged");

    // set sequence string and optionally comments
    void setSequences(vector<string> names, vector<string> sequences,
                 vector<string> comments = nullVector);

    // set abundances
    // for datasets with samples
    void setAbundances(vector<string> names, vector<vector<int>> abunds,
                       string reason = "merged");
    // for datasets without samples
    void setAbundance(vector<string> names, vector<int> abunds,
                       string reason = "merged");

    // for datasets without samples
    void setOtuAbundance(vector<string> otuIDS, vector<int> abunds,
                         string reason = "merged");
    // for datasets with samples
    void setOtuAbundances(vector<string> otuIDS, vector<vector<int>> abunds,
                          string reason = "merged");

    // OTU functions
    // vector string containing sequence names for each otu
    vector<string> getListVector();
    // 2 column dataframe - otu_id, seq_id
    Rcpp::DataFrame getList();
    // names of OTUs
    vector<string> getOtuIds();
    // 2 column dataframe - otu_id, abundance
    Rcpp::DataFrame getRAbund();
    // vector of total abundances for each OTU
    vector<int> getRAbundVector();
    // abundances for each OTU broken down by sample
    vector<vector<int> > getSharedVector();
    // 3 column dataframe - otu_id, abundance, sample
    Rcpp::DataFrame getShared();

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

    // count table data
    AbundTable* count;

    // otu table, shared / rabund
    OtuTable* otuTable;
    // otuName -> vector of sequence indexes
    // "otu1" -> 1,3,5
    // "otu1" -> "seq1,seq3,seq5"
    map<string, vector<int> > list;
    vector<string> seqOtus;

    // if unaligned, returns -1
    int getAlignedLength();
    vector<int> getIncludedNamesIndexes();
    vector<int> getIndexes(vector<string>&);
    // don't update totals when merging
    void removeSequence(int index, string reasons, bool update = true);
};

/******************************************************************************/

#endif

