#ifndef RCPP_rdataset_H_GEN_
#define RCPP_rdataset_H_GEN_

// io libraries
#include <string.h>
#include <iostream>
#include <sstream> // For in-memory serialization

// containers
#include <vector>
#include <map>
#include <string>

// Rcpp
#include <Rcpp.h>
#include <RcppThread.h>
#include <cli/progress.h>

// serialization
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp> // Include for std::string
#include <cereal/types/vector.hpp> // Include for std::vector serialization
#include <cereal/types/set.hpp> // Include for std::set serialization
#include <cereal/types/map.hpp>   // Include for std::map serialization

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
// [[Rcpp::depends(Rcereal)]]

using namespace std;

const vector<string> nullVector;  // used to pass blank vector
const vector< vector<string> > null2DVector;  // used to pass blank vector
const vector<int> nullIntVector;  // used to pass blank ints
const vector< vector<int> > null2DIntVector;  // used to pass blank vector
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
 * This structure is also used to store bin abundance data for samples
 *
 * bin1 <- c(10, 0, 0, 250, 1) meaning bin1 has abundance of 10 in sample1
 *                                     bin1 has abundance of 0 in sample2
 *                                     bin1 has abundance of 0 in sample3
 *                                     bin1 has abundance of 250 in sample4
 *                                     bin1 has abundance of 1 in sample5
 *
 */
struct sampleAbunds {
    vector<int> sampleIndex;
    vector<int> abunds;

    sampleAbunds() {}
    ~sampleAbunds() {}

    sampleAbunds(const sampleAbunds& sa) {
        sampleIndex = sa.sampleIndex;
        abunds = sa.abunds;
    }

    // no samples constructor
    sampleAbunds(int i, int a) {
        sampleIndex.push_back(i);
        abunds.push_back(a);
    }

    // sparse format constructor
    sampleAbunds(vector<int> ind, vector<int> a, string mode)  {
        if (mode == "sparse") {
            sampleIndex = ind;
            abunds = a;
        }else{
            for (int i = 0; i < a.size(); i++){
                if (a[i] != 0) {
                    sampleIndex.push_back(ind[i]);
                    abunds.push_back(a[i]);
                }
            }
        }
    }

    // This method lets Cereal know how to serialize.
    template<class Archive>
    void serialize(Archive & archive) {
        archive(sampleIndex, abunds);
    }
};

/******************************************************************************/
/*
 * The 'AbundTable' class will store abundance data for samples.
 *
 * It is used by 'Dataset' to store sequence abundances by sample / treatment.
 * It is used by 'BinTable' to store bin abundances by sample.
 *
 * The "names" of the sequences or bins are stored as indexes to save space.
 */

class AbundTable {

public:

    AbundTable();
    AbundTable(const AbundTable& abundTable);
    ~AbundTable();

    void clear();
    void clone(const AbundTable& abundTable);

    // 2, 3 or 4 columns: id, abundance, sample (optional -
    //                                   added when table includes sample data),
    //                                   treatment (optional -
    //                                   added when table includes treatment data),
    // used to export AbundTable
    Rcpp::DataFrame getAbundanceTable(vector<string> outputNames,
                                              vector<int> names,
                                              string tag = "sequence");

    // names, sets abundance to 1
    double add(vector<int>& names);

    double assignTreatments(vector<string> samples,
                          vector<string> treatments);

    // names, abundances, samples (optional), treatments (optional)
    double assignAbundance(vector<int> names,
                         vector<int> abunds,
                         vector<string> samples = nullVector,
                         vector<string> treatments = nullVector);

    // set abundance parsed by sample
    void setAbundance(int name, vector<int> abunds);
    // set abundance, assumes no samples if sample is blank
    void setAbundance(int name, int abund, string sample = "");

    // removes id, returns abund. Be sure to run updateTotals after.
    // totals are not updated in function for time savings when removing multiple
    // ids. Only calc totals once rather than after each removal.
    int remove(int name);
    void removeSamples(vector<string> samples);
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
    bool hasSamples() { return hasSampleData; }

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
    vector<int> getSampleIndexes();
    void addSamples(vector<string> samples);
    void addTreatments(vector<string> treatments);
    void updateSampleTotals(vector<int> diffAbunds);

    friend class cereal::access; // Grants Cereal access

    template<class Archive>
    void serialize(Archive& ar) {
        ar(counts, total, numSamples, numTreatments,
           sampleIndex, sampleNames, sampleTotals, tableSamples,
           treatmentIndex, treatmentTotals, tableTreatments,
           sampleTreatment, hasSampleData, hasTreatments);
    }
};

/******************************************************************************/
/*
 * The 'BinTable' class will store asv / otu / phylotype data.
 */

class BinTable {

public:

    BinTable();
    BinTable(string label);
    BinTable(const BinTable& binTable);
    ~BinTable();

    int numBins;
    string label;
    bool hasListAssignments, hasBinTaxonomy, runClassify;

    // ids, abundances, samples(optional)
    double assignAbundance(vector<string> ids, vector<int> abundance,
             vector<string> samples, vector<int> seqIds, AbundTable& count,
             bool update = false);

    double assignTaxonomy(vector<string> ids, vector<string> taxonomies);

    double assignTreatments(vector<string> samples,
                          vector<string> treatments);

    void clear(string tag = "");
    void clone(const BinTable& binTable);

    // string containing seqs in bin, comma separated
    string get(string binName, vector<string>& seqNames);
    // total abundance for a given binId, optional sample
    int getAbundance(string binId, string sample = "");
    // abundances for given binId broken down by sample
    vector<int> getAbundances(string binId);
    vector<string> getListVector(vector<string>& seqNames);
    // 2 column dataframe - bin_id, seq_id
    Rcpp::DataFrame getList(vector<string>& seqNames);
    // names of bins
    vector<string> getIds();
    // 2 column dataframe - bin_id, abundance
    Rcpp::DataFrame getRAbund();
    // vector of total abundances for each binId
    vector<int> getRAbundVector();
    // 3 column dataframe - bin_id, abundance, sample
    Rcpp::DataFrame getShared();
    // abundances for each bin broken down by sample
    vector<vector<int> > getSharedVector();
    // classifications of bins
    vector<string> getTaxonomies(vector<string>& tax, AbundTable& count);

    // sample functions
    int getNumSamples();
    int getNumTreatments();
    vector<string> getSamples();
    vector<int> getSampleTotals();
    Rcpp::DataFrame getScrapReport();
    // trashCode, binCount, abundanceCount
    Rcpp::DataFrame getScrapSummary();
    vector<string> getTreatments();
    // vector containing total abundance for each treatment
    vector<int> getTreatmentTotals();
    // total number of sequences
    int getTotal(string sample = "");
    bool hasSample(string sample);

    void merge(vector<string> binIds, string reason = "merged");
    // returns false if seqs are in different bins
    bool okToMerge(vector<int> seqIds);
    void remove(int seqId, AbundTable& count, string reason, bool update = true);
    // remove given binId, return seqs removed
    vector<int> remove(string binId, string reason, bool update = true);
    void removeSamples(vector<string> samples);

    // for datasets without samples, return seqs removed by binRemoval
    vector<int> setAbundance(vector<string> binIds, vector<int> abunds,
                      string reason = "merged");
    // for datasets with samples, return seqs removed by binRemoval
    vector<int> setAbundances(vector<string> binIds,
                                 vector<vector<int>> abunds,
                       string reason = "merged");
    void updateTotals();

private:

    // maps bin name to index
    map<string, int> binIndex;
    // filter for "good" bins
    vector<bool> tableBins;
    vector<string> binNames, trashCodes, taxonomies;

    // binList[binIndex] -> vector of sequence indexes
    // binList[1] (aka "bin1") -> 1,3,5
    // "bin1" -> "seq1,seq3,seq5"
    vector<set<int> > binList;
    // seqIndex -> binIndex
    // 2 -> 1
    // "seq2" -> "bin1"
    map<int, int> seqBins;

    // map reason for deletion to vector containing the number of bins removed
    // for that reason, and the total number of sequences removed by those bins
    // example: "cons_taxonomy" ->  c(10,  230) means cons_taxonomy merged
    // 10 bins that represented 230 sequences.
    map<string, vector<int> > badAccnos;
    int uniqueBad;

    // count table data
    // bin abundances, parsed by sample
    AbundTable binCount;

    vector<int> getGoodIndexes();
    vector<int> getIndexes(vector<string>&);
    void classify(vector<string>& taxs, AbundTable& count);
    string classifyBin(int binIndex, vector<string>& tax, AbundTable& count);
    void updateBinAbunds(map<int, map<string, int>>& binAbunds,
                         AbundTable& count, int bIndex, int seqIndex,
                         vector<string>& allSamples, bool firstTimeSeq = true);
    double updateBins(map<int, map<string, int>>& binAbunds, AbundTable& count,
                    bool hasSamples);

    friend class cereal::access; // Grants Cereal access

    template<class Archive>
    void serialize(Archive& ar) {
        ar(numBins, label, hasListAssignments, hasBinTaxonomy,
           binIndex, tableBins, binNames, trashCodes, taxonomies,
           runClassify, binList, seqBins, badAccnos, uniqueBad, binCount);
    }

};
/******************************************************************************/
/*
 * The 'Dataset' class will store sequence data for DNA analysis.
 */


class Dataset {

public:

    Dataset();
    Dataset(string name, int processors);
    Dataset(const Dataset& dataset);
    ~Dataset();

    // public fields exposed through RCPP_MODULE
    string datasetName;
    int alignmentLength; // -1 if unaligned
    bool isAligned, hasSequenceData, hasSequenceTaxonomy;
    long long numUnique;
    int processors;

    // ********** public functions exposed through RCPP_MODULE ********** //
    void clear(vector<string> tags = nullVector);
    Rcpp::List exportDataset();

    // add seqs
    double addSequences(vector<string> n, vector<string> s = nullVector,
                 vector<string> c = nullVector);

    // names, abundances, samples(optional), treatments(optional)
    double assignSequenceAbundance(vector<string> names,
                               vector<int> abunds,
                               vector<string> samples = nullVector,
                               vector<string> treatments = nullVector);
    // binIds, abundances, samples(optional), seqIDs(optional), label (optional)
    double assignBins(vector<string> binIds,
                                 vector<int> abunds = nullIntVector,
                                 vector<string> samples = nullVector,
                                 vector<string> seqIDs = nullVector,
                                 string type = "otu");

    double assignSequenceTaxonomy(vector<string> names, vector<string> taxonomies);
    double assignBinTaxonomy(vector<string> binIds, vector<string> taxonomies,
                           string type = "otu");

    double assignTreatments(vector<string> samples,
                          vector<string> treatments);

    int getAbundance(string name);
    // abundances for seq broken down by sample
    vector<int> getAbundances(string name);
    // total abundance for a given outID, optional sample
    int getBinAbundance(string binID, string type = "otu");
    // abundances for given binID broken down by sample
    vector<int> getBinAbundances(string binID, string type = "otu");
    // string containing sequence names for given binID
    string getBin(string binID, string type = "otu");
    // names of bins
    vector<string> getBinIds(string type = "otu");
    // n columns: id, taxonomy split by level
    Rcpp::DataFrame getBinTaxonomyReport(string type = "otu");
    // 2 column dataframe - bin_id, seq_id
    Rcpp::DataFrame getList(string type = "otu");
    vector<string> getListVector(string type = "otu");
    int getNumBins(string type = "otu");
    int getNumSamples();
    int getNumTreatments();
    // 2 column dataframe - bin_id, abundance
    Rcpp::DataFrame getRAbund(string type = "otu");
    // vector of total abundances for each bin
    vector<int> getRAbundVector(string type = "otu");
    vector<string> getSamples();
    vector<int> getSampleTotals();
    Rcpp::DataFrame getScrapReport(string mode = "sequence");
    // trashCode, uniqueCount, totalCount
    Rcpp::List getScrapSummary();
    // vector[5][1] contains the abundance of seq5 in sample1
    vector<vector<int>> getSeqsAbundsBySample();
    // total abundance for each sequence
    vector<int> getSequenceAbundances();
    // 3 columns: id, sample, abundance
    Rcpp::DataFrame getSequenceAbundanceTable();
    // sequence report: starts, ends, lengths, ambigs, polymers, numns
    vector<string> getSequenceNames(string sample = "");
    vector<vector<string> > getSequenceNamesBySample(vector<string> samples = nullVector);
    Rcpp::DataFrame getSequenceReport();
    vector<string> getSequences(string sample = "");
    vector<vector<string> > getSequencesBySample(vector<string> samples);
    // sequence summary summarizes sequence, and scrap reports
    Rcpp::List getSequenceSummary();
    // n columns: id, taxonomy split by level
    Rcpp::DataFrame getSequenceTaxonomyReport();

    // abundances for each bin broken down by sample
    vector<vector<int> > getSharedVector(string type = "otu");
    // 3 column dataframe - bin_id, abundance, sample
    Rcpp::DataFrame getShared(string type = "otu");
    long long getTotal(string sample = "");
    vector<string> getTreatments();
    vector<int> getTreatmentTotals();
    long long getUniqueTotal(string sample = "");

    bool hasSample(string sample);
    bool hasListAssignments(string type = "otu");
    bool hasSeqs();
    void mergeBins(vector<string> binIDS, string reason = "merged",
                   string type = "otu");
    void mergeSequences(vector<string>, string reason = "merged");

    void removeBins(vector<string> binIDs, vector<string> trashTags,
                    string type = "otu");
    void removeLineages(vector<string> taxonomies,
                        string trashTag = "contaminant");
    void removeSamples(vector<string> samples);
    void removeSequences(vector<string> names, vector<string> trashTags);

    // for datasets without samples
    void setAbundance(vector<string> names, vector<int> abunds,
                      string reason = "update");
    // for datasets with samples
    void setAbundances(vector<string> names, vector<vector<int>> abunds,
                       string reason = "update");
    // for datasets without samples
    void setBinAbundance(vector<string> binIDS, vector<int> abunds,
                         string reason = "update", string type = "otu");
    // for datasets with samples
    void setBinAbundances(vector<string> binIDS, vector<vector<int>> abunds,
                          string reason = "update", string type = "otu");

    // set sequence string and optionally comments
    void setSequences(vector<string> names, vector<string> sequences,
                 vector<string> comments = nullVector);

    Rcpp::RawVector serializeDataset();
    void loadFromSerialized(Rcpp::RawVector);

private:
    // fasta data
    vector<string> names, seqs, comments, trashCodes;

    // fasta summary data
    vector<int> starts, ends, lengths, ambigs, polymers, numns;

    // sequence taxonomy assignments
    vector<string> taxonomies;

    // boolean indicating if sequences is "good"
    vector<bool> tableSeqs;

    // maps sequence name to index in vectors above ^
    map<string, int> seqIndex;

    // map reason for deletion to vector containing unique and total counts
    // example: "pre_cluster" ->  c(10,  230) means precluster removed 10 unique
    // sequences that represented 230 total sequences.
    map<string, vector<int> > badAccnos;
    int uniqueBad;

    // count table data
    AbundTable count;

    // bin table, shared / rabund
    // tables label stores type
    // ie. otu, asv, phylotype
    vector<BinTable> binTables;

    // if unaligned, returns -1
    int getAlignedLength();
    vector<int> getIncludedNamesIndexes();
    vector<int> getIndexes(vector<string>&);
    bool hasBinTable(string type);
    void removeSequence(int index, string reasons, bool update = true,
                        bool removeFromBin = true);
    Rcpp::DataFrame fillTaxReport(string mode);
    int getBinTableIndex(string type);

    friend class cereal::access; // Grants Cereal access

    template<class Archive>
    void serialize(Archive& ar) {
        ar(names, seqs, comments, trashCodes,
           starts, ends, lengths, ambigs, polymers, numns, taxonomies,
           tableSeqs, seqIndex, badAccnos, uniqueBad, count, binTables,
           datasetName, alignmentLength, isAligned,
           hasSequenceData, hasSequenceTaxonomy, numUnique, processors);
    }
};

/******************************************************************************/

#endif

