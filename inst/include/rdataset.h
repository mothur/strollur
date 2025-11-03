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
const vector<float> nullFloatVector;  // used to pass blank floats
const vector< vector<float> > null2DFloatVector;  // used to pass blank vector
const vector< vector<int> > null2DIntVector;  // used to pass blank vector
const vector<bool> nullBoolVector;  // used to pass blank ints
const vector<double> nullDoubleVector;  // used to pass blank double
const vector< vector<double> > null2DDoubleVector;  // used to pass blank vector

/******************************************************************************/
struct Reference {
    string name, version, note, url, usage;

    Reference() : name(""), version(""), note(""), url(""), usage("") {}
    ~Reference() {}

    Reference(string n, string v = "", string u = "",
              string no = "", string ur = "") {
        name = n;
        version = v;
        note = no;
        url = ur;
        usage = u;
    }

    // This method lets Cereal know how to serialize.
    template<class Archive>
    void serialize(Archive & archive) {
        archive(name, version, usage, note, url);
    }
};

const Reference nullReference;
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
    vector<float> abunds;

    sampleAbunds() {}
    ~sampleAbunds() {}

    sampleAbunds(const sampleAbunds& sa) {
        sampleIndex = sa.sampleIndex;
        abunds = sa.abunds;
    }

    // no samples constructor
    sampleAbunds(const int i, const float a) {
        sampleIndex.push_back(i);
        abunds.push_back(a);
    }

    // sparse format constructor
    sampleAbunds(const vector<int> ind,
                 const vector<float> a,
                 const string mode)  {

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
    ~AbundTable();

    void clear();
    void clone(const AbundTable& abundTable);

    // 2, 3 or 4 columns: id, abundance, sample (optional -
    //                                   added when table includes sample data),
    //                                   treatment (optional -
    //                                   added when table includes treatment data),
    // used to export AbundTable
    const Rcpp::DataFrame getAbundanceTable(const vector<string>& outputNames,
                                      const vector<int>& names,
                                      const string tag = "sequence",
                                      const bool useNames = true);

    // names, sets abundance to 1
    double add(const vector<int>& names);

    double assignTreatments(const vector<string>& samples,
                          const vector<string>& treatments);

    // names, abundances, samples (optional), treatments (optional)
    double assignAbundance(vector<int>& names,
                         const vector<float>& abunds,
                         const vector<string> samples = nullVector,
                         const vector<string> treatments = nullVector);

    // set abundance parsed by sample
    void setAbundance(const int name, const vector<float>& abunds);
    // set abundance, assumes no samples if sample is blank
    void setAbundance(const int name, const float abund, const string sample = "");

    // removes id, returns abund. Be sure to run updateTotals after.
    // totals are not updated in function for time savings when removing multiple
    // ids. Only calc totals once rather than after each removal.
    int remove(const int name);
    void removeSamples(const vector<string>& samples);
    void updateTotals();
    // adds counts of idsToMerge[1-n] into idsToMerge[0]
    void merge(const vector<int>& idsToMerge);

    // vector containing total abundance for each id
    const vector<float> getTotalAbundances(const vector<int>& names);
    // total abundance for sequence, if sample is provided then abundance for
    // that sequence in that sample
    const float getAbundance(const int name, const string sample = "");
    // abundances by sample for id, (in the same order as the samples)
    const vector<float> getAbundances(const int id);
    // abundances by sample for ids
    const vector<vector<float>> getAbundances(const vector<int>& ids);
    // total number of sequences
    const double getTotal(const string sample = "");

    const int getNumSamples();
    const int getNumTreatments();
    // vector containing total abundance for each sample
    const vector<double> getSampleTotals();
    // vector containing total abundance for each treatment
    const vector<double> getTreatmentTotals();
    // vector containing names of samples
    const vector<string> getSamples();
    const vector<string> getTreatments();
    // maps sampleName to treatmentName
    const map<string, string> getSampleTreatmentAssignments();

    // does the table contain a sample
    // if name provided, does the sequence have this sample
    const bool hasSample(const string sample, const int name = -1);
    // does the table have sample information
    const bool hasSamples() { return hasSampleData; }

private:

    vector<sampleAbunds> counts;
    // numSamples is 1 for datasets without samples
    // numSamples equals the number of "good" samples in dataset
    int numSamples, numTreatments;
    double total;

    // sample name to index.
    map<string, int> sampleIndex;
    vector<string> sampleNames;
    // total abundance for each sample
    vector<double> sampleTotals;
    // are samples "present" in table
    vector<bool> tableSamples;

    // sample name to index.
    map<string, int> treatmentIndex;
    // total abundance for each sample
    vector<double> treatmentTotals;
    // are samples "present" in table
    vector<bool> tableTreatments;
    // sample index to treatment
    map<int, string> sampleTreatment;

    bool hasSampleData, hasTreatments;

    const int getSparseIndex(const int, const int);
    const vector<int> getSampleIndexes();
    void addSamples(vector<string> samples);
    void addTreatments(vector<string> treatments);
    void updateSampleTotals(const vector<float>& diffAbunds);

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
    BinTable(const string label);
    BinTable(const BinTable& binTable);
    ~BinTable();

    string label;
    bool hasListAssignments, hasBinTaxonomy, runClassify, hasBinReps;

    // ids, abundances, samples(optional)
    double assignAbundance(vector<string> ids,
                           const vector<float>& abundance,
                           const vector<string>& samples,
                           const vector<int>& seqIds,
                           AbundTable& count, bool update = false);
    double assignRepresentativeSequences(const vector<string>& binNames,
                                            const vector<int>& repNames);

    double assignTaxonomy(const vector<string>& ids,
                          const vector<string>& taxonomies);

    double assignTreatments(const vector<string>& samples,
                            const vector<string>& treatments);

    void clear(string tag = "");
    void clone(const BinTable& binTable);

    // id, bin_name, seq_id, abundance, sample, treatment, taxonomy, trashCode
    Rcpp::List exportBinTable();

    // string containing seqs in bin, comma separated
    const string get(const string binName, const vector<string>& seqNames);
    // total abundance for a given binId, optional sample
    const float getAbundance(const string binId, const string sample = "");
    // abundances for given binId broken down by sample
    const vector<float> getAbundances(const string binId);
    const vector<string> getListVector(const vector<string>& seqNames);
    // 2 column dataframe - bin_id, seq_id
    const Rcpp::DataFrame getList(const vector<string>& seqNames);
    // names of bins
    const vector<string> getIds();
    const int getNumBins();
    // 2 column dataframe - bin_id, abundance
    const Rcpp::DataFrame getRAbund();
    // vector of total abundances for each binId
    const vector<float> getRAbundVector();
    // 3 column dataframe - bin_id, abundance, sample
    const Rcpp::DataFrame getRepresentativeSequences(const vector<string>& seqNames,
                                                     const vector<string>& seqs);
    const Rcpp::DataFrame getShared();
    // abundances for each bin broken down by sample
    const vector<vector<float> > getSharedVector();
    // classifications of bins
    vector<string> getTaxonomies(const vector<string>& tax,
                                 AbundTable& count);

    // sample functions
    const int getNumSamples();
    const int getNumTreatments();
    const vector<string> getSamples();
    const vector<double> getSampleTotals();
    const map<string, string> getSampleTreatmentAssignments();
    const Rcpp::DataFrame getScrapReport();
    // trashCode, binCount, abundanceCount
    const Rcpp::DataFrame getScrapSummary();
    const vector<string> getTreatments();
    // vector containing total abundance for each treatment
    const vector<double> getTreatmentTotals();
    // total number of sequences
    const double getTotal(const string sample = "");
    const bool hasSample(const string sample);

    void merge(vector<string> binIds, string reason = "merged");
    // returns false if seqs are in different bins
    const bool okToMerge(const vector<int>& seqIds);
    void remove(const int seqId,
                AbundTable& count,
                const string reason, bool update = true);

    // remove given binId, return seqs removed
    vector<int> remove(const string binId, const string reason, bool update = true);
    void removeSamples(const vector<string>& samples);

    // for datasets without samples, return seqs removed by binRemoval
    vector<int> setAbundance(const vector<string>& binIds,
                             const vector<float>& abunds,
                             string reason = "merged");

    // for datasets with samples, return seqs removed by binRemoval
    vector<int> setAbundances(const vector<string>& binIds,
                              const vector<vector<float>>& abunds,
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
    vector<int> repSequences;

    // seqIndex -> binIndex
    // 2 -> 1
    // "seq2" -> "bin1"
    map<int, int> seqBins;

    // map reason for deletion to vector containing the number of bins removed
    // for that reason, and the total number of sequences removed by those bins
    // example: "cons_taxonomy" ->  c(10,  230) means cons_taxonomy merged
    // 10 bins that represented 230 sequences.
    map<string, vector<double> > badAccnos;
    double uniqueBad;

    // count table data
    // bin abundances, parsed by sample
    AbundTable binCount;

    const vector<int> getGoodIndexes();
    const vector<int> getIndexes(vector<string>&);
    void classify(const vector<string>& taxs,
                  AbundTable& count);
    string classifyBin(const int binIndex,
                       const vector<string>& tax,
                       AbundTable& count);
    void updateBinAbunds(map<int, map<string, float>>& binAbunds,
                         AbundTable& count,
                         const int bIndex, const int seqIndex,
                         const vector<string>& allSamples,
                         bool firstTimeSeq = true);
    double updateBins(map<int, map<string, float>>& binAbunds,
                      AbundTable& count,
                      bool hasSamples);

    friend class cereal::access; // Grants Cereal access

    template<class Archive>
    void serialize(Archive& ar) {
        ar(label, hasListAssignments, hasBinTaxonomy, repSequences, hasBinReps,
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
    Dataset(const string name, const int processors);
    Dataset(const Dataset& dataset);
    ~Dataset();

    string datasetName;
    int processors, alignmentLength; // -1 if unaligned
    bool isAligned, hasSequenceData, hasSequenceTaxonomy;
    double numUnique;

    void clear(const vector<string> tags = nullVector);
    Rcpp::List exportDataset(const vector<string> tags = nullVector);

    // add seqs
    double addSequences(const vector<string>& n,
                        vector<string> s = nullVector,
                        vector<string> c = nullVector,
                        Reference reference = nullReference);
    double addReferences(const vector<Reference>& refs);

    // names, abundances, samples(optional), treatments(optional)
    double assignSequenceAbundance(vector<string>& names,
                                   const vector<float>& abunds,
                                   const vector<string> samples = nullVector,
                                   const vector<string> treatments = nullVector);
    // binIds, abundances, samples(optional), seqIDs(optional), label (optional)
    double assignBins(const vector<string>& binIds,
                      vector<float> abunds = nullFloatVector,
                      vector<string> samples = nullVector,
                      vector<string> seqIDs = nullVector,
                      const string type = "otu");

    double assignSequenceTaxonomy(const vector<string>& names,
                                  const vector<string>& taxonomies);
    double assignBinRepresentativeSequences(const vector<string>& binNames,
                                  const vector<string>& repNames,
                                  const string type = "otu");
    double assignBinTaxonomy(const vector<string>& binIds,
                             const vector<string>& taxonomies,
                             const string type = "otu");

    double assignTreatments(const vector<string>& samples,
                            const vector<string>& treatments);

    const float getAbundance(const string name);
    // abundances for seq broken down by sample
    const vector<float> getAbundances(const string name);
    // total abundance for a given outID, optional sample
    const float getBinAbundance(const string binID, string type = "otu");
    // abundances for given binID broken down by sample
    const vector<float> getBinAbundances(const string binID, string type = "otu");
    // string containing sequence names for given binID
    const string getBin(const string binID, string type = "otu");
    // names of bins
    const vector<string> getBinIds(string type = "otu");
    // n columns: id, taxonomy split by level
    Rcpp::DataFrame getBinTaxonomyReport(string type = "otu");
    // 3 columns: bin_names, representative_names, representative_sequences
    const Rcpp::DataFrame getBinRepresentativeSequences(string type = "otu");
    // 2 column dataframe - bin_id, seq_id
    const vector<string> getBinTypes();
    const Rcpp::DataFrame getList(string type = "otu");
    const vector<string> getListVector(string type = "otu");
    const int getNumBins(string type = "otu");
    const int getNumSamples();
    const int getNumTreatments();
    // 2 column dataframe - bin_id, abundance
    const Rcpp::DataFrame getRAbund(string type = "otu");
    // vector of total abundances for each bin
    const vector<float> getRAbundVector(string type = "otu");
    const Rcpp::DataFrame getReferences();
    const vector<string> getSamples();
    const vector<double> getSampleTotals();
    const Rcpp::DataFrame getSampleTreatmentAssignments();
    const Rcpp::DataFrame getScrapReport(string mode = "sequence");
    // trashCode, uniqueCount, totalCount
    const Rcpp::List getScrapSummary();
    // vector[5][1] contains the abundance of seq5 in sample1
    const vector<vector<float>> getSeqsAbundsBySample();
    // total abundance for each sequence
    const vector<float> getSequenceAbundances();
    // 3 columns: id, sample, abundance
    const Rcpp::DataFrame getSequenceAbundanceTable();
    // sequence report: starts, ends, lengths, ambigs, polymers, numns
    const vector<string> getSequenceNames(string sample = "");
    const vector<vector<string> > getSequenceNamesBySample(vector<string> samples = nullVector);
    const Rcpp::DataFrame getSequenceReport();
    const vector<string> getSequences(string sample = "");
    const vector<vector<string> > getSequencesBySample(const vector<string> samples);
    // sequence summary summarizes sequence, and scrap reports
    const Rcpp::List getSequenceSummary();
    // n columns: id, taxonomy split by level
    Rcpp::DataFrame getSequenceTaxonomyReport();

    // abundances for each bin broken down by sample
    const vector<vector<float> > getSharedVector(string type = "otu");
    // 3 column dataframe - bin_id, abundance, sample
    const Rcpp::DataFrame getShared(string type = "otu");
    const double getTotal(string sample = "");
    const vector<string> getTreatments();
    const vector<double> getTreatmentTotals();
    const double getUniqueTotal(string sample = "");

    const bool hasSample(string sample);
    const bool hasListAssignments(string type = "otu");
    const bool hasSeqs();
    void mergeBins(const vector<string>& binIDS, string reason = "merged",
                   string type = "otu");
    void mergeSequences(const vector<string>&, string reason = "merged");

    void removeBins(const vector<string>& binIDs,
                    const vector<string>& trashTags,
                    string type = "otu");
    void removeLineages(const vector<string>& taxonomies,
                        string trashTag = "contaminant");
    void removeSamples(const vector<string>& samples);
    void removeSequences(const vector<string>& names,
                         const vector<string>& trashTags);

    // for datasets without samples
    void setAbundance(const vector<string>& names,
                      const vector<float>& abunds,
                      string reason = "update");
    // for datasets with samples
    void setAbundances(const vector<string>& names,
                       const vector<vector<float>>& abunds,
                       string reason = "update");
    // for datasets without samples
    void setBinAbundance(const vector<string>& binIDS,
                         const vector<float>& abunds,
                         string reason = "update", string type = "otu");
    // for datasets with samples
    void setBinAbundances(const vector<string>& binIDS,
                          const vector<vector<float>>& abunds,
                          string reason = "update", string type = "otu");

    // set sequence string and optionally comments
    void setSequences(const vector<string>& names,
                      const vector<string>& sequences,
                      vector<string> comments = nullVector);

    const Rcpp::RawVector serializeDataset();
    void loadFromSerialized(const Rcpp::RawVector);

private:
    // fasta data
    vector<string> names, seqs, comments, trashCodes;

    // fasta summary data
    vector<int> starts, ends, lengths, ambigs, polymers, numns;

    // sequence taxonomy assignments
    vector<string> taxonomies;

    vector<Reference> references;

    // boolean indicating if sequences is "good"
    vector<bool> tableSeqs;

    // maps sequence name to index in vectors above ^
    map<string, int> seqIndex;

    // map reason for deletion to vector containing unique and total counts
    // example: "pre_cluster" ->  c(10,  230) means precluster removed 10 unique
    // sequences that represented 230 total sequences.
    map<string, vector<double> > badAccnos;
    double uniqueBad;

    // count table data
    AbundTable count;

    // bin table, shared / rabund
    // tables label stores type
    // ie. otu, asv, phylotype
    vector<BinTable> binTables;

    // if unaligned, returns -1
    int getAlignedLength();
    const vector<int> getIncludedNamesIndexes();
    const vector<int> getIndexes(const vector<string>&);
    const bool hasBinTable(string type);
    void removeSequence(const int index,
                        const string reasons,
                        bool update = true,
                        bool removeFromBin = true);
    Rcpp::DataFrame fillTaxReport(string mode);
    const int getBinTableIndex(string type);

    friend class cereal::access; // Grants Cereal access

    template<class Archive>
    void serialize(Archive& ar) {
        ar(names, seqs, comments, trashCodes, references,
           starts, ends, lengths, ambigs, polymers, numns, taxonomies,
           tableSeqs, seqIndex, badAccnos, uniqueBad, count, binTables,
           datasetName, alignmentLength, isAligned,
           hasSequenceData, hasSequenceTaxonomy, numUnique, processors);
    }
};

/******************************************************************************/

#endif

