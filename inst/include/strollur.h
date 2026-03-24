#ifndef RCPP_strollur_H_GEN_
#define RCPP_strollur_H_GEN_

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
#include <cereal/types/list.hpp>

//[[Rcpp::plugins(cpp11)]]
//[[Rcpp::depends(RcppThread)]]
//[[Rcpp::depends(Rcereal)]]

using namespace std;

const set<string> nullSet;  // used to pass blank set
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
            for (size_t i = 0; i < a.size(); i++){
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
    Rcpp::DataFrame getAbundanceTable(const vector<string>& outputNames,
                                      const vector<int>& names,
                                      const string& tag = "sequence",
                                      const bool useNames = true) const;

    // names, sets abundance to 1
    double add(const vector<int>& names);

    double assignTreatments(const vector<string>& samples,
                          const vector<string>& treatments);

    // names, abundances, samples (optional), treatments (optional)
    double assignAbundance(vector<int>& names,
                         const vector<float>& abunds,
                         const vector<string>& samples = nullVector,
                         const vector<string>& treatments = nullVector);

    // set abundance parsed by sample
    void setAbundance(const int name, const vector<float>& abunds);
    // set abundance, assumes no samples if sample is blank
    void setAbundance(const int name, const float abund, const string& sample = "");
    // index -> (sampleName -> abundance)
    void setAbundances(const map<int, map<string, float>>& binAbunds);

    // removes id, returns abund. Be sure to run updateTotals after.
    // totals are not updated in function for time savings when removing multiple
    // ids. Only calc totals once rather than after each removal.
    int remove(const int name);
    void removeSamples(const vector<string>& samples);
    void updateTotals();
    // adds counts of idsToMerge[1-n] into idsToMerge[0]
    void merge(const vector<int>& idsToMerge);

    // vector containing total abundance for each id
    vector<float> getTotalAbundances(const vector<int>& names) const;
    // total abundance for sequence, if sample is provided then abundance for
    // that sequence in that sample
   float getAbundance(const int name, const vector<string>& samples = nullVector) const;
    // abundances by sample for id, (in the same order as the samples)
    vector<float> getAbundances(const int id) const;
    // abundances by sample for ids
    vector<vector<float>> getAbundances(const vector<int>& ids) const;
    // results[0][1:numSeqsInSample0] -> sample 0's abundances for each sequence
    vector<vector<float>> getAbundanceBySample(const vector<int>& ids,
        vector<string> samplesToSelect = nullVector) const;
    // total number of sequences
    double getTotal(const string& sample = "") const;

    int getNumSamples(const int name = -1) const;
    int getNumTreatments() const;
    // vector containing total abundance for each sample
    vector<double> getSampleTotals() const;
    // vector containing total abundance for each treatment
    vector<double> getTreatmentTotals() const;
    // vector containing names of samples
    vector<string> getSamples(const int name = -1) const;
    vector<string> getTreatments() const;
    // maps sampleName to treatmentName
    map<string, string> getSampleTreatmentAssignments() const;

    // does the table contain a sample
    // if name provided, does the sequence have this sample
    bool hasSample(const string& sample, const int name = -1) const;
    bool hasSamples(const vector<string>& samples, const int name = -1) const;
    // does the table have sample information
    bool hasSamplesData() const { return hasSampleData; }

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

    int getSparseIndex(const int, const int) const;
    vector<int> getSampleIndexes() const;
    void addSamples(vector<string> samples);
    void addTreatments(vector<string> treatments);
    void updateSampleTotals(const vector<float>& diffAbunds);

    friend class cereal::access; // Grants Cereal access

    // 14 items
    template<class Archive>
    void serialize(Archive& ar) {
        ar(counts, numSamples, numTreatments, total,
           sampleIndex, tableSamples, sampleNames, sampleTotals,
           sampleTreatment, treatmentIndex, treatmentTotals, tableTreatments,
           hasSampleData, hasTreatments);
    }
};
/******************************************************************************/
/*
 * The 'Report' class store a data.frame as vectors to allow for serialization
 */

class Report {

public:

    Report();
    ~Report() {}

    void addReport(Rcpp::DataFrame& report);
    void clear();
    Rcpp::DataFrame getReport(set<string> datasetNames);

    Rcpp::DataFrame summarizeReport(set<string> datasetNames,
                                    int proc, vector<float> counts);

    bool hasReport;

private:

    map<int, string> columnNames;
    string sequence_name;
    int sequence_name_col, numRows;
    bool hasStr, hasInt, hasNum, hasLog, hasColumnNames;

    // index in df -> df[index] values
    map<int, vector<string>> strColumns;
    map<int, vector<int>> intColumns;
    map<int, vector<double>> numColumns;
    map<int, vector<bool>> logColumns;

    void pruneReport(set<string> names);

    friend class cereal::access; // Grants Cereal access

    // 14
    template<class Archive>
    void serialize(Archive& ar) {
        ar(hasReport, columnNames, sequence_name, sequence_name_col, numRows,
           hasStr, hasInt, hasNum, hasLog, hasColumnNames,
           strColumns, intColumns, numColumns, logColumns);
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

    string label;
    bool hasBinTaxonomy, hasBinReps, runClassify;

    // ids, seqids
    double assignBins(AbundTable& count, vector<string> ids,
                           const vector<int>& seqIds);
    double assignRepresentativeSequences(const vector<string>& binNamesVector,
                                            const vector<int>& repNames);

    double assignTaxonomy(const vector<string>& ids,
                          const vector<string>& taxonomies);

    void clear();
    void clone(const BinTable& binTable);

    // id, bin_name, seq_id, abundance, sample, treatment, taxonomy, trashCode
    Rcpp::List exportBinTable(AbundTable& count);

    // names of bins
    vector<string> getIds(const AbundTable& count,
                                const vector<string>& sample = nullVector,
                                bool distinct = false) const;
    const vector<string> getListVector(const vector<string>& seqNames);
    // 2 column dataframe - bin_id, seq_id
    const Rcpp::DataFrame getList(const vector<string>& seqNames);

    int getNumBins(const AbundTable& count,
                         const vector<string>& samples = nullVector,
                         bool distinct = false) const;
    // 2 column dataframe - bin_id, abundance
    Rcpp::DataFrame getRAbund(const AbundTable& count) const;
    // vector of total abundances for each binId
    vector<float> getRAbundVector(const AbundTable& count) const;
    // 3 column dataframe - bin_id, abundance, sample
    Rcpp::DataFrame getRepresentativeSequences(const vector<string>& seqNames,
                                                const vector<string>& seqs) const;
    Rcpp::DataFrame getScrapReport();
    // type, trashCode, binCount, abundanceCount
    Rcpp::DataFrame getScrapSummary() const;

    Rcpp::DataFrame getShared(const AbundTable& count) const;
    // abundances for each bin broken down by sample
    vector<vector<float>> getSharedVector(const AbundTable& count) const;
    // classifications of bins
    vector<string> getTaxonomies(const vector<string>& tax,
                                 const AbundTable& count);

    void merge(vector<string> binIds, const string& reason = "merged");
    // returns false if seqs are in different bins
    bool okToMerge(const vector<int>& seqIds) const;

    // remove given binId, return seqs removed
    vector<int> remove(const string& binId, const string& reason);
    void removeSeq(const int seqId, const string& reason);

private:

    // maps bin name to index
    map<string, int> binIndex;
    // filter for "good" bins
    vector<bool> tableBins;
    // used for when a bin is removed to indicate how many sequences were in
    // original otu
    vector<float> originalBinAbunds;
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

    void classify(const vector<string>& taxs,
                  const AbundTable& count);
    string classifyBin(const int binIndex,
                       const vector<string>& tax,
                       const AbundTable& count) const;

    // string containing seqs in bin, comma separated
    const string get(const string binName, const vector<string>& seqNames);
    // total abundance for a given bin
    float getAbundance(const AbundTable& count,
                             const string& binName,
                             const vector<string>& samples = nullVector) const;
    // total abundance for a given bin
    float getAbundance(const AbundTable& count,
                             const int binId,
                             const vector<string>& samples = nullVector) const;
    // abundances for given binId broken down by sample
    vector<float> getAbundances(const AbundTable& count, const string& binName) const;
    // abundances for given binId broken down by sample
    vector<float> getAbundances(const AbundTable& count, const int binId) const;
    // total abundance for each bin
    vector<float> getAbundances(const AbundTable& count, bool onlyGood = true) const;
    // 2, 3 or 4 columns: id, abundance, sample (optional -
    //                                   added when table includes sample data),
    //                                   treatment (optional -
    //                                   added when table includes treatment data),
    // used to export BinTable and Shared functions
    Rcpp::DataFrame getAbundanceTable(const AbundTable& count,
                                            const bool useNames = true) const;

    vector<int> getGoodIndexes(const AbundTable& count) const;
    vector<int> getIndexes(vector<string>&) const;
    vector<int> remove(const int binId, const string& reason);

    friend class cereal::access; // Grants Cereal access

    // 15
    template<class Archive>
    void serialize(Archive& ar) {
        ar(label, hasBinTaxonomy, hasBinReps, runClassify,
           binIndex, tableBins, originalBinAbunds, binList,
           binNames, trashCodes, taxonomies, repSequences,
           seqBins, badAccnos, uniqueBad);
    }
};
/******************************************************************************/
/*
 * The 'Dataset' class will store sequence data for DNA analysis.
 */


class Dataset {

public:

    Dataset();
    Dataset(string name, const int processors);
    Dataset(const Dataset& dataset);
    ~Dataset();

    string datasetName;
    int processors, alignmentLength; // -1 if unaligned
    bool isAligned, hasSequenceData, hasSequenceTaxonomy;
    double numUnique;

    void clear();
    Rcpp::List exportDataset(const vector<string> tags = nullVector);

    // add seqs
    double addSequences(const vector<string>& n,
                        vector<string> s = nullVector,
                        vector<string> c = nullVector,
                        Reference reference = nullReference);
    double addReferences(const vector<Reference>& refs);
    void addReport(Rcpp::DataFrame& report, string type);
    void addMetadata(Rcpp::DataFrame& metadata);

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


    const vector<vector<float> > getSequenceAbundanceBySample(vector<string> samples = nullVector);

    // names of bins
    const vector<string> getBinIds(string type = "otu",
                                   vector<string> samples = nullVector,
                                   bool distinct = false);
    // n columns: id, taxonomy split by level
    Rcpp::DataFrame getBinTaxonomyReport(string type = "otu");
    // 3 columns: bin_names, representative_names, representative_sequences
    const Rcpp::DataFrame getBinRepresentativeSequences(string type = "otu");
    // 2 column dataframe - bin_id, seq_id
    const vector<string> getBinTypes();
    // fasta data.frame 2 or 3 columns, sequence_names, sequences, comments
    const Rcpp::DataFrame getFastaReport();

    const Rcpp::DataFrame getList(string type = "otu");
    const vector<string> getListVector(string type = "otu");
    const int getNumBins(string type = "otu",
                         vector<string> samples = nullVector,
                         bool distinct = false);
    const int getNumSamples();
    const int getNumTreatments();
    const int getNumResourceReferences();

    const Rcpp::DataFrame getMetadata();
    const Rcpp::DataFrame getReferences();
    const Rcpp::DataFrame getReports(string type);
    const vector<string> getReportTypes();
    // sequence report: starts, ends, lengths, ambigs, polymers, numns
    const Rcpp::DataFrame getSequenceReport();
    // n columns: id, taxonomy split by level
    Rcpp::DataFrame getSequenceTaxonomyReport();
    const vector<string> getSamples();
    const Rcpp::DataFrame getSampleTreatmentAssignments();
    const Rcpp::DataFrame getScrapReport(string mode = "sequence");
    // type, trashCode, uniqueCount, totalCount
    const Rcpp::DataFrame getScrapSummary();
    // total abundance for each sequence
    const Rcpp::DataFrame getSequenceAbundances(bool bySample = false);
    const Rcpp::DataFrame getBinAbundances(string bin_type = "otu",
                                           bool bySample = false);
    const vector<string> getSequenceNames(vector<string> sample = nullVector,
                                          bool distinct = false);
    const vector<vector<string> > getSequenceNamesBySample(vector<string> samples = nullVector);


    const vector<string> getSequences(string sample = "");
    // 2 columns: sequence names, sequence strings
    const Rcpp::DataFrame getSequenceTable(string sample = "");
    const vector<vector<string> > getSequencesBySample(const vector<string> samples);

    const Rcpp::DataFrame getSummary(string type = "sequences",
                                     string reportType = "");

    const double getTotal(vector<string> samples = nullVector);
    const Rcpp::DataFrame getTotals(string type = "samples");
    const vector<string> getTreatments();
    const double getUniqueTotal(vector<string> samples = nullVector);

    const bool hasSample(string sample);
    const bool hasSamples(vector<string> samples = nullVector);
    const bool hasListAssignments() { return hasList; }
    const bool hasSeqs();

    void mergeBins(const vector<string>& binIDS, string reason = "merged",
                   string type = "otu");
    void mergeSequences(const vector<string>&, string reason = "merged");

    void removeBins(const vector<string>& binIDs,
                    const vector<string>& trashTags,
                    string type = "otu");
    void removeLineages(const vector<string>& taxonomies,
                        string trashTag = "contaminant");
    void removeSamples(const vector<string>& samples,
                       string reason = "remove_samples");
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
    bool hasList;

    // count table data
    AbundTable count;

    // bin table, shared / rabund
    // tables label stores type
    // ie. otu, asv, phylotype
    vector<BinTable> binTables;

    // sequence reports
    map<string, Report> reports;
    Report metadata;

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

    // 28
    template<class Archive>
    void serialize(Archive& ar) {
        ar(datasetName, processors, alignmentLength, isAligned, hasSequenceData,
           hasSequenceTaxonomy, numUnique, taxonomies, references, tableSeqs,
           starts, ends, lengths, ambigs, polymers,
           numns, names, seqs, comments, trashCodes,
           seqIndex, badAccnos, uniqueBad, hasList, count,
           binTables, reports, metadata);
    }
};

/******************************************************************************/

#endif

