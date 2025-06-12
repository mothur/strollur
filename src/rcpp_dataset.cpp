#include <Rcpp.h>
#include "../inst/include/rdataset.h"

RCPP_EXPOSED_AS(Dataset);

/******************************************************************************/
RCPP_MODULE(Dataset) {
    Rcpp::class_< Dataset >("Dataset")
    .constructor<string, int>("Create new dataset")

    // ******* properties ******* //
    .field("dataset_name", &Dataset::datasetName, "Get dataset name")
    .field("label", &Dataset::label, "Get otu label")
    .field("is_aligned", &Dataset::isAligned, "Get dataset alignment status")
    .field("num_samples", &Dataset::numSamples, "Get number of samples in dataset")
    .field("num_treatments", &Dataset::numTreatments, "Get number of treatments in dataset")
    .field("num_otus", &Dataset::numOtus, "Get number of otus in dataset")
    .field("num_unique", &Dataset::numUnique, "Get number of unique sequences")
    .field("has_contigs_data", &Dataset::hasContigsData, "Get contigs data status")
    .field("has_align_data", &Dataset::hasAlignData, "Get align data status")

    // ******* data transfer ******* //
    .method("export", &Dataset::exportDataset, "Export list containing dataset")
    .method("get_pointer", &Dataset::getPointer, "get pointer to self")

    // ******* modifiers ******* //
    .method("add_sequences", &Dataset::addSequences, "Add sequences to dataset")
    .method("add_align_report", &Dataset::addAlignReport,
    "Add alignment report to dataset")
    .method("add_contigs_report", &Dataset::addContigsReport,
    "Add contigs report to dataset")

    .method("assign_otus", &Dataset::assignOTUAbundance, "Add otu assignments to dataset")
    .method("assign_sequence_abundance", &Dataset::assignSequenceAbundance,
    "Set sequence abundance and optionally assign sample and treatment data")
    .method("clear", &Dataset::clear, "Clear dataset")
    .method("merge_sequences", &Dataset::mergeSequences,
    "Merge sequences")
    .method("remove_sequences", &Dataset::removeSequences, "Remove sequences")
    .method("set_sequences", &Dataset::setSequences, "Set sequence strings")
    .method("set_abundances", &Dataset::setAbundances, "Set sequence abundances")

    // ******* getters ******* //
    // abundances for seq broken down by sample
    .method("get_abundances", &Dataset::getAbundances, "Get abundances for sequence")
    .method("get_sequence_abundances", &Dataset::getSequenceAbundances,
    "Get total abundance for each sequence")
    .method("get_count_matrix", &Dataset::getSeqsAbundsBySample,
    "Get dataset's sequence's abundances by sample")
    .method("get_samples", &Dataset::getSamples, "Get samples in dataset")
    .method("get_sample_totals", &Dataset::getSampleTotals,
    "Get sample totals")
    .method("has_sample", &Dataset::hasSample,
    "Determine if dataset contains sample")
    .method("get_treatments", &Dataset::getTreatments, "Get treatments in dataset")
    .method("get_treatment_totals", &Dataset::getTreatmentTotals,
    "Get treatment totals")

    .method("get_sequence_report", &Dataset::getSequenceReport,
    "Get sequence report data: starts, ends, lengths, ambigs, homopolymers, numns")
    .method("get_sequence_summary", &Dataset::getSequenceSummary,
    "Get sequence summary report")
    .method("get_align_report", &Dataset::getAlignReport,
    "Get align report data: search_score, sim_score, longest_insert")
    .method("get_contigs_report", &Dataset::getContigsReport,
    "Get contigs report data: olengths, ostarts, oends, mismatches, ee")
    .method("get_scrap_report", &Dataset::getScrapReport,
    "Get scrap report data: id, trash_code")

    // 3 column table: id, sample, abundance
    .method("get_sequence_abundance_table", &Dataset::getSequenceAbundanceTable,
    "Get sequence abundance table")

    // can pass optional sample name
    .method("get_abundance", &Dataset::getAbundance, "Get abundance of sequence")
    .method("get_names", &Dataset::getNames,
    "Get names of sequences in dataset")
    .method("get_sequences", &Dataset::getSequences, "Get sequences in dataset")
    .method("get_total", &Dataset::getTotal,
    "Get total number of sequences")
    .method("get_unique_total", &Dataset::getUniqueTotal,
    "Get total number of unique sequences")

    // parse by sample options
    .method("get_names_by_sample", &Dataset::getNamesBySample,
    "Get names of sequences parsed by sample")
    .method("get_sequences_by_sample", &Dataset::getSequencesBySample,
    "Get sequences parsed by sample")
    ;
}
/******************************************************************************/

