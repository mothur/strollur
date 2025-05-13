#include <Rcpp.h>
#include "../inst/include/rdataset.h"

RCPP_EXPOSED_AS(Dataset);

/******************************************************************************/
RCPP_MODULE(Dataset) {
    Rcpp::class_< Dataset >("Dataset")
    .constructor<string>("Create new dataset")

    // ******* properties ******* //
    .field("dataset_name", &Dataset::datasetName, "Get dataset name")
    .field("is_aligned", &Dataset::isAligned, "Get dataset alignment status")
    .field("num_groups", &Dataset::numGroups, "Get number of groups in dataset")
    .field("num_unique", &Dataset::numUnique, "Get number of unique sequences")

    // ******* data transfer ******* //
    .method("print", &Dataset::print, "Summary of dataset")
    .method("export", &Dataset::exportDataset, "Export list containing dataset")
    .method("get_pointer", &Dataset::getPointer, "get pointer to self")

    // ******* modifiers ******* //
    .method("add_seqs", &Dataset::addSeqs, "Add sequences to dataset")
    .method("add_align_report", &Dataset::addAlignReport,
    "Add alignment report to dataset")
    .method("add_contigs_report", &Dataset::addContigsReport,
    "Add contigs report to dataset")

    .method("assign_sample_abundance", &Dataset::assignSampleAbundance,
    "Assign sequence sample abundance")
    .method("clear", &Dataset::clear, "Clear dataset")
    .method("merge_seqs", &Dataset::mergeSeqs,
    "Merge sequences")
    .method("reinstate_seqs", &Dataset::reinstateSeqs, "Reinstate sequences")
    .method("remove_seqs", &Dataset::removeSeqs, "Remove sequences")
    .method("set_seqs", &Dataset::setSeqs, "Set sequence strings")
    .method("set_abunds", &Dataset::setAbundances, "Set sequence abundances")

    // ******* getters ******* //
    // abundances for seq broken down by sample
    .method("get_abunds", &Dataset::getAbunds, "Get abundances for sequence")
    .method("get_seqs_abunds", &Dataset::getSeqsAbunds,
    "Get total abundance for each sequence")
    .method("get_count_matrix", &Dataset::getSeqsAbundsBySample,
    "Get dataset's sequence's abundances by sample")
    .method("get_groups", &Dataset::getGroups, "Get groups in dataset")
    .method("get_group_totals", &Dataset::getGroupTotals,
    "Get group totals")
    .method("has_group", &Dataset::hasGroup,
    "Determine if dataset contains sample")

    .method("get_fasta_report", &Dataset::getFastaReport,
    "Get fasta summary data: starts, ends, lengths, ambigs, homopolymers, numns")
    .method("get_align_report", &Dataset::getAlignReport,
    "Get align summary data: search_score, sim_score, longest_insert")
    .method("get_contigs_report", &Dataset::getContigsReport,
    "Get contigs sumary data: olengths, ostarts, oends, mismatches, ee")
    // 3 column table: id, group, abundance
    .method("get_sequence_abundance_table", &Dataset::getSequenceAbundanceTable,
    "Get sequence abundance table")

    // can pass optional group name
    .method("get_abund", &Dataset::getAbund, "Get abundance of sequence")
    .method("get_names", &Dataset::getNames,
    "Get names of sequences in dataset")
    .method("get_seqs", &Dataset::getSeqs, "Get sequences in dataset")
    .method("get_total", &Dataset::getTotal,
    "Get total number of sequences")

    // parse by sample options
    .method("get_names_by_sample", &Dataset::getNamesBySample,
    "Get names of sequences parsed by sample")
    .method("get_seqs_by_sample", &Dataset::getSeqsBySample,
    "Get sequences parsed by sample")
    ;
}
/******************************************************************************/

