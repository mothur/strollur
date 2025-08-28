#include <Rcpp.h>
#include "../inst/include/rdataset.h"

RCPP_EXPOSED_AS(Dataset);

/******************************************************************************/
RCPP_MODULE(Dataset) {
    Rcpp::class_< Dataset >("Dataset")
    .constructor<string, int>("Create new dataset")
    .constructor<const Dataset&>("Exposes a copy constructor")

    // ******* exposed properties ******* //
    .field("dataset_name", &Dataset::datasetName, "Get dataset name")
    .field("processors", &Dataset::processors, "Get number of threads used")
    .field("is_aligned", &Dataset::isAligned, "Get dataset alignment status")
    .field("num_unique", &Dataset::numUnique, "Get number of unique sequences")

    // ******* exposed functions ******* //
    .method("add_sequences", &Dataset::addSequences, "Add sequences to dataset")
    .method("assign_bins", &Dataset::assignBins, "Add bin assignments to dataset")
    .method("assign_bin_taxonomy", &Dataset::assignBinTaxonomy, "Assign bin classification")
    .method("assign_sequence_abundance", &Dataset::assignSequenceAbundance, "Set sequence abundance and optionally assign sample and treatment data")
    .method("assign_sequence_taxonomy", &Dataset::assignSequenceTaxonomy, "Assign sequence classification")
    .method("assign_treatments", &Dataset::assignTreatments, "Assign samples to treatments")

    .method("clear", &Dataset::clear, "Clear dataset")
    .method("load", &Dataset::loadFromSerialized, "Load dataset from serialized data")
    .method("export", &Dataset::exportDataset, "Export list containing dataset")
    .method("get_pointer", &Dataset::getPointer, "Get pointer to dataset c++ class")

    .method("get_abundance", &Dataset::getSequenceAbundances, "Get total abundance for each sequence")
    .method("get_abundances", &Dataset::getSeqsAbundsBySample, "Get abundances for sequences by sample")
    .method("get_list", &Dataset::getList, "Get data.frame containing sequence bin assignments")
    .method("get_names", &Dataset::getNames, "Get names of sequences in dataset")
    .method("get_names_by_sample", &Dataset::getNamesBySample, "Get names of sequences parsed by sample")
    .method("get_num_bins", &Dataset::getNumBins, "Get number of bins")
    .method("get_num_samples", &Dataset::getNumSamples, "Get number of samples in dataset")
    .method("get_num_treatments", &Dataset::getNumTreatments, "Get number of treatments in dataset")
    .method("get_bin", &Dataset::getBin, "Get bin containing sequence names")
    .method("get_bin_abundance", &Dataset::getBinAbundance, "Get abundance of bin")
    .method("get_bin_abundances", &Dataset::getBinAbundances, "Get abundance of bin parsed by sample")
    .method("get_bin_taxonomy_report", &Dataset::getBinTaxonomyReport, "Get data.frame containing bin taxonomy data")
    .method("get_rabund", &Dataset::getRAbund, "Get data.frame containing bin abundance data")
    .method("get_samples", &Dataset::getSamples, "Get samples in dataset")
    .method("get_sample_totals", &Dataset::getSampleTotals, "Get sample totals")
    .method("get_scrap_report", &Dataset::getScrapReport, "Get scrap report data: id, trash_code")
    .method("get_sequences", &Dataset::getSequences, "Get sequences in dataset")
    .method("get_sequences_by_sample", &Dataset::getSequencesBySample, "Get sequences parsed by sample")
    // 4 column table: id, abundance, sample(optional), treatment(optional)
    .method("get_sequence_abundance_table", &Dataset::getSequenceAbundanceTable, "Get sequence abundance table")
    .method("get_sequence_report", &Dataset::getSequenceReport, "Get sequence report data: starts, ends, lengths, ambigs, homopolymers, numns")
    .method("get_sequence_summary", &Dataset::getSequenceSummary, "Get sequence summary report")
    .method("get_sequence_taxonomy_report", &Dataset::getSequenceTaxonomyReport, "Get data.frame containing sequence taxonomy data")
    .method("get_shared", &Dataset::getShared, "Get data.frame containing bin abundance data by sample")
    .method("get_total", &Dataset::getTotal, "Get total number of sequences")
    .method("get_treatments", &Dataset::getTreatments, "Get treatments in dataset")
    .method("get_treatment_totals", &Dataset::getTreatmentTotals, "Get treatment totals")
    .method("get_unique_total", &Dataset::getUniqueTotal, "Get total number of unique sequences")

    .method("has_sample", &Dataset::hasSample, "Determine if dataset contains sample")
    .method("has_sequence_strings", &Dataset::hasSeqs, "Get sequence status")

    .method("merge_bins", &Dataset::mergeBins, "Merge bins")
    .method("merge_sequences", &Dataset::mergeSequences, "Merge sequences")

    .method("remove_lineages", &Dataset::removeLineages, "Get sequence contaminants")
    .method("remove_bins", &Dataset::removeBins, "Remove bins")
    .method("remove_samples", &Dataset::removeSamples, "Remove samples from dataset")
    .method("remove_sequences", &Dataset::removeSequences, "Remove sequences")

    .method("set_abundance", &Dataset::setAbundance, "Set sequence abundances for datasets without sample data")
    .method("set_abundances", &Dataset::setAbundances, "Set sequence abundances for datasets with multiple samples")
    .method("set_bin_abundance", &Dataset::setBinAbundance, "Set abundance of bins without samples")
    .method("set_bin_abundances", &Dataset::setBinAbundances, "Set abundance of bins without samples")
    .method("set_sequences", &Dataset::setSequences, "Set sequence strings")
    .method("serialize_dataset", &Dataset::serializeDataset, "Serialize Dataset class into binary string")
    ;
}
/******************************************************************************/
