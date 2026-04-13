#ifndef INTERNAL_DEVELOPMENT_H_
#define INTERNAL_DEVELOPMENT_H_

#include <Rcpp.h>
#include "../inst/include/strollur.h"
#include "dataset.h"

/******************************************************************************/
SEXP xint_fill_required_parameters(const Rcpp::DataFrame& df,
                                   const string& given_column_name,
                                   const string& type = "string");

SEXP xint_fill_optional_parameters(const Rcpp::DataFrame& df,
                                   const string& default_column_name,
                                   const string& given_column_name,
                                   const string& type = "string");
/******************************************************************************/
//' @title Get a data.frame containing the requested abundance data
//' @name xdev_abundance
//' @rdname xdev_abundance
//' @description
//' Get a table containing the requested abundance data in a
//' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param type, string containing the type of data you want the number of.
//' Options include: "sequences", "bins".
//' Default = "sequences".
//'
//' @param bin_type, string containing the bin type you would like the number of
//' bins for. Default = "otu".
//'
//' @param by_sample, Boolean. When by_sample is TRUE, the abundance data will
//' be parsed by sample. Default = FALSE.
//'
//' @examples
//'
//' miseq <- miseq_sop_example()
//'
//' # To the total abundance for each sequence
//' xdev_abundance(data = miseq, type = "sequences")
//'
//' # To the total abundance for each sequence parsed by sample
//' xdev_abundance(data = miseq, type = "sequences", by_sample = TRUE)
//'
//' # To the total abundance for each "otu" bin
//' xdev_abundance(data = miseq, type = "bins", bin_type = "otu")
//'
//' # To the total abundance for each "otu" bin parsed by sample
//' xdev_abundance(data = miseq, type = "bins", bin_type = "otu", by_sample = TRUE)
//'
//' # To the total abundance for each "asv" bin
//' xdev_abundance(data = miseq, type = "bins", bin_type = "asv")
//'
//' # To the total abundance for each "asv" bin parsed by sample
//' xdev_abundance(data = miseq, type = "bins", bin_type = "asv", by_sample = TRUE)
//'
//' # To the total abundance of each sample
//' xdev_abundance(data = miseq, type = "samples")
//'
//' # To the total abundance of each treatment
//' xdev_abundance(data = miseq, type = "treatments")
//'
//' @return data.frame
//' @export
//[[Rcpp::export]]
Rcpp::DataFrame xdev_abundance(const Rcpp::Environment& data,
                               const string& type = "sequences",
                               const string& bin_type = "otu",
                               bool by_sample = false);
/******************************************************************************/
//' @title Add resource references
//' @name xdev_add_references
//' @rdname xdev_add_references
//' @description
//' Add resource references to a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing reference_names, reference_versions
//' (optional), reference_usages (optional), reference_notes (optional), and
//' reference_urls (optional).
//'
//' @param reference_name, a string containing the name of the column in 'table'
//' that contains the reference names. Default column name is 'reference_names'.
//' @param reference_version, a string containing the name of the column in
//' 'table' that contains the reference versions. Default column name is
//' 'reference_versions'.
//' @param reference_usage, a string containing the name of the column in
//' 'table' that contains the reference usages. Default column name is
//'  reference_usages'.
//' @param reference_note, a string containing the name of the column in
//' 'table' that contains the reference notes. Default column name is
//'  reference_notes'.
//' @param reference_url, a string containing the name of the column in
//' 'table' that contains the reference urls. Default column name is
//'  reference_urls'.
//'
//' @param verbose, a boolean whether or not you want progress messages.
//' Default = TRUE.
//'
//' @examples
//'
//' data <- new_dataset("just for fun", 2)
//' reference_table <- readr::read_csv(strollur_example("references.csv"),
//'                              col_names = TRUE, show_col_types = FALSE)
//' xdev_add_references(data, reference_table)
//'
//' @return double containing the number of references added
//' @export
//[[Rcpp::export]]
double xdev_add_references(const Rcpp::Environment& data,
                       const Rcpp::DataFrame& table,
                       const string& reference_name = "reference_names",
                       const string& reference_version = "reference_versions",
                       const string& reference_usage = "reference_usages",
                       const string& reference_note = "reference_notes",
                       const string& reference_url = "reference_urls",
                       bool verbose = true);
/******************************************************************************/
//' @title Add a report to a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @name xdev_add_report
//' @rdname xdev_add_report
//' @description
//' Add a report to a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing your report.
//'
//' @param type, a string containing the type of report. Options include:
//' "metadata" and custom report tags. Default = "metadata".
//'
//' @param sequence_name, a string containing the name of the column in 'table'
//' that contains the sequence names. This is used for custom reports, metadata
//' does not require a sequence_name column. Default column name is 'sequence_names'.
//'
//' @param verbose, a boolean whether or not you want progress messages.
//' Default = TRUE.
//'
//' @examples
//'
//' # To add a custom report including your contigs assembly data
//'
//' data <- new_dataset("just for fun", 2)
//' contigs_report <- readRDS(strollur_example("miseq_contigs_report.rds"))
//'
//' xdev_add_report(data, contigs_report, "contigs_report", "Name")
//'
//' # To add metadata related to your study
//'
//' metadata <- readRDS(strollur_example("miseq_metadata.rds"))
//'
//' xdev_add_report(data, metadata, "metadata")
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_add_report(const Rcpp::Environment& data,
                 Rcpp::DataFrame table,
                 const string& type = "metadata",
                 const string& sequence_name = "sequence_names",
                 bool verbose = true);
/******************************************************************************/
//' @title xdev_add_sequences
//' @description
//' Add sequence data to a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing names, sequences(optional) and
//' comments(optional).
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//' @param sequence_name, a string containing the name of the column in 'table'
//' that contains the sequence names. Default column name is 'sequence_names'.
//' @param sequence, a string containing the name of the column in 'table' that
//' contains the sequence nucleotide strings. Default column name is
//' 'sequences'.
//' @param comment, a string containing the name of the column in
//' 'table' that contains the sequence comments. Default column name is
//' 'comments'.
//'
//' @param verbose, a boolean whether or not you want progress messages.
//' Default = TRUE.
//'
//' @examples
//'
//'  data <- new_dataset("miseq_sop", 2)
//'  fasta_data <- read_fasta(strollur_example("final.fasta.gz"))
//'  xdev_add_sequences(data, fasta_data)
//'
//' # With the additional parameters to add information about the reference
//'
//'  data <- new_dataset("miseq_sop", 2)
//'  fasta_data <- read_fasta(strollur_example("final.fasta.gz"))
//'
//'  xdev_add_sequences(data, fasta_data,
//'                new_reference("silva.bacteria.fasta",
//'                "1.38.1",
//'                "alignment by mothur2 v1.0 using default options",
//'                "https://mothur.org/wiki/silva_reference_files/"))
//'
//' # You can also add references using the 'add_references' function.
//'
//' @return double containing the number of sequences added
//' @export
//[[Rcpp::export]]
double xdev_add_sequences(const Rcpp::Environment& data,
                      const Rcpp::DataFrame& table,
                      Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                      const string& sequence_name = "sequence_names",
                      const string& sequence = "sequences",
                      const string& comment = "comments",
                      bool verbose = true);
/******************************************************************************/
//' @title xdev_assign_bins
//' @description
//' Add bin assignments to a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing bin_data assignments
//' @param bin_type a string indicating the type of bin assignments. Default "otu".
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//'
//' @param bin_name, a string containing the name of the column in 'table' that
//' contains the bin names. Default column name is 'bin_names'.
//' @param abundance, a string containing the name of the column in 'table'
//' that contains the bin abundances. Default column name is 'abundances'. Note:
//'  You must provide either abundances or sequence_names in the table.
//' @param sample, a string containing the name of the column in 'table' that
//' contains the sample names for datasets where the abundances are broken down
//' by sample. Default column name is 'samples'.
//' @param sequence_name, a string containing the name of the column in 'table'
//' that contains the sequence names. Default column name is 'sequence_names'.
//' Note: You must provide either abundances or sequence_names in the table.
//'
//' @param verbose, a boolean whether or not you want progress messages.
//' Default = TRUE.
//'
//' @examples
//'
//'   # To assign sequences to bins:
//'
//'   data <- new_dataset(dataset_name = "miseq_sop")
//'   otu_data <- read_mothur_list(list = strollur_example("final.opti_mcc.list.gz"))
//'
//'   xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
//'
//'   # To add abundance only bin assignments:
//'
//'   data <- new_dataset(dataset_name = "miseq_sop")
//'   otu_data <- read_mothur_rabund(rabund = strollur_example("final.opti_mcc.rabund"))
//'
//'   xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
//'
//'   # To add abundance bin assignments parsed by sample:
//'
//'   data <- new_dataset(dataset_name = "miseq_sop")
//'   otu_data <- readRDS(strollur_example("miseq_shared_otu.rds"))
//'
//'   xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
//'
//' @return double containing the number of bins assigned
//' @export
//[[Rcpp::export]]
double xdev_assign_bins(const Rcpp::Environment& data,
                    const Rcpp::DataFrame& table,
                    const string& bin_type = "otu",
                    Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                    const string& bin_name = "bin_names",
                    const string& abundance = "abundances",
                    const string& sample = "samples",
                    const string& sequence_name = "sequence_names",
                    bool verbose = true);
/******************************************************************************/
//' @title xdev_assign_bin_representative_sequences
//' @description
//' Assign representative sequences to bins.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing bin representative assignments
//'
//' @param bin_type a string indicating the type of bin assignments.
//' Default "otu".
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//'
//' @param bin_name, a string containing the name of the column in 'table' that
//' contains the bin names. Default column name is 'bin_names'.
//' @param sequence_name a string containing the name of the column in 'table' that
//' contains the bin names. Default column name is 'sequence_names'.
//'
//' @param verbose, a boolean whether or not you want progress messages.
//' Default = TRUE.
//'
//' @examples
//'
//'   miseq <- miseq_sop_example()
//'
//'   bin_reps <- readRDS(strollur_example(
//'          "miseq_representative_sequences.rds"))
//'
//'   xdev_assign_bin_representative_sequences(data = miseq,
//'                                       table = bin_reps,
//'                                       bin_type = "otu")
//'
//' @return double containing the number of representative sequences assigned
//' @export
//[[Rcpp::export]]
double xdev_assign_bin_representative_sequences(const Rcpp::Environment& data,
                                            const Rcpp::DataFrame& table,
                                            const string& bin_type = "otu",
                                            Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                                            const string& bin_name = "bin_names",
                                            const string& sequence_name = "sequence_names",
                                            bool verbose = true);
/******************************************************************************/
//' @title xdev_assign_bin_taxonomy
//' @description
//' Assign bin classifications to a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' Note, if you assign sequence taxonomies and assign bins, 'Dataset' will find
//'  the concensus taxonomy for each bin for you.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing bin taxonomy assignments
//' @param bin_type a string indicating the type of bin assignments. Default "otu".
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//'
//' @param bin_name, a string containing the name of the column in 'table' that
//' contains the bin names. Default column name is 'bin_names'.
//' @param taxonomy, a string containing the name of the column in 'table' that
//' contains the bin taxonomies. Default column name is 'taxonomies'.
//'
//' @param verbose, a boolean whether or not you want progress messages.
//' Default = TRUE.
//'
//' @examples
//'
//' otu_data <- read_mothur_cons_taxonomy(strollur_example(
//'                         "final.cons.taxonomy"))
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//'
//' # assign otu abundances
//' xdev_assign_bins(data = data, table = otu_data, bin_type = "otu")
//'
//' # assign otu classifications
//' xdev_assign_bin_taxonomy(data = data, table = otu_data,
//'                          bin_type = "otu")
//'
//' @return double containing the number of bins assigned
//' @export
//[[Rcpp::export]]
double xdev_assign_bin_taxonomy(const Rcpp::Environment& data,
                            const Rcpp::DataFrame& table,
                            const string& bin_type = "otu",
                            Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                            const string& bin_name = "bin_names",
                            const string& taxonomy = "taxonomies",
                            bool verbose = true);
/******************************************************************************/
//' @title xdev_assign_sequence_taxonomy
//' @description
//' Assign sequence classifications to a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' Note, if you assign sequence taxonomies and assign bins, 'Dataset' will find
//'  the consensus taxonomy for each bin for you.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing sequence taxonomy assignments
//'
//' @param reference, a list created by the function [new_reference]. Optional.
//'
//' @param sequence_name, a string containing the name of the column in 'table'
//' that contains the sequence names. Default column name is 'sequence_names'.
//' @param taxonomy, a string containing the name of the column in 'table' that
//' contains the sequence taxonomies. Default column name is 'taxonomies'.
//'
//' @param verbose, a boolean whether or not you want progress messages.
//' Default = TRUE.
//'
//' @examples
//'
//' sequence_classifications <- read_mothur_taxonomy(strollur_example(
//'                         "final.taxonomy.gz"))
//'
//' data <- new_dataset("my_dataset", 2)
//'
//' xdev_assign_sequence_taxonomy(data, sequence_classifications)
//'
//' # With the reference parameter you can add information about the reference
//' # you used to classify your sequences. You can also add references using the
//' # 'add_references' function.
//'
//' reference <- new_reference("trainset9_032012.pds.zip", "9_032012",
//'               "classification by mothur2 v1.0 using default options", "",
//' "https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip")
//'
//' xdev_assign_sequence_taxonomy(data, sequence_classifications, reference)
//'
//' @return double containing the number of sequence assigned
//' @export
//[[Rcpp::export]]
double xdev_assign_sequence_taxonomy(const Rcpp::Environment& data,
                                 const Rcpp::DataFrame& table,
                                 Rcpp::Nullable<Rcpp::List> reference = R_NilValue,
                                 const string& sequence_name = "sequence_names",
                                 const string& taxonomy = "taxonomies",
                                 bool verbose = true);

/******************************************************************************/
//' @title xdev_assign_sequence_abundance
//' @description
//' Assign sequence abundance and optionally assign sample and treatment data to
//'  a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing sequence abundance assignments
//'
//' @param sequence_name, a string containing the name of the column in 'table'
//'  that contains the sequence names. Default column name is 'sequence_names'.
//' @param abundance, a string containing the name of the column in 'table'
//'  that contains the sequence abundances. Default column name is 'abundances'.
//' @param sample, a string containing the name of the column in 'table'
//' that contains the sequence samples. Default column name is 'samples'.
//' (Optional)
//' @param treatment, a string containing the name of the column in 'table'
//' that contains the sequence treatments. Default column name is 'treatments'.
//'
//' @param verbose, a boolean whether or not you want progress messages.
//' Default = TRUE.
//'
//' @examples
//'
//' data <- new_dataset("my_dataset")
//' sequence_abundance <- readRDS(strollur_example("miseq_abundance_by_sample.rds"))
//'
//' xdev_assign_sequence_abundance(data = data, table = sequence_abundance)
//'
//' @return double containing the number of sequences assigned
//' @export
//[[Rcpp::export]]
double xdev_assign_sequence_abundance(const Rcpp::Environment& data,
                                  const Rcpp::DataFrame& table,
                                  const string& sequence_name = "sequence_names",
                                  const string& abundance = "abundances",
                                  const string& sample = "samples",
                                  const string& treatment = "treatments",
                                  bool verbose = true);
/******************************************************************************/
//' @title xdev_assign_treatments
//' @description
//' Assign samples to treatments in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param table, a data.frame containing sample treatment assignments
//'
//' @param sample, a string containing the name of the column in 'table'
//' that contains the samples. Default column name is 'samples'.
//' @param treatment, a string containing the name of the column in 'table'
//' that contains the treatments. Default column name is 'treatments'.
//'
//' @param verbose, a boolean indicating whether or not you want progress
//'  messages. Default = TRUE.
//'
//' @examples
//'
//' data <- new_dataset("my_dataset")
//' sequence_abundance <- readRDS(strollur_example("miseq_abundance_by_sample.rds"))
//'
//' xdev_assign_sequence_abundance(data, sequence_abundance)
//'
//' sample_assignments <- readRDS(strollur_example("miseq_sample_design.rds"))
//'
//' xdev_assign_treatments(data, sample_assignments)
//'
//' @return double containing the number of samples assigned to treatments
//' @export
//[[Rcpp::export]]
double xdev_assign_treatments(const Rcpp::Environment& data,
                          const Rcpp::DataFrame& table,
                          const string& sample = "samples",
                          const string& treatment = "treatments",
                          bool verbose = true);
/******************************************************************************/
//' @title xdev_count
//' @description
//' Find the number of sequences, samples, treatments or bins of a given type in
//' a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param type, string containing the type of data you want the number of.
//' Options include: "sequences", "samples", "treatments", "bins" and
//' "references". Default = "sequences".
//'
//' @param bin_type, string containing the bin type you would like the number of
//' bins for. Default = "otu".
//'
//' @param samples, vector of strings. samples is only used when 'type' =
//' "sequences" or 'type' = "bins" . samples should contain the names of the
//' samples you want the count for. Default = NULL.
//'
//' @param distinct, Boolean. distinct is used when 'type' =
//' "sequences" or 'type' = "bins". When 'type' = "sequences" and distinct is
//' TRUE the number of unique sequences is returned. When 'type' = "sequences"
//' and distinct is FALSE total number of sequences is returned. This can also
//' be combined with samples to find the number of unique sequences found only
//' in a given set of samples, or to find the total number of sequences in a
//' given set of samples.
//' When 'type' = "bins", you can set distinct = TRUE to return the number of
//' bins that ONLY contain sequences from the given samples. When distinct is
//' FALSE the count returned contains bins with sequences from a given samples,
//' but those bins may also contain other samples.
//' Default = FALSE.
//' @examples
//'
//' miseq <- miseq_sop_example()
//'
//' # To get the total number of sequences
//' xdev_count(data = miseq, type = "sequences")
//'
//' # To get number of unique sequences
//' xdev_count(data = miseq, type = "sequences", distinct = TRUE)
//'
//' # To get number of unique sequences from samples 'F3D0' and 'F3D1'
//' # Note these sequences will be present in both samples but may be
//' # be present in other samples as well
//' xdev_count(data = miseq, type = "sequences", samples = c("F3D0", "F3D1"))
//'
//' # To get number of unique sequences exclusive to samples 'F3D0' and 'F3D1'
//' # Note these sequences are present in both samples and NOT present in
//' # other samples
//' xdev_count(data = miseq, type = "sequences", samples = c("F3D0", "F3D1"),
//' distinct = TRUE)
//'
//' # To get the number of samples in the dataset
//' xdev_count(data = miseq, type = "samples")
//'
//' # To get the number of treatments in the dataset
//' xdev_count(data = miseq, type = "treatments")
//'
//' # To get the number of "otu" bins in the dataset
//' xdev_count(data = miseq, type = "bins", bin_type = "otu")
//'
//' # To get the number of "asv" bins in the dataset
//' xdev_count(data = miseq, type = "bins", bin_type = "asv")
//'
//' # To get the number of "phylotype" bins in the dataset
//' xdev_count(data = miseq, type = "bins", bin_type = "phylotype")
//'
//' # To get number of bins from samples 'F3D0' and 'F3D1'
//' # Note these bins will have sequences from both samples but there may be
//' # other samples present as well
//' xdev_count(data = miseq, type = "bins", samples = c("F3D0", "F3D1"))
//'
//' # To get number of bins unique to samples 'F3D0' and 'F3D1'
//' # Note these bins will have sequences from both samples and NO other samples
//' # will be present in the bins.
//' xdev_count(data = miseq, type = "bins", samples = c("F3D0", "F3D1"),
//' distinct = TRUE)
//'
//' @return double
//' @export
//[[Rcpp::export]]
double xdev_count(const Rcpp::Environment& data,
            const string& type = "sequences",
            const string& bin_type = "otu",
            Rcpp::Nullable<Rcpp::List> samples = R_NilValue,
            bool distinct = false);
/******************************************************************************/
//' @title xdev_get_abundances_by_sample
//' @description
//' Get the sequence abundance data in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object parsed by sample
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param samples a vector of strings containing the names of the samples you
//' would like sequence names for. By default all samples are included.
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' # To get the sequence names parsed by sample
//' abunds <- xdev_get_abundances_by_sample(data)
//'
//' @return 2D vector of float containing data requested parsed by sample.
//' @export
//[[Rcpp::export]]
vector<vector<float> > xdev_get_abundances_by_sample(const Rcpp::Environment& data,
                                            const Rcpp::CharacterVector& samples = Rcpp::CharacterVector::create());
/******************************************************************************/
//' @title xdev_get_list_vector
//' @description
//' Get vector of strings containing the sequences bin data
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param type a string indicating the type of bin assignments. Default "otu".
//' @examples
//'
//' data <- miseq_sop_example()
//' xdev_get_list_vector(data)
//'
//' @return vector of strings containing the names of the sequences in each bin
//' separated by commas
//' @export
//[[Rcpp::export]]
vector<string> xdev_get_list_vector(const Rcpp::Environment& data,
                                    const string& type = "otu");
/******************************************************************************/
//' @title xdev_get_by_sample
//' @description
//' Get the requested data in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object parsed by sample
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param type, string containing the type of data you want the totals of.
//' Options include: "sequence_names", "sequences". Default = "sequence_names".
//'
//' @param samples a vector of strings containing the names of the samples you
//' would like sequence names for. By default all samples are included.
//'
//' @param degap a logical. Default = FALSE. When degap = `TRUE`, all gap
//' characters will be removed from the sequences.
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' # To get the sequence names parsed by sample
//' xdev_get_by_sample(data, "sequence_names")
//'
//' # To get the sequence nucleotide strings parsed by sample
//' parsed_sequences <- xdev_get_by_sample(data, "sequences")
//'
//' @return 2D vector of strings ([num_seqs][num_samples]) containing data
//' requested parsed by sample.
//' @export
//[[Rcpp::export]]
vector<vector<string> > xdev_get_by_sample(const Rcpp::Environment& data,
                                      const string& type = "sequence_names",
                                      const Rcpp::CharacterVector& samples = Rcpp::CharacterVector::create(),
                                      bool degap = false);
/******************************************************************************/
//' @title xdev_get_sequences
//' @description
//' Get the nucleotide strings for each sequence in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param sample, a string containing the name of the sample you
//' would like sequence names for. For all samples in dataset, sample = "".
//' @param degap, a logical. Default = FALSE. When degap = `TRUE`, all gap
//' characters ('-', '.') will be removed from the sequences.
//' @examples
//'
//'  data <- miseq_sop_example()
//'  xdev_get_sequences(data)
//'
//' @return vector of string containing nucleotide strings of the sequences in
//' a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @export
//[[Rcpp::export]]
vector<string> xdev_get_sequences(const Rcpp::Environment& data, const string& sample = "", bool degap = false);
/******************************************************************************/
//' @title xdev_has_sequence_taxonomy
//' @description
//' Determine if a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object has sequence taxonomy assignments
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @examples
//'
//'  data <- miseq_sop_example()
//'  xdev_has_sequence_taxonomy(data)
//'
//' @return boolean
//' @export
//[[Rcpp::export]]
bool xdev_has_sequence_taxonomy(const Rcpp::Environment& data);
/******************************************************************************/
// ******** merging *********

//' @title xdev_merge_bins
//' @description
//' Designed with package integration in mind, the merge bins function allows
//' you to merge bins in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param bin_names, a vector of strings containing the names of the bins you
//' would like merge. The resulting merged bin will be stored in the first
//' bin_id in the vector.
//' @param reason, a string indicating why you are merging bins. Default =
//' "merged".
//' @param bin_type, a string indicating the type of bin clusters.
//'  Default = "otu"
//'
//' @examples
//'
//'  data <- miseq_sop_example()
//'
//'  # to merge otu5 and otu6
//'
//'  bins_to_merge <- c("Otu005", "Otu006")
//'
//'  xdev_merge_bins(data = data, bin_names = bins_to_merge)
//'
//'  # If you look at the scrap report, you will see Otu006 with the trash code
//'  # set to "merged".
//'
//'  report(data = data, type = "bin_scrap")
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_merge_bins(const Rcpp::Environment& data, const vector<string>& bin_names,
                     const string& reason = "merged", const string& bin_type = "otu");
/******************************************************************************/
//' @title xdev_merge_sequences
//' @description
//' Designed with package integration in mind, the merge sequences function
//' allows you to merge sequences in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param sequence_names, a vector of strings containing the names of the
//' sequences you would like merge. The resulting merged sequence will be stored
//' in the first sequence name in the vector.
//' @param reason a string indicating why you are merging sequences.
//' Default = "merged"
//'
//' @examples
//'
//' sequence_names <- c("seq1", "seq2", "seq3", "seq3",
//'                "seq4", "seq4", "seq5", "seq6",
//'                "seq7", "seq8", "seq9", "seq9",
//'                "seq10", "seq10", "seq10", "seq10")
//'
//' samples <- c("sample1", "sample2", "sample4", "sample5",
//'              "sample1", "sample2", "sample1", "sample1",
//'              "sample2", "sample4", "sample4", "sample5",
//'              "sample1", "sample3", "sample5", "sample6")
//'
//' abundances <- c(10, 10, 5, 5, 5, 5,
//'                 10, 10, 10, 10, 5, 5,
//'                 1, 2, 3, 4)
//'
//' data <- new_dataset("my_data")
//'
//'
//' assign(data = data,
//'        table = data.frame(sequence_names = sequence_names,
//'                           abundances = abundances,
//'                           samples = samples),
//'        type = "sequence_abundance")
//'
//' # For the sake of example let's merge the first 3 sequences.
//'
//' seqs_to_merge <- c("seq1", "seq2", "seq3")
//'
//' xdev_merge_sequences(data = data, sequence_names = seqs_to_merge)
//'
//' # If you look at the scrap report, you will see the second two sequence
//' # names, listed with the trash code set to "merged".
//'
//' report(data = data, type = "sequence_scrap")
//'
//' # You can see from the get_num_sequences function that the merged sequence's
//' # abundances are added to the first sequence.
//'
//' count(data = data, type = "sequences")
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_merge_sequences(const Rcpp::Environment& data, const vector<string>& sequence_names,
                          const string& reason = "merged");

/******************************************************************************/
//' @title xdev_names
//' @description
//' Get the names of a given type of data in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param type, string containing the type of data you would like. Options
//' include: "dataset", "sequences", "bins", "samples", "treatments", "reports".
//' Default = "sequences".
//'
//' @param bin_type, string containing the bin type you would like the names
//' for. Default = "otu".
//'
//' @param samples, vector of strings. samples is only used when 'type' =
//' "sequences" or 'type' = "bins" . samples should contain the names of the
//' samples you want names for. Default = NULL.
//'
//' @param distinct, Boolean. distinct is used when 'type' =
//' "sequences" or 'type' = "bins" and the samples parameter is used. The
//' distinct parameter allows you to get the names that are unique to a given
//' set of samples. When distinct is TRUE, the names function will return the
//' names that ONLY contain data from the given samples. When distinct is FALSE
//' the data returned contains data from a given samples, but may ALSO contain
//' data from other samples. Default = FALSE.
//'
//' @examples
//'
//' miseq <- miseq_sop_example()
//'
//' # To get the name of the dataset
//' xdev_names(data = miseq, type = "dataset")
//'
//' # To get the names of the sequences in the dataset
//' xdev_names(data = miseq, type = "sequences")
//'
//' # To get the names of the sequences that are unique to sample 'F3D0'
//' xdev_names(data = miseq, type = "sequences", samples = c("F3D0"), distinct = TRUE)
//'
//' # To get the names of the sequences that include sample 'F3D0'
//' xdev_names(data = miseq, type = "sequences", samples = c("F3D0"))
//'
//' # To get the names of the samples in the dataset
//' xdev_names(data = miseq, type = "samples")
//'
//' # To get the names of the treatments in the dataset
//' xdev_names(data = miseq, type = "treatments")
//'
//' # To get the names of the bins in the dataset
//' xdev_names(data = miseq, type = "bins")
//'
//' # To get the names of the bins in the dataset that are unique to 'F3D0'
//' xdev_names(data = miseq, type = "bins", samples = c("F3D0"), distinct = TRUE)
//'
//' # To get the names of the bins in the dataset that include sequences
//' # from 'F3D0'
//' xdev_names(data = miseq, type = "bins", samples = c("F3D0"), distinct = FALSE)
//'
//' # To get the names of the reports in the dataset
//' xdev_names(data = miseq, type = "reports")
//'
//' @return vector of strings, containing the names requested
//' @export
//[[Rcpp::export]]
vector<string> xdev_names(const Rcpp::Environment& data,
                           const string& type = "sequences",
                           const string& bin_type = "otu",
                           Rcpp::Nullable<Rcpp::List> samples = R_NilValue,
                           bool distinct = false);

// ************** removing ******************

//' @title xdev_remove_bins
//' @description
//' Designed with package integration in mind, the remove bins function allows
//' you to remove bins from a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param bin_names, a vector of strings containing the names of the bins you
//' would like removed.
//' @param trash_tags, a vector of strings containing the reasons you are
//' removing each bin
//' @param bin_type a string indicating the type of clusters.
//' @examples
//'
//'   data <- new_dataset(dataset_name = "my_dataset")
//'
//'   bin_names <- c("bin1", "bin2", "bin3")
//'   abundances <- c(110, 525, 80)
//'
//'   xdev_assign_bins(data = data, table = data.frame(bin_names = bin_names,
//'                                abundances = abundances), bin_type = "otu")
//'
//'   count(data = data, type = "bins", bin_type = "otu")
//'
//'   bins_to_remove <- c("bin1")
//'   trash_tag <- c("bad_bin")
//'
//'   xdev_remove_bins(data = data,
//'                    bin_names = bins_to_remove,
//'                    trash_tags = trash_tag)
//'
//'   count(data = data, type = "bins", bin_type = "otu")
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_remove_bins(const Rcpp::Environment& data, const vector<string>& bin_names,
                      const vector<string>& trash_tags, const string& bin_type = "otu");

//' @title xdev_remove_lineages
//' @description
//' Designed with package integration in mind, the remove lineages function
//' allows you to remove contaminents from a \href{https://mothur.org/strollur/reference/strollur.html}{strollur}
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param contaminants, vector of strings containing the taxonomies you would
//' like to remove
//' @param reason, a string containing reason you are removing the lineages.
//' Default = "contaminant".
//'
//' @examples
//' data <- read_mothur(fasta = strollur_example("final.fasta.gz"),
//'                        count = strollur_example("final.count_table.gz"),
//'                        taxonomy = strollur_example("final.taxonomy.gz"),
//'                        design = strollur_example("mouse.time.design"),
//'                        otu_list = strollur_example("final.opti_mcc.list.gz"),
//'                        dataset_name = "miseq_sop")
//'
//' contaminants <- c("Chloroplast", "Mitochondria", "unknown", "Archaea",
//'  "Eukaryota")
//'
//' xdev_remove_lineages(data = data, contaminants = contaminants)
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_remove_lineages(const Rcpp::Environment& data, const vector<string>& contaminants,
                          const string& reason = "contaminant");

//' @title xdev_remove_samples
//' @description
//' Designed with package integration in mind, the remove samples function allows
//' you to remove samples from a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param samples, vector of strings containing the names of the samples to
//' remove.
//'
//' @param reason, string containing the reason for removal.
//'  Default = "remove_samples".
//'
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' count(data = data, type = "samples")
//'
//' # To remove samples 'F3D0' and 'F3D1'
//'
//' xdev_remove_samples(data, c("F3D0", "F3D1"))
//'
//' count(data = data, type = "samples")
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_remove_samples(const Rcpp::Environment& data, const vector<string>& samples,
                         const string& reason = "remove_samples");

//' @title xdev_remove_sequences
//' @description
//' Designed with package integration in mind, the remove sequences function
//' allows you to remove sequences from a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param sequence_names, vector of strings containing the names of the
//' sequences to remove
//' @param trash_tags vector of strings containing the reasons for the sequences
//' removals
//'
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' count(data = data, type = "sequences")
//'
//' # For the sake of example let's remove the first 3 sequences from
//' # miseq_sop_example:
//'
//' seqs_to_remove <- c("M00967_43_000000000-A3JHG_1_2101_16474_12783",
//'                    "M00967_43_000000000-A3JHG_1_1113_12711_3318",
//'                    "M00967_43_000000000-A3JHG_1_2108_14707_9807")
//' trash_codes <- c("example", "removing", "sequences")
//'
//' xdev_remove_sequences(data = data, sequence_names = seqs_to_remove,
//'                       trash_tags = trash_codes)
//'
//' # If you look at the scrap report, you the sequences names, listed with the
//' # trash codes set to "example", "removing", "sequences".
//'
//' report(data = data, type = "sequence_scrap")
//'
//' # You can see from the get_num_sequences function that the removed
//' # sequence's abundances are removed from the dataset.
//'
//' count(data = data, type = "sequences")
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_remove_sequences(const Rcpp::Environment& data,
                           const vector<string>& sequence_names,
                           const vector<string>& trash_tags) ;

//' @title xdev_report
//' @description
//' Get a data.frame containing the given report in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param type, string containing the type of report you would like. Options
//' include: "fasta", "sequences", "sequence_bin_assignments",
//' "sequence_taxonomy", "bin_taxonomy", "bin_representatives",
//'  "sample_assignments", "metadata", "references", "sequence_scrap",
//' "bin_scrap". If you have added custom reports for alignment,
//' contigs_assembly or chimeras, you can get those as well.
//'  Default = "sequences".
//'
//' @param bin_type, string containing the bin type you would like a bin_taxonomy
//' report for. Default = "otu".
//'
//' @examples
//'
//' # First let's create a dataset from the \href{https://mothur.org/wiki/miseq_sop/}{MiSeq_SOP}
//'
//' miseq <- miseq_sop_example()
//'
//' # To get the FASTA data
//'
//' fasta <- xdev_report(data = miseq, type = "fasta")
//' head(fasta, n = 10)
//'
//' # To get a report about the FASTA data
//'
//' sequence_report <- xdev_report(data = miseq, type = "sequences")
//' head(sequence_report, n = 10)
//'
//' # To get the sequence bin assignments
//'
//' bin_assignments <- xdev_report(data = miseq, type = "sequence_bin_assignments",
//'                           bin_type = "otu")
//' head(bin_assignments, n = 10)
//'
//' # To get the sample treatment assignments
//'
//' xdev_report(data = miseq, type = "sample_assignments")
//'
//' # To get a report about sequence classifications
//'
//' sequence_taxonomy_report <- xdev_report(data = miseq,
//'                                        type = "sequence_taxonomy")
//' head(sequence_taxonomy_report, n = 10)
//'
//' # To get a report about bin classifications for 'otu' data
//'
//' otu_taxonomy_report <- xdev_report(data = miseq,
//'                                        type = "bin_taxonomy",
//'                                        bin_type = "otu")
//' head(otu_taxonomy_report, n = 10)
//'
//' # To get a report about bin classifications for 'asv' data
//'
//' asv_taxonomy_report <- xdev_report(data = miseq, type = "bin_taxonomy",
//'                               bin_type = "asv")
//' head(asv_taxonomy_report, n = 10)
//'
//' # To get a report about bin classifications for 'phylotype' data
//'
//' phylotype_taxonomy_report <- xdev_report(data = miseq, type = "bin_taxonomy",
//'                                     bin_type = "phylotype")
//' head(phylotype_taxonomy_report, n = 10)
//'
//' # To get the 'otu' bin representative sequences
//'
//' otu_bin_reps <- xdev_report(data = miseq, type = "bin_representatives",
//'                        bin_type = "otu")
//' head(otu_bin_reps, n = 10)
//'
//' # To get a report about the sequences removed during your analysis:
//'
//' scrapped_sequence_report <- xdev_report(data = miseq, type = "sequence_scrap")
//'
//' # To get a report about the "otu" bins removed during your analysis:
//'
//' scrapped_otu_report <- xdev_report(data = miseq, type = "bin_scrap",
//'                               bin_type = "otu")
//'
//' # To get a report about the "phylotype" bins removed during your analysis:
//'
//' scrapped_phylotype_report <- xdev_report(data = miseq, type = "bin_scrap",
//'                                     bin_type = "phylotype")
//'
//' # To get the metadata associated with your data:
//'
//' metadata <- xdev_report(data = miseq, type = "metadata")
//'
//' # To get the resource references associated with your data:
//'
//' references <- xdev_report(data = miseq, type = "references")
//'
//' # To get our custom report containing the contigs assembly data:
//'
//' contigs_report <- xdev_report(data = miseq, type = "contigs_report")
//' head(contigs_report, n = 10)
//'
//' @return data.frame
//' @export
//[[Rcpp::export]]
Rcpp::DataFrame xdev_report(const Rcpp::Environment& data, const string& type = "sequences",
                        const string& bin_type = "otu");
// ****************** setting *******************

//' @title xdev_set_abundance
//' @description
//' Designed with package integration in mind, the set abundance function
//' allows you to change the abundances of sequences in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' without samples.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param sequence_names, a vector of strings containing sequence names
//' @param sequence_abundances, vector containing the abundances of each
//' sequence.
//' @param reason, a string containing the trash tag to be applied to any
//' sequences set to 0 abundance. Default = "update".
//'
//' @examples
//'
//' names <- c("seq1", "seq2", "seq3",  "seq4")
//' abunds <- c(1250, 65, 50, 4)
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//'
//' xdev_assign_sequence_abundance(data = data, table = data.frame(sequence_names = names,
//'                                            abundances = abunds))
//' abundance(data = data, type = "sequences")
//'
//' seqs_to_update <- c("seq1", "seq3")
//' new_abunds <- c(1000, 100)
//'
//' xdev_set_abundance(data = data,
//'                    sequence_names = seqs_to_update,
//'                    sequence_abundances = new_abunds)
//'
//' abundance(data = data, type = "sequences")
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_set_abundance(const Rcpp::Environment& data,
                        const vector<string>& sequence_names,
                        const vector<float>& sequence_abundances,
                        const string& reason = "update");

//' @title xdev_set_abundances
//' @description
//' Designed with package integration in mind, the set abundances function
//' allows you to change the abundances of sequences in a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' with samples.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param sequence_names, a vector of strings containing sequence names
//' @param abundances, 2D vector ([num_seqs][num_samples]) containing
//' the abundances of each sequence parsed by sample.
//' @param reason, a string containing the trash tag to be applied to any
//' sequences set to 0 abundance. Default = "update".
//' @examples
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//'
//' sequence_names <- c("seq1", "seq1", "seq1", "seq2", "seq2", "seq2", "seq3",
//'                     "seq3", "seq4")
//' samples <- c("sample2", "sample3", "sample4", "sample2", "sample3",
//'              "sample4", "sample2", "sample3", "sample4")
//' abundances <- c(250, 400, 500, 25, 40, 50, 25, 25, 4)
//'
//' xdev_assign_sequence_abundance(data = data,
//'                           table = data.frame(sequence_names = sequence_names,
//'                                              abundances = abundances,
//'                                              samples = samples))
//'
//' seqs_to_update <- c("seq4")
//' new_abunds <- list(c(20, 10, 4))
//'
//' xdev_set_abundances(data = data,
//'                     sequence_names = seqs_to_update,
//'                     abundances = new_abunds)
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_set_abundances(const Rcpp::Environment& data,
                         const vector<string>& sequence_names,
                         const vector<vector<float>>& abundances,
                         const string& reason = "update");

//' @title xdev_set_sequences
//' @description
//' Designed with package integration in mind, the set sequences function allows
//' you to change the nucleotide strings of sequences in a
//' \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object. For example, set_sequences may be used
//' after alignment to overwrite the unaligned sequences with aligned sequences.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @param sequence_names, a vector of strings containing sequence names
//' @param sequences, a vector of strings containing sequence nucleotide strings
//' @param comments, a vector of strings containing sequence comments.
//' (Optional)
//'
//' @examples
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//'
//' xdev_add_sequences(data = data,
//'               table = data.frame(sequence_names = c("seq1", "seq2",
//'                                                   "seq3", "seq4")))
//'
//' xdev_set_sequences(data = data,
//'                    sequence_names = c("seq1", "seq2","seq3", "seq4"),
//'                    sequences = c("ATTGC", "ACTGC", "AGTGC", "TTTGC"))
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_set_sequences(const Rcpp::Environment& data,
                        const vector<string>& sequence_names,
                        const vector<string>& sequences,
                        const Rcpp::CharacterVector& comments = Rcpp::CharacterVector::create());

//' @title xdev_set_dataset_name
//' @description
//' Designed with package integration in mind, set the name of a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object.
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @param dataset_name, a string containing the desired name
//'
//' @examples
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//' xdev_set_dataset_name(data = data, dataset_name = "new_dataset_name")
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_set_dataset_name(const Rcpp::Environment& data, const string& dataset_name);

//' @title xdev_set_num_processors
//' @description
//' Designed with package integration in mind, set the number of processors used
//'  to summarize a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//' @param processors, a integer containing the desired number of processors
//' @examples
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//' xdev_set_num_processors(data = data, processors = 1)
//'
//' @return No return value, called for side effects.
//' @export
//[[Rcpp::export]]
void xdev_set_num_processors(const Rcpp::Environment& data, int processors);

/******************************************************************************/
//' @title xdev_summarize
//' @description
//' Summarize the sequences data, custom reports, and scrapped data in a
//' \link{strollur} object
//'
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @param type, string containing the type of data you want the number of.
//' Options include: "sequences", "reports" and "scrap". Default = "sequences".
//'
//' @param report_type, string containing the report type you would summarized.
//' For example, the miseq_sop_example includes contigs assembly data and can be
//' accessed with report_type = "contigs_report". Default = NULL.
//'
//' @examples
//'
//'  data <- miseq_sop_example()
//'
//'  # summarize FASTA data
//'  xdev_summarize(data = data, type = "sequences")
//'
//'  # summarize contigs_report
//'  xdev_summarize(data = data, type = "reports",
//'                  report_type = "contigs_report")
//'
//'  # remove sample 'F3D0'
//'  xdev_remove_samples(data = data, samples = c("F3D0"))
//'
//'  # summarize FASTA data after removal of sample F3D0
//'  xdev_summarize(data = data, type = "sequences")
//'
//'  # summarize scrapped data
//'  xdev_summarize(data = data, type = "scrap")
//'
//' @return data.frame()
//' @export
//[[Rcpp::export]]
Rcpp::DataFrame xdev_summarize(const Rcpp::Environment& data,
                               const string& type = "sequences",
                               Rcpp::Nullable<Rcpp::CharacterVector> report_type = R_NilValue);

// ***************** internal ******************

//' @title xint_copy_pointer
//' @name xint_copy_pointer
//' @description
//' For internal use only, copy an instance of the C++ 'Dataset' class.
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//'
//' @examples
//'
//' data <- read_mothur(fasta = strollur_example("final.fasta.gz"))
//'
//' copy_data <- new_dataset("copy")
//' copy_data
//'
//' copy_data$data <- xint_copy_pointer(data)
//' copy_data
//'
//' @return pointer to an instance of the C++ 'Dataset' class.
//' @keywords internal
//[[Rcpp::export]]
Rcpp::XPtr<Dataset> xint_copy_pointer(const Rcpp::Environment& data);

//' @title xint_new_pointer
//' @name xint_new_pointer
//' @description
//' For internal use only, create an instance of the C++ 'Dataset' class.
//' @param dataset_name, string containing dataset name
//' @param processors, number of processors to use
//' @return pointer to an instance of the C++ 'Dataset' class.
//' @keywords internal
//[[Rcpp::export]]
Rcpp::XPtr<Dataset> xint_new_pointer(const string& dataset_name, int processors);

//' @title xint_deserialize_dobject
//' @name xint_deserialize_dobject
//' @description
//' For internal use only, deserialize_dobject an instance of the C++ 'Dataset'
//'  class.
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//[[Rcpp::export]]
void xint_deserialize_dobject(Rcpp::Environment data);

//' @title xint_serialize_dobject
//' @name xint_serialize_dobject
//' @description
//' For internal use only, xint_serialize_dobject an instance of the C++ 'Dataset'
//' class.
//' @param data, a \href{https://mothur.org/strollur/reference/strollur.html}{strollur} object
//[[Rcpp::export]]
void xint_serialize_dobject(Rcpp::Environment data);

/******************************************************************************/

#endif
