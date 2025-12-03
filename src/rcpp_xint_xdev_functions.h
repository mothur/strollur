#ifndef INTERNAL_DEVELOPMENT_H_
#define INTERNAL_DEVELOPMENT_H_

#include <Rcpp.h>
#include "../inst/include/rdataset.h"
#include "dataset.h"

/******************************************************************************/
SEXP xint_fill_required_parameters(const Rcpp::DataFrame df,
                                   const string& given_column_name,
                                   string type = "string");

SEXP xint_fill_optional_parameters(const Rcpp::DataFrame df,
                                   const string& default_column_name,
                                   const string& given_column_name,
                                   string type = "string");


/******************************************************************************/
//' @title xdev_get_by_sample
//' @description
//' Get the requested data in a \link{dataset} object parsed by sample
//'
//' @param data, a \link{dataset} object
//'
//' @param type, string containing the type of data you want the totals of.
//' Options include: "sequence_names", "sequences". Default = "sequence_names".
//'
//' @param samples a vector of strings containing the names of the samples you
//' would like sequence names for. By default all samples are included.
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' # To get the sequence names parsed by sample
//' xdev_get_by_sample(data, "sequence_names")
//'
//' # To get the sequence nucleotide strings parsed by sample
//' xdev_get_by_sample(data, "sequences")
//'
//' @return 2D vector of strings ([num_seqs][num_samples]) containing data
//' requested parsed by sample.
//[[Rcpp::export]]
vector<vector<string> > xdev_get_by_sample(Rcpp::Environment data,
                                      string type = "sequence_names",
                                      Rcpp::CharacterVector samples = Rcpp::CharacterVector::create());

// ******** merging *********

//' @title xdev_merge_bins
//' @description
//' Designed with package integration in mind, the merge bins function allows
//' you to merge bins in a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//'
//' @param bin_names, a vector of strings containing the names of the bins you
//' would like merge. The resulting merged bin will be stored in the first
//' bin_id in the vector.
//' @param reason, a string indicating why you are merging bins. Default =
//' "merged".
//' @param type, a string indicating the type of bin clusters. Default = "otu"
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
//[[Rcpp::export]]
void xdev_merge_bins(Rcpp::Environment data, vector<string> bin_names,
                     string reason = "merged", string type = "otu");

//' @title xdev_merge_sequences
//' @description
//' Designed with package integration in mind, the merge sequences function
//' allows you to merge sequences in a \link{dataset} object.
//'
//' @param data, a \link{dataset} object.
//'
//' @param sequence_names, a vector of strings containing the names of the
//' sequences you would like merge. The resulting merged sequence will be stored
//' in the first sequence name in the vector.
//' @param reason a string indicating why you are merging sequences.
//' Default = "merged"
//'
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' num(data = data, type = "sequences")
//'
//' # For the sake of example let's merge the first 3 sequences from
//' # miseq_sop_example:
//'
//' seqs_to_merge <- c("M00967_43_000000000-A3JHG_1_2101_16474_12783",
//'                    "M00967_43_000000000-A3JHG_1_1113_12711_3318",
//'                    "M00967_43_000000000-A3JHG_1_2108_14707_9807")
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
//' num(data = data, type = "sequences")
//'
//[[Rcpp::export]]
void xdev_merge_sequences(Rcpp::Environment data, vector<string> sequence_names,
                          string reason = "merged");

/******************************************************************************/
//' @title xdev_names
//' @description
//' Get the names of a given type of data in a \link{dataset} object
//'
//' @param data, a \link{dataset} object
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
//[[Rcpp::export]]
const vector<string> xdev_names(Rcpp::Environment data,
                           string type = "sequences",
                           string bin_type = "otu",
                           Rcpp::Nullable<Rcpp::List> samples = R_NilValue,
                           bool distinct = false);

// ************** removing ******************

//' @title xdev_remove_bins
//' @description
//' Designed with package integration in mind, the remove bins function allows
//' you to remove bins from a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//'
//' @param bin_names, a vector of strings containing the names of the bins you
//' would like removed.
//' @param trash_tags, a vector of strings containing the reasons you are
//' removing each bin
//' @param type a string indicating the type of clusters.
//' @examples
//'
//'   data <- new_dataset(dataset_name = "my_dataset")
//'
//'   bin_names <- c("bin1", "bin2", "bin3")
//'   abundances <- c(110, 525, 80)
//'
//'   assign_bins(data = data, table = data.frame(bin_names = bin_names,
//'                                abundances = abundances), bin_type = "otu")
//'
//'   num(data = data, type = "bins", bin_type = "otu")
//'
//'   bins_to_remove <- c("bin1")
//'   trash_tag <- c("bad_bin")
//'
//'   xdev_remove_bins(data = data,
//'                    bin_names = bins_to_remove,
//'                    trash_tags = trash_tag)
//'
//'   num(data = data, type = "bins", bin_type = "otu")
//'
//[[Rcpp::export]]
void xdev_remove_bins(Rcpp::Environment data, vector<string> bin_names,
                      vector<string> trash_tags, string type = "otu");

//' @title xdev_remove_lineages
//' @description
//' Designed with package integration in mind, the remove lineages function
//' allows you to remove contaminents from a \link{dataset}
//'
//' @param data, a \link{dataset} object.
//'
//' @param contaminants, vector of strings containing the taxonomies you would
//' like to remove
//' @param trash_tag, a string containing reason you are removing the lineages.
//' Default = "contaminant".
//'
//' @examples
//' data <- read_mothur(fasta = rdataset_example("final.fasta"),
//'                        count = rdataset_example("final.count_table"),
//'                        taxonomy = rdataset_example("final.taxonomy"),
//'                        design = rdataset_example("mouse.time.design"),
//'                        otu_list = rdataset_example("final.opti_mcc.list"),
//'                        dataset_name = "miseq_sop")
//'
//' contaminants <- c("Chloroplast", "Mitochondria", "unknown", "Archaea",
//'  "Eukaryota")
//'
//' xdev_remove_lineages(data = data, contaminants = contaminants)
//'
//[[Rcpp::export]]
void xdev_remove_lineages(Rcpp::Environment data, vector<string> contaminants,
                          string trash_tag = "contaminant");

//' @title xdev_remove_samples
//' @description
//' Designed with package integration in mind, the remove samples function allows
//' you to remove samples from a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
//'
//' @param samples, vector of strings containing the names of the samples to
//' remove.
//'
//' @examples
//'
//' data <- miseq_sop_example()
//'
//' num(data = data, type = "samples")
//'
//' # To remove samples 'F3D0' and 'F3D1'
//'
//' xdev_remove_samples(data, c("F3D0", "F3D1"))
//'
//' num(data = data, type = "samples")
//'
//[[Rcpp::export]]
void xdev_remove_samples(Rcpp::Environment data, vector<string> samples);

//' @title xdev_remove_sequences
//' @description
//' Designed with package integration in mind, the remove sequences function
//' allows you to remove sequences from a \link{dataset} object
//'
//' @param data, a \link{dataset} object.
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
//' num(data = data, type = "sequences")
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
//' num(data = data, type = "sequences")
//'
//[[Rcpp::export]]
void xdev_remove_sequences(Rcpp::Environment data,
                           vector<string> sequence_names,
                           vector<string> trash_tags) ;

// ****************** setting *******************

//' @title xdev_set_abundance
//' @description
//' Designed with package integration in mind, the set abundance function
//' allows you to change the abundances of sequences in a \link{dataset} object
//' without samples.
//'
//' @param data, a \link{dataset} object
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
//' assign_sequence_abundance(data = data, table = data.frame(sequence_names = names,
//'                                            abundances = abunds))
//' get_sequence_abundances(data = data)
//'
//' seqs_to_update <- c("seq1", "seq3")
//' new_abunds <- c(1000, 100)
//'
//' xdev_set_abundance(data = data,
//'                    sequence_names = seqs_to_update,
//'                    sequence_abundances = new_abunds)
//'
//' get_sequence_abundances(data = data)
//'
//[[Rcpp::export]]
void xdev_set_abundance(Rcpp::Environment data,
                        vector<string> sequence_names,
                        vector<float> sequence_abundances,
                        string reason = "update");

//' @title xdev_set_abundances
//' @description
//' Designed with package integration in mind, the set abundances function
//' allows you to change the abundances of sequences in a \link{dataset} object
//' with samples.
//'
//' @param data, a \link{dataset} object
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
//' assign_sequence_abundance(data = data,
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
//[[Rcpp::export]]
void xdev_set_abundances(Rcpp::Environment data,
                         vector<string> sequence_names,
                         vector<vector<float>> abundances,
                         string reason = "update");

//' @title xdev_set_bin_abundance
//' @description
//' Designed with package integration in mind, the set bin abundance function
//' allows you to change the abundances of bins in a \link{dataset} object
//' without sample data.
//'
//' @param data, a \link{dataset} object
//'
//' @param bin_names, a vector strings containing of bin names to set the
//' abundances for.
//' @param abundances, vector containing the abundances of each bin.
//' @param type, a string indicating the type of clusters. Default = "otu".
//' @param reason, a string containing the trash tag to be applied to any bins
//'  set to 0 abundance. Default = "update".
//'
//' @examples
//'   # For example sake, let's create a dataset with 3 bins:
//'
//'   data <- new_dataset(dataset_name = "my_dataset")
//'
//'   bin_ids <- c("bin1", "bin2", "bin3")
//'   abundances <- c(110, 525, 80)
//'
//'   assign_bins(data = data, table = data.frame(bin_names = bin_ids,
//'                                               abundances = abundances))
//'
//'   get_bin_abundance(data, "bin1")
//'   get_bin_abundance(data, "bin2")
//'
//'   # Now we can use set_bin_abundance to change the abundances of bin1 and
//'   # bin2
//'
//'   bins <- c("bin1", "bin2")
//'   new_abunds <- c(300, 250)
//'
//'   xdev_set_bin_abundance(data = data,
//'                          bin_names = bins,
//'                          abundances = new_abunds)
//'
//'   get_bin_abundance(data, "bin1")
//'   get_bin_abundance(data, "bin2")
//'
//[[Rcpp::export]]
void xdev_set_bin_abundance(Rcpp::Environment data,
                            vector<string> bin_names,
                            vector<float> abundances,
                            string type = "otu",
                            string reason = "update");

//' @title xdev_set_bin_abundances
//' @description
//' Designed with package integration in mind, the set bin abundances function
//' allows you to change the abundances of bins in a \link{dataset} object
//' with sample data.
//'
//' @param data, a \link{dataset} object
//'
//' @param bin_names, a vector strings containing of bin names to set the
//' abundances for.
//' @param abundances, 2D vector ([num_seqs][num_samples]) containing the
//' abundances of each bin parsed by sample.
//' @param type a string indicating the type of clusters. Default = "otu".
//' @param reason, a string containing the trash tag to be applied to any bins
//'  set to 0 abundance. Default = "update".
//'
//' @examples
//'
//'   # For example sake, let's create a dataset with 3 bins:
//'
//'   data <- new_dataset(dataset_name = "my_dataset")
//'
//'   bin_ids <- c("bin1", "bin1", "bin1", "bin2", "bin2", "bin3")
//'   samples <- c("sample1", "sample2", "sample5", "sample1", "sample3",
//'                "sample1")
//'   sample_abundances <- c(10, 100, 1, 500, 25, 80)
//'
//'   assign_bins(data = data, table = data.frame(bin_names = bin_ids,
//'                                               abundances = sample_abundances,
//'                                               samples = samples))
//'
//'   # You can see bin1's abundances parsed by sample using get_bin_abundances:
//'   get_bin_abundances(data = data, bin_name = "bin1")
//'
//'   # You can change bin1's abundances as follows:
//'
//'   new_bin1_abunds <- list(c(10,50,0,0))
//'   bins <- c("bin1")
//'
//'   xdev_set_bin_abundances(data = data,
//'                           bin_names = bins,
//'                           abundances = new_bin1_abunds)
//'
//'   get_bin_abundances(data = data, bin_name = "bin1")
//'
//[[Rcpp::export]]
void xdev_set_bin_abundances(Rcpp::Environment data,
                             vector<string> bin_names,
                             vector<vector<float>> abundances,
                             string type = "otu", string reason = "update");

//' @title xdev_set_sequences
//' @description
//' Designed with package integration in mind, the set sequences function
//' allows you to change the nucleotide strings of sequences in a \link{dataset}
//' object. For example, set_sequences may be used after alignment to overwrite
//' the unaligned sequences with aligned sequences.
//'
//' @param data, a \link{dataset} object
//' @param sequence_names, a vector of strings containing sequence names
//' @param sequences, a vector of strings containing sequence nucleotide strings
//' @param comments, a vector of strings containing sequence comments.
//' (Optional)
//'
//' @examples
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//'
//' add_sequences(data = data,
//'               table = data.frame(sequence_names = c("seq1", "seq2",
//'                                                   "seq3", "seq4")))
//'
//' xdev_set_sequences(data = data,
//'                    sequence_names = c("seq1", "seq2","seq3", "seq4"),
//'                    sequences = c("ATTGC", "ACTGC", "AGTGC", "TTTGC"))
//'
//[[Rcpp::export]]
void xdev_set_sequences(Rcpp::Environment data,
                        vector<string> sequence_names,
                        vector<string> sequences,
                        Rcpp::CharacterVector comments = Rcpp::CharacterVector::create());

//' @title xdev_set_dataset_name
//' @description
//' Designed with package integration in mind, set the name of a \link{dataset} object.
//'
//' @param data, a \link{dataset} object
//' @param dataset_name, a string containing the desired name
//'
//' @examples
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//' xdev_set_dataset_name(data = data, dataset_name = "new_dataset_name")
//'
//[[Rcpp::export]]
void xdev_set_dataset_name(Rcpp::Environment data, string dataset_name);

//' @title xdev_set_num_processors
//' @description
//' Designed with package integration in mind, set the number of processors used
//'  to summarize a \link{dataset} object
//'
//' @param data, a \link{dataset} object
//' @param processors, a integer containing the desired number of processors
//' @examples
//'
//' data <- new_dataset(dataset_name = "my_dataset")
//' xdev_set_num_processors(data = data, processors = 1)
//'
//[[Rcpp::export]]
void xdev_set_num_processors(Rcpp::Environment data, int processors);

// ***************** internal ******************

//' @title xint_copy_pointer
//' @name xint_copy_pointer
//' @description
//' For internal use only, copy an instance of the C++ 'Dataset' class.
//' @param data, a \link{dataset} object
//' @return pointer to an instance of the C++ 'Dataset' class.
//' @keywords internal
//[[Rcpp::export]]
Rcpp::XPtr<Dataset> xint_copy_pointer(Rcpp::Environment data);

//' @title xint_new_pointer
//' @name xint_new_pointer
//' @description
//' For internal use only, create an instance of the C++ 'Dataset' class.
//' @param dataset_name, string containing dataset name
//' @param processors, number of processors to use
//' @return pointer to an instance of the C++ 'Dataset' class.
//' @keywords internal
//[[Rcpp::export]]
Rcpp::XPtr<Dataset> xint_new_pointer(string dataset_name, int processors);

//' @title xint_deserialize_dobject
//' @name xint_deserialize_dobject
//' @description
//' For internal use only, deserialize_dobject an instance of the C++ 'Dataset'
//'  class.
//' @param data, a \link{dataset} object
//[[Rcpp::export]]
void xint_deserialize_dobject(Rcpp::Environment data);

//' @title xint_serialize_dobject
//' @name xint_serialize_dobject
//' @description
//' For internal use only, xint_serialize_dobject an instance of the C++ 'Dataset'
//' class.
//' @param data, a \link{dataset} object
//[[Rcpp::export]]
void xint_serialize_dobject(Rcpp::Environment data);

/******************************************************************************/

#endif
