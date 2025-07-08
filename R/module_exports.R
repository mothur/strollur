#' @name Dataset
#' @title Dataset
#' @description The 'Dataset' class is the backend C++ implementation of the R6
#' 'sequence_data' object. Along with the rcpp_module this class allows
#' package developers access to additional functionality. 'Dataset' stores
#' nucleotide sequences, abundance, sample and treatment assignments,
#' taxonomic classifications, asv / otu clusters. It creates various reports
#' and summaries. It is designed to facilitate data transfer and access across
#' multiple R packages.
#'
#' @field new Constructs a new Dataset object \itemize{
#' \item \emph{name}, a string containing the dataset name
#' \item \emph{processors}, an integer indicating the number of processors
#' \item Returns: 'Dataset' object
#' \item \code{data <- new(Dataset, "my_dataset", 4)}
#' }
#'
#' @field add_align_report Add alignment data
#' \itemize{
#' \item \emph{names}, a vector of strings containing the sequence names
#' \item \emph{search_scores}, a vector of doubles containing kmer search score
#' used to find the closest reference sequence match
#' \item \emph{sim_scores}, a vector of doubles containing the similarity score
#' between the reference sequence and the query sequence
#' \item \emph{longest_inserts}, a vector of integers containing the longest
#' insert added to the query sequence
#' \item \code{data$add_align_report(names, search_scores, sim_scores,
#'  longest_inserts)}
#' }
#'
#' @field add_contigs_report Add paired read assembly data
#' \itemize{
#' \item \emph{names}, a vector of strings containing the sequence names
#' \item \emph{overlap_lengths}, a vector of integers containing length of the
#' overlap between the paired reads
#' \item \emph{overlap_starts}, a vector of integers containing the start
#' position of the paired reads overlap
#' \item \emph{overlap_ends}, a vector of integers containing the end
#' position of the paired reads overlap
#' \item \emph{mismatches}, a vector of integers containing the number of
#' mismatches between the paired reads
#' \item \emph{expected_errors}, a vector of doubles containing the number of
#' expected errors based on the quality scores of the paired reads
#' \item \code{data$add_contigs_report(names, overlap_lengths, overlap_starts,
#'  overlap_ends, mismatches, expected_errors)}
#' }
#'
#' @field add_sequences Add sequences to dataset
#' \itemize{
#' \item \emph{names}, a vector of strings containing the sequence names
#' \item \emph{sequences}, a vector of strings containing the sequence
#' nucleotides
#' \item \emph{comments}, a vector of strings containing the sequence comments
#' \item \code{data$add_sequences(names, sequences, comments)}
#' }
#'
#'
#' @field dataset_name a string containing the dataset name
#' @field is_aligned boolean indicating the dataset's alignment status
#' @field num_samples integer containing the number of samples in the dataset
#' @field num_treatments integer containing the number of treatments in the
#' dataset
#' @field num_otus integer containing the number of otus in the dataset
#' @field num_unique integer containing the number of distinct sequences in the
#'  dataset
#' @field has_align_data boolean indicating the dataset includes an alignment
#'  report
#' @field has_contigs_data boolean indicating the dataset includes a contigs
#'  assembly report
#'

loadModule(module = "Dataset", TRUE)
