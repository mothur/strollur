# Package index

## Getting Started

Adding and using data

- [`add()`](https://mothur.org/strollur/reference/add.md) : Add
  sequences, reports, metadata or resource references to a strollur
  object

- [`assign()`](https://mothur.org/strollur/reference/assign.md) :

  Assign sequence abundances, sequence classifications, bins, bin
  representative sequences, bin classifications or treatments to a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object

- [`names()`](https://mothur.org/strollur/reference/names.md) :

  Get the names sequences, bins, samples, treatments, and reports data
  in a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- [`count()`](https://mothur.org/strollur/reference/count.md) :

  Find the number of sequences, samples, treatments or bins of a given
  type in a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object

- [`abundance()`](https://mothur.org/strollur/reference/abundance.md) :

  Get the abundance data for sequences, bins, samples, and treatments in
  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- [`report()`](https://mothur.org/strollur/reference/report.md) :

  Get a data.frame containing the given report in a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object

- [`summary()`](https://mothur.org/strollur/reference/summary.md) :

  Summarize the sequences data, custom reports, and scrapped data in a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object

- [`strollur`](https://mothur.org/strollur/reference/strollur.md) : The
  \`strollur\` object stores the data associated with your microbial
  analysis.

- [`strollur-package`](https://mothur.org/strollur/reference/strollur-package.md)
  : strollur: Store and transfer data used in the analysis of microbial
  RNA

## Importing Data

Importing from other software tools namely: mothur, qiime2, phyloseq and
DADA2

- [`read_dada2()`](https://mothur.org/strollur/reference/read_dada2.md)
  :

  Create a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object
  from dada2 outputs

- [`read_fasta()`](https://mothur.org/strollur/reference/read_fasta.md)
  : read_fasta

- [`read_mothur()`](https://mothur.org/strollur/reference/read_mothur.md)
  :

  Create a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object
  from mothur outputs

- [`read_mothur_cons_taxonomy()`](https://mothur.org/strollur/reference/read_mothur_cons_taxonomy.md)
  : read_mothur_cons_taxonomy

- [`read_mothur_count()`](https://mothur.org/strollur/reference/read_mothur_count.md)
  : read_mothur_count

- [`read_mothur_list()`](https://mothur.org/strollur/reference/read_mothur_list.md)
  : read_mothur_list

- [`read_mothur_rabund()`](https://mothur.org/strollur/reference/read_mothur_rabund.md)
  : read_mothur_rabund

- [`read_mothur_shared()`](https://mothur.org/strollur/reference/read_mothur_shared.md)
  : read_mothur_shared

- [`read_mothur_taxonomy()`](https://mothur.org/strollur/reference/read_mothur_taxonomy.md)
  : read_mothur_taxonomy

- [`read_phyloseq()`](https://mothur.org/strollur/reference/read_phyloseq.md)
  :

  Create a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object
  from a phyloseq object

- [`read_qiime2()`](https://mothur.org/strollur/reference/read_qiime2.md)
  :

  Create a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object
  from a qiime2 outputs

- [`read_qiime2_feature_table()`](https://mothur.org/strollur/reference/read_qiime2_feature_table.md)
  : read_qiime2_feature_table

- [`read_qiime2_metadata()`](https://mothur.org/strollur/reference/read_qiime2_metadata.md)
  : read_qiime2_metadata

- [`read_qiime2_taxonomy()`](https://mothur.org/strollur/reference/read_qiime2_taxonomy.md)
  : read_qiime2_taxonomy

- [`unpack_qiime2_artifact()`](https://mothur.org/strollur/reference/unpack_qiime2_artifact.md)
  : unpack_qiime2_artifact

## General

Create datasets and references, clear data

- [`miseq_sop_example()`](https://mothur.org/strollur/reference/miseq_sop_example.md)
  : miseq_sop_example
- [`new_dataset()`](https://mothur.org/strollur/reference/new_dataset.md)
  : new_dataset
- [`new_reference()`](https://mothur.org/strollur/reference/new_reference.md)
  : new_reference
- [`clear()`](https://mothur.org/strollur/reference/clear.md) : clear
- [`is_aligned()`](https://mothur.org/strollur/reference/is_aligned.md)
  : is_aligned
- [`has_sample()`](https://mothur.org/strollur/reference/has_sample.md)
  : has_sample
- [`get_bin_types()`](https://mothur.org/strollur/reference/get_bin_types.md)
  : get_bin_types

## Data Transfers

Save, copy, load, import, export and write to file

- [`save_dataset()`](https://mothur.org/strollur/reference/save_dataset.md)
  : save_dataset
- [`load_dataset()`](https://mothur.org/strollur/reference/load_dataset.md)
  : load_dataset
- [`export_dataset()`](https://mothur.org/strollur/reference/export_dataset.md)
  : export_dataset
- [`import_dataset()`](https://mothur.org/strollur/reference/import_dataset.md)
  : import_dataset
- [`copy_dataset()`](https://mothur.org/strollur/reference/copy_dataset.md)
  : copy_dataset
- [`write_fasta()`](https://mothur.org/strollur/reference/write_fasta.md)
  : write_fasta
- [`write_mothur()`](https://mothur.org/strollur/reference/write_mothur.md)
  : write_mothur
- [`write_mothur_cons_taxonomy()`](https://mothur.org/strollur/reference/write_mothur_cons_taxonomy.md)
  : write_mothur_cons_taxonomy
- [`write_mothur_count()`](https://mothur.org/strollur/reference/write_mothur_count.md)
  : write_mothur_count
- [`write_mothur_design()`](https://mothur.org/strollur/reference/write_mothur_design.md)
  : write_mothur_design
- [`write_mothur_list()`](https://mothur.org/strollur/reference/write_mothur_list.md)
  : write_mothur_list
- [`write_mothur_rabund()`](https://mothur.org/strollur/reference/write_mothur_rabund.md)
  : write_mothur_rabund
- [`write_mothur_shared()`](https://mothur.org/strollur/reference/write_mothur_shared.md)
  : write_mothur_shared
- [`write_phyloseq()`](https://mothur.org/strollur/reference/write_phyloseq.md)
  : write_phyloseq
- [`write_taxonomy()`](https://mothur.org/strollur/reference/write_taxonomy.md)
  : write_taxonomy

## Functions for Package Developers

Want to create and modify strollur objects from your package? Check out
these functions

- [`strollur_example()`](https://mothur.org/strollur/reference/strollur_example.md)
  : strollur_example

- [`get_available_processors()`](https://mothur.org/strollur/reference/get_available_processors.md)
  : get_available_processors

- [`has_sequence_strings()`](https://mothur.org/strollur/reference/has_sequence_strings.md)
  : has_sequence_strings

- [`added_message()`](https://mothur.org/strollur/reference/added_message.md)
  : added_message

- [`assigned_message()`](https://mothur.org/strollur/reference/assigned_message.md)
  : assigned_message

- [`remove_file()`](https://mothur.org/strollur/reference/remove_file.md)
  : remove_file

- [`sort_dataframe()`](https://mothur.org/strollur/reference/sort_dataframe.md)
  : sort_dataframe

- [`xdev_abundance()`](https://mothur.org/strollur/reference/xdev_abundance.md)
  : Get a data.frame containing the requested abundance data

- [`xdev_add_references()`](https://mothur.org/strollur/reference/xdev_add_references.md)
  : Add resource references

- [`xdev_add_report()`](https://mothur.org/strollur/reference/xdev_add_report.md)
  :

  Add a report to a
  [strollur](https://mothur.org/strollur/reference/strollur.html) object

- [`xdev_add_sequences()`](https://mothur.org/strollur/reference/xdev_add_sequences.md)
  : xdev_add_sequences

- [`xdev_assign_bin_representative_sequences()`](https://mothur.org/strollur/reference/xdev_assign_bin_representative_sequences.md)
  : xdev_assign_bin_representative_sequences

- [`xdev_assign_bin_taxonomy()`](https://mothur.org/strollur/reference/xdev_assign_bin_taxonomy.md)
  : xdev_assign_bin_taxonomy

- [`xdev_assign_bins()`](https://mothur.org/strollur/reference/xdev_assign_bins.md)
  : xdev_assign_bins

- [`xdev_assign_sequence_abundance()`](https://mothur.org/strollur/reference/xdev_assign_sequence_abundance.md)
  : xdev_assign_sequence_abundance

- [`xdev_assign_sequence_taxonomy()`](https://mothur.org/strollur/reference/xdev_assign_sequence_taxonomy.md)
  : xdev_assign_sequence_taxonomy

- [`xdev_assign_treatments()`](https://mothur.org/strollur/reference/xdev_assign_treatments.md)
  : xdev_assign_treatments

- [`xdev_count()`](https://mothur.org/strollur/reference/xdev_count.md)
  : xdev_count

- [`xdev_get_abundances_by_sample()`](https://mothur.org/strollur/reference/xdev_get_abundances_by_sample.md)
  : xdev_get_abundances_by_sample

- [`xdev_get_by_sample()`](https://mothur.org/strollur/reference/xdev_get_by_sample.md)
  : xdev_get_by_sample

- [`xdev_get_list_vector()`](https://mothur.org/strollur/reference/xdev_get_list_vector.md)
  : xdev_get_list_vector

- [`xdev_get_sequences()`](https://mothur.org/strollur/reference/xdev_get_sequences.md)
  : xdev_get_sequences

- [`xdev_has_sequence_taxonomy()`](https://mothur.org/strollur/reference/xdev_has_sequence_taxonomy.md)
  : xdev_has_sequence_taxonomy

- [`xdev_merge_bins()`](https://mothur.org/strollur/reference/xdev_merge_bins.md)
  : xdev_merge_bins

- [`xdev_merge_sequences()`](https://mothur.org/strollur/reference/xdev_merge_sequences.md)
  : xdev_merge_sequences

- [`xdev_names()`](https://mothur.org/strollur/reference/xdev_names.md)
  : xdev_names

- [`xdev_remove_bins()`](https://mothur.org/strollur/reference/xdev_remove_bins.md)
  : xdev_remove_bins

- [`xdev_remove_lineages()`](https://mothur.org/strollur/reference/xdev_remove_lineages.md)
  : xdev_remove_lineages

- [`xdev_remove_samples()`](https://mothur.org/strollur/reference/xdev_remove_samples.md)
  : xdev_remove_samples

- [`xdev_remove_sequences()`](https://mothur.org/strollur/reference/xdev_remove_sequences.md)
  : xdev_remove_sequences

- [`xdev_report()`](https://mothur.org/strollur/reference/xdev_report.md)
  : xdev_report

- [`xdev_set_abundance()`](https://mothur.org/strollur/reference/xdev_set_abundance.md)
  : xdev_set_abundance

- [`xdev_set_abundances()`](https://mothur.org/strollur/reference/xdev_set_abundances.md)
  : xdev_set_abundances

- [`xdev_set_dataset_name()`](https://mothur.org/strollur/reference/xdev_set_dataset_name.md)
  : xdev_set_dataset_name

- [`xdev_set_num_processors()`](https://mothur.org/strollur/reference/xdev_set_num_processors.md)
  : xdev_set_num_processors

- [`xdev_set_sequences()`](https://mothur.org/strollur/reference/xdev_set_sequences.md)
  : xdev_set_sequences

- [`xdev_summarize()`](https://mothur.org/strollur/reference/xdev_summarize.md)
  : xdev_summarize

- [`xint_deserialize_dobject()`](https://mothur.org/strollur/reference/xint_deserialize_dobject.md)
  : xint_deserialize_dobject

- [`xint_serialize_dobject()`](https://mothur.org/strollur/reference/xint_serialize_dobject.md)
  : xint_serialize_dobject
