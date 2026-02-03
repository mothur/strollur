# Package index

## Getting Started

Adding and accessing data

- [`add()`](add.md) : add
- [`assign()`](assign.md) : assign
- [`names()`](names.md) : names
- [`count()`](count.md) : count
- [`abundance()`](abundance.md) : abundance
- [`report()`](report.md) : report
- [`summary()`](summary.md) : summary
- [`dataset`](dataset.md) : dataset
- [`rdataset-package`](rdataset-package.md)
  [`rdataset`](rdataset-package.md) : Store and transfer data used in
  the analysis of microbial DNA
- [`get_bin_types()`](get_bin_types.md) : get_bin_types

## Importing Data

Importing from other software tools namely: mothur, qiime2, phyloseq and
DADA2

- [`read_dada2()`](read_dada2.md) : read_dada2
- [`read_fasta()`](read_fasta.md) : read_fasta
- [`read_mothur()`](read_mothur.md) : read_mothur
- [`read_mothur_cons_taxonomy()`](read_mothur_cons_taxonomy.md) :
  read_mothur_cons_taxonomy
- [`read_mothur_count()`](read_mothur_count.md) : read_mothur_count
- [`read_mothur_list()`](read_mothur_list.md) : read_mothur_list
- [`read_mothur_rabund()`](read_mothur_rabund.md) : read_mothur_rabund
- [`read_mothur_shared()`](read_mothur_shared.md) : read_mothur_shared
- [`read_mothur_taxonomy()`](read_mothur_taxonomy.md) :
  read_mothur_taxonomy
- [`read_qiime2()`](read_qiime2.md) : read_qiime2
- [`read_qiime2_feature_table()`](read_qiime2_feature_table.md) :
  read_qiime2_feature_table
- [`read_qiime2_metadata()`](read_qiime2_metadata.md) :
  read_qiime2_metadata
- [`read_qiime2_taxonomy()`](read_qiime2_taxonomy.md) :
  read_qiime2_taxonomy
- [`unpack_qiime2_artifact()`](unpack_qiime2_artifact.md) :
  unpack_qiime2_artifact

## General

Create datasets and references, clear data

- [`miseq_sop_example()`](miseq_sop_example.md) : miseq_sop_example
- [`new_dataset()`](new_dataset.md) : new_dataset
- [`new_reference()`](new_reference.md) : new_reference
- [`clear()`](clear.md) : clear
- [`is_aligned()`](is_aligned.md) : is_aligned
- [`has_sample()`](has_sample.md) : has_sample

## Data Transfers

Save, copy, load, import, export and write to file

- [`save_dataset()`](save_dataset.md) : save_dataset
- [`load_dataset()`](load_dataset.md) : load_dataset
- [`export_dataset()`](export_dataset.md) : export_dataset
- [`import_dataset()`](import_dataset.md) : import_dataset
- [`copy_dataset()`](copy_dataset.md) : copy_dataset
- [`write_fasta()`](write_fasta.md) : write_fasta
- [`write_mothur()`](write_mothur.md) : write_mothur
- [`write_mothur_cons_taxonomy()`](write_mothur_cons_taxonomy.md) :
  write_mothur_cons_taxonomy
- [`write_mothur_count()`](write_mothur_count.md) : write_mothur_count
- [`write_mothur_design()`](write_mothur_design.md) :
  write_mothur_design
- [`write_mothur_list()`](write_mothur_list.md) : write_mothur_list
- [`write_mothur_rabund()`](write_mothur_rabund.md) :
  write_mothur_rabund
- [`write_mothur_shared()`](write_mothur_shared.md) :
  write_mothur_shared
- [`write_taxonomy()`](write_taxonomy.md) : write_taxonomy

## Functions for Package Developers

Want to create and modify dataset objects from your package? Check out
these functions

- [`rdataset_example()`](rdataset_example.md) : rdataset_example
- [`get_available_processors()`](get_available_processors.md) :
  get_available_processors
- [`has_sequence_strings()`](has_sequence_strings.md) :
  has_sequence_strings
- [`added_message()`](added_message.md) : added_message
- [`assigned_message()`](assigned_message.md) : assigned_message
- [`remove_file()`](remove_file.md) : remove_file
- [`sort_dataframe()`](sort_dataframe.md) : sort_dataframe
- [`xdev_abundance()`](xdev_abundance.md) : xdev_abundance
- [`xdev_add_references()`](xdev_add_references.md) :
  xdev_add_references
- [`xdev_add_report()`](xdev_add_report.md) : xdev_add_report
- [`xdev_add_sequences()`](xdev_add_sequences.md) : xdev_add_sequences
- [`xdev_assign_bin_representative_sequences()`](xdev_assign_bin_representative_sequences.md)
  : xdev_assign_bin_representative_sequences
- [`xdev_assign_bin_taxonomy()`](xdev_assign_bin_taxonomy.md) :
  xdev_assign_bin_taxonomy
- [`xdev_assign_bins()`](xdev_assign_bins.md) : xdev_assign_bins
- [`xdev_assign_sequence_abundance()`](xdev_assign_sequence_abundance.md)
  : xdev_assign_sequence_abundance
- [`xdev_assign_sequence_taxonomy()`](xdev_assign_sequence_taxonomy.md)
  : xdev_assign_sequence_taxonomy
- [`xdev_assign_treatments()`](xdev_assign_treatments.md) :
  xdev_assign_treatments
- [`xdev_count()`](xdev_count.md) : xdev_count
- [`xdev_get_by_sample()`](xdev_get_by_sample.md) : xdev_get_by_sample
- [`xdev_get_list_vector()`](xdev_get_list_vector.md) :
  xdev_get_list_vector
- [`xdev_get_sequences()`](xdev_get_sequences.md) : xdev_get_sequences
- [`xdev_merge_bins()`](xdev_merge_bins.md) : xdev_merge_bins
- [`xdev_merge_sequences()`](xdev_merge_sequences.md) :
  xdev_merge_sequences
- [`xdev_names()`](xdev_names.md) : xdev_names
- [`xdev_remove_bins()`](xdev_remove_bins.md) : xdev_remove_bins
- [`xdev_remove_lineages()`](xdev_remove_lineages.md) :
  xdev_remove_lineages
- [`xdev_remove_samples()`](xdev_remove_samples.md) :
  xdev_remove_samples
- [`xdev_remove_sequences()`](xdev_remove_sequences.md) :
  xdev_remove_sequences
- [`xdev_report()`](xdev_report.md) : xdev_report
- [`xdev_set_abundance()`](xdev_set_abundance.md) : xdev_set_abundance
- [`xdev_set_abundances()`](xdev_set_abundances.md) :
  xdev_set_abundances
- [`xdev_set_dataset_name()`](xdev_set_dataset_name.md) :
  xdev_set_dataset_name
- [`xdev_set_num_processors()`](xdev_set_num_processors.md) :
  xdev_set_num_processors
- [`xdev_set_sequences()`](xdev_set_sequences.md) : xdev_set_sequences
- [`xdev_summarize()`](xdev_summarize.md) : xdev_summarize
- [`xint_deserialize_dobject()`](xint_deserialize_dobject.md) :
  xint_deserialize_dobject
- [`xint_serialize_dobject()`](xint_serialize_dobject.md) :
  xint_serialize_dobject
