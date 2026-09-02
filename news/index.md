# Changelog

## strollur 0.1.4

New Features:

- Expands built in read / write functions.
  - [`strollur::read_fastq()`](https://mothur.org/strollur/reference/read_fastq.md)
    returns data.frame with fastq data.
  - [`strollur::write_fastq()`](https://mothur.org/strollur/reference/write_fastq.md)
    writes fastq data to file.
  - [`strollur::read_quality()`](https://mothur.org/strollur/reference/read_quality.md)
    returns data.frame with sequence quality data.
  - [`strollur::write_quality()`](https://mothur.org/strollur/reference/write_quality.md)
    writes sequence quality data to file.
  - [`strollur::read_biom()`](https://mothur.org/strollur/reference/read_biom.md)
    returns a
    [strollur_object](https://mothur.org/strollur/reference/strollur.html).
  - [`strollur::write_biom()`](https://mothur.org/strollur/reference/write_biom.md)
    writes biom files containing the
    [strollur_object](https://mothur.org/strollur/reference/strollur.html).
- Expands member functions of the
  [strollur_object](https://mothur.org/strollur/reference/strollur.html).
  - `strollur_object$alignment_length()` returns the sequences alignment
    length or -1 if the sequences are unaligned.
  - `strollur_object$has_bin_taxonomy()` returns `TRUE` if your dataset
    includes bin classifications.
  - `strollur_object$has_sequence_taxonomy()` returns `TRUE` if your
    dataset includes sequence classifications.
  - `strollur_object$rename()` rename the strollur object.

Enhancements:

- Adds `fastq`, `quality`, `sequence_tree`, `sample_tree`,
  `sample_distance` types to:
  - [`strollur::report()`](https://mothur.org/strollur/reference/report.md)
  - `strollur_object$report()`
- Adds `fastq` type to:
  - [`strollur::add()`](https://mothur.org/strollur/reference/add.md)
  - `strollur_object$add()`
- Adds `quality` and `sample_distance` type to:
  - [`strollur::assign()`](https://mothur.org/strollur/reference/assign.md)
  - `strollur_object$assign()`

Functions for Package developers

- Adds
  [`strollur::xdev_add_sequence_fastq_scores()`](https://mothur.org/strollur/reference/xdev_add_sequence_fastq_scores.md)
  function.
- Adds `strollur::xdev_alignment_length()` function
- Adds
  [`strollur::xdev_assign_bin_taxonomy_tidy()`](https://mothur.org/strollur/reference/xdev_assign_bin_taxonomy_tidy.md)
  function.
- Adds
  [`strollur::xdev_assign_sample_distances()`](https://mothur.org/strollur/reference/xdev_assign_sample_distances.md)
  function.
- Adds
  [`strollur::xdev_assign_sequence_quality_scores`](https://mothur.org/strollur/reference/xdev_assign_sequence_quality_scores.md)
  function.
- Adds
  [`strollur::xdev_has_bin_taxonomy()`](https://mothur.org/strollur/reference/xdev_has_bin_taxonomy.md)
  function.
- Adds
  [`strollur::xdev_get_sample_distances()`](https://mothur.org/strollur/reference/xdev_get_sample_distances.md)
  function.
- Adds
  [`strollur::xdev_get_sequence_indexes_by_sample()`](https://mothur.org/strollur/reference/xdev_get_sequence_indexes_by_sample.md)
  function.

Bug Fixes: -
[`strollur::summary()`](https://mothur.org/strollur/reference/summary.md)
will now display results that accurately represent the data.

Deprecated:

- Simplifies member functions of the
  [strollur_object](https://mothur.org/strollur/reference/strollur.html).
  - Removes `strollur_object$add_sequence_tree()` : Functionality now
    included in`strollur_object$add()`.
  - Removes `strollur_object$add_sample_tree()` : Functionality now
    included in`strollur_object$add()`.
  - Removes `strollur_object$get_sequence_tree()` : Functionality now
    included in`strollur_object$report()`.
  - Removes `strollur_object$get_sample_tree()` : Functionality now
    included in`strollur_object$report()`.

## strollur 0.1.3

CRAN release: 2026-07-02

- Adds
  [`strollur::read_mothur_oligos()`](https://mothur.org/strollur/reference/read_mothur_oligos.md)
  function for the import of oligo data from
  [mothur](https://mothur.org).
- Makes report format more flexible. You can now add reports that are
  not tied to sequence data. These can be general metadata, oligo
  information, paired fastq file tables and much more.

## strollur 0.1.2

CRAN release: 2026-06-24

- Initial CRAN release.
