# strollur 0.1.4

New Features:

- Expands built in read / write functions.
     - `strollur::read_fastq()` returns data.frame with fastq data.
     - `strollur::write_fastq()` writes fastq data to file.
     - `strollur::read_quality()` returns data.frame with sequence quality data.
     - `strollur::write_quality()` writes sequence quality data to file.
     - `strollur::read_biom()` returns a [strollur_object](https://mothur.org/strollur/reference/strollur.html).
     - `strollur::write_biom()` writes biom files containing the [strollur_object](https://mothur.org/strollur/reference/strollur.html).
- Expands member functions of the [strollur_object](https://mothur.org/strollur/reference/strollur.html).
     - `strollur_object$alignment_length()` returns the sequences alignment length or -1 if the sequences are unaligned.
     - `strollur_object$has_bin_taxonomy()` returns `TRUE` if your dataset includes bin classifications.
     - `strollur_object$has_sequence_taxonomy()` returns `TRUE` if your dataset includes sequence classifications.
     - `strollur_object$rename()` rename the strollur object.

Enhancements:

- Adds `fastq`, `quality`, `sequence_tree`, `sample_tree`, `sample_distance` types to:
     - `strollur::report()`
     - `strollur_object$report()`
- Adds `fastq` type to:
     - `strollur::add()`
     - `strollur_object$add()`
- Adds `quality` and `sample_distance` type to:
     - `strollur::assign()`
     - `strollur_object$assign()`

Functions for Package developers

- Adds `strollur::xdev_alignment_length()` function 
- Adds `strollur::xdev_assign_bin_taxonomy_tidy()` function.
- Adds `strollur::xdev_has_bin_taxonomy()` function.
- Adds `strollur::xdev_get_sequence_indexes_by_sample()` function.
- Adds `strollur::xdev_add_sequence_fastq_scores()` function.
- Adds `strollur::xdev_get_sample_distances()` function.
- Adds `strollur::xdev_assign_sample_distances()` function.

Bug Fixes:
 - `strollur::summary()` will now display results that accurately represent the data.

Deprecated:

 - Simplifies member functions of the [strollur_object](https://mothur.org/strollur/reference/strollur.html). 
     - Removes `strollur_object$add_sequence_tree()` :  Functionality now included in`strollur_object$add()`.
     - Removes `strollur_object$add_sample_tree()`  :  Functionality now included in`strollur_object$add()`. 
     - Removes `strollur_object$get_sequence_tree()` :  Functionality now included in`strollur_object$report()`.
     - Removes `strollur_object$get_sample_tree()` :  Functionality now included in`strollur_object$report()`. 

# strollur 0.1.3

* Adds `strollur::read_mothur_oligos()` function for the import of oligo data from [mothur](https://mothur.org). 
* Makes report format more flexible. You can now add reports that are not tied to sequence data. These can be general metadata, oligo information, paired fastq file tables and much more.

# strollur 0.1.2

* Initial CRAN release.
