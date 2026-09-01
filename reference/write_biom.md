# Write a [BIOM](https://biom-format.org) formatted file containing your [strollur](https://mothur.org/strollur/reference/strollur.html) objects data.

Write a [BIOM](https://biom-format.org) formatted file containing your
strollur::strollur objects data.

## Usage

``` r
write_biom(
  data,
  file_root = NULL,
  path = NULL,
  taxonomy = TRUE,
  sequence = TRUE
)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- file_root:

  a string containing the root name of the output file. Default =
  `{dataset_name}`.`{bin_type}`.biom

- path:

  string containing the name of directory where the files should be
  written. Default = current working directory.

- taxonomy:

  logical, when `TRUE` if consensus taxonomies are available include
  them in the biom file. Default = `TRUE`.

- sequence:

  logical, when `TRUE` if bins have a representative sequence or all
  bins contain a single feature sequence include the sequences in the
  biom file. Default = `TRUE`.

## Value

vector of BIOM file names. One for each bin type.

## Examples

``` r

if (requireNamespace("h5lite", quietly = TRUE)) {
  miseq <- strollur::miseq_sop_example()
  strollur::write_biom(miseq, tempfile())
} else {
  message(paste(
    "To use this functionality you have to install the",
    "h5lite package."
  ))
}
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 171 samples distances.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.
#> [1] "/tmp/RtmpUVSXMn/file3efd5c4ce97e.otu.biom"      
#> [2] "/tmp/RtmpUVSXMn/file3efd5c4ce97e.asv.biom"      
#> [3] "/tmp/RtmpUVSXMn/file3efd5c4ce97e.phylotype.biom"
```
