# clear

Clear data from a [dataset](dataset.md) object

## Usage

``` r
clear(data, tags = as.character(c()))
```

## Arguments

- data, :

  a [dataset](dataset.md) object

- tags:

  a vector of strings containing the items you wish to clear. Options
  are 'sequence_data', 'bin_data', 'metadata', 'references',
  'sequence_tree', 'sample_tree' and 'reports'. By default, everything
  is cleared.

## Examples

``` r
data <- miseq_sop_example()
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 19 samples to treatments.
#> ℹ Assigned 531 otu bin taxonomies.
#> ℹ Assigned 531 otu bin representative sequences.
#> ℹ Added metadata.
#> ℹ Added 2 resource references.
#> ℹ Added a contigs_report.
clear(data)
#> → `` is not a valid item to clear, ignoring.
```
