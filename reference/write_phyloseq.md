# write_phyloseq

The \`write_phyloseq()\` function will take any strollur object and
return it as a "phyloseq" object.

## Usage

``` r
write_phyloseq(data)
```

## Arguments

- data:

  the strollur object you created using one of the many read functions
  in this package.

## Value

returns a "phyloseq" object.

## Examples

``` r
miseq <- miseq_sop_example()
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
phylo_obj <- write_phyloseq(miseq)
```
