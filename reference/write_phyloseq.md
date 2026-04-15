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
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added metadata.
#> Added 2 resource references.
#> Added a contigs_report.
phylo_obj <- write_phyloseq(miseq)
```
