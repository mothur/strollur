# Create a phyloseq object from your [strollur object](https://mothur.org/strollur/reference/strollur.html)

The `strollur::write_phyloseq()` function will take any strollur object
and return it as a `phyloseq` object.

## Usage

``` r
write_phyloseq(data)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

## Value

returns a `phyloseq` object.

## References

McMurdie,P.J. and Holmes,S. (2013), phyloseq: An R Package for
Reproducible Interactive Analysis and Graphics of Microbiome Census
Data. PLoS ONE 8:e61217. <doi:10.1371/journal.pone.0061217>

## Examples

``` r
miseq <- strollur::miseq_sop_example()
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
if (requireNamespace("phyloseq", quietly = TRUE)) {
  phylo_obj <- strollur::write_phyloseq(miseq)
} else {
  message(paste(
    "To use this functionality you have to install the",
    "phyloseq package."
  ))
}
```
