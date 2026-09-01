# sort_dataframe

Sort dataframe

## Usage

``` r
sort_dataframe(data, order, named_col)
```

## Arguments

- data:

  the data.frame to be sorted

- order:

  vector containing the order desired

- named_col:

  name of column in data.frame to match order

## Value

sorted data.frame

## Examples

``` r

# sort results alphabetically

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

sequence_names <- strollur::names(miseq)

fasta <- strollur::report(miseq, type = fasta)

sorted_fasta <- strollur::sort_dataframe(fasta,
  order = sort(sequence_names),
  named_col = "sequence_name"
)
```
