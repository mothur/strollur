# sort_dataframe

Sort dataframe

## Usage

``` r
sort_dataframe(data, order, named_col)
```

## Arguments

- data, :

  the data.frame to be sorted

- order, :

  vector containing the order desired

- named_col, :

  name of column in data.frame to match order

  \# sort results alphabetically

  miseq \<- miseq_sop_example()

  sequence_names \<- names(miseq)

  fasta \<- report(miseq, type = fasta)

  sorted_fasta \<- sort_dataframe(fasta, order = sort(sequence_names),
  named_col = "sequence_names")

## Value

sorted data.frame
