# Find the number of sequences, samples, treatments or bins of a given type in a [strollur](https://mothur.org/strollur/reference/strollur.html) object

Find the number of sequences, samples, treatments or bins of a given
type in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
count(
  data,
  type = "sequences",
  bin_type = "otu",
  samples = NULL,
  distinct = FALSE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- type, :

  string containing the type of data you want the number of. Options
  include: "sequences", "samples", "treatments", "bins", and
  "references". Default = "sequences".

- bin_type, :

  string containing the bin type you would like the number of bins for.
  Default = "otu".

- samples, :

  vector of strings. samples is only used when 'type' = "sequences" or
  'type' = "bins" . samples should contain the names of the samples you
  want the count for. Default = NULL.

- distinct, :

  Boolean. distinct is used when 'type' = "sequences" or 'type' =
  "bins". When 'type' = "sequences" and distinct is TRUE the number of
  unique sequences is returned. When 'type' = "sequences" and distinct
  is FALSE the total number of sequences is returned. This can also be
  combined with samples to find the number of unique sequences found
  ONLY in a given set of samples, or to find the number of unique
  sequences in given set of samples that may also be present in other
  samples. When 'type' = "bins", you can set distinct = TRUE to return
  the number of bins that ONLY contain sequences from the given samples.
  When distinct is FALSE the count returned contains bins with sequences
  from a given samples, but those bins may also contain other samples.
  Default = FALSE.

## Value

double

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

# To get the total number of sequences
count(data = miseq, type = "sequences")
#> [1] 113963

# To get number of unique sequences
count(data = miseq, type = "sequences", distinct = TRUE)
#> [1] 2425

# To get number of unique sequences from samples 'F3D0' and 'F3D1'
# Note these sequences will be present in both samples but may be
# be present in other samples as well
count(data = miseq, type = "sequences", samples = c("F3D0", "F3D1"))
#> [1] 9385

# To get number of unique sequences exclusive to samples 'F3D0' and 'F3D1'
# Note sequences are present in both samples and NOT present in any other
# samples.
count(
  data = miseq, type = "sequences", samples = c("F3D0", "F3D1"),
  distinct = TRUE
)
#> [1] 2

# To get the number of samples in the dataset
count(data = miseq, type = "samples")
#> [1] 19

# To get the number of treatments in the dataset
count(data = miseq, type = "treatments")
#> [1] 2

# To get the number of "otu" bins in the dataset
count(data = miseq, type = "bins", bin_type = "otu")
#> [1] 531

# To get the number of "asv" bins in the dataset
count(data = miseq, type = "bins", bin_type = "asv")
#> [1] 2425

# To get the number of "phylotype" bins in the dataset
count(data = miseq, type = "bins", bin_type = "phylotype")
#> [1] 63

# To get number of "otu" bins from samples 'F3D0' and 'F3D1'
# Note these bins will have sequences from both samples but there may be
# other samples present as well
count(
  data = miseq,
  type = "bins", bin_type = "otu", samples = c("F3D0", "F3D1")
)
#> [1] 125

# To get number of "otu" bins unique to samples 'F3D0' and 'F3D1'
# Note these bins will have sequences from both samples and NO other samples
# will be present in the bins.
count(
  data = miseq, type = "bins", bin_type = "otu",
  samples = c("F3D0", "F3D1"), distinct = TRUE
)
#> [1] 1
```
