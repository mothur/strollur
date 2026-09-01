# Get indexes of sequences parsed by sample

As a developer you may want to process data by sample. To save space you
can request the indexes parsed by sample, and then get a single copy of
the sequences and sequence names.

## Usage

``` r
xdev_get_sequence_indexes_by_sample(data, samples = as.character(c()))
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- samples:

  a vector of strings containing the names of the samples you would like
  sequence names for. By default all samples are included.

## Value

2D vector of strings indexes for use with xdev_get_sequences and
xdev_get_names. requested parsed by sample. (Indexes start at 1)

## Examples

``` r

data <- strollur::miseq_sop_example()
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

sequences <- strollur::xdev_get_sequences(data)
names <- strollur::xdev_names(data)

# To get the indexes of the names and sequences by sample
indexes <- strollur::xdev_get_sequence_indexes_by_sample(data)

# First sequence in first sample
names[indexes[[1]][1]]
#> [1] "M00967_43_000000000-A3JHG_1_2103_25452_6018"
sequences[indexes[[1]][1]]
#> [1] "TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-T-T--GG-GT--TT-A-AA-GG-GTGGC-G-CA-GGC-G-G-G-AT-G-C-C--A-G-T-C-A-G-C-G-G--TC-A-AA-TT-T-C-GG-GG--CT-C-AA-C-C-C-C-G-A-C--CT-G-C-CGTT-GAAAC-TG-G-TGTCC-TAGA-GT-GG-GC-GA-G-A---AG-T-A-TGCGGAATGCGTGGTGT-AGCGGT-GAAATGCATAG-AT-A-TC-AC-GC-AG-AACTCCGAT-TGCGAAGGCA------GCATA-CCG-G-CG-CC-C-G-ACTGACG-CTCA-TGCA-CGAAA-GCG-TGGGT-ATC-GAACAGG"
```
