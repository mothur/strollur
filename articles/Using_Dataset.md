# Using_Dataset

## Overview

The *rdataset* package stores the data associated with your microbial
DNA analysis. This tutorial will familiarize you with some of the
functions available in the *rdataset* package. If you haven’t reviewed
the “Getting Started” tuturial, we recommend you start there.

First let’s load the rdataset package.

``` r
library(rdataset)
#> 
#> Attaching package: 'rdataset'
#> The following objects are masked from 'package:base':
#> 
#>     assign, names, summary
```

## Loading the example dataset

We can use the **miseq_sop_example()** function to create a dataset
object from the [Miseq SOP Example](https://mothur.org/wiki/miseq_sop/).

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
```

To print a summary of the miseq dataset, run the following:

``` r
miseq
#> miseq_sop:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0   2850.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28491.750
#> Median:          1  375 252.0000      0 4.000000     0  56982.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85473.250
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 111114.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.4472      0 4.368699     0      0.000
#> Unique seqs:  2425 
#> Total seqs:   113963 
#> 
#> contigs_report :
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000   0.000000      0    248.0000
#> 2.5%-tile:     0.0016101600 252.0000   0.000000      0    251.0000
#> 25%-tile:      0.0028177700 252.0000   0.000000      0    251.0000
#> Median:        0.0062948698 252.0000   2.000000      0    251.0000
#> 75%-tile:      0.0264780000 253.0000   4.000000      0    251.0000
#> 97.5%-tile:    1.1412299871 253.0000  76.000000      0    252.0000
#> Maximum:       3.0126200000 270.0000 120.000000      0    256.0000
#> Mean:          0.0984788569 252.5128   7.534147      0    251.0762
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000      0.000000
#> 2.5%-tile:        248.0000      1.000000
#> 25%-tile:         249.0000      2.000000
#> Median:           249.0000      2.000000
#> 75%-tile:         249.0000      2.000000
#> 97.5%-tile:       251.0000      3.000000
#> Maximum:          255.0000     22.000000
#> Mean:             249.2136      1.862692
#> Unique seqs:  2425 
#> Total seqs:   113963 
#> 
#> 
#> Sample   Total:
#> F3D0 6191 
#> F3D1 4652 
#> F3D141   4656 
#> F3D142   2423 
#> F3D143   2403 
#> F3D144   3449 
#> F3D145   5532 
#> F3D146   3831 
#> F3D147   12430 
#> F3D148   9465 
#> F3D149   10014 
#> F3D150   4126 
#> F3D2 15686 
#> F3D3 5199 
#> F3D5 3469 
#> F3D6 6394 
#> F3D7 4055 
#> F3D8 4253 
#> F3D9 5735 
#> 
#> Treatment   Total:
#> Early    55634 
#> Late 58329 
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> Total number of otus: 531 
#> Total number of asvs: 2425 
#> Total number of phylotypes: 63
```

## Accessing Data

The *rdataset* package includes several functions to access your data.

**names()** - The names() function is an extension of base R’s names
function. It allows you to get the name of your dataset or the names of
the sequences, bins, samples, treatments and reports in your dataset.

**count()** - The count() function allows you to get the number of
sequences, samples, treatments or bins in your dataset.

**abundance()** - The abundance() function can be used to access the
abundance data for sequences, bins, samples, and treatments.

**report()** - The report() function can be used to access the various
data reports associated with your dataset. If you have added custom
reports for alignment, contigs_assembly or chimeras, you can get those
as well using the report_type name you provided in the add() function.

**summary()** - The summary() function allows you to summarize
sequences, your custom reports, and scrapped data.

### names() and count()

Let’s take a closer look at the names() and count() functions.

- **Parameters:**
  - *data* - a dataset object
  - *type* - string containing the type of data you would like.
  - *bin_type* - string containing the bin type you would like data for.
  - *samples* - vector of strings containing sample names. samples is
    only used when ‘type’ = “sequences” or ‘type’ = “bins” .
  - *distinct* - boolean. ‘distinct’ can be used when ‘type’ is
    ‘sequences’ or ‘type’ is ‘bins’. When ‘type’ = ‘sequences’ and
    distinct is TRUE the number of unique sequences is returned. When
    ‘type’ = ‘sequences’ and distinct is FALSE the total number of
    sequences is returned. When ‘type’ = “bins”, you can set distinct =
    TRUE to return the number of bins that ONLY contain sequences from
    the given samples. When distinct is FALSE the count returned
    contains bins with sequences from a given samples, but those bins
    may also contain other samples.

To get the name of the dataset, set the type parameter to ‘dataset’:

``` r
names(data = miseq, type = "dataset")
#> [1] "miseq_sop"
```

To get the names of the sequences in your dataset, set the type
parameter to ‘sequences’:

``` r
names(data = miseq, type = "sequences") |> head(n = 5)
#> [1] "M00967_43_000000000-A3JHG_1_2101_16474_12783"
#> [2] "M00967_43_000000000-A3JHG_1_1113_12711_3318" 
#> [3] "M00967_43_000000000-A3JHG_1_2108_14707_9807" 
#> [4] "M00967_43_000000000-A3JHG_1_1110_4126_16552" 
#> [5] "M00967_43_000000000-A3JHG_1_2102_8408_13436"
```

To find the number of sequences in your dataset, you can use the count()
function. Set the type to ‘sequences’. Take special note of the
‘distinct’ parameter.

To find the total number of sequences in the dataset, run the following:

``` r
count(data = miseq, type = "sequences", distinct = FALSE)
#> [1] 113963
```

To find the number of unique sequences in the dataset, set ‘distinct’ to
TRUE:

``` r
count(data = miseq, type = "sequences", distinct = TRUE)
#> [1] 2425
```

To get the names of the unique sequences, run the following:

``` r
names(data = miseq, type = "sequences") |> head(n = 5)
#> [1] "M00967_43_000000000-A3JHG_1_2101_16474_12783"
#> [2] "M00967_43_000000000-A3JHG_1_1113_12711_3318" 
#> [3] "M00967_43_000000000-A3JHG_1_2108_14707_9807" 
#> [4] "M00967_43_000000000-A3JHG_1_1110_4126_16552" 
#> [5] "M00967_43_000000000-A3JHG_1_2102_8408_13436"
```

To get number of sequences from sample ‘F3D0’, you can set the samples
parameter. Note, these sequences will be present in the sample but may
be be present in other samples as well.

``` r
count(data = miseq, type = "sequences", samples = c("F3D0"))
#> [1] 6191
```

To get the names of the sequences present in sample ‘F3D0’

``` r
names(data = miseq, type = "sequences", samples = c("F3D0")) |> head(n = 5)
#> [1] "M00967_43_000000000-A3JHG_1_2103_25452_6018" 
#> [2] "M00967_43_000000000-A3JHG_1_1109_13330_21597"
#> [3] "M00967_43_000000000-A3JHG_1_1110_5315_13833" 
#> [4] "M00967_43_000000000-A3JHG_1_2104_26311_10309"
#> [5] "M00967_43_000000000-A3JHG_1_1101_9620_19745"
```

To get number of unique sequences exclusive to sample ‘F3D0’ *Note
sequences are present in the sample and NOT present in any other
samples.*

``` r
count(data = miseq, type = "sequences", samples = c("F3D0"), distinct = TRUE)
#> [1] 101
```

To get the names of the sequences unique to sample ‘F3D0’

``` r
names(data = miseq, type = "sequences", samples = c("F3D0"), distinct = TRUE) |>
  head(n = 5)
#> [1] "M00967_43_000000000-A3JHG_1_2103_25452_6018" 
#> [2] "M00967_43_000000000-A3JHG_1_1101_9620_19745" 
#> [3] "M00967_43_000000000-A3JHG_1_2109_17345_6668" 
#> [4] "M00967_43_000000000-A3JHG_1_1114_14431_2336" 
#> [5] "M00967_43_000000000-A3JHG_1_2106_14305_11884"
```

To get the number of “otu” bins in the dataset

``` r
count(data = miseq, type = "bins", bin_type = "otu")
#> [1] 531
```

To get the names of the “otu” bins

``` r
names(data = miseq, type = "bins", bin_type = "otu") |> head(n = 5)
#> [1] "Otu001" "Otu002" "Otu003" "Otu004" "Otu005"
```

To get the number of “asv” bins in the dataset

``` r
count(data = miseq, type = "bins", bin_type = "asv")
#> [1] 2425
```

To get the names of the “asv” bins

``` r
names(data = miseq, type = "bins", bin_type = "asv") |> head(n = 5)
#> [1] "Asv0001" "Asv0002" "Asv0003" "Asv0004" "Asv0005"
```

To get the number of “phylotype” bins in the dataset

``` r
count(data = miseq, type = "bins", bin_type = "phylotype")
#> [1] 63
```

To get the names of the “phylotype” bins

``` r
names(data = miseq, type = "bins", bin_type = "phylotype") |> head(n = 5)
#> [1] "Phylo01" "Phylo02" "Phylo03" "Phylo04" "Phylo05"
```

To get number of “otu” bins from sample ‘F3D0’ *Note these bins will
have sequences from the sample but there may be other* *samples present
in the bins as well*

``` r
count(data = miseq, type = "bins", bin_type = "otu", samples = c("F3D0"))
#> [1] 191
```

To get the names of the “otu” bins that include sequences from ‘F3D0’

``` r
names(
  data = miseq,
  type = "bins", samples = c("F3D0"), distinct = FALSE
) |> head(n = 5)
#> [1] "Otu001" "Otu002" "Otu003" "Otu004" "Otu005"
```

To get number of “otu” bins unique to samples ‘F3D0’ *Note these bins
will have sequences from the sample and NO other samples will* *be
present in the bins.*

``` r
count(
  data = miseq, type = "bins",
  bin_type = "otu", samples = c("F3D0"), distinct = TRUE
)
#> [1] 14
```

To get the names of the “otu” bins that are unique to ‘F3D0’

``` r
names(
  data = miseq,
  type = "bins", samples = c("F3D0"), distinct = TRUE
) |> head(n = 5)
#> [1] "Otu330" "Otu339" "Otu341" "Otu345" "Otu347"
```

To get the number of samples in the dataset

``` r
count(data = miseq, type = "samples")
#> [1] 19
```

To get the names of the samples

``` r
names(data = miseq, type = "samples")
#>  [1] "F3D0"   "F3D1"   "F3D141" "F3D142" "F3D143" "F3D144" "F3D145" "F3D146"
#>  [9] "F3D147" "F3D148" "F3D149" "F3D150" "F3D2"   "F3D3"   "F3D5"   "F3D6"  
#> [17] "F3D7"   "F3D8"   "F3D9"
```

To get the number of treatments in the dataset

``` r
count(data = miseq, type = "treatments")
#> [1] 2
```

To get the names of the treatments

``` r
names(data = miseq, type = "treatments")
#> [1] "Early" "Late"
```

To get the names of the reports

``` r
names(data = miseq, type = "reports")
#> [1] "contigs_report" "sequence_data"
```

Now that we are familar with the names() and count() functions. Let’s
take a closer look at how we can use the **abundance()** function.

``` r
# To the total abundance for each sequence
abundance(data = miseq, type = "sequences") |> head(n = 5)
#>                                 sequence_names abundances
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783          1
#> 2  M00967_43_000000000-A3JHG_1_1113_12711_3318          1
#> 3  M00967_43_000000000-A3JHG_1_2108_14707_9807          1
#> 4  M00967_43_000000000-A3JHG_1_1110_4126_16552          1
#> 5  M00967_43_000000000-A3JHG_1_2102_8408_13436          1

# To the total abundance for each sequence parsed by sample
abundance(data = miseq, type = "sequences", by_sample = TRUE) |> head(n = 5)
#>                                 sequence_names abundances samples treatments
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783          1  F3D150       Late
#> 2  M00967_43_000000000-A3JHG_1_1113_12711_3318          1  F3D142       Late
#> 3  M00967_43_000000000-A3JHG_1_2108_14707_9807          1    F3D3      Early
#> 4  M00967_43_000000000-A3JHG_1_1110_4126_16552          1    F3D8      Early
#> 5  M00967_43_000000000-A3JHG_1_2102_8408_13436          1    F3D7      Early

# To the total abundance for each "otu" bin
abundance(data = miseq, type = "bins", bin_type = "otu") |> head(n = 5)
#>   otu_id abundance
#> 1 Otu001     12288
#> 2 Otu002      8892
#> 3 Otu003      7794
#> 4 Otu004      7476
#> 5 Otu005      7450

# To the total abundance for each "otu" bin parsed by sample
abundance(
  data = miseq,
  type = "bins", bin_type = "otu", by_sample = TRUE
) |> head(n = 5)
#>   bin_names abundances samples treatments
#> 1    Otu001        499    F3D0      Early
#> 2    Otu001        351    F3D1      Early
#> 3    Otu001        388  F3D141       Late
#> 4    Otu001        244  F3D142       Late
#> 5    Otu001        189  F3D143       Late

# To the total abundance for each "asv" bin
abundance(data = miseq, type = "bins", bin_type = "asv") |> head(n = 5)
#>    asv_id abundance
#> 1 Asv0001     12196
#> 2 Asv0002      8829
#> 3 Asv0003      7698
#> 4 Asv0004      7436
#> 5 Asv0005      7307

# To the total abundance for each "asv" bin parsed by sample
abundance(
  data = miseq,
  type = "bins", bin_type = "asv", by_sample = TRUE
) |> head(n = 5)
#>   bin_names abundances samples treatments
#> 1   Asv0001        495    F3D0      Early
#> 2   Asv0001        340    F3D1      Early
#> 3   Asv0001        386  F3D141       Late
#> 4   Asv0001        242  F3D142       Late
#> 5   Asv0001        188  F3D143       Late

# To the total abundance for each sample
abundance(data = miseq, type = "samples") |> head(n = 5)
#>   samples abundances
#> 1    F3D0       6191
#> 2    F3D1       4652
#> 3  F3D141       4656
#> 4  F3D142       2423
#> 5  F3D143       2403

# To the total abundance for each treatment
abundance(data = miseq, type = "treatments") |> head(n = 5)
#>   treatments abundances
#> 1      Early      55634
#> 2       Late      58329
```

The **report()** and **summary()** functions allow you to access
sequence and classification reports, metadata, resource references,
scrapped data reports, sequence data summaries, custom report summaries
and scrapped data summaries. Let’s take a closer look using the miseq
example.

``` r
# To get a report about the FASTA sequence data
report(data = miseq, type = "sequence_data") |> head(n = 5)
#> 
#> Your dataset does not include a report named sequence_data, ignoring request.
#> data frame with 0 columns and 0 rows

# To summarize the FASTA sequence data
summary(data = miseq, type = "sequences", verbose = FALSE)
#>             starts ends   nbases ambigs polymers numns    numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.000
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0   2850.075
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28491.750
#> Median:          1  375 252.0000      0 4.000000     0  56982.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85473.250
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 111114.925
#> Maximum:         1  375 256.0000      0 6.000000     0 113963.000
#> Mean:            1  375 252.4472      0 4.368699     0      0.000

# To get a classification report about your sequence data
report(data = miseq, type = "sequence_taxonomy") |> head(n = 7)
#>                                             id level
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783     1
#> 2 M00967_43_000000000-A3JHG_1_2101_16474_12783     2
#> 3 M00967_43_000000000-A3JHG_1_2101_16474_12783     3
#> 4 M00967_43_000000000-A3JHG_1_2101_16474_12783     4
#> 5 M00967_43_000000000-A3JHG_1_2101_16474_12783     5
#> 6 M00967_43_000000000-A3JHG_1_2101_16474_12783     6
#> 7  M00967_43_000000000-A3JHG_1_1113_12711_3318     1
#>                               taxon confidence
#> 1                          Bacteria        100
#> 2                   "Bacteroidetes"        100
#> 3                     "Bacteroidia"         99
#> 4                   "Bacteroidales"         99
#> 5              "Porphyromonadaceae"         88
#> 6 "Porphyromonadaceae"_unclassified         88
#> 7                          Bacteria        100

# To get the consensus taxonomies assigned to your 'otu' bins
report(data = miseq, type = "bin_taxonomy", bin_type = "otu") |> head(n = 7)
#>       id level                             taxon confidence
#> 1 Otu001     1                          Bacteria        100
#> 2 Otu001     2                   "Bacteroidetes"        100
#> 3 Otu001     3                     "Bacteroidia"        100
#> 4 Otu001     4                   "Bacteroidales"        100
#> 5 Otu001     5              "Porphyromonadaceae"        100
#> 6 Otu001     6 "Porphyromonadaceae"_unclassified        100
#> 7 Otu002     1                          Bacteria        100

# To get the consensus taxonomies assigned to your 'asv' bins
report(data = miseq, type = "bin_taxonomy", bin_type = "asv") |> head(n = 7)
#>        id level                             taxon confidence
#> 1 Asv0001     1                          Bacteria        100
#> 2 Asv0001     2                   "Bacteroidetes"        100
#> 3 Asv0001     3                     "Bacteroidia"        100
#> 4 Asv0001     4                   "Bacteroidales"        100
#> 5 Asv0001     5              "Porphyromonadaceae"        100
#> 6 Asv0001     6 "Porphyromonadaceae"_unclassified        100
#> 7 Asv0002     1                          Bacteria        100

# To get the consensus taxonomies assigned to your 'phylotype' bins
report(data = miseq, type = "bin_taxonomy", bin_type = "phylotype") |>
  head(n = 7)
#>        id level                        taxon confidence
#> 1 Phylo01     1                     Bacteria        100
#> 2 Phylo01     2                   Firmicutes        100
#> 3 Phylo01     3                   Clostridia        100
#> 4 Phylo01     4                Clostridiales        100
#> 5 Phylo01     5              Lachnospiraceae        100
#> 6 Phylo01     6 Lachnospiraceae_unclassified        100
#> 7 Phylo02     1                     Bacteria        100

# To get the scrapped sequences report
# Note, the miseq example does not include scrapped data
report(data = miseq, type = "sequence_scrap")
#> data frame with 0 columns and 0 rows

# To get the "otu" scrapped bins report
# Note, the miseq example does not include scrapped data
report(data = miseq, type = "bin_scrap", bin_type = "otu")
#> data frame with 0 columns and 0 rows

# To get a summary of the data scrapped during your analysis
# Note, the miseq example does not include scrapped data
summary(data = miseq, type = "scrap", verbose = FALSE)
#> data frame with 0 columns and 0 rows

# To get the metadata associated with your dataset
report(data = miseq, type = "metadata") |> head(n = 7)
#>   sample days_post_wean
#> 1   F3D0              0
#> 2   F3D1              1
#> 3 F3D141            141
#> 4 F3D142            142
#> 5 F3D143            143
#> 6 F3D144            144
#> 7 F3D145            145

# To get the resource references associated with your dataset
report(data = miseq, type = "references")
#>            reference_names reference_versions         reference_usages
#> 1 trainset9_032012.pds.zip                 NA classification by mothur
#> 2           silva.v4.fasta             1.38.1                alignment
#>                                                              reference_notes
#> 1                                                                         NA
#> 2 custom reference created by trimming silva.bacteria.fasta to the V4 region
#>                                                            reference_urls
#> 1 https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip
#> 2                          https://mothur.org/wiki/silva_reference_files/

# To access the custom reports, first let's find the names of the reports
names(data = miseq, type = "reports")
#> [1] "contigs_report" "sequence_data"

# To access our contigs assembly report
report(data = miseq, type = "contigs_report") |> head(n = 7)
#>                                           Name Length Overlap_Length
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783    253            250
#> 2  M00967_43_000000000-A3JHG_1_1113_12711_3318    253            249
#> 3  M00967_43_000000000-A3JHG_1_2108_14707_9807    253            249
#> 4  M00967_43_000000000-A3JHG_1_1110_4126_16552    252            249
#> 5  M00967_43_000000000-A3JHG_1_2102_8408_13436    253            249
#> 6 M00967_43_000000000-A3JHG_1_1107_22580_21773    252            250
#> 7 M00967_43_000000000-A3JHG_1_1108_14299_17220    253            249
#>   Overlap_Start Overlap_End MisMatches Num_Ns Expected_Errors
#> 1             2         252         19      0      0.29461400
#> 2             2         251          0      0      0.00183396
#> 3             2         251          0      0      0.00196774
#> 4             2         251          4      0      0.05629750
#> 5             2         251          0      0      0.00259554
#> 6             1         251          8      0      0.05068300
#> 7             2         251          0      0      0.00215398

# To get a summary of your contigs assembly report
summary(
  data = miseq, type = "reports",
  report_type = "contigs_report", verbose = FALSE
)
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000   0.000000      0    248.0000
#> 2.5%-tile:     0.0016101600 252.0000   0.000000      0    251.0000
#> 25%-tile:      0.0028177700 252.0000   0.000000      0    251.0000
#> Median:        0.0062948698 252.0000   2.000000      0    251.0000
#> 75%-tile:      0.0264780000 253.0000   4.000000      0    251.0000
#> 97.5%-tile:    1.1412299871 253.0000  76.000000      0    252.0000
#> Maximum:       3.0126200000 270.0000 120.000000      0    256.0000
#> Mean:          0.0984788569 252.5128   7.534147      0    251.0762
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000      0.000000
#> 2.5%-tile:        248.0000      1.000000
#> 25%-tile:         249.0000      2.000000
#> Median:           249.0000      2.000000
#> 75%-tile:         249.0000      2.000000
#> 97.5%-tile:       251.0000      3.000000
#> Maximum:          255.0000     22.000000
#> Mean:             249.2136      1.862692
```

## General

**copy_dataset()** - Create a new dataset from an existing dataset.

**save_dataset()** - Save a dataset as a .rds file.

**load_dataset()** - Create a new dataset from a .rds file.

**export_dataset()** - Create a human readable list containing the
dataset data.

**import_dataset()** - Create a new dataset from an exported dataset.

**get_bin_types()** - Returns the types of bins in your dataset

**is_aligned()** - Returns TRUE if sequence data is aligned

``` r
# To create a new dataset that is a copy of miseq
copy_of_miseq <- copy_dataset(data = miseq)

# To save the dataset as a .rds file
save_dataset(data = miseq, file = "miseq.rds")
#> [1] "miseq.rds"

# To create a dataset from a .rds file
data <- load_dataset(file = "miseq.rds")

# To export the data into a human readable format
table <- export_dataset(data = miseq)

# For an overview of the exported data.frames
names(table)
#>  [1] "sequence_data"                      "sequence_report"                   
#>  [3] "sequence_abundance_table"           "otu_bin_data"                      
#>  [5] "otu_sequence_bin_assignments"       "otu_bin_representative_sequences"  
#>  [7] "asv_bin_data"                       "asv_sequence_bin_assignments"      
#>  [9] "phylotype_bin_data"                 "phylotype_sequence_bin_assignments"
#> [11] "references"                         "metadata"                          
#> [13] "contigs_report"                     "sequence_tree"                     
#> [15] "sample_tree"

# To create a new dataset from the exported table
data <- import_dataset(table = table)
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 531 otu bin representative sequences.
#> ℹ Assigned 531 otu bin taxonomies.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 2425 asv bin taxonomies.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 63 phylotype bin taxonomies.
#> ℹ Added metadata.
#> ℹ Added 2 resource references.

# To find the types of bins in your dataset
get_bin_types(data = miseq)
#> [1] "otu"       "asv"       "phylotype"

# To determine if the FASTA data is aligned
is_aligned(data = miseq)
#> [1] TRUE
```

## Sample Trees and Sequence Trees

\#TODO
