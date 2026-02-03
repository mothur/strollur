# Functions for package developers

*rdataset* includes several functions designed for package developers.
These functions provide additional access to the back end data and allow
developers to merge, remove and set data directly. The function names
begin with *‘xdev\_’* or *‘xint\_’* to clearly indicate that they are
designed for development or for internal use, and not for the general
user.

## Setting functions

Over the course of an analysis your package may want to change the
abundances of sequences, modify sequence nucleotide strings due to
alignment, screening or filtering, and change the data set name. Let’s
take a close look at functions you will need to do so using the
[`miseq_sop_example()`](../reference/miseq_sop_example.md) data set.

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

### Setting the dataset name

To change the data set name you can use the
[`xdev_set_dataset_name()`](../reference/xdev_set_dataset_name.md)
function.

``` r
names(data = miseq, type = "dataset")
#> [1] "miseq_sop"
xdev_set_dataset_name(data = miseq, dataset_name = "modified_miseq")
names(data = miseq, type = "dataset")
#> [1] "modified_miseq"
```

### Setting nucleotide sequence strings

The miseq example data set is already aligned, filtered and screened,
but for the sake of example let’s set the nucleotide strings of
sequences *exclusive* to sample ‘F3D0’ to ‘NNNN’. First we will get the
names of the sequences only present in sample ‘F3D0’.

``` r
f3d0_names <- names(
  data = miseq,
  type = "sequences",
  samples = c("F3D0"),
  distinct = TRUE
)
length(f3d0_names)
#> [1] 101
```

Now we will assign nucleotide strings of the sequences *exclusive* to
sample ‘F3D0’ to ‘NNNN’. Then we can use the
[`xdev_get_by_sample()`](../reference/xdev_get_by_sample.md) get all the
sequences in present in sample ‘F3D0’ for closer inspection.

``` r
nnnn_neucleotides <- rep("NNNN", length(f3d0_names))
comments <- rep("example_set_sequences", length(f3d0_names))

xdev_set_sequences(
  data = miseq, sequence_names = f3d0_names,
  sequences = nnnn_neucleotides, comments = comments
)

f3d0_sequences <- xdev_get_by_sample(
  data = miseq,
  type = "sequences", samples = c("F3D0")
)
f3d0_sequences[[1]][300:305]
#> [1] "TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-T-TA-T-G-C-AA-G-T-C-A-G-A-A-G--TG-A-AA-GC-C-C-AA-AG--CT-C-AA-C-T-T-C-G-G-G-ACT-G-C-TTTT-GAAAC-TG-T-GTAAC-TAGA-GT-GC-AG-GA-G-G---GG-T-A-AGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCGGC-GGCGAAGGCG------GCTTA-CTG-G-AC-TG-T-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG"
#> [2] "NNNN"                                                                                                                                                                                                                                                                                                                                                                                   
#> [3] "TAC--GG-AG-GGT--GCA-A-G-C-G-T-T--AA-T-CGG-AA--TT-A-C-T--GG-GC--GT-A-AA-GC-GC-AC-G-CA-GGC-G-G-T-TT-G-T-T-AA-G-T-C-A-G-A-T-G--TG-A-AA-TC-C-C-CG-GG--CT-C-AA-C-C-T-G-G-G-A-ACT-G-C-ATCT-GATAC-TG-G-CAAGC-TTGA-GT-CT-CG-TA-G-A---GG-G-G-GGTAGAATTCCAGGTGT-AGCGGT-GAAATGCGTAG-AG-A-TC-TG-GA-GG-AATACCGGT-GGCGAAGGCG------GCCCC-CTG-G-AC-GA-A-G-ACTGACG-CTCA-GGTG-CGAAA-GCG-TGGGG-AGC-AAACAGG"
#> [4] "NNNN"                                                                                                                                                                                                                                                                                                                                                                                   
#> [5] "TAC--GT-AG-GGA--GCG-A-G-C-G-T-T--GT-C-CGG-AA--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GGC-G-G-G-AT-T-G-C-AA-G-T-T-G-G-A-T-G--TG-A-AA-AC-T-G-CG-GG--CT-C-AA-C-C-C-G-G-A-G-AGT-G-C-ATTC-AAAAC-TG-C-GATTC-TTGA-GT-GA-AG-TA-G-A---GG-C-A-AGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCTTG-CTG-G-GC-TT-T-T-ACTGACG-CTGA-GGCT-CGAAA-GTG-TGGGG-AGC-AAACAGG"
#> [6] "TAC--GT-AG-GGG--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-C-GT-C-G-C-AA-G-T-C-T-G-A-A-G--TG-A-AA-GC-C-C-GT-GG--CT-C-AA-C-C-G-C-G-G-A-ACC-G-C-TTTG-GAAAC-TG-C-GAGGC-TGGA-GT-GC-TG-GA-G-A---GG-T-A-AGCGGAATTCCTGGTGT-AGCGGT-GAAATGCGTAG-AT-A-TC-AG-GA-GG-AACACCGGT-GGCGAAGGCG------GCTTA-CTG-G-AC-AG-T-G-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-GAACAGG"
```

### Setting sequence abundances

Now that we know how to set sequence strings, let’s learn how to set
sequence abundances. In most cases you can use the
[`assign()`](../reference/assign.md) function to set sequence
abundances, but there may be cases where you want more flexibility. The
assign function requires abundances to be provided for each sequence. If
you only need to update a subset of the sequences, there are two
functions to that allow you to do that,
[`xdev_set_abundances()`](../reference/xdev_set_abundances.md) and
[`xdev_set_abundance()`](../reference/xdev_set_abundance.md).
*xdev_set_abundances* is used with data sets that include samples,
*xdev_set_abundance* is used for data sets without samples. Let’s set
the abundances of all sequences *exclusive* to sample F3D0 to 0, in
essence removing them. First we will get the abundance table for the
whole data set and create the inputs for
[`xdev_set_abundances()`](../reference/xdev_set_abundances.md).

``` r
abundance_table <- abundance(
  data = miseq,
  type = "sequences", by_sample = TRUE
)
head(abundance_table, n = 10)
#>                                  sequence_names abundances samples treatments
#> 1  M00967_43_000000000-A3JHG_1_2101_16474_12783          1  F3D150       Late
#> 2   M00967_43_000000000-A3JHG_1_1113_12711_3318          1  F3D142       Late
#> 3   M00967_43_000000000-A3JHG_1_2108_14707_9807          1    F3D3      Early
#> 4   M00967_43_000000000-A3JHG_1_1110_4126_16552          1    F3D8      Early
#> 5   M00967_43_000000000-A3JHG_1_2102_8408_13436          1    F3D7      Early
#> 6  M00967_43_000000000-A3JHG_1_1107_22580_21773          1    F3D3      Early
#> 7  M00967_43_000000000-A3JHG_1_1108_14299_17220         22  F3D146       Late
#> 8  M00967_43_000000000-A3JHG_1_1108_14299_17220         19  F3D147       Late
#> 9  M00967_43_000000000-A3JHG_1_1108_14299_17220         12  F3D148       Late
#> 10 M00967_43_000000000-A3JHG_1_1108_14299_17220          9  F3D149       Late

num_samples <- count(data = miseq, type = "samples")
new_abunds <- rep(list(rep(0, num_samples)), length(f3d0_names))

xdev_set_abundances(
  data = miseq,
  sequence_names = f3d0_names,
  abundances = new_abunds,
  reason = "F3D0_exclusive_sequences"
)
miseq
#> modified_miseq:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns   numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.00
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0   2847.35
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28464.50
#> Median:          1  375 252.0000      0 4.000000     0  56928.00
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85391.50
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 111008.65
#> Maximum:         1  375 256.0000      0 6.000000     0 113854.00
#> Mean:            1  375 252.4468      0 4.368533     0      0.00
#> Unique seqs:  2324 
#> Total seqs:   113854 
#> 
#> contigs_report :
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000   0.000000      0    248.0000
#> 2.5%-tile:     0.0016101600 252.0000   0.000000      0    251.0000
#> 25%-tile:      0.0028177700 252.0000   0.000000      0    251.0000
#> Median:        0.0062948698 252.0000   2.000000      0    251.0000
#> 75%-tile:      0.0261963997 253.0000   4.000000      0    251.0000
#> 97.5%-tile:    1.1412299871 253.0000  76.000000      0    252.0000
#> Maximum:       3.0126200000 270.0000 120.000000      0    256.0000
#> Mean:          0.0985028379 252.5124   7.536301      0    251.0762
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000      0.000000
#> 2.5%-tile:        248.0000      1.000000
#> 25%-tile:         249.0000      2.000000
#> Median:           249.0000      2.000000
#> 75%-tile:         249.0000      2.000000
#> 97.5%-tile:       251.0000      3.000000
#> Maximum:          255.0000     22.000000
#> Mean:             249.2138      1.862447
#> Unique seqs:  2324 
#> Total seqs:   113854 
#> 
#> scrap_summary:
#>        type                trash_code unique total
#> 1  sequence  F3D0_exclusive_sequences    101   109
#> 2       otu F3D0_exclusive_sequences,     14    14
#> 3       asv F3D0_exclusive_sequences,    101   109
#> 4 phylotype F3D0_exclusive_sequences,      2     2
#> 
#> Sample   Total:
#> F3D0 6082 
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
#> Early    55525 
#> Late 58329 
#> 
#> Number of unique seqs: 2324 
#> Total number of seqs: 113854 
#> Total number of otus: 517 
#> Total number of asvs: 2324 
#> Total number of phylotypes: 61
```

You can see that not all of the sequences in ‘F3D0’ are removed. That is
because the majority of sequences in the data set are present in
multiple samples. mothur2 uses xdev_set_abundances when running the
pre_cluster() function, sequences are clustered by sample, and the
abundances are set accordingly.

## Merging Functions

In the course of analyzing your data, there may be times when you want
to merge sequence abundances. When sequence nucleotide string are
identical you can save time processing by merging the identical read’s
abundances. mothur2 merges sequences abundances in the unique_seqs
function. The miseq example includes merged sequences that have been
assigned to bins. In practice, you would would not merge sequences after
assigning them to bins. Note, the code below is set to ‘eval=FALSE’
because you can only merge sequences assigned to the same bin and it
will throw errors accordingly. Let’s take a look at how to use the
[`xdev_merge_sequences()`](../reference/xdev_merge_sequences.md)
function for your reference.

``` r
random_sequence_names <- sample(names(data = miseq),
  size = 100, replace = FALSE
)

xdev_merge_sequences(
  data = miseq,
  sequence_names = random_sequence_names,
  reason = "merge_sequences_example"
)
```

Similarly, you may want to merge bins and the
[`xdev_merge_bins()`](../reference/xdev_merge_bins.md) function will
allow you to do so. For this example we will randomly select 100 ‘otu’
bins to merge.

``` r
random_bin_names <- sample(
  names(
    data = miseq,
    type = "bins",
    bin_type = "otu"
  ),
  size = 100, replace = FALSE
)

xdev_merge_bins(
  data = miseq,
  bin_names = random_bin_names,
  reason = "merge_bins_example",
  bin_type = "otu"
)
miseq
#> modified_miseq:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns   numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.00
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0   2847.35
#> 25%-tile:        1  375 252.0000      0 4.000000     0  28464.50
#> Median:          1  375 252.0000      0 4.000000     0  56928.00
#> 75%-tile:        1  375 253.0000      0 5.000000     0  85391.50
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 111008.65
#> Maximum:         1  375 256.0000      0 6.000000     0 113854.00
#> Mean:            1  375 252.4468      0 4.368533     0      0.00
#> Unique seqs:  2324 
#> Total seqs:   113854 
#> 
#> contigs_report :
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000   0.000000      0    248.0000
#> 2.5%-tile:     0.0016101600 252.0000   0.000000      0    251.0000
#> 25%-tile:      0.0028177700 252.0000   0.000000      0    251.0000
#> Median:        0.0062948698 252.0000   2.000000      0    251.0000
#> 75%-tile:      0.0261963997 253.0000   4.000000      0    251.0000
#> 97.5%-tile:    1.1412299871 253.0000  76.000000      0    252.0000
#> Maximum:       3.0126200000 270.0000 120.000000      0    256.0000
#> Mean:          0.0985028379 252.5124   7.536301      0    251.0762
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000      0.000000
#> 2.5%-tile:        248.0000      1.000000
#> 25%-tile:         249.0000      2.000000
#> Median:           249.0000      2.000000
#> 75%-tile:         249.0000      2.000000
#> 97.5%-tile:       251.0000      3.000000
#> Maximum:          255.0000     22.000000
#> Mean:             249.2138      1.862447
#> Unique seqs:  2324 
#> Total seqs:   113854 
#> 
#> scrap_summary:
#>        type                trash_code unique total
#> 1  sequence  F3D0_exclusive_sequences    101   109
#> 2       otu F3D0_exclusive_sequences,     14    14
#> 3       otu        merge_bins_example     99  6674
#> 4       asv F3D0_exclusive_sequences,    101   109
#> 5 phylotype F3D0_exclusive_sequences,      2     2
#> 
#> Sample   Total:
#> F3D0 6082 
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
#> Early    55525 
#> Late 58329 
#> 
#> Number of unique seqs: 2324 
#> Total number of seqs: 113854 
#> Total number of otus: 418 
#> Total number of asvs: 2324 
#> Total number of phylotypes: 61
```

You can see that the sample totals, total sequence and unique sequences
remain the same but the number of otu bins is reduced by 99.

## Removing Functions

Over the course of your analysis, you may want to remove sequences,
bins, samples or contaminants. rdataset has several functions to help
with that, namely
[`xdev_remove_sequences()`](../reference/xdev_remove_sequences.md),
[`xdev_remove_bins()`](../reference/xdev_remove_bins.md),
[`xdev_remove_samples()`](../reference/xdev_remove_samples.md), and
[`xdev_remove_lineages()`](../reference/xdev_remove_lineages.md).

### Removing Sequences

There are several reasons you may want to remove sequences including
removing chimeras, removing sequence without good overlap and removing
sequences with ambiguous bases. The
[`xdev_remove_sequences()`](../reference/xdev_remove_sequences.md)
function allows you to easily do that. The
[`miseq_sop_example()`](../reference/miseq_sop_example.md) has already
been screened for chimeras, overlap, sequence length and ambiguous
bases, so let’s randomly select 100 sequences to remove.

``` r
random_sequence_names <- sample(names(data = miseq),
  size = 100, replace = FALSE
)
xdev_remove_sequences(
  data = miseq,
  sequence_names = random_sequence_names,
  trash_tags = rep("remove_sequences_example", 100)
)
miseq
#> modified_miseq:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns  numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0      1.0
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0   2529.6
#> 25%-tile:        1  375 252.0000      0 4.000000     0  25287.0
#> Median:          1  375 252.0000      0 4.000000     0  50573.0
#> 75%-tile:        1  375 253.0000      0 5.000000     0  75859.0
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0  98616.4
#> Maximum:         1  375 256.0000      0 6.000000     0 101144.0
#> Mean:            1  375 252.4739      0 4.400399     0      0.0
#> Unique seqs:  2224 
#> Total seqs:   101144 
#> 
#> contigs_report :
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000   0.000000      0    248.0000
#> 2.5%-tile:     0.0016101600 252.0000   0.000000      0    251.0000
#> 25%-tile:      0.0025611699 252.0000   0.000000      0    251.0000
#> Median:        0.0077774799 253.0000   2.000000      0    251.0000
#> 75%-tile:      0.0276404992 253.0000   4.000000      0    251.0000
#> 97.5%-tile:    1.1412299871 253.0000  76.000000      0    252.0000
#> Maximum:       3.0126200000 270.0000 120.000000      0    256.0000
#> Mean:          0.1093153473 252.5478   8.078799      0    251.0886
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000      0.000000
#> 2.5%-tile:        248.0000      1.000000
#> 25%-tile:         249.0000      2.000000
#> Median:           249.0000      2.000000
#> 75%-tile:         249.0000      2.000000
#> 97.5%-tile:       251.0000      3.000000
#> Maximum:          255.0000     22.000000
#> Mean:             249.1741      1.914419
#> Unique seqs:  2224 
#> Total seqs:   101144 
#> 
#> scrap_summary:
#>        type                trash_code unique total
#> 1  sequence  F3D0_exclusive_sequences    101   109
#> 2  sequence  remove_sequences_example    100 12710
#> 3       otu F3D0_exclusive_sequences,     14    14
#> 4       otu        merge_bins_example     99  6674
#> 5       otu  remove_sequences_example      9   157
#> 6       asv F3D0_exclusive_sequences,    101   109
#> 7       asv  remove_sequences_example    100 12710
#> 8 phylotype F3D0_exclusive_sequences,      2     2
#> 9 phylotype  remove_sequences_example      1     1
#> 
#> Sample   Total:
#> F3D0 5576 
#> F3D1 4094 
#> F3D141   4114 
#> F3D142   2067 
#> F3D143   2174 
#> F3D144   3103 
#> F3D145   4860 
#> F3D146   3472 
#> F3D147   10785 
#> F3D148   8454 
#> F3D149   8937 
#> F3D150   3729 
#> F3D2 13925 
#> F3D3 4543 
#> F3D5 3093 
#> F3D6 5655 
#> F3D7 3522 
#> F3D8 3841 
#> F3D9 5200 
#> 
#> Treatment   Total:
#> Early    49449 
#> Late 51695 
#> 
#> Number of unique seqs: 2224 
#> Total number of seqs: 101144 
#> Total number of otus: 409 
#> Total number of asvs: 2224 
#> Total number of phylotypes: 60
```

Note, sequences can also be removed by setting their abundance to 0.

### Removing Bins

Similarly, lets randomly select and remove 10 phylotype bins with the
[`xdev_remove_bins()`](../reference/xdev_remove_bins.md) function. Note,
removing the bins also removes sequences. The removal of the sequences
from the data set also effects the bins in the ‘otu’ and ‘asv’ clusters.

``` r
random_bin_names <- sample(
  names(
    data = miseq,
    type = "bins", bin_type = "phylotype"
  ),
  size = 10, replace = FALSE
)
xdev_remove_bins(
  data = miseq,
  bin_names = random_bin_names,
  trash_tags = rep("remove_bins_example", 10),
  bin_type = "phylotype"
)
miseq
#> modified_miseq:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0     1.0
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0  1980.9
#> 25%-tile:        1  375 252.0000      0 4.000000     0 19800.0
#> Median:          1  375 252.0000      0 4.000000     0 39599.0
#> 75%-tile:        1  375 253.0000      0 4.000000     0 59398.0
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 77217.1
#> Maximum:         1  375 256.0000      0 6.000000     0 79196.0
#> Mean:            1  375 252.3263      0 4.272842     0     0.0
#> Unique seqs:  1373 
#> Total seqs:   79196 
#> 
#> contigs_report :
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000   0.000000      0     248.000
#> 2.5%-tile:     0.0014405299 252.0000   0.000000      0     251.000
#> 25%-tile:      0.0028177700 252.0000   0.000000      0     251.000
#> Median:        0.0077774799 252.0000   2.000000      0     251.000
#> 75%-tile:      0.0276404992 253.0000   7.000000      0     251.000
#> 97.5%-tile:    1.1412299871 253.0000  76.000000      0     252.000
#> Maximum:       3.0126200000 270.0000 120.000000      0     256.000
#> Mean:          0.1339895157 252.4207   9.832542      0     251.097
#>             Overlap_Length Overlap_Start
#> Minimum:           232.000      0.000000
#> 2.5%-tile:         249.000      1.000000
#> 25%-tile:          249.000      2.000000
#> Median:            249.000      2.000000
#> 75%-tile:          249.000      2.000000
#> 97.5%-tile:        251.000      2.000000
#> Maximum:           255.000     22.000000
#> Mean:              249.252      1.845043
#> Unique seqs:  1373 
#> Total seqs:   79196 
#> 
#> scrap_summary:
#>         type                trash_code unique total
#> 1   sequence  F3D0_exclusive_sequences    101   109
#> 2   sequence       remove_bins_example    851 21948
#> 3   sequence  remove_sequences_example    100 12710
#> 4        otu F3D0_exclusive_sequences,     14    14
#> 5        otu        merge_bins_example     99  6674
#> 6        otu       remove_bins_example    134 15766
#> 7        otu  remove_sequences_example      9   157
#> 8        asv F3D0_exclusive_sequences,    101   109
#> 9        asv       remove_bins_example    851 21948
#> 10       asv  remove_sequences_example    100 12710
#> 11 phylotype F3D0_exclusive_sequences,      2     2
#> 12 phylotype       remove_bins_example     10 23530
#> 13 phylotype  remove_sequences_example      1     1
#> 
#> Sample   Total:
#> F3D0 3984 
#> F3D1 2323 
#> F3D141   3284 
#> F3D142   1784 
#> F3D143   1741 
#> F3D144   2718 
#> F3D145   4410 
#> F3D146   2619 
#> F3D147   9428 
#> F3D148   7186 
#> F3D149   7076 
#> F3D150   2827 
#> F3D2 10106 
#> F3D3 4083 
#> F3D5 1986 
#> F3D6 4391 
#> F3D7 3072 
#> F3D8 2643 
#> F3D9 3535 
#> 
#> Treatment   Total:
#> Early    36123 
#> Late 43073 
#> 
#> Number of unique seqs: 1373 
#> Total number of seqs: 79196 
#> Total number of otus: 275 
#> Total number of asvs: 1373 
#> Total number of phylotypes: 50
```

Looking closer at the scrap summary we can see that removing 10
phylotype bins, removed 851 unique sequences that represented 21948
reads. After the removal of the 851 unique sequences the ‘otu’ cluster
removed 134 bins and the ‘asv’ cluster removed 851 bins.

### Removing Samples

If you included a mock community in your data set, you will want to
remove it after assessing your error rates in preparation for the rest
of your analysis. The
[`miseq_sop_example()`](../reference/miseq_sop_example.md) already has
the mock community removed so for the sake of example we will remove
sample ‘F3D142’.

``` r
xdev_remove_samples(
  data = miseq,
  samples = c("F3D142"), reason = "remove_samples_example"
)
miseq
#> modified_miseq:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns numseqs
#> Minimum:         1  375 249.0000      0  3.00000     0     1.0
#> 2.5%-tile:       1  375 252.0000      0  3.00000     0  1936.3
#> 25%-tile:        1  375 252.0000      0  4.00000     0 19354.0
#> Median:          1  375 252.0000      0  4.00000     0 38707.0
#> 75%-tile:        1  375 253.0000      0  4.00000     0 58060.0
#> 97.5%-tile:      1  375 253.0000      0  6.00000     0 75477.7
#> Maximum:         1  375 256.0000      0  6.00000     0 77412.0
#> Mean:            1  375 252.3263      0  4.27196     0     0.0
#> Unique seqs:  1335 
#> Total seqs:   77412 
#> 
#> contigs_report :
#>             Expected_Errors  Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.000   0.000000      0    248.0000
#> 2.5%-tile:     0.0014405299 252.000   0.000000      0    251.0000
#> 25%-tile:      0.0028177700 252.000   0.000000      0    251.0000
#> Median:        0.0077774799 252.000   2.000000      0    251.0000
#> 75%-tile:      0.0276404992 253.000   7.000000      0    251.0000
#> 97.5%-tile:    1.1412299871 253.000  76.000000      0    252.0000
#> Maximum:       3.0126200000 270.000 120.000000      0    256.0000
#> Mean:          0.1343094619 252.421   9.853434      0    251.0969
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000       0.00000
#> 2.5%-tile:        249.0000       1.00000
#> 25%-tile:         249.0000       2.00000
#> Median:           249.0000       2.00000
#> 75%-tile:         249.0000       2.00000
#> 97.5%-tile:       251.0000       2.00000
#> Maximum:          255.0000      22.00000
#> Mean:             249.2524       1.84452
#> Unique seqs:  1335 
#> Total seqs:   77412 
#> 
#> scrap_summary:
#>         type                trash_code unique total
#> 1   sequence  F3D0_exclusive_sequences    101   109
#> 2   sequence       remove_bins_example    851 21948
#> 3   sequence    remove_samples_example     38    42
#> 4   sequence  remove_sequences_example    100 12710
#> 5        otu F3D0_exclusive_sequences,     14    14
#> 6        otu        merge_bins_example     99  6674
#> 7        otu       remove_bins_example    134 15766
#> 8        otu    remove_samples_example      9    15
#> 9        otu  remove_sequences_example      9   157
#> 10       asv F3D0_exclusive_sequences,    101   109
#> 11       asv       remove_bins_example    851 21948
#> 12       asv    remove_samples_example     38    42
#> 13       asv  remove_sequences_example    100 12710
#> 14 phylotype F3D0_exclusive_sequences,      2     2
#> 15 phylotype       remove_bins_example     10 23530
#> 16 phylotype    remove_samples_example      1     1
#> 17 phylotype  remove_sequences_example      1     1
#> 
#> Sample   Total:
#> F3D0 3984 
#> F3D1 2323 
#> F3D141   3284 
#> F3D143   1741 
#> F3D144   2718 
#> F3D145   4410 
#> F3D146   2619 
#> F3D147   9428 
#> F3D148   7186 
#> F3D149   7076 
#> F3D150   2827 
#> F3D2 10106 
#> F3D3 4083 
#> F3D5 1986 
#> F3D6 4391 
#> F3D7 3072 
#> F3D8 2643 
#> F3D9 3535 
#> 
#> Treatment   Total:
#> Early    36123 
#> Late 41289 
#> 
#> Number of unique seqs: 1335 
#> Total number of seqs: 77412 
#> Total number of otus: 266 
#> Total number of asvs: 1335 
#> Total number of phylotypes: 49
```

Lastly, we can remove contaminants from the data set using the
[`xdev_remove_lineages()`](../reference/xdev_remove_lineages.md)
function. The miseq example has already had the contaminants removed
after classification and before bin assignment. For the sake of example
we will remove all sequences assigned to
‘Bacteria;Firmicutes;Clostridia;Clostridiales;Lachnospiraceae;’.

``` r
bad_tax <- paste0(
  "Bacteria;Firmicutes;Clostridia;Clostridiales;",
  "Clostridiales_unclassified;"
)

xdev_remove_lineages(
  data = miseq,
  contaminants = c(bad_tax),
  reason = "remove_contaminants_example"
)
miseq
#> modified_miseq:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns  numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0     1.00
#> 2.5%-tile:       1  375 252.0000      0 3.000000     0  1893.25
#> 25%-tile:        1  375 252.0000      0 4.000000     0 18923.50
#> Median:          1  375 252.0000      0 4.000000     0 37846.00
#> 75%-tile:        1  375 253.0000      0 4.000000     0 56768.50
#> 97.5%-tile:      1  375 253.0000      0 6.000000     0 73798.75
#> Maximum:         1  375 256.0000      0 6.000000     0 75690.00
#> Mean:            1  375 252.3107      0 4.258555     0     0.00
#> Unique seqs:  1221 
#> Total seqs:   75690 
#> 
#> contigs_report :
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000    0.00000      0    248.0000
#> 2.5%-tile:     0.0014405299 252.0000    0.00000      0    251.0000
#> 25%-tile:      0.0028177700 252.0000    0.00000      0    251.0000
#> Median:        0.0077774799 252.0000    2.00000      0    251.0000
#> 75%-tile:      0.0276404992 253.0000    7.00000      0    251.0000
#> 97.5%-tile:    1.1412299871 253.0000   76.00000      0    252.0000
#> Maximum:       3.0126200000 270.0000  120.00000      0    256.0000
#> Mean:          0.1372073502 252.4075   10.05614      0    251.1045
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000      0.000000
#> 2.5%-tile:        249.0000      1.000000
#> 25%-tile:         249.0000      2.000000
#> Median:           249.0000      2.000000
#> 75%-tile:         249.0000      2.000000
#> 97.5%-tile:       251.0000      2.000000
#> Maximum:          255.0000     22.000000
#> Mean:             249.2658      1.838711
#> Unique seqs:  1221 
#> Total seqs:   75690 
#> 
#> scrap_summary:
#>         type                  trash_code unique total
#> 1   sequence    F3D0_exclusive_sequences    101   109
#> 2   sequence         remove_bins_example    851 21948
#> 3   sequence remove_contaminants_example    114  1722
#> 4   sequence      remove_samples_example     38    42
#> 5   sequence    remove_sequences_example    100 12710
#> 6        otu   F3D0_exclusive_sequences,     14    14
#> 7        otu          merge_bins_example     99  6674
#> 8        otu         remove_bins_example    134 15766
#> 9        otu remove_contaminants_example     37  4428
#> 10       otu      remove_samples_example      9    15
#> 11       otu    remove_sequences_example      9   157
#> 12       asv   F3D0_exclusive_sequences,    101   109
#> 13       asv         remove_bins_example    851 21948
#> 14       asv remove_contaminants_example    114  1752
#> 15       asv      remove_samples_example     38    42
#> 16       asv    remove_sequences_example    100 12710
#> 17 phylotype   F3D0_exclusive_sequences,      2     2
#> 18 phylotype         remove_bins_example     10 23530
#> 19 phylotype remove_contaminants_example      1  1773
#> 20 phylotype      remove_samples_example      1     1
#> 21 phylotype    remove_sequences_example      1     1
#> 
#> Sample   Total:
#> F3D0 3859 
#> F3D1 2224 
#> F3D141   3191 
#> F3D143   1719 
#> F3D144   2683 
#> F3D145   4355 
#> F3D146   2551 
#> F3D147   9266 
#> F3D148   6977 
#> F3D149   6877 
#> F3D150   2725 
#> F3D2 9935 
#> F3D3 4018 
#> F3D5 1926 
#> F3D6 4331 
#> F3D7 3019 
#> F3D8 2586 
#> F3D9 3448 
#> 
#> Treatment   Total:
#> Early    35346 
#> Late 40344 
#> 
#> Number of unique seqs: 1221 
#> Total number of seqs: 75690 
#> Total number of otus: 229 
#> Total number of asvs: 1221 
#> Total number of phylotypes: 48
```

Thanks for following along. To explore more *xdev\_* functions you can
look at the rcpp_xint_xdev_functions.h and rcpp_xint_xdev_functions.cpp
files located in the src folder of the package.
