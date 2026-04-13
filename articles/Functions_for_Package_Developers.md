# Functions for package developers

*strollur* includes several functions designed for package developers.
These functions provide additional access to the back end data and allow
developers to merge, remove and set data directly. The function names
begin with *‘xdev\_’* or *‘xint\_’* to clearly indicate that they are
designed for development or for internal use, and not for the general
user.

### Under The Hood

The `dataset` object is an R6 object to keep the memory usage low. An R6
class is inherently passed by reference rather than by value. You can
make a deep copy of your dataset using the
[`copy_dataset()`](https://mothur.org/strollur/reference/copy_dataset.md)
function. Note, if you use an assignment operator to copy it’s a shallow
copy.

The `dataset` object has several public fields, but I’d like to bring
special attention to the *data* field. *data* is an Rcpp external
pointer (safe pointer) to ‘Dataset’ c++ class (class definitions found
in stroller.h). This format allows package developers an easy access
point to the underlying C++ code with additional functionality. While
many of the most frequently used functions have exported Rcpp functions,
you can also write your own functions to access any public function in
the ‘Dataset’ class.

### Setting functions

Over the course of an analysis your package may want to change the
abundances of sequences, modify sequence nucleotide strings due to
alignment, screening or filtering, and change the data set name. Let’s
take a close look at functions you will need to do so using the
[`miseq_sop_example()`](https://mothur.org/strollur/reference/miseq_sop_example.md)
data set.

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
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        3     0   2850.08
#> 25%-tile:        1  375    252      0        4     0  28491.75
#> Median:          1  375    252      0        4     0  56982.50
#> 75%-tile:        1  375    253      0        5     0  85473.25
#> 97.5%-tile:      1  375    253      0        6     0 111114.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0      0.00
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 531 
#> Total number of otu bin classifications: 531 
#> Total number of asvs: 2425 
#> Total number of asv bin classifications: 2425 
#> Total number of phylotypes: 63 
#> Total number of phylotype bin classifications: 63 
#> Total number of sequence classifications: 2425 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
```

#### Setting the dataset name

To change the data set name you can use the
[`xdev_set_dataset_name()`](https://mothur.org/strollur/reference/xdev_set_dataset_name.md)
function.

``` r
names(miseq, type = "dataset")
#> [1] "miseq_sop"
xdev_set_dataset_name(miseq, dataset_name = "modified_miseq")
names(miseq, type = "dataset")
#> [1] "modified_miseq"
```

#### Setting nucleotide sequence strings

The miseq example data set is already aligned, filtered and screened,
but for the sake of example let’s set the nucleotide strings of
sequences *exclusive* to sample ‘F3D0’ to ‘NNNN’. First we will get the
names of the sequences only present in sample ‘F3D0’.

``` r
f3d0_names <- names(
  miseq,
  type = "sequences",
  samples = c("F3D0"),
  distinct = TRUE
)
length(f3d0_names)
#> [1] 101
```

Now we will assign nucleotide strings of the sequences *exclusive* to
sample ‘F3D0’ to ‘NNNN’. Then we can use the
[`xdev_get_by_sample()`](https://mothur.org/strollur/reference/xdev_get_by_sample.md)
get all the sequences in present in sample ‘F3D0’ for closer inspection.

``` r
nnnn_neucleotides <- rep("NNNN", length(f3d0_names))
comments <- rep("example_set_sequences", length(f3d0_names))

xdev_set_sequences(
  miseq,
  sequence_names = f3d0_names,
  sequences = nnnn_neucleotides, comments = comments
)

f3d0_sequences <- xdev_get_by_sample(
  miseq,
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

#### Setting sequence abundances

Now that we know how to set sequence strings, let’s learn how to set
sequence abundances. In most cases you can use the
[`assign()`](https://mothur.org/strollur/reference/assign.md) function
to set sequence abundances, but there may be cases where you want more
flexibility. The assign function requires abundances to be provided for
each sequence. If you only need to update a subset of the sequences,
there are two functions to that allow you to do that,
[`xdev_set_abundances()`](https://mothur.org/strollur/reference/xdev_set_abundances.md)
and
[`xdev_set_abundance()`](https://mothur.org/strollur/reference/xdev_set_abundance.md).
*xdev_set_abundances* is used with data sets that include samples,
*xdev_set_abundance* is used for data sets without samples. Let’s set
the abundances of all sequences *exclusive* to sample F3D0 to 0, in
essence removing them. First we will get the abundance table for the
whole data set and create the inputs for
[`xdev_set_abundances()`](https://mothur.org/strollur/reference/xdev_set_abundances.md).

``` r
abundance_table <- abundance(
  miseq,
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

num_samples <- count(miseq, type = "samples")
new_abunds <- rep(list(rep(0, num_samples)), length(f3d0_names))

xdev_set_abundances(
  miseq,
  sequence_names = f3d0_names,
  abundances = new_abunds,
  reason = "F3D0_exclusive_sequences"
)
miseq
#> modified_miseq:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        3     0   2847.35
#> 25%-tile:        1  375    252      0        4     0  28464.50
#> Median:          1  375    252      0        4     0  56928.00
#> 75%-tile:        1  375    253      0        5     0  85391.50
#> 97.5%-tile:      1  375    253      0        6     0 111008.65
#> Maximum:         1  375    256      0        6     0 113854.00
#> Mean:            1  375    252      0        4     0      0.00
#> scrap_summary:
#>        type               trash_code unique total
#> 1  sequence F3D0_exclusive_sequences    101   109
#> 2       otu F3D0_exclusive_sequences     14    14
#> 3       asv F3D0_exclusive_sequences    101   109
#> 4 phylotype F3D0_exclusive_sequences      2     2
#> 
#> Number of unique seqs: 2324 
#> Total number of seqs: 113854 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 517 
#> Total number of otu bin classifications: 517 
#> Total number of asvs: 2324 
#> Total number of asv bin classifications: 2324 
#> Total number of phylotypes: 61 
#> Total number of phylotype bin classifications: 61 
#> Total number of sequence classifications: 2324 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
```

You can see that not all of the sequences in ‘F3D0’ are removed. That is
because the majority of sequences in the data set are present in
multiple samples. mothur2 uses xdev_set_abundances when running the
pre_cluster() function, sequences are clustered by sample, and the
abundances are set accordingly.

### Merging Functions

In the course of analyzing your data, there may be times when you want
to merge sequence abundances. When sequence nucleotide string are
identical you can save time processing by merging the identical read’s
abundances. mothur2 merges sequences abundances in the unique_seqs
function. The miseq example includes merged sequences that have been
assigned to bins. In practice, you would would not merge sequences after
assigning them to bins. Note, the code below is set to ‘eval=FALSE’
because you can only merge sequences assigned to the same bin and it
will throw errors accordingly. Let’s take a look at how to use the
[`xdev_merge_sequences()`](https://mothur.org/strollur/reference/xdev_merge_sequences.md)
function for your reference.

``` r
random_sequence_names <- sample(names(miseq),
  size = 100, replace = FALSE
)

xdev_merge_sequences(
  miseq,
  sequence_names = random_sequence_names,
  reason = "merge_sequences_example"
)
```

Similarly, you may want to merge bins and the
[`xdev_merge_bins()`](https://mothur.org/strollur/reference/xdev_merge_bins.md)
function will allow you to do so. For this example we will randomly
select 100 ‘otu’ bins to merge.

``` r
random_bin_names <- sample(
  names(
    miseq,
    type = "bins",
    bin_type = "otu"
  ),
  size = 100, replace = FALSE
)

xdev_merge_bins(
  miseq,
  bin_names = random_bin_names,
  reason = "merge_bins_example",
  bin_type = "otu"
)
miseq
#> modified_miseq:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        3     0   2847.35
#> 25%-tile:        1  375    252      0        4     0  28464.50
#> Median:          1  375    252      0        4     0  56928.00
#> 75%-tile:        1  375    253      0        5     0  85391.50
#> 97.5%-tile:      1  375    253      0        6     0 111008.65
#> Maximum:         1  375    256      0        6     0 113854.00
#> Mean:            1  375    252      0        4     0      0.00
#> scrap_summary:
#>        type               trash_code unique total
#> 1  sequence F3D0_exclusive_sequences    101   109
#> 2       otu F3D0_exclusive_sequences     14    14
#> 3       otu       merge_bins_example     99  6674
#> 4       asv F3D0_exclusive_sequences    101   109
#> 5 phylotype F3D0_exclusive_sequences      2     2
#> 
#> Number of unique seqs: 2324 
#> Total number of seqs: 113854 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 418 
#> Total number of otu bin classifications: 418 
#> Total number of asvs: 2324 
#> Total number of asv bin classifications: 2324 
#> Total number of phylotypes: 61 
#> Total number of phylotype bin classifications: 61 
#> Total number of sequence classifications: 2324 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
```

You can see that the sample totals, total sequence and unique sequences
remain the same but the number of otu bins is reduced by 99.

### Removing Functions

Over the course of your analysis, you may want to remove sequences,
bins, samples or contaminants. strollur has several functions to help
with that, namely
[`xdev_remove_sequences()`](https://mothur.org/strollur/reference/xdev_remove_sequences.md),
[`xdev_remove_bins()`](https://mothur.org/strollur/reference/xdev_remove_bins.md),
[`xdev_remove_samples()`](https://mothur.org/strollur/reference/xdev_remove_samples.md),
and
[`xdev_remove_lineages()`](https://mothur.org/strollur/reference/xdev_remove_lineages.md).

#### Removing Sequences

There are several reasons you may want to remove sequences including
removing chimeras, removing sequence without good overlap and removing
sequences with ambiguous bases. The
[`xdev_remove_sequences()`](https://mothur.org/strollur/reference/xdev_remove_sequences.md)
function allows you to easily do that. The
[`miseq_sop_example()`](https://mothur.org/strollur/reference/miseq_sop_example.md)
has already been screened for chimeras, overlap, sequence length and
ambiguous bases, so let’s randomly select 100 sequences to remove.

``` r
random_sequence_names <- sample(names(miseq),
  size = 100, replace = FALSE
)
xdev_remove_sequences(
  miseq,
  sequence_names = random_sequence_names,
  trash_tags = rep("remove_sequences_example", 100)
)
miseq
#> modified_miseq:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        3     0   2529.60
#> 25%-tile:        1  375    252      0        4     0  25287.00
#> Median:          1  375    252      0        4     0  50573.00
#> 75%-tile:        1  375    253      0        5     0  75859.00
#> 97.5%-tile:      1  375    253      0        6     0  98616.40
#> Maximum:         1  375    256      0        6     0 101144.00
#> Mean:            1  375    252      0        4     0      0.00
#> scrap_summary:
#>        type               trash_code unique total
#> 1  sequence F3D0_exclusive_sequences    101   109
#> 2  sequence remove_sequences_example    100 12710
#> 3       otu F3D0_exclusive_sequences     14    14
#> 4       otu       merge_bins_example     99  6674
#> 5       otu remove_sequences_example      9   157
#> 6       asv F3D0_exclusive_sequences    101   109
#> 7       asv remove_sequences_example    100 12710
#> 8 phylotype F3D0_exclusive_sequences      2     2
#> 9 phylotype remove_sequences_example      1     1
#> 
#> Number of unique seqs: 2224 
#> Total number of seqs: 101144 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 409 
#> Total number of otu bin classifications: 409 
#> Total number of asvs: 2224 
#> Total number of asv bin classifications: 2224 
#> Total number of phylotypes: 60 
#> Total number of phylotype bin classifications: 60 
#> Total number of sequence classifications: 2224 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
```

Note, sequences can also be removed by setting their abundance to 0.

#### Removing Bins

Similarly, lets randomly select and remove 10 phylotype bins with the
[`xdev_remove_bins()`](https://mothur.org/strollur/reference/xdev_remove_bins.md)
function. Note, removing the bins also removes sequences. The removal of
the sequences from the data set also effects the bins in the ‘otu’ and
‘asv’ clusters.

``` r
random_bin_names <- sample(
  names(
    miseq,
    type = "bins", bin_type = "phylotype"
  ),
  size = 10, replace = FALSE
)
xdev_remove_bins(
  miseq,
  bin_names = random_bin_names,
  trash_tags = rep("remove_bins_example", 10),
  bin_type = "phylotype"
)
miseq
#> modified_miseq:
#> 
#>             starts ends nbases ambigs polymers numns  numseqs
#> Minimum:         1  375    249      0        3     0     1.00
#> 2.5%-tile:       1  375    252      0        3     0  1980.90
#> 25%-tile:        1  375    252      0        4     0 19800.00
#> Median:          1  375    252      0        4     0 39599.00
#> 75%-tile:        1  375    253      0        4     0 59398.00
#> 97.5%-tile:      1  375    253      0        6     0 77217.10
#> Maximum:         1  375    256      0        6     0 79196.00
#> Mean:            1  375    252      0        4     0     0.00
#> scrap_summary:
#>         type               trash_code unique total
#> 1   sequence F3D0_exclusive_sequences    101   109
#> 2   sequence      remove_bins_example    851 21948
#> 3   sequence remove_sequences_example    100 12710
#> 4        otu F3D0_exclusive_sequences     14    14
#> 5        otu       merge_bins_example     99  6674
#> 6        otu      remove_bins_example    134 15766
#> 7        otu remove_sequences_example      9   157
#> 8        asv F3D0_exclusive_sequences    101   109
#> 9        asv      remove_bins_example    851 21948
#> 10       asv remove_sequences_example    100 12710
#> 11 phylotype F3D0_exclusive_sequences      2     2
#> 12 phylotype      remove_bins_example     10 23530
#> 13 phylotype remove_sequences_example      1     1
#> 
#> Number of unique seqs: 1373 
#> Total number of seqs: 79196 
#> 
#> Total number of samples: 19 
#> Total number of treatments: 2 
#> Total number of otus: 275 
#> Total number of otu bin classifications: 275 
#> Total number of asvs: 1373 
#> Total number of asv bin classifications: 1373 
#> Total number of phylotypes: 50 
#> Total number of phylotype bin classifications: 50 
#> Total number of sequence classifications: 1373 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
```

Looking closer at the scrap summary we can see that removing 10
phylotype bins, removed 851 unique sequences that represented 21948
reads. After the removal of the 851 unique sequences the ‘otu’ cluster
removed 134 bins and the ‘asv’ cluster removed 851 bins.

#### Removing Samples

If you included a mock community in your data set, you will want to
remove it after assessing your error rates in preparation for the rest
of your analysis. The
[`miseq_sop_example()`](https://mothur.org/strollur/reference/miseq_sop_example.md)
already has the mock community removed so for the sake of example we
will remove sample ‘F3D142’.

``` r
xdev_remove_samples(
  miseq,
  samples = c("F3D142"), reason = "remove_samples_example"
)
miseq
#> modified_miseq:
#> 
#>             starts ends nbases ambigs polymers numns  numseqs
#> Minimum:         1  375    249      0        3     0     1.00
#> 2.5%-tile:       1  375    252      0        3     0  1936.30
#> 25%-tile:        1  375    252      0        4     0 19354.00
#> Median:          1  375    252      0        4     0 38707.00
#> 75%-tile:        1  375    253      0        4     0 58060.00
#> 97.5%-tile:      1  375    253      0        6     0 75477.70
#> Maximum:         1  375    256      0        6     0 77412.00
#> Mean:            1  375    252      0        4     0     0.00
#> scrap_summary:
#>         type               trash_code unique total
#> 1   sequence F3D0_exclusive_sequences    101   109
#> 2   sequence      remove_bins_example    851 21948
#> 3   sequence   remove_samples_example     38    42
#> 4   sequence remove_sequences_example    100 12710
#> 5        otu F3D0_exclusive_sequences     14    14
#> 6        otu       merge_bins_example     99  6674
#> 7        otu      remove_bins_example    134 15766
#> 8        otu   remove_samples_example      9    15
#> 9        otu remove_sequences_example      9   157
#> 10       asv F3D0_exclusive_sequences    101   109
#> 11       asv      remove_bins_example    851 21948
#> 12       asv   remove_samples_example     38    42
#> 13       asv remove_sequences_example    100 12710
#> 14 phylotype F3D0_exclusive_sequences      2     2
#> 15 phylotype      remove_bins_example     10 23530
#> 16 phylotype   remove_samples_example      1     1
#> 17 phylotype remove_sequences_example      1     1
#> 
#> Number of unique seqs: 1335 
#> Total number of seqs: 77412 
#> 
#> Total number of samples: 18 
#> Total number of treatments: 2 
#> Total number of otus: 266 
#> Total number of otu bin classifications: 266 
#> Total number of asvs: 1335 
#> Total number of asv bin classifications: 1335 
#> Total number of phylotypes: 49 
#> Total number of phylotype bin classifications: 49 
#> Total number of sequence classifications: 1335 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
```

Lastly, we can remove contaminants from the data set using the
[`xdev_remove_lineages()`](https://mothur.org/strollur/reference/xdev_remove_lineages.md)
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
  miseq,
  contaminants = c(bad_tax),
  reason = "remove_contaminants_example"
)
miseq
#> modified_miseq:
#> 
#>             starts ends nbases ambigs polymers numns  numseqs
#> Minimum:         1  375    249      0        3     0     1.00
#> 2.5%-tile:       1  375    252      0        3     0  1893.25
#> 25%-tile:        1  375    252      0        4     0 18923.50
#> Median:          1  375    252      0        4     0 37846.00
#> 75%-tile:        1  375    253      0        4     0 56768.50
#> 97.5%-tile:      1  375    253      0        6     0 73798.75
#> Maximum:         1  375    256      0        6     0 75690.00
#> Mean:            1  375    252      0        4     0     0.00
#> scrap_summary:
#>         type                  trash_code unique total
#> 1   sequence    F3D0_exclusive_sequences    101   109
#> 2   sequence         remove_bins_example    851 21948
#> 3   sequence remove_contaminants_example    114  1722
#> 4   sequence      remove_samples_example     38    42
#> 5   sequence    remove_sequences_example    100 12710
#> 6        otu    F3D0_exclusive_sequences     14    14
#> 7        otu          merge_bins_example     99  6674
#> 8        otu         remove_bins_example    134 15766
#> 9        otu remove_contaminants_example     37  4428
#> 10       otu      remove_samples_example      9    15
#> 11       otu    remove_sequences_example      9   157
#> 12       asv    F3D0_exclusive_sequences    101   109
#> 13       asv         remove_bins_example    851 21948
#> 14       asv remove_contaminants_example    114  1752
#> 15       asv      remove_samples_example     38    42
#> 16       asv    remove_sequences_example    100 12710
#> 17 phylotype    F3D0_exclusive_sequences      2     2
#> 18 phylotype         remove_bins_example     10 23530
#> 19 phylotype remove_contaminants_example      1  1773
#> 20 phylotype      remove_samples_example      1     1
#> 21 phylotype    remove_sequences_example      1     1
#> 
#> Number of unique seqs: 1221 
#> Total number of seqs: 75690 
#> 
#> Total number of samples: 18 
#> Total number of treatments: 2 
#> Total number of otus: 229 
#> Total number of otu bin classifications: 229 
#> Total number of asvs: 1221 
#> Total number of asv bin classifications: 1221 
#> Total number of phylotypes: 48 
#> Total number of phylotype bin classifications: 48 
#> Total number of sequence classifications: 1221 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
```

Thanks for following along. To explore more *xdev\_* functions you can
look at the rcpp_xint_xdev_functions.h and rcpp_xint_xdev_functions.cpp
files located in the src folder of the package.
