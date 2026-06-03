# Data Transfers

The *strollur* package stores the data associated with your Amplicon
Sequence Analysis. This tutorial will explain how to save, load, copy,
export, and import your `strollur` object. If you haven’t reviewed the
[Getting
Started](https://mothur.org/strollur/articles/Getting-Started.md)
tuturial, we recommend you start there.

Let’s use the
[`miseq_sop_example()`](https://mothur.org/strollur/reference/miseq_sop_example.md)
function to create a strollur object from the [Miseq SOP
Example](https://mothur.org/wiki/miseq_sop/).

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
miseq
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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

## Saving and Loading

The strollur package has a function to save a dataset object as an
*.rds* file,
[`save_dataset()`](https://mothur.org/strollur/reference/save_dataset.md),
and a function to create a dataset from an *.rds* file,
[`load_dataset()`](https://mothur.org/strollur/reference/load_dataset.md).
Let’s use the miseq data object to learn how to do that.

``` r

file_name <- file.path(tempdir(), "miseq_sop.rds")
save_dataset(miseq, file = file_name)
#> [1] "/tmp/RtmpQ0v2SK/miseq_sop.rds"

miseq_from_rds <- load_dataset(file = file_name)
miseq_from_rds
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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
unlink(file_name)
```

We can see that the summaries of miseq and miseq_from_rds are identical.
Let’s modify miseq_from_rds to verify they are not referring to the same
object. We will add clusters created by [mothur](https://mothur.org)
using [vsearch’s](https://github.com/torognes/vsearch) distance-based
greedy clustering (dgc) algorithm.

``` r

dgc_data <- read_mothur_list(list = strollur_example("final.dgc.list.gz"))

assign(miseq_from_rds, table = dgc_data, bin_type = "dgc")
#> Assigned 361 dgc bins.
miseq_from_rds
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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
#> Total number of dgcs: 361 
#> Total number of dgc bin classifications: 361 
#> Total number of sequence classifications: 2425 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
miseq
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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

We can see from the summary that 361 ‘dgc’ bins were added to
miseq_from_rds and not to miseq.

## Export and Import

The *.rds* file is in binary format and is not human readable. You can
use the
[`export_dataset()`](https://mothur.org/strollur/reference/export_dataset.md)
to see a human readable form of the raw data stored in the dataset.
Let’s export *miseq* and look at the table created.

``` r

table <- export_dataset(miseq)
str(table)
#> List of 15
#>  $ sequence_data                    :'data.frame':   2425 obs. of  5 variables:
#>   ..$ sequence_id     : int [1:2425] 0 1 2 3 4 5 6 7 8 9 ...
#>   ..$ sequence_name   : chr [1:2425] "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10331_23332" "M00967_43_000000000-A3JHG_1_1101_10382_22128" "M00967_43_000000000-A3JHG_1_1101_11035_15765" ...
#>   ..$ sequence        : chr [1:2425] "TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GGC-G-G-C-CA-T-G-C-AA-G-T-"| __truncated__ "TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-T-T--GG-GT--TT-A-AA-GG-GA-GC-G-CA-GGC-G-G-C-AT-G-G-C-AA-G-T-"| __truncated__ "TAC--GT-AG-GTA--GCA-A-G-C-G-T-T--GT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GC-GT-G-TA-GCC-G-G-G-CT-T-A-C-AA-G-T-"| __truncated__ "TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GG-GC-G-TA-GAC-G-G-C-AG-T-G-C-AA-G-T-"| __truncated__ ...
#>   ..$ taxonomy        : chr [1:2425] "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(85);Lachnospiraceae_unclassified(85);" "Bacteria(100);Firmicutes(99);Clostridia(98);Clostridiales(98);Lachnospiraceae(95);Lachnospiraceae_unclassified(95);" "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(98);Ruminococcaceae_unclassified(98);" "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);Lachnospiraceae_unclassified(100);" ...
#>   ..$ include_sequence: logi [1:2425] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ sequence_report                  :'data.frame':   2425 obs. of  7 variables:
#>   ..$ sequence_id        : int [1:2425] 0 1 2 3 4 5 6 7 8 9 ...
#>   ..$ start              : int [1:2425] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ end                : int [1:2425] 375 375 375 375 375 375 375 375 375 375 ...
#>   ..$ length             : int [1:2425] 253 253 253 252 253 252 253 253 252 252 ...
#>   ..$ ambig              : int [1:2425] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ longest_homopolymer: int [1:2425] 5 4 5 5 4 5 5 5 4 4 ...
#>   ..$ num_n              : int [1:2425] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ sequence_abundance_table         :'data.frame':   5539 obs. of  4 variables:
#>   ..$ sequence_id: int [1:5539] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ abundance  : num [1:5539] 32 127 1 1 1 222 5 13 20 17 ...
#>   ..$ sample     : chr [1:5539] "F3D0" "F3D1" "F3D146" "F3D149" ...
#>   ..$ treatment  : chr [1:5539] "Early" "Early" "Late" "Late" ...
#>  $ otu_bin_data                     :'data.frame':   531 obs. of  5 variables:
#>   ..$ bin_id     : int [1:531] 0 1 2 3 4 5 6 7 8 9 ...
#>   ..$ bin_name   : chr [1:531] "Otu001" "Otu002" "Otu003" "Otu004" ...
#>   ..$ abundance  : num [1:531] 12288 8892 7794 7476 7450 ...
#>   ..$ taxonomy   : chr [1:531] "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);\""| __truncated__ "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);\""| __truncated__ "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);\""| __truncated__ "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);Barnesiella(100);" ...
#>   ..$ include_bin: logi [1:531] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ otu_sequence_bin_assignment      :'data.frame':   2425 obs. of  2 variables:
#>   ..$ bin_id     : int [1:2425] 28 274 34 32 41 3 123 333 220 0 ...
#>   ..$ sequence_id: int [1:2425] 0 1 2 3 4 5 6 7 8 9 ...
#>  $ otu_bin_representative_sequence  :'data.frame':   531 obs. of  2 variables:
#>   ..$ bin_id     : int [1:531] 0 1 2 3 4 5 6 7 8 9 ...
#>   ..$ sequence_id: int [1:531] 592 439 21 373 1501 419 666 859 2027 87 ...
#>  $ asv_bin_data                     :'data.frame':   2425 obs. of  5 variables:
#>   ..$ bin_id     : int [1:2425] 0 1 2 3 4 5 6 7 8 9 ...
#>   ..$ bin_name   : chr [1:2425] "Asv0001" "Asv0002" "Asv0003" "Asv0004" ...
#>   ..$ abundance  : num [1:2425] 12196 8829 7698 7436 7307 ...
#>   ..$ taxonomy   : chr [1:2425] "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);\""| __truncated__ "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);\""| __truncated__ "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);\""| __truncated__ "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);Barnesiella(100);" ...
#>   ..$ include_bin: logi [1:2425] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ asv_sequence_bin_assignment      :'data.frame':   2425 obs. of  2 variables:
#>   ..$ bin_id     : int [1:2425] 27 2419 2404 2401 518 2382 2377 2372 2363 283 ...
#>   ..$ sequence_id: int [1:2425] 0 1 2 3 4 5 6 7 8 9 ...
#>  $ phylotype_bin_data               :'data.frame':   63 obs. of  5 variables:
#>   ..$ bin_id     : int [1:63] 0 1 2 3 4 5 6 7 8 9 ...
#>   ..$ bin_name   : chr [1:63] "Phylo01" "Phylo02" "Phylo03" "Phylo04" ...
#>   ..$ abundance  : num [1:63] 21639 53147 2805 1773 5337 ...
#>   ..$ taxonomy   : chr [1:63] "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Lachnospiraceae(100);Lachnospiraceae_unclassified(100);" "Bacteria(100);\"Bacteroidetes\"(100);\"Bacteroidia\"(100);\"Bacteroidales\"(100);\"Porphyromonadaceae\"(100);\""| __truncated__ "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Ruminococcaceae(100);Ruminococcaceae_unclassified(100);" "Bacteria(100);Firmicutes(100);Clostridia(100);Clostridiales(100);Clostridiales_unclassified(100);Clostridiales_"| __truncated__ ...
#>   ..$ include_bin: logi [1:63] TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ phylotype_sequence_bin_assignment:'data.frame':   2425 obs. of  2 variables:
#>   ..$ bin_id     : int [1:2425] 0 0 2 0 0 1 3 0 24 1 ...
#>   ..$ sequence_id: int [1:2425] 0 1 2 3 4 5 6 7 8 9 ...
#>  $ resource_reference               :'data.frame':   2 obs. of  10 variables:
#>   ..$ vendor           : chr [1:2] "Schloss Lab - University of Michigan" "SILVA"
#>   ..$ name             : chr [1:2] "R phylotypr package" "silva.bacteria.fasta"
#>   ..$ version          : chr [1:2] "0.1.1" "1.38.1"
#>   ..$ usage            : chr [1:2] "classification of sequences" "alignment of sequences"
#>   ..$ note             : chr [1:2] "classification using Bayesian method" "alignment reference trimmed to V4 region"
#>   ..$ method_url       : chr [1:2] "doi:10.1128/mra.01144-24" "NA"
#>   ..$ documentation_url: chr [1:2] "https://mothur.org/phylotypr/" "https://mothur.org/wiki/silva_reference_files/"
#>   ..$ parameter        : chr [1:2] "kmer_size=8,num_bootstraps=100,min_confidence=80" "NA"
#>   ..$ citation         : chr [1:2] "@article{doi:10.1128/AEM.00062-07, author = {Qiong Wang and George M. Garrity and James M. Tiedje and James R. "| __truncated__ "NA"
#>   ..$ creation_date    : chr [1:2] "2026-06-03" "2026-06-03"
#>  $ metadata                         :'data.frame':   19 obs. of  2 variables:
#>   ..$ sample        : chr [1:19] "F3D0" "F3D1" "F3D141" "F3D142" ...
#>   ..$ days_post_wean: num [1:19] 0 1 141 142 143 144 145 146 147 148 ...
#>  $ contigs_report                   :'data.frame':   2425 obs. of  8 variables:
#>   ..$ Name           : chr [1:2425] "M00967_43_000000000-A3JHG_1_1101_10133_8460" "M00967_43_000000000-A3JHG_1_1101_10331_23332" "M00967_43_000000000-A3JHG_1_1101_10382_22128" "M00967_43_000000000-A3JHG_1_1101_11035_15765" ...
#>   ..$ Length         : num [1:2425] 253 253 253 252 253 252 253 253 252 252 ...
#>   ..$ Overlap_Length : num [1:2425] 249 249 249 250 249 249 248 249 250 248 ...
#>   ..$ Overlap_Start  : num [1:2425] 2 2 2 1 2 2 3 2 1 3 ...
#>   ..$ Overlap_End    : num [1:2425] 251 251 251 251 251 251 251 251 251 251 ...
#>   ..$ MisMatches     : num [1:2425] 0 0 2 19 2 0 4 1 2 8 ...
#>   ..$ Num_Ns         : num [1:2425] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ Expected_Errors: num [1:2425] 0.00207 0.00231 0.01135 0.08721 0.00674 ...
#>   ..- attr(*, "sequence_name")= chr "Name"
#>  $ sequence_tree                    :List of 4
#>   ..$ edge       : int [1:4848, 1:2] 2426 2427 2427 2426 2428 2429 2429 2428 2430 2431 ...
#>   ..$ edge.length: num [1:4848] NaN 0.00395 0.00395 0 0.00198 ...
#>   ..$ Nnode      : int 2424
#>   ..$ tip.label  : chr [1:2425] "M00967_43_000000000-A3JHG_1_1114_15727_25995" "M00967_43_000000000-A3JHG_1_2109_19976_22044" "M00967_43_000000000-A3JHG_1_1102_9244_9305" "M00967_43_000000000-A3JHG_1_2101_14159_9619" ...
#>   ..- attr(*, "class")= chr "phylo"
#>   ..- attr(*, "order")= chr "cladewise"
#>  $ sample_tree                      :List of 5
#>   ..$ edge       : int [1:36, 1:2] 20 21 22 23 24 25 26 27 27 26 ...
#>   ..$ edge.length: num [1:36] 0.03913 0.01741 0.02565 0.00379 0.01732 ...
#>   ..$ Nnode      : int 18
#>   ..$ tip.label  : chr [1:19] "F3D9" "F3D8" "F3D6" "F3D5" ...
#>   ..$ root.edge  : num 0.221
#>   ..- attr(*, "class")= chr "phylo"
#>   ..- attr(*, "order")= chr "cladewise"
#>  - attr(*, "strollur_version")= chr "0.1.0"
#>  - attr(*, "dataset_name")= chr "miseq_sop"
```

Similarly to
[`load_dataset()`](https://mothur.org/strollur/reference/load_dataset.md),
you can use the
[`import_dataset()`](https://mothur.org/strollur/reference/import_dataset.md)
function to create a new dataset object from the exported table.

``` r

miseq_import <- import_dataset(table = table)
#> Added 2425 sequences.
#> Assigned 2425 sequence taxonomies.
#> Assigned 2425 sequence abundances.
#> Assigned 531 otu bins.
#> Assigned 531 otu bin representative sequences.
#> Assigned 531 otu bin taxonomies.
#> Assigned 2425 asv bins.
#> Assigned 2425 asv bin taxonomies.
#> Assigned 63 phylotype bins.
#> Assigned 63 phylotype bin taxonomies.
#> Added metadata.
#> Added 2 resource references.
#> Added a contigs_report.
miseq_import
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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

Again, we can see that the summary of miseq_import is identical to the
summary of miseq.

## Copy

Lastly, you can make a deep copy of your dataset using the
[`copy_dataset()`](https://mothur.org/strollur/reference/copy_dataset.md)
function. Note, if you use an assignment operator to copy it’s a shallow
copy. The dataset object is an R6 object to keep the memory usage low.
First let’s learn how to use the
[`copy_dataset()`](https://mothur.org/strollur/reference/copy_dataset.md)
function, then we will take a closer look at how deep and shallow
copying differ.

``` r

miseq_deep_copy <- copy_dataset(miseq)

miseq_shallow_copy <- miseq
```

Let’s add the dgc_data to miseq_shallow_copy and then compare miseq,
miseq_deep_copy, and mise_shallow_copy.

``` r

assign(miseq_shallow_copy, table = dgc_data, bin_type = "dgc")
#> Assigned 361 dgc bins.

miseq
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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
#> Total number of dgcs: 361 
#> Total number of dgc bin classifications: 361 
#> Total number of sequence classifications: 2425 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata

miseq_shallow_copy
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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
#> Total number of dgcs: 361 
#> Total number of dgc bin classifications: 361 
#> Total number of sequence classifications: 2425 
#> Total number of resource references: 2 
#> Total number of custom reports: 1 
#> Your dataset includes metadata

miseq_deep_copy
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns   numseqs
#> Minimum:         1  375    249      0        3     0      1.00
#> 2.5%-tile:       1  375    252      0        4     0   2849.08
#> 25%-tile:        1  375    252      0        4     0  28490.75
#> Median:          1  375    253      0        4     0  56981.50
#> 75%-tile:        1  375    253      0        5     0  85472.25
#> 97.5%-tile:      1  375    254      0        6     0 111113.93
#> Maximum:         1  375    256      0        6     0 113963.00
#> Mean:            1  375    252      0        4     0  56981.64
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

You can see from the summaries that the dgc_data was added to both miseq
and miseq_shallow_copy because they actually reference the same object,
but miseq_deep_copy was not modified.
