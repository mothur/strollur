# Importing from mothur

*rdataset* includes the function
[`read_mothur()`](../reference/read_mothur.md) as well as several
functions to read [mothur](https://mothur.org) output files
individually. To create a data set from the outputs of the [Miseq SOP
Example](https://mothur.org/wiki/miseq_sop/), run the following:

## Using `read_mothur()`

``` r
data <- read_mothur(
  fasta = rdataset_example("final.fasta"),
  count = rdataset_example("final.count_table"),
  taxonomy = rdataset_example("final.taxonomy"),
  design = rdataset_example("mouse.time.design"),
  otu_list = rdataset_example("final.opti_mcc.list"),
  asv_list = rdataset_example("final.asv.list"),
  phylo_list = rdataset_example("final.tx.list"),
  sample_tree = rdataset_example("final.opti_mcc.jclass.ave.tre"),
  dataset_name = "miseq_sop"
)
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 19 samples to treatments.
```

To view a summary of data:

``` r
data
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

## Importing Individual Files

- [`read_fasta()`](../reference/read_fasta.md) read a
  [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) formatted
  sequence file
- [`read_mothur_count()`](../reference/read_mothur_count.md) read a
  mothur formatted [count file](https://mothur.org/wiki/count_file/)
- [`read_mothur_taxonomy()`](../reference/read_mothur_taxonomy.md) read
  a mothur formatted [taxonomy
  file](https://mothur.org/wiki/taxonomy_file/)
- [`read_mothur_cons_taxonomy()`](../reference/read_mothur_cons_taxonomy.md)
  read a mothur formatted [cons_taxonomy
  file](https://mothur.org/wiki/constaxonomy_file/)
- [`read_mothur_list()`](../reference/read_mothur_list.md) read a mothur
  formatted [list file](https://mothur.org/wiki/list_file/)
- [`read_mothur_shared()`](../reference/read_mothur_shared.md) read a
  mothur formatted [shared file](https://mothur.org/wiki/shared_file/)
- [`read_mothur_rabund()`](../reference/read_mothur_rabund.md) read a
  mothur formatted [rabund file](https://mothur.org/wiki/rabund_file/)

To create a data set and read the individual file types, you can use the
functions below. First let’s create a data set named my_data.

``` r
my_data <- new_dataset(dataset_name = "my_data")
```

To add [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) data
to your data set you can use the
[`read_fasta()`](../reference/read_fasta.md) function:

``` r
fasta_data <- read_fasta(fasta = rdataset_example("final.fasta"))
```

fasta_data is a data.frame containing sequence names, sequence
nucleotide strings, and comments if provided. You can add the FASTA
sequences to your data set using the [`add()`](../reference/add.md)
function:

``` r
add(data = my_data, table = fasta_data, type = "sequences")
#> ℹ Added 2425 sequences.
#> [1] 2425
my_data
#> my_data:
#> 
#> sequence_summary:
#>             starts ends   nbases ambigs polymers numns  numseqs
#> Minimum:         1  375 249.0000      0 3.000000     0    1.000
#> 2.5%-tile:       1  375 252.0000      0 4.000000     0   61.625
#> 25%-tile:        1  375 252.0000      0 4.000000     0  607.250
#> Median:          1  375 253.0000      0 4.000000     0 1213.500
#> 75%-tile:        1  375 253.0000      0 5.000000     0 1819.750
#> 97.5%-tile:      1  375 254.0000      0 6.000000     0 2365.375
#> Maximum:         1  375 256.0000      0 6.000000     0 2425.000
#> Mean:            1  375 252.7406      0 4.496082     0    0.000
#> Unique seqs:  2425 
#> Total seqs:   2425
#> → Your dataset does not include sample data, ignoring.
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425
```

To add your sequence abundance data, you can read a [mothur count
file](https://mothur.org/wiki/count_file/) file using the
[`read_mothur_count()`](../reference/read_mothur_count.md) function:

``` r
sample_table <- read_mothur_count(
  filename = rdataset_example("final.count_table")
)
```

sample_table is a data.frame containing sequence_names, samples, and
abundances. You can add the sequence abundance data to your data set
using the [`assign()`](../reference/assign.md) function:

``` r
assign(data = my_data, table = sample_table, type = "sequence_abundance")
#> ℹ Assigned 2425 sequence abundances.
#> [1] 2425
my_data
#> my_data:
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
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963
```

To add sequence taxonomy assignments, you can read a [taxonomy
file](https://mothur.org/wiki/taxonomy_file/) file using the
[`read_mothur_taxonomy()`](../reference/read_mothur_taxonomy.md)
function:

``` r
classification_data <- read_mothur_taxonomy(
  taxonomy = rdataset_example("final.taxonomy")
)
```

classification_data is a data.frame containing sequence names and
taxonomies. You can add the sequence classification data to your data
set as follows:

``` r
assign(data = my_data, table = classification_data, type = "sequence_taxonomy")
#> ℹ Assigned 2425 sequence taxonomies.
#> [1] 2425
```

To assign sequences to bins, you can read a [mothur list
file](https://mothur.org/wiki/list_file/) file using the
[`read_mothur_list()`](../reference/read_mothur_list.md) function:

``` r
otu_data <- read_mothur_list(list = rdataset_example("final.opti_mcc.list"))
asv_data <- read_mothur_list(list = rdataset_example("final.asv.list"))
phylotype_data <- read_mothur_list(list = rdataset_example("final.tx.list"))
```

otu_data, asv_data and phylotype_data are data.frames containing bin
names and sequence names. You can add the bin data to your data set as
follows:

``` r
assign(data = my_data, table = otu_data, type = "bins", bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531
assign(data = my_data, table = asv_data, type = "bins", bin_type = "asv")
#> ℹ Assigned 2425 asv bins.
#> [1] 2425
assign(
  data = my_data, table = phylotype_data,
  type = "bins", bin_type = "phylotype"
)
#> ℹ Assigned 63 phylotype bins.
#> [1] 63
my_data
#> my_data:
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
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> Total number of otus: 531 
#> Total number of asvs: 2425 
#> Total number of phylotypes: 63
```

When you assign bins to sequences with taxonomic assignments the data
set object will find the consensus taxonomy of the bins automatically.
If you wish to assign the bin taxonomy separately, you can read a
[mothur cons_taxonomy file](https://mothur.org/wiki/constaxonomy_file/)
file using the
[`read_mothur_cons_taxonomy()`](../reference/read_mothur_cons_taxonomy.md)
function:

``` r
otu_taxonomy_data <- read_mothur_cons_taxonomy(
  taxonomy =
    rdataset_example("final.cons.taxonomy")
)
```

otu_taxonomy_data is a data.frame containing bin names, abundances and
taxonomies. You can add the bin taxonomic data to your data set as
follows:

``` r
assign(data = my_data, table = otu_taxonomy_data, type = "bin_taxonomy")
#> ℹ Assigned 531 otu bin taxonomies.
#> [1] 531
```

## Writing mothur formatted file types

- [`write_mothur()`](../reference/write_mothur.md) write mothur
  formatted files for all data
- [`write_fasta()`](../reference/write_fasta.md) read a
  [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) formatted
  sequence file
- [`write_mothur_count()`](../reference/write_mothur_count.md) write a
  mothur formatted [count file](https://mothur.org/wiki/count_file/)
- [`write_mothur_design()`](../reference/write_mothur_design.md) write a
  mothur formatted [design file](https://mothur.org/wiki/design_file/)
- [`write_taxonomy()`](../reference/write_taxonomy.md) write a mothur
  formatted [taxonomy file](https://mothur.org/wiki/taxonomy_file/)
- [`write_mothur_cons_taxonomy()`](../reference/write_mothur_cons_taxonomy.md)
  write a mothur formatted [cons_taxonomy
  file](https://mothur.org/wiki/constaxonomy_file/)
- [`write_mothur_list()`](../reference/write_mothur_list.md) write a
  mothur formatted [list file](https://mothur.org/wiki/list_file/)
- [`write_mothur_shared()`](../reference/write_mothur_shared.md) write a
  mothur formatted [shared file](https://mothur.org/wiki/shared_file/)
- [`write_mothur_rabund()`](../reference/write_mothur_rabund.md) write a
  mothur formatted [rabund file](https://mothur.org/wiki/rabund_file/)
