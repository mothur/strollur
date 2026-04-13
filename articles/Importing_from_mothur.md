# Importing from mothur

*strollur* includes the function
[`read_mothur()`](https://mothur.org/strollur/reference/read_mothur.md)
as well as several functions to read [mothur](https://mothur.org) output
files individually. To create a data set from the outputs of the [Miseq
SOP Example](https://mothur.org/wiki/miseq_sop/), run the following:

## Using `read_mothur()`

``` r
data <- read_mothur(
  fasta = strollur_example("final.fasta.gz"),
  count = strollur_example("final.count_table.gz"),
  taxonomy = strollur_example("final.taxonomy.gz"),
  design = strollur_example("mouse.time.design"),
  otu_list = strollur_example("final.opti_mcc.list.gz"),
  asv_list = strollur_example("final.asv.list.gz"),
  phylo_list = strollur_example("final.tx.list.gz"),
  sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
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
```

## Importing Individual Files

- [`read_fasta()`](https://mothur.org/strollur/reference/read_fasta.md)
  read a [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
  formatted sequence file
- [`read_mothur_count()`](https://mothur.org/strollur/reference/read_mothur_count.md)
  read a mothur formatted [count
  file](https://mothur.org/wiki/count_file/)
- [`read_mothur_taxonomy()`](https://mothur.org/strollur/reference/read_mothur_taxonomy.md)
  read a mothur formatted [taxonomy
  file](https://mothur.org/wiki/taxonomy_file/)
- [`read_mothur_cons_taxonomy()`](https://mothur.org/strollur/reference/read_mothur_cons_taxonomy.md)
  read a mothur formatted [cons_taxonomy
  file](https://mothur.org/wiki/constaxonomy_file/)
- [`read_mothur_list()`](https://mothur.org/strollur/reference/read_mothur_list.md)
  read a mothur formatted [list
  file](https://mothur.org/wiki/list_file/)
- [`read_mothur_shared()`](https://mothur.org/strollur/reference/read_mothur_shared.md)
  read a mothur formatted [shared
  file](https://mothur.org/wiki/shared_file/)
- [`read_mothur_rabund()`](https://mothur.org/strollur/reference/read_mothur_rabund.md)
  read a mothur formatted [rabund
  file](https://mothur.org/wiki/rabund_file/)

To create a data set and read the individual file types, you can use the
functions below. First let’s create a data set named my_data.

``` r
my_data <- new_dataset(dataset_name = "my_data")
```

To add [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) data
to your data set you can use the
[`read_fasta()`](https://mothur.org/strollur/reference/read_fasta.md)
function:

``` r
fasta_data <- read_fasta(fasta = strollur_example("final.fasta.gz"))
```

fasta_data is a data.frame containing sequence names, sequence
nucleotide strings, and comments if provided. You can add the FASTA
sequences to your data set using the
[`add()`](https://mothur.org/strollur/reference/add.md) function:

``` r
add(my_data, table = fasta_data, type = "sequences")
#> ℹ Added 2425 sequences.
#> [1] 2425
my_data
#> my_data:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  375    249      0        3     0    1.00
#> 2.5%-tile:       1  375    252      0        4     0   61.62
#> 25%-tile:        1  375    252      0        4     0  607.25
#> Median:          1  375    253      0        4     0 1213.50
#> 75%-tile:        1  375    253      0        5     0 1819.75
#> 97.5%-tile:      1  375    254      0        6     0 2365.38
#> Maximum:         1  375    256      0        6     0 2425.00
#> Mean:            1  375    252      0        4     0    0.00
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425
```

To add your sequence abundance data, you can read a [mothur count
file](https://mothur.org/wiki/count_file/) file using the
[`read_mothur_count()`](https://mothur.org/strollur/reference/read_mothur_count.md)
function:

``` r
sample_table <- read_mothur_count(
  filename = strollur_example("final.count_table.gz")
)
```

sample_table is a data.frame containing sequence_names, samples, and
abundances. You can add the sequence abundance data to your data set
using the [`assign()`](https://mothur.org/strollur/reference/assign.md)
function:

``` r
assign(my_data, table = sample_table, type = "sequence_abundance")
#> ℹ Assigned 2425 sequence abundances.
#> [1] 2425
my_data
#> my_data:
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
```

To add sequence taxonomy assignments, you can read a [taxonomy
file](https://mothur.org/wiki/taxonomy_file/) file using the
[`read_mothur_taxonomy()`](https://mothur.org/strollur/reference/read_mothur_taxonomy.md)
function:

``` r
classification_data <- read_mothur_taxonomy(
  taxonomy = strollur_example("final.taxonomy.gz")
)
```

classification_data is a data.frame containing sequence names and
taxonomies. You can add the sequence classification data to your data
set as follows:

``` r
assign(my_data, table = classification_data, type = "sequence_taxonomy")
#> ℹ Assigned 2425 sequence taxonomies.
#> [1] 2425
```

To assign sequences to bins, you can read a [mothur list
file](https://mothur.org/wiki/list_file/) file using the
[`read_mothur_list()`](https://mothur.org/strollur/reference/read_mothur_list.md)
function:

``` r
otu_data <- read_mothur_list(list = strollur_example("final.opti_mcc.list.gz"))
asv_data <- read_mothur_list(list = strollur_example("final.asv.list.gz"))
phylotype_data <- read_mothur_list(list = strollur_example("final.tx.list.gz"))
```

otu_data, asv_data and phylotype_data are data.frames containing bin
names and sequence names. You can add the bin data to your data set as
follows:

``` r
assign(my_data, table = otu_data, type = "bins", bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531
assign(my_data, table = asv_data, type = "bins", bin_type = "asv")
#> ℹ Assigned 2425 asv bins.
#> [1] 2425
assign(
  my_data,
  table = phylotype_data,
  type = "bins", bin_type = "phylotype"
)
#> ℹ Assigned 63 phylotype bins.
#> [1] 63
my_data
#> my_data:
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
#> Total number of otus: 531 
#> Total number of otu bin classifications: 531 
#> Total number of asvs: 2425 
#> Total number of asv bin classifications: 2425 
#> Total number of phylotypes: 63 
#> Total number of phylotype bin classifications: 63 
#> Total number of sequence classifications: 2425
```

When you assign bins to sequences with taxonomic assignments the data
set object will find the consensus taxonomy of the bins automatically.
If you wish to assign the bin taxonomy separately, you can read a
[mothur cons_taxonomy file](https://mothur.org/wiki/constaxonomy_file/)
file using the
[`read_mothur_cons_taxonomy()`](https://mothur.org/strollur/reference/read_mothur_cons_taxonomy.md)
function:

``` r
otu_taxonomy_data <- read_mothur_cons_taxonomy(
  taxonomy =
    strollur_example("final.cons.taxonomy")
)
```

otu_taxonomy_data is a data.frame containing bin names, abundances and
taxonomies. You can add the bin taxonomic data to your data set as
follows:

``` r
assign(my_data, table = otu_taxonomy_data, type = "bin_taxonomy")
#> ℹ Assigned 531 otu bin taxonomies.
#> [1] 531
```

## Writing mothur formatted file types

- [`write_mothur()`](https://mothur.org/strollur/reference/write_mothur.md)
  write mothur formatted files for all data
- [`write_fasta()`](https://mothur.org/strollur/reference/write_fasta.md)
  read a [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
  formatted sequence file
- [`write_mothur_count()`](https://mothur.org/strollur/reference/write_mothur_count.md)
  write a mothur formatted [count
  file](https://mothur.org/wiki/count_file/)
- [`write_mothur_design()`](https://mothur.org/strollur/reference/write_mothur_design.md)
  write a mothur formatted [design
  file](https://mothur.org/wiki/design_file/)
- [`write_taxonomy()`](https://mothur.org/strollur/reference/write_taxonomy.md)
  write a mothur formatted [taxonomy
  file](https://mothur.org/wiki/taxonomy_file/)
- [`write_mothur_cons_taxonomy()`](https://mothur.org/strollur/reference/write_mothur_cons_taxonomy.md)
  write a mothur formatted [cons_taxonomy
  file](https://mothur.org/wiki/constaxonomy_file/)
- [`write_mothur_list()`](https://mothur.org/strollur/reference/write_mothur_list.md)
  write a mothur formatted [list
  file](https://mothur.org/wiki/list_file/)
- [`write_mothur_shared()`](https://mothur.org/strollur/reference/write_mothur_shared.md)
  write a mothur formatted [shared
  file](https://mothur.org/wiki/shared_file/)
- [`write_mothur_rabund()`](https://mothur.org/strollur/reference/write_mothur_rabund.md)
  write a mothur formatted [rabund
  file](https://mothur.org/wiki/rabund_file/)
