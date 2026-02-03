# General Importing

The rdataset package stores data associated with your microbial DNA
analysis. This tutorial will familiarize you some of with the functions
available in the rdataset package. If you haven’t reviewed the “Getting
Started” tutorial, we recommend you start there.

## Creating a new dataset

First let’s create an empty data set named my_data.

``` r
data <- new_dataset(dataset_name = "my_data")
```

## Importing Data

The *rdataset* package includes two functions to allow you to add
microbial data.

[`add()`](../reference/add.md) - The add function allows you to add
sequences, reports, metadata, and resource references to your data set.

[`assign()`](../reference/assign.md) - The assign function allows assign
sequence abundances, sequence classifications, bins, bin representative
sequences, bin classifications, samples and treatments to your data set.

### add

The add function allows you to add sequences, reports, metadata, and
resource references to your data set.

- **Parameters:**
  - *data* - a data set object
  - *table* - a data.frame containing the data you wish to add
  - *type* - string containing the type of data you are adding
  - *report_type* - string containing the report type you would like to
    add
  - *table_names* - a named list used to indicate the names of the
    columns in the table
  - *reference* - a list created by the
    [`new_reference()`](../reference/new_reference.md) function
  - *verbose* - Boolean, indicating whether you want outputs about what
    is being added

#### Adding FASTA sequences

First, let’s add some
[FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) data.
rdataset has a function for reading FASTA files named
[`read_fasta()`](../reference/read_fasta.md). We will use it to read the
sequence data into a data.frame.

``` r
fasta_data <- read_fasta(rdataset_example("final.fasta"))

add(
  data = data,
  table = fasta_data,
  type = "sequences"
)
#> ℹ Added 2425 sequences.
#> [1] 2425
data
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

If you want to include a resource reference about your fasta data you
can use the `new_reference` function and the reference parameter.
rdataset does not allow you to add sequences with the same name, so
let’s use the [`clear()`](../reference/clear.md) function to remove all
data from our data set.

``` r
clear(data)

resource_url <- "https://mothur.org/wiki/silva_reference_files/"

resource_reference <- new_reference(
  reference_name = "silva.bacteria.fasta",
  reference_version = "1.38.1",
  reference_usage = "alignment by mothur2 v1.0",
  reference_note = "aligned with default options",
  reference_url = resource_url
)

add(
  data = data,
  table = fasta_data,
  type = "sequences",
  reference = resource_reference
)
#> ℹ Added 2425 sequences.
#> [1] 2425
data
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

#### Adding Custom Reports

You may want to add custom reports to your data set such as an contigs
assembly report, chimera report or alignment report. You can do so by
setting type = “reports”. You must also provide a report_type.

This is also a good time to explain what the table_names parameter does
for you. rdataset expects the columns in custom reports to have specific
names. If your table’s names differ from what rdataset is expecting, you
will see an error like that below.

``` r
contigs_report <- readr::read_tsv(rdataset_example("final.contigs_report"),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data = data,
  table = contigs_report,
  type = "reports",
  report_type = "contigs_report"
)

# Error: The report must include a column containing sequence names.
# sequence_names is not a named column in your report.

# Called from: xdev_add_report(data = data, table = table, type = report_type,
#    sequence_name = table_names[["sequence_name"]], verbose)
```

You can use the table_names parameter to tell rdataset what the specific
column is called in your custom report table. In the contigs_report the
*sequence_name* column is called *Name*, so we will add one more line to
the add function.

``` r
contigs_report <- readr::read_tsv(rdataset_example("final.contigs_report"),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data = data,
  table = contigs_report,
  type = "reports",
  report_type = "contigs_report",
  table_names = list(sequence_name = "Name")
)
#> ℹ Added a contigs_report.
#> [1] 1
data
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
#> 
#> contigs_report :
#>             Expected_Errors   Length MisMatches Num_Ns Overlap_End
#> Minimum:       0.0000452496 250.0000   0.000000      0    248.0000
#> 2.5%-tile:     0.0010250499 252.0000   0.000000      0    250.0000
#> 25%-tile:      0.0022657500 252.0000   0.000000      0    251.0000
#> Median:        0.0092338603 253.0000   1.000000      0    251.0000
#> 75%-tile:      0.0559640005 253.0000   5.000000      0    251.0000
#> 97.5%-tile:    0.4990670085 254.0000  26.000000      0    253.0000
#> Maximum:       3.0126200000 270.0000 120.000000      0    256.0000
#> Mean:          0.0738509483 252.7575   5.162474      0    251.1555
#>             Overlap_Length Overlap_Start
#> Minimum:          232.0000      0.000000
#> 2.5%-tile:        246.0000      1.000000
#> 25%-tile:         249.0000      2.000000
#> Median:           249.0000      2.000000
#> 75%-tile:         250.0000      2.000000
#> 97.5%-tile:       251.0000      4.000000
#> Maximum:          255.0000     22.000000
#> Mean:             249.1501      2.005361
#> Unique seqs:  2425 
#> Total seqs:   2425
#> → Your dataset does not include sample data, ignoring.
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425
```

#### Adding Metadata

Now that we have added our custom contigs assembly report, let’s learn
how to add metadata. We can add metadata to our data set by setting the
type = “metadata”.

``` r
metadata <- readr::read_tsv(rdataset_example("mouse.dpw.metadata"),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data = data,
  table = metadata,
  type = "metadata"
)
#> ℹ Added metadata.
#> [1] 1
```

#### Adding Resource References

We can add additional resource references to our data set by setting the
type = “references”.

``` r
reference <- readr::read_csv(rdataset_example("references.csv"),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data = data,
  table = reference,
  type = "references"
)
#> ℹ Added 2 resource references.
#> [1] 2
```

### assign

The [`assign()`](../reference/assign.md) function allows assign sequence
abundances, sequence classifications, bins, bin representative
sequences, bin classifications, samples and treatments to your data set.

- **Parameters:**
  - *data* - a data set object
  - *table* - a data.frame containing the data you wish to assign
  - *type* - string containing the type of data you are assigning
  - *bin_type* - string containing the bin type you would like to assign
  - *table_names* - a named list used to indicate the names of the
    columns in the table
  - *reference* - a list created by the
    [`new_reference()`](../reference/new_reference.md) function
  - *verbose* - boolean, indicating whether you want outputs about what
    is being added

#### Assigning Abundances

After adding your FASTA sequences, you can assign abundance and sample
data using the assign function with the type = “sequence_abundance”.

``` r
abundance_table <- readr::read_tsv(rdataset_example("mothur2_count_table.tsv"),
  show_col_types = FALSE
)
assign(
  data = data,
  table = abundance_table,
  type = "sequence_abundance",
  table_names = list(sequence_name = "names")
)
#> ℹ Assigned 2425 sequence abundances.
#> [1] 2425
data
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
```

#### Assigning Bins

As you can see we now have abundances, samples and treatments added to
the data set. Next, let’s assign the sequences to bins using type =
“bins”. When you assign sequences to bins you must provide a *bin_type*.
The bin_type is a tag of your choosing used to reference the bin
clusters you are adding. Let’s add some [Operational Taxonomic
Unit](https://en.wikipedia.org/wiki/Operational_taxonomic_unit) clusters
and set the bin_type = “otu”.

``` r
bin_table <- readr::read_tsv(
  rdataset_example(
    "mothur2_bin_assignments_list.tsv"
  ),
  show_col_types = FALSE
)

assign(
  data = data,
  table = bin_table,
  type = "bins",
  bin_type = "otu",
  table_names = list(bin_name = "otu_id", sequence_name = "seq_id")
)
#> ℹ Assigned 531 otu bins.
#> [1] 531
data
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
```

You can see from the summary, we now have 531 otus in our dataset.

Note, if you are importing data from packages that preprocess the
microbial data into features, you can assign the feature table
abundances as sequence abundances and then assign the features to
[Amplicon Sequence
Variant](https://en.wikipedia.org/wiki/Amplicon_sequence_variant)
clusters or *asv* bins.

#### Assigning Taxonomic Classifications

Now that we have assigned our sequences to bins, let’s assign taxonomy
to our sequences.

``` r
sequence_classification_data <- read_mothur_taxonomy(
  taxonomy = rdataset_example("final.taxonomy")
)

assign(
  data = data,
  table = sequence_classification_data,
  type = "sequence_taxonomy"
)
#> ℹ Assigned 2425 sequence taxonomies.
#> [1] 2425
```

Note, when you assign taxonomy to sequences that are assigned to bins,
rdataset will automatically assign the bin taxonomies to be the
consensus taxonomy of the sequences in the bins. You can also set bin
taxonomies independently by setting the type = “bin_taxonomy”.

``` r
otu_data <- read_mothur_cons_taxonomy(rdataset_example(
  "final.cons.taxonomy"
))

assign(
  data = data,
  table = otu_data,
  type = "bin_taxonomy",
  bin_type = "otu"
)
#> ℹ Assigned 531 otu bin taxonomies.
#> [1] 531
```

#### Assigning Bin Representatives

rdataset allows you to assign a bin representative sequences to the bins
in your clusters. Let’s assign bin representatives to our *otu* bins.

``` r
bin_representatives <- readr::read_tsv(
  rdataset_example(
    "otu_representative_sequences.tsv"
  ),
  show_col_types = FALSE
)

assign(
  data = data,
  table = bin_representatives,
  type = "bin_representatives"
)
#> ℹ Assigned 531 otu bin representative sequences.
#> [1] 531
```

#### Assigning Treatments

In our case the abundance_table included treatment assignments, but you
can also assign samples to treatments by setting type = “treatments”.

``` r
sample_assignments <- readr::read_table(
  rdataset_example("mouse.time.design"),
  col_names = TRUE, show_col_types = FALSE
)

assign(
  data = data,
  table = sample_assignments,
  type = "treatments"
)
#> ℹ Assigned 19 samples to treatments.
#> [1] 19
```

## Sample Trees and Sequence Trees

Lastly, rdataset allows you to add tree that relate your samples or
sequences. Let’s look at some examples together.

``` r
sample_tree <- ape::read.tree(rdataset_example("final.opti_mcc.jclass.ave.tre"))
sequence_tree <- ape::read.tree(rdataset_example("final.phylip.tre"))

data$add_sample_tree(sample_tree)
data$add_sequence_tree(sequence_tree)

#| fig.alt: >
#|   Plot of Miseq_SOP's sample relationship tree
#|
plot(data$get_sample_tree(),
  no.margin = TRUE,
  cex = 0.5, edge.color = "maroon", tip.color = "navy"
)
```

![](General_Importing_files/figure-html/unnamed-chunk-15-1.png)

Thanks for following along. To learn more about the functions used to
access the data in your data set, take a look at the [Accessing
Data](vignettes/Using_Dataset.md) tutorial.
