# General Importing

The strollur package stores data associated with your microbial DNA
analysis. This tutorial will familiarize you some of with the functions
available in the strollur package. If you haven’t reviewed the “Getting
Started” tutorial, we recommend you start there.

## Creating a new dataset

First let’s create an empty data set named my_data.

``` r
data <- new_dataset(dataset_name = "my_data")
```

## Importing Data

The *strollur* package includes two functions to allow you to add
microbial data.

[`add()`](https://mothur.org/strollur/reference/add.md) - The add
function allows you to add sequences, reports, metadata, and resource
references to your data set.

[`assign()`](https://mothur.org/strollur/reference/assign.md) - The
assign function allows assign sequence abundances, sequence
classifications, bins, bin representative sequences, bin
classifications, samples and treatments to your data set.

### add

The add function allows you to add sequences, reports, metadata, and
resource references to your data set.

#### Adding FASTA sequences

First, let’s add some
[FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) data.
strollur has a function for reading FASTA files named
[`read_fasta()`](https://mothur.org/strollur/reference/read_fasta.md).
We will use it to read the sequence data into a data.frame.

``` r
fasta_data <- read_fasta(strollur_example("final.fasta.gz"))

add(
  data,
  table = fasta_data,
  type = "sequences"
)
#> ℹ Added 2425 sequences.
#> [1] 2425
data
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

If you want to include a resource reference about your fasta data you
can use the `new_reference` function and the reference parameter.
strollur does not allow you to add sequences with the same name, so
let’s use the
[`clear()`](https://mothur.org/strollur/reference/clear.md) function to
remove all data from our data set.

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
  data,
  table = fasta_data,
  type = "sequences",
  reference = resource_reference
)
#> ℹ Added 2425 sequences.
#> [1] 2425
data
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
#> 
#> Total number of resource references: 1
```

#### Adding Custom Reports

You may want to add custom reports to your data set such as an contigs
assembly report, chimera report or alignment report. You can do so by
setting type = “reports”. You must also provide a report_type.

This is also a good time to explain what the table_names parameter does
for you. strollur expects the columns in custom reports to have specific
names. If your table’s names differ from what strollur is expecting, you
will see an error like that below.

``` r
contigs_report <- readr::read_tsv(strollur_example("final.contigs_report.gz"),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data,
  table = contigs_report,
  type = "reports",
  report_type = "contigs_report"
)

# Error: The report must include a column containing sequence names.
# sequence_names is not a named column in your report.

# Called from: xdev_add_report(data, table = table, type = report_type,
#    sequence_name = table_names[["sequence_name"]], verbose)
```

You can use the table_names parameter to tell strollur what the specific
column is called in your custom report table. In the contigs_report the
*sequence_name* column is called *Name*, so we will add one more line to
the add function.

``` r
contigs_report <- readr::read_tsv(strollur_example("final.contigs_report.gz"),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data,
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
#> 
#> Total number of resource references: 1 
#> Total number of custom reports: 1
```

#### Adding Metadata

Now that we have added our custom contigs assembly report, let’s learn
how to add metadata. We can add metadata to our data set by setting the
type = “metadata”.

``` r
metadata <- readr::read_tsv(strollur_example("mouse.dpw.metadata"),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data,
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
reference <- readr::read_csv(strollur_example("references.csv"),
  col_names = TRUE, show_col_types = FALSE
)

add(
  data,
  table = reference,
  type = "references"
)
#> ℹ Added 2 resource references.
#> [1] 2
```

### assign

The [`assign()`](https://mothur.org/strollur/reference/assign.md)
function allows assign sequence abundances, sequence classifications,
bins, bin representative sequences, bin classifications, samples and
treatments to your data set.

#### Assigning Abundances

After adding your FASTA sequences, you can assign abundance and sample
data using the assign function with the type = “sequence_abundance”.

``` r
abundance_table <- readr::read_tsv(
  strollur_example(
    "mothur2_count_table.tsv.gz"
  ),
  show_col_types = FALSE
)
assign(
  data,
  table = abundance_table,
  type = "sequence_abundance",
  table_names = list(sequence_name = "names")
)
#> ℹ Assigned 2425 sequence abundances.
#> [1] 2425
data
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
#> Total number of treatments: 2 
#> Total number of resource references: 3 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
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
  strollur_example(
    "mothur2_bin_assignments_list.tsv.gz"
  ),
  show_col_types = FALSE
)

assign(
  data,
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
#> Total number of resource references: 3 
#> Total number of custom reports: 1 
#> Your dataset includes metadata
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
  taxonomy = strollur_example("final.taxonomy.gz")
)

assign(
  data,
  table = sequence_classification_data,
  type = "sequence_taxonomy"
)
#> ℹ Assigned 2425 sequence taxonomies.
#> [1] 2425
```

Note, when you assign taxonomy to sequences that are assigned to bins,
strollur will automatically assign the bin taxonomies to be the
consensus taxonomy of the sequences in the bins. You can also set bin
taxonomies independently by setting the type = “bin_taxonomy”.

``` r
otu_data <- read_mothur_cons_taxonomy(strollur_example(
  "final.cons.taxonomy"
))

assign(
  data,
  table = otu_data,
  type = "bin_taxonomy",
  bin_type = "otu"
)
#> ℹ Assigned 531 otu bin taxonomies.
#> [1] 531
```

#### Assigning Bin Representatives

strollur allows you to assign a bin representative sequences to the bins
in your clusters. Let’s assign bin representatives to our *otu* bins.

``` r
bin_representatives <- readr::read_tsv(
  strollur_example(
    "otu_representative_sequences.tsv"
  ),
  show_col_types = FALSE
)

assign(
  data,
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
  strollur_example("mouse.time.design"),
  col_names = TRUE, show_col_types = FALSE
)

assign(
  data,
  table = sample_assignments,
  type = "treatments"
)
#> ℹ Assigned 19 samples to treatments.
#> [1] 19
```

## Sample Trees and Sequence Trees

Lastly, strollur allows you to add tree that relate your samples or
sequences. Let’s look at some examples together.

``` r
sample_tree <- ape::read.tree(strollur_example("final.opti_mcc.jclass.ave.tre"))
sequence_tree <- ape::read.tree(strollur_example("final.phylip.tre.gz"))

data$add_sample_tree(sample_tree)
data$add_sequence_tree(sequence_tree)

#| fig.alt: >
#|   Plot of Miseq_SOP's sample relationship tree
par(bg = "white")
ape::plot.phylo(data$get_sample_tree(),
  no.margin = TRUE,
  cex = 0.5, edge.color = "maroon", tip.color = "navy"
)
```

![](General_Importing_files/figure-html/unnamed-chunk-15-1.png)

Thanks for following along. To learn more about the functions used to
access the data in your data set, take a look at the [Accessing
Data](https://mothur.org/strollur/articles/vignettes/Accessing_Dataset.md)
tutorial.
