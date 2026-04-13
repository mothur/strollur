# Create a [strollur](https://mothur.org/strollur/reference/strollur.html) object from dada2 outputs

This function reads a dada2 sequence table and creates a \`strollur\`
object. The dada2 sequence table is a 2D matrix containing the abundance
counts by sample for each ASV. The sample names are stored as row names
and the sequence nucleotide strings are stored as column names.

To generate the dada2 sequence table from your own files you can follow
[this dada2 tutorial](https://benjjneb.github.io/dada2/tutorial.html).

## Usage

``` r
read_dada2(sequence_table, dataset_name = "")
```

## Arguments

- sequence_table:

  A dada2 sequence table

- dataset_name:

  A string containing a name for your dataset.

## Value

A \`strollur\` object

## Examples

``` r
seqtab <- readRDS(strollur_example("dada2.rds"))
dim(seqtab)
#> [1]  19 279

data <- read_dada2(sequence_table = seqtab, dataset_name = "dada2 example")
#> ℹ Added 279 sequences.
#> ℹ Assigned 279 sequence abundances.
```
