# read_dada2

This function reads a dada2 sequence table and creates a \`strollur\`
object. The dada2 sequence table is a 2D matrix containing the abundance
counts by sample for each ASV. The sample names are stored as row names
and the sequence nucleotide strings are stored as column names.

To generate the dada2 sequence table you can follow [this dada2
tutorial](https://benjjneb.github.io/dada2/tutorial.html).

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
if (FALSE) { # \dontrun{
# Follow the dada2 tutorial above to generate seqtab from your fastq files.

# dim(seqtab)
# [1]  20 293

# data <- read_dada2(seqtab, "dada2")
} # }
```
