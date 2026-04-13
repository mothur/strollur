# strollur

## Overview

The **`strollur`** package stores the data associated with your
microbial DNA analysis. This includes nucleotide sequences, abundance,
sample and treatment assignments, taxonomic classifications, sequence
bin assignments, metadata, trees and various reports. It is designed to
facilitate data analysis across multiple R packages with utility
functions to import from [mothur](https://mothur.org),
[qiime2](https://qiime2.org), [dada2](https://benjjneb.github.io/dada2/)
and
[phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html).

- [`add()`](https://mothur.org/strollur/reference/add.md) adds
  sequences, reports, metadata, and resource references
- [`assign()`](https://mothur.org/strollur/reference/assign.md) assigns
  abundances, classifications, bins, samples and treatments and more
- [`names()`](https://mothur.org/strollur/reference/names.md) gets the
  names of sequences, bins, samples, treatments and reports
- [`count()`](https://mothur.org/strollur/reference/count.md) gets the
  number of sequences, bins, samples and treatments
- [`abundance()`](https://mothur.org/strollur/reference/abundance.md)
  gets the abundances for sequences, bins, samples, and treatments
- [`report()`](https://mothur.org/strollur/reference/report.md) gets
  [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) sequences,
  sequence and classification reports, bin assignments, sample
  assignments, metadata, sequence data reports, custom reports, resource
  references and scrapped data reports.
- [`summary()`](https://mothur.org/strollur/reference/summary.md)
  summarizes sequences, your custom reports, and scrapped data

## Installation

You can install the CRAN version with:

``` r
install.packages("strollur")
```

### Development version

You can install the development version of strollur from
[GitHub](https://github.com/mothur/strollur) with:

``` r
# install.packages("devtools")
devtools::install_github("mothur/strollur")
```

## Usage

The example below adds
[FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) sequence
data, assigns sequence abundance, samples and treatments, as well as
assigning bins and taxonomic data to a [strollur
object](https://mothur.org/strollur/reference/strollur.html).

``` r
fasta_data <- read_fasta(strollur_example("final.fasta.gz"))
abundance_table <- readRDS(strollur_example("miseq_abundance_by_sample.rds"))
bin_table <- readRDS(strollur_example("miseq_list_otu.rds"))
classification_data <- read_mothur_taxonomy(taxonomy = strollur_example("final.taxonomy.gz"))

data <- new_dataset(dataset_name = "microbial RNA example")

add(data, table = fasta_data, type = "sequences")
#> ℹ Added 2425 sequences.
#> [1] 2425
assign(data, table = abundance_table, type = "sequence_abundance")
#> ℹ Assigned 2425 sequence abundances.
#> [1] 2425
assign(data, table = bin_table, type = "bins", bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531
assign(data, table = classification_data, type = "sequence_taxonomy")
#> ℹ Assigned 2425 sequence taxonomies.
#> [1] 2425

data
#> microbial RNA example:
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
#> Total number of sequence classifications: 2425
```

## Getting help

If you encounter an issue, please file an issue on
[GitHub](https://github.com/mothur/strollur/issues). Please include a
minimal reproducible example with your issue.

## Contributing

Is there a feature you’d like to see included, please let us know! Pull
requests are welcome on
[GitHub](https://github.com/mothur/strollur/pulls).

## Code of Conduct

Please note that the strollur project is released with a [Contributor
Code of Conduct](https://mothur.org/strollur/CODE_OF_CONDUCT.html). By
contributing to this project, you agree to abide by its terms.
