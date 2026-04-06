
<!-- README.md is generated from README.Rmd. Please edit that file -->

# strollur <a href="https://mothur.org/strollur/"><img src="man/figures/logo.png" align="right" height="139" alt="strollur website" /></a>

<!-- # strollur <a href="https://mothur.org/strollur/"><img src="man/figures/logo.png" align="right" height="200" style="float:right; height:200px;" alt="strollur website"> -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/mothur/strollur/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mothur/strollur/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mothur/strollur/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mothur/strollur?branch=main)
[![pkgdown](https://github.com/mothur/strollur/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/mothur/strollur/actions/workflows/pkgdown.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/strollur)](https://CRAN.R-project.org/package=strollur)
<!-- badges: end -->

## Overview

The **strollur** package stores the data associated with your microbial
DNA analysis. This includes nucleotide sequences, abundance, sample and
treatment assignments, taxonomic classifications, sequence bin
assignments, metadata, trees and various reports. It is designed to
facilitate data analysis across multiple R packages with utility
functions to import from [mothur](https://mothur.org),
[qiime2](https://qiime2.org), [dada2](https://benjjneb.github.io/dada2/)
and
[phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html).

- `add()` adds sequences, reports, metadata, and resource references
- `assign()` assigns abundances, classifications, bins, samples and
  treatments and more
- `names()` gets the names of sequences, bins, samples, treatments and
  reports
- `count()` gets the number of sequences, bins, samples and treatments
- `abundance()` gets the abundances for sequences, bins, samples, and
  treatments
- `report()` gets
  [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) sequences,
  sequence and classification reports, bin assignments, sample
  assignments, metadata, sequence data reports, custom reports, resource
  references and scrapped data reports.
- `summary()` summarizes sequences, your custom reports, and scrapped
  data

## Installation

You can install the CRAN version with:

``` r
library(strollur)
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
```

## Getting help

If you encounter an issue, please file an issue on
[GitHub](https://github.com/mothur/strollur/issues). Please include a
minimal reproducible example with your issue.

## Contributing

Is there a feature you’d like to see included, please let us know! Pull
requests are welcome on
[GitHub](https://github.com/mothur/strollur/pulls).
