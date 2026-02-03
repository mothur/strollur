
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rdataset <img src="man/figures/logo.png" align="right" height="136" alt="rdataset website" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/mothur/rdataset/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mothur/rdataset/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/SchlossLab/rdataset/branch/main/graph/badge.svg)](https://app.codecov.io/gh/SchlossLab/rdataset?branch=main)
[![pkgdown](https://github.com/mothur/rdataset/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/mothur/rdataset/actions/workflows/pkgdown.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/rdataset)](https://CRAN.R-project.org/package=rdataset)
<!-- badges: end -->

## Overview

The **rdataset** package stores the data associated with your microbial
DNA analysis. This includes nucleotide sequences, abundance, sample and
treatment assignments, taxonomic classifications, sequence bin
assignments, metadata, trees and various reports. It is designed to
facilitate data analysis across multiple R packages with utility
functions to import from [mothur](https://mothur.org),
[qiime2](https://qiime2.org), [dada2](https://benjjneb.github.io/dada2/)
and
[phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html).

- `add()` add sequences, reports, metadata, and resource references
- `assign()` assign abundances, classifications, bins, samples and
  treatments and more
- `names()` get the names of sequences, bins, samples, treatments and
  reports
- `count()` get the number of sequences, bins, samples and treatments
- `abundance()` get the abundances for sequences, bins, samples, and
  treatments
- `report()` get
  [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) sequences,
  sequence and classification reports, bin assignments, sample
  assignments, metadata, sequence data reports, custom reports, resource
  references and scrapped data reports.
- `summary()` summarize sequences, your custom reports, and scrapped
  data

## Installation

You can install the CRAN version with:

``` r
install.packages("rdataset")
```

### Development version

You can install the development version of rdataset from
[GitHub](https://github.com/mothur/rdataset) with:

``` r
# install.packages("devtools")
devtools::install_github("SchlossLab/rdataset")
```

## Usage

The example below adds
[FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) sequence
data, assigns sequence abundance, samples and treatments, as well as
assigning bins and taxonomic data.

``` r
library(rdataset)

data <- new_dataset(dataset_name = "microbial RNA example")

fasta_data <- read_fasta(rdataset_example("final.fasta.gz"))

add(
  data = data,
  table = fasta_data,
  type = "sequences"
)
#> ℹ Added 2425 sequences.
#> [1] 2425

abundance_table <- readr::read_tsv(rdataset_example("mothur2_count_table.tsv.gz"),
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

bin_table <- readr::read_tsv(
  rdataset_example(
    "mothur2_bin_assignments_list.tsv.gz"
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

sequence_classification_data <- read_mothur_taxonomy(
  taxonomy = rdataset_example("final.taxonomy.gz")
)

assign(
  data = data,
  table = sequence_classification_data,
  type = "sequence_taxonomy"
)
#> ℹ Assigned 2425 sequence taxonomies.
#> [1] 2425
data
#> microbial RNA example:
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
```

## Getting help

If you encounter an issue, please file an issue on
[GitHub](https://github.com/mothur/rdataset/issues). Please include a
minimal reproducible example with your issue.

## Contributing

Is there a feature you’d like to see included, please let us know! Pull
requests are welcome on
[GitHub](https://github.com/mothur/rdataset/pulls).
