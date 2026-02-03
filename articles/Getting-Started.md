# Getting Started

The *rdataset* package stores the data associated with your microbial
DNA analysis. This includes nucleotide sequences, abundance, sample and
treatment assignments, taxonomic classifications, asv, otu and phylotype
clusters, metadata, trees and various reports. It is designed to
facilitate data analysis across multiple R packages with utility
functions to read from [mothur](https://mothur.org),
[qiime2](https://qiime2.org), [dada2](https://benjjneb.github.io/dada2/)
and
[phyloseq](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html).

## Installing rdataset

You can install the CRAN version with:

``` r
install.packages("rdataset")
```

You can install the development version of *rdataset* from
[GitHub](https://github.com/SchlossLab/rdataset) with:

``` r
# install.packages("devtools")
devtools::install_github("SchlossLab/rdataset")
```

## Importing Data

- [`vignette("General_Importing")`](../articles/General_Importing.md)
- [`vignette("Importing_from_mothur")`](../articles/Importing_from_mothur.md)
- `vignette("Importing_from_qiime")`

## Accessing Data

- `vignette("Using_Dataset")`

## Transfering Data

- [`vignette("Data_Transfers")`](../articles/Data_Transfers.md)
