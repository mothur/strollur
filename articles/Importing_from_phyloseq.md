# Importing from phyloseq

``` r
library(strollur)
#> 
#> Attaching package: 'strollur'
#> The following objects are masked from 'package:base':
#> 
#>     assign, names, summary
library(phyloseq)
```

`strollur` includes functionality for the reading and writing of
phyloseq objects. To convert a phyloseq object to a strollur object, you
need to run the
[`read_phyloseq()`](https://mothur.org/strollur/reference/read_phyloseq.md)
function.

``` r
# Using the phyloseq example data
phylo_object <- readRDS(strollur_example("GlobalPatterns.RDS"))
rdata_object <- read_phyloseq(phylo_object)
#> Added 19216 sequences.
#> Assigned 19216 sequence abundances.
#> Assigned 19216 sequence taxonomies.
#> Added metadata.
rdata_object
#> 
#> Number of unique seqs: 19216 
#> Total number of seqs: 28216678 
#> 
#> Total number of samples: 26 
#> Total number of sequence classifications: 19216 
#> Your dataset includes metadata
```

Now that are phyloseq object is converted into a strollur object, we can
utilize functions like
[`count()`](https://mothur.org/strollur/reference/count.md),
[`abundance()`](https://mothur.org/strollur/reference/abundance.md), and
[`names()`](https://mothur.org/strollur/reference/names.md) to inspect
the data.

``` r
count(rdata_object, "samples")
#> [1] 26
head(names(rdata_object, "sequences"))
#> [1] "549322" "522457" "951"    "244423" "586076" "246140"
head(abundance(rdata_object, "sequences"))
#>   sequence_names abundances
#> 1         549322        259
#> 2         522457          8
#> 3            951          1
#> 4         244423         51
#> 5         586076          3
#> 6         246140          4
```

Furthermore, we can output strollur objects as phyloseq objects using
the
[`write_phyloseq()`](https://mothur.org/strollur/reference/write_phyloseq.md)
function.

``` r
phyloseq_object <- write_phyloseq(rdata_object)
phyloseq_object
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 19216 taxa and 26 samples ]
#> sample_data() Sample Data:       [ 26 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 19216 taxa by 7 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 19216 tips and 19215 internal nodes ]


# With the miseq example data
miseq <- miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added metadata.
#> Added 2 resource references.
#> Added a contigs_report.
miseq_phyloseq <- write_phyloseq(miseq)
miseq_phyloseq
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 2425 taxa and 19 samples ]
#> sample_data() Sample Data:       [ 19 samples by 2 sample variables ]
#> tax_table()   Taxonomy Table:    [ 2425 taxa by 6 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 2425 tips and 2424 internal nodes ]
```
