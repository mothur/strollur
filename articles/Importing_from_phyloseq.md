# Importing_from_phyloseq

``` r
library(strollur)
#> 
#> Attaching package: 'strollur'
#> The following objects are masked from 'package:base':
#> 
#>     assign, names, summary
library(phyloseq)
```

rdataset includes functionality for the reading and writing of phyloseq
objects. To convert a phyloseq object to a rdataset object, you need to
run the
[`read_phyloseq()`](https://mothur.org/strollur/reference/read_phyloseq.md)
function.

``` r
# Using the phyloseq example data
phylo_object <- readRDS(strollur_example("GlobalPatterns.RDS"))
rdata_object <- read_phyloseq(phylo_object)
#> ℹ Added 19216 sequences.
#> ℹ Assigned 19216 sequence abundances.
#> ℹ Assigned 19216 sequence taxonomies.
#> ℹ Added metadata.
rdata_object
#> 
#> Number of unique seqs: 19216 
#> Total number of seqs: 28216678 
#> 
#> Total number of samples: 26 
#> Total number of sequence classifications: 19216 
#> Your dataset includes metadata
```

Now that are phyloseq object is converted into a rdataset object, we can
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

Furthermore, we can output rdataset objects as phyloseq objects using
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
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 19 samples to treatments.
#> ℹ Assigned 531 otu bin taxonomies.
#> ℹ Assigned 531 otu bin representative sequences.
#> ℹ Added metadata.
#> ℹ Added 2 resource references.
#> ℹ Added a contigs_report.
miseq_phyloseq <- write_phyloseq(miseq)
miseq_phyloseq
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 2425 taxa and 19 samples ]
#> sample_data() Sample Data:       [ 19 samples by 2 sample variables ]
#> tax_table()   Taxonomy Table:    [ 2425 taxa by 6 taxonomic ranks ]
#> phy_tree()    Phylogenetic Tree: [ 2425 tips and 2424 internal nodes ]
```
