# Create a [strollur](https://mothur.org/strollur/reference/strollur.html) object from a phyloseq object

The \`read_phyloseq()\` function reads phyloseq objects created from the
phyloseq package
(https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
and converts it into a strollur object.

## Usage

``` r
read_phyloseq(phyloseq_object, treatment_column_name = NULL, dataset_name = "")
```

## Arguments

- phyloseq_object:

  the phyloseq object that is returned when using any read function in
  the phyloseq package. It has to be of type "phyloseq"

- treatment_column_name:

  the column name inside your phyloseq object within your sample data
  that is used to descrbe treatments. It must be a character. Defaults
  to NULL.

- dataset_name:

  A string containing a name for your dataset.

## Value

a strollur object.

## Examples

``` r
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
phylo_obj <- write_phyloseq(miseq)
miseq_re_read <- read_phyloseq(phylo_obj)
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Added metadata.
```
