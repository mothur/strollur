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
phylo_obj <- write_phyloseq(miseq)
miseq_re_read <- read_phyloseq(phylo_obj)
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Added metadata.
```
