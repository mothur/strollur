# Create a [strollur](https://mothur.org/strollur/reference/strollur.html) object from mothur outputs

The read_mothur function reads various [file
types](https://mothur.org/wiki/tags/#file_types) created by mothur, and
creates a \`strollur\` object.

To generate the various input files you can follow Pat's [Miseq example
analysis](https://mothur.org/wiki/miseq_sop/).

## Usage

``` r
read_mothur(
  fasta = NULL,
  count = NULL,
  taxonomy = NULL,
  otu_list = NULL,
  asv_list = NULL,
  phylo_list = NULL,
  design = NULL,
  cons_taxonomy = NULL,
  otu_shared = NULL,
  asv_shared = NULL,
  phylo_shared = NULL,
  sample_tree = NULL,
  sequence_tree = NULL,
  dataset_name = ""
)
```

## Arguments

- fasta:

  filename, a FASTA formatted file containing sequence strings. [fasta
  file](https://mothur.org/wiki/fasta_file/)

- count:

  filename, a mothur [count file](https://mothur.org/wiki/count_file/)

- taxonomy:

  filename, a mothur [taxonomy
  file](https://mothur.org/wiki/taxonomy_file/), created by
  [classify.seqs](https://mothur.org/wiki/classify.seqs/)

- otu_list:

  filename, a mothur [list file](https://mothur.org/wiki/list_file/)
  containing otu bin assignments. The otu_list file is created by
  [cluster](https://mothur.org/wiki/cluster/),
  [cluster.split](https://mothur.org/wiki/cluster.split/), and
  [cluster.fit](https://mothur.org/wiki/cluster.fit/)

- asv_list:

  filename, a mothur [list file](https://mothur.org/wiki/list_file/)
  containing asv bin assignments. The asv_list file is created by
  [cluster](https://mothur.org/wiki/cluster/) using the 'unique' method.

- phylo_list:

  filename, a mothur [list file](https://mothur.org/wiki/list_file/)
  containing phylotype bin assignments. The phylo_list file is created
  by [phylotype](https://mothur.org/wiki/phylotype/).

- design:

  filename, a mothur [design file](https://mothur.org/wiki/design_file/)

- cons_taxonomy:

  filename, a mothur consensus taxonomy file [constaxonomy
  file](https://mothur.org/wiki/constaxonomy_file/). The cons_taxonomy
  file is created by
  [classify.otu](https://mothur.org/wiki/classify.otu/).

- otu_shared:

  filename, a mothur [shared file](https://mothur.org/wiki/shared_file/)
  containing otu bin sample abundance assignments.

- asv_shared:

  filename, a mothur [shared file](https://mothur.org/wiki/shared_file/)
  containing asv bin sample abundance assignments.

- phylo_shared:

  filename, a mothur [shared file](https://mothur.org/wiki/shared_file/)
  containing phylotype bin sample abundance assignments.

- sample_tree:

  filename, a tree that relates samples. The sample tree is created by
  [tree.shared](https://mothur.org/wiki/tree.shared/). We recommend
  running tree.shared with subsample = true, and using the 'ave.tre'
  output for best results.

- sequence_tree:

  filename, a tree that relates sequences. The sequence tree is created
  by [clearcut](https://mothur.org/wiki/clearcut/). We DO NOT recommend
  using sequence trees. With the ever growing size of modern datasets,
  sequence tree can be difficult / impossible to build without hitting a
  memory limitation.

- dataset_name:

  A string containing a name for your dataset.

## Value

A [strollur](https://mothur.org/strollur/reference/strollur.html) object

## Note

- *consensus taxonomy*, The \`strollur\` object will generate consensus
  taxonomies for you based on the sequence taxonomy assignment. You only
  need to provide the ".cons.taxonomy" file if you are not providing
  sequence taxonomy assignments.

- *shared / rabund file*, The \`strollur\` object will generate shared
  and rabund data for you based on the otu assignment in the list file
  and the count data. You only need to provide the ".shared" file if you
  are not providing the list and count files.

## Examples

``` r
# For dataset's including sequence data:

data <- read_mothur(
  fasta = strollur_example("final.fasta.gz"),
  count = strollur_example("final.count_table.gz"),
  taxonomy = strollur_example("final.taxonomy.gz"),
  design = strollur_example("mouse.time.design"),
  otu_list = strollur_example("final.opti_mcc.list.gz"),
  asv_list = strollur_example("final.asv.list.gz"),
  phylo_list = strollur_example("final.tx.list.gz"),
  sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
  dataset_name = "miseq_sop"
)
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.

# For dataset's with only otu data:

data <- read_mothur(
  otu_shared = strollur_example("final.opti_mcc.shared"),
  cons_taxonomy = strollur_example(
    "final.cons.taxonomy"
  ),
  design = strollur_example("mouse.time.design"),
  sample_tree = strollur_example("final.opti_mcc.jclass.ave.tre"),
  dataset_name = "miseq_sop"
)
#> Assigned 531 otu bins.
#> Assigned 19 samples to treatments.
#> Assigned 531 otu bin taxonomies.
```
