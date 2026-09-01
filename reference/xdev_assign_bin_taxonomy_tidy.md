# Assign bin classifications to a [strollur object](https://mothur.org/strollur/reference/strollur.html)

Assign bin classifications to a [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Usage

``` r
xdev_assign_bin_taxonomy_tidy(
  data,
  table,
  bin_type = "otu",
  reference = NULL,
  bin_name = "bin_name",
  level = "level",
  taxonomy = "taxonomy",
  confidence = "confidence",
  verbose = TRUE
)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- table:

  a data.frame containing bin taxonomy assignments

- bin_type:

  a string indicating the type of bin assignments. Default "otu".

- reference:

  a list created by the function
  [new_reference](https://mothur.org/strollur/reference/new_reference.md).
  Optional.

- bin_name:

  a string containing the name of the column in 'table' that contains
  the bin names. Default column name is 'bin_name'.

- level:

  a string containing the name of the column in 'table' that contains
  the taxonomy levels. Default column name is 'level'.

- taxonomy:

  a string containing the name of the column in 'table' that contains
  the bin taxonomies. Default column name is 'taxonomy'.

- confidence:

  a string containing the name of the column in 'table' that contains
  the taxonomies confidence. Default column name is 'confidence'.

- verbose:

  a logical whether or not you want progress messages. Default = TRUE.

## Value

an updated [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Examples

``` r

bin_classifications <- readRDS(strollur_example("miseq_tidy_bin_taxonomy.rds"))
str(bin_classifications)
#> 'data.frame':    3186 obs. of  4 variables:
#>  $ bin_name  : chr  "Otu001" "Otu001" "Otu001" "Otu001" ...
#>  $ level     : int  1 2 3 4 5 6 1 2 3 4 ...
#>  $ taxonomy  : chr  "Bacteria" "\"Bacteroidetes\"" "\"Bacteroidia\"" "\"Bacteroidales\"" ...
#>  $ confidence: int  100 100 100 100 100 100 100 100 100 100 ...

data <- strollur::read_mothur(otu_shared =
                                strollur_example("final.opti_mcc.shared"))
#> Assigned 531 otu bins.

xdev_assign_bin_taxonomy_tidy(data, bin_classifications)
#> Assigned 531 otu bin taxonomies.
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of otus: 531 
#> Total number of otu bin classifications: 531 
#> 

# With the reference parameter you can add information about the reference
# you used to classify your bins. You can also add references using the
# 'strollur::xdev_add_references' function.

reference <- new_reference("trainset9_032012.pds.zip", "9_032012",
              "classification_otu by mothur2 v1.0.0 using default options",
               "",
"https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip")

xdev_assign_bin_taxonomy_tidy(data, bin_classifications,
                                 reference = reference)
#> Assigned 531 otu bin taxonomies.
#> Added 1 resource references.
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of otus: 531 
#> Total number of otu bin classifications: 531 
#> Total number of resource references: 1 
#> 
```
