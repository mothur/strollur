# xdev_assign_sequence_taxonomy

Assign sequence classifications to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

Note, if you assign sequence taxonomies and assign bins, strollur will
find the consensus taxonomy for each bin for you.

## Usage

``` r
xdev_assign_sequence_taxonomy(
  data,
  table,
  reference = NULL,
  sequence_name = "sequence_name",
  taxonomy = "taxonomy",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing sequence taxonomy assignments

- reference, :

  a list created by the function \[new_reference\]. Optional.

- sequence_name, :

  a string containing the name of the column in 'table' that contains
  the sequence names. Default column name is 'sequence_name'.

- taxonomy, :

  a string containing the name of the column in 'table' that contains
  the sequence taxonomies. Default column name is 'taxonomy'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

sequence_classifications <- read_mothur_taxonomy(strollur_example(
                        "final.taxonomy.gz"))

data <- new_dataset("my_dataset", 2)

xdev_assign_sequence_taxonomy(data, sequence_classifications)
#> Assigned 2425 sequence taxonomies.
#> my_dataset:
#> 
#> data frame with 0 columns and 0 rows
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> Total number of sequence classifications: 2425 
#> 

# With the reference parameter you can add information about the reference
# you used to classify your sequences. You can also add references using the
# 'add_references' function.

reference <- new_reference("trainset9_032012.pds.zip", "9_032012",
              "classification by mothur2 v1.0 using default options", "",
"https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip")

xdev_assign_sequence_taxonomy(data, sequence_classifications, reference)
#> Assigned 2425 sequence taxonomies.
#> Added 1 resource references.
#> my_dataset:
#> 
#> data frame with 0 columns and 0 rows
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> Total number of sequence classifications: 2425 
#> Total number of resource references: 1 
#> 
```
