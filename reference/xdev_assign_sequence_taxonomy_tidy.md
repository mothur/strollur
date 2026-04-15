# xdev_assign_sequence_taxonomy_tidy

Assign sequence classifications to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

Note, if you assign sequence taxonomies and assign bins, 'Dataset' will
find the consensus taxonomy for each bin for you.

## Usage

``` r
xdev_assign_sequence_taxonomy_tidy(
  data,
  table,
  reference = NULL,
  sequence_name = "sequence_names",
  level = "levels",
  taxonomy = "taxonomies",
  confidence = "confidences",
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
  the sequence names. Default column name is 'sequence_names'.

- level, :

  a string containing the name of the column in 'table' that contains
  the taxonomy levels. Default column name is 'levels'.

- taxonomy, :

  a string containing the name of the column in 'table' that contains
  the sequence taxonomies. Default column name is 'taxonomies'.

- confidence, :

  a string containing the name of the column in 'table' that contains
  the taxonomies confidence. Default column name is 'confidences'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

double containing the number of sequence taxonomies assigned

## Examples

``` r
sequence_classifications <- readRDS(strollur_example("miseq_tidy_taxonomy.rds"))

data <- new_dataset("my_dataset")

xdev_assign_sequence_taxonomy_tidy(data, sequence_classifications)
#> Assigned 2425 sequence taxonomies.
#> [1] 2425

# With the reference parameter you can add information about the reference
# you used to classify your sequences. You can also add references using the
# 'add_references' function.

reference <- new_reference("trainset9_032012.pds.zip", "9_032012",
              "classification by mothur2 v1.0 using default options", "",
"https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip")

xdev_assign_sequence_taxonomy_tidy(data, sequence_classifications, reference)
#> Assigned 2425 sequence taxonomies.
#> [1] 2425
```
