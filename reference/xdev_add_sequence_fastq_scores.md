# Add FASTQ data to a [strollur object](https://mothur.org/strollur/reference/strollur.html)

Add FASTQ data to a [strollur
object](https://mothur.org/strollur/reference/strollur.html).

## Usage

``` r
xdev_add_sequence_fastq_scores(
  data,
  table,
  reference = NULL,
  sequence_name = "sequence_name",
  sequence = "sequence",
  quality_score = "quality_score",
  verbose = TRUE
)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- table:

  a data.frame containing FASTQ data you are trying to add

- reference:

  a list created by the function
  [new_reference](https://mothur.org/strollur/reference/new_reference.md).
  Optional.

- sequence_name:

  a string containing the name of the column in 'table' that contains
  the sequence names. Default column name is 'sequence_name'.

- sequence:

  a string containing the name of the column in 'table' that contains
  the sequence nucleotide strings. Default column name is 'sequence'.
  (Optional)

- quality_score:

  a string containing the name of the column in 'table' that contains
  the quality scores. Default column name is 'quality_score'.

- verbose:

  a logical whether or not you want progress messages. Default = TRUE.

## Value

an updated [strollur
object](https://mothur.org/strollur/reference/strollur.html)

## Examples

``` r

table <- strollur::read_fastq(strollur_example("tiny.fastq.gz"))
data <- strollur::new_dataset("example")
strollur::xdev_add_sequence_fastq_scores(data, table)
#> Added 3 sequences.
#> Assigned 3 quality scores.
#> example:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  251    251      0        4     0    1.00
#> 2.5%-tile:       1  251    251      0        4     0    1.05
#> 25%-tile:        1  251    251      0        4     0    1.50
#> Median:          1  251    251      0        4     0    2.00
#> 75%-tile:        1  251    251      0        4     0    2.50
#> 97.5%-tile:      1  251    251      0        4     0    2.95
#> Maximum:         1  251    251      1        4     1    3.00
#> Mean:            1  251    251      0        4     0    2.00
#> 
#> Number of unique seqs: 3 
#> Total number of seqs: 3 
#> 
#> 
data
#> example:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  251    251      0        4     0    1.00
#> 2.5%-tile:       1  251    251      0        4     0    1.05
#> 25%-tile:        1  251    251      0        4     0    1.50
#> Median:          1  251    251      0        4     0    2.00
#> 75%-tile:        1  251    251      0        4     0    2.50
#> 97.5%-tile:      1  251    251      0        4     0    2.95
#> Maximum:         1  251    251      1        4     1    3.00
#> Mean:            1  251    251      0        4     0    2.00
#> 
#> Number of unique seqs: 3 
#> Total number of seqs: 3 
#> 
#> 
```
