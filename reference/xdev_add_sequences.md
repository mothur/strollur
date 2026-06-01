# xdev_add_sequences

Add sequence data to a
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Usage

``` r
xdev_add_sequences(
  data,
  table,
  reference = NULL,
  sequence_name = "sequence_name",
  sequence = "sequence",
  comment = "comment",
  verbose = TRUE
)
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- table, :

  a data.frame containing names, sequences(optional) and
  comments(optional).

- reference, :

  a list created by the function \[new_reference\]. Optional.

- sequence_name, :

  a string containing the name of the column in 'table' that contains
  the sequence names. Default column name is 'sequence_name'.

- sequence, :

  a string containing the name of the column in 'table' that contains
  the sequence nucleotide strings. Default column name is 'sequence'.

- comment, :

  a string containing the name of the column in 'table' that contains
  the sequence comments. Default column name is 'comment'.

- verbose, :

  a boolean whether or not you want progress messages. Default = TRUE.

## Value

an updated
[strollur](https://mothur.org/strollur/reference/strollur.html) object

## Examples

``` r

 data <- new_dataset("miseq_sop")
 fasta_data <- read_fasta(strollur_example("final.fasta.gz"))
 xdev_add_sequences(data, fasta_data)
#> Added 2425 sequences.
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  375    249      0        3     0    1.00
#> 2.5%-tile:       1  375    252      0        4     0   61.62
#> 25%-tile:        1  375    252      0        4     0  607.25
#> Median:          1  375    253      0        4     0 1213.50
#> 75%-tile:        1  375    253      0        5     0 1819.75
#> 97.5%-tile:      1  375    254      0        6     0 2365.38
#> Maximum:         1  375    256      0        6     0 2426.00
#> Mean:            1  375    252      0        4     0 1213.50
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> 

# With the additional parameters to add information about the reference

 data <- new_dataset("miseq_sop")
 fasta_data <- read_fasta(strollur_example("final.fasta.gz"))

 xdev_add_sequences(data, fasta_data,
               new_reference("silva.bacteria.fasta",
               "1.38.1",
               "alignment by mothur2 v1.0 using default options",
               "https://mothur.org/wiki/silva_reference_files/"))
#> Added 2425 sequences.
#> Added 1 resource references.
#> miseq_sop:
#> 
#>             starts ends nbases ambigs polymers numns numseqs
#> Minimum:         1  375    249      0        3     0    1.00
#> 2.5%-tile:       1  375    252      0        4     0   61.62
#> 25%-tile:        1  375    252      0        4     0  607.25
#> Median:          1  375    253      0        4     0 1213.50
#> 75%-tile:        1  375    253      0        5     0 1819.75
#> 97.5%-tile:      1  375    254      0        6     0 2365.38
#> Maximum:         1  375    256      0        6     0 2426.00
#> Mean:            1  375    252      0        4     0 1213.50
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 2425 
#> 
#> Total number of resource references: 1 
#> 

# You can also add references using the 'add_references' function.
```
