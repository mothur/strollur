# read_mothur_count

Read a mothur formatted [count
file](https://mothur.org/wiki/count_file/)

## Usage

``` r
read_mothur_count(filename)
```

## Arguments

- filename:

  count file name (required)

## Value

data.frame

## Examples

``` r
# mothur count file
# Representative_Sequence     total   sample2  sample3  sample4
# seq1  1150  250  400  500
# seq2  115  25  40  50
# seq3  50  25  25  0
# seq4  4  0  0  4

# returns
# id   sample abundance
# <char>  <char>     <int>
#  1:   seq1 sample2       250
#  2:   seq1 sample3       400
#  3:   seq1 sample4       500
#  4:   seq2 sample2        25
#  5:   seq2 sample3        40
#  6:   seq2 sample4        50
#  7:   seq3 sample2        25
#  8:   seq3 sample3        25
#  9:   seq4 sample4         4

# read a count file with samples
sample_table <- read_mothur_count(strollur_example("final.count_table.gz"))

# You can add your sequence abundance data to your dataset as follows:

# create a new empty dataset
data <- new_dataset()

# assign sequence abundances parsed by sample
assign(data = data, table = sample_table, type = "sequence_abundance")
#> ℹ Assigned 2425 sequence abundances.
#> [1] 2425

# print summary of data
data
#> 
#> Sample   Total:
#> F3D0 6191 
#> F3D1 4652 
#> F3D141   4656 
#> F3D142   2423 
#> F3D143   2403 
#> F3D144   3449 
#> F3D145   5532 
#> F3D146   3831 
#> F3D147   12430 
#> F3D148   9465 
#> F3D149   10014 
#> F3D150   4126 
#> F3D2 15686 
#> F3D3 5199 
#> F3D5 3469 
#> F3D6 6394 
#> F3D7 4055 
#> F3D8 4253 
#> F3D9 5735 
#> 
#> Number of unique seqs: 2425 
#> Total number of seqs: 113963 
#> 
```
