# read_mothur_shared

Read a mothur formatted [shared
file](https://mothur.org/wiki/shared_file/)

## Usage

``` r
read_mothur_shared(shared)
```

## Arguments

- shared:

  file name (required)

## Value

A data.frame containing the sequence otu assignments

## Examples

``` r

# You can add your otu assignments to the your data set using the following:

# read mothur shared file into data.frame
otu_data <- read_mothur_shared(strollur_example("final.opti_mcc.shared"))

# create a new empty `strollur` object
data <- new_dataset()

# assign abundance only 'otu' bins parsed by sample
assign(data = data, table = otu_data, type = "bin", bin_type = "otu")
#> Assigned 531 otu bins.
#> 
#> Number of unique seqs: 531 
#> Total number of seqs: 113963 
#> 
#> Total number of samples: 19 
#> Total number of otus: 531 
#> 
```
