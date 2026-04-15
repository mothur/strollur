# read_mothur_rabund

Read a mothur formatted [rabund
file](https://mothur.org/wiki/rabund_file/)

## Usage

``` r
read_mothur_rabund(rabund)
```

## Arguments

- rabund:

  file name (required)

## Value

A data.frame containing the sequence otu assignments

## Examples

``` r
# You can add your otu assignments to the your data set using the following:

# read rabund file into data.frame
otu_data <- read_mothur_rabund(
  rabund =
    strollur_example("final.opti_mcc.rabund")
)

data <- new_dataset()

# assign abundance only 'otu' bins
assign(data = data, table = otu_data, type = "bins", bin_type = "otu")
#> Assigned 531 otu bins.
#> [1] 531
```
