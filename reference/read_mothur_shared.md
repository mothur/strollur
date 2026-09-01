# Read a mothur formatted [shared file](https://mothur.org/wiki/shared_file/)

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
otu_data <-
     strollur::read_mothur_shared(strollur_example("final.opti_mcc.shared"))

# create a new empty `strollur::strollur` object
data <- strollur::new_dataset()

# assign abundance only 'otu' bins parsed by sample
strollur::assign(data = data, table = otu_data,
                  type = "bin", bin_type = "otu")
#> Assigned 531 otu bins.
```
