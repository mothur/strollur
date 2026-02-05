# read_mothur_list

Read a mothur formatted [list file](https://mothur.org/wiki/list_file/)

## Usage

``` r
read_mothur_list(list)
```

## Arguments

- list:

  file name. The [list file](https://mothur.org/wiki/list_file/) can be
  created using several of mothur's commands.
  [cluster](https://mothur.org/wiki/cluster/),
  [cluster.split](https://mothur.org/wiki/cluster.split/),
  [cluster.fit](https://mothur.org/wiki/cluster.fit/) and
  [phylotype](https://mothur.org/wiki/phylotype/).

## Value

A data.frame containing the sequence otu assignments

## Examples

``` r
# You can add your otu assignments to the your data set using the following:

# read mothur's list file into data.frame
otu_data <- read_mothur_list(strollur_example("final.opti_mcc.list.gz"))

# create a new empty dataset
data <- new_dataset()

# assign sequences to 'otu' bins
assign(data = data, table = otu_data, type = "bins", bin_type = "otu")
#> ℹ Assigned 531 otu bins.
#> [1] 531
```
