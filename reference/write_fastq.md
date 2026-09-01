# Write a [FASTQ](https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/) formatted file

Write a [FASTQ](https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/)
formatted file

## Usage

``` r
write_fastq(data, filename = NULL)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- filename:

  a string containing the name of the output file. Default =
  'dataset_name'.fastq

## Value

name of FASTQ file

## Examples

``` r

table <- strollur::read_fastq(strollur_example("tiny.fastq.gz"))

# Create `strollur::strollur` object
data <- strollur::new_dataset("example")

# Add FASTQ data `strollur::strollur` object
strollur::add(data, table, type = "fastq")
#> Added 3 sequences.
#> Assigned 3 quality scores.

# Write `strollur::strollur` objects FASTQ data to file
strollur::write_fastq(data, tempfile())
#> [1] "/tmp/RtmpUVSXMn/file3efd38ec90d5"
```
