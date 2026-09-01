# Write a file containing sequence quality scores

Write a file containing sequence quality scores

## Usage

``` r
write_quality(data, filename = NULL)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- filename:

  a string containing the name of the output file. Default =
  'dataset_name'.qual

## Value

name of sequence quality file

## Examples

``` r

table <- strollur::read_fastq(strollur_example("tiny.fastq.gz"))

# Create `strollur::strollur` object
data <- strollur::new_dataset("example")

# Add FASTQ data `strollur::strollur` object
strollur::add(data, table, type = "fastq")
#> Added 3 sequences.
#> Assigned 3 quality scores.

# Write `strollur::strollur` objects sequence quality data to file
strollur::write_quality(data, tempfile())
#> [1] "/tmp/RtmptQFCni/file1d8b287376ae"
```
