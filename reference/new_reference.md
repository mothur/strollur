# new_reference

Create a reference you can add to your dataset

## Usage

``` r
new_reference(
  reference_name,
  reference_version = "",
  reference_usage = "",
  reference_note = "",
  reference_url = ""
)
```

## Arguments

- reference_name, :

  a string containing the name of the reference used in the preparing of
  the sequences. For example: 'silva.bacteria.fasta'.

- reference_version, :

  a string containing the version of the reference used in the preparing
  of the sequences. For example: '1.38.1'. Default = "".

- reference_usage, :

  a string containing the usage of the reference in your analysis. For
  example: 'alignment using mothur2' Default = NULL.

- reference_note, :

  a string containing the any additional notes about the reference.
  Default = "".

- reference_url, :

  a string containing a web address where the reference may be
  downloaded. Default = "".

## Value

a list

## See also

\[add()\]

## Examples

``` r
reference <- new_reference("silva.bacteria.fasta",
                           "1.38.1",
                           "alignment by mothur2 v1.0 using default options",
                           "",
                           "https://mothur.org/wiki/silva_reference_files/")
```
