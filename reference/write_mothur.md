# Write the various [file types](https://mothur.org/wiki/tags/#file_types) from your [strollur object](https://mothur.org/strollur/reference/strollur.html) for use with mothur.

The write_mothur function will write various [file
types](https://mothur.org/wiki/tags/#file_types) for use with mothur.

## Usage

``` r
write_mothur(data, dir_path = NULL, compress = TRUE, tags = NULL)
```

## Arguments

- data:

  a [strollur
  object](https://mothur.org/strollur/reference/strollur.html)

- dir_path:

  a string containing the name of directory where the files should be
  written. Default = current working directory.

- compress:

  logical, Default = TRUE.

- tags:

  a vector of strings containing the items you wish to write Options are
  'sequence_data', 'bin_data', 'metadata', 'resource_reference',
  'sequence_tree', 'sample_tree' and 'report'. By default, everything is
  written to files.

## Value

a vector of file names

## References

Schloss,P.D., Westcott,S.L., Ryabin,T., Hall,J.R., Hartmann,M.,
Hollister,E.B., Lesniewski,R.A., Oakley,B.B., Parks,D.H., Robinson,C.J.,
Sahl,J.W., Stres,B., Thallinger,G.G., Van Horn,D.J. and Weber,C.F.
(2009), Introducing mothur: Open-source, platform-independent,
community-supported software for describing and comparing microbial
communities. Applied and Environmental Microbiology 75:7537-7541.
<doi:10.1128/AEM.01541-09>

## Examples

``` r

miseq <- strollur::miseq_sop_example()
#> Added 2425 sequences.
#> Assigned 2425 sequence abundances.
#> Assigned 2425 sequence taxonomies.
#> Assigned 531 otu bins.
#> Assigned 2425 asv bins.
#> Assigned 63 phylotype bins.
#> Assigned 19 samples to treatments.
#> Assigned 171 samples distances.
#> Assigned 531 otu bin taxonomies.
#> Assigned 531 otu bin representative sequences.
#> Added a metadata report.
#> Added 2 resource references.
#> Added a contigs_report report.
files <- strollur::write_mothur(miseq, tempdir(), compress = FALSE)
```
