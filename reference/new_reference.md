# new_reference

Create a resource reference for your
[strollur](https://mothur.org/strollur/reference/strollur.html) object
to aid in reproducibility.

## Usage

``` r
new_reference(
  name,
  vendor = "",
  version = "",
  usage = "",
  note = "",
  documentation_url = "",
  method_url = "",
  parameter = "",
  citation = ""
)
```

## Arguments

- name, :

  a string containing the name of the resource used. For example:
  'silva.bacteria.fasta' or 'R package phylotypr'.

- vendor:

  a string containing name of entity that created original resource.
  example: "Silva" or "Schloss Lab - University of Michigan"

- version, :

  a string containing the version of the reference resource. For
  example: '1.38.1' or '0.1.1'. Default = "".

- usage, :

  a string containing the usage of the resource reference in your
  analysis. For example: 'alignment of sequences' or 'classification of
  sequences'. Default = "".

- note, :

  a string containing additional notes about the resource reference in
  your analysis. For example: 'alignment reference trimmed to V4 region'
  or 'classification of sequences using Bayesian method'. Default = "".

- documentation_url, :

  a string containing a web address where the reference may be
  downloaded or documentation may be found. Default = "".

- method_url, :

  a string containing any publications describing the methods used by
  the resource reference. For example: 'doi:10.1128/mra.01144-24'.
  Default = "".

- parameter, :

  a string containing the any specific parameters used by the resource.
  For example: 'kmer_size = 8, num_bootstraps = 100, min_confidence =
  80' Default = "".

- citation, :

  a string containing the citation information for the resource
  reference. For example: "citation_key = "doi:10.1128/AEM.00062-07",
  author = "Qiong Wang and George M. Garrity and James M. Tiedje and
  James R. Cole", title = "Naïve Bayesian Classifier for Rapid
  Assignment of rRNA Sequences into the New Bacterial Taxonomy", journal
  = "Applied and Environmental Microbiology", volume = "73", number =
  "16", pages = "5261-5267", year = "2007", doi =
  "10.1128/AEM.00062-07"". Default = "".

## Value

a list

## Examples

``` r

silva_resource <- new_reference( vendor = "SILVA", name =
"silva.bacteria.fasta", version = "1.38.1", usage = "alignment of sequences",
note = "alignment reference trimmed to V4 region", documentation_url =
"https://mothur.org/wiki/silva_reference_files/", method_url =
"https://mothur.org/blog/2024/SILVA-v138_2-reference-files/" )

phylotypr_resource <- new_reference( vendor = "Schloss Lab - University of
Michigan", name = "R phylotypr package", version = "0.1.1", usage =
"classification of sequences", note = "classification using Bayesian method",
parameter = "kmer_size = 8, num_bootstraps = 100, min_confidence = 80",
documentation_url = "https://mothur.org/phylotypr/", method_url =
"doi:10.1128/mra.01144-24", citation = "@article{doi:10.1128/AEM.00062-07,
author = {Qiong Wang and George M. Garrity and James M. Tiedje and James R.
Cole}, title = {Naïve Bayesian Classifier for Rapid Assignment of rRNA
Sequences into the New Bacterial Taxonomy}, journal = {Applied and
Environmental Microbiology}, volume = {73}, number = {16}, pages =
{5261-5267}, year = {2007}, doi = {10.1128/AEM.00062-07}, URL =
{https://journals.asm.org/doi/abs/10.1128/aem.00062-07}, eprint =
{https://journals.asm.org/doi/pdf/10.1128/aem.00062-07}}"
)
```
