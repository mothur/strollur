# Get a data.frame containing the given report in a [strollur](https://mothur.org/strollur/reference/strollur.html) object

Get a data.frame containing the report. Reports include FASTA format,
sequences reports, sequence_bin_assignments, sequence_taxonomy,
bin_taxonomy, bin_representatives, sample_assignments, metadata,
references, sequence_scrap, and bin_scrap in a
[strollur](https://mothur.org/strollur/reference/strollur.html) object.

## Usage

``` r
report(data, type = "sequences", bin_type = "otu")
```

## Arguments

- data, :

  a [strollur](https://mothur.org/strollur/reference/strollur.html)
  object

- type, :

  string containing the type of report you would like. Options include:
  "fasta", "sequences", "sequence_bin_assignments", "sequence_taxonomy",
  "bin_taxonomy", "bin_representatives", "sample_assignments",
  "metadata", "references", "sequence_scrap", "bin_scrap". If you have
  added custom reports for alignment, contigs_assembly or chimeras, you
  can get those as well. Default = "sequences".

- bin_type, :

  string containing the bin type you would like a bin_taxonomy report
  for. Default = "otu".

## Value

data.frame

## Examples

``` r
miseq <- miseq_sop_example()
#> ℹ Added 2425 sequences.
#> ℹ Assigned 2425 sequence abundances.
#> ℹ Assigned 2425 sequence taxonomies.
#> ℹ Assigned 531 otu bins.
#> ℹ Assigned 2425 asv bins.
#> ℹ Assigned 63 phylotype bins.
#> ℹ Assigned 19 samples to treatments.
#> ℹ Assigned 531 otu bin taxonomies.
#> ℹ Assigned 531 otu bin representative sequences.
#> ℹ Added metadata.
#> ℹ Added 2 resource references.
#> ℹ Added a contigs_report.

# To get the FASTA data

report(data = miseq, type = "fasta") |> head(n = 5)
#>                                 sequence_names
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783
#> 2  M00967_43_000000000-A3JHG_1_1113_12711_3318
#> 3  M00967_43_000000000-A3JHG_1_2108_14707_9807
#> 4  M00967_43_000000000-A3JHG_1_1110_4126_16552
#> 5  M00967_43_000000000-A3JHG_1_2102_8408_13436
#>                                                                                                                                                                                                                                                                                                                                                                                 sequences
#> 1 TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGTGAT--TT-A-T-T--GG-GT--TT-A-AA-GG-GT-GC-G-TA-GGC-G-G-A-CA-G-T-T-AA-G-T-C-A-G-C-G-G--TA-A-AA-TT-G-A-GA-GG--CT-C-AA-C-C-T-C-T-T-C--CC-G-C-CGTT-GAAAC-TG-A-TTGTC-TTGA-GT-GG-GC-GA-G-A---AG-T-A-TGTGGAATGCGTGGTGT-AGCGGT-GAAATGCATAG-AT-A-TC-AC-GC-AG-AACTCCGAT-TGCGAAGGCA------GCATA-CCG-G-CG-CC-C-A-ACTGACG-CTGA-AGCA-CGAAA-GCG-TGGGT-ATC-GAACAGG
#> 2 TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GGC-G-G-C-CA-T-G-C-AA-G-T-C-A-G-A-A-G--TG-A-AA-AC-C-C-GG-GG--CT-C-AA-C---C-C-TGG-G-AGT-G-C-TTTT-GAAAC-TG-T-GCGGC-TAGA-GT-GT-CG-GA-G-G---GG-T-A-AGTGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCTTA-CTG-G-AC-GA-T-G-ACTGACG-CTGA-GGCT-CGAAA-GCG-TGGGT-ATC-GAACAGG
#> 3 TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-C-GG-C-G-C-AA-G-T-C-T-G-A-A-G--TG-A-AA-GC-C-C-GT-GG--CT-C-AA-C-C-G-C-G-G-A-ACC-G-C-TTTG-GAAAC-TG-C-GAGGC-TGGA-GT-GC-TG-GA-G-A---GG-T-A-AGCGGAATTCCTGGTGT-AGCGGT-GAAATGCGTAG-AT-A-TC-AG-GA-GG-AACACCGGT-GGCGAAGGCG------GCTTA-CTG-G-AC-AG-T-G-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-GAACAGG
#> 4 TAC--GG-AG-GAT--TCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-T-T--GG-GT--TT-A-AA-GG-GT-GC-G-TA-GGC-G-G-G-CT-G-T-T-AA-G-T-C-A-G-C-G-G--TC-A-AA-TG-T-C-GG-GG--CT-C-AA-C-C-C-C-G-G-C--CT-G-C-CGTT-GAAAC-TG-G-CGGCC-TCGA-GT-GG-GC-GA-G-A---AG-T-A-TGCGGAATGCGTGGTGT-AGCGGT-GAAATGCATAG-AT-A-TC-AC-GC-AG-AACTCCGAT-TGCGAAGGCA------GCATA-CCG-G-CG-CC-C-T-ACTGACG-CTGA-GGCA-CGAAA-GCG-TGGGT-ATC-GAACAGG
#> 5 TAC--GG-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GGC-G-G-C-AG-T-G-C-AA-G-T-C-A-G-A-A-G--TG-A-AA-GC-C-C-AA-GG--CT-C-AA-C---C-A-TGG-G-ACT-G-C-TTTT-GAAAC-TG-T-ACAGC-TAGA-TT-GC-AG-GA-G-A---GG-T-A-AGTGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCTTA-CTG-G-AC-TG-T-A-AATGACG-CTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG

# To get a report about the FASTA data

report(data = miseq, type = "sequences") |> head(n = 5)
#>                                             id start end length ambig
#> 1 M00967_43_000000000-A3JHG_1_2101_16474_12783     1 375    253     0
#> 2  M00967_43_000000000-A3JHG_1_1113_12711_3318     1 375    253     0
#> 3  M00967_43_000000000-A3JHG_1_2108_14707_9807     1 375    253     0
#> 4  M00967_43_000000000-A3JHG_1_1110_4126_16552     1 375    252     0
#> 5  M00967_43_000000000-A3JHG_1_2102_8408_13436     1 375    253     0
#>   longest_homopolymer num_n
#> 1                   4     0
#> 2                   5     0
#> 3                   4     0
#> 4                   4     0
#> 5                   5     0

# To get the sequence bin assignments

report(data = miseq, type = "sequence_bin_assignments", bin_type = "otu") |>
head(n = 5)
#>   otu_id                                       seq_id
#> 1 Otu001  M00967_43_000000000-A3JHG_1_1111_20933_6700
#> 2 Otu001  M00967_43_000000000-A3JHG_1_1113_17095_9759
#> 3 Otu001 M00967_43_000000000-A3JHG_1_1114_22144_24942
#> 4 Otu001   M00967_43_000000000-A3JHG_1_1112_5981_8948
#> 5 Otu001  M00967_43_000000000-A3JHG_1_2106_5509_18056

# To get the sample treatment assignments

report(data = miseq, type = "sample_assignments")
#>    samples treatments
#> 1     F3D0      Early
#> 2     F3D1      Early
#> 3   F3D141       Late
#> 4   F3D142       Late
#> 5   F3D143       Late
#> 6   F3D144       Late
#> 7   F3D145       Late
#> 8   F3D146       Late
#> 9   F3D147       Late
#> 10  F3D148       Late
#> 11  F3D149       Late
#> 12  F3D150       Late
#> 13    F3D2      Early
#> 14    F3D3      Early
#> 15    F3D5      Early
#> 16    F3D6      Early
#> 17    F3D7      Early
#> 18    F3D8      Early
#> 19    F3D9      Early

# To get a report about sequence classifications

report(data = miseq, type = "sequence_taxonomy") |> head(n = 10)
#>                                              id level
#> 1  M00967_43_000000000-A3JHG_1_2101_16474_12783     1
#> 2  M00967_43_000000000-A3JHG_1_2101_16474_12783     2
#> 3  M00967_43_000000000-A3JHG_1_2101_16474_12783     3
#> 4  M00967_43_000000000-A3JHG_1_2101_16474_12783     4
#> 5  M00967_43_000000000-A3JHG_1_2101_16474_12783     5
#> 6  M00967_43_000000000-A3JHG_1_2101_16474_12783     6
#> 7   M00967_43_000000000-A3JHG_1_1113_12711_3318     1
#> 8   M00967_43_000000000-A3JHG_1_1113_12711_3318     2
#> 9   M00967_43_000000000-A3JHG_1_1113_12711_3318     3
#> 10  M00967_43_000000000-A3JHG_1_1113_12711_3318     4
#>                                taxon confidence
#> 1                           Bacteria        100
#> 2                    "Bacteroidetes"        100
#> 3                      "Bacteroidia"         99
#> 4                    "Bacteroidales"         99
#> 5               "Porphyromonadaceae"         88
#> 6  "Porphyromonadaceae"_unclassified         88
#> 7                           Bacteria        100
#> 8                         Firmicutes        100
#> 9                         Clostridia        100
#> 10                     Clostridiales        100

# To get a report about bin classifications for 'otu' data

report(data = miseq, type = "bin_taxonomy", bin_type = "otu") |> head(n = 10)
#>        id level                             taxon confidence
#> 1  Otu001     1                          Bacteria        100
#> 2  Otu001     2                   "Bacteroidetes"        100
#> 3  Otu001     3                     "Bacteroidia"        100
#> 4  Otu001     4                   "Bacteroidales"        100
#> 5  Otu001     5              "Porphyromonadaceae"        100
#> 6  Otu001     6 "Porphyromonadaceae"_unclassified        100
#> 7  Otu002     1                          Bacteria        100
#> 8  Otu002     2                   "Bacteroidetes"        100
#> 9  Otu002     3                     "Bacteroidia"        100
#> 10 Otu002     4                   "Bacteroidales"        100

# To get the 'otu' bin representative sequences

report(
  data = miseq, type = "bin_representatives",
  bin_type = "otu"
) |> head(n = 5)
#>   otu_names                         representative_names
#> 1    Otu001 M00967_43_000000000-A3JHG_1_1108_14299_17220
#> 2    Otu002  M00967_43_000000000-A3JHG_1_1106_22705_6123
#> 3    Otu003  M00967_43_000000000-A3JHG_1_1101_15533_5293
#> 4    Otu004 M00967_43_000000000-A3JHG_1_1105_25642_17588
#> 5    Otu005  M00967_43_000000000-A3JHG_1_2102_7041_13746
#>                                                                                                                                                                                                                                                                                                                                                                  representative_sequences
#> 1 TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-C-TG-T-G-C-AA-G-T-C-T-G-A-A-G--TG-A-AA-TG-C-C-GG-GG--CT-C-AA-C-C-C-C-G-G-A-ACT-G-C-TTTG-GAAAC-TG-T-ACAGC-TAGA-GT-GC-AG-GA-G-G---GG-T-G-AGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCGGT-GGCGAAGGCG------GCTCA-CTG-G-AC-TG-T-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 2 TAC--GT-AG-GGG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-CA-GAC-G-G-C-TG-T-G-C-AA-G-T-C-T-G-G-A-G--TG-A-AA-GG-C-G-GG-GG--CC-C-AA-C-C-C-C-C-G-G-ACT-G-C-TCTG-GAAAC-TG-T-AAAGC-TGGA-GT-GC-AG-GA-G-A---GG-T-A-AGCGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCTTA-CTG-G-AC-TG-C-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGT-ATC-GAACAGG
#> 3 TAC--GG-AG-GAT--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GA-GC-G-TA-GAC-G-G-C-GA-T-G-C-AA-G-T-C-T-G-A-A-G--TG-A-AA-GG-C-G-GG-GG--CC-C-AA-C-C-C-C-C-G-G-ACT-G-C-TTTG-GAAAC-TG-T-ATAGC-TGGA-GT-GC-AG-GA-G-A---GG-T-A-AGTGGAATTCCTAGTGT-AGCGGT-GAAATGCGTAG-AT-A-TT-AG-GA-GG-AACACCAGT-GGCGAAGGCG------GCTTA-CTG-G-AC-TG-T-A-ACTGACG-TTGA-GGCT-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 4 TAC--GT-AG-GTG--GCA-A-G-C-G-T-T--AT-C-CGG-AT--TT-A-C-T--GG-GT--GT-A-AA-GG-GC-GT-G-TA-GGC-G-G-G-AC-T-G-C-AA-G-T-C-A-G-A-T-G--TG-A-AA-CC-C-A-TG-GG--CT-C-AA-C-C-C-A-T-G-G-CCT-G-C-ATTT-GAAAC-TG-T-AGTTC-TTGA-GT-GA-TG-GA-G-A---GG-C-A-GGCGGAATTCCGTGTGT-AGCGGT-GAAATGCGTAG-AT-A-TA-CG-GA-GG-AACACCAGT-GGCGAAGGCG------GCCTG-CTG-G-AC-AT-T-A-ACTGACG-CTGA-GGCG-CGAAA-GCG-TGGGG-AGC-AAACAGG
#> 5 TAC--GT-AG-GGG--GCG-A-G-C-G-T-T--AT-C-CGG-AT--TC-A-T-T--GG-GC--GT-A-AA-GC-GC-GC-G-CA-GGC-G-G-A-CT-C-A-T-AA-G-C-G-G-A-G-C-C--TT-T-AA-TC-T-T-GG-GG--CT-T-AA-C-C-T-C-A-A-G-T-C-G-G-GCCC-CGAAC-TG-T-GAGTC-TCGA-GT-GT-GG-TA-G-G---GG-A-A-GGCGGAATTCCCGGTGT-AGCGGT-GGAATGCGCAG-AT-A-TC-GG-GA-AG-AACACCGAT-GGCGAAGGCA------GCCTT-CTG-G-GC-CA-T-C-ACTGACG-CTGA-GGCG-CGAAA-GCT-AGGGG-AGC-AAACAGG

# To get a report about the sequences removed during your analysis:

report(data = miseq, type = "sequence_scrap")
#> data frame with 0 columns and 0 rows

# To get a report about the "otu" bins removed during your analysis:

report(data = miseq, type = "bin_scrap", bin_type = "otu")
#> data frame with 0 columns and 0 rows

# To get the metadata associated with your data:

metadata <- report(data = miseq, type = "metadata")

# To get the resource references associated with your data:

references <- report(data = miseq, type = "references")

# To get our custom report containing the contigs assembly data:

report(data = miseq, type = "contigs_report") |> head(n = 10)
#>                                            Name Length Overlap_Length
#> 1  M00967_43_000000000-A3JHG_1_2101_16474_12783    253            250
#> 2   M00967_43_000000000-A3JHG_1_1113_12711_3318    253            249
#> 3   M00967_43_000000000-A3JHG_1_2108_14707_9807    253            249
#> 4   M00967_43_000000000-A3JHG_1_1110_4126_16552    252            249
#> 5   M00967_43_000000000-A3JHG_1_2102_8408_13436    253            249
#> 6  M00967_43_000000000-A3JHG_1_1107_22580_21773    252            250
#> 7  M00967_43_000000000-A3JHG_1_1108_14299_17220    253            249
#> 8   M00967_43_000000000-A3JHG_1_1114_8059_18290    253            249
#> 9    M00967_43_000000000-A3JHG_1_2112_9811_9982    253            249
#> 10  M00967_43_000000000-A3JHG_1_2103_25452_6018    252            249
#>    Overlap_Start Overlap_End MisMatches Num_Ns Expected_Errors
#> 1              2         252         19      0      0.29461400
#> 2              2         251          0      0      0.00183396
#> 3              2         251          0      0      0.00196774
#> 4              2         251          4      0      0.05629750
#> 5              2         251          0      0      0.00259554
#> 6              1         251          8      0      0.05068300
#> 7              2         251          0      0      0.00215398
#> 8              3         252         11      0      0.05362360
#> 9              2         251          0      0      0.00184906
#> 10             2         251         18      0      0.54816100
```
