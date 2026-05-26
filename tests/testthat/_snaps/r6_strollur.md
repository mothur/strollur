# dataset - intialize from read_mothur / print

    Code
      waldo::compare(dataset_t$print(), dataset_t$print())
    Output
      miseq_sop:
      
                  starts ends   nbases ambigs polymers numns numseqs
      Minimum:         1  375 249.0000      0 3.000000     0     1.0
      2.5%-tile:       1  375 252.0000      0 4.000000     0  2436.7
      25%-tile:        1  375 252.0000      0 4.000000     0 24358.0
      Median:          1  375 253.0000      0 4.000000     0 48715.0
      75%-tile:        1  375 253.0000      0 5.000000     0 73072.0
      97.5%-tile:      1  375 254.0000      0 6.000000     0 94993.3
      Maximum:         1  375 256.0000      0 6.000000     0 97429.0
      Mean:            1  375 252.7296      0 4.502876     0 48715.0
          Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
      1 250.0000       232.0000      0.000000    248.0000   0.000000      0
      2 252.0000       247.0000      1.000000    250.0000   0.000000      0
      3 252.0000       249.0000      2.000000    251.0000   0.000000      0
      4 253.0000       249.0000      2.000000    251.0000   1.000000      0
      5 253.0000       250.0000      2.000000    251.0000   5.000000      0
      6 254.0000       251.0000      4.000000    253.0000  26.000000      0
      7 270.0000       255.0000     22.000000    256.0000 120.000000      0
      8 252.7488       249.1649      1.995206    251.1601   5.209012      0
        Expected_Errors
      1       1.0000000
      2       1.0000000
      3       1.0000000
      4       1.0000000
      5       1.0000000
      6       1.0000000
      7       4.0000000
      8       0.0740319
              type     trash_code unique total
      1   sequence    contaminant     38    39
      2   sequence remove_samples    200   211
      3   sequence           test    101  6052
      4        otu    contaminant      7     8
      5        otu remove_samples     30    31
      6        otu           test     19  6052
      7        asv    contaminant     38    39
      8        asv remove_samples    200   211
      9        asv           test    101  6052
      10 phylotype remove_samples      4     4
      11 phylotype           test      2  6052
                  starts ends nbases ambigs polymers numns  numseqs
      Minimum:         1  375    249      0        3     0     1.00
      2.5%-tile:       1  375    252      0        4     0  2436.70
      25%-tile:        1  375    252      0        4     0 24358.00
      Median:          1  375    253      0        4     0 48715.00
      75%-tile:        1  375    253      0        5     0 73072.00
      97.5%-tile:      1  375    254      0        6     0 94993.30
      Maximum:         1  375    256      0        6     0 97429.00
      Mean:            1  375    252      0        4     0 48715.00
      scrap_summary:
              type     trash_code unique total
      1   sequence    contaminant     38    39
      2   sequence remove_samples    200   211
      3   sequence           test    101  6052
      4        otu    contaminant      7     8
      5        otu remove_samples     30    31
      6        otu           test     19  6052
      7        asv    contaminant     38    39
      8        asv remove_samples    200   211
      9        asv           test    101  6052
      10 phylotype remove_samples      4     4
      11 phylotype           test      2  6052
      
      Number of unique seqs: 2086 
      Total number of seqs: 97428 
      
      Total number of samples: 17 
      Total number of treatments: 2 
      Total number of otus: 475 
      Total number of otu bin classifications: 475 
      Total number of asvs: 2086 
      Total number of asv bin classifications: 2086 
      Total number of phylotypes: 57 
      Total number of phylotype bin classifications: 57 
      Total number of sequence classifications: 2086 
      Total number of resource references: 2 
      Total number of custom reports: 1 
      Your dataset includes metadata 
      
      miseq_sop:
      
                  starts ends   nbases ambigs polymers numns numseqs
      Minimum:         1  375 249.0000      0 3.000000     0     1.0
      2.5%-tile:       1  375 252.0000      0 4.000000     0  2436.7
      25%-tile:        1  375 252.0000      0 4.000000     0 24358.0
      Median:          1  375 253.0000      0 4.000000     0 48715.0
      75%-tile:        1  375 253.0000      0 5.000000     0 73072.0
      97.5%-tile:      1  375 254.0000      0 6.000000     0 94993.3
      Maximum:         1  375 256.0000      0 6.000000     0 97429.0
      Mean:            1  375 252.7296      0 4.502876     0 48715.0
          Length Overlap_Length Overlap_Start Overlap_End MisMatches Num_Ns
      1 250.0000       232.0000      0.000000    248.0000   0.000000      0
      2 252.0000       247.0000      1.000000    250.0000   0.000000      0
      3 252.0000       249.0000      2.000000    251.0000   0.000000      0
      4 253.0000       249.0000      2.000000    251.0000   1.000000      0
      5 253.0000       250.0000      2.000000    251.0000   5.000000      0
      6 254.0000       251.0000      4.000000    253.0000  26.000000      0
      7 270.0000       255.0000     22.000000    256.0000 120.000000      0
      8 252.7488       249.1649      1.995206    251.1601   5.209012      0
        Expected_Errors
      1       1.0000000
      2       1.0000000
      3       1.0000000
      4       1.0000000
      5       1.0000000
      6       1.0000000
      7       4.0000000
      8       0.0740319
              type     trash_code unique total
      1   sequence    contaminant     38    39
      2   sequence remove_samples    200   211
      3   sequence           test    101  6052
      4        otu    contaminant      7     8
      5        otu remove_samples     30    31
      6        otu           test     19  6052
      7        asv    contaminant     38    39
      8        asv remove_samples    200   211
      9        asv           test    101  6052
      10 phylotype remove_samples      4     4
      11 phylotype           test      2  6052
                  starts ends nbases ambigs polymers numns  numseqs
      Minimum:         1  375    249      0        3     0     1.00
      2.5%-tile:       1  375    252      0        4     0  2436.70
      25%-tile:        1  375    252      0        4     0 24358.00
      Median:          1  375    253      0        4     0 48715.00
      75%-tile:        1  375    253      0        5     0 73072.00
      97.5%-tile:      1  375    254      0        6     0 94993.30
      Maximum:         1  375    256      0        6     0 97429.00
      Mean:            1  375    252      0        4     0 48715.00
      scrap_summary:
              type     trash_code unique total
      1   sequence    contaminant     38    39
      2   sequence remove_samples    200   211
      3   sequence           test    101  6052
      4        otu    contaminant      7     8
      5        otu remove_samples     30    31
      6        otu           test     19  6052
      7        asv    contaminant     38    39
      8        asv remove_samples    200   211
      9        asv           test    101  6052
      10 phylotype remove_samples      4     4
      11 phylotype           test      2  6052
      
      Number of unique seqs: 2086 
      Total number of seqs: 97428 
      
      Total number of samples: 17 
      Total number of treatments: 2 
      Total number of otus: 475 
      Total number of otu bin classifications: 475 
      Total number of asvs: 2086 
      Total number of asv bin classifications: 2086 
      Total number of phylotypes: 57 
      Total number of phylotype bin classifications: 57 
      Total number of sequence classifications: 2086 
      Total number of resource references: 2 
      Total number of custom reports: 1 
      Your dataset includes metadata 
      
      v No differences

