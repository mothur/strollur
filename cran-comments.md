## Resubmission
* Updates link to 'Getting_Started' in 'Data_Transfers' vignette to remove note.
    * checking CRAN incoming feasibility ... [4s/6s] NOTE
      - Found the following (possibly) invalid file URI. 
      - URI: Getting-Started.html From: inst/doc/Data_Transfers.html
* Removes unused childNumber variable from TaxNode in phylotree.h.
    * checking whether package ‘strollur’ can be installed ... [57s/57s] WARNING
      - Found the following significant warnings:
      - phylotree.h:8:8: warning: ‘<unnamed>.TaxNode::childNumber’ is used uninitialized [-Wuninitialized]
      
## Test environments
* local OS package environment: macOS Tahoe, R 4.6.0
* GitHub Actions: 
    * R CMD check of Ubuntu-latest (R-release, R-devel), Windows-latest, macOS-latest.
    * Test Coverage
    * Linting
    * Pkgdown build and deployment
    * Dependency only build and check
    
## R CMD check results
0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Notes for Reviewer
This is a new package. All automated checks on GitHub Actions pass cleanly.

