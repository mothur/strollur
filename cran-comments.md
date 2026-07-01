## Test environments

- local OS package environment: macOS Tahoe, R 4.6.0
- GitHub Actions:
  - R CMD check of Ubuntu-latest (R-release, R-devel), Windows-latest, macOS-latest.
  - Test Coverage
  - Linting
  - Pkgdown build and deployment
  - Dependency only build and check
  - rhub::rhub_check() - clang-asan, gcc-asan
    * gcc-asan: PASS
    * clang-asan: 1 warning (ODR-violation in upstream dependency RcppParallel)

## R CMD check results

0 errors ✔ \| 0 warnings ✔ \| 0 notes ✔

## Notes for Reviewer

- All automated checks on GitHub Actions pass cleanly. 
- Version 0.1.3 corrects a `clang-ASAN` and `gcc-ASAN` issue detected in `strollur` version 0.1.2.
    - The clang-asan build notes a duplicate symbol (`__itt_detach_ptr__3_0`) present simultaneously in `libtbbmalloc.so.2` and `libtbb.so.2`. This is a known issue tracked under RcppParallel (Issue #152). This is outside the scope of the strollur package codebase, and all package-specific memory leaks and heap overflows have been resolved.

## Version - 0.1.3

Initial release of v0.1.3

## Resubmission - 0.1.2

- Updates Accessing_Data and General_Importing vignettes replacing on.exit(par(...)) with par(...).
  - Please always make sure to reset to user's options(), working directory or par() after you changed it in examples and vignettes and demos.

## Resubmission - 0.1.2

- Removes references from description in the DESCRIPTION file and adds references in the read / write functions for the software packages.
  - Please write references in the description of the DESCRIPTION file in the form authors (year) \<doi:...\> authors (year, ISBN:...) or if those are not available: authors (year) \<https:...\> with no space after 'doi:', 'https:' and angle brackets for auto-linking.
- Modifies get_full_name() function in helper.R for test outputs. Previously wrote to the testthat::test_path(), now writes to tempdir(). Also, updates other tests to write to tempdir(), instead of test_path().
  - Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace.
- Restores par() values in Accessing_Data and General_Importing vignettes.
  - Please always make sure to reset to user's options(), working directory or par() after you changed it in examples and vignettes and demos.
- Updates strollur version from 0.1.1 to 0.1.2.
- Updates maintainer to Pat Schloss [pschloss\@umich.edu](mailto:pschloss@umich.edu){.email}

## Resubmission - 0.1.1

- Updates link to 'Getting_Started' in 'Data_Transfers' vignette to remove note.
  - checking CRAN incoming feasibility ... [4s/6s] NOTE
    - Found the following (possibly) invalid file URI.
    - URI: Getting-Started.html From: inst/doc/Data_Transfers.html
- Removes unused childNumber variable from TaxNode in phylotree.h.
  - checking whether package ‘strollur’ can be installed ... [57s/57s] WARNING
    - Found the following significant warnings:
    - phylotree.h:8:8: warning: ‘<unnamed>.TaxNode::childNumber’ is used uninitialized [-Wuninitialized]
- Updates strollur version from 0.1.0 to 0.1.1.

