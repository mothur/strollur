# require_namespace

a wrapper for \`requireNamespace\`. Allowing us to more easily create
mock test for this functionality.

## Usage

``` r
require_namespace(package_name)
```

## Arguments

- package_name:

  the name of the package.

## Value

A logical TRUE or FALSE boolean depending on whether the package was
added to the search path.
