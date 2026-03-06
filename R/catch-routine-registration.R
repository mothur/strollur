# This dummy function definition is included with the package to ensure that
# 'tools::package_native_routine_registration_skeleton()' generates the required
# registration info for the 'run_testthat_tests' symbol.

(function() {
  #' @name run_testthat_tests
  #' @title run_testthat_tests
  #' @description This dummy function definition is included with the package
  #'  to ensure that 'tools::package_native_routine_registration_skeleton()'
  #'  generates the required registration info for the 'run_testthat_tests'
  #'  symbol. It allows for testing using c++ unit tests.
  #'
  #' @keywords internal
  .Call("run_testthat_tests", FALSE, PACKAGE = "strollur")
})
