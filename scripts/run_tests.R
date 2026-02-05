# Script that runs all test cases.

source("tests/probit_tests.R")
source("tests/data_tests.R")
source("tests/prediction_tests.R")
source("tests/eval_tests.R")
source("tests/persistence_tests.R")

run_all_data_tests()
run_all_persistence_tests()
run_all_probit_tests()
run_all_prediction_tests()
run_all_eval_tests()

cat("All tests passed.")