# Script that runs all test cases.

rm(list = ls())
devtools::load_all()

source("tests/probit_tests.R")
source("tests/data_tests.R")
source("tests/prediction_tests.R")
source("tests/eval_tests.R")
source("tests/persistence_tests.R")
source("tests/pt_tests.R")
source("tests/tune_tests.R")


run_all_probit_tests()
run_all_pt_tests()
run_all_tune_tests()
run_all_data_tests()
run_all_prediction_tests()
run_all_eval_tests()
#run_all_persistence_tests()


cat("All tests passed.")