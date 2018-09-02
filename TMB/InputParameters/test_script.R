setwd("C:/Work/Projects/PhD/Logistic_normal/TMB/InputParameters")
library(TMB)

## this script tests the functionality of the parameter configuration that I want to take from stockassessment.org
## and add it to future general programs that I will construct.
compile("test_model.cpp")
dyn.load(dynlib("model"))