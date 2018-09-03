## (03) Compare_comp_likelihood.R
## @Author C Marsh
## @Date 3/9/2018
## @Description
## This script uses actual observations from the HAK assessment and compares fits in a 'fixed' effects model to look at residuals to convince me that
## all this trouble was worth it.

## Add Libraries
library(TMBdebug) # to stop windows from crashing
library(casal)

## Add other dependency funcitons
setwd("C:/Work/Projects/PhD/Logistic_normal/Rstudio");
source("Initialisation.r");

Hak  = extract.fits(path = DIR$'csl_files', file = "HAK.log");





