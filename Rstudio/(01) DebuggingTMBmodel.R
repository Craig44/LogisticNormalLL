## (01) DebuggingTMBmodel.R
## @Author C Marsh
## @Date 29/08/2018
## @Description
## This script is for debugging the TMB models while I try to translate my C++ code -> TMB .hpp files for use in 
## the assessments I want to explore
## This shouldn't really be used after I am happt the model performs as expected.

## Add Libraries
library(TMB)
library(TMBdebug) # to stop windows from crashing
library(casal)
## Add other dependency funcitons
setwd("C:/Work/Projects/PhD/Logistic_normal/Rstudio");
source("Initialisation.r");

## Bring some data to debug this model
load(file = make.filename(file = "tmb.RData", path = DIR$TMB))


setwd(DIR$TMB)
setwd("AddLogisticNormal")
compile("orig.cpp")
dyn.load(dynlib("orig"))


## Bring in the HAK data to test
Hak  = extract.fits(path = DIR$'csl_files', file = "HAK.log");
compdat = Hak$subaTANageDEC
## restructure for Chris's Code
compdat$obs =  as.matrix(Hak$subaTANageDEC$obs)
compdat$exp =  round(as.matrix(Hak$subaTANageDEC$fits),5)
compdat$N = Hak$subaTANageDEC$error.value[,1]

Obs_data = Hak$subaOBSage
Obs_data$obs = as.matrix(Hak$subaOBSage$obs)
Obs_data$exp = round(as.matrix(Hak$subaOBSage$fits),5)
Obs_data$N = Hak$subaOBSage$error.value[,1]

## add the new switches
data$survey_at_age_obs = Obs_data$obs[1:15,];
data$survey_at_age_error = Obs_data$N[1:15];
data$survey_years = as.numeric(rownames(data$survey_at_age_obs))

data$ages = 1:19
data$ageing_error = matrix(0,length(data$ages),length(data$ages))
diag(data$ageing_error) = 1
data$fishery_at_age_obs = compdat$obs;
data$fishery_years = as.numeric(rownames(compdat$obs))
data$fishery_at_age_exp_LN_test = compdat$exp;
data$fishery_at_age_error_LN_test = compdat$N;
data$use_logistic_normal = 1;
data$ARMA = 1;

unit_YCS = rep(1,length(data$years))
#YCS_log_simp = Simplex_transform(unit_YCS)$transformed_vars
YCS_start = rep(0,length(data$years) - 1)
pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma = 0.1, log_norm_phi = c(-0.5,-0.5))

#gdbsource("model.R", TRUE)
compile("orig.cpp")
dyn.load(dynlib("orig"))
obj_ran <- MakeADFun(data=data,parameters=pars, DLL = "orig")

output = obj_ran$report()
output$covar
output$log_det_please
output$log_det
output$V_mat_inv
solve(output$V_mat)
log(det(output$V_mat))
output$neg_ll_LN
output$log_obs_tot


output$neg_ll_survey_age
output$neg_ll_fishery_age


