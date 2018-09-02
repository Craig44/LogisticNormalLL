## Add Libraries
library(TMB)
#library(TMBdebug) # to stop windows from crashing
library(casal)
## Add other dependency funcitons
setwd("C:/Work/Projects/PhD/Logistic_normal/Rstudio");
source("Initialisation.r");

## Bring some data to debug this model
load(file = make.filename(file = "tmb.RData", path = DIR$TMB))

## add the new switches
data$use_logistic_normal = 0;
pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma = 0.6, log_norm_phi = 0.4)

setwd(DIR$TMB)
setwd("AddLogisticNormal")
compile("model.cpp", "-O1 -g3", DLLFLAGS="")
dyn.load(dynlib("model"))


## Bring in the HAK data to test
Hak  = extract.fits(path = DIR$'csl_files', file = "HAK.log");
compdat = Hak$subaTANageDEC
## restructure for Chris's Code
compdat$obs =  as.matrix(Hak$subaTANageDEC$obs)
compdat$exp =  round(as.matrix(Hak$subaTANageDEC$fits),5)
compdat$N = Hak$subaTANageDEC$error.value[,1]

data$fishery_at_age_obs = compdat$obs;
data$fishery_at_age_exp_LN_test = compdat$exp;
data$fishery_at_age_error_LN_test = compdat$N;
data$deviations = 0
YCS_start = rep(1,length(data$years))
pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma = 0.6, log_norm_phi = 0.4)

obj_ran <- MakeADFun(data=data,parameters=pars, random = "YCS", DLL = "model",hessian = TRUE)

