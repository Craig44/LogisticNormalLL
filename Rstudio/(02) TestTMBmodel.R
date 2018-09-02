## (02) TestTMBmodel.R
## @Author C Marsh
## @Date 2/9/2018
## @Description
## Compare my TMB model against the original script for unsexed comp data, as my TMB model is unsexed at them moment


## Add Libraries
library(TMB)
library(TMBdebug) # to stop windows from crashing
library(casal)

## Add other dependency funcitons
setwd("C:/Work/Projects/PhD/Logistic_normal/Rstudio");
source("Chris_original_code.r")
source("Initialisation.r");


Hak  = extract.fits(path = DIR$'csl_files', file = "HAK.log");


#########
## HAK
#########
compdat = Hak$subaTANageDEC
## restructure for Chris's Code
compdat$obs =  as.matrix(Hak$subaTANageDEC$obs[1:15,])
compdat$exp =  round(as.matrix(Hak$subaTANageDEC$fits[1:15,]),5)
compdat$N = Hak$subaTANageDEC$error.value[1:15,1]

# save this information to text files to compare with C++ code
write.table(x = compdat$obs, file = make.filename(file = "observed_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = compdat$exp, file = make.filename(file = "expected_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = compdat$N, file = make.filename(file = "error_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)

#### ------------------- C++ section ---------------
## run the program
setwd(DIR$`C++_work`);
if (file.exists("results.txt")) file.remove("results.txt")
system("LogNorm.exe notsexed 1",ignore.stdout = TRUE,ignore.stderr = TRUE)
setwd(DIR$R)
C_results = as.vector(read.table(make.filename(file = "results.txt", path = DIR$`C++_work`), header = F))

#### ------------------- Read in TMB model and load data needed --------
setwd(DIR$TMB)
setwd("AddLogisticNormal")
compile("orig.cpp")
dyn.load(dynlib("orig"))

load(file = make.filename(file = "tmb.RData", path = DIR$TMB))

## add the new switches
data$survey_at_age_obs = Obs_data$obs;
data$survey_at_age_error = Obs_data$N;
data$survey_years = as.numeric(rownames(data$survey_at_age_obs))

data$ages = 1:19
data$ageing_error = matrix(0,length(data$ages),length(data$ages))
diag(data$ageing_error) = 1;
data$fishery_at_age_obs = compdat$obs;
data$fishery_years = as.numeric(rownames(compdat$obs))
data$fishery_at_age_exp_LN_test = compdat$exp;
data$fishery_at_age_error_LN_test = compdat$N;
data$use_logistic_normal = 1;
unit_YCS = rep(1,length(data$years))
YCS_start = rep(0,length(data$years) - 1)
pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma = 0.1, log_norm_phi = c(-0.5,-0.5))


### ------------------- Run through R and TMB code to test that it does what we hope
sepbysex=F;
sexlag=F;
robust=F;
ARMA=F
sigma = 0.9;
phi = c(0.2)

pars$log_norm_sigma = log(sigma);
pars$log_norm_phi = log(phi);
data$ARMA = 0;
obj <- MakeADFun(data=data,parameters=pars, DLL = "orig")
obj_run = obj$report()
true_results = vector();
tmb_results = vector();
true_results[1] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
tmb_results[1] = obj_run$neg_ll_fishery_age

sigma = 0.2;
pars$log_norm_sigma = log(sigma);
obj <- MakeADFun(data=data,parameters=pars, DLL = "orig")
obj_run = obj$report()
true_results[2] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
tmb_results[2] = obj_run$neg_ll_fishery_age

sigma = 0.4;
phi = c(0.2,0.5)
pars$log_norm_sigma = log(sigma);
pars$log_norm_phi = log(phi);
obj <- MakeADFun(data=data,parameters=pars, DLL = "orig")
obj_run = obj$report()
true_results[3] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
tmb_results[3] = obj_run$neg_ll_fishery_age

sigma = 0.6;
phi = c(0.4,0.1)
pars$log_norm_sigma = log(sigma);
pars$log_norm_phi = log(phi);
obj <- MakeADFun(data=data,parameters=pars, DLL = "orig")
obj_run = obj$report()
true_results[4] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
tmb_results[4] = obj_run$neg_ll_fishery_age

ARMA = TRUE
data$ARMA = 0;
obj <- MakeADFun(data=data,parameters=pars, DLL = "orig")
obj_run = obj$report()
true_results[5] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
tmb_results[5] = obj_run$neg_ll_fishery_age

true_results - C_results
cbind(true_results, t(C_results),tmb_results,true_results -t(C_results))
