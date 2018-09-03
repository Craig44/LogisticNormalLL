## (03) Compare_comp_likelihood.R
## @Author C Marsh
## @Date 3/9/2018
## @Description
## This script uses actual observations from the HAK assessment and compares fits in a 'fixed' effects model to look at residuals to convince me that
## all this trouble was worth it.

## Add Libraries
library(TMBdebug) # to stop windows from crashing
library(TMBhelper) 
library(casal)

## Add other dependency funcitons
setwd("C:/Work/Projects/PhD/Logistic_normal/Rstudio");
source("Initialisation.r");
source("Chris_original_code.r")

Hak  = extract.fits(path = DIR$'csl_files', file = "HAK.log");

survey_compdat = Hak$subaTANageDEC
survey_bio = Hak$subaTANbiomassDEC
length(survey_bio$year); length(survey_compdat$year)
surv_ndx = survey_compdat$year %in% survey_bio$year

fishery_compdat = Hak$subaTANageDEC

load(file = make.filename(file = "tmb.RData", path = DIR$TMB))
## add the new switches
data$survey_at_age_obs = as.matrix(survey_compdat$obs[surv_ndx,]);
data$survey_at_age_error = survey_compdat$error.value[surv_ndx,1];
data$survey_years = survey_compdat$year[surv_ndx]

data$survey_biomass_obs = survey_bio$obs
data$survey_biomass_error = survey_bio$error.value


data$ages = 1:19
data$years = 1989:2013
data$ageing_error = matrix(0,length(data$ages),length(data$ages))
diag(data$ageing_error) = 1;
data$fishery_at_age_obs = as.matrix(fishery_compdat$obs);
data$fishery_at_age_error = fishery_compdat$error.value[,1];
data$fishery_years = fishery_compdat$year
data$use_logistic_normal = 1;
#data$LN_resids_centered = 0;
data$ARMA = 1;
ARMA = data$ARMA;
unit_YCS = rep(1,length(data$years))
YCS_start = rep(0,length(data$years) - 1)

pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma_fishery = 0.1, log_norm_phi_fishery = c(-0.5,-0.5), log_norm_sigma_survey = 0.1, log_norm_phi_survey = c(-0.5,-0.5))

## check it works
setwd(DIR$TMB)
compile("model.cpp")
dyn.load(dynlib("model"))
         
obj <- MakeADFun(data=data,parameters=pars, DLL = "model");
#obj_run = obj$report();
rm(obj_run);

obj$fn()
obj$gr() ## if a parameter is zero something is not right, means it is not-identifiable

lower = c(16, 0.01,1,1,1,1,-2.995732,rep(-10, length(YCS_start)), -3,-9,-9, -3,-9,-9)
upper = c(19, 1,20 ,20 ,20 ,20, 0.9162907,rep(10, length(YCS_start)),0.7,-0.01,-0.01,0.7,-0.01,-0.01)

obj$fn()
obj$gr()

obj$env$tracepar = TRUE

names(lower) = names(obj$par)
names(upper) = names(obj$par)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, lower=lower, upper=upper, control = list(abs.tol = 0.0001, rel.tol = 0.0001, eval.max = 20000, iter.max = 20000))
opt_obj_run = obj$report();

## look at the likelihood components
opt_obj_run$neg_ll
opt_obj_run$neg_ll_survey_bio
opt_obj_run$neg_ll_survey_age
opt_obj_run$neg_ll_fishery_age

obj$env$last.par
opt_obj_run$true_ycs

## look at pearson residuals 
opt_obj_run$survey_age_pred
opt_obj_run$fishery_age_pred


## Double check the Survey age comp
surv_comp_Dat = list();
surv_comp_Dat$obs = data$survey_at_age_obs
colnames(surv_comp_Dat$obs) = paste0("X",1:19)
surv_comp_Dat$exp = opt_obj_run$survey_age_expectations
colnames(surv_comp_Dat$exp) = paste0("X",1:19)

surv_comp_Dat$N = data$survey_at_age_error
sigma = exp(obj$env$last.par["log_norm_sigma_survey"])
#phi =c(exp(obj$env$last.par[c(33,34)]))
phi =c(exp(obj$env$last.par[c(36,37)]))

sepbysex=F;
sexlag=F;
robust=F;
#ARMA=F

value =  NLLlogistnorm(surv_comp_Dat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
resids_LN =  Sres.logistnorm(surv_comp_Dat,sigma = sigma,phi,centred=T,sepbysex=F,sexlag=F,
                             ARMA=T)

resids_LN1 = opt_obj_run$survey_age_pred
rbind(resids_LN[1,],resids_LN1[1,])

value;opt_obj_run$neg_ll_survey_age

## Check the fishery
## Double check the Survey age comp
fish_comp_Dat = list();
fish_comp_Dat$obs = data$fishery_at_age_obs
colnames(fish_comp_Dat$obs) = paste0("X",1:19)
fish_comp_Dat$exp = opt_obj_run$fishery_age_expectations
colnames(fish_comp_Dat$exp) = paste0("X",1:19)

fish_comp_Dat$N = data$fishery_at_age_error
sigma = exp(obj$env$last.par["log_norm_sigma_fishery"])
phi =c(exp(obj$env$last.par[c(33,34)]))


value =  NLLlogistnorm(fish_comp_Dat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
resids_LN =  Sres.logistnorm(fish_comp_Dat,sigma = sigma,phi,centred=T,sepbysex=F,sexlag=F,
                             ARMA=T)
resids_LN1 = opt_obj_run$fishery_age_pred
rbind(resids_LN[1,],resids_LN1[1,])

value;opt_obj_run$neg_ll_fishery_age


## Lets look at the correlation structure of the residuals
Plcor(opt_obj_run$fishery_age_pred)
Plcor(opt_obj_run$survey_age_pred)



## simulating a MVN
library(MASS)
covar = covmat.logistnorm(sigma,phi,binnam = colnames(fish_comp_Dat$obs),sepbysex=F,sexlag=F,ARMA=F)
N = 10000
error = mvrnorm(N,mu = rep(0,ncol(covar)), Sigma = covar)
chris_dist = rmultnorm(N,log(fish_comp_Dat$exp[1,]),covar)
## offset 
vals = sweep(error,MARGIN=2,log(fish_comp_Dat$exp[1,]),'+')
vals[1,]
error[1,] + log(fish_comp_Dat$exp[1,])

rbind(apply(vals,2,mean),
log(fish_comp_Dat$exp[1,]),apply(chris_dist,1,mean))

diag(cov(vals))
diag(covar)

