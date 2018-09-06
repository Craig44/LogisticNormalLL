## (04) Compare_comp_likelihood.R
## @Author C Marsh
## @Date 3/9/2018
## @Description
## This script uses actual observations from the HAK assessment and compares fits in a 'fixed' effects model to look at residuals to convince me that
## all this trouble was worth it.

## Add Libraries
library(TMB)
library(TMBdebug)
library(casal)
library(cyrils)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(lemon) ## grid_arrange_shared_legend()
library(grid)
## Add other dependency funcitons
setwd("C:/Work/Projects/PhD/Logistic_normal/Rstudio");
source("Initialisation.r");
source("Chris_original_code.r")
source("C:/Work/Software/TMB/MyUtilityFuns/fix_pars.R")
Hak  = extract.fits(path = DIR$'csl_files', file = "HAK.log");

survey_compdat = Hak$subaTANageDEC
survey_bio = Hak$subaTANbiomassDEC
length(survey_bio$year); length(survey_compdat$year)
surv_ndx = survey_compdat$year %in% survey_bio$year
fishery_compdat = Hak$subaOBSage

load(file = make.filename(file = "tmb.RData", path = DIR$TMB))
## add the new switches
data$survey_at_age_obs = as.matrix(survey_compdat$obs[surv_ndx,]);
data$survey_at_age_error = survey_compdat$error.value[surv_ndx,1];
data$survey_years = survey_compdat$year[surv_ndx]

data$survey_biomass_obs = survey_bio$obs
data$survey_biomass_error = survey_bio$error.value

data$fishery_at_age_obs = as.matrix(fishery_compdat$obs);
data$fishery_at_age_error = fishery_compdat$error.value[,1];
data$fishery_years = fishery_compdat$year

data$ages = 1:19
data$years = 1989:2013
data$ageing_error = matrix(0,length(data$ages),length(data$ages))
diag(data$ageing_error) = 1;

data$use_logistic_normal = 1;

#data$LN_resids_centered = 0;
data$LN_AR_structure = 4;
ARMA = T;
unit_YCS = rep(1,length(data$years))
YCS_start = rep(0,length(data$years) - 1)

pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma_fishery = 0.1, log_norm_phi_fishery = c(0,0), log_norm_sigma_survey = 0.1, log_norm_phi_survey = c(0,0))
#parm_labels_to_fix = c("log_norm_sigma_fishery","log_norm_phi_survey","log_norm_sigma_fishery","log_norm_phi_fishery")
#mapped_pars = fix_pars(pars,parm_labels_to_fix);

data$steepness = data$h
data$a_weight = data$a
data$b_weight = data$b
data$L_inf_vb = data$L_inf
data$k_vb = data$k
data$t0_vb = data$t0
data$M_mort = data$M

## load model
setwd(DIR$TMB)
compile("model.cpp")
dyn.load(dynlib("model"))

# if not estimatng LN parameters
#obj <- MakeADFun(data=data,parameters=pars, DLL = "model", map= mapped_pars);

### -------------- LN3m ---------------
pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma_fishery = 0.1, norm_phi_fishery = c(0.1,0.1), log_norm_sigma_survey = 0.1, norm_phi_survey = c(0.1,0.1))
data$LN_AR_structure = 4;
data$use_logistic_normal = 1;
obj_LN3m <- MakeADFun(data=data,parameters=pars, DLL = "model");

lower = c(16, 0.01,1,1,1,1,-2.995732,rep(-10, length(YCS_start)), -3,-0.99,-Inf, -3,-0.99,-Inf)
upper = c(19, 1,20 ,20 ,20 ,20, 0.9162907,rep(10, length(YCS_start)),0.7,0.99,Inf,0.7,0.99,Inf)

names(upper) = names(pars)
names(lower) = names(lower)

obj_LN3m$fn()
obj_LN3m$gr()

opt_LN3m <- nlminb(start = obj_LN3m$par, objective = obj_LN3m$fn, gradient = obj_LN3m$gr, lower=lower, upper=upper, control = list(abs.tol = 0.0001, rel.tol = 0.0001, eval.max = 20000, iter.max = 20000))
opt_LN3m$convergence ## needs to be 0
obj_LN3m_report = obj_LN3m$report();
obj_LN3m_simulate = obj_LN3m$simulate(complete = T);

### -------------- LN3 ---------------
data$LN_AR_structure = 3
pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma_fishery = 0.1, norm_phi_fishery = c(0,0.5), log_norm_sigma_survey = 0.1, norm_phi_survey = c(0,0.5))
lower = c(16, 0.01,1,1,1,1,-2.995732,rep(-10, length(YCS_start)), -3,-1,0.001, -3,-1,0.001)
upper = c(19, 1,20 ,20 ,20 ,20, 0.9162907,rep(10, length(YCS_start)),0.7,1,0.99,0.7,1,0.99)

names(upper) = names(pars)
names(lower) = names(lower)

obj_LN3 <- MakeADFun(data=data,parameters=pars, DLL = "model");
obj_LN3$fn()
obj_LN3$gr()


opt_LN3 <- nlminb(start = obj_LN3$par, objective = obj_LN3$fn, gradient = obj_LN3$gr, lower=lower, upper=upper, control = list(abs.tol = 0.0001, rel.tol = 0.0001, eval.max = 20000, iter.max = 20000))
opt_LN3$convergence
obj_LN3_report = obj_LN3$report();

### -------------- LN2 ---------------
pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma_fishery = 0.1, norm_phi_fishery = c(0.5), log_norm_sigma_survey = 0.1, norm_phi_survey = c(0.5))
data$LN_AR_structure = 2

lower = c(16, 0.01,1,1,1,1,-2.995732,rep(-10, length(YCS_start)), -3,-0.99, -3,-0.99)
upper = c(19, 1,20 ,20 ,20 ,20, 0.9162907,rep(10, length(YCS_start)),0.7,0.99,0.7,0.99)
names(upper) = names(pars)
names(lower) = names(pars)


obj_LN2 <- MakeADFun(data=data,parameters=pars, DLL = "model");
obj_LN2$gr()
obj_LN2$fn()

pars2 = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma_fishery = 0.1, norm_phi_fishery = c(0.5), log_norm_sigma_survey = 0.1, norm_phi_survey = c(0.5))
obj2_LN2 <- MakeADFun(data=data,parameters=pars2, DLL = "model");

pars3 = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma_fishery = 0.1, norm_phi_fishery = c(-0.5), log_norm_sigma_survey = 0.1, norm_phi_survey = c(-0.5))
obj3_LN2 <- MakeADFun(data=data,parameters=pars3, DLL = "model");
c(obj_LN2$fn(),obj2_LN2$fn(), obj3_LN2$fn())

opt_LN2 <- nlminb(start = obj_LN2$par, objective = obj_LN2$fn, gradient = obj_LN2$gr, lower=lower, upper=upper, control = list(abs.tol = 0.0001, rel.tol = 0.0001, eval.max = 20000, iter.max = 20000))
opt2_LN2 <- nlminb(start = obj2_LN2$par, objective = obj2_LN2$fn, gradient = obj2_LN2$gr, lower=lower, upper=upper, control = list(abs.tol = 0.0001, rel.tol = 0.0001, eval.max = 20000, iter.max = 20000))
opt3_LN2 <- nlminb(start = obj3_LN2$par, objective = obj3_LN2$fn, gradient = obj3_LN2$gr, lower=lower, upper=upper, control = list(abs.tol = 0.0001, rel.tol = 0.0001, eval.max = 20000, iter.max = 20000))

opt_LN2$convergence
rbind(opt$par,opt2_LN2$par, opt3_LN2$par)

obj_LN2_report = obj_LN2$report();
obj_LN2_simulate = obj_LN2$simulate(complete = T);

### -------------- LN1 ---------------
pars = list(log_R0 = 18, q = 0.2, s_a50 = 3, s_ato95 = 3, f_a50 = 3, f_ato95 = 3, log_sigma_r = log(0.5), YCS = YCS_start, log_norm_sigma_fishery = 0.1, norm_phi_fishery = c(0), log_norm_sigma_survey = 0.1, norm_phi_survey = c(0))
data$LN_AR_structure = 1

lower = c(16, 0.01,1,1,1,1,-2.995732,rep(-10, length(YCS_start)), -3,-1, -3,-1)
upper = c(19, 1,20 ,20 ,20 ,20, 0.9162907,rep(10, length(YCS_start)),0.7,1,0.7,1)
names(upper) = names(pars)
names(lower) = names(lower)

parm_labels_to_fix = c("norm_phi_survey", "norm_phi_fishery")
mapped_pars = fix_pars(pars,parm_labels_to_fix);

obj_LN1 <- MakeADFun(data=data,parameters=pars, DLL = "model", map = mapped_pars);
obj_LN1$fn()
obj_LN1$gr()
opt_LN1 <- nlminb(start = obj_LN1$par, objective = obj_LN1$fn, gradient = obj_LN1$gr, lower=lower, upper=upper, control = list(abs.tol = 0.0001, rel.tol = 0.0001, eval.max = 20000, iter.max = 20000))
opt_LN1$convergence
obj_LN1_report = obj_LN1$report();
obj_LN1_simulate = obj_LN1$simulate(complete = T);

### --------------- Multinomial ------------
data$use_logistic_normal = 0
parm_labels_to_fix = c("norm_phi_survey","log_norm_sigma_survey","log_norm_sigma_fishery", "norm_phi_fishery")

obj_multi <- MakeADFun(data=data,parameters=pars, DLL = "model", map = mapped_pars);
obj_multi$fn();
obj_multi$gr();
opt_multi <- nlminb(start = obj_multi$par, objective = obj_multi$fn, gradient = obj_multi$gr, lower=lower, upper=upper, control = list(abs.tol = 0.0001, rel.tol = 0.0001, eval.max = 20000, iter.max = 20000))
opt_multi$convergence

obj_multi_report = obj_multi$report();
obj_multi_simulate = obj_multi$simulate(complete = T);

###############################
## To double triple check 
##############################
compdat = list()
compdat$obs = data$survey_at_age_obs
compdat$N = data$survey_at_age_error
colnames(compdat$obs) = paste0("X",1:19)
compdat$exp = obj_LN1_report$survey_age_expectations
colnames(compdat$exp) = colnames(compdat$obs)
sigma = exp(opt_LN1$par["log_norm_sigma_survey"])
NLLlogistnorm(compdat,sigma = sigma,phi = 0,covmat=NULL,sepbysex=F, sexlag=F, robust=F, ARMA=F)
obj_LN1_report$neg_ll_survey_age

compdat$exp = obj_LN2_report$survey_age_expectations
colnames(compdat$exp) = colnames(compdat$obs)
sigma = exp(opt_LN2$par["log_norm_sigma_survey"])
phi = (opt_LN2$par["norm_phi_survey"])
NLLlogistnorm(compdat,sigma = sigma,phi = phi,covmat=NULL,sepbysex=F, sexlag=F, robust=F, ARMA=F)
obj_LN2_report$neg_ll_survey_age

compdat$exp = obj_LN3_report$survey_age_expectations
colnames(compdat$exp) = colnames(compdat$obs)
sigma = exp(opt_LN3$par["log_norm_sigma_survey"])
phi = (opt_LN3$par[names(opt_LN3$par) %in% "norm_phi_survey"])
ph2 = -1 + (2 - abs(phi[1])) * phi[2];
phi = c(phi[1], ph2)
NLLlogistnorm(compdat,sigma = sigma,phi = phi,covmat=NULL,sepbysex=F, sexlag=F, robust=F, ARMA=F)
obj_LN3_report$neg_ll_survey_age

compdat$exp = obj_LN3m_report$survey_age_expectations
colnames(compdat$exp) = colnames(compdat$obs)
sigma = exp(opt_LN3m$par["log_norm_sigma_survey"])
phi = (opt_LN3m$par[names(opt_LN3m$par) %in% "norm_phi_survey"])
NLLlogistnorm(compdat,sigma = sigma,phi = phi,covmat=NULL,sepbysex=F, sexlag=F, robust=F, ARMA=T)
obj_LN3m_report$neg_ll_survey_age


#######################################################################
## Compare fits to Biomass and frequency data
#######################################################################

survey_bio_ci  = lognormal_CI(data$survey_biomass_obs, data$survey_biomass_error, CI = 0.95)

png(filename = make.filename(file = "biomass_across_Obs_models.png",path = DIR$figures),units ="cm", width = 27, height = 17, res = 1200)
plot(data$survey_years, data$survey_biomass_obs, pch = 20, col = "blue", xlab = "years", ylab = "Biomass (t)", ylim = c(0,7000), main = "Survey Biomass")
segments(data$survey_years,survey_bio_ci$upper, data$survey_years,survey_bio_ci$lower)
## fits
lines(data$survey_years, obj_LN3m_report$survey_biomass_expectations, lty = 2, col = "purple",lwd = 2)
lines(data$survey_years, obj_LN3_report$survey_biomass_expectations, lty = 2, col = "orange",lwd = 2)
lines(data$survey_years, obj_LN2_report$survey_biomass_expectations, lty = 2, col = "red",lwd = 2)
lines(data$survey_years, obj_LN1_report$survey_biomass_expectations, lty = 2, col = "darkgreen",lwd = 2)
lines(data$survey_years, obj_multi_report$survey_biomass_expectations, lty = 2, col = "cyan",lwd = 2)
legend('topright', legend = c("LN3m", "LN3", "LN2","LN1","Multinom"), col = c("purple","orange","red","darkgreen","cyan"), lty = 2, lwd = 2, cex = 0.8)
dev.off()
#### Age frequency look at standardised residuals
## Fishery

p1 = plot_resids(obj_LN3m_report$survey_age_pred, bin_labels = data$ages,years = data$survey_years)
p2 = plot_resids(obj_LN3m_report$fishery_age_pred, bin_labels = data$ages,years = data$fishery_years)
## add margins
p1 = p1 + theme(plot.margin = unit(c(2,0,0,0), "cm"))
p2 = p2 + theme(plot.margin = unit(c(2,0,0,0), "cm"))
png(filename = make.filename(file = "LN3m_resids.png",path = DIR$figures),units ="cm", width = 27, height = 17, res = 120)
grid_arrange_shared_legend(p1, p2, nrow = 1)
grid.text(label = "Survey", x = 0.25,y = 0.95, gp=gpar(fontsize=20))
grid.text(label = "Fishery", x = 0.75,y = 0.95, gp=gpar(fontsize=20))
dev.off();

## look at the correlation of residuals
par(mfrow = c(3,2), mar = c(2,2,1,1))
Plcor(cor(data$survey_at_age_obs -  obj_multi_report$survey_age_expectations))
Plcor(cor(data$survey_at_age_obs -  obj_LN1_report$survey_age_expectations))
Plcor(cor(data$survey_at_age_obs -  obj_LN2_report$survey_age_expectations))
Plcor(cor(data$survey_at_age_obs -  obj_LN3_report$survey_age_expectations))
Plcor(cor(data$survey_at_age_obs -  obj_LN3m_report$survey_age_expectations))

par(mfrow = c(3,2), mar = c(2,2,1,1))
Plcor(cor(data$fishery_at_age_obs -  obj_multi_report$fishery_age_expectations))
Plcor(cor(data$fishery_at_age_obs -  obj_LN1_report$fishery_age_expectations))
Plcor(cor(data$fishery_at_age_obs -  obj_LN2_report$fishery_age_expectations))
Plcor(cor(data$fishery_at_age_obs -  obj_LN3_report$fishery_age_expectations))
Plcor(cor(data$fishery_at_age_obs -  obj_LN3m_report$fishery_age_expectations))

### plot SSBs
plot(data$years, obj_LN3m_report$SSB, lty = 2, col = "purple",lwd = 2, xlab = "years", ylab = "SSB (t)", type = "l", ylim = c(0,200000))
lines(data$years, obj_LN3_report$SSB, lty = 2, col = "orange",lwd = 2)
lines(data$years, obj_LN2_report$SSB, lty = 2, col = "red",lwd = 2)
lines(data$years, obj_LN1_report$SSB, lty = 2, col = "darkgreen",lwd = 2)
lines(data$years, obj_multi_report$SSB, lty = 2, col = "cyan",lwd = 2)
legend('topright', legend = c("LN3m", "LN3", "LN2","LN1","Multinom"), col = c("purple","orange","red","darkgreen","cyan"), lty = 2, lwd = 2)


## Perhaps we should weight the Multinomial data.

## look at negative log-likelihood
obj_LN3m_report$neg_ll
obj_LN3_report$neg_ll
obj_LN2_report$neg_ll
obj_LN1_report$neg_ll

