# test_C_code.r
#@author C.Marsh
#@date 15/10/16
#@description
## This script will be the basis for a unit test in Casal2, it will check all the negative likelihood scores for a range of inputs
## it will also test the Random number generator.

## Add Libraries
library(casal)

## Add other dependency funcitons
## Will be using the HAK1 example for unisex test case.
source("C:/Craig/projects/2016/DEE2015-02(TrawlSurveySimulation)/R/Initialisation.r");
source(make.filename(file = "auxiliary.functions.R", path = DIR$'General functions'));
source(make.filename(file = "Logistic_normal\\Chris's_r_code.r", path = DIR$'Base'))

## Import example Casal file as this is 
casal_fits  = extract.fits(path = DIR$'HAK1_csl', file = "out.log");

###################################################
####  Single sex comp data
###################################################
compdat = casal_fits$subaTANageDEC
## restructure for Chris's Code
compdat$obs =  as.matrix(casal_fits$subaTANageDEC$obs)
compdat$exp =  round(as.matrix(casal_fits$subaTANageDEC$fits),5)
compdat$N = casal_fits$subaTANageDEC$error.value[,1]

sigma = 0.243
phi=0.0;
covmat=NULL;
sepbysex=F;
sexlag=F;
robust=F;
ARMA=F

#########
## Test 1
#########
## unisex data, sigma and no phi
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex=F, sexlag=F, robust=F, ARMA=F)
## C++ answer
## 936.652
#########
## Test 2
#########
## unisex data, sigma and single phi
phi= 0.435;
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex=F, sexlag=F, robust=F, ARMA=F)
## C++ answer
## 1509.43
#########
## Test 3
#########
## unisex data, sigma and two phi and ARMA = F
phi= c(0.235,-0.284);
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex=F, sexlag=F, robust=F, ARMA=F)
## C++ answer
## 1280.26
#########
## Test 4
#########
## unisex data, sigma and two phi and ARMA = T
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex=F, sexlag=F, robust=F, ARMA=T)
## C++ answer
## 982.383

#########
## Test 5
#########
## unisex data, sigma and two phi and ARMA = T,robust=T
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex=F, sexlag=F, robust=T, ARMA=T)
## C++ answer
## -2925.01


###################################################
####  Sexed comp data
###################################################
## Use the ling example now.
casal_fits  = extract.fits(path = DIR$'LIN3_csl', file = "out.log");
compdat = casal_fits$Tangaroa_propn_at_age_Jan
## restructure for Chris's Code
compdat$obs =  as.matrix(casal_fits$Tangaroa_propn_at_age_Jan$obs)
compdat$exp =  round(as.matrix(casal_fits$Tangaroa_propn_at_age_Jan$fits),5)
compdat$N = casal_fits$Tangaroa_propn_at_age_Jan$error.value[,1]

#########
## Test 6
#########
## Sexed data, sepsex = T
sigma = 0.243
phi= c(0.235,-0.284);
sepbysex=T
robust = F
ARMA=T
sepbysex=T
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex, sexlag=F, robust, ARMA)

## C++ answer
## -2416.6
#########
## Test 7
#########
## Sexed data, sepsex = T
sigma = 0.243
phi= c(0.235,-0.284);
sepbysex=T
robust = F
ARMA=T
sepbysex=F
sexlag=T
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex, sexlag, robust, ARMA)
## C++ answer
## -2428.48

## now just check that the robustification works with the sexes as well
#########
## Test 8
#########
## Sexed data, sepsex = T
sepbysex=T
robust = T
ARMA=T
sepbysex=T
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex, sexlag=F, robust, ARMA)
## C++ answer
## -4926.83

#########
## Test 9
#########
sepbysex=T
robust = T
ARMA=T
sepbysex=F
sexlag=T
NLLlogistnorm(compdat,sigma ,phi,covmat=NULL,sepbysex, sexlag, robust, ARMA)
## C++ answer
## -4933.73

## Now think about implementing this in Casal2, the algorithm.
## 
## DoReset()
## Calculate the covar matrix
## do the inversions

## add a null 2d matrix on the comparisons for a covariace matrix.
## in the likelihood
## Calculate score()



###############################################################
##### Now test for the simulating component of the likelihood
###############################################################
## read in HAK data and format it for Chris's functions
setwd("C:/Craig/projects/2016/DEE2015-02(TrawlSurveySimulation)/Logistic_normal/")
source(file = "Chris_original_code.R")

HAK1 = read.csv("Test_code/Casal2_example1.csv")

## start with Observer data
N = length(unique(HAK1$year))
obs_bins = length(unique(HAK1$age))
error = vector();
Obs = Exp = matrix(NA, nrow = N,ncol =  obs_bins)
## convert to a structure that Chris code will accept
year_iter = 1;
for (year in unique(HAK1$year)) {
  row_index = HAK1$year ==  year;
  Obs[year_iter,] = HAK1[row_index,"observed"];
  Exp[year_iter,] = HAK1[row_index,"expected"];
  error[year_iter] = unique(HAK1[row_index,"error_value"])
  year_iter = year_iter + 1;
}
colnames(Obs) = paste("X",unique(HAK1$age),sep = "")
compdat_obs = list();
## restructure for Chris's Code
compdat_obs$obs =  Obs
compdat_obs$exp =  Exp
compdat_obs$N = error
## re-estimate
Observer = Estlogistnorm(compdat_obs)

## Now simulate data off the expectations many times and calculate teh mean and variance and compare with 
## my own R code that will be used to convert to C++ code
N = 50 # number of simulations
Sim_dat = resid_dat = array(0,dim=c(nrow(compdat_obs$exp),ncol(compdat_obs$exp),N))
for(n in 1:N) {
  for(i in 1:nrow(compdat_obs$exp))
    Sim_dat[i,,n] = rlogistnorm(n = 1,compdat_obs$exp[i,],sigma = Observer["sigma"], phi = c(Observer["phi1"],Observer["phi2"]))
  ## for each simualation calculate the raw residual between simualated and expected and save it
  resid_dat[,,n] = Sim_dat[,,n] - compdat_obs$exp;
}

## Now look at the observered residual correlations and compare with the simualted residual correlations
Resid = compdat_obs$obs - compdat_obs$exp
Resid_cor = cor(Resid)

## plot the patter
plot(1, type="n", xlab="Lag", ylab="Residual correlatons", xlim = c(0,nrow(Resid_cor)), ylim = c(-1,1))
for (i in 1:(ncol(Resid_cor) - 1)) {
 for(j in (i + 1):ncol(Resid_cor)) {
  points(x = (j - i), y = Resid_cor[i,j], col = "black", pch = 20)
 }
}
## Now summarise and add simualated data to see if it also has this correlation structure.
MAT = apply(array(apply(resid_dat,3,cor),c(obs_bins ,obs_bins,N)),c(1,2),mean) 
for (i in 1:(ncol(MAT) - 1)) {
 for(j in (i + 1):ncol(MAT)) {
  points(x = (j - i), y = MAT[i,j], col = "red", pch = 20)
 }
}
legend('topright',col = c("black","red"),c("Observed","Simulated"),pch = c(20,20))

##################################
## Now do the survey observation #
##################################
HAK1_survey = read.csv("Test_code/Casal2_subantartic.csv")
N = length(unique(HAK1_survey$year))
bins = length(unique(HAK1_survey$age))
error = vector();
Obs = Exp = matrix(NA, nrow = N,ncol =  bins)
## convert to a structure that Chris code will accept
year_iter = 1;
for (year in unique(HAK1_survey$year)) {
  row_index = HAK1_survey$year ==  year;
  Obs[year_iter,] = HAK1_survey[row_index,"observed"];
  Exp[year_iter,] = HAK1_survey[row_index,"expected"];
  error[year_iter] = unique(HAK1_survey[row_index,"error_value"])
  year_iter = year_iter + 1;
}
colnames(Obs) = paste("X",unique(HAK1_survey$age),sep = "")
compdat = list();
## restructure for Chris's Code
compdat$obs =  Obs
compdat$exp =  Exp
compdat$N = error

## re-estimate
Survey = Estlogistnorm(compdat)

######################################################################
### test multivariate normal generator that will be used in the logistic normal 
### simulator in Casal2. The algorithm is base on chloskey decomposition of the 
### covariance matrix
######################################################################
covar = matrix(c(1.4,0.6,0.6,1.54), nrow = 2, byrow = TRUE)
mu = c(3.1,1.35)
cho_covar = chol(covar)
nrow(covar)

rmultinorm = function(n = 100, covar) {
  store_mat = matrix(0.0,nrow = nrow(covar), ncol = n)
  cho_covar = (chol(covar))
  for( k in 1:n) {
  rng_norm = rnorm(n = nrow(covar))
    for (i in 1:nrow(covar)) {
      for(j in 1:nrow(covar)) {
        store_mat[i,k] = store_mat[i,k] +  cho_covar[j,i] * rng_norm[j]
      }
    }
  }
  t(store_mat)
}

rngs = rmultinorm(n = 10000,covar)
apply(rngs,2,mean)
cov(rngs)
plot(rngs,pch = 20)
library(mvtnorm)
actual = rmvnorm(n = 10000, sigma = covar)
cov(actual)
apply(actual,2,mean)

points(actual, col = "red",pch = 20)

## take into account a mean
rmultinorm = function(n = 100, covar,mu = rep(0,nrow(covar))) {
  store_mat = matrix(0.0,nrow = nrow(covar), ncol = n)
  cho_covar = (chol(covar))
  for( k in 1:n) {
  rng_norm = rnorm(n = nrow(covar))
  #result = mu + t(cho_covar) %*% (rng_norm)
  ##store_mat[,k] = result  
    for (i in 1:nrow(covar)) {
      for(j in 1:nrow(covar)) {
        store_mat[i,k] = store_mat[i,k] + cho_covar[j,i] * rng_norm[j]
      }
      store_mat[i,k] = store_mat[i,k] + mu[i]
    }
  }
  t(store_mat)
}

rngs = rmultinorm(n = 10000,covar, mu = c(0.5,0.6))
apply(rngs,MARGIN = 2, FUN = mean)

cov(rngs)
## compare against an R library for unit testing.
library(mvtnorm)
actual = rmvnorm(n = 10000, mean = c(0.5,0.6),sigma = covar)
points(actual, col = "red")
apply(actual,MARGIN = 2, FUN = mean)
cov(rngs)

##########################################################
## now Fully test the simulator conforms to Chris's rlogistnorm simulator
## Now simulate data and numerically calculate the expectations and variances from Chris simulator.
## do the same for the one used in CASAL2 for testing purposes.
###########################################################
N = 10000 # number of simulations
#Sim_dat = resid_dat = array(0,dim=c(nrow(compdat_obs$exp),ncol(compdat_obs$exp),N))
Sim_dat = matrix(0,ncol = ncol(compdat_obs$exp),nrow = N)

for(n in 1:N) { 
    Sim_dat[n,] = rlogistnorm(n = 1,compdat_obs$exp[3,],sigma = Observer["sigma"], phi = c(Observer["phi1"],Observer["phi2"]))
  ## for each simualation calculate the raw residual between simualated and expected and save it
  #resid_dat[,,n] = Sim_dat[,,n] - compdat_obs$exp;
}

## Write out exactly line for line what the CASAL2 simulator will do.
## calculate nbin standard normal random variables

#expprop vector of proportions
RlogisNormal = function(expprop,covar) {
  nbin = length(expprop)
  Covar = chol(covar);
  observed = vector(); ## or simulated
  normals = rnorm(nbin,0,1)
  Total = 0;
  for(i in 1:nbin) {
    ## loop over all bins
    row_sum = 0.0
    for(j in 1:nbin) {
      ## loop over all bins to take into effect all elements of the covariance matrix.
      row_sum = row_sum + Covar[j,i] * normals[j]; #Covar is chol(covar)
    }
  observed[i] = exp(row_sum + log(expprop[i]));
  Total = Total + observed[i];
  }
  ## now do the transformation
  for(i in 1:nbin) {
    observed[i] = observed[i] / Total;
  }
  return(observed);
}


## calculate the covariance
covar = covmat.logistnorm(sigma = Observer["sigma"],phi = c(Observer["phi1"],Observer["phi2"]),ARMA = F,binnam = colnames(compdat_obs$obs))## get an expectation
expprop = compdat_obs$exp[3,]
Sim = matrix(0,ncol = length(expprop), nrow= N)
## Run Casals simulator
for(n in 1:N) { 
  Sim[n,]=RlogisNormal(expprop,covar)
}

## now calculate the mean age for both sets of simulated data.
## Calculate the mean and variance in each bin 
apply(Sim_dat_casal,2,median)
apply(Sim_dat,2,median)
apply(Sim,2,median)
round(apply(Sim_dat_casal,2,var),4)
round(apply(Sim_dat,2,var),4)
round(apply(Sim,2,var),4)

## slightly differene but close enough for me to be happy. Iassume the difference comes down to the different
## multivariate random number generators
## plot the median simualted data over the range
Chris = apply(Sim_dat,2,median)
CASAL2 = apply(Sim,2,median)

plot(density(Chris),xlab = "ages", yaxt = "n", col = "black",type = "l")
lines(density(CASAL2), col = "red")
## I am satisfied that these are drawing the same conclusion.

