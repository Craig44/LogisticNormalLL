#
# test Stan model
#
#
library("rstan")
library(casal)
source("../Chris_original_code.r")

## change some settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("C:/Work/Projects/PhD/Logistic_normal/Stan")
source("../Rstudio/Initialisation.r");

stan_model = stan_model("SimpleStan.stan")

## get some data
Hak  = extract.fits(path = DIR$'csl_files', file = "HAK.log");

survey_compdat = Hak$subaTANageDEC

stan_data = list()
stan_data$Y = nrow(survey_compdat$obs)
stan_data$A = ncol(survey_compdat$obs)
stan_data$Observations = as.matrix(survey_compdat$obs)
stan_data$expectations = as.matrix(survey_compdat$fits)
covmat = covmat.logistnorm(sigma = 1.11, phi = c(0.66), binnam = paste0("X",1:stan_data$A))
stan_data$covariance = covmat
stan_data$sample_size_by_year = as.numeric(survey_compdat$error.value[,1])
stan_data$centered = 0
stan_data$covariance_structure = 1
stan_data$sigma = 1.11
stan_data$rho = 0.66
covmat = covmat.logistnorm(sigma = stan_data$sigma, phi = stan_data$rho, binnam = paste0("X",1:stan_data$A))

init_pars = list()
init_pars$p = 0.5

## run model
stan_out= optimizing(stan_model,data=stan_data,init=init_pars,
                      iter=1,algorithm="LBFGS",
                      verbose=F,as_vector = FALSE)

covmat - stan_out$par$covariance_matrix
## run model
stan_data$covariance_structure = 2
stan_out1= optimizing(stan_model,data=stan_data,init=init_pars,
                     iter=1,algorithm="LBFGS",
                     verbose=F,as_vector = FALSE)

negloglik

NLLlogistnorm(compdat = survey_compdat, covmat = covmat) 
stan_out$par$neg_ll
stan_out$par$gmean_exp
stan_out$par$gmean_obs


stan_out1$par$covariance_matrix

isSymmetric(stan_out$par$covariance_matrix)
diag(stan_out$par$covariance_matrix)
stan_data$sigma * stan_data$sigma 

stan_out$par$cholesky_cor


stan_out$par$cholesky_cor


rowSums(stan_out$par$cholesky_cor)
colSums(stan_out$par$cholesky_cor)

stan_out$par$covariance_matrix %*% t(stan_out$par$covariance_matrix)

stan_out$par$cholesky_cor %*% t(stan_out$par$cholesky_cor)[2,2]
t(stan_out$par$cholesky_cor) %*% (stan_out$par$cholesky_cor)

diag(stan_out$par$cholesky_cor %*% t(stan_out$par$cholesky_cor))

## what it should be
rhovec <- getrho(stan_data$rho,stan_data$A,F)
corr_mat = matrix(0,stan_data$A,stan_data$A)
for(i in 1:(stan_data$A-1))
  corr_mat[abs(row(corr_mat)-col(corr_mat))==i] <- rhovec[i]
diag(corr_mat) = 1
chol_ar = chol(corr_mat)
sum(chol_ar)
rowSums(chol_ar)
colSums(chol_ar)
t(chol_ar) %*% (chol_ar)

## compare R with Stan
t(chol_ar)[1:5,1:5]
stan_out$par$cholesky_cor[1:5,1:5]



diag(stan_out$par$covariance_matrix)

stan_out$par$log_total_obs
stan_out$par$log_det
stan_out$par$ww - ww

stan_res = stan_out$par$sres
survey_compdat$obs = as.matrix(survey_compdat$obs)
survey_compdat$exp = as.matrix(survey_compdat$fits)
survey_compdat$N = as.numeric(survey_compdat$error.value[,1])

chris_red = Sres.logistnorm(compdat = survey_compdat,sigma = 1.11, phi = c(0.7, -0.4))
chris_res = Sres.logistnorm(compdat = survey_compdat,centred = FALSE, sigma = 1.11, phi = c(0.7, -0.4))
stan_res_non = stan_out$par$sres_non_centered

chris_res[1,]
stan_res_non[1,]

chris_red[1,]
stan_res[1,]


stan_out$par$gmean_exp
stan_out$par$gmean_obs


ar = stan_data$rho
ma = stan_data$ma

A = stan_data$A
mat = matrix(0,A,A)
gamma = vector() 
diag(mat) = 1;
for (i in 2:A) { 
  gamma[i - 1] = pow(ar, i - 1); 
  for (j in 1:(i - 1)) { 
    mat[i, j] = gamma[i - j]; 
    mat[j, i] = gamma[i - j]; 
  } 
} 
rhovec
mat[1,]



mat = matrix(0,A,A)
gamma = vector() 
diag(mat) = 1 + ma^2 + 2 * ar * ma;
gamma[1] = (1 + ar * ma) * (ar + ma); 
for (i in 2:A) { 
  gamma[i] = gamma[1] * pow(ar, i - 1); 
  for (j in 1:(i - 1)) { 
    mat[i, j] = gamma[i - j]; 
    mat[j, i] = gamma[i - j]; 
  } 
} 
final = (1 / (1 - ar^2) * mat)

gamma * (1 / (1 - ar^2))
getrho(c(ar,ma),A,T)


