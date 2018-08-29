###############################################################
##### Now test for the simulating component of the likelihood
###############################################################
setwd("C:/Craig/projects/2016/DEE2015-02(TrawlSurveySimulation)/Logistic_normal/")
## read in Chris's functions we need Estlogistnorm() and rlogistnorm()
source(file = "Chris_original_code.R") 

## read in HAK 1 data and format it for Chris's functions
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
## Estimate the sigma and phi parameters
Observer = Estlogistnorm(compdat_obs,ARMA = F) ## try different ARMA values

## Now simulate data off the expectations and see if this simulator recreates the 
## correlation pattern
N = 5000 # number of simulations
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
par(mfrow=c(1,1))
## plot the patter
plot(1, type="n", xlab="Lag", ylab="Residual correlatons",main = "Observed", xlim = c(0,nrow(Resid_cor)), ylim = c(-1,1))
for (i in 1:(ncol(Resid_cor) - 1)) {
 for(j in (i + 1):ncol(Resid_cor)) {
  points(x = (j - i), y = Resid_cor[i,j], col = "black", pch = 20)
 }
}
## Now summarise and add simualated data to see if it also has this correlation structure.
MAT = apply(array(apply(resid_dat,3,cor),c(obs_bins ,obs_bins,N)),c(1,2),mean) 
#plot(1, type="n", xlab="Lag", ylab="Raw residual correlatons",main = "Simulated", xlim = c(0,nrow(Resid_cor)), ylim = c(-1,1))
for (i in 1:(ncol(MAT) - 1)) {
 for(j in (i + 1):ncol(MAT)) {
  points(x = (j - i), y = MAT[i,j], col = "red", pch = 20)
 }
}
legend('topright',col = c("black","red"),c("Observed","Simulated"),pch = c(20,20))


### Now look at my R and C++ code to see if I get the same mean and expectation
