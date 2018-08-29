## Casal2_spec.r
## @Author C Marsh
## @Date 8/10/2016
## @Description
## This script is setting up the logistic normal likelihood for composition data

## Add Libraries
library(casal)

## Add other dependency funcitons
setwd("C:/Work/Projects/PhD/Logistic_normal/Rstudio");
source("Chris_original_code.r")
source("Initialisation.r");

## Import example Casal file as this is 
Hak  = extract.fits(path = DIR$'csl_files', file = "HAK.log");
Hok  = extract.fits(path = DIR$'csl_files', file = "HOKI.log");
Lin  = extract.fits(path = DIR$'csl_files', file = "LIN.log");


## Run through Chris code manually with all stocks

#########
## HAK
#########
compdat = Hak$subaTANageDEC
## restructure for Chris's Code
compdat$obs =  as.matrix(Hak$subaTANageDEC$obs)
compdat$exp =  round(as.matrix(Hak$subaTANageDEC$fits),5)
compdat$N = Hak$subaTANageDEC$error.value[,1]

# save this information to text files to compare with C++ code
write.table(x = compdat$obs, file = make.filename(file = "observed_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = compdat$exp, file = make.filename(file = "expected_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = compdat$N, file = make.filename(file = "error_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)

## run the program
setwd(DIR$`C++_work`);
system("LogNorm.exe notsexed 1",ignore.stdout = TRUE,ignore.stderr = TRUE)
setwd(DIR$R)
C_results = as.vector(read.table(make.filename(file = "results.txt", path = DIR$`C++_work`), header = F))
## Check against Chris's R-Code


sepbysex=F;
sexlag=F;
robust=F;
ARMA=F
sigma = 0.9;
phi = c(0.2)

true_results = vector();
true_results[1] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
sigma = 0.2;
true_results[2] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
sigma = 0.4;
phi = c(0.2,0.5)
true_results[3] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
sigma = 0.6;
phi = c(0.4,0.1)
true_results[4] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
ARMA = TRUE
true_results[5] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)


true_results - C_results
cbind(true_results, t(C_results))

#########
## LIN
#########
compdat = Lin$Tangaroa_propn_at_age_summer
## restructure for Chris's Code
compdat$obs =  as.matrix(Lin$Tangaroa_propn_at_age_summer$obs)
compdat$exp =  round(as.matrix(Lin$Tangaroa_propn_at_age_summer$fits),5)
compdat$N = Lin$Tangaroa_propn_at_age_summer$error.value[,1]

# save this information to text files to compare with C++ code
write.table(x = compdat$obs, file = make.filename(file = "observed_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = compdat$exp, file = make.filename(file = "expected_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(x = compdat$N, file = make.filename(file = "error_data.txt", path = DIR$`C++_work`),quote = FALSE, row.names = FALSE, col.names = FALSE)

## run the program
setwd(DIR$`C++_work`);
system("LogNorm.exe sex 3",ignore.stdout = TRUE,ignore.stderr = TRUE)
setwd(DIR$R)
C_results = as.vector(read.table(make.filename(file = "results.txt", path = DIR$`C++_work`), header = F))
## Check against Chris's R-Code

sepbysex=F;
sexlag=F;
robust=F;
ARMA=F
sigma = 0.9;
phi = c(0.2)

true_results = vector();
true_results[1] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
sigma = 0.2;
true_results[2] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
sigma = 0.4;
phi = c(0.2,0.5)
true_results[3] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
sigma = 0.6;
phi = c(0.4,0.1)
true_results[4] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
ARMA = TRUE
true_results[5] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
sepbysex=T;
ARMA = TRUE
sigma = 0.9;
phi = c(0.2,0.3)
true_results[6] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)

sepbysex=F;
sexlag = T;
true_results[7] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)
sexlag = F;
true_results[8] = NLLlogistnorm(compdat,sigma = sigma,phi,covmat=NULL,sepbysex=sepbysex, sexlag=sexlag, robust=robust, ARMA=ARMA)

true_results - C_results
cbind(true_results, t(C_results))


