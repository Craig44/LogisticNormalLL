## This file provides R functions used in analyses in the
## Fisheries Research paper titled 'Replacing the Multinomial in Stock
## Assessment Models: a First Step' by R.I.C.C. Francis.
##
## To access these R functions type the following line at the R prompt
## source('thisfile.txt') where 'thisfile.txt' is the name of this file.
##
NLLlogistnorm <-
function(compdat,sigma,phi=0,covmat=NULL,sepbysex=F,
                          sexlag=F,robust=F,ARMA=F)
    ## Given observed and expected composition data (single- or
    ## multi-year) that have a logistic-normal distribution, this
    ## function calculates and returns the negative log-likelihood of
    ## the parameters given the data.
    ##
    ## If Ob and Eb are the observed and expected proportions in bin b
    ## it is assumed that Ob = exp(Xb)/sumc(exp(Xc)) where X is a
    ## multivariate normal vector with E(Xb) = log(Eb) and whose
    ## covariance matrix is either defined by the parameters sigma and
    ## phi or provided (as covmat).
    ##
    ## Between-year weighting uses component N of compdat:
    ## sigma[y] = sigma * sqrt(mean(N)/N[y])) or, if covmat is provided
    ## the covariance matrix for year y is covmat * mean(N)/N[y].  The
    ## correlation parameter, phi, is not modified by between-year weighting.
    ##
    ## compdat - the composition data in a list with components
    ##   obs - a Nbin-vector, or a Nyear x Nbin matrix of Nyear observed
    ##     Nbin-vectors.
    ##   exp - the expected compositions in an Nbin-vector or an
    ##     Nyear x Nbin matrix (can use an Nbin-vector if the observations
    ##     have the same expected values in all years)
    ##   N - optional Nyear-vector for between-year weighting
    ## sigma - specifies the s.d.s of the Xb (either a Nbin-vector or a single
    ##         number).
    ## phi - an n-vector (n = 1 or 2), for the special case that X has the
    ##       correlation matrix of an AR(n) process, containing the
    ##       parameters of this process.
    ##       If length(phi)==1, then cor(Xb,Xc) = phi^|b-c|).
    ## covmat - covariance matrix of X
    ## sepbysex - if T and the data are sexed, then the correlation structure
    ##          is separated by sex (i.e., there is no correlation between
    ##          sexes and the within-sex correlation is the same as for an
    ##          unsexed data set).
    ## sexlag - if T and data are sexed, then the AR(n) correlation structure
    ##          ignores sex and sets lag = |i-j|+1, where i & j index the age
    ##          or length classes in the data.  Ignored if data are not sexed.
    ## robust - if T the robustification described in CompLiklMs.doc is
    ##          implemented
    ## ARMA - if T and if phi is of length 2 than X is assumed to have
    ##          the correlation matrix of an ARMA(1,1) process, where
    ##          phi[1] and phi[2] are the AR and MA parameters, respectively
{
    obsmat <- compdat$obs
    if(!is.matrix(obsmat))obsmat <- matrix(obsmat,nrow=1,
                                           dimnames=list(NULL,names(obsmat)))
    Nbin <- ncol(obsmat)
    Nyear <- nrow(obsmat)
    if(is.in('N',names(compdat))){
        if(length(compdat$N)!=Nyear)stop('Wrong length for compdat$N')
        wts <- sqrt(mean(compdat$N)/compdat$N)
    } else wts <- rep(1,Nyear)
    
    if(length(wts)==1)wts <- rep(1,Nyear)
    expmat <- compdat$exp
    if(!is.matrix(expmat)){
        if(length(expmat)!=Nbin)stop('Wrong size for compdat$exp')
        expmat <- matrix(rep(expmat,Nyear),Nyear,byrow=T,
                         dimnames=list(NULL,names(expmat)))
    } else if(!all(dim(expmat)==dim(obsmat)))stop('Wrong size for compdat$exp')
    if(is.null(covmat))
        covmat <- covmat.logistnorm(sigma,phi,colnames(obsmat),sepbysex,sexlag,
                                    ARMA)
    Kmat <- cbind(diag(Nbin-1),-1)
    Vmat <- Kmat %*% (covmat %*% t(Kmat))
    Vinv <- solve(Vmat)
    negloglik <- 0.5*Nyear*(Nbin-1)*log(2*pi)+sum(log(obsmat))+
        0.5*Nyear*log(det(Vmat))
    if(!all(wts==1))negloglik <- negloglik + (Nbin-1)*sum(log(wts))
     ww <- log(sweep(obsmat[,-Nbin,drop=F],1,obsmat[,Nbin],'/'))- log(sweep(expmat[,-Nbin,drop=F],1,expmat[,Nbin],'/'))
    Vinvdiag <- diag(Vinv)
    for(i in 1:Nyear){
        if(robust){
            ## original proposed robustification
#            uu <- ww[i,] %*% Vinv
#            negloglik <- negloglik-(0.5/(wts[i]^2))*
#                sum(log(exp(-uu*ww[i,])+0.01))
            ## modified robustification
            tmp <- (ww[i,]^2)*Vinvdiag
            negloglik <- negloglik+(0.5/(wts[i]^2))*
                ((ww[i,] %*% Vinv) %*% ww[i,] - sum(tmp) -
                 sum(log(exp(-tmp)+0.01)))
        }
        else negloglik <- negloglik+(0.5/(wts[i]^2))*
            (ww[i,] %*% Vinv) %*% ww[i,]
    }
    return(as.vector(negloglik))
    tmp <- rep(0,Nyear)
    for(i in 1:Nyear){
        tmp[i] <- (0.5/(wts[i]^2))*(ww[i,] %*% Vinv) %*% ww[i,]
    }
    tmp
}
covmat.logistnorm <-
function (sigma,phi,binnam,sepbysex=F,sexlag=F,ARMA=F)
    ## Constructs the covariance matrix of the multivariate normal variate X
    ## from which a logistic-normal variate O is created using the
    ## transformation Ob = exp(Xb)/sumc(exp(Xc)).
    ## For more details, see my function NLLlogistnorm
    ##
    ## sigma - either a single number or a vector of the same length as
    ##         binnam, defining the s.d. of X
    ## phi - a 1- or 2-vector defining the correlation structure of X
    ## binnam - character vector containing the bin names for the composition O.
    ##          The first character of each bin name identifies the sex
    ##          for a sexed composition; for an unsexed composition it
    ##           must be 'X' or a numeral
    ## sepbysex, sexlag - options for sexed compositions (see function
    ##            NLLlogistnorm for details)
    ## ARMA - if T and if phi is of length 2 than X is assumed to have
    ##          the correlation matrix of an ARMA(1,1) process, where
    ##          phi[1] and phi[2] are the AR and MA parameters, respectively
    ##
{
    getrho <- function(phi,kk,ARMA) {
      if(length(phi)==1) {
        phi^(1:(kk-1)) 
      } else if(ARMA) {
        ARMAacf(ar=phi[1],ma=phi[2],lag.max=kk)[-1]
      } else {
        ARMAacf(ar=phi,lag.max=kk)[-1]
      }
    }
    Nbin <- length(binnam)
    if(length(sigma)==1)
      sigma <- rep(sigma,Nbin)
    if(length(sigma)!=Nbin)
      stop('Wrong length for argument sigma')
    if(!is.in(length(phi),1:2))
      stop('Wrong length for phi')
    if(length(phi)==2 & !ARMA){
        if(phi[2]<=(-1) | phi[2]>=(1-abs(phi[1])))
            stop('Invalid value for phi')
    }
    if(all(phi==0)){
        covmat <- diag(sigma^2)
    } else{
        sexed <- !is.in(substring(binnam[1],1,1),c('X',paste(0:9)))
        if(sexed & sepbysex) {
            covmat <- diag(Nbin)
            sexlab <- substring(binnam,1,1)
            for(sx in unique(sexlab)){
                sel <- sexlab==sx
                ik <- sum(sel)
                rhovec <- getrho(phi,ik,ARMA)
                subcov <- diag(ik)
                for(i in 1:(ik-1))
                    subcov[abs(row(subcov)-col(subcov))==i] <- rhovec[i]
                subcov <- sweep(subcov,1,sigma[sel],'*')
                subcov <- sweep(subcov,2,sigma[sel],'*')
                covmat[sel,sel] <- subcov
            }
        }
        else if(sexed & sexlag){
            binlab <- as.numeric(substring(binnam,2))
            ubinlab <- sort(unique(binlab))
            binindex <- match(binlab,ubinlab)
            covmat <- diag(Nbin)
            rhovec <- getrho(phi,length(ubinlab)+1,ARMA)
            binlag <- abs(binlab[row(covmat)]-binlab[col(covmat)])+1
            for(lg in unique(binlag))covmat[binlag==lg] <- rhovec[lg]
            covmat[row(covmat)==col(covmat)] <- 1
            covmat <- sweep(covmat,1,sigma,'*')
            covmat <- sweep(covmat,2,sigma,'*')
        }
        else{
            rhovec <- getrho(phi,Nbin,ARMA)
            covmat <- diag(Nbin)
            for(i in 1:(Nbin-1))
                covmat[abs(row(covmat)-col(covmat))==i] <- rhovec[i]
            covmat <- sweep(covmat,1,sigma,'*')
            covmat <- sweep(covmat,2,sigma,'*')
        }
    }
    covmat
}
Estlogistnorm <-
function(compdat,phi=c(NA,NA),sigbnd=c(0.1,50),
             twosig=F,sepbysex=F,add2comp=0,sexlag=F,robust=F,ARMA=F)
    ## Fits a logistic-normal distribution to some composition data
    ## (observed and expected proportions at age or length from a
    ## stock assessment model).
    ##
    ## The multivariate normal variate X =(X1, X2,..., Xn) that
    ## underlies the logistic-normal always has mean log(compdat$exp),
    ## s.d.(Xi) = sigmai and a correlation matrix defined by lag,
    ## i.e., Cor(Xi,Xj) = rho(|i-j|) for some function rho.
    ##
    ## There are four parameterisations for the sigmai, depending on twosig & N.
    ## sigmai = sigline (if N==NULL) or sigline*sqrt(mean(N)/N) otherwise.
    ## If twosig==F the function estimates a single parameter, sigma and sets
    ## sigline = sigma for all i.  If twosig==T the function estimates
    ## sigma1 and sigma2 and sets
    ## sigline(i) = sigma1+(sigma2-sigma1)*(A(i)-minA)/(maxA-minA),
    ## where A(i) is the age/length associated with the ith bin and
    ## minA & maxA are the minimum & maximum age/length in the data set.
    ##
    ## The parameterisation of rho treats phi (a 2-vector) as the
    ## parameter of an AR(2) process, of which rho is the
    ## auto-correlation function.  To keep the AR(2) process stable
    ## the parameters are constrained to the triangle defined by -1 <=
    ## phi2 <= 1-abs(phi1), which means that -2 < phi1 < 2 [see
    ## Wikipedia article on Autoregressive model]
    ##
    ## The function estimates only those components of phi that are
    ## NA.  Thus, various submodels can be selected by fixing one or
    ## more components of phi.  For example, with phi = c(NA,0) an
    ## AR(1) submodel is fitted; with either phi=c(0,0) the the
    ## one-parameter (sigma) version described by Schnute & Haigh
    ## (2007) is fitted.
    ##
    ## If either compdat$obs or compdat$exp includes any zeroes (not allowed
    ## by the logistic-normal likelihood) both are modified by compressing
    ## the tails to remove all zeroes.
    ##
    ## compdat - the composition data in a list with components
    ##   obs - a Nbin-vector, or a Nyear x Nbin matrix of Nyear observed
    ##     Nbin-vectors.
    ##   exp - the expected compositions in an Nbin-vector or an
    ##     Nyear x Nbin matrix (can use an Nbin-vector if the observations
    ##     have the same expected values in all years)
    ##   N - optional Nyear-vector for between-year weighting
    ## phi - a 2-vector for the AR(2) parameterisation.
    ## sigbnd - to avoid overflow, estimates of sigma
    ##      are constrained lie between sigbnd[1] & sigbnd[2]
    ## twosig - if true, sigma is a linear function of age/length
    ## sepbysex - if T and data set is sexed, then the correlation structure
    ##          is separated by sex (i.e., there is no correlation between
    ##          sexes and the within-sex correlation is the same as for an
    ##          unsexed data set).
    ## add2comp - a small robustifying constant; after zero suppression
    ##        (if any) both compdat$obs and $exp are 'adjusted', which involves
    ##        adding add2comp to all proportions and then normalising
    ##        the proportions to sum to 1 (set to 0 for no robustification)
    ## sexlag - if T and data set is sexed, then the AR(n) correlation structure
    ##          ignores sex and sets lag = |i-j|+1, where i & j index the age
    ##          or length classes in the data.  Ignored if data set not sexed.
    ## robust - if T the robustification described in CompLiklMs.doc is
    ##          implemented
    ## ARMA - if T and if phi is of length 2 than X is assumed to have
    ##          the correlation matrix of an ARMA(1,1) process, where
    ##          phi[1] and phi[2] are the AR and MA parameters, respectively
    ##
{
    Robustify <- function(cdat,const){
        for(nm in c('obs','exp')){
            dat <- cdat[[nm]]
            cdat[[nm]] <-
                if(is.matrix(dat))t(apply(dat,1,function(x)(x+const)/
                                          sum(x+const))) else (
                                              (dat+const)/sum(dat+const))
        }
        cdat
    }
    phi1bnd <- c(-2,2)
    if(!is.matrix(compdat$exp))
      compdat$exp = matrix(compdat$exp ,nrow = 1)
    compdat <- Suppress.zeroes(compdat)
    if(nrow(compdat$exp) == 1)
      compdat$exp = as.numeric(compdat$exp)
      
    ## craig added this, for some reason compress.zeroes was converting my matrix to a data.frame
    if(class(compdat$exp) == "data.frame") {
      compdat$exp = data.matrix(compdat$exp);
    }
    if(class(compdat$obs) == "data.frame") {
      compdat$obs = data.matrix(compdat$obs);
    }    
    if(add2comp>0)compdat <- Robustify(compdat,add2comp)
    Nbin <- ncol(compdat$obs)
    bins <- if(is.in(substring(colnames(compdat$obs)[1],1,1),paste(0:9)))
        as.numeric(colnames(compdat$obs)) else as.numeric(substring(
            colnames(compdat$obs),2))
    minbin <- min(bins)
    maxbin <- max(bins)
    getsig <- function(pars,sigbnd,fullvect){
        if(twosig){
            sigma1 <- inf2posbnd(pars[1],sigbnd)
            sigma2 <- inf2posbnd(pars[2],sigbnd)
            if(fullvect){
                sigma <- sigma1+(sigma2-sigma1)*(bins-minbin)/(maxbin-minbin)
            }
            else sigma <- c(sigma1,sigma2)
        }
        else sigma <- inf2posbnd(pars[1],sigbnd)
        sigma
    }
    if(length(phi)!=2)stop('Illegal length for phi')
    if(all(is.na(phi))){ # 3/4-parameter model: sigma, phi[1], phi[2]
        if(!ARMA){ # correlations are AR(2)
            negloglik <- function(pars,compdat,sigbnd,phi1bnd){
                sigma <- getsig(pars,sigbnd,T)
                phipars <- pars[if(twosig)3:4 else 2:3]
                phi1 <- inf2bnd(phipars[1],phi1bnd)
                absphi1 <- sqrt(phi1^2)
                phi2 <- inf2bnd(phipars[2],c(-1,1-absphi1))
                nll <- NLLlogistnorm(compdat,sigma,c(phi1,phi2),robust=robust,
                                     sepbysex=sepbysex,sexlag=sexlag,ARMA=ARMA)
                nll
            }
            startsig <- rep(posbnd2inf(5,sigbnd),ifelse(twosig,2,1))
            fit <- nlm(negloglik,c(startsig,bnd2inf(0,phi1bnd),0),compdat,sigbnd,phi1bnd)
            sigma <- getsig(fit$estimate,sigbnd,F)
            phipars <- fit$estimate[if(twosig)3:4 else 2:3]
            phi1 <- inf2bnd(phipars[1],phi1bnd)
            absphi1 <- sqrt(phi1^2)
            phi2 <- inf2bnd(phipars[2],c(-1,1-absphi1))
            out <- c(sigma=sigma,phi=c(phi1,phi2),
                     negloglik=fit$minimum)
        }
        else{ # correlations are AR(1)MA(1)
            negloglik <- function(pars,compdat,sigbnd,phi1bnd){
                sigma <- getsig(pars,sigbnd,T)
                phipars <- pars[if(twosig)3:4 else 2:3]
                phi1 <- inf2bnd(phipars[1],phi1bnd)
                absphi1 <- sqrt(phi1^2)
                phi2 <- phipars[2]
                nll <- NLLlogistnorm(compdat,sigma,c(phi1,phi2),robust=robust,
                                     sepbysex=sepbysex,sexlag=sexlag,ARMA=ARMA)
                nll
            }
            phi1bnd <- c(-1,1)
            startsig <- rep(posbnd2inf(5,sigbnd),ifelse(twosig,2,1))
            fit <- nlm(negloglik,c(startsig,bnd2inf(0,phi1bnd),0),
                       compdat,sigbnd,phi1bnd)
            sigma <- getsig(fit$estimate,sigbnd,F)
            phipars <- fit$estimate[if(twosig)3:4 else 2:3]
            phi1 <- inf2bnd(phipars[1],phi1bnd)
            absphi1 <- sqrt(phi1^2)
            phi2 <- phipars[2]
            out <- c(sigma=sigma,phi=c(phi1,phi2),
                     negloglik=fit$minimum)
        }

    }
    else if(all(!is.na(phi))){ # 1/2-parameter model: sigma
        negloglik <- function(pars,compdat,sigbnd){
            sigma <- getsig(pars,sigbnd,T)
            nll <- NLLlogistnorm(compdat,sigma,phi,robust=robust,
                                 sepbysex=sepbysex,sexlag=sexlag,ARMA=ARMA)
            nll
        }
        startsig <- rep(posbnd2inf(5,sigbnd),ifelse(twosig,2,1))
        fit <- nlm(negloglik,startsig,
                   compdat,sigbnd)
        sigma <- getsig(fit$estimate,sigbnd,F)
        out <- c(sigma=sigma,negloglik=fit$minimum)
    }
    else if(is.na(phi[1])){ # 2/3-parameter model: sigma, phi[1]
        if(abs(phi[2])>=1)stop('Illegal value for phi[2]')
        negloglik <- function(pars,compdat,sigbnd){
            sigma <- getsig(pars,sigbnd,T)
            phipar <- pars[ifelse(twosig,3,2)]
            phi1 <- inf2bnd(phipar,c(phi[2]-1,1-phi[2]))
            nll <- NLLlogistnorm(compdat,sigma,c(phi1,phi[2]),robust=robust,
                                 sepbysex=sepbysex,sexlag=sexlag,ARMA=ARMA)
            nll
        }
        startsig <- rep(posbnd2inf(5,sigbnd),ifelse(twosig,2,1))
        fit <- nlm(negloglik,c(startsig,
                               bnd2inf(0,c(phi[2]-1,1-phi[2])),0),
                   compdat,sigbnd)
        sigma <- getsig(fit$estimate,sigbnd,F)
        phipar <- fit$estimate[ifelse(twosig,3,2)]
        phi1 <- inf2bnd(phipar,c(phi[2]-1,1-phi[2]))
        out <- c(sigma=sigma,phi1=phi1,
                 negloglik=fit$minimum)
    }
    else{ # 2/3-parameter model (sigma, phi[2])
        if(abs(phi[1])>=2)stop('Illegal value for phi[1]')
        negloglik <- function(pars,compdat,sigbnd){
            sigma <- getsig(pars,sigbnd,T)
            absphi1 <- sqrt(phi[1]^2)
            phi2 <- -0.5*absphi1+atan(pars[3])*2*(1-absphi1)/pi
            nll <- NLLlogistnorm(compdat,sigma,c(phi[1],phi2),robust=robust,
                                 sepbysex=sepbysex,sexlag=sexlag,ARMA=ARMA)
            nll
        }
        startsig <- rep(posbnd2inf(5,sigbnd),ifelse(twosig,2,1))
        fit <- nlm(negloglik,c(startsig,0,0),compdat,
                   sigbnd)
        sigma <- getsig(fit$estimate,sigbnd,F)
        absphi1 <- sqrt(phi[1]^2)
        phi2 <- -0.5*absphi1+atan(fit$estimate[3])*2*(1-absphi1)/pi
        out <- c(sigma=sigma,phi2=phi2,
                 negloglik=fit$minimum)
    }
    out
}

rlogistnorm <-
function(n,expprop,sigma,phi=0,covmat=NULL,ARMA=F)
    ## Generates n random vectors - in a k x n matrix - from a
    ## logistic-normal distribution with parameters (expprop,sigma,phi) or
    ## (expprop,covmat).
    ##
    ## This function follows Schnute & Haigh (ICES Journal of Marine
    ## Science, 64: 218-233, 2007) in defining the logistic-normal variate O
    ## as being created by applying the logistic transformation
    ## Ob = exp(Xb)/sumc(exp(Xc)) to a multivariate normal variate X, but
    ## it allows X to have any covariance matrix.
    ##
    ## n - number of random vectors to generate
    ## expprop - a k-vector; the expected value of X is log(expprop)
    ## sigma - specifies the s.d.s of the Xb (either a k-vector or a single
    ##         number).
    ## phi - an m-vector (m = 1 or 2), for the special case that X has the
    ##       correlation matrix of an AR(m) process, containing the parameters
    ##       of this process.
    ##       Ignored if covmat is not null.
    ##       If length(phi)==1, then phi = cor(Xb,X[b+1]) and
    ##       cor(Xb,Xc) = rho^|b-c|).
    ## covmat - covariance matrix of X
    ## ARMA - if T and if phi is of length 2 than X is assumed to have
    ##          the correlation matrix of an ARMA(1,1) process, where
    ##          phi[1] and phi[2] are the AR and MA parameters, respectively
    ##
{
    kk <- length(expprop)
    if(is.null(covmat)){
        if(length(sigma)==1)sigma <- rep(sigma,kk)
        if(length(sigma)!=kk)stop('Wrong length for argument sigma')
    }
    if(is.null(covmat) & all(phi==0))
            Xmat <- matrix(rnorm(n*kk)*sigma,kk)+log(expprop)
    else{
        if(is.null(covmat)){
            if(length(phi)==1){
                if(abs(phi)>=1)stop('Invalid value for phi')
            }
            else if(length(phi)==2){
                if(!ARMA & (phi[2]<=(-1) | phi[2]>=(1-abs(phi[1]))))
                    stop('Invalid value for phi')
            }
            else stop('Wrong length for phi')
            rhovec <- (if(!ARMA)ARMAacf(ar=phi,lag.max=kk) else
                       ARMAacf(ar=phi[1],ma=phi[2],lag.max=kk))
            covmat <- diag(kk)
            for(i in 1:(kk-1))
                covmat[abs(row(covmat)-col(covmat))==i] <- rhovec[paste(i)]
            covmat <- sweep(covmat,1,sigma,'*')
            covmat <- sweep(covmat,2,sigma,'*')
        }
        Xmat <- rmultnorm(n,log(expprop),covmat)
    }
    tmp <- exp(Xmat)
    YY <- sweep(tmp,2,apply(tmp,2,sum),'/') # logistic transformation
    YY
}

Sres.logistnorm <-
function(compdat,sigma,phi,centred=T,sepbysex=F,sexlag=F,
                            ARMA=F)
    ## Calculates standardised residuals for composition data
    ## (observed and expected proportions at age or length) on the
    ## assumption that these are logistic-normal with parameters
    ## sigma, phi, & N
    ##
    ## compdat - the composition data in a list with components
    ##   obs - a Nbin-vector, or a Nyear x Nbin matrix of Nyear observed
    ##     Nbin-vectors.
    ##   exp - the expected compositions in an Nbin-vector or an
    ##     Nyear x Nbin matrix (can use an Nbin-vector if the observations
    ##     have the same expected values in all years)
    ##   N - optional Nyear-vector used to weight the sigma
    ##      (sigmay, the sigma for year y is sigma*sqrt(mean(N)/Ny)
    ##      (if there is no component N it is set to N = rep(1,Nyear))
    ## sigma, phi - parameters of the logistic-normal.  sigma may be either
    ##      a single number or a vector of length ncol(obsprop).  phi may be
    ##      either a 1- or 2-vector  (for the AR(1) or AR(2) parameterisations
    ##      of the logistic-normal.
    ## centered - if T, standardised residuals are based on the MVN
    ##      variate Z, where Zi = log(Xi/GM(X)); otherwise they are based
    ##      on the MVN variate Y, where Yi = log(Xi/XD), where D = ncol(obsprop)
    ## sepbysex - if T and data set is sexed, then the correlation structure
    ##          is separated by sex (i.e., there is no correlation between
    ##          sexes and the within-sex correlation is the same as for an
    ##          unsexed data set).
    ## sexlag - if T and data set is sexed, then the AR(n) correlation structure
    ##          ignores sex and sets lag = |i-j|+1, where i & j index the age
    ##          or length classes in the data.  Ignored if data set not sexed.
    ## ARMA - if T and if phi is of length 2 than X is assumed to have
    ##          the correlation matrix of an ARMA(1,1) process, where
    ##          phi[1] and phi[2] are the AR and MA parameters, respectively
    ##
{
    obsprop <- compdat$obs
    expprop <- compdat$exp
    if(!is.matrix(obsprop))obsprop <- matrix(obsprop,1)
    Nbin <- ncol(obsprop)
    Nyear <- nrow(obsprop)
    if(!is.matrix(expprop))expprop <- matrix(rep(expprop,Nyear),
                                             Nyear, byrow=T,
                                             dimnames=list(NULL,names(expprop)))
    if(ncol(expprop)!=Nbin)stop('Wrong size for expprop')
    covmat <- covmat.logistnorm(sigma,phi,colnames(expprop),
                                sepbysex,sexlag,ARMA)
    wts <- if(length(compdat$N)!=0) {
      sqrt(mean(compdat$N)/compdat$N) 
    } else {
      rep(1,Nyear)
    }
    if(length(wts)==1)
      wts <- rep(1,Nyear)
    
    Kmat <- cbind(diag(Nbin-1),-1)
    Vmat <- Kmat %*% (covmat %*% t(Kmat))
    if(centred){
        FF <- cbind(diag(Nbin-1),1)
        Hinv <- solve(diag(Nbin-1)+1)
        FHinv <- t(FF) %*% Hinv
        Gam <- FHinv %*% (Vmat %*% t(FHinv))
        sdmat <- matrix(rep(sqrt(diag(Gam)),Nyear),Nyear,byrow=T)*wts
        gmean <- function(x)prod(x)^(1/length(x))
        sres <- (log(obsprop/apply(obsprop,1,gmean))-
            log(expprop/apply(expprop,1,gmean)))/sdmat
    }
    else{
        sdmat <- matrix(rep(sqrt(diag(Vmat)),Nyear),Nyear,byrow=T)*wts
        sres <- (log(obsprop[,-Nbin]/obsprop[,Nbin])-
            log(expprop[,-Nbin]/expprop[,Nbin]))/sdmat
    }
    sres
}
Suppress.zeroes <-
function(compdat,compress=T)
    ## Removes zeroes in composition data (matrices of observed and expected
    ## proportions at age/length) either by
    ## - compressing the tails (making plus and/or minus groups) or
    ## - combining adjacent columns (if compress==F) or
    ##
    ## If compress==F the name of each column in the new matrices is the
    ## average of the names of the columns from the old matrices that
    ## were summed to make it.  If input matrices have no column names
    ## they are given names paste(1:ncol(compdat$obs)) before any columns
    ## are combined.
    ##
{
    obsprop <- compdat$obs
    expprop <- compdat$exp
    if(!is.matrix(obsprop)){
        obsprop <- matrix(obsprop,1,dimnames=list(NULL,names(obsprop)))
        expprop <- matrix(expprop,1,dimnames=list(NULL,names(expprop)))
        cvt.to.matrix <- T
    } else cvt.to.matrix <- F
    if(is.null(colnames(obsprop))){
        colnames(obsprop) <- paste(1:ncol(obsprop))
    }
    if(all(obsprop!=0) & all(expprop!=0)){
        if(cvt.to.matrix){
            compdat$obs <- obsprop[1,]
            compdat$exp <- expprop[1,]
        }
        return(compdat)
    }
    if(is.in(substring(colnames(obsprop)[1],1,1),c('M','F')))
        stop("Function Suppress.zeroes can't yet deal with sexed data")
    if(compress){ # Remove zeroes using plus/minus groups
        Nbin <- ncol(obsprop)
        has.zero <- apply(obsprop,2,function(y)ifelse(any(y==0),1,0))
        zerocol <- grep(1,c(1,has.zero,1))-1
        indx <- match(max(diff(zerocol)),diff(zerocol))
        col1 <- ifelse(zerocol[indx]>0,zerocol[indx],1)
        col2 <- ifelse(zerocol[indx+1]<=Nbin,zerocol[indx+1],Nbin)
        obsprop[,col1] <- apply(obsprop[,1:col1,drop=F],1,sum)
        if(any(obsprop[,col1]==0)){
            col1 <- col1+1
            obsprop[,col1] <- obsprop[,col1-1]+obsprop[,col1]
        }
        obsprop[,col2] <- apply(obsprop[,col2:Nbin,drop=F],1,sum)
        if(any(obsprop[,col2]==0)){
            col2 <- col2-1
            obsprop[,col2] <- obsprop[,col2+1]+obsprop[,col2]
        }
        expprop[,col1] <- apply(expprop[,1:col1,drop=F],1,sum)
        expprop[,col2] <- apply(expprop[,col2:Nbin,drop=F],1,sum)
        obsprop <- obsprop[,col1:col2,drop=F]
        expprop <- expprop[,col1:col2,drop=F]
        if(any(expprop==0))stop('Bug: expprop has zeroes!!')
    }
    else{ ## remove zeroes by combining adjacent columns
        if(substring(colnames(obsprop)[1],1,1)=='X')
            colnames(obsprop) <- colnames(expprop) <-
                substring(colnames(obsprop),2)
        Nyr <- nrow(obsprop)
        oldbin <- as.numeric(colnames(obsprop))
        Noldbin <- length(oldbin)
        newbin <- oldbin*0
        oldindx <- newindx <- 1
        oldmat <- rbind(obsprop,expprop)
        newmat <- oldmat*0
        while(oldindx<=Noldbin){
            if(all(oldmat[,oldindx]!=0)){ # transfer one bin to new matrix
                newbin[newindx] <- oldbin[oldindx]
                newmat[,newindx] <- oldmat[,oldindx]
                Nbin.in.new.bin <- 1
                newindx <- newindx+1
                oldindx <- oldindx+1
            }
            else{
                if(oldindx<Noldbin){ # start to combine bins
                    bin1 <- oldindx
                    newcol <- oldmat[,oldindx]
                    oldindx <- oldindx+1
                    done <- F
                    while(!done){
                        newcol <- newcol + oldmat[,oldindx]
                        oldindx <- oldindx+1
                        done <- all(newcol!=0) | (oldindx>Noldbin)
                    }
                    bin2 <- oldindx-1
                    if(all(newcol!=0)){ # transfer combined bins to new matrix
                        newbin[newindx] <- round(mean(oldbin[bin1:bin2]),1)
                        newmat[,newindx] <- newcol
                        newindx <- newindx+1
                        Nbin.in.new.bin <- bin2 - bin1 + 1

                    }
                    else{ # combine newcol with previous new bin
                        newmat[,newindx-1] <- newmat[,newindx-1] + newcol
                        newbin[newindx-1] <- round(
                            (Nbin.in.new.bin*newbin[newindx-1]+
                             sum(oldbin[bin1:bin2]))/
                            (Nbin.in.new.bin + bin2 - bin1 + 1),1)
                    }
                }
                else{# combine last old bin with last new bin
                    newmat[,newindx-1] <- newmat[,newindx-1] + oldmat[,Noldbin]
                    newbin[newindx-1] <- round(
                        (Nbin.in.new.bin*newbin[newindx-1]+ oldbin[Noldbin])/
                        (Nbin.in.new.bin + 1),1)
                    oldindx <- oldindx+1
                }
            }
        }
        Nnewbin <- newindx-1
        newmat <- newmat[,1:Nnewbin,drop=F]
        colnames(newmat) <- paste(newbin[1:Nnewbin])
        obsprop <- newmat[1:Nyr,,drop=F]
        expprop <- newmat[-(1:Nyr),,drop=F]
    }
    if(cvt.to.matrix){
        compdat$obs <- obsprop[1,]
        compdat$exp <- expprop[1,]
    } else{
        compdat$obs <- obsprop
        compdat$exp <- expprop
    }
    return(compdat)
}
bnd2inf <-
function(transpar,bnd)
    ## Inverse transformation of inf2bnd
{
    tan((transpar-0.5*(bnd[1]+bnd[2]))*pi/(bnd[2]-bnd[1]))
}
inf2bnd <-
function(par,bnd)
    ## Transforms par from (-inf,inf) to (bnd[1],bnd[2])
{
    atan(par)*(bnd[2]-bnd[1])/pi + 0.5*(bnd[1]+bnd[2])
}
inf2posbnd <-
function(par,bnds)
    ## Transforms par from (-inf,inf) to (bnds[1],bnds[2])
{
    exp(atan(par)*(log(bnds[2]/bnds[1])/pi) + 0.5*log(bnds[1]*bnds[2]))
}
Plcor <-
function(cormat,axis.var='bin',size=1,col.neg='red')
    ## Makes a bubble plot of a correlation matrix, with bubble areas
    ## proportional to the absolute correlation and the sign of the correlation
    ## indicated by whether the bubble is empty (+ve) or filled (-ve).
    ## Default bubble size is such that adjacent bubbles just touch if they
    ## are for correlations 1 or -1.
    ##
    ## axis.var - if 'bin', then the variable shown on the x- and y-axes
    ##            is the bin number, otherwise it is that age or length
    ##            (assuming that the column names of cormat are either
    ##             '3','4', ... or 'X3','X4', ... for ages 3,4,...)
    ## size - default bubble radii are multiplied by size
    ## col.neg - colour of bubbles for -ve correlations
{
   if(axis.var=='bin') x <- 1:ncol(cormat)
   else x <- (if(substring(colnames(cormat)[1],1,1)=='X')
              as.numeric(substring(colnames(cormat),2)) else
              as.numeric(colnames(cormat)))
   delx <- mean(diff(x))
   plrng <- range(x)+c(-0.5,0.5)*delx
   plot(0,0,type='n',xlim=plrng,xaxs='i',ylim=plrng,yaxs='i',xlab='',ylab='')
   radii <- size*delx*sqrt(abs(cormat))/2
   sel <- row(cormat)!=col(cormat) & cormat>0
   if(any(sel))
        symbols(x[row(cormat)][sel],x[col(cormat)][sel],
                circles=radii[sel],add=T,inches=F)
    sel <- row(cormat)!=col(cormat) & cormat<0
    if(any(sel))
        symbols(x[row(cormat)][sel],x[col(cormat)][sel],
                circles=radii[sel],add=T,
                fg=col.neg,bg=col.neg,inches=F)
    invisible()
}
posbnd2inf <-
function(transpar,bnds)
    ## Inverse transformation of inf2posbnd
{
    tan((log(transpar)-0.5*log(bnds[1]*bnds[2]))*pi/log(bnds[2]/bnds[1]))
}
unlapply <-
function (lst, fun, ...)
unlist(lapply(lst, fun, ...))
is.in <-
function (x, y)
!is.na(match(x, y))
rmultnorm <-
function (n, mu, sigma, tol = 1e-07)
{
    p <- ncol(sigma)
    if (length(mu) != p)
        stop("mu vector is the wrong length")
    if (max(abs(sigma - t(sigma))) > tol)
        stop("sigma not symmetric")
    vs <- svd(sigma)
    vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
    ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
    ans <- sweep(ans, 2, mu, "+")
    dimnames(ans) <- list(NULL, dimnames(sigma)[[2]])
    t(ans)
}
Plcor <-
function(cormat,axis.var='bin',size=1,col.neg='red')
    ## Makes a bubble plot of a correlation matrix, with bubble areas
    ## proportional to the absolute correlation and the sign of the correlation
    ## indicated by whether the bubble is empty (+ve) or filled (-ve).
    ## Default bubble size is such that adjacent bubbles just touch if they
    ## are for correlations 1 or -1.
    ##
    ## axis.var - if 'bin', then the variable shown on the x- and y-axes
    ##            is the bin number, otherwise it is that age or length
    ##            (assuming that the column names of cormat are either
    ##             '3','4', ... or 'X3','X4', ... for ages 3,4,...)
    ## size - default bubble radii are multiplied by size
    ## col.neg - colour of bubbles for -ve correlations
{
   if(axis.var=='bin') x <- 1:ncol(cormat)
   else x <- (if(substring(colnames(cormat)[1],1,1)=='X')
              as.numeric(substring(colnames(cormat),2)) else
              as.numeric(colnames(cormat)))
   delx <- mean(diff(x))
   plrng <- range(x)+c(-0.5,0.5)*delx
   plot(0,0,type='n',xlim=plrng,xaxs='i',ylim=plrng,yaxs='i',xlab='',ylab='')
   radii <- size*delx*sqrt(abs(cormat))/2
   sel <- row(cormat)!=col(cormat) & cormat>0
   if(any(sel))
        symbols(x[row(cormat)][sel],x[col(cormat)][sel],
                circles=radii[sel],add=T,inches=F)
    sel <- row(cormat)!=col(cormat) & cormat<0
    if(any(sel))
        symbols(x[row(cormat)][sel],x[col(cormat)][sel],
                circles=radii[sel],add=T,
                fg=col.neg,bg=col.neg,inches=F)
    invisible()
}
Simcomp <-
function(compdat,dist,parvec,nsim=1000)
    ## Simulates multiple copies of a composition data set according to
    ## a specified (logistic-normal or Dirichlet) distributions with
    ## parameters in compdat$exp and parvec
    ##
    ## The simulated data is in the same format as compdat, but will have
    ## nsim x Nyr years of data, if compdat has Nyr years of data.
    ##
    ## compdat - the composition data in a list with components obs, exp, & N
    ## dist - the assumed distribution (one of 'Dir','LN1','LN2','LN3','LN3a')
    ## parvec - a named vector containing the weighting parameters, which
    ##          are named as in bestfits, viz,
    ##          for Dir: 'Dirn'
    ##          for LN1: 'sig1'
    ##          for LN2: 'sig2','phi'
    ##          for LN3: 'sig3','phi1','phi2'
    ##          for LN3a: 'sigma','phi1','phi2'
    ## nsim - number of copies to simulate
    ##
{
    if(!is.in(dist,c('Dir','LN1','LN2','LN3','LN3a')))
        stop('Invalid value for argument dist')
    Nyr <- nrow(compdat$obs)
    bin.names <- colnames(compdat$obs)
    Nbin <- length(bin.names)
    N <- compdat$N
    if(length(N)==1)N <- rep(N,Nyr)
    expmat <- obsmat <- matrix(0,nsim*Nyr,Nbin,dimnames=list(NULL,bin.names))
    for(j in 1:ncol(expmat))expmat[,j] <- rep(compdat$exp[,j],rep(nsim,Nyr))
    if(dist=='Dir'){
        wts <- N/mean(N)
        Dirn <- parvec['Dirn']
    }
    else{
        wts <- sqrt(mean(N)/N)
        siglab <- c(LN1='sig1',LN2='sig2',LN3='sig3',LN3a='sigma')[dist]
        sigma <- parvec[siglab]
        phi <- if(dist=='LN1')0 else if(dist=='LN2')parvec['phi'] else (
                                        parvec[c('phi1','phi2')])
        ARMA <- dist=='LN3a'
    }
    for(j in 1:Nyr){
        obsmat[(1:nsim)+nsim*(j-1),] <-
            if(dist=='Dir')t(rdirichlet(nsim,compdat$exp[j,],Dirn*wts[j]))
            else t(rlogistnorm(nsim,compdat$exp[j,],sigma*wts[j],phi,
                               ARMA=ARMA))
    }
    out <- list(obs=obsmat,exp=expmat,N=rep(N,rep(nsim,Nyr)))
}
