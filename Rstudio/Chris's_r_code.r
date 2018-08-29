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
  ww <- log(sweep(obsmat[,-Nbin,drop=F],1,obsmat[,Nbin],'/')) - log(sweep(expmat[,-Nbin,drop=F],1,expmat[,Nbin],'/'))
  Vinvdiag <- diag(Vinv)
  for(i in 1:Nyear){
    if(robust){
      ## original proposed robustification
      ##            uu <- ww[i,] %*% Vinv
      ##            negloglik <- negloglik-(0.5/(wts[i]^2))*
      ##                sum(log(exp(-uu*ww[i,])+0.01))
      ## modified robustification
      tmp <- (ww[i,]^2)*Vinvdiag
      negloglik <- negloglik+(0.5/(wts[i]^2))*
        ((ww[i,] %*% Vinv) %*% ww[i,] - sum(tmp) - sum(log(exp(-tmp)+0.01)))
    }
    else negloglik <- negloglik+(0.5/(wts[i]^2))* (ww[i,] %*% Vinv) %*% ww[i,]
  }
  return(as.vector(negloglik))
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
        Xmat <- rmvnorm(n,log(expprop),covmat)
    }
    tmp <- exp(Xmat)
    YY <- sweep(tmp,2,apply(tmp,1,sum),'/') # logistic transformation
    YY
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
    getrho <- function(phi,kk,ARMA)if(length(phi)==1)
        phi^(1:(kk-1)) else if(ARMA) { ARMAacf(ar=phi[1],ma=phi[2],lag.max=kk)[-1]
        } else ARMAacf(ar=phi,lag.max=kk)[-1]
    Nbin <- length(binnam)
    if(length(sigma)==1)sigma <- rep(sigma,Nbin)
    if(length(sigma)!=Nbin)stop('Wrong length for argument sigma')
    if(!is.in(length(phi),1:2))stop('Wrong length for phi')
    if(length(phi)==2 & !ARMA){
        if(phi[2]<=(-1) | phi[2]>=(1-abs(phi[1])))
            stop('Invalid value for phi')
    }
    if(all(phi==0)){
        covmat <- diag(sigma^2)
    }
    else{
        sexed <- !is.in(substring(binnam[1],1,1),c('X',paste(0:9)))
        if(sexed & sepbysex){
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
is.in <-
structure(function (x, y) 
!is.na(match(x, y)), source = c("function(x, y)", "!is.na(match(x, y))"
))
