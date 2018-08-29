    method <- match.arg(method)
    x <- as.ts(x)
    storage.mode(x) <- "double"
    xtsp <- tsp(x)
    n <- as.integer(NROW(x))
    if (is.na(n)) 
        stop("invalid value of nrow(x)", domain = NA)
    nser <- NCOL(x)
    filter <- as.double(filter)
    nfilt <- as.integer(length(filter))
    if (is.na(n)) 
        stop("invalid value of length(filter)", domain = NA)
    if (anyNA(filter)) 
        stop("missing values in 'filter'")
  ## always rec ursive

        if (missing(init)) {
            init <- matrix(0, nfilt, nser)
        } else {
            ni <- NROW(init)
            if (ni != nfilt) 
                stop("length of 'init' must equal length of 'filter'")
            if (NCOL(init) != 1L && NCOL(init) != nser) {
                stop(sprintf(ngettext(nser, "'init' must have %d column", 
                  "'init' must have 1 or %d columns", domain = "R-stats"), 
                  nser), domain = NA)
            }
            if (!is.matrix(init)) 
                dim(init) <- c(nfilt, nser)
        }
        ind <- seq_len(nfilt)
        if (is.matrix(x)) {
            y <- matrix(NA, n, nser)
            for (i in seq_len(nser)) 
            y[, i] <- .Call(C_rfilter,  x[, i], filter, c(rev(init[, i]), double(n)))[-ind]
        }
        else y <- .Call(C_rfilter, x, filter, c(rev(init[, 1L]), double(n)))[-ind]
   }
    tsp(y) <- xtsp
    class(y) <- if (nser > 1L) 
        c("mts", "ts")
    else "ts"
    y

 

ar = c(0.3508772, 0.5001754)
Acf = c(0.7020007, 0.7464914)

nBin = 16
xx = rep(0.0, nBin)

c(Acf,filter(xx, ar, "recursive", init = rev(Acf)))
c(Acf,Rec_filter(ar, xx ,Acf))
Cor = c(Acf,Rec_filter(ar, xx ,Acf))
## The result = 
ARMAacf(ar=phi,lag.max=kk)[-1]
P = vector()

## Recreating the uni_pacf issue.
uni_pcf = function(Cor, nBin) {
  P = vector();
  w = v = rep(0.0, nBin);
  w[1] = P[1] = Cor[2];
  for (ll in 2:(nBin - 1)) {
    a = Cor[ll + 1];
    b = 1.0;
    for (i in 1:(ll - 1)) {
      a = a - w[i] * Cor[ll - i];
      b = b - w[i] * Cor[i + 1];
    }
    P[ll] = c = a / b;
       if (ll == nlag) {
         next;
    }
    w[ll] = c;
    for (i in 1:(ll - 1)) {
      v[ll-i] = w[i];
    }
    
    for (i in 1:(ll - 1)) {
      w[i] = w[i] - c * v[i];
    }    
  }
  P
}

uni_pcf(Cor,nBin)


ndx = 1:3
rev_ndx = rev(3:1)
store_mat = matrix(0.0, nrow = 3, ncol = 3)

for ( i in ndx) {
  for ( j in ndx) {
    store_mat[rev_ndx[i],rev_ndx[j]] = A[rev_ndx[i],rev_ndx[j]]
  }
}


## now lets look at the when we have a ma coeffiecient
phi = c(0.3,0.5)
kk = 15
round(ARMAacf(ar=phi[1],ma=phi[2],lag.max=kk)[-1],4)

ar=phi[1]
ma=phi[2]
lag.max=kk



## test code for solving the inverse with known B
x = c(1 , 2 , 3, 4 , 5 , 6 ,7 , 8, 10)

A_mat = matrix(x, nrow = 3, byrow = T)
solve(A_mat, c(3,3,4))

