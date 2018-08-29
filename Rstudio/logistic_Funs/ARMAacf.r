    p <- length(ar)
    q <- length(ma)
    if (!p && !q) 
        stop("empty model supplied")
    r <- max(p, q + 1)
        if (r > 1) {
        ## if we have a ma coeffecient
            if (r > p) {
            
                ar <- c(ar, rep(0, r - p))
                p <- r
            }
            A <- matrix(0, p + 1L, 2L * p + 1L)
            ind <- as.matrix(expand.grid(1L:(p + 1), 1L:(p + 
                1)))[, 2L:1L]
            ind[, 2] <- ind[, 1L] + ind[, 2L] - 1L
            A[ind] <- c(1, -ar)
            A[, 1:p] <- A[, 1:p] + A[, (2 * p + 1):(p + 2)]
            rhs <- c(1, rep(0, p))
            
            if (q > 0) {
                psi <- c(1, ARMAtoMA(ar, ma, q))
                theta <- c(1, ma, rep(0, q + 1L))
                for (k in 1 + 0:q) rhs[k] <- sum(psi * theta[k + 0:q])
            }
            ind <- (p + 1):1
            Acf <- solve(A[ind, ind], rhs)
            Acf <- Acf[-1L]/Acf[1L]
        } else  ## we have a single ar coeffiecient
          Acf <- ar
        
        if (lag.max > p) {
            xx <- rep(0, lag.max - p)
            Acf <- c(Acf, filter(xx, ar, "recursive", init = rev(Acf)))
        }
        Acf <- c(1, Acf[1L:lag.max])
    }
    names(Acf) <- 0:lag.max
    if (pacf) 
        drop(.Call(C_pacf1, Acf, lag.max))
    else Acf
}


