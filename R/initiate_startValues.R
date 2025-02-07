initiate_startValues <-
function(ppp, K, model="MLGCP", gamma=NULL, m=NULL, sigSq=NULL, phi=NULL, gamma0=NULL, alpha=NULL, phi0=NULL)
{
    ## MLGCP
    
    if(model == "MLGCP")
    {
        cat(paste("Starting values are initiated for ", model, sep = ""), cat("\n"))
        
        J = ppp$J
        J = ppp$J
        nC.ext = ppp$nC.ext
        n.ji.ext = ppp$n.ji.ext
        obs.inx <- ppp$obs.inx
        A <- ppp$A
        
        ini <- PLN(n.ji.ext[obs.inx,]~ 1 + offset(log(rep(A, length(obs.inx)))))
        
        if(is.null(gamma))
        {
            gamma <- matrix(rnorm(nC.ext * J, 0, 1), nC.ext, J)
        }
        if(is.null(m))
        {
            m <- as.vector(ini$model_par$Theta) * runif(1, 0.99, 1.01)
        }
        if(is.null(sigSq))
        {
            sigSq <- as.vector(diag(ini$model_par$Sigma)) * runif(1, 0.99, 1.01)
        }
        
        if(is.null(phi))
        {
            phi <- rep(0.05 * min(diff(ppp$W$xrange), diff(ppp$W$yrange)), J) * runif(1, 0.99, 1.01)
        }
        
        if(K == 0)
        {
            gamma0 <- NULL
            alpha <- NULL
            phi0 <- NULL
        }else
        {
            gamma0 <- matrix(rnorm(nC.ext * K, 0, 0.0001), nC.ext, K)
            alpha <- matrix(NA, J, K)
            for(j in 1:J)
            {
                alpha[j,] <- rep(0.0001, K) * runif(1, 0.99, 1.01)
            }
            phi0 <- rep(0.05 * min(diff(ppp$W$xrange), diff(ppp$W$yrange)), K) * runif(1, 0.99, 1.01)
        }
        
        ret <- list(gamma=gamma, m=m, sigSq=sigSq, phi=phi, gamma0=gamma0, alpha=alpha, phi0=phi0, model=model)        
    }
    
    ##
    return(ret)
}
