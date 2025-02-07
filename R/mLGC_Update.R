
###################################
######     updating gammmas using HMC
###################################
Update_GAM <- function(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C.eigen, norm.gam, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, L.gam, eps.gam, M.gam.var, accept.gam, accept.gam.100, n.GAM, n.nan.GAM, MM, numReps)
{
    J <- length(m)
    if(is.null(phi0))
    {
        K <- 0
    }else
    {
        K <- length(phi0)
    }
    M.e <- dim(D.e.base)[1]
    N.e <- dim(D.e.base)[2]
    nC.e <- M.e * N.e
    
    nan_gam <- rep(0, J)
    
    for(jj in 1:J)
    {
        gamma.ini <- gamma.star <- gamma[,jj]
        p.star = rnorm(nC.e, 0, sqrt(M.gam.var))
        
        Delta.star <- sqrt(sigSq[jj])*as.vector(Re(fft(sqrt.C.eigen[[jj]] * norm.aVec[,jj], TRUE)/sqrt(nC.e))) - gamma.star
        
        p.ini <- p.star + 0.5 * eps.gam[jj] * Delta.star
        for(l in 1:L.gam)
        {
            gamma.prop <- gamma.ini + eps.gam[jj] * rep(1/M.gam.var, nC.e) * p.ini
            
            norm.gam.prop <- fft(matrix(gamma.prop, M.e, N.e))/sqrt(nC.e)
            sqrt.C_gam.prop <- as.vector(Re(fft(sqrt.C.eigen[[jj]] * norm.gam.prop, TRUE)/sqrt(nC.e)))
            
            Y.e.prop <- m[jj] * obs.ind + sqrt(sigSq[jj])*sqrt.C_gam.prop
            if(K >0)
            {
                for(kk in 1:K)
                {
                    Y.e.prop <- Y.e.prop + alpha[jj, kk]*sqrt.C0_gam0[[kk]]
                }
            }
            Y.e.prop[Y.e.prop> 700] <- 700
            aVec.prop <- n.e[,jj] -  A * obs.ind * exp(Y.e.prop)
            norm.aVec.prop <- fft(matrix(aVec.prop, M.e, N.e))/sqrt(nC.e)
            Delta.prop <- sqrt(sigSq[jj]) * as.vector(Re(fft(sqrt.C.eigen[[jj]] * norm.aVec.prop, TRUE)/sqrt(nC.e))) - gamma.prop
            
            if(l < L.gam)
            {
                p.prop <- p.ini + eps.gam[jj] * Delta.prop
                if(sum(is.nan(p.prop))==0 & sum(is.infinite(p.prop)) == 0)
                {
                    gamma.ini <- gamma.prop
                    p.ini <- p.prop
                    l <- L.gam
                }
            }else
            {
                p.prop <- p.ini + 0.5 * eps.gam[jj] * Delta.prop
                if(sum(is.nan(p.prop)) > 0 | sum(is.infinite(p.prop)) > 0)
                {
                    p.prop <- p.ini
                    gamma.prop <- gamma.ini
                }
            }
        }
        
        gamma.prop.mat <- gamma
        gamma.prop.mat[,jj] <- gamma.prop
        Y.e.prop.mat <- Y.e
        Y.e.prop.mat[,jj] <- Y.e.prop
        
        # accept/reject
        
        U.prop <- -sum(n.e[,jj] * sqrt(sigSq[jj])*sqrt.C_gam.prop - A * obs.ind * exp(Y.e.prop)) - sum(dnorm(gamma.prop, 0, 1, log=T))
        K.prop <- 0.5*(sum(as.vector(p.prop)^2/M.gam.var))
        
        if(is.nan(U.prop) | is.nan(K.prop) | is.infinite(U.prop) | is.infinite(K.prop))
        {
            nan_gam[jj] <- 1
            
            accept.gam.100[[jj]] <- c(accept.gam.100[[jj]], 0)
            if(n.GAM %% 10 == 0 & n.GAM >= 100)
            {
                if(sum(accept.gam.100[[jj]])/100 < 0.60)
                {
                    eps.gam[jj] <- eps.gam[jj] * 0.9
                }else if(sum(accept.gam.100[[jj]])/100 > 0.70)
                {
                    eps.gam[jj] <- eps.gam[jj] * 1.1
                }
                accept.gam.100[[jj]] <- accept.gam.100[[jj]][-c(1:10)]
            }
            
        }else
        {
            U.star <- -sum(n.e[,jj] * sqrt(sigSq[jj])*sqrt.C_gam[[jj]] - A * obs.ind * exp(Y.e[,jj])) - sum(dnorm(gamma[,jj], 0, 1, log=T))
            K.star <- 0.5*(sum(as.vector(p.star)^2/M.gam.var))
            
            logR = -(U.prop + K.prop) + U.star + K.star
            
            u = log(runif(1)) < logR
            
            if(u == 1)
            {
                gamma <- gamma.prop.mat
                norm.gam[[jj]] <- norm.gam.prop
                sqrt.C_gam[[jj]] <-  sqrt.C_gam.prop
                
                Y.e <- Y.e.prop.mat
                aVec[,jj] <- aVec.prop
                norm.aVec[,jj] <- norm.aVec.prop
            }
            
            accept.gam[jj] <- accept.gam[jj] + u
            
            accept.gam.100[[jj]] <- c(accept.gam.100[[jj]], u)
            if(n.GAM %% 10 == 0 & n.GAM >= 100)
            {
                if(sum(accept.gam.100[[jj]])/100 < 0.60)
                {
                    eps.gam[jj] <- eps.gam[jj] * 0.9
                }else if(sum(accept.gam.100[[jj]])/100 > 0.70)
                {
                    eps.gam[jj] <- eps.gam[jj] * 1.1
                }
                accept.gam.100[[jj]] <- accept.gam.100[[jj]][-c(1:10)]
            }
        }
    }
    
    n.nan.GAM <- n.nan.GAM + nan_gam
    
    return(list(gamma=gamma, norm.gam=norm.gam, sqrt.C_gam=sqrt.C_gam, Y.e=Y.e, aVec=aVec, norm.aVec=norm.aVec, accept.gam=accept.gam, accept.gam.100 =accept.gam.100, eps.gam=eps.gam, n.nan.GAM=n.nan.GAM))
}


#########################
######     updating phi using RW-MH
#########################
Update_PHI <- function(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, C.e.base, sqrt.C.eigen, norm.gam, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, l.phi, u.phi, propVar.phi, accept.phi, n.PHI, n.nan.PHI)
{
    J <- length(phi)
    if(is.null(phi0))
    {
        K <- 0
    }else
    {
        K <- length(phi0)
    }
    
    M.e <- dim(D.e.base)[1]
    N.e <- dim(D.e.base)[2]
    nC.e <- M.e * N.e
    
    nan_phi <- rep(0, J)
    
    for(jj in 1:J)
    {
        theta <- log(phi[jj])
        
        log.lh.ini <- sum(n.e[, jj] * sqrt(sigSq[jj]) * sqrt.C_gam[[jj]] -  A * obs.ind * exp(Y.e[,jj])) + theta
        
        # var = wd^2 / 12 -> wd = sqrt(12 * var)
        wd.phi = sqrt(12*propVar.phi)
        
        if(abs(theta-l.phi) < wd.phi/2)
        {
            prop.lb <- l.phi
            prop.ub <- l.phi + wd.phi
        }else if(abs(theta-u.phi) < wd.phi/2)
        {
            prop.lb <- u.phi - wd.phi
            prop.ub <- u.phi
        }else
        {
            prop.lb <- theta - wd.phi/2
            prop.ub <- theta + wd.phi/2
        }
        
        theta.prop <- runif(1, prop.lb, prop.ub)
        phi.prop <- exp(theta.prop)
        
        if(theta.prop > l.phi & theta.prop < u.phi)
        {
            C.e.base.prop <- corfunc(D.e.base, phi.prop)
            sqrt.C.eigen.prop <-  sqrt(Re(fft(C.e.base.prop, TRUE)))
            
            if(sum(is.nan(sqrt.C.eigen.prop)) != 0)
            {
                nan_phi[jj] <- 1
                logR <- -Inf
            }else
            {
                sqrt.C_gam.prop <- as.vector(Re(fft(sqrt.C.eigen.prop * norm.gam[[jj]], TRUE)/sqrt(nC.e)))
                
                Y.e.prop <- Y.e[,jj]
                
                Y.e.prop <- m[jj] * obs.ind + sqrt(sigSq[jj])*sqrt.C_gam.prop
                if(K >0)
                {
                    for(kk in 1:K)
                    {
                        Y.e.prop <- Y.e.prop + alpha[jj, kk]*sqrt.C0_gam0[[kk]]
                    }
                }
                Y.e.prop[Y.e.prop> 700] <- 700
                
                log.lh.prop <- sum(n.e[, jj] * sqrt(sigSq[jj]) * sqrt.C_gam.prop -  A * obs.ind * exp(Y.e.prop)) + theta.prop
                
                if(log.lh.prop == -Inf)
                {
                    logR <- -Inf
                }else
                {
                    log.prior.ini  <- 0
                    log.prior.prop  <- 0
                    
                    log.prop.prop  <- dunif(theta.prop, prop.lb, prop.ub, log = T)
                    log.prop.ini   <- dunif(theta, prop.lb, prop.ub, log = T)
                }
                logR  <- log.lh.prop - log.lh.ini + log.prior.prop - log.prior.ini + log.prop.ini - log.prop.prop
            }
            
            u = log(runif(1)) < logR
            
            if(u == 1)
            {
                phi[jj] <- phi.prop
                C.e.base[[jj]] <- C.e.base.prop
                sqrt.C.eigen[[jj]] <- sqrt.C.eigen.prop
                sqrt.C_gam[[jj]] <-  sqrt.C_gam.prop
                
                Y.e[,jj] <- Y.e.prop
                aVec[,jj] <- n.e[,jj] -  A * obs.ind * exp(Y.e.prop)
                norm.aVec[,jj] <- fft(matrix(aVec[,jj], M.e, N.e))/sqrt(nC.e)
            }
            accept.phi[jj] <- accept.phi[jj] + u
            
            
            
        }
    }
    n.nan.PHI <- n.nan.PHI + nan_phi
    
    return(list(phi=phi, Y.e=Y.e, C.e.base=C.e.base, sqrt.C.eigen=sqrt.C.eigen, sqrt.C_gam=sqrt.C_gam, aVec=aVec, norm.aVec=norm.aVec, accept.phi=accept.phi, n.nan.PHI=n.nan.PHI))
}



###################################
######     updating gammma0... (reparameterized log-intensity)
###################################
Update_GAM0 <- function(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C_gam, C0.e.base, sqrt.C0.eigen, norm.gam0, sqrt.C0_gam0, Y.e, aVec, norm.aVec, L.gam0, eps.gam0, M.gam0.var, accept.gam0, accept.gam0.100, n.GAM0, n.nan.GAM0, MM, numReps)
{
    J <- length(m)
    if(is.null(phi0))
    {
        K <- 0
    }else
    {
        K <- length(phi0)
    }
    M.e <- dim(D.e.base)[1]
    N.e <- dim(D.e.base)[2]
    nC.e <- M.e * N.e
    
    nan_gam0 <- rep(0, K)
    
    for(kk in 1:K)
    {
        gamma0.ini <- gamma0.star <- gamma0[,kk]
        p.star = rnorm(nC.e, 0, sqrt(M.gam0.var))
        Delta.star <- rep(0, nC.e)
        for(jj in 1:J)
        {
            Delta.star <- Delta.star + alpha[jj,kk] * as.vector(Re(fft(sqrt.C0.eigen[[kk]] * norm.aVec[,jj], TRUE)/sqrt(nC.e)))
        }
        Delta.star <- Delta.star - gamma0.star
        
        p.ini <- p.star + 0.5 * eps.gam0[kk] * Delta.star
        
        Y.e.prop <- matrix(NA, nC.e, J)
        aVec.prop <- matrix(NA, nC.e, J)
        norm.aVec.prop <- matrix(NA, nC.e, J)
        
        for(l in 1:L.gam0)
        {
            gamma0.prop <- gamma0.ini + eps.gam0[kk] * rep(1/M.gam0.var, nC.e) * p.ini
            norm.gam0.prop <- fft(matrix(gamma0.prop, M.e, N.e))/sqrt(nC.e)
            sqrt.C0_gam0.prop <- as.vector(Re(fft(sqrt.C0.eigen[[kk]] * norm.gam0.prop, TRUE)/sqrt(nC.e)))
            
            for(jj in 1:J)
            {
                Y.e.prop[,jj] <- m[jj] * obs.ind + sqrt(sigSq[jj])*sqrt.C_gam[[jj]]
                for(ll in 1:K)
                {
                    if(ll == kk)
                    {
                        Y.e.prop[,jj] <- Y.e.prop[,jj] + alpha[jj, ll]*sqrt.C0_gam0.prop
                    }else
                    {
                        Y.e.prop[,jj] <- Y.e.prop[,jj] + alpha[jj, ll]*sqrt.C0_gam0[[ll]]
                    }
                }
                Y.e.prop[Y.e.prop[,jj]> 700, jj] <- 700
                aVec.prop[,jj] <- n.e[,jj] - A * obs.ind * exp(Y.e.prop[,jj])
                norm.aVec.prop[,jj] <- fft(matrix(aVec.prop[,jj], M.e, N.e))/sqrt(nC.e)
            }
            
            Delta.prop <- rep(0, nC.e)
            for(jj in 1:J)
            {
                Delta.prop <- Delta.prop + alpha[jj,kk] * as.vector(Re(fft(sqrt.C0.eigen[[kk]] * norm.aVec.prop[,jj], TRUE)/sqrt(nC.e)))
            }
            Delta.prop <- Delta.prop - gamma0.prop
            
            if(l < L.gam0)
            {
                p.prop <- p.ini + eps.gam0[kk] * Delta.prop
                
                if(sum(is.nan(p.prop))==0 & sum(is.infinite(p.prop)) == 0)
                {
                    gamma0.ini <- gamma0.prop
                    p.ini <- p.prop
                    l <- L.gam0
                }
            }else
            {
                p.prop <- p.ini + 0.5 * eps.gam0[kk] * Delta.prop
                if(sum(is.nan(p.prop)) > 0 | sum(is.infinite(p.prop)) > 0)
                {
                    p.prop <- p.ini
                    gamma0.prop <- gamma0.ini
                }
            }
        }
        
        gamma0.prop.mat <- gamma0
        gamma0.prop.mat[,kk] <- gamma0.prop
        
        # accept/reject
        U.prop <- -sum(dnorm(gamma0.prop, 0, 1, log=T))
        
        for(jj in 1:J)
        {
            U.prop <- U.prop - sum(n.e[, jj] * alpha[jj, kk] * sqrt.C0_gam0.prop -  A * obs.ind * exp(Y.e.prop[,jj]))
        }
        
        K.prop <- 0.5*(sum(as.vector(p.prop)^2/M.gam0.var))
        
        if(is.nan(U.prop) | is.nan(K.prop) | is.infinite(U.prop) | is.infinite(K.prop))
        {
            nan_gam0[kk] <- 1
            
            accept.gam0.100[[kk]] <- c(accept.gam0.100[[kk]], 0)
            if(n.GAM0 %% 10 == 0 & n.GAM0 >= 100)
            {
                if(sum(accept.gam0.100[[kk]])/100 < 0.60)
                {
                    eps.gam0[kk] <- eps.gam0[kk] * 0.9
                }else if(sum(accept.gam0.100[[kk]])/100 > 0.70)
                {
                    eps.gam0[kk] <- eps.gam0[kk] * 1.1
                }
                accept.gam0.100[[kk]] <- accept.gam0.100[[kk]][-c(1:10)]
            }
        }else
        {
            U.star <- -sum(dnorm(gamma0[,kk], 0, 1, log=T))
            for(jj in 1:J)
            {
                U.star <- U.star - sum(n.e[, jj] * alpha[jj, kk] * sqrt.C0_gam0[[kk]] -  A * obs.ind * exp(Y.e[,jj]))
            }
            
            K.star <- 0.5*(sum(as.vector(p.star)^2/M.gam0.var))
            
            logR = -(U.prop + K.prop) + U.star + K.star
            
            u = log(runif(1)) < logR
            
            if(u == 1)
            {
                gamma0 <- gamma0.prop.mat
                norm.gam0[[kk]] <- norm.gam0.prop
                sqrt.C0_gam0[[kk]] <-  sqrt.C0_gam0.prop
                
                Y.e <- Y.e.prop
                aVec <- aVec.prop
                norm.aVec <- norm.aVec.prop
            }
            
            accept.gam0[kk] <- accept.gam0[kk] + u
            
            accept.gam0.100[[kk]] <- c(accept.gam0.100[[kk]], u)
            if(n.GAM0 %% 10 == 0 & n.GAM0 >= 100)
            {
                if(sum(accept.gam0.100[[kk]])/100 < 0.60)
                {
                    eps.gam0[kk] <- eps.gam0[kk] * 0.9
                }else if(sum(accept.gam0.100[[kk]])/100 > 0.70)
                {
                    eps.gam0[kk] <- eps.gam0[kk] * 1.1
                }
                accept.gam0.100[[kk]] <- accept.gam0.100[[kk]][-c(1:10)]
            }
        }
        
    }
    
    n.nan.GAM0 <- n.nan.GAM0 + nan_gam0
    
    return(list(gamma0=gamma0, norm.gam0=norm.gam0, sqrt.C0_gam0=sqrt.C0_gam0, Y.e=Y.e, aVec=aVec, norm.aVec=norm.aVec, accept.gam0=accept.gam0, accept.gam0.100=accept.gam0.100, eps.gam0=eps.gam0, n.nan.GAM0=n.nan.GAM0))
}



#########################
######     updating phi0 -- common cor param
#########################
Update_PHI0 <- function(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, C0.e.base, sqrt.C0.eigen, norm.gam0, sqrt.C0_gam0, Y.e, aVec, norm.aVec, l.phi0, u.phi0, propVar.phi0, accept.phi0, n.PHI0, n.nan.PHI0)
{
    J <- dim(alpha)[1]
    K <- dim(alpha)[2]
    
    M.e <- dim(D.e.base)[1]
    N.e <- dim(D.e.base)[2]
    nC.e <- M.e * N.e
    
    nan_phi0 <- rep(0, K)
    
    for(kk in 1:K)
    {
        theta <- log(phi0[kk])
        
        log.lh.ini <- 0
        for(jj in 1:J)
        {
            log.lh.ini <- log.lh.ini + sum(n.e[, jj] * alpha[jj, kk] * sqrt.C0_gam0[[kk]] -  A * obs.ind * exp(Y.e[,jj]))
        }
        log.lh.ini <- log.lh.ini + theta
        
        # var = wd^2 / 12 -> wd = sqrt(12 * var)
        wd.phi0 = sqrt(12*propVar.phi0)
        
        if(abs(theta-l.phi0) < wd.phi0/2)
        {
            prop.lb <- l.phi0
            prop.ub <- l.phi0 + wd.phi0
        }else if(abs(theta-u.phi0) < wd.phi0/2)
        {
            prop.lb <- u.phi0 - wd.phi0
            prop.ub <- u.phi0
        }else
        {
            prop.lb <- theta - wd.phi0/2
            prop.ub <- theta + wd.phi0/2
        }
        
        theta.prop <- runif(1, prop.lb, prop.ub)
        phi0.prop <- exp(theta.prop)
        
        if(theta.prop > l.phi0 & theta.prop < u.phi0)
        {
            C0.e.base.prop <- corfunc(D.e.base, phi0.prop)
            sqrt.C0.eigen.prop <-  sqrt(Re(fft(C0.e.base.prop, TRUE)))
            
            if(sum(is.nan(sqrt.C0.eigen.prop)) != 0)
            {
                nan_phi0[kk] <- 1
                logR <- -Inf
            }else
            {
                sqrt.C0_gam0.prop <- as.vector(Re(fft(sqrt.C0.eigen.prop * norm.gam0[[kk]], TRUE)/sqrt(nC.e)))
                
                Y.e.prop <- Y.e
                for(jj in 1:J)
                {
                    Y.e.prop[,jj] <- Y.e.prop[,jj] - alpha[jj,kk] * sqrt.C0_gam0[[kk]]
                    Y.e.prop[,jj] <- Y.e.prop[,jj] + alpha[jj,kk] * sqrt.C0_gam0.prop
                }
                
                log.lh.prop <- 0
                for(jj in 1:J)
                {
                    log.lh.prop <- log.lh.prop + sum(n.e[, jj] * alpha[jj, kk] * sqrt.C0_gam0.prop -  A * obs.ind * exp(Y.e.prop[,jj]))
                }
                log.lh.prop <- log.lh.prop + theta.prop
                
                if(log.lh.prop == -Inf)
                {
                    logR <- -Inf
                }else
                {
                    log.prior.ini  <- 0
                    log.prior.prop  <- 0
                    
                    log.prop.prop  <- dunif(theta.prop, prop.lb, prop.ub, log = T)
                    log.prop.ini   <- dunif(theta, prop.lb, prop.ub, log = T)
                }
                logR  <- log.lh.prop - log.lh.ini + log.prior.prop - log.prior.ini + log.prop.ini - log.prop.prop
            }
            
            u = log(runif(1)) < logR
            
            if(u == 1)
            {
                phi0[kk] <- phi0.prop
                Y.e <- Y.e.prop
                
                C0.e.base[[kk]] <- C0.e.base.prop
                sqrt.C0.eigen[[kk]] <- sqrt.C0.eigen.prop
                sqrt.C0_gam0[[kk]] <-  sqrt.C0_gam0.prop
                
                for(jj in 1:J)
                {
                    aVec[,jj] <- n.e[,jj] -  A * obs.ind * exp(Y.e.prop[,jj])
                    norm.aVec[,jj] <- fft(matrix(aVec[,jj], M.e, N.e))/sqrt(nC.e)
                }
            }
            accept.phi0[kk] <- accept.phi0[kk] + u
        }
        
    }
    
    n.nan.PHI0 <- n.nan.PHI0 + nan_phi0
    
    return(list(phi0=phi0, Y.e=Y.e, C0.e.base=C0.e.base, sqrt.C0.eigen=sqrt.C0.eigen, sqrt.C0_gam0=sqrt.C0_gam0, aVec=aVec, norm.aVec=norm.aVec, accept.phi0=accept.phi0, n.nan.PHI0=n.nan.PHI0))
}




#####################
######     updating m -- adaptive M-H
#####################
Update_M <- function(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, mu.m, sigSq.m, accept.m, n.M)
{
    J <- length(m)
    if(is.null(phi0))
    {
        K <- 0
    }else
    {
        K <- length(phi0)
    }
    M.e <- dim(D.e.base)[1]
    N.e <- dim(D.e.base)[2]
    nC.e <- M.e * N.e
    
    for(jj in 1:J)
    {
        log.lh.ini <- sum(n.e[, jj] * m[jj] * obs.ind -  A * obs.ind * exp(Y.e[,jj]))
        
        D1 <- sum(aVec[obs.inx, jj])  - (m[jj]-mu.m)/sigSq.m
        D2 <- sum( - A * obs.ind *exp(Y.e[,jj])) - 1/sigSq.m
        
        m.prop.me <- m[jj] - D1/D2
        m.prop.var <- -2.4^2/D2
        
        m.prop <- rnorm(1, mean=m.prop.me, sd=sqrt(m.prop.var))
        
        Y.e.prop <- Y.e[,jj]
        Y.e.prop <- m.prop * obs.ind + sqrt(sigSq[jj])*sqrt.C_gam[[jj]]
        if(K >0)
        {
            for(kk in 1:K)
            {
                Y.e.prop <- Y.e.prop + alpha[jj, kk]*sqrt.C0_gam0[[kk]]
            }
        }
        Y.e.prop[Y.e.prop> 700] <- 700
        aVec.prop <- n.e[,jj] -  A * obs.ind * exp(Y.e.prop)
        
        log.lh.prop <- sum((n.e[, jj] * m.prop * obs.ind -  A * obs.ind * exp(Y.e.prop))[obs.inx])
        
        if(log.lh.prop == -Inf)
        {
            logR <- -Inf
        }else
        {
            D1.prop <- sum(aVec.prop[obs.inx])  - (m.prop-mu.m)/sigSq.m
            D2.prop <- sum((- A * obs.ind *exp(Y.e.prop))[obs.inx]) - 1/sigSq.m
            
            if(D1.prop == -Inf | D2.prop == -Inf)
            {
                logR <- -Inf
            }else
            {
                m.prop.me.ini <- m.prop - D1.prop/D2.prop
                m.prop.var.ini <- -2.4^2/D2.prop
                
                log.prior.ini  <- dnorm(m[jj], mean = mu.m, sd = sqrt(sigSq.m), log = TRUE)
                log.prior.prop  <- dnorm(m.prop, mean = mu.m, sd = sqrt(sigSq.m), log = TRUE)
                
                log.prop.prop  <- dnorm(m.prop, mean = m.prop.me, sd = sqrt(m.prop.var), log = TRUE)
                log.prop.ini   <- dnorm(m[jj], mean = m.prop.me.ini, sd = sqrt(m.prop.var.ini), log = TRUE)
                
                logR  <- log.lh.prop - log.lh.ini + log.prior.prop - log.prior.ini + log.prop.ini - log.prop.prop
            }
        }
        
        u = log(runif(1)) < logR
        
        if(u == 1)
        {
            m[jj] <- m.prop
            
            Y.e[,jj] <- Y.e.prop
            aVec[,jj] <- aVec.prop
            norm.aVec[,jj] <- fft(matrix(aVec[,jj], M.e, N.e))/sqrt(nC.e)
        }
        
        accept.m[jj] <- accept.m[jj] + u
    }
    
    return(list(m=m, Y.e=Y.e, aVec=aVec, norm.aVec=norm.aVec, accept.m=accept.m))
}




#########################
######     updating sigSq -- adaptive M-H
#########################
Update_SIG <- function(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, a.sig, b.sig, accept.sig, n.SIG)
{
    J <- dim(n.e)[2]
    if(is.null(phi0))
    {
        K <- 0
    }else
    {
        K <- length(phi0)
    }
    M.e <- dim(D.e.base)[1]
    N.e <- dim(D.e.base)[2]
    nC.e <- M.e * N.e
    
    for(jj in 1:J)
    {
        theta <- sqrt(sigSq[jj])
        log.lh.ini <- sum(n.e[, jj] * theta * sqrt.C_gam[[jj]] -  A * obs.ind * exp(Y.e[,jj])) + log(theta)
        
        D1 <- sum(aVec[,jj] * sqrt.C_gam[[jj]]) + (a.sig-1)/theta - b.sig + 1/theta
        D2 <- -sum(A * obs.ind * (sqrt.C_gam[[jj]])^2 * exp(Y.e[,jj])) - (a.sig-1)/theta^2 - 1/theta^2
        
        theta.prop.me <- theta - D1/D2
        theta.prop.var <- -2.4^2/D2
        
        theta.prop <- rtnorm(1, mean = theta.prop.me, sd = sqrt(theta.prop.var), lower = 0)
        sigSq.prop <- theta.prop^2
        
        Y.e.prop <- Y.e[,jj]
        Y.e.prop <- m[jj] * obs.ind + sqrt(sigSq.prop)*sqrt.C_gam[[jj]]
        if(K >0)
        {
            for(kk in 1:K)
            {
                Y.e.prop <- Y.e.prop + alpha[jj, kk]*sqrt.C0_gam0[[kk]]
            }
        }
        Y.e.prop[Y.e.prop> 700] <- 700
        aVec.prop <- n.e[,jj] -  A * obs.ind * exp(Y.e.prop)
        
        log.lh.prop <- sum((n.e[, jj] * theta.prop*sqrt.C_gam[[jj]] -  A * obs.ind * exp(Y.e.prop))[obs.inx]) + log(theta.prop)
        
        if(log.lh.prop == -Inf)
        {
            logR <- -Inf
        }else
        {
            D1.prop <- sum((aVec.prop * sqrt.C_gam[[jj]])[obs.inx]) + (a.sig-1)/theta.prop - b.sig + 1/theta.prop
            D2.prop <- -sum((A * obs.ind * (sqrt.C_gam[[jj]])^2 * exp(Y.e.prop))[obs.inx]) - (a.sig-1)/theta.prop^2 - 1/theta.prop^2
            
            if(D1.prop == -Inf | D2.prop == -Inf)
            {
                logR <- -Inf
            }else
            {
                theta.prop.me.ini <- theta.prop - D1.prop/D2.prop
                theta.prop.var.ini <- -2.4^2/D2.prop
                
                log.prior.ini  <- dgamma(theta, shape = a.sig, rate = b.sig, log = TRUE)
                log.prior.prop  <- dgamma(theta.prop, shape = a.sig, rate = b.sig, log = TRUE)
                
                log.prop.prop  <- dnorm(theta.prop, mean=theta.prop.me, sd=sqrt(theta.prop.var), log = TRUE)-pnorm(0, mean=theta.prop.me, sd=sqrt(theta.prop.var), log.p = TRUE, lower.tail = FALSE)
                
                log.prop.ini  <- dnorm(theta, mean=theta.prop.me.ini, sd=sqrt(theta.prop.var.ini), log = TRUE)-pnorm(0, mean=theta.prop.me.ini, sd=sqrt(theta.prop.var.ini), log.p = TRUE, lower.tail = FALSE)
                
                logR  <- log.lh.prop - log.lh.ini + log.prior.prop - log.prior.ini + log.prop.ini - log.prop.prop
            }
        }
        
        u = log(runif(1)) < logR
        
        if(u == 1)
        {
            sigSq[jj] <- sigSq.prop
            
            Y.e[,jj] <- Y.e.prop
            aVec[,jj] <- aVec.prop
            norm.aVec[,jj] <- fft(matrix(aVec[,jj], M.e, N.e))/sqrt(nC.e)
        }
        accept.sig[jj] <- accept.sig[jj] + u
    }
    
    return(list(sigSq = sigSq, Y.e = Y.e, aVec = aVec, norm.aVec = norm.aVec, accept.sig=accept.sig))
}





#########################
######     updating alpha -- species-specific weights for common GF
#########################
Update_ALP <- function(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, sigSq.alp, tauSq, accept.alp, n.ALP)
{
    J <- length(phi)
    if(is.null(phi0))
    {
        K <- 0
    }else
    {
        K <- length(phi0)
    }
    M.e <- dim(D.e.base)[1]
    N.e <- dim(D.e.base)[2]
    nC.e <- M.e * N.e
    
    p_add = 1/3
    p_del = 1/3
    
    for(jj in 1:J)
    {
        for(kk in 1:K)
        {
            log.lh.ini <- sum(n.e[, jj] * alpha[jj, kk] * sqrt.C0_gam0[[kk]] -  A * obs.ind * exp(Y.e[,jj]))
            
            D1 <- sum(aVec[,jj] * sqrt.C0_gam0[[kk]]) - (alpha[jj, kk])/sigSq.alp
            D2 <- -sum(A * obs.ind * (sqrt.C0_gam0[[kk]])^2 * exp(Y.e[,jj])) - 1/sigSq.alp
            
            alpha.prop.me <- alpha[jj,kk] - D1/D2
            alpha.prop.var <- -2.4^2/D2
            
            alpha.prop <- rnorm(1, mean=alpha.prop.me, sd=sqrt(alpha.prop.var))
            
            Y.e.prop <- Y.e[,jj]
            
            alpha.prop.vec <- alpha
            alpha.prop.vec[jj,kk] <- alpha.prop
            Y.e.prop <- m[jj] * obs.ind + sqrt(sigSq[jj])*sqrt.C_gam[[jj]]
            if(K >0)
            {
                for(ll in 1:K)
                {
                    Y.e.prop <- Y.e.prop + alpha.prop.vec[jj,ll]*sqrt.C0_gam0[[ll]]
                }
            }
            Y.e.prop[Y.e.prop> 700] <- 700
            aVec.prop <- n.e[,jj] -  A * obs.ind * exp(Y.e.prop)
            
            log.lh.prop <- sum((n.e[, jj] * alpha.prop * sqrt.C0_gam0[[kk]] -  A * obs.ind * exp(Y.e.prop))[obs.inx])
            
            if(log.lh.prop == -Inf)
            {
                logR <- -Inf
            }else
            {
                D1.prop <- sum((aVec.prop * sqrt.C0_gam0[[kk]])[obs.inx]) - (alpha.prop)/sigSq.alp
                D2.prop <- -sum((A * obs.ind * (sqrt.C0_gam0[[kk]])^2 * exp(Y.e.prop))[obs.inx]) - 1/sigSq.alp
                
                if(D1.prop == -Inf | D2.prop == -Inf)
                {
                    logR <- -Inf
                }else
                {
                    alpha.prop.me.ini <- alpha.prop - D1.prop/D2.prop
                    alpha.prop.var.ini <- -2.4^2/D2.prop
                    
                    log.prior.ini  <- dnorm(alpha[jj,kk], mean=0, sd=sqrt(sigSq.alp*tauSq[jj,kk]), log = TRUE)
                    log.prior.prop  <- dnorm(alpha.prop, mean=0, sd=sqrt(sigSq.alp*tauSq[jj,kk]), log = TRUE)
                    
                    log.prop.prop  <- dnorm(alpha.prop, mean = alpha.prop.me, sd = sqrt(alpha.prop.var), log = TRUE)
                    log.prop.ini   <- dnorm(alpha[jj,kk], mean = alpha.prop.me.ini, sd = sqrt(alpha.prop.var.ini), log = TRUE)
                    
                    logR  <- log.lh.prop - log.lh.ini + log.prior.prop - log.prior.ini + log.prop.ini - log.prop.prop
                }
            }
            
            u = log(runif(1)) < logR
            
            if(u == 1)
            {
                alpha[jj,kk] <- alpha.prop
                
                Y.e[,jj] <- Y.e.prop
                aVec[,jj] <- aVec.prop
                norm.aVec[,jj] <- fft(matrix(aVec[,jj], M.e, N.e))/sqrt(nC.e)
            }
            
            accept.alp[jj, kk] <- accept.alp[jj, kk] + u
        }
    }
    
    return(list(alpha=alpha, Y.e=Y.e, aVec=aVec, norm.aVec=norm.aVec, accept.alp=accept.alp))
}

Update_Tau <- function(lambdaSq, sigSq.alp, alpha, tauSq){
    
    J <- dim(alpha)[1]
    K <- dim(alpha)[2]
    
    nu.ind<-NULL
    nu=sqrt(lambdaSq * sigSq.alp/(alpha^2))
    nu.ind <- which(nu == Inf)
    if(length(nu.ind) > 0){nu[nu.ind] <- max(nu[-nu.ind]) + 10}
    
    gam <- matrix(NA, J, K)
    
    for(jj in 1:J)
    {
        for(kk in 1:K)
        {
            repeat
            {
                gam[jj, kk]  <- rinvGauss(1, nu = nu[jj, kk], lambda = lambdaSq)
                if (gam[jj, kk] > 0) break
            }
            tauSq[jj, kk] <- 1/gam[jj, kk]
        }
    }
    
    return(tauSq)
    
}


Update_Lambda<- function(tauSq, a.lam, b.lam)
{
    
    J <- dim(tauSq)[1]
    K <- dim(tauSq)[2]
    
    lambdaSq    <- rgamma(1, shape = (J+K+1)/2 + a.lam, rate = b.lam + sum(tauSq)/2)
    
    return(lambdaSq)
    
}




Llam <- function(n.e, A, obs.ind, D.e.base, gamma, m, sigSq, phi, gamma0, phi0, alpha)
{
    J <- length(phi)
    if(is.null(phi0))
    {
        K <- 0
    }else
    {
        K <- length(phi0)
    }
    
    M.e <- dim(D.e.base)[1]
    N.e <- dim(D.e.base)[2]
    nC.e <- M.e * N.e
    
    C.e.base <- vector("list", J)
    sqrt.C.eigen <- vector("list", J)
    norm.gam <- vector("list", J)
    sqrt.C_gam <- vector("list", J)
    
    for(jj in 1:J)
    {
        C.e.base[[jj]] <- corfunc(D.e.base, phi[jj])
        sqrt.C.eigen[[jj]] <-  sqrt(Re(fft(C.e.base[[jj]], TRUE)))
        norm.gam[[jj]] <- fft(matrix(gamma[,jj], M.e, N.e))/sqrt(nC.e)
        sqrt.C_gam[[jj]] <- as.vector(Re(fft(sqrt.C.eigen[[jj]] * norm.gam[[jj]], TRUE)/sqrt(nC.e)))
    }
    
    if(K >0)
    {
        C0.e.base <- vector("list", K)
        sqrt.C0.eigen <- vector("list", K)
        norm.gam0 <- vector("list", K)
        sqrt.C0_gam0 <- vector("list", K)
        
        for(kk in 1:K)
        {
            C0.e.base[[kk]] <- corfunc(D.e.base, phi0[kk])
            sqrt.C0.eigen[[kk]] <-  sqrt(Re(fft(C0.e.base[[kk]], TRUE)))
            norm.gam0[[kk]] <- fft(matrix(gamma0[,kk], M.e, N.e))/sqrt(nC.e)
            sqrt.C0_gam0[[kk]] <- as.vector(Re(fft(sqrt.C0.eigen[[kk]] * norm.gam0[[kk]], TRUE)/sqrt(nC.e)))
        }
    }
    
    Y.e <- matrix(NA, nC.e, J)
    aVec <- matrix(NA, nC.e, J)
    norm.aVec <- matrix(NA, nC.e, J)
    
    for(jj in 1:J)
    {
        Y.e[,jj] <- m[jj] * obs.ind + sqrt(sigSq[jj])*sqrt.C_gam[[jj]]
        if(K >0)
        {
            for(kk in 1:K)
            {
                Y.e[,jj] <- Y.e[,jj] + alpha[jj, kk]*sqrt.C0_gam0[[kk]]
            }
        }
        Y.e[Y.e[,jj]> 700, jj] <- 700
        aVec[,jj] <- n.e[,jj] -  A * obs.ind * exp(Y.e[,jj])
        norm.aVec[,jj] <- fft(matrix(aVec[,jj], M.e, N.e))/sqrt(nC.e)
    }
    
    if(K > 0)
    {
        return(list(C.e.base=C.e.base, sqrt.C.eigen=sqrt.C.eigen, norm.gam=norm.gam, sqrt.C_gam=sqrt.C_gam, C0.e.base=C0.e.base, sqrt.C0.eigen=sqrt.C0.eigen, norm.gam0=norm.gam0, sqrt.C0_gam0=sqrt.C0_gam0, Y.e=Y.e, aVec=aVec, norm.aVec=norm.aVec))
    }else
    {
        return(list(C.e.base=C.e.base, sqrt.C.eigen=sqrt.C.eigen, norm.gam=norm.gam, sqrt.C_gam=sqrt.C_gam, Y.e=Y.e, aVec=aVec, norm.aVec=norm.aVec))
    }
    
}



#########################
######     calculating log joint posterior probability
#########################
Update_logPost <- function(n.e, A, obs.ind, obs.inx, Y.e, gamma, m, sigSq, phi, gamma0, alpha, phi0, mu.m, sigSq.m, a.sig, b.sig, l.phi, u.phi, l.phi0, u.phi0, sigSq.alp, tauSq, lambdaSq, a.lam, b.lam)
{
    N <- dim(n.e)[1]
    J <- dim(n.e)[2]
    K <- length(phi0)
    
    logLH_ij <- matrix(NA, N, J)
    val <- 0
    pri.val2 <- rep(0, J)
    for(jj in 1:J)
    {
        logLH_ij[,jj] <- n.e[, jj] * (Y.e[,jj]+log(A)) -  A *  obs.ind * exp(Y.e[,jj]) - lfactorial(n.e[,jj])
        val <- val + sum(logLH_ij[,jj])
        pri.val2[jj] <- pri.val2[jj] + sum(dnorm(gamma[obs.inx,jj], mean=0, sd=1, log=T))
    }
    val <- val + sum(log(phi))
    val <- val + sum(log(sqrt(sigSq))) #final
    
    pri.val <- 0
    pri.val <- pri.val + sum(dunif(log(phi), min=l.phi, max=u.phi, log=T))
    pri.val <- pri.val + sum(dnorm(m, mean=mu.m, sd=sqrt(sigSq.m), log=T)) #final
    pri.val <- pri.val + sum(dgamma(sqrt(sigSq), shape = a.sig, rate = b.sig, log = TRUE)) #final
    
    pri.val3 <- 0
    if(K >0)
    {
        pri.val3 <- rep(0, K)
        for(kk in 1:K)
        {
            pri.val3[kk] <- pri.val3[kk] + sum(dnorm(gamma0[obs.inx,kk], mean=0, sd=1, log=T))
        }
        pri.val <- pri.val + sum(dunif(log(phi0), min=l.phi0, max=u.phi0, log=T))
        pri.val <- pri.val + sum(dnorm(alpha, mean=0, sd=sqrt(sigSq.alp*tauSq), log=T))
        
        pri.val <- pri.val + sum(dgamma(tauSq, shape = 1, rate = lambdaSq/2, log = TRUE)) #final
        pri.val <- pri.val + sum(dgamma(lambdaSq, shape = a.lam, rate = b.sig, log = TRUE)) #final
        
        val <- val + sum(log(phi0))
    }
    
    list(post=val+pri.val+sum(pri.val2)+sum(pri.val3), lh=val, lh.other=val+pri.val, lh.gam=val+sum(pri.val2), lh.gam0=val+sum(pri.val3), lh.gamgam0=val+sum(pri.val2)+sum(pri.val3), pri.gam=pri.val2, pri.gam0=pri.val3, logLH_ij=logLH_ij)
}



corfunc <- function(distance, param, model="exp")
{
    if(model == "exp")
    {
        phi <- param[1]
        val <- exp(-abs(distance)/phi)
    }
    return(val)
}


## alpha1 = alpha[j, ]: a vector of length K
## alpha2 = alpha[j', ]: a vector of length K
## param.corfunc: a vector of length K (e.g. phi0 under exponential corr)

ccf <- function(distance, alpha1, alpha2, param.corfunc, model = "exp")
{
    K <- length(alpha1)
    
    num <- 0
    for(k in 1:K)
    {
        num <- num + alpha1[k]*alpha2[k]*corfunc(distance, param.corfunc[k], model)
    }
    den <- sqrt(sum(alpha1^2) * sum(alpha2^2))
    
    if(num[1] == 0 & den[1] == 0)
    {
        val <- rep(0, length(distance))
    }else
    {
        val <- num/sqrt(sum(alpha1^2) * sum(alpha2^2))
    }
    
    return(val)
}



ccf2 <- function(distance, alpha1, alpha2, sigSq1, sigSq2, param.corfunc, model = "exp")
{
    K <- length(alpha1)
    
    num <- 0
    for(k in 1:K)
    {
        num <- num + alpha1[k]*alpha2[k]*corfunc(distance, param.corfunc[k], model)
    }
    den <- sqrt((sum(alpha1^2) + sigSq1) * (sum(alpha2^2) + sigSq2))
    
    if(num[1] == 0 & den[1] == 0)
    {
        val <- rep(0, length(distance))
    }else
    {
        val <- num/sqrt((sum(alpha1^2) + sigSq1) * (sum(alpha2^2) + sigSq2))
    }
    
    return(val)
}



PV <- function(distance, alpha_j, sigSq_j, param.corfunc_j, param.corfunc_0, model = "exp")
{
    K <- length(alpha_j)
    num <- 0
    for (k in 1:K) {
        num <- num + alpha_j[k]^2 * corfunc(distance, param.corfunc_0[k], model)
    }
    den <- num + sigSq_j * corfunc(distance, param.corfunc_j, model)
    
    val <- num/den
    
    return(val)
}




post_ccf <- function(fit, distance_vec, p.ind=NULL)
{
    J <- dim(fit$phi.p)[2]
    if(!is.null(fit$phi0.p))
    {
        K <- dim(fit$phi0.p)[2]
    }else
    {
        K <- 0
    }
    if(is.null(p.ind))
    {
        p.ind <- 1:dim(fit$phi.p)[1]
    }
    n.p <- length(p.ind)
    
    col1 <- NULL
    col2 <- NULL
    for(jj in 1:(J-1))
    {
        col1 <- c(col1, rep(jj, J-jj))
        col2 <- c(col2, (jj+1):J)
    }
    
    pair.ind <- cbind(col1, col2)
    n.pair <- dim(pair.ind)[1]
    n.dist <- length(distance_vec)
    
    ccf.p <- array(NA, c(n.pair, n.p, n.dist))
    
    for(rr in 1:n.dist)
    {
        for(ii in 1:n.pair)
        {
            type1.j <- pair.ind[ii, 1]
            type2.j <- pair.ind[ii, 2]
            
            post <- NULL
            for(MM in p.ind)
            {
                post <- c(post, ccf(distance_vec[rr], fit$alpha.p[type1.j,, MM], fit$alpha.p[type2.j,, MM], fit$phi0.p[MM,]))
            }
            ccf.p[ii, ,rr] <- post
        }
    }
    
    list(ccf.p=ccf.p, distance=distance_vec)
}


post_ccf2 <- function(fit, distance_vec, p.ind=NULL)
{
    J <- dim(fit$phi.p)[2]
    if(!is.null(fit$phi0.p))
    {
        K <- dim(fit$phi0.p)[2]
    }else
    {
        K <- 0
    }
    if(is.null(p.ind))
    {
        p.ind <- 1:dim(fit$phi.p)[1]
    }
    n.p <- length(p.ind)
    
    col1 <- NULL
    col2 <- NULL
    for(jj in 1:(J-1))
    {
        col1 <- c(col1, rep(jj, J-jj))
        col2 <- c(col2, (jj+1):J)
    }
    
    pair.ind <- cbind(col1, col2)
    n.pair <- dim(pair.ind)[1]
    n.dist <- length(distance_vec)
    
    ccf.p <- array(NA, c(n.pair, n.p, n.dist))
    
    for(rr in 1:n.dist)
    {
        for(ii in 1:n.pair)
        {
            type1.j <- pair.ind[ii, 1]
            type2.j <- pair.ind[ii, 2]
            
            post <- NULL
            for(MM in p.ind)
            {
                post <- c(post, ccf2(distance_vec[rr], fit$alpha.p[type1.j,, MM], fit$alpha.p[type2.j,, MM], fit$sigSq.p[MM,type1.j], fit$sigSq.p[MM,type2.j], fit$phi0.p[MM,]))
            }
            ccf.p[ii, ,rr] <- post
        }
    }
    
    list(ccf.p=ccf.p, distance=distance_vec)
}




post_PV <- function(fit, distance_vec, p.ind = NULL)
{
    J <- dim(fit$phi.p)[2]
    if (!is.null(fit$phi0.p)) {
        K <- dim(fit$phi0.p)[2]
    }
    else {
        K <- 0
    }
    if (is.null(p.ind)) {
        p.ind <- 1:dim(fit$phi.p)[1]
    }
    n.p <- length(p.ind)
    n.dist <- length(distance_vec)
    
    pv.p <- array(NA, c(J, n.p, n.dist))
    for (rr in 1:n.dist) {
        for (jj in 1:J) {
            post <- NULL
            for (MM in p.ind) {
                post <- c(post, PV(distance_vec[rr], fit$alpha.p[jj,,MM], fit$sigSq.p[MM,jj], fit$phi.p[MM,jj], fit$phi0.p[MM,]))
            }
            pv.p[jj, , rr] <- post
        }
    }
    
    list(pv.p = pv.p, distance = distance_vec)
}


# Bayesian multilevel modeling approach for meta-analysis
# z: transformed point estimates
# var: transformed uncertainties (e.g. posterior variances)
# cluster: cluster identifier for each estimate
# measure of interet: "ccf" or "pv"

BMM_meta <- function(z, var, cluster, measure="ccf", startValues=NULL, hyperParams=NULL, numReps, burnin, thin)
{
    N <- length(z)
    cluster.list <- unique(cluster)
    J <- length(cluster.list)
    nj <- as.vector(table(cluster))
    
    if(is.null(hyperParams))
    {
        mu0 <- 0
        sigSq0 <- 10^6
        a.b <- 0.001
        b.b <- 0.001
        a.w <- 0.001
        b.w <- 0.001
    }else
    {
        mu0 <- hyperParams$mu0
        sigSq0 <- hyperParams$sigSq0
        a.b <- hyperParams$a.b
        b.b <- hyperParams$b.b
        a.w <- hyperParams$a.w
        b.w <- hyperParams$b.w
    }
    
    if(is.null(startValues))
    {
        mu_j <- rep(NA, J)
        for(j in 1:J) mu_j[j] <- mean(z[which(cluster == j)])
        sigSq.b <- var(mu_j)
        mu <- mean(z)
        theta_ij <- z
        sigSq.w_j <- rep(NA, J)
        for(j in 1:J) sigSq.w_j[j] <- var(z[which(cluster == j)])
    }else
    {
        mu <- startValues$mu
        sigSq.b <- startValues$sigSq.b
        mu_j <- startValues$mu_j
        theta_ij <- startValues$theta_ij
        sigSq.w_j <- startValues$sigSq.w_j
    }
    
    nStore <- numReps/thin*(1-burnin)
    
    mu.p <- rep(NA, nStore)
    mu_j.p <- matrix(NA, nStore, J)
    
    for(MM in 1:numReps)
    {
        # update mu
        mu.mean <- (sum(mu_j)/sigSq.b+mu0/sigSq0)/(J/sigSq.b + 1/sigSq0)
        mu.var <- 1/(J/sigSq.b + 1/sigSq0)
        
        mu <- rnorm(1, mean=mu.mean, sd = sqrt(mu.var))
        
        # update sigSq.b
        shape.sigSq.b <- a.b + 0.5*J
        rate.sigSq.b <- b.b + 0.5*sum((mu_j - mu)^2)
        
        sigSq.b <- rigamma(1, a=shape.sigSq.b, b=rate.sigSq.b)
        
        # update mu_j
        for(j in 1:J)
        {
            mu_j.mean <- (sum(theta_ij[which(cluster == j)])/sigSq.w_j[j]+mu/sigSq.b)/(nj[j]/sigSq.w_j[j]+1/sigSq.b)
            mu_j.var <- 1/(nj[j]/sigSq.w_j[j] + 1/sigSq.b)
            
            mu_j[j] <- rnorm(1, mean=mu_j.mean, sd=sqrt(mu_j.var))
        }
        
        # update theta_ij
        for(i in 1:N)
        {
            j = which(cluster.list == cluster[i])
            theta_ij.mean <- (z[i]/var[i] + mu_j[j]/sigSq.w_j[j])/(1/var[i]+1/sigSq.w_j[j])
            theta_ij.var <- 1/(1/var[i]+1/sigSq.w_j[j])
            
            theta_ij[i] <- rnorm(1, mean=theta_ij.mean, sd=sqrt(theta_ij.var))
        }
        
        # update sigSq.w_j
        
        for(j in 1:J)
        {
            shape.sigSq.w_j <- a.w + 0.5*nj[j]
            rate.sigSq.w_j <- b.w + 0.5*sum((theta_ij[which(cluster == j)] - mu_j[j])^2)
            
            sigSq.w_j[j] <- rigamma(1, a=shape.sigSq.w_j, b=rate.sigSq.w_j)
        }
        
        if( (MM %% thin == 0) & (MM > (numReps * burnin)))
        {
            StoreInx <- MM/thin - numReps * burnin/thin
            
            mu.p[StoreInx] <- mu
            mu_j.p[StoreInx,] <- mu_j
        }
    }
    
    if(measure == "ccf")
    {
        mu.new <- (exp(2*mu.p)-1)/(exp(2*mu.p)+1)
        mu_j.new <- (exp(2*mu_j.p)-1)/(exp(2*mu_j.p)+1)
    }else if(measure == "pv")
    {
        mu.new <- exp(mu.p)/(1+exp(mu.p))
        mu_j.new <- exp(mu_j.p)/(1+exp(mu_j.p))
    }
    
    mu.med <- median(mu.new)
    mu.lb <- quantile(mu.new, prob = 0.025)
    mu.ub <- quantile(mu.new, prob = 0.975)
    mu.mean <- mean(mu.new)
    mu.sd <- sd(mu.new)
    
    mu_j.med <- apply(mu_j.new, 2, quantile, prob = 0.5)
    mu_j.lb <- apply(mu_j.new, 2, quantile, prob = 0.025)
    mu_j.ub <- apply(mu_j.new, 2, quantile, prob = 0.975)
    mu_j.mean <- apply(mu_j.new, 2, mean)
    mu_j.sd <- apply(mu_j.new, 2, sd)
    
    list(mu.med=mu.med, mu.lb=mu.lb, mu.ub=mu.ub, mu.mean=mu.mean, mu.sd=mu.sd, mu_j.med=mu_j.med, mu_j.lb=mu_j.lb, mu_j.ub=mu_j.ub, mu_j.mean=mu_j.mean, mu_j.sd=mu_j.sd)
}








