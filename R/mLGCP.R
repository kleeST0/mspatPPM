mLGCP <- function(data, K, hyperParams, startValues, mcmcParams, distance_vec, model="exp")
{
    #########
    ### Data
    dat <- list()
    J <- dat$J <- data$J
    ext <- dat$ext <- data$ext
    chain <- mcmcParams$run$chain
    dataName <- mcmcParams$run$dataName
    
    D.e.base <- dat$D.e.base <- data$D.ext.row1
    n.e <- dat$n.e <- data$n.ji.ext
    obs.ind <- dat$obs.ind <- data$obs.ind
    obs.inx <- dat$obs.inx <- data$obs.inx
    A <- dat$A <- data$A
    dat$scaleF <- data$scaleF
    
    M.e <- dat$M.e <- dim(D.e.base)[1]
    N.e <- dat$N.e <- dim(D.e.base)[2]
    nC.e <- dat$nC.e <- M.e * N.e
    
    #########
    ### Prior
    mu.m <- hyperParams$mu.m
    sigSq.m <- hyperParams$sigSq.m
    a.sig <- hyperParams$a.sig
    b.sig <- hyperParams$b.sig
    l.phi <- hyperParams$l.phi
    u.phi <- hyperParams$u.phi
    
    sigSq.alp <- hyperParams$sigSq.alp
    l.phi0 <- hyperParams$l.phi0
    u.phi0 <- hyperParams$u.phi0
    
    a.lam <- hyperParams$a.lam
    b.lam <- hyperParams$b.lam
        
    #########
    ### Move probabilites
    p.GAM <- mcmcParams$pr.move$p.GAM
    p.GAM0 <- mcmcParams$pr.move$p.GAM0
    p.PHI <- mcmcParams$pr.move$p.PHI
    p.PHI0 <- mcmcParams$pr.move$p.PHI0
    p.M <- mcmcParams$pr.move$p.M
    p.SIG <- mcmcParams$pr.move$p.SIG
    p.ALP <- mcmcParams$pr.move$p.ALP
    
    #########
    ### MCMC parameter
    
    propVar.phi <- mcmcParams$propVar.phi
    propVar.phi0 <- mcmcParams$propVar.phi0
    
    L.gam <- mcmcParams$L.gam
    eps.gam <- mcmcParams$eps.gam
    M.gam.var <- mcmcParams$M.gam.var
    
    L.gam0 <- mcmcParams$L.gam0
    eps.gam0 <- mcmcParams$eps.gam0
    M.gam0.var <- mcmcParams$M.gam0.var
    
    numReps     <- mcmcParams$run$numReps
    thin        <- mcmcParams$run$thin
    burninPerc  <- mcmcParams$run$burninPerc
    nSave  		<- mcmcParams$run$nSave
    nStore      <- numReps/thin * (1 - burninPerc)
    
    #########
    ### Starting values
    gamma <- startValues$gamma
    m <- startValues$m
    phi <- startValues$phi
    sigSq <- startValues$sigSq
    gamma0 <- startValues$gamma0
    alpha <- startValues$alpha
    phi0 <- startValues$phi0
    
    if(is.null(phi0))
    {
        K <- 0
        tauSq <- NULL
        lambdaSq <- NULL
    }else
    {
        K <- length(phi0)
        tauSq <- matrix(0.1, J, K)
        lambdaSq <- 0.01
    }
    
    #########
    ### Posterior samples
    
    saveGamInx <- c(sample(obs.inx, min(1, length(obs.inx), floor(nC.e/2))), sample(setdiff(1:nC.e, obs.inx), min(1, floor(nC.e/2))))
    
    gamma.p <- array(NA, c(length(saveGamInx), J, nStore))
    m.p <- matrix(NA, nStore, J)
    sigSq.p <- matrix(NA, nStore, J)
    phi.p <- matrix(NA, nStore, J)
    gamma0.p <- array(NA, c(length(saveGamInx), K, nStore))
    alpha.p <- array(NA, c(J, K, nStore))
    phi0.p <- matrix(NA, nStore, K)
    lambdaSq.p <- rep(NA, nStore)
    
    accept.gam <- rep(0, J)
    accept.gam.100 <- vector("list", J)
    accept.m <- rep(0, J)
    accept.sig <- rep(0, J)
    accept.phi <- rep(0, J)
    accept.gam0 <- rep(0, K)
    accept.gam0.100 <- vector("list", K)
    accept.alp <- matrix(0, J, K)
    accept.phi0 <- rep(0, K)
    
    logpost.p <- rep(NA, nStore)
    loglh.p <- rep(NA, nStore)
    loglh.other.p <- rep(NA, nStore)
    loglh.gam.p <- rep(NA, nStore)
    loglh.gam0.p <- rep(NA, nStore)
    loglh.gamgam0.p <- rep(NA, nStore)
    logpri.gam.p <- matrix(NA, nStore, J)
    logpri.gam0.p <- matrix(NA, nStore, K)
    
    n.GAM <- 0
    n.M <- 0
    n.SIG <- 0
    n.PHI <- 0
    n.GAM0 <- 0
    n.ALP <- 0
    n.PHI0 <- 0
    
    n.nan.GAM <- rep(0, J)
    n.nan.GAM0 <- rep(0, K)
    n.nan.PHI <- rep(0, J)
    n.nan.PHI0 <- rep(0, K)
    
    Lam.mean.all <- matrix(0, nC.e, J)
    Lam.var.all <- matrix(0, nC.e, J)
    logLam.mean.all <- matrix(0, nC.e, J)
    logLam.var.all <- matrix(0, nC.e, J)
    
    Lam.mean <- matrix(0, nC.e, J)
    Lam.var <- matrix(0, nC.e, J)
    logLam.mean <- matrix(0, nC.e, J)
    logLam.var <- matrix(0, nC.e, J)
    
    gamma.mean <- matrix(0, nC.e, J)
    gamma.var <- matrix(0, nC.e, J)
    
    gamma0.mean <- matrix(0, nC.e, K)
    gamma0.var <- matrix(0, nC.e, K)
    
    LH_ij.mean <- matrix(0, nC.e, J)
    
    #########
    ### Initialize the intensity functions
    
    res <- Llam(n.e, A, obs.ind, D.e.base, gamma, m, sigSq, phi, gamma0, phi0, alpha)
    
    C.e.base <- res$C.e.base
    sqrt.C.eigen <- res$sqrt.C.eigen
    norm.gam <- res$norm.gam
    sqrt.C_gam <- res$sqrt.C_gam
    C0.e.base <- res$C0.e.base
    sqrt.C0.eigen <- res$sqrt.C0.eigen
    norm.gam0 <- res$norm.gam0
    sqrt.C0_gam0 <- res$sqrt.C0_gam0
    Y.e <- res$Y.e
    aVec <- res$aVec
    norm.aVec <- res$norm.aVec
    
    logpost <- Update_logPost(n.e, A, obs.ind, obs.inx, Y.e, gamma, m, sigSq, phi, gamma0, alpha, phi0, mu.m, sigSq.m, a.sig, b.sig, l.phi, u.phi, l.phi0, u.phi0,  sigSq.alp, tauSq, lambdaSq, a.lam, b.lam)$post
    
    for(MM in 1:numReps)
    {
        if(MM %% 10000 == 0)
        {
            cat("iteration: ", MM, " out of ", numReps, fill = T)
        }
        
        move <- sample(c("GAM", "M", "SIG", "PHI", "GAM0", "PHI0", "ALP"), 1, prob=c(p.GAM, p.M, p.SIG, p.PHI, p.GAM0, p.PHI0, p.ALP))
        
        # gamma
        if(move == "GAM")
        {
            n.GAM <- n.GAM + 1
            res <- Update_GAM(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C.eigen, norm.gam, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, L.gam, eps.gam, M.gam.var, accept.gam, accept.gam.100, n.GAM, n.nan.GAM, MM, numReps)
            
            gamma <- res$gamma
            norm.gam <- res$norm.gam
            sqrt.C_gam <- res$sqrt.C_gam
            
            Y.e <- res$Y.e
            aVec <- res$aVec
            norm.aVec <- res$norm.aVec
            
            n.nan.GAM <- res$n.nan.GAM
            accept.gam <- res$accept.gam
            accept.gam.100 <- res$accept.gam.100
            eps.gam <- res$eps.gam
        }
        
        # m
        if(move == "M")
        {
            n.M <- n.M + 1
            res <- Update_M(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, mu.m, sigSq.m, accept.m, n.M)
            
            m <- res$m
            Y.e <- res$Y.e
            aVec <- res$aVec
            norm.aVec <- res$norm.aVec
            accept.m <- res$accept.m
        }
        
        # sigma
        if(move == "SIG")
        {
            n.SIG <- n.SIG + 1
            res <- Update_SIG(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, a.sig, b.sig, accept.sig, n.SIG)
            
            sigSq<- res$sigSq
            Y.e <- res$Y.e
            aVec <- res$aVec
            norm.aVec <- res$norm.aVec
            accept.sig <- res$accept.sig
        }
        
        # phi
        if(move == "PHI")
        {
            n.PHI <- n.PHI + 1
            res <- Update_PHI(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, C.e.base, sqrt.C.eigen, norm.gam, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, l.phi, u.phi, propVar.phi, accept.phi, n.PHI, n.nan.PHI)
            
            phi <- res$phi
            C.e.base <- res$C.e.base
            sqrt.C.eigen <- res$sqrt.C.eigen
            sqrt.C_gam <- res$sqrt.C_gam
            Y.e <- res$Y.e
            aVec <- res$aVec
            norm.aVec <- res$norm.aVec
            n.nan.PHI <- res$n.nan.PHI
            accept.phi <- res$accept.phi
        }
        
        # gamma0
        if(move == "GAM0")
        {
            n.GAM0 <- n.GAM0 + 1
            res <- Update_GAM0(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C_gam, C0.e.base, sqrt.C0.eigen, norm.gam0, sqrt.C0_gam0, Y.e, aVec, norm.aVec, L.gam0, eps.gam0, M.gam0.var, accept.gam0, accept.gam0.100, n.GAM0, n.nan.GAM0, MM, numReps)
            
            gamma0 <- res$gamma0
            norm.gam0 <- res$norm.gam0
            sqrt.C0_gam0 <- res$sqrt.C0_gam0
            Y.e <- res$Y.e
            aVec <- res$aVec
            norm.aVec <- res$norm.aVec
            
            n.nan.GAM0 <- res$n.nan.GAM0
            accept.gam0 <- res$accept.gam0
            accept.gam0.100 <- res$accept.gam0.100
            eps.gam0 <- res$eps.gam0
        }
        
        # alpha
        if(move == "ALP")
        {
            n.ALP <- n.ALP + 1
            res <- Update_ALP(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, sqrt.C_gam, sqrt.C0_gam0, Y.e, aVec, norm.aVec, sigSq.alp, tauSq, accept.alp, n.ALP)
            
            alpha   <- res$alpha
            Y.e     <- res$Y.e
            aVec    <- res$aVec
            norm.aVec <- res$norm.aVec
            accept.alp <- res$accept.alp
            
            tauSq <- Update_Tau(lambdaSq, sigSq.alp, alpha, tauSq)
            lambdaSq <- Update_Lambda(tauSq, a.lam, b.lam)
        }
        
        
        # phi0
        if(move == "PHI0")
        {
            n.PHI0 <- n.PHI0 + 1
            res <- Update_PHI0(n.e, A, obs.ind, obs.inx, D.e.base, gamma, m, sigSq, phi, gamma0, alpha, phi0, C0.e.base, sqrt.C0.eigen, norm.gam0, sqrt.C0_gam0, Y.e, aVec, norm.aVec, l.phi0, u.phi0, propVar.phi0, accept.phi0, n.PHI0, n.nan.PHI0)
            
            phi0 <- res$phi0
            C0.e.base <- res$C0.e.base
            sqrt.C0.eigen <- res$sqrt.C0.eigen
            sqrt.C0_gam0 <- res$sqrt.C0_gam0
            Y.e <- res$Y.e
            aVec <- res$aVec
            norm.aVec <- res$norm.aVec
            n.nan.PHI0 <- res$n.nan.PHI0
            accept.phi0 <- res$accept.phi0
        }
  
        eval <- Update_logPost(n.e, A, obs.ind, obs.inx, Y.e, gamma, m, sigSq, phi, gamma0, alpha, phi0, mu.m, sigSq.m, a.sig, b.sig, l.phi, u.phi, l.phi0, u.phi0, sigSq.alp, tauSq, lambdaSq, a.lam, b.lam)
        
        logpost <- eval$post
        loglh <- eval$lh
        loglh.other <- eval$lh.other
        loglh.gam <- eval$lh.gam
        loglh.gam0 <- eval$lh.gam0
        loglh.gamgam0 <- eval$lh.gamgam0
        logpri.gam <- eval$pri.gam
        logpri.gam0 <- eval$pri.gam0
        LH_ij <- exp(eval$logLH_ij)
        
        ###################################
        ######     storing posterior samples
        ###################################
        
        if( (MM %% thin == 0) & (MM > (numReps * burninPerc)))
        {
            StoreInx <- MM/thin - numReps * burninPerc/thin
            
            gamma.p[,,StoreInx] <- gamma[saveGamInx,]
            m.p[StoreInx, ] <- m
            phi.p[StoreInx,] <- phi
            sigSq.p[StoreInx,] <- sigSq
            
            if(K > 0)
            {
                gamma0.p[,,StoreInx] <- gamma0[saveGamInx,]
                phi0.p[StoreInx,] <- phi0
                alpha.p[,,StoreInx] <- alpha
                lambdaSq.p[StoreInx] <- lambdaSq
            }
            
            if(StoreInx > 1)
            {
                Lam.var <- (StoreInx-1)/StoreInx * Lam.var + (StoreInx-1)/StoreInx^2 * (Lam.mean - exp(Y.e))^2
                logLam.var <- (StoreInx-1)/StoreInx * logLam.var + (StoreInx-1)/StoreInx^2 * (logLam.mean - Y.e)^2
                
                gamma.var <- (StoreInx-1)/StoreInx * gamma.var + (StoreInx-1)/StoreInx^2 * (gamma.mean - gamma)^2
                gamma0.var <- (StoreInx-1)/StoreInx * gamma0.var + (StoreInx-1)/StoreInx^2 * (gamma0.mean - gamma0)^2
            }
            
            Lam.mean <- ((StoreInx-1)* Lam.mean + exp(Y.e))/StoreInx
            logLam.mean <- ((StoreInx-1)* logLam.mean + Y.e)/StoreInx
            
            gamma.mean <- ((StoreInx-1)* gamma.mean + gamma)/StoreInx
            gamma0.mean <- ((StoreInx-1)* gamma0.mean + gamma0)/StoreInx
            
            LH_ij.mean <- ((StoreInx-1)* LH_ij.mean + LH_ij)/StoreInx
            
            logpost.p[StoreInx] <- logpost
            loglh.p[StoreInx] <- loglh
            loglh.other.p[StoreInx] <- loglh.other
            loglh.gam.p[StoreInx] <- loglh.gam
            loglh.gam0.p[StoreInx] <- loglh.gam0
            loglh.gamgam0.p[StoreInx] <- loglh.gamgam0
            logpri.gam.p[StoreInx,] <- logpri.gam
            logpri.gam0.p[StoreInx,] <- logpri.gam0
        }
        
        if((MM %% nSave == 0) & (MM > (numReps * burninPerc)))
        {
            last <- list()
            last$gamma <- gamma
            last$m <- m
            last$sigSq <- sigSq
            last$phi <- phi
            
            if(K > 0)
            {
                last$gamma0 <- gamma0
                last$alpha <- alpha
                last$phi0 <- phi0
            }
            
            m.pm <- apply(m.p, 2, mean, na.rm = T)
            sigSq.pm <- apply(sigSq.p, 2, mean, na.rm = T)
            phi.pm <- apply(phi.p, 2, mean, na.rm = T)
            
            if(K > 0)
            {
                alpha.pm <- apply(alpha.p, c(1,2), mean, na.rm = T)
                phi0.pm <- apply(phi0.p, 2, mean, na.rm = T)
            }else
            {
                alpha.pm <- NULL
                phi0.pm <- NULL
            }
            
            #########################

            loglh.pm <- Update_logPost(n.e, A, obs.ind, obs.inx, Y.e, gamma.mean, m.pm, sigSq.pm, phi.pm, gamma0.mean, alpha.pm, phi0.pm, mu.m, sigSq.m, a.sig, b.sig, l.phi, u.phi, l.phi0, u.phi0, sigSq.alp, tauSq, lambdaSq, a.lam, b.lam)$lh
            
            DIC <- -4 * mean(loglh.p, na.rm=T) + 2 * loglh.pm
            
            DIC2 <- -4 * mean(loglh.p, na.rm=T) + 2 * sum((n.e*(logLam.mean+log(A)) - A*Lam.mean - lfactorial(n.e))[obs.inx,])
            
            DIC3 <- -4 * mean(loglh.p, na.rm=T) + 2 * sum(log(LH_ij.mean))
            
            
            fit <- list(Lam.mean.all=Lam.mean.all, Lam.var.all=Lam.var.all, logLam.mean.all=logLam.mean.all, logLam.var.all=logLam.var.all, Lam.mean=Lam.mean, Lam.var=Lam.var, logLam.mean=logLam.mean, logLam.var=logLam.var, m.p=m.p, phi0.p=phi0.p, phi.p=phi.p, alpha.p = alpha.p, sigSq.p=sigSq.p, lambdaSq.p=lambdaSq.p, accept.alp=accept.alp/n.ALP, accept.gam=accept.gam/n.GAM, accept.gam0=accept.gam0/n.GAM0, accept.phi=accept.phi/n.PHI,accept.phi0=accept.phi0/n.PHI0, accept.m=accept.m/n.M, accept.sig=accept.sig/n.SIG, startValues=startValues, hyperParams=hyperParams, mcmcParams=mcmcParams, logpost.p=logpost.p, loglh.p=loglh.p, n.nan.PHI=n.nan.PHI, n.nan.PHI0=n.nan.PHI0, n.nan.GAM=n.nan.GAM, n.nan.GAM0=n.nan.GAM0, eps.gam=eps.gam, eps.gam0=eps.gam0, gamma.mean=gamma.mean, gamma0.mean=gamma0.mean, gamma.var=gamma.var, gamma0.var=gamma0.var, saveGamInx=saveGamInx, last=last, dat=dat, DIC=DIC3)
            
            if(MM == numReps)
            {
                results <- list()
                results$DIC <- DIC
                results$DIC2 <- DIC2
                results$DIC3 <- DIC3
                
                if(K > 0)
                {
                    ## cross-correlation functions
                    res.ccf <- post_ccf(fit, distance_vec)
                    ccf.p <- res.ccf$ccf.p
                    suppressWarnings(z.ccf.p <- 0.5 * log((1+ccf.p)/(1-ccf.p)))
                    
                    med.ccf <- apply(ccf.p, c(1,3), quantile, prob = 0.5)
                    ub.ccf <- apply(ccf.p, c(1,3), quantile, prob = 0.975)
                    lb.ccf <- apply(ccf.p, c(1,3), quantile, prob = 0.025)
                    
                    mean.ccf <- apply(ccf.p, c(1,3), mean)
                    sd.ccf <- apply(ccf.p, c(1,3), sd)
                    
                    mean.z.ccf <- apply(z.ccf.p, c(1,3), mean)
                    var.z.ccf <- apply(z.ccf.p, c(1,3), var)
                    
                    ## cross-correlation functions (ccf2: based on all intensity components)
                    res.ccf2 <- post_ccf2(fit, distance_vec)
                    ccf2.p <- res.ccf2$ccf.p
                    z.ccf2.p <- 0.5 * log((1+ccf2.p)/(1-ccf2.p))
                    
                    med.ccf2 <- apply(ccf2.p, c(1,3), quantile, prob = 0.5)
                    ub.ccf2 <- apply(ccf2.p, c(1,3), quantile, prob = 0.975)
                    lb.ccf2 <- apply(ccf2.p, c(1,3), quantile, prob = 0.025)
                    
                    mean.ccf2 <- apply(ccf2.p, c(1,3), mean)
                    sd.ccf2 <- apply(ccf2.p, c(1,3), sd)
                    
                    mean.z.ccf2 <- apply(z.ccf2.p, c(1,3), mean)
                    var.z.ccf2 <- apply(z.ccf2.p, c(1,3), var)
                    
                    ## proportion of varince functions
                    res.pv <- post_PV(fit, distance_vec)
                    pv.p <- pv.p.new <- res.pv$pv.p
                    pv.p.new[pv.p.new == 1] <- 0.999999
                    suppressWarnings(z.pv.p <- log(pv.p.new/(1-pv.p.new)))
                    
                    med.pv <- apply(pv.p, c(1,3), quantile, prob = 0.5)
                    ub.pv <- apply(pv.p, c(1,3), quantile, prob = 0.975)
                    lb.pv <- apply(pv.p, c(1,3), quantile, prob = 0.025)
                    
                    mean.pv <- apply(pv.p, c(1,3), mean)
                    sd.pv <- apply(pv.p, c(1,3), sd)
                    
                    mean.z.pv <- apply(z.pv.p, c(1,3), mean)
                    var.z.pv <- apply(z.pv.p, c(1,3), var)
                    
                    fit$distance_vec <- distance_vec
                    fit$ccf.p <- ccf.p
                    fit$ccf2.p <- ccf2.p
                    fit$pv.p <- pv.p
                    
                    results$ccf.med <- med.ccf
                    results$ccf.ub <- ub.ccf
                    results$ccf.lb <- lb.ccf
                    results$ccf.mean <- mean.ccf
                    results$ccf.sd <- sd.ccf
                    
                    results$ccf2.med <- med.ccf2
                    results$ccf2.ub <- ub.ccf2
                    results$ccf2.lb <- lb.ccf2
                    results$ccf2.mean <- mean.ccf2
                    results$ccf2.sd <- sd.ccf2
                    
                    results$pv.med <- med.pv
                    results$pv.ub <- ub.pv
                    results$pv.lb <- lb.pv
                    results$pv.mean <- mean.pv
                    results$pv.sd <- sd.pv
                    
                    results$mean.z.ccf <- mean.z.ccf
                    results$var.z.ccf <- var.z.ccf
                    results$mean.z.ccf2 <- mean.z.ccf2
                    results$var.z.ccf2 <- var.z.ccf2
                    
                    results$mean.z.pv <- mean.z.pv
                    results$var.z.pv <- var.z.pv
                }
            }
        }
    }    ### end of MCMC
    
    if(K >0)
    {
        list(Lam.mean.all=Lam.mean.all, Lam.var.all=Lam.var.all, logLam.mean.all=logLam.mean.all, logLam.var.all=logLam.var.all, Lam.mean=Lam.mean, Lam.var=Lam.var, logLam.mean=logLam.mean, logLam.var=logLam.var, m.p=m.p, phi0.p=phi0.p, phi.p=phi.p, alpha.p = alpha.p, sigSq.p=sigSq.p, lambdaSq.p=lambdaSq.p, accept.alp=accept.alp/n.ALP, accept.gam=accept.gam/n.GAM, accept.gam0=accept.gam0/n.GAM0, accept.phi=accept.phi/n.PHI,accept.phi0=accept.phi0/n.PHI0, accept.m=accept.m/n.M, accept.sig=accept.sig/n.SIG, startValues=startValues, hyperParams=hyperParams, mcmcParams=mcmcParams, logpost.p=logpost.p, loglh.p=loglh.p, n.nan.PHI=n.nan.PHI, n.nan.PHI0=n.nan.PHI0, n.nan.GAM=n.nan.GAM, n.nan.GAM0=n.nan.GAM0, eps.gam=eps.gam, eps.gam0=eps.gam0, gamma.mean=gamma.mean, gamma0.mean=gamma0.mean, gamma.var=gamma.var, gamma0.var=gamma0.var, saveGamInx=saveGamInx, last=last, dat=dat, DIC=DIC3, distance_vec=distance_vec, ccf.p=ccf.p,ccf2.p=ccf2.p, pv.p=pv.p, results=results)
    }else
    {
        list(Lam.mean.all=Lam.mean.all, Lam.var.all=Lam.var.all, logLam.mean.all=logLam.mean.all, logLam.var.all=logLam.var.all, Lam.mean=Lam.mean, Lam.var=Lam.var, logLam.mean=logLam.mean, logLam.var=logLam.var, m.p=m.p, phi0.p=phi0.p, phi.p=phi.p, alpha.p = alpha.p, sigSq.p=sigSq.p, lambdaSq.p=lambdaSq.p, accept.alp=accept.alp/n.ALP, accept.gam=accept.gam/n.GAM, accept.gam0=accept.gam0/n.GAM0, accept.phi=accept.phi/n.PHI,accept.phi0=accept.phi0/n.PHI0, accept.m=accept.m/n.M, accept.sig=accept.sig/n.SIG, startValues=startValues, hyperParams=hyperParams, mcmcParams=mcmcParams, logpost.p=logpost.p, loglh.p=loglh.p, n.nan.PHI=n.nan.PHI, n.nan.PHI0=n.nan.PHI0, n.nan.GAM=n.nan.GAM, n.nan.GAM0=n.nan.GAM0, eps.gam=eps.gam, eps.gam0=eps.gam0, gamma.mean=gamma.mean, gamma0.mean=gamma0.mean, gamma.var=gamma.var, gamma0.var=gamma0.var, saveGamInx=saveGamInx, last=last, dat=dat, DIC=DIC3, distance_vec=distance_vec, results=results)
    }

}
