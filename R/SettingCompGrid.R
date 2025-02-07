

## ppp: ppp object
## MN: A two-element vector (M, N) specifying the dimensions of the computational grid, where M is the number of cells along the horizontal axis and N is the number of cells along the vertical axis.
## ext: A value to determine extended computational grid for circulant embedding
##      c.f.) ext*(M x N)
## object.only: If TRUE, the analysis will focus on the object specified within its esimtated boundary; otherwise, the entire window will be considered for the intensity function

## exclude: An optional argument when object.only=TRUE. If TRUE, the last type (mark) in the "ppp" object will be considered as coordinates to be excluded


SettingCompGrid <- function(ppp, MN, ext, object.only = F, exclude=F)
{
    if(object.only == F & exclude==T)
    {
        cat("\n Analysis for the entire image window")
        cat("\n 'exclude' argument is ignored")
    }else if(object.only == F & exclude==F)
    {
        cat("\n Analysis for the entire image window")
        cat("\n 'exclude' argument is ignored")
    }else if(object.only == T & exclude==T)
    {
        cat("\n Analysis focusing on the object")
        cat("\n exclude=TRUE: the last type in the 'ppp' object will be used for exclusion")
    }else if(object.only == T & exclude==F)
    {
        cat("\n Analysis focusing on the object")
        cat("\n exclude=FALSE: only the estimated boundary of the object will be considered")
    }
    
    type <- ppp$type <- names((table(ppp$marks)))
    J <- ppp$J <- length(type)
    
    marksNum <- as.vector(ppp$marks, mode = "character")
    for(j in 1:J)
    {
        marksNum[marksNum == type[j]] <- j
    }
    marksNum <- ppp$marksNum <- as.vector(marksNum, mode = "numeric")
    
    if(ppp$window$xrange[2] >= ppp$window$yrange[2])
    {
        if(ppp$window$xrange[2]/2 >= ppp$window$yrange[2])
        {
            M <- ppp$M <- MN[1]
            N <- ppp$N <- M/2
            W <- ppp$W <- owin(c(0, ppp$window$xrange[2]), c(0, ppp$window$xrange[2]/2))
        }else
        {
            M <- ppp$M <- MN[1]
            N <- ppp$N <- MN[2]
            W <- ppp$W <- owin(c(0, ppp$window$xrange[2]), c(0, ppp$window$xrange[2]))
        }
    }else
    {
        if(ppp$window$yrange[2]/2 >= ppp$window$xrange[2])
        {
            N <- ppp$N <- MN[1]
            M <- ppp$M <- N/2
            W <- ppp$W <- owin(c(0, ppp$window$yrange[2]/2), c(0, ppp$window$yrange[2]))
        }else
        {
            N <- ppp$N <- MN[1]
            M <- ppp$M <- MN[2]
            W <- ppp$W <- owin(c(0, ppp$window$yrange[2]), c(0, ppp$window$yrange[2]))
        }
    }
    nC <- ppp$nC <- M * N
    
    mygrid <- ppp$mygrid <- grid.prep(W, M, N, ext)
    
    M.ext <- ppp$M.ext <- mygrid$M.ext
    N.ext <- ppp$N.ext <- mygrid$N.ext
    nC.ext <- ppp$nC.ext <- M.ext * N.ext
    cell.width <- ppp$cell.width <- mygrid$cell.width
    cell.height <- ppp$cell.height <- mygrid$cell.height
    mcens <- ppp$mcens <- mygrid$mcens
    ncens <- ppp$ncens <- mygrid$ncens
    mcens.ext <- ppp$mcens.ext <- mygrid$mcens.ext
    ncens.ext <- ppp$ncens.ext <- mygrid$ncens.ext
    
    W.im <- matrix(NA, N, M)
    for (i in 1:ncol(W.im)) {
        W.im[, i] <- as.numeric(inside.owin(x = rep(mygrid$mcens[i], nrow(W.im)), y = mygrid$ncens,
        w = W))
    }
    W.im[W.im == 0] <- NA
    W.im <- ppp$W.im <-  im(W.im, xcol = mygrid$mcens, yrow = mygrid$ncens)
    
    W.im.ext <- matrix(NA, N.ext, M.ext)
    for (i in 1:ncol(W.im.ext)) {
        W.im.ext[, i] <- as.numeric(inside.owin(x = rep(mygrid$mcens.ext[i], nrow(W.im.ext)),
        y = mygrid$ncens.ext, w = W))
    }
    W.im.ext[W.im.ext == 0] <- NA
    W.im.ext <- ppp$W.im.ext <- im(W.im.ext, xcol = mygrid$mcens.ext, yrow = mygrid$ncens.ext)
    
    ulx.ext <- ppp$ulx.ext <- rep(mygrid$mgrid[-1], N)
    uly.ext <- ppp$uly.ext <- rep(mygrid$ngrid[-1], each=M)
    
    llx.ext <- ppp$llx.ext <- rep(mygrid$mgrid[-(M+1)], N)
    lly.ext <- ppp$lly.ext <- rep(mygrid$ngrid[-(N+1)], each=M)
    
    CalNji     <- .C("Cal_nji_loop",
    nn = as.integer(ppp$n),
    nc = as.integer(nC),
    Jval = as.integer(J),
    typeVec = as.double(1:J),
    marksVec = as.double(marksNum),
    xVec = as.double(ppp$x),
    yVec = as.double(ppp$y),
    ulxExt = as.double(ulx.ext),
    llxExt = as.double(llx.ext),
    ulyExt = as.double(uly.ext),
    llyExt = as.double(lly.ext),
    nji = as.double(rep(0, nC*J)))
    
    n.ji <- ppp$n.ji <- matrix(as.vector(CalNji$nji[1:(nC*J)]), nrow = nC, byrow = FALSE)
    rm(CalNji)
    
    Rx <- ppp$Rx <- M.ext * cell.width
    Ry <- ppp$Ry <- N.ext * cell.height
    
    xidx <- rep(1:M.ext,N.ext)
    yidx <- rep(1:N.ext,each=M.ext)
    dxidx <- pmin(abs(xidx-xidx[1]),M.ext-abs(xidx-xidx[1]))
    dyidx <- pmin(abs(yidx-yidx[1]),N.ext-abs(yidx-yidx[1]))
    D.ext.row1 <- ppp$D.ext.row1 <-  matrix(sqrt(((mcens.ext[2]-mcens.ext[1])*dxidx)^2+((ncens.ext[2]-ncens.ext[1])*dyidx)^2), M.ext, N.ext)
    
    A <- ppp$A <- cell.width * cell.height
    
    n.ji.ext <- matrix(NA, M.ext*N.ext, J)
    for(jj in 1:J)
    {
        for(j in 1:N.ext)
        {
            for(i in 1:M.ext)
            {
                ii = (j-1)*M.ext+i
                ii.old = (j-1)*M+i
                if(i <= M & j <=N)
                {
                    n.ji.ext[ii, jj] <- n.ji[ii.old, jj]
                }else
                {
                    n.ji.ext[ii, jj] <- 0
                }
            }
        }
    }
    ppp$n.ji.ext <- n.ji.ext
    
    ppp$xmin.ind <- xmin.ind <- min(which(ppp$window$xrange[1] <= mygrid$mgrid[-1]))
    ppp$xmax.ind <- xmax.ind <- min(which(ppp$window$xrange[2] <= mygrid$mgrid[-1]))
    ppp$ymin.ind <- ymin.ind <- min(which(ppp$window$yrange[1] <= mygrid$ngrid[-1]))
    ppp$ymax.ind <- ymax.ind <- min(which(ppp$window$yrange[2] <= mygrid$ngrid[-1]))
    
    obs.ind <- rep(0, nC.ext)
    
    if(object.only == TRUE)
    {
        if(exclude == T)
        {
            j.range <- 1:(J-1)
        }else
        {
            j.range <- 1:J
        }
        
        for(i in 1:M)
        {
            temp.ind <- NULL
            
            for(j in j.range)
            {
                temp.ind <- unique(c(temp.ind, which(n.ji.ext[c(0:(N-1))*M.ext+i, j] != 0)))
            }
            
            if(length(temp.ind) >0)
            {
                obs.ind[c(min(temp.ind):max(temp.ind)-1)*M.ext+i] <- 1
            }
            
            if(exclude == T)
            {
                exc.ind <- which(n.ji.ext[c(0:(N-1))*M.ext+i, J] != 0)
                ppp$n.ji.ext[(exc.ind-1)*M.ext+i, ] <- 0
                obs.ind[(exc.ind-1)*M.ext+i] <- 0
            }
        }
    }else
    {
        for(j in ymin.ind:ymax.ind)
        {
            obs.ind[(j-1)*M.ext+(xmin.ind:xmax.ind)] <- 1
        }
    }
    ppp$obs.ind <- obs.ind
    ppp$obs.inx <- which(obs.ind == 1)
    
    if(object.only == T & exclude==T)
    {
        ppp$type <- type[-J]
        ppp$J <- J-1
        ppp$n.ji <- ppp$n.ji[,-J]
        ppp$n.ji.ext <- ppp$n.ji.ext[,-J]
    }
    
    return(ppp)
}


## Internal function to extend the computional grid based on W, M, N, ext
grid.prep <- function(W, M, N, ext = 2)
{
    cell.width <- diff(W$xrange)/M
    cell.height <- diff(W$yrange)/N
    
    mgrid <- seq(W$xrange[1], W$xrange[2], by = cell.width)
    ngrid <- seq(W$yrange[1], W$yrange[2], by = cell.height)
    mcens <- (mgrid + 0.5 * cell.width)[-(M + 1)]
    ncens <- (ngrid + 0.5 * cell.height)[-(N + 1)]
    
    if (ext <= 1)
    mgrid.ext <- ngrid.ext <- mcens.ext <- ncens.ext <- M.ext <- N.ext <- NULL else {
        M.ext <- ext * M
        N.ext <- ext * N
        mgrid.ext <- seq(W$xrange[1], W$xrange[2] + (ext - 1) * diff(W$xrange), by = cell.width)
        ngrid.ext <- seq(W$yrange[1], W$yrange[2] + (ext - 1) * diff(W$yrange), by = cell.height)
        mcens.ext <- (mgrid.ext + 0.5 * cell.width)[-(M.ext + 1)]
        ncens.ext <- (ngrid.ext + 0.5 * cell.height)[-(N.ext + 1)]
    }
    
    return(list(M = M, N = N, mgrid = mgrid, ngrid = ngrid, mcens = mcens, ncens = ncens, cell.width = cell.width, cell.height = cell.height, M.ext = M.ext, N.ext = N.ext, mgrid.ext = mgrid.ext, ngrid.ext = ngrid.ext, mcens.ext = mcens.ext, ncens.ext = ncens.ext))
}






plot_data <- function(ppp, type)
{
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
    
    ow.mat <- matrix(ppp$obs.ind, length(ppp$mcens.ext), length(ppp$ncens.ext), byrow = F)[1:ppp$M, 1:ppp$N]
    
    jj = which(ppp$type == type)
    
    n.ji.mat.jj <- matrix(ppp$n.ji.ext[,jj], length(ppp$mygrid$mcens.ext), length(ppp$mygrid$ncens.ext), byrow = F)[1:ppp$M, 1:ppp$N]
    
    image.plot(ppp$mygrid$mcens, ppp$mygrid$ncens, ow.mat, ylim = c(0, max(ppp$mygrid$ncens)+ppp$cell.width/2), col = c("black", "white"), xlim = c(0, max(ppp$mygrid$mcens)+ppp$cell.width/2), main = type, zlim=c(0,1), xlab=NA, ylab=NA)
    
    M <- ppp$M
    N <- ppp$N
    
    for(j in 1:M) abline(v = ppp$mygrid$mgrid[j+1], lty = 1, col = "blue")
    for(j in 1:N) abline(h = ppp$mygrid$ngrid[j+1], lty = 1, col = "blue")
    
    points(ppp$x[ppp$marksNum == jj], ppp$y[ppp$marksNum == jj], cex = 0.2, col = jj, pch=10)
#    for(i in 1:ppp$M)
#    {
#        for(j in 1:ppp$N)
#        {
#            text(ppp$mygrid$mcens[i], ppp$mygrid$ncens[j], paste(n.ji.mat.jj[i,j]), cex = 0.5, col = "black")
#        }
#    }
}


