.generate_random <- function(x, nRows, numGEN)
{
    return(numGEN(nRows))
}

.generate_perm <- function(x, nRows)
{
    return(sample(x, size=nRows, replace=FALSE))
}

.parallel_random <- function(loadMat, emptyMat, nRows, numGEN, fn, corfn, ...)
{
    # loadMat is used to iterate with apply().
    emptyMat <- apply(emptyMat, 2, .generate_random, nRows, numGEN)
    eigenValues <- fn(corfn(emptyMat),...)$values
    return(eigenValues)
}

.parallel_perm <- function(loadMat, x, nRows, fn, corfn, ...)
{
    # loadMat is used to iterate with apply().
    x[,] <- apply(x, 2, .generate_perm, nRows)
    eigenValues <- fn(corfn(x),...)$values
    return(eigenValues)
}

parallel <- function(x, iter=1000, ordinal=FALSE, method=c("perm","random"), alpha=0.05, standard=FALSE, plot=TRUE, fn=eigen, ...)
{
    x <- as.matrix(x)
    nRows <- nrow(x)
    nComp <- ncol(x)
    method <- match.arg(method)
    if(!ordinal) {
        numGEN <- rnorm
        corfn <- cor
        correlation <- "pearson"
    } else {
        maxCateg <- max(x)-min(x)
        numGEN <- function(n) return(rbinom(n,maxCateg,0.5))
        if(maxCateg == 1) {
            corfn <- function(m) return(psych::tetrachoric(m)$rho)
            correlation <- "tetrachoric"
        } else {
            corfn <- function(m) return(psych::polychoric(m)$rho)
            correlation <- "polychoric"
        }
    }
    if(!ordinal & standard)
        x[,] <- apply(x, 2, function(x) as.numeric(scale(x)))
    compLabels <- paste0("c",1:nComp)
    randLoad <- matrix(nrow=nComp,ncol=iter)
    rownames(randLoad) <- compLabels
    eigenVal <- fn(corfn(x),...)$values
    names(eigenVal) <- compLabels
    if(method=="perm") {
        randLoad <- apply(randLoad, 2, .parallel_perm, x, nRows, fn, corfn, ...)
    } else {
        emptyMat <- matrix(nrow=nRows,ncol=nComp)
        randLoad <- apply(randLoad, 2, .parallel_random, emptyMat, nRows, numGEN, fn, corfn, ...)
    }
    randLoad <- t(randLoad)
    z <- qnorm(1-alpha/2)
    rM <- colMeans(randLoad)
    rSE <- apply(randLoad,2,mean)/sqrt(nRows)
    rQ <- apply(randLoad,2,quantile,probs=c(1-alpha,0.5))
    rCI <- rbind(rM+rSE,rM,rM-rSE)
    rownames(rCI) <- c("CI Sup","Mean","CI Inf")
    result <- list(
        correlation = correlation,
        method = method,
        synthetic.eigen = randLoad,
        pca.eigen = eigenVal,
        parallel.CI = rCI,
        parallel.quantiles = rQ,
        suggest.ncomp = sum(eigenVal > rQ[1,])
    )
    class(result) <- "parallel"
    if(plot) plot(result)
    return(result)
}

print.parallel <- function(x,...)
{
    x$method <- ifelse(x$method=="random","Random Data","Permutations")
    cat("Parallel Analysis by",x$method,"\n")
    cat("Correlation matrix:",x$correlation,"\n\n")
    cat("Observed values:\n")
    print(round(rbind("eigen"=x$pca.eigen),2))
    cat("\nParallel Confidence Intervals:\n")
    print(round(x$parallel.CI,2))
    cat("\nParallel Quantiles:\n")
    print(round(x$parallel.quantiles,2))
}

plot.parallel <- function(x,...)
{
    main <- list(...)$main
    xlab <- list(...)$xlab
    ylab <- list(...)$ylab
    ylim <- list(...)$ylim
    nComp <- list(...)$n.comp
    if(is.null(main))
        main <- "Parallel Analysis"
    if(is.null(xlab))
        xlab <- "Component number"
    if(is.null(ylab))
        ylab <- "Eigenvalue"
    if(is.null(nComp))
        nComp <- length(x$pca.eigen)
    component <- 1:nComp
    x$pca.eigen <- x$pca.eigen[component]
    x$parallel.CI <- x$parallel.CI[,component]
    x$parallel.quantiles <- x$parallel.quantiles[,component]
    Q <- x$parallel.quantiles[1,]
    xlim <- c(1,nComp)
    if(is.null(ylim)) {
        v <- c(x$pca.eigen,Q)
        ylim <- c(decimal.floor(min(v),digits=1),decimal.ceiling(max(v),digits=1))
    }
    col <- c("gray56","black")
    # Simulated data
    stripchart(Q~component,vertical=TRUE,xlim=xlim,ylim=ylim,col=col[1],pch=20,cex=1.7,cex.lab=1.1,main=main,xlab=xlab,ylab=ylab)
    segments(component[-nComp],Q[-nComp],component[-1],Q[-1],col=col[1],lwd=1.5)
    par(new=TRUE)
    # Observed data
    stripchart(x$pca.eigen~component,vertical=TRUE,xlim=xlim,ylim=ylim,xlab="",ylab="",pch=20,cex=1.7,col=col[2])
    segments(component[-nComp],x$pca.eigen[-nComp],component[-1],x$pca.eigen[-1],lwd=1.5)
    # Segment for eigen=1
    segments(0,1,nComp+1,1,lty=2)
}
