# ------------------------------------------------------------------------------------------
# PARALLEL PRINCIPAL COMPONENT ANALYSIS
# Script version: 1.0
# Author: Davide Massidda
# e-mail: davide.massidda@humandata.it
# Date: June 11, 2013
# URL: http://www.insular.it, http://www.humandata.it
# License: GPLv3
# Description: this script provides the function parallelPCA, which performs a parallel
# principal component analysis on the correlation matrix from continuous and ordinal data.
# Dependencies: psych
# ------------------------------------------------------------------------------------------
# USAGE
# parallelPCA(x, iter = 1000, ordinal = FALSE, method = "random",
#             alpha = 0.05, standard = FALSE, plot = TRUE)
# ARGUMENTS
# x: row data matrix.
# iter: the number of iterations for parallel analysis.
# ordinal: specifies if the data matrix contains ordinal data. If ordinal data are
#          present, tetrachoric or polycoric correlations are used instead of the
#          Pearson's indices. The number of categories for ordinal data is calculated
#          as max(x)-(min(x)-1).
# method: specifies if the parallel analysis must be performed using random data
#         (method="random") or permutations of observed data (method="perm").
# alpha: alpha level to calculate quantiles and confidence intervals.
# standard: specifies if the data must be standardized.
# plot: plots the scree test.
# VALUE
# The function returns and object of class 'parpca', containing the slots:
# @correlation: the type of correlation index used.
# @method: the type of data generation (random or permutations).
# @synthetic.eigen: the matrix of eigen values from each iteration of parallel analysis.
# @pca.eigen: estimated eigenvalues from observed data.
# @parallel.CI = averages and confidence intevals of eigenvalues estimated by parallel analysis.
# @parallel.quantiles = quantiles of eigenvalues estimated by parallel analysis.
# The plot method is implemented on the class 'parpca'. The function plot() applied
# to class 'parpca' returns the scree plot overlapping parallel quantiles (gray
# points and bars).
# ------------------------------------------------------------------------------------------
# REFERENCES
# [1] Buja A., Eyuboglu N. (1992). Remarks on parallel analysis. Multivariate Behavioral
#       Research, 27(4), 509–540.
# [2] Crawford A.V., Green S.B., Levy R., Lo W.J. Scott L., Svetina D., Thompson M.S. (2010)
#       Evaluation of parallel analysis methods for determining the number of factors.
#       Educational and Psychological Measurement, 70(6), 885–901.
# [3] Horn J.L. (1965). A rationale and test for the number of factors in factor analysis.
#       Psychometrika, 30, 179–185.
# [4] Franklin S.B., Gibson D.J., Robertson P.A., Pohlmann J.T., Fralish, J.S. (1995).
#       Parallel analysis: a method for determining significant principal components.
#       Journal of Vegetation Science, 6(1), 99-106.
# [5] Peres-Neto P.R., Jackson D.A., Somers K.M. (2005). How many principal components?
#       Stopping rules for determining the number of non-trivial axes revisited.
#       Computational Statistics and Data Analysis, 49, 974-997.
# [6] Weng L.J., Cheng C.P. (2005). Parallel analysis with unidimensional binary data.
#       Educational and Psychological Measurement, 65(5), 697-716.
# ------------------------------------------------------------------------------------------

require(psych)

setClass("parpca",
    representation(
        correlation = "character",
        method = "character",
        synthetic.eigen = "matrix",
        pca.eigen = "vector",
        parallel.CI = "matrix",
        parallel.quantiles = "matrix"
    )
)

setMethod("show","parpca",
    function(object) {
        object@method <- ifelse(object@method=="random","Random Data","Permutations")
        cat("Parallel Analysis by",object@method,"\n")
        cat("Correlation matrix:",object@correlation,"\n\n")
        cat("Observed values:\n")
        print(round(rbind("eigen"=object@pca.eigen),2))
        cat("\nParallel Confidence Intervals:\n")
        print(round(object@parallel.CI,2))
        cat("\nParallel Quantiles:\n")
        print(round(object@parallel.quantiles,2))
    }
)

setMethod("plot","parpca",
    function(object,x="parpca",y="parpca",main="Parallel Analysis",xlab="Component number",ylab="Eigenvalue")
    {
        nComp <- length(object@pca.eigen)
        component <- 1:nComp
        Q <- object@parallel.quantiles
        xlim <- c(1,nComp)
        ylim <- c(0,ceiling(max(c(object@pca.eigen,Q[1,]))))
        col <- c("gray56","black")    
        stripchart(object@pca.eigen~component,vertical=TRUE,xlim=xlim,ylim=ylim,xlab="",ylab="",pch=20,cex=1.5,col=col[2])
        segments(0,1,nComp+1,1,lty=2)
        segments(component[-nComp],object@pca.eigen[-nComp],component[-1],object@pca.eigen[-1])
        #textxy(1:nComp, object@pca.eigen, round(object@pca.eigen,2), cx=0.9)
        par(new=TRUE)
        stripchart(Q[2,]~component,vertical=TRUE,xlim=xlim,ylim=ylim,col=col[1],pch=20,cex=1.5,cex.lab=1.1,main=main,xlab=xlab,ylab=ylab)
        arrows(component, Q[1,], component, Q[3,], angle=90, code=3, col=col[1])
    }
)

parallelPCA <- function(x,iter=1000,ordinal=FALSE,method="random",alpha=0.05,standard=FALSE,plot=TRUE)
{
    x <- as.matrix(x)
    nRows <- dim(x)[1]
    nComp <- ncol(x)
    if(!ordinal) {
        numGEN <- rnorm
        corFUN <- cor
        correlation = "pearson"
    } else {
        maxCateg <- max(x)-min(x)
        numGEN <- function(n) return(rbinom(n,maxCateg,0.5))
        if(maxCateg == 1) {
            corFUN <- function(m,N) return(tetrachoric(m)$rho)
            correlation = "tetrachoric"
        } else {
            corFUN <- function(m,N) return(polychoric(m)$rho)
            correlation = "polychoric"
        }
    }
    if(!ordinal & standard) {
        m <- apply(x,2,mean,na.rm=TRUE)
        s <- apply(x,2,sd,na.rm=TRUE)
        for(i in 1:nComp)
            x[,i] <- (x[,i]-m[i])/s[i]
    }
    rowIndex <- 1:nRows
    xRand <- vector("list",iter)
    randLoad <- matrix(NA,ncol=nComp,nrow=iter)
    eigenVal <- eigen(corFUN(x),symmetric=TRUE,only.values=TRUE)$values
    colnames(randLoad) <- names(eigenVal) <- paste("c",1:nComp,sep="")
    if(method=="random") {
        for(i in 1:iter) {
            xRand[[i]] <- matrix(NA,nrow=nRows,ncol=nComp)
            for(j in 1:nComp) xRand[[i]][,j] <- numGEN(nRows)
            randLoad[i,] <- eigen(corFUN(xRand[[i]]),symmetric=TRUE,only.values=TRUE)$values
        }
    } else {
        if(method=="perm") {
            N <- length(x)
            K <- 1:N
            V <- c(x)
            for(i in 1:iter) {
                xRand[[i]] <- matrix(V[sample(K,size=N,replace=FALSE)],nrow=nRows,ncol=nComp)
                randLoad[i,] <- eigen(corFUN(xRand[[i]]),symmetric=TRUE,only.values=TRUE)$values
            }
        }
    }
    z <- abs(qnorm(alpha/2))
    rM <- colMeans(randLoad)
    rSE <- apply(randLoad,2,mean)/sqrt(nRows)
    rQ <- apply(randLoad,2,quantile,probs=c(1-alpha/2,0.5,alpha/2))
    rCI <- rbind(rM+rSE,rM,rM-rSE)
    rownames(rCI) <- c("CI Sup","Mean","CI Inf")
    result <- new("parpca",
        correlation = correlation,
        method = method,
        synthetic.eigen = randLoad,
        pca.eigen = eigenVal,
        parallel.CI = rCI,
        parallel.quantiles = rQ
    )
    if(plot) plot(result)
    return(result)
}
