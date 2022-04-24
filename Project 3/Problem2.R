
source("Project 3/functions.R")
library(ggplot2)
library(rgdal)
library(spdep)
library(Matrix)

admin1 <- read.table("Project 3/data/DirectEstimates.txt", header = TRUE)
table1 <- read.table("Project 3/data/Admin1Graph.txt", header=TRUE, sep=" ")


# neighborhood matrices
options(max.print=1000000)
N1 <- as.matrix(table1[,-1])
colnames(N1)<-NULL
isSymmetric(N1)
# precision matrices
Q_besag <- diag(rowSums(N1)) - N1
Q_besag <- Matrix(Q_besag,sparse=TRUE)

# loads the map data into variable of name 'nigeriaAdm1'
load("Project 3/data/Admin1Geography.RData")

inverseLogit <- function(x) 1/(1+exp(-x))

obsProportions <- inverseLogit(admin1$Observation)

head(admin1)
summary(admin1)


# Function for getting conditional expectiations and precision matrix
getParam <- function(Q1,Q2, M,Y) {
    mu <- solve(Q1+t(M)%*%Q2%*%M, t(M)%*%Q2%*%Y)
    Q <- Q1 + t(M)%*%Q2%*%M
    list(mu=mu, Q=Q)
}

# Function for sampling from GMRF
sampleGMRF <- function(mu,Q, n_samp) {
    L <- chol(Q)
    z <- matrix( rnorm(37*n_samp,mean=0,sd=1), 37, n_samp) 
    v <- solve(L,z)
    x <- v+mu
}

# Function for calculating median, cv
calcStat <- function(x) {
    med_x <- apply(x,1,median)
    mean_x <- apply(x, 1, mean)
    sd_x <- apply(x, 1, sd)
    cv_x <- sd_x/mean_x
    low_cv <- floor(min(cv_x))
    high_cv <- ceiling(max(cv_x))
    list(med = med_x, cv = cv_x, low_cv = low_cv, high_cv= high_cv )
}

## Problem 2 a)

plotAreaCol(fName = "Project 3/figures/obsProp.pdf", width = 5, height = 4, estVal = obsProportions, 
            geoMap = nigeriaAdm1, leg = expression(~~hat(p)), colLim = c(0,1))


## Problem 2b)

#Sampler for x|y

sigma2 <- 100^2

Q2 <- Diagonal(x = 1/admin1$StdDev^2)  # ^2 or not?
Q1 <- (1/sigma2)*Diagonal(37)
Q <- Q1+Q2
mu <-  solve(Q1+Q2,Q2 %*% admin1$Observation) 
M <- diag(37)
paramsB <- getParam(Q1,Q2,M,admin1$Observation)

xB <- sampleGMRF(paramsB$mu,paramsB$Q, 100)
x_pB <- inverseLogit(xB)

statB <- calcStat(x_pB)

plotAreaCol(fName = "Project 3/figures/medianPost.pdf", width = 5, height = 4, estVal = statB$med, 
            geoMap = nigeriaAdm1, leg = "Median", colLim = c(0,1))
plotAreaCol(fName = "Project 3/figures/cvPost.pdf", width = 5, height = 4, estVal = statB$cv, 
            geoMap = nigeriaAdm1, leg = "CV", colLim = c(statB$low_cv,statB$high_cv))

## Problem 2c)

samp2C <- function(tau) {
Q_besag <- tau*Q_besag
paramsC <- getParam(Q_besag,Q2,diag(37),admin1$Observation)

rankMatrix(paramsC$Q)

xC <- sampleGMRF(paramsC$mu,paramsC$Q, 100)
x_pC <- inverseLogit(xC)

statC <- calcStat(x_pC)

plotAreaCol(fName = paste("Project 3/figures/medianPostC", tau, ".pdf", sep=""), width = 5, height = 4, estVal = statC$med, 
            geoMap = nigeriaAdm1, leg = "Median", colLim = c(0,1))
plotAreaCol(fName = paste("Project 3/figures/cvPostC", tau, ".pdf", sep=""), width = 5, height = 4, estVal = statC$cv, 
            geoMap = nigeriaAdm1, leg = "CV", colLim = c(statC$low_cv,statC$high_cv))

}

samp2C(tau=1)
samp2C(tau=0.1)
samp2C(tau=10)

## Problem 2d)

# Create M matrix:
ind_Kaduna <- grep("Kaduna", colnames(table1))
y_38 <- 0.5
M <- matrix(0, nrow = 38,ncol=37)
M[1:37,1:37] <- diag(37) 
M[38,ind_Kaduna] = 1
M <- Matrix(M, sparse=TRUE)

Q2D <- Diagonal(x = 1/c(admin1$StdDev^2,0.1^2))
paramsD <- getParam(Q_besag,Q2D,M,c(admin1$Observation,y_38))

xD <- sampleGMRF(paramsD$mu,paramsD$Q, 100)
x_pD <- inverseLogit(xD)

statD <- calcStat(x_pD)

plotAreaCol(fName = "Project 3/figures/medianPostD.pdf", width = 5, height = 4, estVal = statD$med, 
            geoMap = nigeriaAdm1, leg = "Median", colLim = c(0,1))
plotAreaCol(fName = "Project 3/figures/cvPostD.pdf", width = 5, height = 4, estVal = statD$cv, 
            geoMap = nigeriaAdm1, leg = "CV", colLim = c(statD$low_cv,statD$high_cv))

# Compare changes in neighbours
#statC$med[c(5,18,19,21)]
#statD$med[c(5,18,19,21)]

## Problem 2 f)


loglike <- function(tau) {
    paramf <- getParam(tau*Q_besag,Q2,diag(37),admin1$Observation)
    x <- rep(100, 37)
    -(((37-1)/2)*log(tau) -(tau/2)*t(x)%*%Q_besag%*%x -(1/2)*t(admin1$Observation -x)%*%Q2%*%(admin1$Observation -x)-
        (1/2)*log(det(paramf$Q))+(1/2)*t(x-paramf$mu)%*%paramf$Q%*%(x-paramf$mu))[1,1]
}

loglike2 <- function(tau) {
    paramf <- getParam(tau*Q_besag,Q2,diag(37),admin1$Observation)

    -(((37-1)/2)*log(tau) -(1/2)*t(admin1$Observation)%*%Q2%*%(admin1$Observation)-
            (1/2)*log(det(paramf$Q))+(1/2)*t(-paramf$mu)%*%paramf$Q%*%(-paramf$mu))[1,1]
}

tau_est <- optimize(loglike2,c(0.001,550)) 
tau_est

samp2C(tau_est$minimum)


