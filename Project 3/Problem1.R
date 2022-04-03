# Project 3 in TMA4250 Spatial Statistics

#library(MASS)
#library(fields)
#library(akima)
#library(geoR)
#library(gridExtra)
#library(latex2exp)
#library(tidyverse)
#library(reshape2)

library(Matrix)
library(ggplot2)
library(rgdal)
library(spData)
library(spdep)
library(MASS)

source('functions.R')

# seed
set.seed(4250)


### Problem 1

## a)

# read the data
load("data/Admin1Geography.RData")
# --> contains an object nigeriaAdm1 that contains the borders of the 37 admin1 areas

load("data/Admin2Geography.RData")
# --> contains an object nigeriaAdm2 that contains the borders of the 775 admin2 areas

# read the border graphs (note that a space " " was added in the beginning of each file)
admin1 <- read.table("data/Admin1Graph.txt", header=TRUE, sep=" ")
admin2 <- read.table("data/Admin2Graph.txt", header=TRUE, sep=" ")

# neighborhood matrices
N1 <- as.matrix(admin1[,-1])
N2 <- as.matrix(admin2[,-1])

# precision matrices
tau1 <- 1
tau2 <- 1
Q1 <- tau1 * (diag(rowSums(N1)) - N1)
Q2 <- tau2 * (diag(rowSums(N2)) - N2)

# find ranks
rankMatrix(Q1)
rankMatrix(Q2)

# sparsity
print(paste("Matrix Q1 has sparsity", round(100*sum(Q1 == 0)/length(Q1), 2), "%." ))
print(paste("Matrix Q2 has sparsity", round(100*sum(Q2 == 0)/length(Q2), 2), "%." ))

plt.sparsity <- function(Q){
  x <- 1:dim(Q)[1]
  y <- 1:dim(Q)[2]-1
  
  df <- data.frame(Nonzero = array(Q!=0), coor = expand.grid(x,y))
  
  plot <- ggplot(df, aes(x=coor.Var1, y=coor.Var2, fill=Nonzero)) +
    geom_raster(hjust = 0, vjust = 0) + xlab("Column") + ylab("Row") +
    scale_fill_manual(values=c("white", "black")) + scale_y_reverse() +
    coord_fixed() + theme_classic() + theme(legend.position="none")
  
  return(plot)
}

ggsave("./figures/sparsityQ1.pdf", plot = plt.sparsity(Q1), width = 4, height = 4)
ggsave("./figures/sparsityQ2.pdf", plot = plt.sparsity(Q2), width = 4, height = 4)


## b)
rGMRF <- function(n, Q, eps = 1e-8){
  nr <- nrow(Q)
  L = chol(Q + eps*diag(nr))
  
  if(n==1){
    z <- rnorm(nr)
    v <- solve(L, z)
    return(v - mean(v))
  } else {
    #V <- matrix(rep(0, nr*n), nrow=nr, ncol=n)
    Z <- matrix(rnorm(nr*n), nrow=nr)
    V <- solve(L, Z)
    return( V - matrix(rep(colMeans(V), nr), nrow=nr, byrow=TRUE) )
  }
}

set.seed(4250)
admin1.besag1 <- rGMRF(1, Q1, eps=1e-8)
admin1.besag2 <- rGMRF(1, Q1, eps=1e-8)
admin1.rnorm1 <- rnorm(nrow(Q1))
admin1.rnorm2 <- rnorm(nrow(Q1))

# same scale to all plots
lim <- ceiling( max(abs(c(admin1.besag1, admin1.besag2, admin1.rnorm1, admin1.rnorm2))) )

# display two realizations from beasg model
plotAreaCol(fName="./figures/admin1_besag_realization1.pdf", width=5, height=4, estVal=admin1.besag1,
            geoMap=nigeriaAdm1, leg=expression(x), colLim = c(-lim,lim))
plotAreaCol(fName="./figures/admin1_besag_realization2.pdf", width=5, height=4, estVal=admin1.besag2,
            geoMap=nigeriaAdm1, leg=expression(x), colLim = c(-lim,lim))

# display two realizations from independent model
plotAreaCol(fName="./figures/admin1_rnorm_realization1.pdf", width=5, height=4, estVal=admin1.rnorm1,
            geoMap=nigeriaAdm1, leg=expression(x), colLim = c(-lim,lim))
plotAreaCol(fName="./figures/admin1_rnorm_realization2.pdf", width=5, height=4, estVal=admin1.rnorm2,
            geoMap=nigeriaAdm1, leg=expression(x), colLim = c(-lim,lim))


## c)

set.seed(4250)
admin2.besag1 <- rGMRF(1, Q2, eps=1e-8)
admin2.besag2 <- rGMRF(1, Q2, eps=1e-8)
admin2.rnorm1 <- rnorm(nrow(Q2))
admin2.rnorm2 <- rnorm(nrow(Q2))

# same scale to all plots
lim <- ceiling( max(abs(c(admin2.besag1, admin2.besag2, admin2.rnorm1, admin2.rnorm2))) )

# display two realizations from beasg model
plotAreaCol(fName="./figures/admin2_besag_realization1.pdf", width=5, height=4, estVal=admin2.besag1,
            geoMap=nigeriaAdm2, leg=expression(x), colLim = c(-lim,lim))
plotAreaCol(fName="./figures/admin2_besag_realization2.pdf", width=5, height=4, estVal=admin2.besag2,
            geoMap=nigeriaAdm2, leg=expression(x), colLim = c(-lim,lim))

# display two realizations from independent model
plotAreaCol(fName="./figures/admin2_rnorm_realization1.pdf", width=5, height=4, estVal=admin2.rnorm1,
            geoMap=nigeriaAdm2, leg=expression(x), colLim = c(-lim,lim))
plotAreaCol(fName="./figures/admin2_rnorm_realization2.pdf", width=5, height=4, estVal=admin2.rnorm2,
            geoMap=nigeriaAdm2, leg=expression(x), colLim = c(-lim,lim))


## d)
set.seed(4250)
admin2.100 <- rGMRF(100, Q2, eps=1e-8)

# empirical variance
admin2.var <- sapply(as.data.frame(t(admin2.100)), var)

max_var <- ceiling(max(admin2.var))

plotAreaCol(fName="./figures/admin2_variance.pdf", width=5, height=4, estVal=admin2.var,
            geoMap=nigeriaAdm2, leg="var", colLim = c(0,max_var))

ne <- diag(Q2/tau2)

# plot empirical variance with number of neighbors
median <- rep(0, nrow(Q2))
for (u in unique(ne)){
  median[ne == u] <- median(admin2.var[ne==u])
}

ne.var <- ggplot(data.frame(ne=ne, var=admin2.var, median_var=median), aes(x=ne, y=var, group=ne, fill=median)) +
  geom_boxplot() + labs(x="Neighbors", y="Empirical variance") +
  scale_fill_viridis_c(direction=1, begin=1, end=0, limit=c(0,2)) + theme_classic()
ne.var

ggsave("./figures/admin2_variance_boxplot.pdf", plot = ne.var, width = 5, height = 4)


# empirical correlation with Gubio
cor.Gubio <- function(x) cor(admin2.100[150,], x)

admin2.corr <- sapply(as.data.frame(t(admin2.100)), cor.Gubio)

plotAreaCol(fName="./figures/admin2_correlation_Gubio.pdf", width=5, height=4, estVal=admin2.corr,
            geoMap=nigeriaAdm2, leg="corr", colLim = c(-1,1))

