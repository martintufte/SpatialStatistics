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
table1 <- read.table("data/Admin1Graph.txt", header=TRUE, sep=" ")
table2 <- read.table("data/Admin2Graph.txt", header=TRUE, sep=" ")

# neighborhood matrices
M1 <- as.matrix(table1[,-1])
M2 <- as.matrix(table2[,-1])

# precision matrices
Q1 <- diag(rowSums(M1)) - M1
Q2 <- diag(rowSums(M2)) - M2

# find ranks
rankMatrix(Q1)
rankMatrix(Q2)

# sparsity
print(paste("Matrix Q1 has sparsity", round(100*sum(Q1 == 0)/length(Q1), 2), "%." ))
print(paste("Matrix Q1 has sparsity", round(100*sum(Q2 == 0)/length(Q2), 2), "%." ))


## b)

plt.sparsity <- function(Q){
  x <- 1:dim(Q)[1]
  y <- 1:dim(Q)[2]
  
  Q.flipped <- apply(Q, 1, rev)
  return(image(x, y, Q.flipped!=0, xaxt='n', ann=FALSE, asp=1, yaxt='n', xlab="Column", ylab="Row", col=c("white", "black")))
}

plt.sparsity(Q1)
plt.sparsity(Q2)

install.packages("ggmap")
library(ggmap)

ggimage(Q1)
?ggimage
