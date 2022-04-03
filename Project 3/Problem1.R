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
N1 <- as.matrix(table1[,-1])
N2 <- as.matrix(table2[,-1])

# precision matrices
Q1 <- diag(rowSums(N1)) - N1
Q2 <- diag(rowSums(N2)) - N2

# find ranks
rankMatrix(Q1)
rankMatrix(Q2)

# sparsity
print(paste("Matrix Q1 has sparsity", round(100*sum(Q1 == 0)/length(Q1), 2), "%." ))
print(paste("Matrix Q2 has sparsity", round(100*sum(Q2 == 0)/length(Q2), 2), "%." ))


## b)

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

sparsity1 <- plt.sparsity(Q1)
sparsity2 <- plt.sparsity(Q2)

ggsave("./figures/sparsityQ1.pdf", plot = sparsity1, width = 5, height = 5)
ggsave("./figures/sparsityQ2.pdf", plot = sparsity2, width = 5, height = 5)












