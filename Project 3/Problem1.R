# Project 3 in TMA4250 Spatial Statistics

#library(MASS)
#library(fields)
#library(akima)
#library(geoR)
#library(gridExtra)
#library(latex2exp)
#library(tidyverse)
#library(reshape2)

library(ggplot2)
library(rgdal)
library(spData)
library(spdep)
source('functions.R')


# seed
set.seed(4250)


### Problem 1

# read the data
load("data/Admin1Geography.RData")
# --> contains an object nigeriaAdm1 that contains the borders of the 37 admin1 areas

load("data/Admin2Geography.RData")
# --> contains an object nigeriaAdm2 that contains the borders of the 775 admin2 areas

# read the border graphs (note that a space " " was added in the beginning of each file)
table1 <- read.table("data/Admin1Graph.txt", header=TRUE, sep=" ")
table2 <- read.table("data/Admin2Graph.txt", header=TRUE, sep=" ")

