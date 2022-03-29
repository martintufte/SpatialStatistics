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

# seed
set.seed(4250)


### Problem 1

# read the data
load("data/Admin1Geography.RData")
# --> contains an object nigeriaAdm1 that contains the borders of the 37 admin1 areas

load("data/Admin2Geography.RData")
# --> contains an object nigeriaAdm2 that contains the borders of the 775 admin2 areas



