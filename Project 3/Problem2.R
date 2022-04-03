
source("Project 3/functions.R")
library(ggplot2)
library(rgdal)
library(spdep)

admin1 <- read.table("Project 3/data/DirectEstimates.txt", header = TRUE)
#admin1Map <- read.table("Project 3/data/admin1Graph.txt", header=TRUE, sep=" ")

# loads the map data into variable of name 'nigeriaAdm1'
load("Project 3/data/Admin1Geography.RData")

inverseLogit <- function(x) 1/(1+exp(-x))

obsProportions <- inverseLogit(admin1$Observation)

head(admin1)
summary(admin1)

lowLim <- floor(min(obsProportions))
uppLim <- ceiling(max(obsProportions))

## Problem 2 a)

plotAreaCol(fName = "Project 3/figures/obsProp.pdf", width = 5, height = 4, estVal = obsProportions, 
            geoMap = nigeriaAdm1, leg = expression(~~hat(p)), colLim = c(lowLim,uppLim))


## Problem 2 b)






