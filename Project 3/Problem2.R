
#source("Project 3/functions.R")
library(ggplot2)
library(rgdal)
library(spdep)

admin1 <- read.table("Project 3/data/DirectEstimates.txt", header = TRUE)
#admin1Map <- read.table("Project 3/data/admin1Graph.txt", header=TRUE, sep=" ")
load("Project 3/data/Admin1Geography.RData")
head(admin1)
summary(admin1)

plotAreaCol(fName = "Project 3/figures/obsProp.pdf", estVal = admin1$Observation,geoMap = nigeriaAdm1)
