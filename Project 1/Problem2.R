library(MASS)
library(geoR)
library(akima)
library(fields)
library(ggplot2)
library(dplyr)
library(viridis)
library(plotly)
library(plot3D)

# plotting parameters
theme_set(theme_bw())
theme_update(text=element_text(size=14))

# loading data from the file topo.dat, and examining it
topo_data <- read.table("topo.dat")
topo_data <- as.data.frame(topo_data)
topo <- as.geodata(topo_data)
head(topo_data)
summary(topo_data)

#interpolate data for plotting 
interp_topo_data <- interp(topo_data$x,topo_data$y,topo_data$z)
interp_topo_data_df <- as.data.frame(cbind(interp_topo_data$x,interp_topo_data$y, interp_topo_data$z))
# construct different plots of the data 
pdf(file="imPlotTopo.pdf")
filled.contour(x = interp_topo_data$x, y=interp_topo_data$y, z=interp_topo_data$z, color.palette = turbo, xlab="s1", ylab="s2", cex.lab=1.5)
dev.off()

ggplot(topo_data, aes(x=x,y=y)) +
    geom_point(aes(color = z), size=4) +
    xlab("s1")+ylab("s2")+
    scale_color_viridis(name="X",option = "H") +
    theme(text = element_text(size = 14))   
ggsave("scatterPlotTopo.pdf", width = 5, height = 4.5)




# Construct grid for Kriging
grid <- expand.grid(1:315,1:315)

# Perform ordinary Kriging 
ordKrigPred <- krige.conv(topo, locations=grid, krige=krige.control(type.krige = "ok", trend.d = "cte", trend.l = "cte",
                                                                    , cov.model = "powered.exponential", cov.pars = c(2500,100), kappa=1.5)) 
interp_ordKrigPred <- interp(grid$Var1,grid$Var2,ordKrigPred$predict)
interp_ordKrigPredVar <- interp(grid$Var1,grid$Var2,ordKrigPred$krige.var)
# Display the prediction results and corresponding variance

pdf(file="contPlotKrig.pdf")
filled.contour(x = interp_ordKrigPred$x, y=interp_ordKrigPred$y, z=interp_ordKrigPred$z, color.palette = turbo, xlab="s1", ylab="s2", cex.lab=1.5)
dev.off()

pdf(file="contPlotKrigVar.pdf")
filled.contour(x = interp_ordKrigPredVar$x, y=interp_ordKrigPredVar$y, z=interp_ordKrigPredVar$z, color.palette = turbo, xlab="s1", ylab="s2", cex.lab=1.5)
dev.off()

pdf(file="3DKrig.pdf")
persp3D(x=interp_ordKrigPred$x,y=interp_ordKrigPred$y,z =interp_ordKrigPred$z, phi=40, theta=130, xlab="s1", ylab="s2", zlab="X",cex.lab=1.5)
dev.off()

pdf(file="3DKrigVar.pdf")
persp3D(x=interp_ordKrigPredVar$x,y=interp_ordKrigPredVar$y,z =interp_ordKrigPredVar$z, phi=40, theta=130, xlab="s1", ylab="s2", zlab="Var", cex.lab=1.5)
dev.off()

# Find specific values for task e)
ordKrigPred$predict[which(grid[,1]==100 & grid[,2]==100)]
ordKrigPred$krige.var[which(grid[,1]==100 & grid[,2]==100)]

# Perform universal Kriging
uniKrigPred <- krige.conv(topo, locations=grid, krige=krige.control(type.krige = "ok", trend.d = "2nd", trend.l = "2nd",
                                                                    , cov.model = "powered.exponential", cov.pars = c(2500,100), kappa=1.5)) 

interp_uniKrigPred <- interp(grid$Var1,grid$Var2,uniKrigPred$predict)
interp_uniKrigPredVar <- interp(grid$Var1,grid$Var2,uniKrigPred$krige.var)
#Display prediction and variance results

pdf(file="contPlotUniKrig.pdf")
filled.contour(x = interp_uniKrigPred$x, y=interp_uniKrigPred$y, z=interp_uniKrigPred$z, color.palette = turbo, xlab="s1", ylab="s2", cex.lab=1.5)
dev.off()

pdf(file="contPlotUniKrigVar.pdf")
filled.contour(x = interp_uniKrigPredVar$x, y=interp_uniKrigPredVar$y, z=interp_uniKrigPredVar$z, color.palette = turbo, xlab="s1", ylab="s2", cex.lab=1.5)
dev.off()

pdf(file="3DUniKrig.pdf")
persp3D(x=interp_uniKrigPred$x,y=interp_uniKrigPred$y,z =interp_uniKrigPred$z, phi=40, theta=130, xlab="s1", ylab="s2", zlab="X",cex.lab=1.5)
dev.off()

pdf(file="3DUniKrigVar.pdf")
persp3D(x=interp_uniKrigPredVar$x,y=interp_uniKrigPredVar$y,z =interp_uniKrigPredVar$z, phi=40, theta=130, xlab="s1", ylab="s2", zlab="Var", cex.lab=1.5)
dev.off()
