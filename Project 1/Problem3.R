library(MASS)
library(geoR)
library(akima)
library(fields)
library(ggplot2)
library(dplyr)

# Plotting parameters
theme_set(theme_bw())
theme_update(text=element_text(size=14))

# Set parameters for the GRF
sigma2 <- 2
a <- 3

# Construct a grid of distances to evaluate the semi-variogram fnc. on, and calculate the true s-variogram
gridh <- 0:43
true_svariogram <- sigma2*(1-exp(-gridh/a))
true_svariogram_df <- as.data.frame(cbind(gridh,true_svariogram))
grid <- expand.grid(1:30,1:30)

# Genereate three realizations of X on the grid, plot it and then estimate the variogram and plot
for (i in 1:3) {
    set.seed(2022+i)
    x_realization <- grf(grid=grid, cov.model = "exponential", cov.pars = c(sigma2,a))
    svariogram_est <- variog(x_realization)
    df1 <- as.data.frame(cbind(grid, x_realization$data))
    colnames(df1) <- c("s1","s2", "realization")
    df2 <- as.data.frame(cbind(svariogram_est$u, svariogram_est$v)) %>% setNames(c("u","v"))
    
    #plot realization
    ggplot(data = df1, aes(s1,s2)) +
        geom_raster(aes(fill = realization)) +
        scale_fill_viridis_c(name="X") +
        coord_fixed() + xlab("s1") + ylab("s2")+
        ggtitle(paste("Realization", i))
    ggsave(filename = paste("Real",i,".pdf", sep=""), scale=0.66, device="pdf")
    
    #plot estimated semivariogram and true semivariogram
    ggplot() +
        geom_line(data=df2, aes(x=u,y=v, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
        geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
        scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
        scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
        ggtitle(paste("Realization", i)) +
        theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
              legend.background = element_rect(fill = "white", colour = "gray30")) + 
        xlab("Distance") +
        ylab("Semi variogram")
        
    ggsave(filename = paste("SVario",i,".pdf",sep=""), scale=0.66, device="pdf")
}

set.seed(2060)

#Generate a realization for use in the following tasks
x_realization <- grf(grid=grid, cov.model = "exponential", cov.pars = c(sigma2,a))

# Case 36 samples (Identical method for subsequent cases)

# Draw random a subgrid
sub36_grid_ind <- sample(1:900, 36, replace=FALSE)
sub36_grid <- grid[sub36_grid_ind,]

# Estimate the s-variogram on the subgrid
svariogram_sub36est <- variog(coords = sub36_grid, data=x_realization$data[sub36_grid_ind])
df32 <- as.data.frame(cbind(svariogram_sub36est$u, svariogram_sub36est$v)) %>% setNames(c("u","v"))
#plot the semivariogram and the true semivariogram
ggplot() +
    geom_line(data=df32, aes(x=u,y=v, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
    geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
    scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
    scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
    ggtitle(paste("Subgrid, 36 locations")) +
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    xlab("Distance") +
    ylab("Semi variogram")

ggsave(filename = "SVarioSub36.pdf", scale=0.66, device="pdf")

#Now we want to estimate sigma2 and a based on full grid and subgrid
full_param_est <- likfit(x_realization, cov.model = "exponential", ini.cov.pars =c(1,1.5))
full_param_est$phi
full_param_est$sigmasq

sub36_param_est <- likfit(coords = sub36_grid, data=x_realization$data[sub36_grid_ind], cov.model="exponential", ini.cov.pars = c(1,1.5))
sub36_param_est$phi
sub36_param_est$sigmasq

#Case 9 samples

sub9_grid_ind <- sample(1:900, 9, replace=FALSE)
sub9_grid <- grid[sub9_grid_ind,]
svariogram_sub9est <- variog(coords = sub9_grid, data=x_realization$data[sub9_grid_ind])
df9 <- as.data.frame(cbind(svariogram_sub9est$u, svariogram_sub9est$v)) %>% setNames(c("u","v"))
ggplot() +
    geom_line(data=df9, aes(x=u,y=v, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
    geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
    scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
    scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
    ggtitle(paste("Subgrid, 9 locations")) +
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    xlab("Distance") +
    ylab("Semi variogram")

ggsave(filename = "SVarioSub9.pdf", scale=0.66, device="pdf")

sub9_param_est <- likfit(coords = sub9_grid, data=x_realization$data[sub9_grid_ind], cov.model="exponential", ini.cov.pars = c(1,1.5))
sub9_param_est$phi
sub9_param_est$sigmasq

#Case 64 samples

sub64_grid_ind <- sample(1:900, 64, replace=FALSE)
sub64_grid <- grid[sub64_grid_ind,]
svariogram_sub64est <- variog(coords = sub64_grid, data=x_realization$data[sub64_grid_ind])
df64 <- as.data.frame(cbind(svariogram_sub64est$u, svariogram_sub64est$v)) %>% setNames(c("u","v"))
ggplot() +
    geom_line(data=df64, aes(x=u,y=v, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
    geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
    scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
    scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
    ggtitle(paste("Subgrid, 64 locations")) +
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    xlab("Distance") +
    ylab("Semi variogram")

ggsave(filename = "SVarioSub64.pdf", scale=0.66, device="pdf")

sub64_param_est <- likfit(coords = sub64_grid, data=x_realization$data[sub64_grid_ind], cov.model="exponential", ini.cov.pars = c(1,1.5))
sub64_param_est$phi
sub64_param_est$sigmasq

#Case 100 samples

sub100_grid_ind <- sample(1:900, 100, replace=FALSE)
sub100_grid <- grid[sub100_grid_ind,]
svariogram_sub100est <- variog(coords = sub100_grid, data=x_realization$data[sub100_grid_ind])
df100 <- as.data.frame(cbind(svariogram_sub100est$u, svariogram_sub100est$v)) %>% setNames(c("u","v"))
ggplot() +
    geom_line(data=df100, aes(x=u,y=v, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
    geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
    scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
    scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
    ggtitle(paste("Subgrid, 100 locations")) +
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    xlab("Distance") +
    ylab("Semi variogram")

ggsave(filename = "SVarioSub100.pdf", scale=0.66, device="pdf")

sub100_param_est <- likfit(coords = sub100_grid, data=x_realization$data[sub100_grid_ind], cov.model="exponential", ini.cov.pars = c(1,1.5))
sub100_param_est$phi
sub100_param_est$sigmasq

semVarFnc <- function(sigma2,a) {
    sigma2*(1-exp(-gridh/a))
}

#lastly we plot the semi-variogram from estimated parameters compared to true semivariogram
sv9 <- semVarFnc(sub9_param_est$sigmasq,sub9_param_est$phi)
df9 <- data.frame(x = gridh, y=sv9)
ggplot() +
    geom_line(data=df9, aes(x=x,y=y, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
    geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
    scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
    scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
    ggtitle(paste("Subgrid, 9 locations")) +
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    xlab("Distance") +
    ylab("Semi variogram")
ggsave(filename = "SVarioSub9ParEst.pdf", scale=0.66, device="pdf")

sv36 <- semVarFnc(sub36_param_est$sigmasq,sub36_param_est$phi)
df36 <- data.frame(x = gridh, y=sv36)
ggplot() +
    geom_line(data=df36, aes(x=x,y=y, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
    geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
    scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
    scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
    ggtitle(paste("Subgrid, 36 locations")) +
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    xlab("Distance") +
    ylab("Semi variogram")
ggsave(filename = "SVarioSub36ParEst.pdf", scale=0.66, device="pdf")

sv64 <- semVarFnc(sub64_param_est$sigmasq,sub64_param_est$phi)
df64 <- data.frame(x = gridh, y=sv64)
ggplot() +
    geom_line(data=df64, aes(x=x,y=y, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
    geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
    scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
    scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
    ggtitle(paste("Subgrid, 64 locations")) +
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    xlab("Distance") +
    ylab("Semi variogram")
ggsave(filename = "SVarioSub64ParEst.pdf", scale=0.66, device="pdf")

sv100 <- semVarFnc(sub100_param_est$sigmasq,sub100_param_est$phi)
df100 <- data.frame(x = gridh, y=sv100)
ggplot() +
    geom_line(data=df100, aes(x=x,y=y, color="Estimated semivariogram", linetype="Estimated semivariogram")) +
    geom_line(data=true_svariogram_df, aes(x=gridh, y=true_svariogram, color="True semivariogram", linetype="True semivariogram")) +
    scale_color_manual(name="", values = c("Estimated semivariogram" = "blue", "True semivariogram" = "black")) +
    scale_linetype_manual(name="", values = c("Estimated semivariogram"= "solid", "True semivariogram" = "longdash")) +
    ggtitle(paste("Subgrid, 100 locations")) +
    theme(legend.direction = "horizontal", legend.position = "bottom", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    xlab("Distance") +
    ylab("Semi variogram")
ggsave(filename = "SVarioSub100ParEst.pdf", scale=0.66, device="pdf")