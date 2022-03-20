
library(spatial)
library(sf)
library(ggplot2)
library(viridis)
theme_set(theme_bw())

#load the required data
obsprob <- read.table("Project 2/data/obsprob.txt", header=TRUE) #x,y, alpha
obspines <- read.table("Project 2/data/obspines.txt", header=TRUE) #x,y, N_obs
obs_grid <- data.frame(x=obsprob$x, y=obsprob$y)
alpha <- obsprob$alpha
N_obs <- obspines$N_obs

num_cells <- length(N_obs)

new_pines <- st_as_sf(obspines, coords = c("x","y"))

area_cells <- 100
cell_size <- 10

new_data <- (st_bbox(new_pines) + cell_size/2*c(-1,-1,1,1)) %>%
    st_make_grid(cellsize=c(cell_size, cell_size)) %>% st_sf()
new_data$N_obs = N_obs
new_data$alpha = alpha

### Task a)

ggplot() + 
    geom_sf(data=new_data, aes(geometry=geometry, fill=factor(N_obs))) + 
    xlab("x") +
    ylab("y") +
    scale_fill_viridis_d(name=" M") +
    coord_sf()
ggsave("Project 2/figures/obspines.pdf", width = 5, height = 4.5)

ggplot() + 
    geom_sf(data=new_data, aes(geometry=geometry, fill=alpha)) + 
    xlab("x") +
    ylab("y") +
    scale_fill_viridis(name=expression(~~alpha),option = "inferno") +
    coord_sf()
ggsave("Project 2/figures/obsprob.pdf", width = 5, height = 4.5)

#### Task c)

#Estimation of lambda based on estimator from 2c)
lambda_hat <- (1/(area_cells*sum(alpha)))*sum(N_obs)
lambda_hat

num_simsC <- 3

# Function for simulating placement of points
place_points <- function(x,y,N) {
    n_t <- sum(N)
    x_new <- rep(NA,n_t)
    y_new <- rep(NA,n_t)
    counter <- 1
    if (n_t>0) {
        for (i in (1:length(N))) {
            if (N[i]>0) {
                x_new[counter:(counter + N[i] -1)] = runif(N[i],x[i]-5,x[i]+5)
                y_new[counter:(counter + N[i] -1)] = runif(N[i],y[i]-5,y[i]+5)  
                counter <- counter + N[i]
            }
        }
        data.frame(x=x_new,y=y_new)
    }
}

# Create matrix for storing sim values, and list to hold placement of points
N_true <- matrix(NA,nrow=num_cells, ncol = num_simsC)
placement_points <- list()

for (i in 1:num_simsC) {
    N_true[, i] <- rpois(num_cells,lambda_hat*area_cells)
    placement_points[[i]] <- place_points(obs_grid$x,obs_grid$y,N_true[,i])
}

length(placement_points[[1]]$x)
sum(N_true[,1])

for (i in 1:num_simsC) {
    ggplot() +
            geom_sf(data=new_data, aes(geometry=geometry)) +
            geom_point(data=placement_points[[i]], aes(x=x,y=y))
    ggsave(file=paste("Project 2/figures/N_true", i,".pdf", sep=""), width = 5, height = 5)
}

#### Task d)
#set.seed(1)

num_simsD <- 3

N_post <- matrix(NA,nrow=num_cells, ncol = num_simsD)
placement_points_post <- list()
for (i in 1:num_simsD) {
    N_post[, i] <- rpois(num_cells,(1-alpha)*lambda_hat*area_cells)
    N_post[,i] <- N_post[,i] + N_obs
    placement_points_post[[i]] <- place_points(obs_grid$x,obs_grid$y,N_post[,i])
}
length(placement_points_post[[1]]$x)
sum(N_post[,1])
for (i in 1:num_simsC) {
    ggplot() +
        geom_sf(data=new_data, aes(geometry=geometry)) +
        geom_point(data=placement_points_post[[i]], aes(x=x,y=y))
    ggsave(file=paste("Project 2/figures/N_post", i,".pdf", sep=""), width = 5, height = 5)
}


#### Task e)

num_simsE = 500

N_trueE <- matrix(NA,nrow=num_cells, ncol = num_simsE)
for (i in 1:num_simsE) {
    N_trueE[, i] <- rpois(num_cells,lambda_hat*area_cells)
}


N_postE <- matrix(NA,nrow=num_cells, ncol = num_simsE)
for (i in 1:num_simsE) {
    N_postE[, i] <- rpois(num_cells,(1-alpha)*lambda_hat*area_cells)
    N_postE[,i] <- N_postE[,i] + N_obs
}

N_true_mean <- rowMeans(N_trueE)
N_true_sd <- apply(N_trueE, MARGIN = 1, sd)

N_post_mean <- rowMeans(N_postE)
N_post_sd <- apply(N_postE, MARGIN = 1, sd)

#find maximum value for color bar
max_mean <-ceiling(max(max(N_true_mean),max(N_post_mean)))
max_sd <- ceiling(max(max(N_true_sd),max(N_post_sd)))

new_data$N_post_mean <- N_post_mean
new_data$N_true_mean <- N_true_mean
new_data$N_true_sd <- N_true_sd
new_data$N_post_sd <- N_post_sd

ggplot() + 
    geom_sf(data=new_data, aes(geometry=geometry, fill=(N_post_mean))) + 
    xlab("x") +
    ylab("y") +
    scale_fill_viridis(name="mean", limits=c(0,max_mean),option="inferno") +
    coord_sf()
#ggsave(file=paste("Project 2/figures/N_post_mean.pdf", sep=""), width = 5, height = 5)

ggplot() + 
    geom_sf(data=new_data, aes(geometry=geometry, fill=(N_true_mean))) + 
    xlab("x") +
    ylab("y") +
    scale_fill_viridis(name="mean",limits=c(0,max_mean),option="inferno") +
    coord_sf()
#ggsave(file=paste("Project 2/figures/N_true_mean.pdf", sep=""), width = 5, height = 5)

ggplot() + 
    geom_sf(data=new_data, aes(geometry=geometry, fill=(N_true_sd))) + 
    xlab("x") +
    ylab("y") +
    scale_fill_viridis(name=" sd",limits=c(0,max_sd),option="inferno") +
    coord_sf()
#ggsave(file=paste("Project 2/figures/N_true_sd.pdf", sep=""), width = 5, height = 5)

ggplot() + 
    geom_sf(data=new_data, aes(geometry=geometry, fill=(N_post_sd))) + 
    xlab("x") +
    ylab("y") +
    scale_fill_viridis(name=" sd",limits=c(0,max_sd),option="inferno") +
    coord_sf()
#ggsave(file=paste("Project 2/figures/N_post_sd.pdf", sep=""), width = 5, height = 5)
