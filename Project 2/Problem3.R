# Project 2 in TMA4250 Spatial Statistics

library(MASS)
library(spatial)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(class)
library(factoextra)
##### Problem 3


### (a)

redwood <- read.table("data/redwood.dat", col.names = c("x", "y"))

# try clustering
#redwood.kmeans <- kmeans(redwood, 5)
#redwood$cluster <- redwood.kmeans$cluster
#redwood.cluster <- ggplot(redwood, aes(x=x, y=y)) +
#  geom_point(aes(colour=cluster), size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
#  labs(x=expression(x), y=expression(y)) + theme_bw() + scale_colour_manual(labels=redwood$cluster)
#redwood.cluster


# Neymann-Scott simulation
sim.parent <- function(lambda, add.dist = 0.25) {
  # update region to include additional areas around the domain
  ppregion(xl=-add.dist, xu=1+add.dist, yl=-add.dist, yu=1+add.dist)
  # volume
  volume = (1+2*add.dist)^2
  # get parent locations
  locs <- Psim(rpois(1,lambda*volume))
  # reset region domain
  ppregion(xl=0, xu=1, yl=0, yu=1)
  
  data.frame("x" = locs$x, "y" = locs$y)
}

sim.daughter <- function(parent.loc, lambda, sigmasq){
  # calculate number of daughter points
  n.daughters <- rpois(1, lambda)
  # return empty data frame if there are no daughters
  if(n.daughters == 0){
    return( data.frame("x" = double(), "y" = double()) )
  }
  # get daughter locations
  locs <- mvrnorm(n.daughters, as.numeric(parent.loc), sigmasq*diag(2))
  if(n.daughters == 1){
    return( data.frame("x" = locs[1], "y" = locs[2]) )
  }
  data.frame("x" = locs[,1], "y" = locs[,2])
}

sim.NeymanScott <- function(theta, add.dist = 0.25){
  # parameters
  lambda.p <- theta[1]
  lambda.d <- theta[2]
  sigmasq  <- theta[3]
  # data frame to store points
  points <- data.frame("x" = double(), "y" = double())
  # simulate parent locations
  parent.locs <- sim.parent(lambda.p, add.dist)
  for(i in 1:nrow(parent.locs)){
    # simulate daughter locations
    new.points <- sim.daughter(parent.locs[i,], lambda.d, sigmasq)
    if(nrow(new.points)){
      points <- rbind(points, new.points)
    }
  }
  # remove points outside the boundary
  return(points[!(points$x<0 | points$x>1 | points$y<0 | points$y>1), ])
}


## evaluating the empirical guess
# theta = (lambda.parent, lambda. daughter, sigma^2)
theta = c(8, 7.75, 0.05^2)

sim.NeymannScott.L <- function(n, theta, add.dist=0.25, fs=0.5, k=100) {
  simulations <- matrix(0, nrow=n, ncol=k)
  for (i in 1:n) {
    simulations[i, ] <- Kfn(pp=sim.NeymanScott(theta, add.dist), fs, k)$y
  }
  simulations
}

set.seed(4250)
NeymannScott.L <- sim.NeymannScott.L(100, theta)

redwood.L <- with(Kfn(pp=redwood, fs=0.5, k=100), data.frame("r"=x, "L"=y))
redwood.L$lower = apply(NeymannScott.L, 2, quantile, probs=0.05)
redwood.L$upper = apply(NeymannScott.L, 2, quantile, probs=0.95)

redwood.L.plot <- ggplot(redwood.L, aes(x=r, y=L)) +
  geom_line(aes(y=lower), colour="red") + geom_line(aes(y=upper), colour="red") +
  geom_ribbon(aes(ymin=lower,ymax=upper), fill="red", alpha=0.05) +
  geom_line(colour="blue") +
  labs(x=expression(r), y=expression(L(r))) + theme_bw()

ggsave("figures/redwood_NeymannScott_L.pdf", plot = redwood.L.plot, width = 3.0, height = 3.0)



## iterate to find theta

# estimate the difference between the two functions
theta0 = c(8, 7.75, 0.05^2) # initial value
n.iterations <- 10 # number of iterations
proposals <- 20 # number of proposed iterations of theta
n.sim <- 10 # number of simulations used to estimate L(r)
change <- 0.1 # percentage change in parameters



iterate.theta <- function(theta0, n.iterations, proposals, change, n.sim){
  # Empirical L
  redwood.L <- with(Kfn(pp=redwood, fs=0.5, k=100), data.frame("r"=x, "L"=y))
  
  # store iterations
  theta.vec <- matrix(NA, nrow = n.iterations+1, ncol = 3)
  theta.vec[1,] <- theta0
  
  # iterations
  for(i in 1:n.iterations){
    print(paste("Iteration number", i))
    
    # generate proposals for theta
    theta.proposals = matrix(rep(theta.vec[i,], proposals), nrow = proposals, byrow=TRUE) *
      runif(3*proposals, min=1-change, max=1+change)
    # force theta.d = 62/theta.d
    theta.proposals[,2] = 62/theta.proposals[,1]
    
    # find proposal with lowest objective
    objectives <- numeric(proposals)
    for(j in 1:proposals){
      # simulate for each proposal
      proposal.L <- colMeans( sim.NeymannScott.L(n.sim, theta.proposals[j,]) )
      objectives[j] <- sum((proposal.L-redwood.L)^2)
    }
    
    # find the lowest objective proposal
    theta.vec[i+1, ] <- theta.proposals[which.min(objectives), ]
  }
  return(theta.vec)
}

# iterating to find "better" values of theta 
theta.vec <- iterate.theta(theta0, n.iterations = 8, proposals = 20, change = 0.1, n.sim = 50)

new.theta <- theta.vec[9,]
# normalize such that lambda.d*lambda.p = 62

new.NeymannScott.L <- sim.NeymannScott.L(100, new.theta)

#redwood.L <- with(Kfn(pp=redwood, fs=0.5, k=100), data.frame("r"=x, "L"=y))
redwood.L$new.lower = apply(new.NeymannScott.L, 2, quantile, probs=0.05)
redwood.L$new.upper = apply(new.NeymannScott.L, 2, quantile, probs=0.95)

redwood.L.plot <- ggplot(redwood.L, aes(x=r, y=L)) +
  geom_line(aes(y=new.lower), colour="red") + geom_line(aes(y=new.upper), colour="red") +
  geom_ribbon(aes(ymin=new.lower,ymax=new.upper), fill="red", alpha=0.05) +
  geom_line(colour="blue") +
  labs(x=expression(r), y=expression(L(r))) + theme_bw()
redwood.L.plot

ggsave("figures/redwood_NeymannScott_L_new.pdf", plot = redwood.L.plot, width = 3.0, height = 3.0)



## Three new realizations using the iterated values of theta

set.seed(4250)
sim1 <- sim.NeymanScott(new.theta)
sim2 <- sim.NeymanScott(new.theta)
sim3 <- sim.NeymanScott(new.theta)


redwood.sim1 <- ggplot(data.frame(x=sim1$x, y=sim1$y), aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()
redwood.sim1

redwood.sim2 <- ggplot(data.frame(x=sim2$x, y=sim2$y), aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()
redwood.sim2

redwood.sim3 <- ggplot(data.frame(x=sim3$x, y=sim3$y), aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()
redwood.sim3

ggsave("figures/redwood_sim1.pdf", plot = redwood.sim1, width = 3.0, height = 3.0)
ggsave("figures/redwood_sim2.pdf", plot = redwood.sim2, width = 3.0, height = 3.0)
ggsave("figures/redwood_sim3.pdf", plot = redwood.sim3, width = 3.0, height = 3.0)

       
       
       
       