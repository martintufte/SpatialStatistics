# Project 2 in TMA4250 Spatial Statistics

library(spatial)
library(ggplot2)
library(Rfast)

##### Problem 4


### (a)

# Read data
cells <- read.table("data/cells.dat", col.names = c("x", "y"))
n.points <- nrow(cells)

# Energy in a state using a Strauss process
potential.energy <- function(state, beta, r0){
  beta * sum(dist(state) < r0)
}

# Acceptance probability
acceptance.prob <- function(current.state, proposal.x_k, k, beta, r0){
  current.x_k <- matrix(current.state[k,], nrow = 1)
  
  # number of points within pairwise distance r0
  current.within.r0  <- sum(dista(current.x_k, current.state[-k,]) < r0)
  proposal.within.r0 <- sum(dista(proposal.x_k, current.state[-k,]) < r0)
  
  # return acceptance probability
  min(1, exp(beta*(current.within.r0 - proposal.within.r0)))
}

# Simulate the MCMC
sim.MCMC <- function(n.points, n.iter, beta, r0){
  n.accepted  <- 0
  energies    <- numeric(1+n.iter)
  states      <- array(rep(NA, (n.iter+1)*n.points*2), dim=c(n.iter+1, n.points, 2))
  
  # initial values
  states[1,,] <- matrix(runif(2*n.points), nrow = n.points)
  energies[1] = potential.energy(states[1,,], beta, r0)
  
  # iterations
  for(i in 1:n.iter){
    # set new state as current state
    states[i+1,,] <- states[i,,]
    
    # try updating each point
    for(k in 1:n.points){
      
      # propose new x_k
      proposal.x_k <- matrix(runif(2), nrow=1)
      
      # calculate acceptance probability
      alpha <- acceptance.prob(states[i+1,,], proposal.x_k, k, beta, r0)
      if(runif(1) < alpha){
        states[i+1,k,] <- proposal.x_k
        n.accepted <- n.accepted + 1
      }
    }
    # store energy
    energies[i+1] <- potential.energy(states[i+1,,], beta, r0)
  }
  list("state" = states[n.iter+1,,],
       "acceptance.rate" = n.accepted/(n.iter*n.points),
       "energies" = energies)
}


sim.MCMC.L <- function(n.sim=10, n.iter = 200, ref.L = cells.L, beta=5.0, r0=0.15, return.sim = FALSE) {
  simulations <- matrix(0, nrow=n.sim, ncol=100)
  for (i in 1:n.sim) {
    print(paste("Simulation number", i, "of", n.sim))
    run <- sim.MCMC(42, n.iter, beta, r0)$state
    simulations[i, ] <- Kfn(pp=data.frame(x=run[,1], y=run[,2]), fs=0.5, k=100)$y
  }
  
  if(return.sim){
    return( simulations )
  }
  
  # add prediction interval to reference
  ref.L$lower = apply(simulations, 2, quantile, probs=0.05)
  ref.L$upper = apply(simulations, 2, quantile, probs=0.95)
  
  plot.L <- ggplot(ref.L, aes(x=r, y=L)) +
    geom_line(aes(y=lower), colour="red") + geom_line(aes(y=upper), colour="red") +
    geom_ribbon(aes(ymin=lower,ymax=upper), fill="red", alpha=0.05) +
    geom_line(colour="blue") +
    labs(x=expression(r), y=expression(L(r))) + theme_bw()
  plot.L
}


## Analyse with guessed parameters
cells.L <- with(Kfn(pp=cells, fs=0.5, k=100), data.frame("r"=x, "L"=y))

# Find the burn-in time
set.seed(4250)
run1 <- sim.MCMC(n.points=42, n.iter=500, beta=5.0, r0=0.15)
run2 <- sim.MCMC(n.points=42, n.iter=500, beta=5.0, r0=0.15)
run3 <- sim.MCMC(n.points=42, n.iter=500, beta=5.0, r0=0.15)
run4 <- sim.MCMC(n.points=42, n.iter=500, beta=5.0, r0=0.15)
run5 <- sim.MCMC(n.points=42, n.iter=500, beta=5.0, r0=0.15)

MCMC.energy <- ggplot(data.frame(run1 = run1$energies, run2 = run2$energies, run3 = run3$energies, run4 = run4$energies, run5 = run5$energies), aes(x=0:500)) +
  geom_line(aes(y=run1), colour="red") +
  geom_line(aes(y=run2), colour="magenta") +
  geom_line(aes(y=run3), colour="orange") +
  geom_line(aes(y=run4), colour="green") +
  geom_line(aes(y=run5), colour="blue") +
  labs(x="MCMC Iteration", y="Potential energy") + theme_bw()

ggsave("figures/MCMC_energies.pdf", plot = MCMC.energy, width = 8.0, height = 4.0)

plot(run3$energies, type='l')
# 100 iterations seems enough to ensure convergence

plot.L <- sim.MC.L(n.sim = 100, n.iter = 100, ref.L = cells.L, beta = 5.0, r0 = 0.15)
plot.L
ggsave("figures/cells_Stross_L.pdf", plot = plot.L, width = 4.0, height = 4.0)



## iterate to find theta


iterate.theta <- function(theta0, n.iterations, proposals, change, n.sim){
  # Empirical L
  cells.L <- with(Kfn(pp=cells, fs=0.5, k=100), data.frame("r"=x, "L"=y))
  
  # store iterations
  theta.vec <- matrix(NA, nrow = n.iterations+1, ncol = 2)
  theta.vec[1,] <- theta0
  
  # iterations
  for(i in 1:n.iterations){
    print(paste("Iteration number", i))
    
    # generate proposals for theta
    theta.proposals = matrix(rep(theta.vec[i,], proposals), nrow = proposals, byrow=TRUE) *
      runif(2*proposals, min=1-change, max=1+change)
    
    # find proposal with lowest objective
    objectives <- numeric(proposals)
    for(j in 1:proposals){
      print(paste("Proposal", j))
      # simulate for each proposal
      sim.L <- sim.MCMC.L(n.sim=n.sim, n.iter = 100, ref.L = cells.L, beta=theta.proposals[j,1], r0=theta.proposals[j,2], return.sim = TRUE)
      proposal.L <- colMeans( sim.L )
      objectives[j] <- sum((proposal.L[1:50]-cells.L$L[1:50])^2) + 0.1*sum((proposal.L[51:100]-cells.L$L[51:100])^2)
    }
    
    # find the lowest objective proposal
    theta.vec[i+1, ] <- theta.proposals[which.min(objectives), ]
  }
  return(theta.vec)
}

# estimate the difference between the two functions
theta0 = c(5.0, 0.15) # initial value

# iterating to find "better" values of theta 
theta.vec <- iterate.theta(theta0, n.iterations = 10, proposals = 5, change = 0.1, n.sim = 10)
theta0 <- theta.vec[11,]
theta0


new.beta <- theta0[1]
new.r0 <- theta0[2]

set.seed(4250)
run1 <- sim.MCMC(n.points=42, n.iter=500, beta=new.beta, r0=new.r0)
run2 <- sim.MCMC(n.points=42, n.iter=500, beta=new.beta, r0=new.r0)
run3 <- sim.MCMC(n.points=42, n.iter=500, beta=new.beta, r0=new.r0)
run4 <- sim.MCMC(n.points=42, n.iter=500, beta=new.beta, r0=new.r0)
run5 <- sim.MCMC(n.points=42, n.iter=500, beta=new.beta, r0=new.r0)

MCMC.energy2 <- ggplot(data.frame(run1 = run1$energies, run2 = run2$energies, run3 = run3$energies, run4 = run4$energies, run5 = run5$energies), aes(x=0:500)) +
  geom_line(aes(y=run1), colour="red") +
  geom_line(aes(y=run2), colour="magenta") +
  geom_line(aes(y=run3), colour="orange") +
  geom_line(aes(y=run4), colour="green") +
  geom_line(aes(y=run5), colour="blue") +
  labs(x="MCMC Iteration", y="Potential energy") + theme_bw()
MCMC.energy2
ggsave("figures/MCMC_energies2.pdf", plot = MCMC.energy, width = 8.0, height = 4.0)

#plot(run3$energies, type='l')
# 100 iterations seems enough to ensure convergence

plot.L <- sim.MC.L(n.sim = 100, n.iter = 100, ref.L = cells.L, beta = new.beta, r0 = 0.1)
plot.L
ggsave("figures/cells_Stross_L2.pdf", plot = plot.L, width = 4.0, height = 4.0)



## Three new realizations using the iterated values of theta

set.seed(4250)
sim1 <- sim.MCMC(n.points = 42, n.iter = 100, beta = new.beta, r0 = new.r0)$state
sim2 <- sim.MCMC(n.points = 42, n.iter = 100, beta = new.beta, r0 = new.r0)$state
sim3 <- sim.MCMC(n.points = 42, n.iter = 100, beta = new.beta, r0 = new.r0)$state


cells.sim1 <- ggplot(data.frame(x=sim1[,1], y=sim1[,2]), aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()
cells.sim1

cells.sim2 <- ggplot(data.frame(x=sim2[,1], y=sim2[,2]), aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()
cells.sim2

cells.sim3 <- ggplot(data.frame(x=sim3[,1], y=sim3[,2]), aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()
cells.sim3

ggsave("figures/cells_sim1.pdf", plot = cells.sim1, width = 3.0, height = 3.0)
ggsave("figures/cells_sim2.pdf", plot = cells.sim2, width = 3.0, height = 3.0)
ggsave("figures/cells_sim3.pdf", plot = cells.sim3, width = 3.0, height = 3.0)

