# Project 1 in TMA4250 Spatial Statistics

library(MASS)
library(fields)
library(akima)
library(geoR)
library(gridExtra)
library(latex2exp)
library(tidyverse)
library(reshape2)

# seed
set.seed(4250)


### Problem 1


## Defining the Gaussian random field
D <- 1:50 # Discretization of D
D.dist <- D-1 # Distances in D
mu <- rep(0, 50) # Expected value
sigma2 <- c(1, 5) # Variance

## Powered exponential correlation
powexp.a <- 10 # Range parameter
powexp.alpha <- c(1, 1.9) # Power parameter

## Matérn correlation
matern.a  <- 20 # Range parameter
matern.nu <- c(1,3) # Smoothness


## a)

# Plotting the four correlations functions (use cov.spatial with sigma2 = 1 to get correlation)
corr.pe1 <- cov.spatial(D.dist, cov.model="powered.exponential", cov.pars = c(1, powexp.a), kappa=powexp.alpha[1])
corr.pe2 <- cov.spatial(D.dist, cov.model="powered.exponential", cov.pars = c(1, powexp.a), kappa=powexp.alpha[2])
#corr3 <- cov.spatial(D.dist, cov.model="matern", cov.pars = c(1, matern.a), kappa=matern.nu[1]) # does not work
#corr4 <- cov.spatial(D.dist, cov.model="matern", cov.pars = c(1, matern.a), kappa=matern.nu[2]) # does not work

matern.corr <- function(h, a, nu){
  arg <- sqrt(8*nu)*h / a
  return.arr <- 2^(1-nu) / gamma(nu) * arg^nu * besselK(arg, nu)
  return.arr[h==0] = 1 # fix for h=0
  return.arr
}
corr.m1 <- matern.corr(D.dist, matern.a, matern.nu[1])
corr.m2 <- matern.corr(D.dist, matern.a, matern.nu[2])


## Plot correlations
df <- data.frame("h"=D.dist, "pe1"=corr.pe1, "pe2"=corr.pe2, "m1"=corr.m1, "m2"=corr.m2)
df <- melt(df, id.vars='h', variable.name = 'Correlations')

plot <- ggplot(df, aes(x = h, y = value)) +
  geom_line(aes(colour = Correlations)) +
  xlab(TeX("Distance $h$")) +
  ylab(TeX("Correlation $\\rho(h)$")) +
  theme_bw()
plot
ggsave("figures/correlation.pdf", plot = plot, height = 3.00, width = 6.00)


## Plot variograms
variogram <- function(corr, sigma2) {
  sigma2*(1 - corr)
}

vari.pe1.s1 <- variogram(corr.pe1, sigma2[1])
vari.pe2.s1 <- variogram(corr.pe2, sigma2[1])
vari.m1.s1 <- variogram(corr.m1, sigma2[1])
vari.m2.s1 <- variogram(corr.m2, sigma2[1])

vari.pe1.s2 <- variogram(corr.pe1, sigma2[2])
vari.pe2.s2 <- variogram(corr.pe2, sigma2[2])
vari.m1.s2 <- variogram(corr.m1, sigma2[2])
vari.m2.s2 <- variogram(corr.m2, sigma2[2])


df <- data.frame("h"=D.dist, "pe1.s1"=vari.pe1.s1, "pe2.s1"=vari.pe2.s1, "m1.s1"=vari.m1.s1, "m2.s1"=vari.m2.s1)
df <- melt(df, id.vars='h', variable.name = 'Variograms')

plot <- ggplot(df, aes(x = h, y = value)) +
  geom_line(aes(colour = Variograms)) +
  xlab(TeX("Distance $h$")) +
  ylab(TeX("Variogram $2\\gamma(h)$")) +
  ylim(0, 5) +
  theme_bw()
plot
ggsave("figures/variogram1.pdf", plot = plot, height = 2.00, width = 4.00)


df <- data.frame("h"=D.dist, "pe1.s2"=vari.pe1.s2, "pe2.s2"=vari.pe2.s2, "m1.s2"=vari.m1.s2, "m2.s2"=vari.m2.s2)
df <- melt(df, id.vars='h', variable.name = 'Variograms')

plot <- ggplot(df, aes(x = h, y = value)) +
  geom_line(aes(colour = Variograms)) +
  xlab(TeX("Distance $h$")) +
  ylab(TeX("Variogram $2\\gamma(h)$")) +
  ylim(0, 5) +
  theme_bw()
plot
ggsave("figures/variogram2.pdf", plot = plot, height = 2.00, width = 4.00)




## b)

# Simulating four realizations four each combination

simulate.X <- function(corr, sigma2, mu=NA, nsim=4){
  # Covarinace matrix
  Sigma <- sigma2*toeplitz(corr)
  
  # Set mu = 0 if not specified
  if(is.na(mu)) mu <- rep(0, length(corr))
  
  # Simulate from multivariate normal distribution
  simulations <- mvrnorm(n = nsim, mu=mu, Sigma=Sigma)
}

plot.X <- function(sim){
  df <- data.frame(D=D, x1 = sim[1,], x2 = sim[2,], x3 = sim[3,], x4 = sim[4,])
  plot <- ggplot(data=df, aes(x=D)) +
    geom_line(aes(y=x1, col="1")) +
    geom_line(aes(y=x2, col="2")) + 
    geom_line(aes(y=x3, col="3")) +
    geom_line(aes(y=x4, col="4")) +
    xlab("s")+
    ylab("X(s)") +
    theme_bw()+
    theme(legend.position = "none")
  plot
}


## Simulate the 8 different combinations

# powexp(10, 1), sigma2 = 1
sim <- simulate.X(corr1, sigma2[1])
plot <- plot.X(sim)
ggsave("figures/realisations1.pdf", plot = plot, height = 2.00, width = 4.00)

# powexp(10, 1.9), sigma2 = 1
sim <- simulate.X(corr2, sigma2[1])
plot <- plot.X(sim)
ggsave("figures/realisations2.pdf", plot = plot, height = 2.00, width = 4.00)

# matern(20, 1), sigma2 = 1
sim <- simulate.X(corr3, sigma2[1])
plot <- plot.X(sim)
ggsave("figures/realisations3.pdf", plot = plot, height = 2.00, width = 4.00)

# matern(20, 3), sigma2 = 1
sim <- simulate.X(corr4, sigma2[1])
plot <- plot.X(sim)
ggsave("figures/realisations4.pdf", plot = plot, height = 2.00, width = 4.00)


# powexp(10, 1), sigma2 = 5
sim <- simulate.X(corr1, sigma2[2])
plot <- plot.X(sim)
ggsave("figures/realisations5.pdf", plot = plot, height = 2.00, width = 4.00)

# powexp(10, 1.9), sigma2 = 5
sim <- simulate.X(corr2, sigma2[2])
plot <- plot.X(sim)
ggsave("figures/realisations6.pdf", plot = plot, height = 2.00, width = 4.00)

# matern(20, 1), sigma2 = 5
sim <- simulate.X(corr3, sigma2[2])
plot <- plot.X(sim)
ggsave("figures/realisations7.pdf", plot = plot, height = 2.00, width = 4.00)

# matern(20, 3), sigma2 = 5
sim <- simulate.X(corr4, sigma2[2])
plot <- plot.X(sim)
ggsave("figures/realisations8.pdf", plot = plot, height = 2.00, width = 4.00)




## c) and d)

# Select X as the last simulation; Matern(20, 3), sigma2=5
X <- sim[3,]


# Observation locations
s <- c(10, 25, 30)
m <- length(s)

# Create the measurement matrix
M <- matrix(0, nrow = m, ncol = length(D))
for (i in 1:m){
  M[i, s[i]] = 1
}

# Measurements without measurement noise
y <- M %*% X



## Create prediction intervals
R <- toeplitz(corr4)
Sigma <- sigma2[2] * R
# Upper 95% quantile for the normal distribution
Zq <- qnorm(0.95)


## Prediction plots with no measurement error
sigma2.N <- 0

# Posterior mean and variance
mu.posterior <- Sigma %*% t(M) %*% solve(M %*% Sigma %*% t(M) + sigma2.N*diag(3)) %*% y
Sigma.posterior <- Sigma - Sigma %*% t(M) %*% solve(M %*% Sigma %*% t(M) + sigma2.N*diag(3)) %*% M %*% Sigma

# Lower and upper bounds for 90% PI
low.posterior <- mu.posterior - Zq * sqrt(diag(Sigma.posterior))
high.posterior <- mu.posterior + Zq * sqrt(diag(Sigma.posterior))


df <- data.frame("s"=D, "X.mean"=mu.posterior, "X.low"=low.posterior, "X.high"=high.posterior)

plot <- ggplot(data=df, aes(x = s)) +
  scale_color_manual(values=c("black", "red")) +
  geom_line(aes(y=X.mean, col="1")) +
  geom_line(aes(y=X.low, col="2")) + 
  geom_line(aes(y=X.high, col="2")) +
  xlab(TeX("Location $s$")) +
  ylab(TeX("Prediction of $X$ given $y$")) +
  theme_bw() +
  theme(legend.position = "none")
plot
ggsave("figures/prediction1.pdf", plot = plot, height = 2.00, width = 4.00)


## 100 realizations from the posterior distribution
sim.posterior1 <- mvrnorm(n=100, mu.posterior, Sigma.posterior)


empirical.mean <- apply(sim.posterior1, MARGIN = 2, mean)
empirical.var <- apply(sim.posterior1, MARGIN = 2, var)

# Lower and upper bounds for 90% PI
low.empirical <- empirical.mean - Zq * sqrt(empirical.var)
high.empirical <- empirical.mean + Zq * sqrt(empirical.var)


df <- as.data.frame(t(sim.posterior))
df[,"s"] = D
df[,"mean"] = empirical.mean
df[,"low"] = low.empirical
df[,"high"] = high.empirical


plot <- ggplot(data = df, aes(x = s))
for (i in 1:100) {
  plot <- plot + geom_line(aes_string(y = paste0("V", i)), col="pink")
}

plot <- plot +
  scale_color_manual(values=c("black", "red")) +
  geom_line(aes(y=mean, col="1")) +
  geom_line(aes(y=low, col="2")) + 
  geom_line(aes(y=high, col="2")) +
  xlab(TeX("Location $s$")) +
  ylab(TeX("Estimation of $X$ given $y$")) +
  theme_bw() +
  theme(legend.position = "none")
plot
ggsave("figures/100realizations1.pdf", plot = plot, height = 2.00, width = 4.00)



## Prediction with measurement error
sigma2.N <- 0.25

# Posterior mean and variance
mu.posterior <- Sigma %*% t(M) %*% solve(M %*% Sigma %*% t(M) + sigma2.N*diag(3)) %*% y
Sigma.posterior <- Sigma - Sigma %*% t(M) %*% solve(M %*% Sigma %*% t(M) + sigma2.N*diag(3)) %*% M %*% Sigma

# Lower and upper bounds for 90% CI
low.posterior <- mu.posterior - Zq * sqrt(diag(Sigma.posterior))
high.posterior <- mu.posterior + Zq * sqrt(diag(Sigma.posterior))


df <- data.frame("s"=D, "X.mean"=mu.posterior, "X.low"=low.posterior, "X.high"=high.posterior)

plot <- ggplot(data=df, aes(x = s)) +
  scale_color_manual(values=c("black", "red")) +
  geom_line(aes(y=X.mean, col="1")) +
  geom_line(aes(y=X.low, col="2")) + 
  geom_line(aes(y=X.high, col="2")) +
  xlab(TeX("Location $s$")) +
  ylab(TeX("Prediction of $X$ given $y$")) +
  theme_bw() +
  theme(legend.position = "none")
plot
ggsave("figures/prediction2.pdf", plot = plot, height = 2.00, width = 4.00)


## 100 realizations from the posterior distribution
sim.posterior2 <- mvrnorm(n=100, mu.posterior, Sigma.posterior)

empirical.mean <- apply(sim.posterior2, MARGIN = 2, mean)
empirical.var <- apply(sim.posterior2, MARGIN = 2, var)

# Lower and upper bounds for 90% PI
low.empirical <- empirical.mean - Zq * sqrt(empirical.var)
high.empirical <- empirical.mean + Zq * sqrt(empirical.var)


df <- as.data.frame(t(sim.posterior))
df[,"s"] = D
df[,"mean"] = empirical.mean
df[,"low"] = low.empirical
df[,"high"] = high.empirical


plot <- ggplot(data = df, aes(x = s))
for (i in 1:100) {
  plot <- plot + geom_line(aes_string(y = paste0("V", i)), col="pink")
}

plot <- plot +
  scale_color_manual(values=c("black", "red")) +
  geom_line(aes(y=mean, col="1")) +
  geom_line(aes(y=low, col="2")) + 
  geom_line(aes(y=high, col="2")) +
  xlab(TeX("Location $s$")) +
  ylab(TeX("Estimation of $X$ given $y$")) +
  theme_bw() +
  theme(legend.position = "none")
plot
ggsave("figures/100realizations2.pdf", plot = plot, height = 2.00, width = 4.00)



## f)
# Calculate the area A using the 100 realizations from e)

sim.A <- apply( (sim.posterior1 > 2)*(sim.posterior1 - 2), MARGIN = 1, sum)

empirical.A.mean <- mean(sim.A)
empirical.A.sd <- sd(sim.A)
empirical.A.sd

# Alternative predictor.
alternative.A <- sum((mu.posterior>2)*(mu.posterior - 2))
alternative.A








