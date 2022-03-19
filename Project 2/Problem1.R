# Project 2 in TMA4250 Spatial Statistics

library(MASS)
library(spatial)
library(ggplot2)
library(tidyverse)
library(dplyr)


### (a)

## Read data
cells <- read.table("data/cells.dat", col.names = c("x", "y"))
redwood <- read.table("data/redwood.dat", col.names = c("x", "y"))
pines <- read.table("data/pines.dat", col.names = c("x", "y"), skip=3)


## Point plots
cells.plot <- ggplot(cells, aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()

redwood.plot <- ggplot(redwood, aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()

pines.plot <- ggplot(pines, aes(x=x, y=y)) +
  geom_point(size=0.5) + coord_fixed(xlim=c(0, 1), ylim=c(0, 1)) +
  labs(x=expression(x), y=expression(y)) + theme_bw()

ggsave("figures/cells.pdf", plot = cells.plot, width = 3.0, height = 3.0)
ggsave("figures/redwood.pdf", plot = redwood.plot, width = 3.0, height = 3.0)
ggsave("figures/pines.pdf", plot = pines.plot, width = 3.0, height = 3.0)



### (b)

# Set the domain for spatial point pattern analyses
ppregion(xl=0, xu=1, yl=0, yu=1)

# ?Kfn
# pp: a list such as a pp object, including components x and y
# fs: full-scale of plot
# k:  number of regularly spaced distances

cells.L <- with(Kfn(pp=cells, fs=0.5, k=100), data.frame("r"=x, "L"=y))
redwood.L <- with(Kfn(pp=redwood, fs=0.5, k=100), data.frame("r"=x, "L"=y))
pines.L <- with(Kfn(pp=pines, fs=0.5, k=100), data.frame("r"=x, "L"=y))

cells.L.plot <- ggplot(cells.L, aes(x=r, y=L)) +
  geom_line(colour="blue") + geom_line(aes(y=r), linetype="dashed") +
  labs(x=expression(r), y=expression(L(r))) + theme_bw()

redwood.L.plot <- ggplot(redwood.L, aes(x=r, y=L)) +
  geom_line(colour="blue") + geom_line(aes(y=r), linetype="dashed") +
  labs(x=expression(r), y=expression(L(r))) + theme_bw()

pines.L.plot <- ggplot(pines.L, aes(x=r, y=L)) +
  geom_line(colour="blue") + geom_line(aes(y=r), linetype="dashed") +
  labs(x=expression(r), y=expression(L(r))) + theme_bw()

ggsave("figures/cells_L.pdf", plot = cells.L.plot, width = 3.0, height = 3.0)
ggsave("figures/redwood_L.pdf", plot = redwood.L.plot, width = 3.0, height = 3.0)
ggsave("figures/pines_L.pdf", plot = pines.L.plot, width = 3.0, height = 3.0)



### c)

# function to simulate L(r) n times
sim.L <- function(n, data, fs=0.5, k=100) {
  count <- nrow(data)
  simulations <- matrix(0, nrow=n, ncol=k)
  for (i in 1:n) {
    simulations[i, ] <- Kfn(pp=Psim(count), fs, k)$y
  }
  simulations
}

set.seed(4250)

cells.sim.L = sim.L(100, cells)
cells.L$lower = apply(cells.sim.L, 2, quantile, probs=0.05)
cells.L$upper = apply(cells.sim.L, 2, quantile, probs=0.95)
cells.L.plot2 <- ggplot(cells.L, aes(x=r, y=L)) +
  geom_line(aes(y=lower), colour="red") + geom_line(aes(y=upper), colour="red") +
  geom_ribbon(aes(ymin=lower,ymax=upper), fill="red", alpha=0.05) +
  geom_line(colour="blue") +
  labs(x=expression(r), y=expression(L(r))) + theme_bw()

redwood.sim.L = sim.L(100, redwood)
redwood.L$lower = apply(redwood.sim.L, 2, quantile, probs=0.05)
redwood.L$upper = apply(redwood.sim.L, 2, quantile, probs=0.95)
redwood.L.plot2 <- ggplot(redwood.L, aes(x=r, y=L)) +
  geom_line(aes(y=lower), colour="red") + geom_line(aes(y=upper), colour="red") +
  geom_ribbon(aes(ymin=lower,ymax=upper), fill="red", alpha=0.05) +
  geom_line(colour="blue") +
  labs(x=expression(r), y=expression(L(r))) + theme_bw()

pines.sim.L = sim.L(100, pines)
pines.L$lower = apply(pines.sim.L, 2, quantile, probs=0.05)
pines.L$upper = apply(pines.sim.L, 2, quantile, probs=0.95)
pines.L.plot2 <- ggplot(pines.L, aes(x=r, y=L)) +
  geom_line(aes(y=lower), colour="red") + geom_line(aes(y=upper), colour="red") +
  geom_ribbon(aes(ymin=lower,ymax=upper), fill="red", alpha=0.05) +
  geom_line(colour="blue") +
  labs(x=expression(r), y=expression(L(r))) + theme_bw()

ggsave("figures/cells_L2.pdf", plot = cells.L.plot2, width = 3.0, height = 3.0)
ggsave("figures/redwood_L2.pdf", plot = redwood.L.plot2, width = 3.0, height = 3.0)
ggsave("figures/pines_L2.pdf", plot = pines.L.plot2, width = 3.0, height = 3.0)


