#odesolver2
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(matlib)

a1 = 0.02
a2 = 0.5
a3 = 0.9
b = 0.975
e = 0.4
d = 2/18.5
g = 0.061

h = 0.1
S = 100

param1 <- c(a = 0.02, b = 0.975, e = 0.4, d = 2/18.5, g = 0.061)
tf <- 30
times <- seq(0, tf, by = 1)
t0 <- c(S = 100, L = 0, I = 0, D = 0, U = 0, R = 0, F = 0)

totalPop <- function(t, x, parms = NULL) {
  dS = -h*(I+U)*S
  dL = h*(I+U)*S - (e*L)
  dI = -(e*I) + (e*L)
  dD = (a*e)*I + b*d*U - g*d*D - (1-g)*d*D
  dU = (1-a)*e*I - b*d*U - (1-b)*d*U
  dR = (1-g)*d*D + (1-b)*d*U
  dF = g*d*D
  list(c(dS, dL, dI, dD, dR, dF))
}

#transition matrix
V <-matrix(c(
  -e, 0, 0, 0,
  e, -e, 0, 0,
  0, a1*e, -d, b*d,
  0, (1-a1)*e, 0, -d),
  nrow = 4, ncol = 4, byrow = TRUE)

w = matrix(c(1, 0, 0, 0), ncol =1, byrow = FALSE)
beta = matrix(c(0, h, 0, h), nrow =1, byrow = TRUE)
Z = (w %*% beta)*S
x = matrix(c(L, I, D, U), ncol = 1, byrow = FALSE)

susceptibles <- function(t, parms = NULL) {
  with(as.list(c(x, parms)))
  dS = -h*(I+U)*S
  dx = Z*x + V*x
  list(c(dS, dx))
}

#says 'unused argument (parms)' ?
susceptODE <- ode(y = t0, times, func = susceptibles, parms = param1)
susceptODE.df <- as.data.frame(susceptODE)

#says V is singular
inverseV <- solve(V)
Ro = beta %*% inverseV %*% w %*% S
