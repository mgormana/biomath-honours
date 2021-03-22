#odesolver2
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(matlib)


#new model
totalPop <- function(t, x, parms = NULL) {
  with(as.list(c(x, parms)), {
  dS = -h*(I+U)*S
  dL = h*(I+U)*S - (e*L)
  dI = -(e*I) + (e*L)
  dD = (a1*e)*I + b*d*U - g*d*D - (1-g)*d*D
  dU = (1-a1)*e*I - b*d*U - (1-b)*d*U
  dR = (1-g)*d*D + (1-b)*d*U
  dF = g*d*D
  list(c(dS, dL, dI, dD, dR, dU, dF))
 })
}

params1a = c(a1=0.02, b= 0.975, e=0.4, d = 0.12, g= 0.061, h = 0.2)
params1b = c(a1=0.5, b= 0.975, e=0.4, d = 0.12, g= 0.061, h = 0.2)
params1c = c(a1 = 0.9, b= 0.975, e = 0.4, d = 0.12, g = 0.061, h = 0.2)
tf <- 30
times <- seq(0,tf, by = 1)
t0.1 <- c(S = 1, L = 0, I = 10, D = 0, U = 10, R = 300, F = 5)


#test the 1 susceptible within first 2
susceptODE1a <- ode(y = t0.1, times, func = totalPop, parms = params1a)
susceptODE1a.df <- melt(as.data.frame(susceptODE1a), id = "time")
susceptODE1a.df <- rename(susceptODE1a.df, state = variable, scenarioA = value)

#test the 1 susceptible within first 5 days
susceptODE1b <- ode(y=t0.1, times, func = totalPop, parms = params1b)
susceptODE1b.df <- melt(as.data.frame(susceptODE1b), id = "time")
susceptODE1b.df <- rename(susceptODE1b.df, state = variable, scenarioB = value)

#test the 1 susceiptible within days 5-10
susceptODE1c <- ode(y=t0.1, times, func = totalPop, parms = params1c)
susceptODE1c.df <- melt(as.data.frame(susceptODE1c), id = "time")
susceptODE1c.df <- rename(susceptODE1c.df, state = variable, scenarioC = value)

#merge to one data frame
susceptODE1.df <- merge(susceptODE1a.df, susceptODE1b.df, by = c("time", "state"))
susceptODE1.df <- merge(susceptODE1.df, susceptODE1c.df, by = c("time", "state"))

#plotting
ggplot(susceptODE1.df, aes(x=time, y = scenarioA, color = state)) + 
  geom_line()

ggplot(susceptODE1.df, aes(x=time, y = scenarioB, color = state)) + 
  geom_line()

ggplot(susceptODE1.df, aes(x=time, y = scenarioC, color = state)) + 
  geom_line()

#t0.2 <- c(new initial conditions given by output of ODE)
susceptODE2a <- ode(y = t0.2, times, func = totalPop, parms = params2a)
susceptODE2a.df <- melt(as.data.frame(susceptODE2a), id = "time")
susceptODE2a.df <- rename(susceptODE2a.df, state = variable, scenarioA = value)


##Finding appropriate h value using the transition matrix V
a = 0.02
b = 0.975
d = 0.12
e = 0.4
g = 0.061

V <- matrix(c(
  -e,       0,  0,    0,
  e,      -e,  0,    0,
  0,     a*e, -d,  b*d,
  0, (1-a)*e,  0,   -d),
  nrow = 4, ncol = 4, byrow = TRUE)

# V <- matrix(c(
#   -e, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#   e, -e, 0, 0, 0, 0, 0, 0, 0, 0,
#   0, e, -e, 0, 0, 0, 0, 0, 0, 0,
#   0, 0, (a*e), -d, 0, 0, 0, 0, 0,0, 
#   0, 0, 0, d, -d, b*d, 0, 0, 0, 0,
#   0, 0, (1-a)*e, 0, 0, -d, 0, 0, 0,0, 
#   0, 0, 0, 0, 0, (1-b)*d, -d, 0, 0, 0,
#   0, 0, 0, 0, (1-g)*d, 0, 0, 0, 0, 0,
#   0, 0, 0, 0, 0, 0, d, 0, 0, 0,
#   0, 0, 0, 0, g*d, 0, 0, 0, 0,0), 
#   nrow = 10, ncol = 10, byrow = TRUE)

inverseV <- solve(-V)

h = 0.2
w = matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1, byrow = FALSE)
j = matrix(c(0, h, 0, h), nrow = 1, ncol = 4, byrow = TRUE)
S = 1

Ro = j %*% inverseV %*% w %*% S
