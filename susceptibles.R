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
  list(c(dS, dL, dI, dD, dU, dR, dF))
 })
}

params1a = c(a1=0.2, b= 0.975, e=0.4, d = 0.12, g= 0.061, h = 0.2)
params1b = c(a1=0.15, b= 0.975, e=0.4, d = 0.12, g= 0.061, h = 0.2)
params1c = c(a1 = 0.1, b= 0.975, e = 0.4, d = 0.12, g = 0.061, h = 0.2)
tf <- 30
times <- seq(0,tf, by = 1)
t0.1 <- c(S = .9, L = 0, I = 0.1, D = 0, U = 0, R = 0, F = 0)


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
susceptODE1.df <- melt(susceptODE1.df, id = c("time", "state"))
susceptODE1.df <- rename(susceptODE1.df, scenario = variable, stateProb = value)

#plotting
ggplot(susceptODE1.df, aes(x=time, y = stateProb, color = state)) + 
  geom_line() + 
  facet_wrap(~scenario)

tail(susceptODE1a)
tail(susceptODE1b)
tail(susceptODE1c)
t1.a <- c(S = 0.42, L = 0.03, I = 0.04, D = 0.12, U = 0.10, R = 0.27, F = 0.02)
t1.b <- c(S = 0.39, L = 0.03, I = 0.04, D = 0.12, U = 0.12, R = 0.28, F = 0.17)
t1.c <- c(S = 0.37, L = 0.04, I = 0.04, D = 0.13, U = 0.13, R = 0.28, F = 0.02)
paramsX = c(a1 = 0.02, b = 0.975, d = 0.12, e = 0.4, g = 0.061, h = 0.2)
paramsY = c(a1 = 0.5, b = 0.975, d = 0.12, e = 0.4, g = 0.061, h = 0.2)
paramsZ = c(a1 = 0.9, b = 0.975, d = 0.12, e = 0.4, g = 0.061, h = 0.2)


susceptODE2ax <- ode(y = t1.a, times, func = totalPop, parms = paramsX)
susceptODE2ax.df <- melt(as.data.frame(susceptODE2ax), id = "time")
susceptODE2ax.df <- rename(susceptODE2ax.df, state = variable, scenarioX = value)

susceptODE2ay <- ode(y=t1.a, times, func = totalPop, parms = paramsY)
susceptODE2ay.df <- melt(as.data.frame(susceptODE2ay), id = "time")
susceptODE2ay.df <- rename(susceptODE2ay.df, state = variable, scenarioY = value)

susceptODE2az <- ode(y=t1.a, times, func = totalPop, parms = paramsZ)
susceptODE2az.df <- melt(as.data.frame(susceptODE2az), id = "time")
susceptODE2az.df <- rename(susceptODE2az.df, state = variable, scenarioZ = value)

susceptODE2a.df <- merge(susceptODE2ax.df, susceptODE2ay.df, by = c("time", "state"))
susceptODE2a.df <- merge(susceptODE2a.df, susceptODE2az.df, by = c("time", "state"))
susceptODE2a.df <- melt(susceptODE2a.df, id = c("time", "state"))
susceptODE2a.df <- rename(susceptODE2a.df, stateProb = value, scenarioA = variable)

ggplot(susceptODE2a.df, aes(x=time, y = stateProb, color = state)) + 
  geom_line() +
  facet_wrap(~scenarioA)

#NOW FOR OTHER SCNEARIO
susceptODE2bx <- ode(y = t1.b, times, func = totalPop, parms = paramsX)
susceptODE2bx.df <- melt(as.data.frame(susceptODE2bx), id = "time")
susceptODE2bx.df <- rename(susceptODE2bx.df, state = variable, scenarioX = value)

susceptODE2by <- ode(y=t1.b, times, func = totalPop, parms = paramsY)
susceptODE2by.df <- melt(as.data.frame(susceptODE2by), id = "time")
susceptODE2by.df <- rename(susceptODE2by.df, state = variable, scenarioY = value)

susceptODE2bz <- ode(y=t1.b, times, func = totalPop, parms = paramsZ)
susceptODE2bz.df <- melt(as.data.frame(susceptODE2bz), id = "time")
susceptODE2bz.df <- rename(susceptODE2bz.df, state = variable, scenarioZ = value)

susceptODE2b.df <- merge(susceptODE2bx.df, susceptODE2by.df, by = c("time", "state"))
susceptODE2b.df <- merge(susceptODE2b.df, susceptODE2bz.df, by = c("time", "state"))
susceptODE2b.df <- melt(susceptODE2b.df, id = c("time", "state"))
susceptODE2b.df <- rename(susceptODE2b.df, scenarioB = variable, stateProb = value)

ggplot(susceptODE2b.df, aes(x=time, y = stateProb, color = state)) + 
  geom_line() +
  facet_wrap(~scenarioB)

##
susceptODE2cx <- ode(y = t1.c, times, func = totalPop, parms = paramsX)
susceptODE2cx.df <- melt(as.data.frame(susceptODE2cx), id = "time")
susceptODE2cx.df <- rename(susceptODE2cx.df, state = variable, scenarioX = value)

susceptODE2cy <- ode(y=t1.c, times, func = totalPop, parms = paramsY)
susceptODE2cy.df <- melt(as.data.frame(susceptODE2cy), id = "time")
susceptODE2cy.df <- rename(susceptODE2cy.df, state = variable, scenarioY = value)

susceptODE2cz <- ode(y=t1.c, times, func = totalPop, parms = paramsZ)
susceptODE2cz.df <- melt(as.data.frame(susceptODE2cz), id = "time")
susceptODE2cz.df <- rename(susceptODE2cz.df, state = variable, scenarioZ = value)

susceptODE2c.df <- merge(susceptODE2cx.df, susceptODE2cy.df, by = c("time", "state"))
susceptODE2c.df <- merge(susceptODE2c.df, susceptODE2cz.df, by = c("time", "state"))
susceptODE2c.df <- melt(susceptODE2c.df, id = c("time", "state"))
susceptODE2c.df <- rename(susceptODE2c.df, scenarioC = variable, stateProb = value)

ggplot(susceptODE2c.df, aes(x=time, y = stateProb, color = state)) + 
  geom_line() +
  facet_wrap(~scenarioC)


##Finding appropriate h value using the transition matrix V
a = 0.9
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


inverseV <- solve(-V)

h = 0.2
w = matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1, byrow = FALSE)
j = matrix(c(0, h, 0, h), nrow = 1, ncol = 4, byrow = TRUE)
S = 1

Ro = j %*% inverseV %*% w %*% S
