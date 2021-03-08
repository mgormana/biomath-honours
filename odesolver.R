library("deSolve")
library(reshape2)
library(ggplot2)
library(dplyr)
asymptomatics <- function(t, x, parms = NULL) {
  with(as.list(c(x, parms)), {
    
    dL = -(e*L)
    dI = -(e*I) + (e*L)
    dD = (a*e)*I + b*d*U - g*d*D - (1-g)*d*D
    dU = (1-a)*e*I - b*d*U - (1-b)*d*U
    dR = (1-g)*d*D + (1-b)*d*U
    dF = g*d*D
    list(c(dL, dI, dD, dU, dR, dF))
    
  })
}

#mean duration of infection
param1 <- c(a = 0.02, b = 0.975, e = 2/5, d = 2/(21-(1/.4)), g = 0.061)
param2 <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(21-(1/.4)), g = 0.061)
param3 <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(21-(1/.4)), g = 0.061)
tf <- 30
times <- seq(0,tf, by = 1)
t0 <- c(L = 1, I = 0, D = 0, U = 0, R = 0, F = 0)

#95% bounds of infectious period 
param1b <- c(a = 0.02, b = 0.975, e = 2/5, d = 2/(2.6-(1/.4)), g = 0.061)
param1c <- c(a = 0.02, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

param2b <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(2.6-(1/.4)), g = 0.061)
param2c <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

param3b <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(2.6-(1/.4)), g = 0.061)
param3c <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

res_ode1 <- ode(y = t0, times, func = asymptomatics, param1)
res_ode1.df <- melt(as.data.frame(res_ode1), id = "time")
res_ode1.df <- rename(res_ode1.df, state = variable, scenario1 = value)

res_ode1b <- ode(y = t0, times, func = asymptomatics, param1b)
res_ode1b.df <- melt(as.data.frame(res_ode1b), id = "time")
res_ode1b.df <- rename(res_ode1b.df, state = variable, scenario1b = value) 

res_ode1c <- ode(y = t0, times, func = asymptomatics, param1c)
res_ode1c.df <- melt(as.data.frame(res_ode1c), id = "time")
res_ode1c.df <- rename(res_ode1c.df, state = variable, scenario1c = value)

res_ode2 <- ode(y = t0, times, func = asymptomatics, param2)
res_ode2.df <- melt(as.data.frame(res_ode2), id = c("time"))
res_ode2.df <- rename(res_ode2.df, state = variable, scenario2 = value)

res_ode2b <- ode(y = t0, times, func = asymptomatics, param2b)
res_ode2b.df <- melt(as.data.frame(res_ode2b), id = c("time"))
res_ode2b.df <- rename(res_ode2b.df, state = variable, scenario2b = value)

res_ode2c <- ode(y = t0, times, func = asymptomatics, param2c)
res_ode2c.df <- melt(as.data.frame(res_ode2c), id = c("time"))
res_ode2c.df <- rename(res_ode2c.df, state = variable, scenario2c = value)

res_ode3 <- ode(y = t0, times, func = asymptomatics, param3)
res_ode3.df <- melt(as.data.frame(res_ode3), id= c("time"))
res_ode3.df <- rename(res_ode3.df, state = variable, scenario3 = value)

res_ode3b <- ode(y = t0, times, func = asymptomatics, param3b)
res_ode3b.df <- melt(as.data.frame(res_ode3b), id= c("time"))
res_ode3b.df <- rename(res_ode3b.df, state = variable, scenario3b = value)

res_ode3c <- ode(y = t0, times, func = asymptomatics, param3c)
res_ode3c.df <- melt(as.data.frame(res_ode3c), id= c("time"))
res_ode3c.df <- rename(res_ode3c.df, state = variable, scenario3c = value)

res_ode.df <- merge(res_ode1.df, res_ode2.df, by = c("state", "time"))
res_ode.df <- merge(res_ode.df, res_ode3.df, by = c("state", "time"))

res_odetest.df <- merge(res_ode.df, res_ode1b.df, by = c("state", "time"))
res_odetest.df <- merge(res_odetest.df, res_ode1c.df, by = c("state", "time"))
res_odetest.df <- merge(res_odetest.df, res_ode2b.df, by = c("state", "time"))
res_odetest.df <- merge(res_odetest.df, res_ode2c.df, by = c("state", "time"))
res_odetest.df <- merge(res_odetest.df, res_ode3b.df, by = c("state", "time"))
res_odetest.df <- merge(res_odetest.df, res_ode3c.df, by = c("state", "time"))

res_ode.df <- melt(res_ode.df, id = c("time", "state"))
res_ode.df <- rename(res_ode.df, scenario = variable, stateProb = value)

res_odetest.df <- melt(res_odetest.df, id = c("time", "state"))
res_odetest.df <- rename(res_odetest.df, scenario = variable, stateProb = value)


res_ode1test.df <- merge(res_ode1.df, res_ode1b.df, by = c("state", "time"))
res_ode1test.df <- merge(res_ode1test.df, res_ode1c.df, by = c("state", "time"))
res_ode1test.df <- melt(res_ode1test.df, id = c("time", "state"))
res_ode1test.df <- rename(res_ode1test.df, scenario = variable, stateProb = value)

res_ode2test.df <- merge(res_ode2.df, res_ode2b.df, by = c("state", "time"))
res_ode2test.df <- merge(res_ode2test.df, res_ode2c.df, by = c("state", "time"))
res_ode2test.df <- melt(res_ode2test.df, id = c("time", "state"))
res_ode2test.df <- rename(res_ode2test.df, scenario = variable, stateProb = value)

res_ode3test.df <- merge(res_ode3.df, res_ode3b.df, by = c("state", "time"))
res_ode3test.df <- merge(res_ode3test.df, res_ode3c.df, by = c("state", "time"))
res_ode3test.df <- melt(res_ode3test.df, id = c("time", "state"))
res_ode3test.df <- rename(res_ode3test.df, scenario = variable, stateProb = value)

#plotting 
ggplot(res_ode.df, aes(x=time, y = stateProb, color = state)) +
  geom_line() + 
  facet_wrap(~scenario)

ggplot(res_odetest.df, aes(x = time, y = stateProb, color = state)) + 
  geom_line() + 
  facet_wrap(~scenario)

ggplot(res_ode1test.df, aes(x=time, y = stateProb, color = state)) + 
  geom_line()  + 
  facet_wrap(~scenario)

ggplot(res_ode2test.df, aes(x=time, y = stateProb, color = state)) + 
  geom_line() + 
  facet_wrap(~scenario)

ggplot(res_ode3test.df, aes(x = time, y = stateProb, color = state)) +
  geom_line() + 
  facet_wrap(~scenario)
  

#matrix exp function
probDetectedFn <- matrix(c(
  exp(-(2*t)/5), 0, 0, 0, 0, 0, 0, 0, 0, 
  (2*t*exp(-(2*t)/5))/5,exp(-(2*t)/5), 0,0,0,0, 0, 0, 0, 
  (13689*exp(-(4*t)/39))/8410 - (13689*exp(-(2*t)/5))/8410 - (351*t*exp(-(2*t)/5))/725,(351*exp(-(4*t)/39))/290 - (351*exp(-(2*t)/5))/290,exp(-(4*t)/39),0,0, 0, 0, 0, 0, 
  (606879*exp(-(2*t)/5))/487780 - (606879*exp(-(4*t)/39))/487780 + (15561*t*exp(-(2*t)/5))/84100 + (15561*t*exp(-(4*t)/39))/84100,(15561*exp(-(2*t)/5))/33640 - (15561*exp(-(4*t)/39))/33640 + (399*t*exp(-(4*t)/39))/2900, (4*t*exp(-(4*t)/39))/39, exp(-(4*t)/39), (t*exp(-(4*t)/39))/10, 0, 0, 0, 0,
  (1521*exp(-(4*t)/39))/8410 - (1521*exp(-(2*t)/5))/8410 - (39*t*exp(-(2*t)/5))/725, (39*exp(-(4*t)/39))/290 - (39*exp(-(2*t)/5))/290, 0,0, exp(-(4*t)/39),0, 0, 0, 0,
  (1521*exp(-(2*t)/5))/487780 - (1521*exp(-(4*t)/39))/487780 + (39*t*exp(-(2*t)/5))/84100 + (39*t*exp(-(4*t)/39))/84100, (39*exp(-(2*t)/5))/33640 - (39*exp(-(4*t)/39))/33640 + (t*exp(-(4*t)/39))/2900, 0, 0,(t*exp(-(4*t)/39))/390,exp(-(4*t)/39), 0, 0, 0,
  374661/400000 - (5128734429*exp(-(4*t)/39))/9755600000 - (374661*t*exp(-(2*t)/5))/8410000 - (14611779*t*exp(-(4*t)/39))/84100000 - (40088727*exp(-(2*t)/5))/97556000, 374661/400000 - (277623801*exp(-(4*t)/39))/336400000 - (374661*t*exp(-(4*t)/39))/2900000 - (374661*exp(-(2*t)/5))/3364000, 939/1000 - (313*t*exp(-(4*t)/39))/3250 - (939*exp(-(4*t)/39))/1000, 939/1000 - (939*exp(-(4*t)/39))/1000, 36621/40000 - (939*t*exp(-(4*t)/39))/10000 - (36621*exp(-(4*t)/39))/40000, 0, 1, 0, 0,
  1/400 - (13689*exp(-(4*t)/39))/9755600 - (t*exp(-(2*t)/5))/8410 - (39*t*exp(-(4*t)/39))/84100 - (107*exp(-(2*t)/5))/97556,1/400 - (741*exp(-(4*t)/39))/336400 - (t*exp(-(4*t)/39))/2900 - exp(-(2*t)/5)/3364, 0, 0,1/40 - (t*exp(-(4*t)/39))/390 - exp(-(4*t)/39)/40, 1 - exp(-(4*t)/39), 0, 1, 0,
  24339/400000 - (333176571*exp(-(4*t)/39))/9755600000 - (24339*t*exp(-(2*t)/5))/8410000 - (949221*t*exp(-(4*t)/39))/84100000 - (2604273*exp(-(2*t)/5))/97556000,24339/400000 - (18035199*exp(-(4*t)/39))/336400000 - (24339*t*exp(-(4*t)/39))/2900000 - (24339*exp(-(2*t)/5))/3364000, 61/1000 - (61*t*exp(-(4*t)/39))/9750 - (61*exp(-(4*t)/39))/1000,   61/1000 - (61*exp(-(4*t)/39))/1000, 2379/40000 - (61*t*exp(-(4*t)/39))/10000 - (2379*exp(-(4*t)/39))/40000, 0, 0, 0, 1
  ), nrow = 9, ncol = 9, byrow = TRUE)




