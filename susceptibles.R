#odesolver2
library(deSolve)
library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)

#new model
totalPop <- function(t, x, parms = NULL) {
  with(as.list(c(x, parms)), {
  dS = -h*(I+U1+U2)*S
  dL = h*(I+U1+U2)*S - (e*L)
  dI = -(e*I) + (e*L)
  dD = (a1*e)*I + b*d*U1 - g*d*D - (1-g)*d*D
  dU1 = (1-a1)*e*I - d*U1
  dU2 = (1-b)*U1 - d*U2
  dRd = (1-g)*d*D  
  dRu = d*U2
  dF = g*d*D
  list(c(dS, dL, dI, dD, dU1, dU2, dRd, dRu, dF))
 })
}

a1 = 0.25
a2 = 0.15
a3 = 0.1
b = 0.975
d = 0.12
e = 0.4
g = 0.061

V1 <- matrix(c(
  -e,    0,      0,    0,
   e,   -e,      0,    0,
   0,   a1*e,   -d,   b*d,
   0, (1-a1)*e,  0,   -d),
  nrow = 4, ncol = 4, byrow = TRUE)

V2 <- matrix(c(
  -e,       0,  0,    0,
  e,      -e,  0,    0,
  0,     a2*e, -d,  b*d,
  0, (1-a2)*e,  0,   -d),
  nrow = 4, ncol = 4, byrow = TRUE)

V3 <- matrix(c(
  -e,       0,  0,    0,
  e,      -e,  0,    0,
  0,     a3*e, -d,  b*d,
  0, (1-a3)*e,  0,   -d),
  nrow = 4, ncol = 4, byrow = TRUE)

inverseV1 <- solve(-V1)
inverseV2 <- solve(-V2)
inverseV3 <- solve(-V3)

h1 = 0.23 #for testing within 2 days
h2 = 0.21 #for testing within 5 days
h3 = 0.2 #for testing after 5 days
w = matrix(c(1, 0, 0, 0), nrow = 4, ncol = 1, byrow = FALSE)
j1 = matrix(c(0, h1, 0, h1), nrow = 1, ncol = 4, byrow = TRUE)
j2 = matrix(c(0, h2, 0, h2), nrow = 1, ncol = 4, byrow = TRUE)
j3 = matrix(c(0, h3, 0, h3), nrow = 1, ncol = 4, byrow = TRUE)
S = 1

Ro1 = j1 %*% inverseV1 %*% w %*% S
Ro2 = j2 %*% inverseV2 %*% w %*% S
Ro3 = j3 %*% inverseV3 %*% w %*% S
rm(h1, h2, h3, w, j1, j2, j3, S, inverseV1, inverseV2, inverseV3, V1, V2, V3, a1, a2, a3, b, e, d, g, Ro1, Ro2, Ro3)


params1a = c(a1=0.25, b= 0.975, e=0.4, d = 0.12, g= 0.061, h = 0.23)
params1b = c(a1=0.1, b= 0.975, e=0.4, d = 0.12, g= 0.061, h = 0.21)
params1c = c(a1 = 0.0, b= 0.975, e = 0.4, d = 0.12, g = 0.061, h = 0.2)
tf <- 30
times <- seq(0,tf, by = 1)
t0.1 <- c(S = .97, L = 0.004, I = 0.006, D = 0.09, U1 = 0.005, U2 = 0.005, Rd = 0, Ru = 0, F = 0.006) #canada


#test the 1 susceptible within first 2
susceptODE1a.df <- as.data.frame(ode(y = t0.1, times, func = totalPop, parms = params1a))
susceptible1a <- susceptODE1a.df %>% 
  select(time, S) %>% 
  rename(susceptible = S)
undetected1a <- susceptODE1a.df %>% 
  select(time, L, I, U1, U2,  Ru) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected1a <- susceptODE1a.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
high <- merge(susceptible1a, undetected1a) 
high <- merge(high, detected1a) %>% 
  melt(id = "time") %>% 
  rename(high = value, state = variable)

#test the 1 susceptible within first 5 days
susceptODE1b.df <- as.data.frame(ode(y=t0.1, times, func = totalPop, parms = params1b)) 
susceptible1b <- susceptODE1b.df %>% 
  select(time, S) %>% 
  rename(susceptible = S)
undetected1b <- susceptODE1b.df %>% 
  select(time, L, I, U1, U2, Ru) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected1b <- susceptODE1b.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
medium <- merge(susceptible1b, undetected1b) 
medium <- merge(medium, detected1b) %>% 
  melt(id = "time") %>% 
  rename(medium = value, state = variable)

#test the 1 susceiptible within days 5-10
susceptODE1c.df <- as.data.frame(ode(y=t0.1, times, func = totalPop, parms = params1c)) 
susceptible1c <- susceptODE1c.df %>% 
  select(time, S) %>% 
  rename(susceptible = S)
undetected1c <- susceptODE1c.df %>% 
  select(time, L, I, U1, U2, Ru) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected1c <- susceptODE1c.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
low <- merge(susceptible1c, undetected1c) 
low <- merge(low, detected1c) %>% 
  melt(id = "time") %>% 
  rename(low = value, state = variable)

forGraph <- merge(high, medium, by = c('time', 'state'))
forGraph <- merge(forGraph, low, by = c('time', 'state')) %>% 
  melt(id = c('time', 'state')) %>% 
  rename(restrictionLevel = variable, probability = value) %>% 
  filter(time == 30) 

ggplot(forGraph, aes(x=restrictionLevel, y = probability, fill = state)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c('slateblue4', 'orange', "chartreuse4")) + 
  theme_minimal() + 
  geom_text(aes(label = label_percent()(probability)), vjust = 2.5, colour = "white") +
  scale_y_continuous(labels = label_percent())

# ARRIVAL AT NB BORDER
tail(susceptODE1a.df)
tail(susceptODE1b.df)
tail(susceptODE1c.df)



#at day 30, assuming if someone is detected/dead they aren't still entering
