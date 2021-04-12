#odesolver2
library(deSolve)
library(ggplot2)
library(dplyr)
library(reshape2)
library(scales)

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

params1a = c(a1=0.25, b= 0.975, e=0.4, d = 0.12, g= 0.061, h = 0.2)
params1b = c(a1=0.1, b= 0.975, e=0.4, d = 0.12, g= 0.061, h = 0.2)
params1c = c(a1 = 0, b= 0.975, e = 0.4, d = 0.12, g = 0.061, h = 0.2)
tf <- 30
times <- seq(0,tf, by = 1)
t0.1 <- c(S = .97, L = 0.004, I = 0.006, D = 0.09, U = 0.01, R = 0, F = 0.006) #canada


#test the 1 susceptible within first 2
susceptODE1a.df <- as.data.frame(ode(y = t0.1, times, func = totalPop, parms = params1a))
susceptible1a <- susceptODE1a.df %>% 
  select(time, S) %>% 
  rename(susceptible = S)
undetected1a <- susceptODE1a.df %>% 
  select(time, L, I, U) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
high <- merge(susceptible1a, undetected1a) %>% 
  melt(id = "time") %>% 
  rename(high = value, state = variable)

#test the 1 susceptible within first 5 days
susceptODE1b.df <- as.data.frame(ode(y=t0.1, times, func = totalPop, parms = params1b)) 
susceptible1b <- susceptODE1b.df %>% 
  select(time, S) %>% 
  rename(susceptible = S)
undetected1b <- susceptODE1b.df %>% 
  select(time, L, I, U) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
medium <- merge(susceptible1b, undetected1b) %>% 
  melt(id = "time") %>% 
  rename(medium = value, state = variable)

#test the 1 susceiptible within days 5-10
susceptODE1c.df <- as.data.frame(ode(y=t0.1, times, func = totalPop, parms = params1c)) 
susceptible1c <- susceptODE1c.df %>% 
  select(time, S) %>% 
  rename(susceptible = S)
undetected1c <- susceptODE1c.df %>% 
  select(time, L, I, U) %>% 
  melt(id = 'time') %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
low <- merge(susceptible1c, undetected1c) %>% 
  melt(id = "time") %>% 
  rename(low = value, state = variable)

forGraph <- merge(high, medium, by = c('time', 'state'))
forGraph <- merge(forGraph, low, by = c('time', 'state')) %>% 
  melt(id = c('time', 'state')) %>% 
  rename(restrictionLevel = variable, probability = value) %>% 
  filter(time == 30) 

ggplot(forGraph, aes(x=restrictionLevel, y = probability, fill = state)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = c('slateblue4', 'orange')) + 
  theme_minimal() + 
  geom_text(aes(label = label_percent()(probability)), vjust = 2.5, colour = "white") +
  scale_y_continuous(labels = label_percent())

# ARRIVAL AT NB BORDER
tail(susceptODE1a.df)
tail(susceptODE1b.df)
tail(susceptODE1c.df)

#at day 30, assuming if someone is detected/dead they aren't still entering
