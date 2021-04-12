library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

asymptomatics <- function(t, x, parms = NULL) {
  with(as.list(c(x, parms)), {
    
    dL = -(e*L)
    dI = -(e*I) + (e*L)
    dD = (a*e)*I + b*d*U - g*d*D - (1-g)*d*D
    dU = (1-a)*e*I - b*d*U - (1-b)*d*U
    dRd = (1-g)*d*D 
    dRu = (1-b)*d*U
    dF = g*d*D
    list(c(dL, dI, dD, dU, dRd, dRu, dF))
  })
}

#mean duration of infection
param1a <- c(a = 0.02, b = 0.975, e = 2/5, d = 2/(21-(1/.4)), g = 0.061)
param1b <- c(a = 0.02, b = 0.975, e = 2/5, d = 2/(2.6-(1/.4)), g = 0.061)
param1c <- c(a = 0.02, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

param2a <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(21-(1/.4)), g = 0.061)
param2b <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(2.6-(1/.4)), g = 0.061)
param2c <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

param3a <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(21-(1/.4)), g = 0.061)
param3b <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(2.6-(1/.4)), g = 0.061)
param3c <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

times <- seq(0,30, by = 1)
t0 <- c(L = 1, I = 0, D = 0, U = 0, Rd = 0, Ru = 0, F = 0)

res_ode1a.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param1a))
preinfectious1a.df <- res_ode1a.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected1a.df <- res_ode1a.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected1a.df <- res_ode1a.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode1a.df <- merge(preinfectious1a.df, undetected1a.df, by.x = c("time"))
res_ode1a.df <- merge(res_ode1a.df, detected1a.df, by.x = c('time'))
res_ode1a.df <- melt(res_ode1a.df, id = "time")
res_ode1a.df <- rename(res_ode1a.df, state = variable, stateProb = value)

res_ode1b.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param1b))
preinfectious1b.df <- res_ode1b.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected1b.df <- res_ode1b.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected1b.df <- res_ode1b.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode1b.df <- merge(preinfectious1b.df, undetected1b.df, by.x = c("time"))
res_ode1b.df <- merge(res_ode1b.df, detected1b.df, by.x = c('time'))
res_ode1b.df <- melt(res_ode1b.df, id = "time")
res_ode1b.df <- rename(res_ode1b.df, state = variable, stateProb = value)

res_ode1c.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param1c))
preinfectious1c.df <- res_ode1c.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected1c.df <- res_ode1c.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected1c.df <- res_ode1c.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode1c.df <- merge(preinfectious1c.df, undetected1c.df, by.x = c("time"))
res_ode1c.df <- merge(res_ode1c.df, detected1c.df, by.x = c('time'))
res_ode1c.df <- melt(res_ode1c.df, id = "time")
res_ode1c.df <- rename(res_ode1c.df, state = variable, stateProb = value)

res_ode1.df <- merge(res_ode1a.df, res_ode1b.df, by= c("time", "state"))
res_ode1.df <- merge(res_ode1.df, res_ode1c.df, by = c("time", "state"))
res_ode1.df <- res_ode1.df %>% 
  rename(median = stateProb.x, short = stateProb.y, long = stateProb) %>% 
  gather(InfectiousPeriodLength, byDay2, -time, -state)
rm(res_ode1a.df, res_ode1b.df, res_ode1c.df, 
   undetected1a.df, undetected1b.df, undetected1c.df,
   detected1a.df, detected1b.df, detected1c.df,
   preinfectious1a.df, preinfectious1b.df, preinfectious1c.df)

res_ode2a.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param2a))
preinfectious2a.df <- res_ode2a.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected2a.df <- res_ode2a.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected2a.df <- res_ode2a.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode2a.df <- merge(preinfectious2a.df, undetected2a.df, by.x = c("time"))
res_ode2a.df <- merge(res_ode2a.df, detected2a.df, by.x = c('time'))
res_ode2a.df <- melt(res_ode2a.df, id = "time")
res_ode2a.df <- rename(res_ode2a.df, state = variable, stateProb = value)

res_ode2b.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param2b))
preinfectious2b.df <- res_ode2b.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected2b.df <- res_ode2b.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected2b.df <- res_ode2b.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode2b.df <- merge(preinfectious2b.df, undetected2b.df, by.x = c("time"))
res_ode2b.df <- merge(res_ode2b.df, detected2b.df, by.x = c('time'))
res_ode2b.df <- melt(res_ode2b.df, id = "time")
res_ode2b.df <- rename(res_ode2b.df, state = variable, stateProb = value)

res_ode2c.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param2c))
preinfectious2c.df <- res_ode2c.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected2c.df <- res_ode2c.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected2c.df <- res_ode2c.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode2c.df <- merge(preinfectious2c.df, undetected2c.df, by.x = c("time"))
res_ode2c.df <- merge(res_ode2c.df, detected2c.df, by.x = c('time'))
res_ode2c.df <- melt(res_ode2c.df, id = "time")
res_ode2c.df <- rename(res_ode2c.df, state = variable, stateProb = value)

res_ode2.df <- merge(res_ode2a.df, res_ode2b.df, by= c("time", "state"))
res_ode2.df <- merge(res_ode2.df, res_ode2c.df, by = c("time", "state"))
res_ode2.df <- res_ode2.df %>% 
  rename(median = stateProb.x, short = stateProb.y, long = stateProb) %>% 
  gather(InfectiousPeriodLength, byDay5, -time, -state)
rm(detected2a.df, detected2b.df, detected2c.df,
   preinfectious2a.df, preinfectious2b.df, preinfectious2c.df,
   undetected2a.df, undetected2b.df, undetected2c.df,
   res_ode2a.df, res_ode2b.df, res_ode2c.df)

res_ode3a.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param3a))
preinfectious3a.df <- res_ode3a.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected3a.df <- res_ode3a.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected3a.df <- res_ode3a.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode3a.df <- merge(preinfectious3a.df, undetected3a.df, by.x = c("time"))
res_ode3a.df <- merge(res_ode3a.df, detected3a.df, by.x = c('time'))
res_ode3a.df <- melt(res_ode3a.df, id = "time")
res_ode3a.df <- rename(res_ode3a.df, state = variable, stateProb = value)

res_ode3b.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param3b))
preinfectious3b.df <- res_ode3b.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected3b.df <- res_ode3b.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected3b.df <- res_ode3b.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode3b.df <- merge(preinfectious3b.df, undetected3b.df, by.x = c("time"))
res_ode3b.df <- merge(res_ode3b.df, detected3b.df, by.x = c('time'))
res_ode3b.df <- melt(res_ode3b.df, id = "time")
res_ode3b.df <- rename(res_ode3b.df, state = variable, stateProb = value)

res_ode3c.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param3c))
preinfectious3c.df <- res_ode3c.df %>% 
  select(time, L, I) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected3c.df <- res_ode3c.df %>% 
  select(time, U, Ru) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(undetected = sum(value))
detected3c.df <- res_ode3c.df %>% 
  select(time, D, Rd, F) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(detected = sum(value))
res_ode3c.df <- merge(preinfectious3c.df, undetected3c.df, by.x = c("time"))
res_ode3c.df <- merge(res_ode3c.df, detected3c.df, by.x = c('time'))
res_ode3c.df <- melt(res_ode3c.df, id = "time")
res_ode3c.df <- rename(res_ode3c.df, state = variable, stateProb = value)

res_ode3.df <- merge(res_ode3a.df, res_ode3b.df, by= c("time", "state"))
res_ode3.df <- merge(res_ode3.df, res_ode3c.df, by = c("time", "state"))
res_ode3.df <- res_ode3.df %>% 
  rename(median = stateProb.x, short = stateProb.y, long = stateProb) %>% 
  gather(InfectiousPeriodLength, afterDay5, -time, -state)
rm(undetected3a.df, undetected3b.df, undetected3c.df,
   detected3a.df, detected3b.df, detected3c.df,
   preinfectious3a.df, preinfectious3b.df, preinfectious3c.df,
   res_ode3a.df, res_ode3b.df, res_ode3c.df)

res_ode.df <- merge(res_ode1.df, res_ode2.df)
res_ode.df <- merge(res_ode.df, res_ode3.df)
res_ode.df <- res_ode.df %>%
  gather(TestingOption, StateProbability, -time, -state, -InfectiousPeriodLength) %>% 
  mutate(TestingOption = fct_relevel(TestingOption, "byDay2", "byDay5", "afterDay5")) %>% 
  mutate(InfectiousPeriodLength = fct_relevel(InfectiousPeriodLength, "short", "median", "long"))

#plotting 
ggplot(res_ode.df, aes(x = time, y = StateProbability, color = state)) + 
  geom_line() + 
  scale_color_manual(values = c("olivedrab2", "slateblue4", "orange", "magenta4")) +
  facet_wrap(~ InfectiousPeriodLength + TestingOption) + 
  theme_minimal()

stillUndetected1a <- as.data.frame(ode(y = t0, times, func = asymptomatics, param1a)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(median = sum(value))
stillUndetected1b <- as.data.frame(ode(y = t0, times, func = asymptomatics, param1b)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(short = sum(value))
stillUndetected1c <- as.data.frame(ode(y = t0, times, func = asymptomatics, param1a)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(long = sum(value))
byDay2 <- merge(stillUndetected1a, stillUndetected1b, by = 'time')
byDay2 <- merge(byDay2, stillUndetected1c, by = 'time') %>% 
  melt(id = 'time') %>% 
  filter(time == 2) %>% 
  select(-time) %>% 
  rename(day2 = value)
stillUndetected2a <- as.data.frame(ode(y = t0, times, func = asymptomatics, param2a)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(median = sum(value))
stillUndetected2b <- as.data.frame(ode(y = t0, times, func = asymptomatics, param2b)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(short = sum(value))
stillUndetected2c <- as.data.frame(ode(y = t0, times, func = asymptomatics, param2c)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(long = sum(value))
byDay5 <- merge(stillUndetected2a, stillUndetected2b, by = 'time')
byDay5 <- merge(byDay5, stillUndetected2c, by = 'time') %>% 
  melt(id = 'time') %>% 
  filter(time == 5) %>% 
  select(-time) %>% 
  rename(day5 = value)
stillUndetected3a <- as.data.frame(ode(y = t0, times, func = asymptomatics, param3a)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(median = sum(value))
stillUndetected3b <- as.data.frame(ode(y = t0, times, func = asymptomatics, param3b)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(short = sum(value))
stillUndetected3c <- as.data.frame(ode(y = t0, times, func = asymptomatics, param3c)) %>% 
  select(time, L, I, U) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(long = sum(value))
afterDay5 <- merge(stillUndetected3a, stillUndetected3b, by = 'time')
afterDay5 <- merge(afterDay5, stillUndetected2c, by = 'time') %>% 
  melt(id = 'time') %>% 
  filter(time == 7) %>% 
  select(-time) %>% 
  rename(day7 = value)
afterDay5x <- merge(stillUndetected3a, stillUndetected3b, by = 'time')
afterDay5x <- merge(afterDay5x, stillUndetected3c, by = 'time') %>% 
  melt(id = 'time') %>% 
  filter(time == 14) %>% 
  select(-time) %>% 
  rename(day14 = value)
afterDay10 <- merge(stillUndetected3a, stillUndetected3b, by = 'time') 
afterDay10 <- merge(afterDay10, stillUndetected3c, by = 'time') %>% 
  melt(id = 'time') %>% 
  filter(time == 10) %>% 
  select(-time) %>% 
  rename(day10 = value)
stillUndetected <- merge(byDay2, byDay5, by = 'variable')
stillUndetected <- merge(stillUndetected, afterDay5, by = 'variable')
stillUndetected <- merge(stillUndetected, afterDay5x, by = 'variable') 
stillUndetected <- merge(stillUndetected, afterDay10, by = 'variable') %>% 
  rename(infectiousPeriod = variable) %>% 
  melt(id = 'infectiousPeriod') %>% 
  rename(quarantine = variable, probability = value) %>% 
  mutate(infectiousPeriod = fct_relevel(infectiousPeriod, "short", "median", "long"))  %>% 
  mutate(quarantine = fct_relevel(quarantine, "day2", "day5", 'day7', 'day10', 'day14'))


ggplot(stillUndetected, aes(x = quarantine, y = probability, fill = infectiousPeriod)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  scale_fill_manual(values = c("olivedrab2", "slateblue4", "orange", "magenta4")) +
  theme_minimal()

#plot matrix exp
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




