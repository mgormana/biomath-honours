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
param1b <- c(a = 0.02, b = 0.975, e = 2/5, d = 2/(3-(1/.4)), g = 0.061)
param1c <- c(a = 0.02, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

param2a <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(21-(1/.4)), g = 0.061)
param2b <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(3-(1/.4)), g = 0.061)
param2c <- c(a = 0.5, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

param3a <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(21-(1/.4)), g = 0.061)
param3b <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(3-(1/.4)), g = 0.061)
param3c <- c(a = 0.9, b = 0.975, e = 2/5, d = 2/(33-(1/.4)), g = 0.061)

times <- seq(0,14, by = .1)
t0 <- c(L = 1, I = 0, D = 0, U = 0, Rd = 0, Ru = 0, F = 0)

res_ode1a.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param1a))
preinfectious1a.df <- res_ode1a.df %>% 
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected1a.df <- res_ode1a.df %>% 
  select(time, I, U, Ru) %>% 
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
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected1b.df <- res_ode1b.df %>% 
  select(time, I, U, Ru) %>% 
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
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected1c.df <- res_ode1c.df %>% 
  select(time, I, U, Ru) %>% 
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
  rename("22 Days" = stateProb.x,"3 Days" = stateProb.y, "31 Days" = stateProb) %>% 
  gather(InfectiousPeriodLength, byDay2, -time, -state)
rm(res_ode1a.df, res_ode1b.df, res_ode1c.df, 
   undetected1a.df, undetected1b.df, undetected1c.df,
   detected1a.df, detected1b.df, detected1c.df,
   preinfectious1a.df, preinfectious1b.df, preinfectious1c.df)

res_ode2a.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param2a))
preinfectious2a.df <- res_ode2a.df %>% 
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected2a.df <- res_ode2a.df %>% 
  select(time, I, U, Ru) %>% 
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
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected2b.df <- res_ode2b.df %>% 
  select(time, I, U, Ru) %>% 
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
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected2c.df <- res_ode2c.df %>% 
  select(time,I,  U, Ru) %>% 
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
  rename("22 Days" = stateProb.x, "3 Days" = stateProb.y, "31 Days" = stateProb) %>% 
  gather(InfectiousPeriodLength, byDay5, -time, -state)
rm(detected2a.df, detected2b.df, detected2c.df,
   preinfectious2a.df, preinfectious2b.df, preinfectious2c.df,
   undetected2a.df, undetected2b.df, undetected2c.df,
   res_ode2a.df, res_ode2b.df, res_ode2c.df)

res_ode3a.df <- as.data.frame(ode(y = t0, times, func = asymptomatics, param3a))
preinfectious3a.df <- res_ode3a.df %>% 
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected3a.df <- res_ode3a.df %>% 
  select(time, I, U, Ru) %>% 
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
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected3b.df <- res_ode3b.df %>% 
  select(time, I, U, Ru) %>% 
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
  select(time, L) %>% 
  melt(id = "time") %>% 
  group_by(time) %>% 
  summarise(preinfectious = sum(value))
undetected3c.df <- res_ode3c.df %>% 
  select(time,I,  U, Ru) %>% 
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
  rename("22 Days" = stateProb.x, "3 Days" = stateProb.y, "31 Days" = stateProb) %>% 
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
  mutate(InfectiousPeriodLength = fct_relevel(InfectiousPeriodLength, "3 Days", "22 Days", "31 Days"))

labsTest <- c("Day 2", "By Day 5", "Day 8")
names(labsTest) = c("byDay2", "byDay5", "afterDay5")

#plotting 
ggplot(res_ode.df, aes(x = time, y = StateProbability, color = state)) + 
  geom_line() + 
  scale_color_manual(values = c("#f00000", "#ffc709","slateblue4"), 
                     labels = c("Pre-Infectious", "Infectious and Undetected", "Detected")) +
  facet_grid(rows = vars(InfectiousPeriodLength), 
             cols = vars(TestingOption), 
             labeller = labeller(TestingOption = labsTest)) + 
  scale_y_continuous(breaks = seq(0,1,.5), labels = scales::percent) + 
  xlim(0,14) +
  theme_minimal() +
  labs(x = "Day", y = "Probability", color = "State of Infection") +
  theme(axis.text = element_text(size = 14), 
        axis.title =element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 18),
        legend.position = "top",
        panel.spacing = unit(2, "lines")) + 
  guides(color = guide_legend(title.position = "top"))
  

