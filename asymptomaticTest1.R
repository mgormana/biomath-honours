library(expm)
library(ggplot2)
library(tidyverse)

#alpha, sensitivity of test  (I want this to be the range, 0.67 to 0.99 from day 0 to 5)
a1 = 0.02
a2 = 0.5
a3 = 0.9
e = 2/5
d = 2/(22 - (1/e))
b = 0.975
g = 0.061

#probability of being in states
pI1 = 0.25
pI2 = 0.25
pU1 = 0.25
pU2 = 0.25


#transition matrix
V1 <- matrix(c(
  -e, 0, 0, 0, 0, 0, 0, 0, 0, 
  e, -e, 0, 0, 0, 0, 0, 0, 0, 
  0, a1*e, -d, 0, 0, 0, 0, 0, 0, 
  0, 0, d, -d, b*d, 0, 0, 0, 0, 
  0, (1-a1)*e, 0, 0, -d, 0, 0, 0, 0, 
  0, 0, 0, 0, (1-b)*d, -d, 0, 0, 0, 
  0, 0, 0, (1-g)*d, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, d, 0, 0, 0, 
  0, 0, 0, g*d, 0, 0, 0, 0, 0), 
  nrow = 9, ncol = 9, byrow = TRUE)

V2 <- matrix(c(
  -e, 0, 0, 0, 0, 0, 0, 0, 0, 
  e, -e, 0, 0, 0, 0, 0, 0, 0, 
  0, a2*e, -d, 0, 0, 0, 0, 0, 0, 
  0, 0, d, -d, b*d, 0, 0, 0, 0, 
  0, (1-a2)*e, 0, 0, -d, 0, 0, 0, 0, 
  0, 0, 0, 0, (1-b)*d, -d, 0, 0, 0, 
  0, 0, 0, (1-g)*d, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, d, 0, 0, 0, 
  0, 0, 0, g*d, 0, 0, 0, 0, 0), 
  nrow = 9, ncol = 9, byrow = TRUE)

V3 <- matrix(c(
  -e, 0, 0, 0, 0, 0, 0, 0, 0, 
  e, -e, 0, 0, 0, 0, 0, 0, 0, 
  0, a3*e, -d, 0, 0, 0, 0, 0, 0, 
  0, 0, d, -d, b*d, 0, 0, 0, 0, 
  0, (1-a3)*e, 0, 0, -d, 0, 0, 0, 0, 
  0, 0, 0, 0, (1-b)*d, -d, 0, 0, 0, 
  0, 0, 0, (1-g)*d, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, d, 0, 0, 0, 
  0, 0, 0, g*d, 0, 0, 0, 0, 0), 
  nrow = 9, ncol = 9, byrow = TRUE)


#probability vector 
undetProb = matrix(c(pI1, pI2, 0, 0, pU1, pU2, 0, 0, 0), ncol=1, byrow = FALSE)

detectedRec = matrix(c(0,0,1,1,0,0,1,0,1), nrow = 1, byrow = TRUE)
IRec = matrix(c(1,1,0,0,0,0,0,0,0), nrow = 1, byrow = TRUE)
undetbyTestRec = matrix(c(0,0,0,0,1,1,0,0,0), nrow = 1, byrow = TRUE)


detectedRec_prob1 = function(t) {detectedRec%*%expm(V1*t)%*%undetProb}
IRec_prob1 = function(t) {-1*IRec %*% expm(V1*t)%*%undetProb}
undetByTestRec_prob1 = function(t) {-1*undetbyTestRec %*% expm(V1*t) %*% undetProb}

detectedRec_prob2 = function(t) {detectedRec%*%expm(V2*t)%*%undetProb}
IRec_prob2 = function(t) {-1*IRec %*% expm(V2*t)%*%undetProb}
undetByTestRec_prob2 = function(t) {-1*undetbyTestRec %*% expm(V2*t) %*% undetProb}

detectedRec_prob3 = function(t) {detectedRec%*%expm(V3*t)%*%undetProb}
IRec_prob3 = function(t) {-1*IRec %*% expm(V3*t)%*%undetProb}
undetByTestRec_prob3 = function(t) {-1*undetbyTestRec %*% expm(V3*t) %*% undetProb}

t = seq(0,30,by=1)

ProbM <- as.data.frame(t)

ProbM1a <-  as.data.frame(sapply(t, detectedRec_prob1))
ProbM1b <- as.data.frame(sapply(t, IRec_prob1))
ProbM1c <- as.data.frame(sapply(t, undetByTestRec_prob1))
ProbM1 <- cbind(ProbM, ProbM1a, ProbM1b, ProbM1c) 
ProbM1 <- ProbM1 %>% 
  filter( t == 5 |  t == 10 | t == 15 ) %>% 
  rename(Detected = "sapply(t, detectedRec_prob1)", 
         pretest = "sapply(t, IRec_prob1)", 
         Undetected = "sapply(t, undetByTestRec_prob1)") %>% 
  gather(State, byDay2, -t)


ProbM2a <- as.data.frame(sapply(t, detectedRec_prob2))
ProbM2b <- as.data.frame(sapply(t, IRec_prob2))
ProbM2c <- as.data.frame(sapply(t, undetByTestRec_prob2))
ProbM2 <- cbind(ProbM, ProbM2a, ProbM2b, ProbM2c)
ProbM2 <- ProbM2 %>% 
  filter( t == 5 |  t == 10 | t == 15  ) %>% 
  rename(Detected = "sapply(t, detectedRec_prob2)", 
         pretest = "sapply(t, IRec_prob2)", 
         Undetected = "sapply(t, undetByTestRec_prob2)") %>% 
  gather(State, byDay5, -t)

ProbM3a <- as.data.frame(sapply(t, detectedRec_prob3))
ProbM3b <- as.data.frame(sapply(t, IRec_prob3))
ProbM3c <- as.data.frame(sapply(t, undetByTestRec_prob3))
ProbM3 <- cbind(ProbM, ProbM3a, ProbM3b, ProbM3c)
ProbM3 <- ProbM3 %>% 
  filter( t == 5 |t == 10 | t == 15 ) %>% 
  rename(Detected = "sapply(t, detectedRec_prob3)", 
         pretest = "sapply(t, IRec_prob3)", 
         Undetected = "sapply(t, undetByTestRec_prob3)") %>% 
  gather(State, afterDay5, -t)

rm(ProbM1a, ProbM1b, ProbM1c, ProbM2a, ProbM2c, ProbM2b, ProbM3a, ProbM3b, ProbM3c)

ProbM <- merge(ProbM1, ProbM2)
ProbM <- merge(ProbM, ProbM3)
ProbM <- ProbM %>% gather(TestingOption, Probability, -t, -State) %>% 
  mutate(TestingOption = fct_relevel(TestingOption, c("afterDay5", "byDay5", "byDay2"))) %>% 
  filter(State != "pretest") %>% 
  mutate(State = fct_relevel(State, "Undetected"))

labs <- c("Day 5", "Day 10", "Day 15")
names(labs) <- c("5", "10", "15")

ggplot(ProbM, aes(x = TestingOption, y = Probability, fill = State)) + 
  geom_bar(position = "stack", stat = "identity") + 
  facet_wrap(~t, nrow = 3, labeller = labeller(t = labs)) + 
  scale_fill_manual(values = c("slateblue4","#ffc709")) +
  theme_minimal() +
  coord_flip() + 
  scale_x_discrete(name = "Day Tested", labels = c("Day 5-8", "Day 0-5", "Day 0-2")) +
  scale_y_continuous(name = "Probability", breaks = seq(-.75, .75, by = 0.25), labels = scales::percent) + 
  theme(axis.text = element_text(size = 12), 
        axis.title =element_text(size = 18), 
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 18),
        legend.position = "right") 


