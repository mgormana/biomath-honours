library(expm)
library(patchwork)
# slope = speed of catching people 

#alpha, sensitivity of test  (I want this to be the range, 0.67 to 0.99 from day 0 to 5)
a1 = 0.02
a2 = 0.5
a3 = 0.9

#epsilon, incubation
e = 2/5

#delta, infectious
d = 2/(22 - (1/e))
d1 = 2/(2.6 - (1/e))
d2 = 2/(33 - (1/e))

#beta, detected symptoms by day 11
b = 0.975

#gamma, fatality rate
g = 0.061

#probability of being in states
pI1 = 0.25
pI2 = 0.25
pU1 = 0.25
pU2 = 0.25

pD1
  
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

#those who are undetected
undetRec = matrix(c(1,1,0,0,1,1,0,1,0), nrow = 1, byrow = TRUE)
detectedRec = matrix(c(0, 0, 1, 1, 0, 0, 1, 0, 1), nrow = 1, byrow = TRUE)

#
undetRec_prob1 = function(t) {undetRec%*%expm(V1*t)%*%undetProb}
detectedRec_prob1 = function(t) {detectedRec%*%expm(V1*t)%*%undetProb}

undetRec_prob2 = function(t) {undetRec%*%expm(V2*t)%*%undetProb}
detectedRec_prob2 = function(t) {detectedRec%*%expm(V2*t)%*%undetProb}

undetRec_prob3 = function(t) {undetRec%*%expm(V3*t)%*%undetProb}
detectedRec_prob3 = function(t) {detectedRec%*%expm(V3*t)%*%undetProb}

t = seq(0,30,by=.1)
t.df = as.data.frame(t)

undetProbM1 <- as.data.frame(sapply(t, undetRec_prob1))
undetProbM2 <- as.data.frame(sapply(t, undetRec_prob2))
undetProbM3 <- as.data.frame(sapply(t, undetRec_prob3))
undetProbM <- merge(t, undetProbM1)
undetProbM <- merge(undetProbM, undetProbM2)
undetProbM <- merge(undetProbM, undetProbM3)

undetplot1 <- qplot(t,sapply(t,undetRec_prob1),ylim=c(0,1))
ggplot(undetProbM, aes(x = t, y = sapply(t, undetRec_prob2))) + geom_line() + xlim(0,1)

undetplot2 <- qplot(t,sapply(t,undetRec_prob2),ylim=c(0,1))
ggplot(undetProbM, aes(x = t, y = sapply(t, undetRec_prob1))) + geom_line() + xlim(0,1)

undetplot3 <- qplot(t,sapply(t,undetRec_prob3),ylim=c(0,1))
ggplot(undetProbM, aes(x = t, y = sapply(t, undetRec_prob1))) + geom_line() + xlim(0,1)


detplot1 <- qplot(t, sapply(t,detectedRec_prob1), ylim = c(0,1))
detplot2 <- qplot(t, sapply(t,detectedRec_prob2), ylim = c(0,1))
detplot3 <- qplot(t, sapply(t,detectedRec_prob3), ylim = c(0,1))

