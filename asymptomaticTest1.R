Library(Matrix::expm)

#alpha, sensitivity of test  (I want this to be the range, 0.67 to 0.99 from day 0 to 5)
a <- 0.99

#epsilon, incubation
e = 2/1.9

#delta, infectious
d = 2/(8.5-0.2)

#beta, detected symptoms if day Y = 11.5 (**this is off**)
b = 0.975

#gamma, fatality rate
g = 0.025

#probability of being in states, should this be different?
pI1 = 0.25
pI2 = 0.25
pU1 = 0.25
pU2 = 0.25
  
#transition matrix
V <- matrix(c(
  -e, 0, 0, 0, 0, 0, 0, 0, 0, 
   e, -e, 0, 0, 0, 0, 0, 0, 0, 
  0, a*e, -d, 0, 0, 0, 0, 0, 0, 
  0, 0, d, -d, b*d, 0, 0, 0, 0, 
  0, (1-a)*e, 0, 0, -d, 0, 0, 0, 0, 
  0, 0, 0, 0, (1-b)*d, -d, 0, 0, 0, 
  0, 0, 0, (1-g)*d, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, d, 0, 0, 0, 
  0, 0, 0, g*d, 0, 0, 0, 0, 0), 
  nrow = 9, ncol = 9, byrow = TRUE)

#probability vector 
undetProb = matrix(c(pI1, pI2, 0, 0, pU1, pU2, 0, 0, 0), ncol=1, byrow = FALSE)

#those who are undetected
undetRec = matrix(c(1,1,0,0,1,1,0,1,0), nrow = 1, byrow = TRUE)

#
undetRec_prob = function(t) {undetRec%*%expm(V*t)%*%undetProb}
t = seq(0,30,by=.1)
plot(t,sapply(t,undetRec_prob),ylim=c(0,1), main = 'a = 0.67')
