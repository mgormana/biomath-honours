#odesolver2
a = 0.5
b = 0.975
d = 0.1
g = 0.061
e = 0.4
h = 0.8
S = 100

totalPop <- function(t, x, parms = NULL) {
  dS = -h*(I+U)*S
  dL = h*(I+U)*S - (e*L)
  dI = -(e*I) + (e*L)
  dD = (a*e)*I + b*d*U - g*d*D - (1-g)*d*D
  dR = (1-g)*d*D + (1-b)*d*U
  dF = g*d*D
  list(c(dS, dL, dI, dD, dR, dF))
}

#transition matrix
V <-matrix(c(
  -e, 0, 0, 0,
  e, -e, 0, 0,
  0, (1-a)*e, 0, -d,
  nrow = 4, ncol = 4, byrow = TRUE))

w = matrix(c(1, 0, 0, 0), ncol =1, byrow = FALSE)
beta = matrix(c(0, h, 0, h), nrow =1, byrow = TRUE)
F = (w %*% beta)*S

dS = -h*(I+U)*S
dx = F*x + V*x

Ro = beta %*% V^(-1)%*%w%*%S
