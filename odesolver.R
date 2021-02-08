library("deSolve")
asymptomatics <- function(x, t, parms = NULL) {
  with(as.list(c(x, parms)), {
    
    #should i be including delta and episol in these?
    dI = -(e*I)
    dD = (a*e)*I + b*d*U - g*d*D - (1-g)*d*D
    dU = (1-a)*e*I - b*d*U - (1-b)*d*U
    dR = (1-g)*d*D + d*U
    dF = g*d*D
    list(c(dI, dD, dU, dR, dF))
  })
}

param <- list(c(a = 0.99, e = 2/1.9, d = 2/8.3, b = 0.975, g = 0.025))

tf <- 60
times <- seq(0,tf, by = 0.01)
t0 <- c(I = 1, D = 0, U = 0, R = 0, F = 0)

res_ode <- ode(y = t0, times, func = asymptomatics, param)
plot(res_ode[,1], res_ode[,I])