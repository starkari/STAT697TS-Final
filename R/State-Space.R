
# this r script loads the results from scripts that takes awhile to 
# run so when sourcing this file in the final report it compiles quickly

library(tidyverse)
library(MARSS)

load("Data/data.Rdata")
load("Data/data.90.percent.Rdata")



y.full <- data$TOBS
y.full <- data$TOBS
n.full <- nrow(data)
n.90 <- nrow(data.90.percent)

y.marss <- y.full
y.marss[(n.90+1):n.full] <- NA


# Adding covariates of indicators for each month
covariates <- t( model.matrix(~month
                              ,data=data)[, -1])

# covariates in the observation part of the equation
model.covy <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   x0=matrix("mu"), tinitx=1, 
                   D="unconstrained", d=covariates)

fit.covy <- MARSS(c(y.marss), model=model.covy, method = "kem")
fit.covy <- MARSS(c(y.marss), model=model.covy,  method = "BFGS",
                  inits = fit.covy)

save(fit.covy,file="Data/fit.covy.Rdata")


