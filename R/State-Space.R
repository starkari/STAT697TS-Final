
# this r script loads the results from scripts that takes awhile to 
# run so when sourcing this file in the final report it compiles quickly

library(tidyverse)
library(MARSS)

load("Data/data.Rdata")
load("Data/data.90.percent.Rdata")



y.full <- data$TOBS
y.90 <- data.90.percent$TOBS
n.full <- nrow(data)
n.90 <- nrow(data.90.percent)


## State Space 90% no covariates

model.ss90nc <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   x0=matrix("mu"), tinitx=1)

fit.ss90nc <- MARSS(c(y.90), model=model.ss90nc, method = "kem")
fit.ss90nc <- MARSS(c(y.90), model=model.ss90nc,  method = "BFGS",
                  inits = fit.ss90nc)



## State Space 90% month covariates in state equation

covariates.90 <- t( model.matrix(~month
                              ,data=data.90.percent)[, -1])
model.ss90covse <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                     Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                     x0=matrix("mu"), tinitx=1,
                     C="unconstrained", c=covariates.90)

fit.ss90covse <- MARSS(c(y.90), model=model.ss90covse, method = "kem")
fit.ss90covse <- MARSS(c(y.90), model=model.ss90covse,  method = "BFGS",
                    inits = fit.ss90covse)


## State Space 90% month covariates in observation equation

model.ss90covoe <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                        Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                        x0=matrix("mu"), tinitx=1,
                        D="unconstrained", d=covariates.90)

fit.ss90covoe <- MARSS(c(y.90), model=model.ss90covoe, method = "kem")
fit.ss90covoe <- MARSS(c(y.90), model=model.ss90covoe,  method = "BFGS",
                       inits = fit.ss90covoe)


State_Space_AIC_Fit90 <- data.frame("Model"=c("No Covariates","Month Cov. State Eq.",
                                              "Month Cov. Obs. Eq." ),
                                    "AIC"=c(fit.ss90nc$AIC, fit.ss90covse$AIC,
                                            fit.ss90covoe$AIC),
                                    "AICc"=c(fit.ss90nc$AICc, fit.ss90covse$AICc,
                                            fit.ss90covoe$AICc))

save(State_Space_AIC_Fit90,file="Data/State_Space_AIC_Fit90.Rdata")


## State Space Prediciton

y.marss <- y.full
y.marss[(n.90+1):n.full] <- NA


# Adding covariates of indicators for each month
covariates <- t( model.matrix(~month
                              ,data=data)[, -1])

# covariates in the state part of the equation as it was AIC selected
model.state.full <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                   Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                   x0=matrix("mu"), tinitx=1, 
                   C="unconstrained", c=covariates)

fit.state.full <- MARSS(c(y.marss), model=model.state.full, method = "kem")
fit.state.full <- MARSS(c(y.marss), model=model.state.full,  method = "BFGS",
                  inits = fit.state.full)

save(fit.covy,file="Data/fit.state.full.Rdata")


