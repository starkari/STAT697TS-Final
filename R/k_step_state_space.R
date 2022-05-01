# using k_step ahead forecast to find which state_space model minimizes forecast error


library(tidyverse)
library(MARSS)

load("Data/data.90.percent.Rdata")


nsubset <- 10


k.step.state.space.selection <- matrix(NA,nrow = nsubset, ncol = 3)


c.names <- c("No_Cov","Cov_in_SE", "Cov_in_OE")

colnames(k.step.state.space.selection) <- c.names


for (i in 1:nsubset) {
  
  data.subset <- data.90.percent[(((i-1)*390 + 1):((i)*390)),]
  data.test <- data.subset[-(1:360),]
  y.test <- data.test$TOBS
  
  y.marss.k.step <- data.subset$TOBS
  y.marss.k.step[-(1:360)] <- NA
  
  covariates.k.step <- t( model.matrix(~month,data=data.subset)[, -1])
  
  # model definitions
  model.no.cov.k <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                       Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                       x0=matrix("mu"), tinitx=1)
  
  model.cov.se.k <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                          Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                          x0=matrix("mu"), tinitx=1,
                          C="unconstrained", c=covariates.k.step)
  
  model.cov.oe.k <- list(B=matrix("phi"), U=matrix(0), Q=matrix("sig.sq.w"), 
                          Z=matrix("a"), A=matrix(0), R=matrix("sig.sq.v"),
                          x0=matrix("mu"), tinitx=1,
                          D="unconstrained", d=covariates.k.step)
  # fitting/predicting
  
  fit.no.cov.k <- MARSS(c(y.marss.k.step), model=model.no.cov.k, method = "kem")
  fit.no.cov.k <- MARSS(c(y.marss.k.step), model=model.no.cov.k,  method = "BFGS",
                         inits = fit.no.cov.k)
  pred.no.cov.k <- fit.no.cov.k$ytT[-(1:360)]
  
  
  fit.cov.se.k <- MARSS(c(y.marss.k.step), model=model.cov.se.k, method = "kem")
  fit.cov.se.k <- MARSS(c(y.marss.k.step), model=model.cov.se.k,  method = "BFGS",
                        inits = fit.cov.se.k)
  pred.cov.se.k <- fit.cov.se.k$ytT[-(1:360)]
  
  
  fit.cov.oe.k <- MARSS(c(y.marss.k.step), model=model.cov.oe.k, method = "kem")
  fit.cov.oe.k <- MARSS(c(y.marss.k.step), model=model.cov.oe.k,  method = "BFGS",
                        inits = fit.cov.oe.k)
  pred.cov.oe.k <- fit.cov.oe.k$ytT[-(1:360)]
  
  
  #  error measures
  
  mse.no.cov.k <- mean((y.test-pred.no.cov.k)^2)
  mse.cov.se.k <- mean((y.test-pred.cov.se.k)^2)
  mse.cov.oe.k <- mean((y.test-pred.cov.oe.k)^2)
  
  # store values
  k.step.state.space.selection[i,] <-c(mse.no.cov.k,mse.cov.se.k,mse.cov.oe.k)
  
}


save(k.step.state.space.selection,file="Data/k.step.state.space.selection.Rdata")




