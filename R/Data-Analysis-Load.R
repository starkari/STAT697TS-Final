
# this r script loads the results from scripts that takes awhile to 
# run so when sourcing this file in the final report it compiles quickly

library(tidyverse)
library(readr)
library(forecast)
library(R.utils)
library(MARSS)

load("Data/data.Rdata")
load("Data/data.90.percent.Rdata")
load("Data/missing_data.Rdata")
load("Data/k.step.arima.selection.Rdata")
load("Data/fit.covy.Rdata")

## ------------ data overview
y.full <- data$TOBS

# data plot
data.plot <- data %>% 
  ggplot(aes(x=DATE, y= TOBS))+
  geom_line() +
  xlab("Date")+
  ylab("Temperature")

# max date in 90%
max(data.90.percent$DATE)

## ______________ what ARIMA model

# k step model selection from 90% of data

which.min(colMeans(k.step.arima.selection))

# The k=30 step ahead forecasts error minimizing model is 
# ARIMA(4,0,2) with a sinusoidal trend


## -------------------- model fitting ARIMA

xreg.90 <- model.matrix(~sin(2*pi*omega) +
               cos(2*pi*omega)
             ,data=data.90.percent)

y.90 <- data.90.percent$TOBS

arima.fit <- arima(y.90, order = c(4, 0, 2), 
                   xreg = xreg.90, method = "ML", include.mean = FALSE)


arima.fit.coef <- arima.fit$coef

## ---------------------- model prediction ARIMA


# compare 30 day predictions and 1 day predictions for last 10% of data
y.full <- data$TOBS
n.full <- nrow(data)
n.90 <- nrow(data.90.percent)

indices.to.predict.from <- (n.90+1):(n.full-30) # gets indices in held out data minus last 30 days to predict from

prediction.arima <- as.data.frame(matrix(NA,nrow = 2*length(indices.to.predict.from),
                           ncol = 4))

colnames(prediction.arima) <- c("horizon","MSE","MAE","num_in_95CI")

for (i in 1:length(indices.to.predict.from)) {
  index <- indices.to.predict.from[i]
  # fitting
  y.fitting <- y.full[1:index]
  xreg.fit <- model.matrix(~sin(2*pi*omega) +
                            cos(2*pi*omega)
                          ,data=data[1:index,])

  # predict 1
  y.predict.1 <- y.full[index+1]
  data.predict.1 <- data[index+1,]
  xreg.predict.1 <- model.matrix(~sin(2*pi*omega) +
                             cos(2*pi*omega)
                           ,data=data.predict.1)
  ### predicts 1 step ahead using coefficients from fit to first 90%
  predict.1 <- predict(arima(y.fitting, order = c(4, 0, 2), xreg = xreg.fit,
                             method = "ML", fixed = arima.fit.coef,
                             include.mean = FALSE),
                       n.ahead = 1, newxreg = xreg.predict.1)
  pred.1 <- as.numeric(predict.1$pred)
  se.pred.1 <- as.numeric(predict.1$se)
  true.val.predict.1 <- y.full[index+1]
  MSE.predict.1 <- as.numeric((pred.1-true.val.predict.1)^2)
  MAE.predict.1 <- as.numeric(abs(pred.1-true.val.predict.1))
  pred.1.lower.bound <- pred.1 + qnorm(0.025)*se.pred.1
  pred.1.upper.bound <- pred.1 + qnorm(0.975)*se.pred.1
  pred.1.ci <- ifelse((pred.1<=pred.1.upper.bound & pred.1>=pred.1.lower.bound),1,0)


  # predict 30
  y.predict.30 <- y.full[(index+1):(index+30)]
  data.predict.30 <- data[(index+1):(index+30),]
  xreg.predict.30 <- model.matrix(~sin(2*pi*omega) +
                                   cos(2*pi*omega)
                                 ,data=data.predict.30)
  ### predicts 30 steps ahead using coefficients from fit to first 90%
  predict.30 <- predict(arima(y.fitting, order = c(4, 0, 2), xreg = xreg.fit,
                             method = "ML", fixed = arima.fit.coef,
                             include.mean = FALSE),
                       n.ahead = 30, newxreg = xreg.predict.30)

  pred.30 <- as.numeric(predict.30$pred)
  se.pred.30 <- as.numeric(predict.30$se)
  true.val.predict.30 <- y.full[index+1]
  MSE.predict.30 <- as.numeric(mean((pred.30-true.val.predict.30)^2))
  MAE.predict.30 <- as.numeric(mean(abs(pred.30-true.val.predict.30)))
  pred.30.ci <- 0
  # count number of times prediction in CI
  for (j in 1:30) {
    pred.j <- pred.30[j]
    se.pred.j <- se.pred.30[j]
    pred.j.lower.bound <- pred.j + qnorm(0.025)*se.pred.j
    pred.j.upper.bound <- pred.j + qnorm(0.975)*se.pred.j
    pred.j.ci <- ifelse((pred.j<=pred.j.upper.bound & pred.j>=pred.j.lower.bound),
                        1,0)
    pred.30.ci <- pred.30.ci+pred.j.ci
  }
  pred.30.ci.pc <- pred.30.ci/30
  
  prediction.arima[2*(i-1)+1,] <- c(1,MSE.predict.1,
                                    MAE.predict.1,pred.1.ci)
  prediction.arima[2*(i-1)+2,] <- c(30,MSE.predict.30,
                                    MAE.predict.30,pred.30.ci.pc)
  
}


# predict entire 10%
predict.last.10.percent <- predict(arima(y.90, order = c(4, 0, 2), xreg = xreg.90,
                                         method = "ML", fixed = arima.fit.coef,
                                         include.mean = FALSE),
                                   n.ahead = n.full-n.90, 
                                   newxreg = model.matrix(~sin(2*pi*omega) +
                                                            cos(2*pi*omega)
                                                          ,data=data[n.full-n.90,]))
### MSE last 10%
mean(as.vector((predict.last.10.percent$pred-y.full[n.full-n.90])^2))

#### MAE last 10%
mean(as.vector(abs(predict.last.10.percent$pred-y.full[n.full-n.90])))

### ---------------- State Space


forc.covy <- fit.covy$ytT
forc.covy.se <- fit.covy$ytT.se

plot(data$DATE, data$TOBS, type = "l",
     xlab = "Date", ylab = "Temp",
     col = "gray", ylim = c(-10, 100))
lines(data$DATE[(n.90+1):n.full], forc.covy[(n.90+1):n.full], col = "blue")


### MSE last 10%
mean(as.vector((forc.covy-y.full[n.full-n.90])^2))

#### MAE last 10%
mean(as.vector(abs(forc.covy-y.full[n.full-n.90])))



