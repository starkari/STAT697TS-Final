
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
#load("Data/State_Space_AIC_Fit90.Rdata")
load("Data/k.step.state.space.selection.Rdata")
load("Data/fit.state.full.Rdata")

## ------------ data overview
y.full <- data$TOBS
n.full <- nrow(data)
n.90 <- nrow(data.90.percent)
y.90 <- data.90.percent$TOBS

# data plot
data.plot <- data %>% 
  ggplot(aes(x=DATE, y= TOBS))+
  geom_line() +
  xlab("Date")+
  ylab("Temperature")

# max date in 90%
max(data.90.percent$DATE)

# acf
acf.data <- acf(y.full,lag.max = length(y.full)-1,plot = FALSE)

acf.data.plot <- ggplot(mapping=aes(x=acf.data$lag)) +
  geom_segment(aes(y=acf.data$acf),yend=0,xend=acf.data$lag) +
  geom_hline(yintercept= qnorm(.975)/sqrt(acf.data$n.used),color="blue",
             linetype="dashed") +
  geom_hline(yintercept= qnorm(.025)/sqrt(acf.data$n.used),color="blue",
             linetype="dashed")+
  xlab("Lag")+
  ylab("ACF")



## ----------- stationarity

## lag 365 plots
lag.365.plot <- ggplot(mapping=aes(x=y.full[1:(n.full-365)],y=y.full[366:n.full]))+
  geom_point() +
  xlab("365 Day Lag Temp")+
  ylab("Temp") +
  labs(title="365 Day Lag Plot")


# trend
lm.fit<-lm(TOBS~sin(2*pi*omega) +
             cos(2*pi*omega)
           ,data=data.90.percent)
# trend's residuals
resid <- lm.fit$residual

resid.plot <- data.90.percent %>% 
  ggplot(aes(x=DATE, y= resid))+
  geom_line() +
  xlab("Date")+
  ylab("Residuals")

acf.resid <- acf(resid,lag.max = n.90 - 1, plot = FALSE)
acf.resid.rej <- (which(c(acf.resid$acf) < qnorm(0.025, 0, sqrt(1/n.90)) |
                    c(acf.resid$acf) > qnorm(0.975, 0, sqrt(1/n.90))) - 1)[-1]
max.rej.acf.resid <- max(acf.resid.rej)
max.rej.acf.resid

pacf.resid <- pacf(resid,lag.max = n.90 - 1, plot = FALSE)
pacf.resid.rej <- (which(c(pacf.resid$acf) < qnorm(0.025, 0, sqrt(1/n.90)) |
                     c(pacf.resid$acf) > qnorm(0.975, 0, sqrt(1/n.90))))
max.rej.pacf.resid <- max(pacf.resid.rej)
max.rej.pacf.resid
  
# level tests since data has a trend
# results in zero differencing so residuals from sinusoidal regression stationary

# ndiffs(resid, alpha = 0.05, test = "adf", type = "level") # not used
ndiffs(resid, alpha = 0.05, test = "pp", type = "level")
  

## ______________ what ARIMA model

# k step model selection from 90% of data

which.min(colMeans(k.step.arima.selection))

# The k=30 step ahead forecasts error minimizing model is 
# ARIMA(4,0,2) with a sinusoidal trend


## -------------------- model fitting ARIMA

xreg.90 <- model.matrix(~sin(2*pi*omega) +
               cos(2*pi*omega)
             ,data=data.90.percent)


arima.fit <- arima(y.90, order = c(4, 0, 2), 
                   xreg = xreg.90, method = "ML", include.mean = FALSE)


arima.fit.coef <- arima.fit$coef

## ---------------------- model prediction ARIMA

# predict entire 10%
predict.last.10.percent <- predict(arima(y.90, order = c(4, 0, 2), xreg = xreg.90,
                                         method = "ML", fixed = arima.fit.coef,
                                         include.mean = FALSE),
                                   n.ahead = n.full-n.90, 
                                   newxreg = model.matrix(~sin(2*pi*omega) +
                                                            cos(2*pi*omega)
                                                          ,data=data[(n.90+1):n.full,]))
### MSE last 10%
mean(as.vector((predict.last.10.percent$pred-y.full[n.full-n.90:n.full])^2))

#### MAE last 10%
mean(as.vector(abs(predict.last.10.percent$pred-y.full[n.full-n.90:n.full])))

#### Prediction Interval Coverage

# count number of times prediction in CI
PI_cov_10_percent_prediction <- 0
for (j in 1: length(predict.last.10.percent$pred)){
  index_in_data <- n.90+j
  tv.j <- y.full[index_in_data]
  pred.j <- predict.last.10.percent$pred[j]
  se.pred.j <- predict.last.10.percent$se[j]
  pred.j.lower.bound <- pred.j + qnorm(0.025)*se.pred.j
  pred.j.upper.bound <- pred.j + qnorm(0.975)*se.pred.j
  pred.j.ci <- ifelse((tv.j<=pred.j.upper.bound & tv.j>=pred.j.lower.bound),
                      1,0)
  
  PI_cov_10_percent_prediction <- PI_cov_10_percent_prediction + pred.j.ci
  
}

PI_cov_10_percent_prediction_final <- 
  (PI_cov_10_percent_prediction/length(predict.last.10.percent$pred))*100







# compare 30 day predictions and 1,7,14 day predictions for last 10% of data


indices.to.predict.from <- (n.90+1):(n.full-30) # gets indices in held out data minus last 30 days to predict from

prediction.arima <- as.data.frame(matrix(NA,nrow = 4*length(indices.to.predict.from),
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
  pred.1.ci <- ifelse((true.val.predict.1<=pred.1.upper.bound & 
                         true.val.predict.1>=pred.1.lower.bound),1,0)

  # predict 7
  y.predict.7 <- y.full[(index+1):(index+7)]
  data.predict.7 <- data[(index+1):(index+7),]
  xreg.predict.7 <- model.matrix(~sin(2*pi*omega) +
                                    cos(2*pi*omega)
                                  ,data=data.predict.7)
  ### predicts 7 steps ahead using coefficients from fit to first 90%
  predict.7 <- predict(arima(y.fitting, order = c(4, 0, 2), xreg = xreg.fit,
                              method = "ML", fixed = arima.fit.coef,
                              include.mean = FALSE),
                        n.ahead = 7, newxreg = xreg.predict.7)
  
  pred.7 <- as.numeric(predict.7$pred)
  se.pred.7 <- as.numeric(predict.7$se)
  true.val.predict.7 <- y.full[(index+1):(index+7)]
  MSE.predict.7 <- as.numeric(mean((pred.7-true.val.predict.7)^2))
  MAE.predict.7 <- as.numeric(mean(abs(pred.7-true.val.predict.7)))
  pred.7.ci <- 0
  # count number of times prediction in CI
  for (j in 1:7) {
    tv.j <- true.val.predict.7[j]
    pred.j <- pred.7[j]
    se.pred.j <- se.pred.7[j]
    pred.j.lower.bound <- pred.j + qnorm(0.025)*se.pred.j
    pred.j.upper.bound <- pred.j + qnorm(0.975)*se.pred.j
    pred.j.ci <- ifelse((tv.j<=pred.j.upper.bound & tv.j>=pred.j.lower.bound),
                        1,0)
    pred.7.ci <- pred.7.ci+pred.j.ci
  }
  pred.7.ci.pc <- pred.7.ci/7
  
  
  # predict 14
  y.predict.14 <- y.full[(index+1):(index+14)]
  data.predict.14 <- data[(index+1):(index+14),]
  xreg.predict.14 <- model.matrix(~sin(2*pi*omega) +
                                    cos(2*pi*omega)
                                  ,data=data.predict.14)
  ### predicts 14 steps ahead using coefficients from fit to first 90%
  predict.14 <- predict(arima(y.fitting, order = c(4, 0, 2), xreg = xreg.fit,
                              method = "ML", fixed = arima.fit.coef,
                              include.mean = FALSE),
                        n.ahead = 14, newxreg = xreg.predict.14)
  
  pred.14 <- as.numeric(predict.14$pred)
  se.pred.14 <- as.numeric(predict.14$se)
  true.val.predict.14 <- y.full[(index+1):(index+14)]
  MSE.predict.14 <- as.numeric(mean((pred.14-true.val.predict.14)^2))
  MAE.predict.14 <- as.numeric(mean(abs(pred.14-true.val.predict.14)))
  pred.14.ci <- 0
  # count number of times prediction in CI
  for (j in 1:14) {
    tv.j <- true.val.predict.14[j]
    pred.j <- pred.14[j]
    se.pred.j <- se.pred.14[j]
    pred.j.lower.bound <- pred.j + qnorm(0.025)*se.pred.j
    pred.j.upper.bound <- pred.j + qnorm(0.975)*se.pred.j
    pred.j.ci <- ifelse((tv.j<=pred.j.upper.bound & tv.j>=pred.j.lower.bound),
                        1,0)
    pred.14.ci <- pred.14.ci+pred.j.ci
  }
  pred.14.ci.pc <- pred.14.ci/14
  

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
  true.val.predict.30 <- y.full[(index+1):(index+30)]
  MSE.predict.30 <- as.numeric(mean((pred.30-true.val.predict.30)^2))
  MAE.predict.30 <- as.numeric(mean(abs(pred.30-true.val.predict.30)))
  pred.30.ci <- 0
  # count number of times prediction in CI
  for (j in 1:30) {
    tv.j <- true.val.predict.30[j]
    pred.j <- pred.30[j]
    se.pred.j <- se.pred.30[j]
    pred.j.lower.bound <- pred.j + qnorm(0.025)*se.pred.j
    pred.j.upper.bound <- pred.j + qnorm(0.975)*se.pred.j
    pred.j.ci <- ifelse((tv.j<=pred.j.upper.bound & tv.j>=pred.j.lower.bound),
                        1,0)
    pred.30.ci <- pred.30.ci+pred.j.ci
  }
  pred.30.ci.pc <- pred.30.ci/30
  
  prediction.arima[4*(i-1)+1,] <- c(1,MSE.predict.1,
                                    MAE.predict.1,pred.1.ci)
  prediction.arima[4*(i-1)+2,] <- c(7,MSE.predict.7,
                                    MAE.predict.7,pred.7.ci.pc)
  prediction.arima[4*(i-1)+3,] <- c(14,MSE.predict.14,
                                    MAE.predict.14,pred.14.ci.pc)
  prediction.arima[4*(i-1)+4,] <- c(30,MSE.predict.30,
                                    MAE.predict.30,pred.30.ci.pc)
  
}

prediction_arima_means <- prediction.arima %>% 
  group_by(horizon) %>% 
  summarise_all("mean") %>% 
  mutate(num_in_95CI=round(num_in_95CI*100,2)) %>% 
  rename("95% Pred. Int. Cov."=num_in_95CI)


prediction_arima_lb <-prediction.arima %>% 
  group_by(horizon) %>% 
  summarise_all(function(val) quantile(val,probs=0.025)) 

prediction_arima_ub <-prediction.arima %>% 
  group_by(horizon) %>% 
  summarise_all(function(val) quantile(val,probs=0.975)) 



### ------------------- state space


# selects covariate in observation equation
which.min(colMeans(k.step.state.space.selection))

y.10 <- y.full[(n.90+1):n.full]

full_state_space_predictions <- as.vector(fit.state.full$ytT)[(n.90+1):n.full]

full_state_space_predictions_se <- as.vector(fit.state.full$ytT.se)[(n.90+1):n.full]

mse.state.space <- mean((y.10-full_state_space_predictions)^2)
mae.state.space <- mean(abs(y.10-full_state_space_predictions))

# Prediction Interval Coverage

PI_cov_10_percent_prediction_state_space <- 0
for (j in 1: length(y.10)){
  tv.j <- y.10[j]
  pred.j <- full_state_space_predictions[j]
  se.pred.j <- full_state_space_predictions_se[j]
  pred.j.lower.bound <- pred.j + qnorm(0.025)*se.pred.j
  pred.j.upper.bound <- pred.j + qnorm(0.975)*se.pred.j
  pred.j.ci <- ifelse((tv.j<=pred.j.upper.bound & tv.j>=pred.j.lower.bound),
                      1,0)
  
  PI_cov_10_percent_prediction_state_space <- 
    PI_cov_10_percent_prediction_state_space + pred.j.ci
  
}

pic_state_space <- (PI_cov_10_percent_prediction_state_space/length(y.10))*100


plot.state.space.predictions <- data[(n.90+1):n.full,] %>% 
  ggplot(aes(x=DATE,y=TOBS)) +
  geom_line(aes())+
  geom_line(aes(y=full_state_space_predictions, linetype="prediction"),
            color="blue") +
  geom_line(aes(y=full_state_space_predictions +
                  qnorm(0.025)*full_state_space_predictions_se,
                linetype="Confidence Interval"),color="blue") +
  geom_line(aes(y=full_state_space_predictions +
                  qnorm(0.975)*full_state_space_predictions_se,
                linetype="Confidence Interval"),color="blue") +
  scale_linetype_manual(values=c("dashed","solid"),name=NULL,
                        labels=c("Confidence\n Interval","Forecast"))+
  xlab("Date")+
  ylab("Temperature")






