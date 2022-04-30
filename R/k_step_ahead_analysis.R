library(tidyverse)
library(readr)
library(forecast)

load("Data/data.90.percent.Rdata")

# this is analysis on first 90% of data

# fits regression with sinusoidal terms
lm.fit<-lm(TOBS~sin(2*pi*omega) +
             cos(2*pi*omega)
           ,data=data.90.percent)

# regression residuals
resid <- lm.fit$residuals

# plot of residuals looks stationary
resid.plot <- data.90.percent %>% 
  ggplot(mapping=aes(x=DATE,y=resid))+
  geom_line()+
  xlab("Date")+
  ylab("Residuals")

# level tests since data has a trend
# results in zero differencing so residuals from sinusoidal regression stationary

ndiffs(resid, alpha = 0.05, test = "adf", type = "level")
ndiffs(resid, alpha = 0.05, test = "pp", type = "level")

# max p and q terms and n subsets
pmax <- 20
qmax <- 5
nsubset <- 10


k.step.arima.selection <- matrix(NA,nrow = nsubset, ncol =  (pmax+1)*(qmax+1))

# get columnames
p.val <- rep(0:pmax,each=(qmax+1))
q.val <- rep(0:qmax,(pmax+1))

c.names <- paste0("p=",p.val,", q=",q.val)

colnames(k.step.arima.selection) <- c.names

for (i in 1:nsubset) {
  
  data.subset <- resid[(((i-1)*390 + 1):((i)*390))]
  data.train <- data.subset[1:360]
  data.test <- data.subset[-(1:360)]
  
  for (p in 0:10) {
    
    for (q in 0:qmax) {
      
      train.arima.fit <- arima(data.train, order = c(p, 0, q), 
                               method = "ML")
      predict.arima.fit <- predict(train.arima.fit, n.ahead = 30)
      
      mse.forecast <- mean((data.test-predict.arima.fit$pred)^2)
      
      k.step.arima.selection[i,p*(qmax+1)+q+1] <- mse.forecast
      
    }
  }
}

save(k.step.arima.selection,file="Data/k.step.arima.selection.Rdata")
