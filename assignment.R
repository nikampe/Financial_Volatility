#A uthors: Basadonna, Gutjahr, Kampe
# ----------------------------------------
# Financial Volatility - Group Assignment
# University of St. Gallen | Spring Semester 2022
# Clara Elisa Basadonna | 21-607-395
# Jan Gutjahr | 21-607-544
# Niklas Leander Kampe | 16-611-618
# ----------------------------------------

# 1 Libraries & Data Preparation ----------------------------------------

# Load packages
library(fBasics)
library(zoo)
library(fGarch)
library(tseries)
library(car)
library(systemfit)
library(mvtnorm)
library(quadprog)
library(VGAM)
library(sandwich)
library(magrittr)
library(tidyverse)
library(ggplot2)

# Working directory
getwd()
setwd("~/Documents/GitHub/Financial_Volatility")

# Import data set and numeric transformation
xtrackers_msci <- read.csv(as.matrix("XDWD.DE.csv")) %>%
  mutate_at(vars(Date),as.Date) %>%
  mutate_at(vars(Open, High, Low, Close, Adj.Close, Volume),as.numeric)

# Check of NA values
xtrackers_msci[is.na(xtrackers_msci[,2]),]

# Check column classes and omit NA values
str(xtrackers_msci)
xtrackers_msci %<>% na.omit()

# Calculate log returns from closing prices
xtrackers_msci[,"log_returns"] <- 0
for (i in 2:nrow(xtrackers_msci)) {
  xtrackers_msci[i,"log_returns"] <- (log(xtrackers_msci[i,5]/xtrackers_msci[i-1,5]))*100
}
xtrackers_msci <- xtrackers_msci[-1,]

# Save log returns as vector
xtrack_returns <- xtrackers_msci[,8]

# 2 Descriptive Statistics --------------------------------------------------------------

# Summary statistics
basicStats(xtrack_returns)

# Plot of log returns
ggplot(xtrackers_msci, aes(x=Date, y=log_returns)) +
  geom_line() + 
  labs(title = "Log Returns", 
       subtitle = "Xtrackers MSCI World UCITS ETF | 08/2014 - 05/2022", 
       x = "Time",
       y = "Log Return") +
  theme(panel.grid.minor = element_blank())

# 3 Analysis --------------------------------------------------------------

# Test for zero mean
t.test(xtrack_returns)

# Test for normality (Jarque-Bera Test)
normalTest(xtrack_returns, method = "jb")

# ACF and PACF of log returns
# -> (Partial) Autocorrelations ouside the confidence bounds are significant
par(mfrow=c(1,2))
acf(xtrack_returns, lag.max=30, main = "ACF | Log Returns", ylim=range(-1,1))
pacf(xtrack_returns, lag.max=30, main = "PACF | Log Returns", ylim=range(-1,1))

# Ljung-Box Test
# -> according to ACF/PACF plot: Test autocorrelation for lags 3 and 7
# -> p-value < 0.01 (or 0.05, 0.1) rejects null hypothesis of autocorrelation for given lag
for (i in c(1,3,5,7)) {
  stat <- Box.test(xtrack_returns, lag=i, type="Ljung")
  print(stat$p.value)
}

# First and Second differences of log returns
par(mfrow=c(2,1))
xtrack_returns_fd <- ts(xtrackers_msci[2:1954,8]-xtrackers_msci[1:1953,8])
ts.plot(xtrack_returns_fd, main="1st Differences | Log Returns", ylab = "1st Difference")
xtrack_returns_sd <- ts(xtrackers_msci[3:1954,8]-2*xtrackers_msci[2:1953,8] + xtrackers_msci[1:1952,8])
ts.plot(xtrack_returns_sd, main="2nd Differences | Log Returns", ylab = "2nd Difference")

# ACF and PACF of First and Second differences of log returns
par(mfrow=c(2,2))
acf(xtrack_returns_fd, lag.max=30, main="ACF | 1st Differences | Log Returns")
pacf(xtrack_returns_fd, lag.max=30, main="PACF | 1st Differences | Log Returns")
acf(xtrack_returns_sd, lag.max=30, main="ACF | 2nd Differences | Log Returns")
pacf(xtrack_returns_sd, lag.max=30, main="PACF | 2nd Differences | Log Returns")

# Comparison of simple and squared returns
par(mfrow=c(2,1))
acf(xtrack_returns, main="Xtrackers MSCI World UCITS ETF Log Returns")
acf(xtrack_returns^2, main="Xtrackers MSCI World UCITS ETF Squared Returns")

#compare returns with normal dist
a=density(xtrack_returns)

max(xtrack_returns)
min(xtrack_returns) #here, we check which limits to use in our following plot

plot(a,main="Daily Xtrackers MSCI World UCITS ETF returns kernel density",xlim=c(-11,8))
lines(a$x,dnorm(a$x,mean=mean(xtrack_returns),sd=sqrt(var(xtrack_returns))),lty=2)

#fit the first garch on our data
m0 <- garchFit(xtrack_returns ~ garch(1,0), data = xtrack_returns,trace=F)
summary(m0)
predict(m0,5)

#look at some comparative statistics again for our time series
par(mfrow=c(2,2))
qqplot(qt(ppoints(length(m0@residuals)), df = 6),y=m0@residuals/m0@sigma.t ,main = expression("Q-Q plot for" ~~ {t}[nu == 6]),xlab="qt",ylab="residuals")
acf(m0@residuals/m0@sigma.t,main="ARCH(1) residuals")
acf(m0@residuals^2/m0@sigma.t^2,main="ARCH(1) squared residuals")
acf(abs(m0@residuals/m0@sigma.t),main="ARCH(1) absolute residuals")







# 4 Appendix --------------------------------------------------------------

#nice snippet from Audrino that shows the relation between the models

par(mfrow=c(2,2))
set.seed(66)      # so you can reproduce these results
a = arima.sim(list(order=c(1,0,0), ar=.9), n=100) #AR(1)  
acf(a,lag.max=30,main="AR(1)")
a = arima.sim(list(order=c(1,0,0), ar=-0.8), n=100) #AR(1)  
acf(a,lag.max=30,main="AR(1)")
a = arima.sim(list(order=c(2,0,0), ar=c(1.2,-0.35)), n=100) #AR(2)
acf(a,lag.max=30,main="AR(2)")
a = arima.sim(list(order=c(2,0,0), ar=c(-0.2,0.35)), n=100) #AR(2)
acf(a,lag.max=30,main="AR(2)")

par(mfrow=c(2,2))
set.seed(148)
a = arima.sim(list(order=c(0,0,1), ma=-0.7), n=100) #MA(1)
acf(a,lag.max=30,main="MA(1)")
a = arima.sim(list(order=c(0,0,2), ma=c(-0.7,0.7)), n=100) #MA(2)
acf(a,lag.max=30,main="MA(2)")
a = arima.sim(list(order=c(1,0,1), ma=0.1,ar=0.9), n=100) #ARMA(1,1)
acf(a,lag.max=30,main="ARMA(1,1)")
a = arima.sim(list(order=c(1,0,1), ma=-0.77,ar=0.55), n=100) #ARMA(1,1)
acf(a,lag.max=30,main="ARMA(1,1)")