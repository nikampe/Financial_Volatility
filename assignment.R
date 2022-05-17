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
library(nortsTest)
library(stargazer)
library(forecast)

# Working directory
getwd()
#setwd("~/Documents/GitHub/Financial_Volatility")

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

#THEORY addition in paper: (describe data set)

# Plot of time series of Clsoing prices
ggplot(xtrackers_msci, aes(x=Date, y=Close)) +
  geom_line() + 
  labs(title = "Price development", 
       subtitle = "Xtrackers MSCI World UCITS ETF | 08/2014 - 05/2022", 
       x = "Time",
       y = "Price (in €)") +
  theme(panel.grid.minor = element_blank())


# Plot of time series of log returns
ggplot(xtrackers_msci, aes(x=Date, y=log_returns)) +
  geom_line() + 
  labs(title = "Log Returns", 
       subtitle = "Xtrackers MSCI World UCITS ETF | 08/2014 - 05/2022", 
       x = "Time",
       y = "Log Return") +
  theme(panel.grid.minor = element_blank())

#THEORY addition in paper: (Take reference to high volatility periods à Dot-Com, Global Financial Crisis, Covid-19 …)

# Description of dataset
basicStats(xtrack_returns) #potentially write about observed patterns in paper

# 3 Model Fit and Forecast Comparisons ------------------------------------
# 3.1 Check for stylized facts of the time series -------------------------

#autocorrelation in time series for different lags
# Ljung-Box Test
# -> p-value < 0.01 (or 0.05, 0.1) rejects null hypothesis of autocorrelation for given lag
for (i in c(1,3,5,7)) {
  stat <- Box.test(xtrack_returns, lag=i, type="Ljung")
  print(stat$p.value)
}

#INSIGHT: We have significant autocorrelations for the lags tested


# Lagrange multiplier test
# -> p-value < 0.01 (or 0.05, 0.1) rejects null hypothesis of homoskedasticity for a maximum lag of 5
Lm.test(xtrack_returns,lag.max = 5,alpha = 0.05)

#INSIGHT: our log returns are, as expected, heteroskedastic for a lag of up to 5 time steps

#TBD Asymmetries and tail tests TBD

#Appendix: Additional tests

# Test for zero mean
t.test(xtrack_returns)
# Test for normality (Jarque-Bera Test)
normalTest(xtrack_returns, method = "jb")

#confirm our results from above


# 3.2 Identification of GARCH models --------------------------------------

# Comparison of log and squared log returns
par(mfrow=c(2,1))
acf(xtrack_returns, main="Xtrackers MSCI World UCITS ETF Log Returns")
acf(xtrack_returns^2, main="Xtrackers MSCI World UCITS ETF Squared Log Returns")

# ACF and PACF of log and squared log returns
# -> (Partial) autocorrelations outside the confidence bounds are significant
par(mfrow=c(2,2))
acf(xtrack_returns, lag.max=30, main = "ACF | Log Returns", ylim=range(-1,1))
pacf(xtrack_returns, lag.max=30, main = "PACF | Log Returns", ylim=range(-1,1))
acf(xtrack_returns^2, lag.max=30, main = "ACF | Squared Log Returns", ylim=range(-1,1))
pacf(xtrack_returns^2, lag.max=30, main = "PACF | Squared Log Returns", ylim=range(-1,1))

#THEORY: Infer from the graph the signficant lags

#Tables with sample autocorrelations up to lag 30
acf(xtrack_returns, lag.max=30, pl=F)
pacf(xtrack_returns, lag.max=30, pl=F)
acf(xtrack_returns^2, lag.max=30, pl=F)
pacf(xtrack_returns^2, lag.max=30, pl=F) #potentially add the significance level and enhance the looks of the tables, still missing here

#INSIGHT: We try a GARCH(1,1) model (benchmark) against other GARCH models, particularly a GARCH(1,2),
          #a GARCH(2,2), as well as a GARCH(3,3) due to the partly strongly auto correlated squared log returns.

### further insight on the special GARCHs to be added! ###


# 3.3 Model Fit of GARCH --------------------------------------------------

#Split return series in training and test sample
xtrackers_msci$Split <- rep(x = c("Training", "Test"),
                        times = c(floor(x = 0.7 * nrow(x = xtrackers_msci)), ceiling(x = 0.3 * nrow(x = xtrackers_msci))))
train <- xtrackers_msci[xtrackers_msci$Split %in% c('Training'),]
test <- xtrackers_msci[xtrackers_msci$Split %in% c('Test'),]

#INSIGHT: this split is interesting since we are essentially going to predict the Corona crisis.

#make a plot of the lof returns over the sample period (pre-Corona)
ggplot(train, aes(x=Date, y=log_returns)) +
  geom_line() + 
  labs(title = "Training sample log returns", 
       subtitle = "Xtrackers MSCI World UCITS ETF | 08/2014 - 01/2020", 
       x = "Time",
       y = "log returns") +
  theme(panel.grid.minor = element_blank())

train_returns <- as.vector(train$log_returns)

#estimate model parameters for General GARCH models
#ARCH(1)
m0 <- garch(train_returns, order = c(0, 1))
summary(m0)
plot(m0)
fit0 <- fitted.values(m0)

#GARCH(1,1)
m1 <- garch(train_returns, order = c(1, 1))
summary(m1)
plot(m1)
fit1 <- fitted.values(m1)

#GARCH(1,2)
m2 <- garch(train_returns, order = c(1, 2))
summary(m2)
plot(m2)
fit2 <- fitted.values(m2)

#GARCH(2,2)
m3 <- garch(train_returns, order = c(2, 2))
summary(m3)
plot(m3)
fit3 <- fitted.values(m3)

#GARCH(3,3)
m4 <- garch(train_returns, order = c(3, 3))
summary(m4)
plot(m4)
fit4 <- fitted.values(m4)

### SPECIAL GARCHS STILL MISSING HERE ### STILL TO BE INCLUDED @ NIKLAS


### EVALUATION WITH AIC/BIC STILL MISSING HERE ### @NIKLAS/CLARA FEEL FREE TO COCLUDE THIS HERE


### plot fitted values against sample values ### TBD



# 4 Model Forecast of GARCH -----------------------------------------------

o Split of training (fit) and testing (forecast) samples à Use/Plot testing sample
o Define forecast horizon/lags n
o Calculate volatility forecasts of General GARCH and best Special GARCH

# 5 Model Forecast Comparison with Realized Volatility (RV) ---------------

- Model Forecast Comparison with Realized Volatility (RV)
o Choose forecast horizon/lags n as before
o Calculate volatility forecasts by Realized Volatility (RV)
o Compare RV forecast results with General GARCH and Special GARCH forecasts
(…) = only referring to written part, not to code

# 4 Appendix --------------------------------------------------------------

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
