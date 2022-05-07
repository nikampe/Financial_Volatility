#Authors: Basadonna, Gutjahr, Kampe

# 1 Libraries and Data Preparation ----------------------------------------

#check what the wd is set to (should be the directory where the data is saved)
getwd()

#load packages needed for analysis
library(fBasics) #Load the package fBasics
library(zoo)
library(fGarch) #Load the package for ARCH/GARCH estimation
library(tseries)
library(car) #consistent standard errors
library(systemfit)
library(mvtnorm)
library(quadprog)
library(VGAM)
library(sandwich)
library(magrittr)

#load data and transform to numeric
xtrackers_msci <- read.csv(as.matrix("XDWD.DE.csv")) %>%
  mutate_at(vars(Date),as.Date) %>%
  mutate_at(vars(Open, High, Low, Close, Adj.Close, Volume),as.numeric)

xtrackers_msci[is.na(xtrackers_msci[,2]),] #I'm not exactly sure what the issue with these two dates is, let's think about this

#check the data quality and eliminate the NAs
str(xtrackers_msci)
xtrackers_msci %<>% na.omit() #omit the two NA days where all our values are NA

#transform the closing prices to obtain log returns
xtrackers_msci[,8] <- 0
for (i in 2:nrow(xtrackers_msci)) {
  
  xtrackers_msci[i,8] <- (log(xtrackers_msci[i,5]/xtrackers_msci[i-1,5]))*100
  
}

colnames(xtrackers_msci)[8] <- "log_returns"
xtrackers_msci <- xtrackers_msci[-1,] #get rid of first observation that has no return value

#take out the returns from the dataset
xtrack_returns <- xtrackers_msci[,8]


# 2 Analysis --------------------------------------------------------------

#obtain basic stats
basicStats(xtrack_returns)

#Test for mean being zero
t.test(xtrack_returns)

#Test for normality of returns
normalTest(xtrack_returns,method="jb") #(Jarque-Bera test)

#Test for normality via Chi-Square
par(mfrow=c(1,2))
acf(xtrack_returns,lag.max=30,main="Xtrackers MSCI World UCITS ETF")

Box.test(xtrack_returns,lag=30,type="Ljung") #Ljung-Box statistic Q(5)

pacf(xtrack_returns, main="Xtrackers MSCI World UCITS ETF")

#find first differences
par(mfrow=c(2,1))
xtrack_returns_fd <- ts(xtrackers_msci[2:1954,8]-xtrackers_msci[1:1953,8])
ts.plot(xtrack_returns_fd, main="Xtrackers MSCI World UCITS ETF First Differences", ylab = "First Differences")

xtrack_returns_sd <- ts(xtrackers_msci[3:1954,8]-2*xtrackers_msci[2:1953,8] + xtrackers_msci[1:1952,8])
ts.plot(xtrack_returns_sd, main="Xtrackers MSCI World UCITS ETF Second Differences", ylab = "Second Differences")

par(mfrow=c(2,2))
acf(xtrack_returns_fd, main="Xtrackers MSCI World UCITS ETF First Differences")
pacf(xtrack_returns_fd, main="Xtrackers MSCI World UCITS ETF First Differences")

acf(xtrack_returns_sd, main="Xtrackers MSCI World UCITS ETF Second Differences")
pacf(xtrack_returns_sd, main="Xtrackers MSCI World UCITS ETF Second Differences")

#compare simple and squared returns
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







# 3 Appendix --------------------------------------------------------------

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