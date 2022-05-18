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
library(moments)
library(modelsummary)
library(xtable)
options("modelsummary_format_numeric_latex" = "plain")

# Working directory
getwd()
#setwd("~/Documents/GitHub/Financial_Volatility") #NIK
setwd("C:/Users/jan_g/OneDrive/HSG/2nd semester/Financial volatility/Assignment") #JAN

# 1.0 Data Import & Cleaning --------------------------------------------------------------

# 1.1 Import data set & numeric transformation -----------------------------------
xtrackers_msci <- read.csv(as.matrix("XDWD.DE.csv")) %>%
  mutate_at(vars(Date),as.Date) %>%
  mutate_at(vars(Open, High, Low, Close, Adj.Close, Volume),as.numeric)

# 1.2 Check for  NA values -------------------------------------------------------
xtrackers_msci[is.na(xtrackers_msci[,2]),]

# 1.3 Check column classes & omit NA values --------------------------------------
str(xtrackers_msci)
xtrackers_msci %<>% na.omit()

# 1.4 Calculate log returns from closing prices ----------------------------------
xtrackers_msci[,"log_returns"] <- 0
for (i in 2:nrow(xtrackers_msci)) {
  xtrackers_msci[i,"log_returns"] <- (log(xtrackers_msci[i,5]/xtrackers_msci[i-1,5]))*100
}
xtrackers_msci <- xtrackers_msci[-1,]

# 1.5 Save log returns as vector -------------------------------------------------
xtrack_returns <- xtrackers_msci[,8]

xtrack_returns
# 2 Descriptive Statistics --------------------------------------------------------------

# 2.1 Plot of time series of Closing prices and more -----------------------------
pdf("Figures/Plot_Prices.pdf") 
ggplot(xtrackers_msci, aes(x=Date, y=Close)) +
  geom_line() + 
  labs(title = "Prices", 
       subtitle = "Xtrackers MSCI World UCITS ETF | 08/2014 - 05/2022", 
       x = "Time",
       y = "Prices (in EUR)") +
  theme(panel.grid.minor = element_blank())
dev.off()

# 2.2 Plot of time series of log returns -----------------------------------------
pdf("Figures/Plot_Returns.pdf") 
ggplot(xtrackers_msci, aes(x=Date, y=log_returns)) +
  geom_line() + 
  labs(title = "Log Returns", 
       subtitle = "Xtrackers MSCI World UCITS ETF | 08/2014 - 05/2022", 
       x = "Time",
       y = "Log Returns") +
  theme(panel.grid.minor = element_blank())
dev.off()

# 2.3 Summary Statistics of dataset ----------------------------------------------
summ_stats <- t(as.data.frame(round(basicStats(xtrack_returns),4)))
rownames(summ_stats)[1] = "Log returns"

#display as latex table
sink(file = "Latex/summary_stats.txt")
xtable(summ_stats)
sink(file = NULL)

# 3 Model Fit and Forecast Comparisons -------------------------------------------

# 3.1 Check for stylized facts of the time series --------------------------------

# 3.1.1 Check for autocorrelation for different lags: Ljung-Box Test -------------
for (i in c(1,3,5,7)) {
  stat <- Box.test(xtrack_returns, lag=i, type="Ljung") 
  print(stat$p.value) # -> p-value < 0.01 (or 0.05, 0.1) rejects null hypothesis of autocorrelation for given lag
} # --> INSIGHT: We have significant autocorrelations for the lags tested

# 3.1.2 Check for (G)ARCH effects: LM (Lagrange Multiplier) Test -----------------
Lm.test(xtrack_returns, lag.max = 5, alpha = 0.05) # -> p-value < 0.01 (or 0.05, 0.1) rejects null hypothesis of homoskedasticity for a maximum lag of 5
# --> INSIGHT: our log returns are, as expected, heteroskedastic for a lag of up to 5 time steps

# 3.1.3 Check for asymmetries and tail distributions: Histogram & Statistics -----
binwidth <- 0.2
pdf("Figures/Histogram_Returns.pdf") 
ggplot(xtrackers_msci, aes(x=log_returns)) +
  geom_histogram(aes(y = ..density..), binwidth = binwidth) + 
  stat_function(fun = dnorm,
                args = list(mean = mean(xtrackers_msci$log_returns, na.rm=TRUE),
                            sd = sd(xtrackers_msci$log_returns, na.rm=TRUE)),
                col = "red",
                size = 0.5) + 
  labs(title = "Log Returns", 
       subtitle = "Xtrackers MSCI World UCITS ETF | 08/2014 - 05/2022", 
       x = "Log Returns",
       y = "Density") +
  theme(panel.grid.minor = element_blank()) # --> INSIGHT: Approximately zero mean, Fat tails --> Further diagnositcs below
dev.off()

# 3.1.4 Test for additional distribution properties: Moments, T-Test, JB Test ---- 
t.test(xtrack_returns) # --> INSIGHT: Zero mean with 95% confidence
normalTest(xtrack_returns, method = "jb") # --> INSIGHT: Not normaly distirbuted
skewness(xtrackers_msci$log_returns, na.rm=TRUE) # --> INSIGHT: Slightly negatively skewed
kurtosis(xtrackers_msci$log_returns, na.rm=TRUE) # --> INSIGHT: Fat/Heavy tails (interesting, note to self)

# 3.1.5 Check for asymmetries & leverage effect: ---------------------------------

# 3.1.5.1 Leverage Effect: Autocorrelations --------------------------------------
a <- c()
c <- c()
for (i in 1:length(xtrack_returns))  {
  a <- c(a, max(xtrack_returns[i], 0))
}
for (h in 1:40) {
  c <- c(c, cor(a[(1+h):(length(a))],xtrack_returns[(1):(length(xtrack_returns)-h)]))
}
cor_summary <- data.frame(row.names = c("XDWD"))
lags <- c(1,2,5,10,20,40)
for (i in 1:length(lags)) {
  j = lags[i]
  cor_summary[1,i] <- c[j]
}
colnames(cor_summary) <- lags
sink(file = "Latex/Table_Leverage_Effect_autocorr.txt")
xtable(cor_summary) # --> INISGHT: Leverage effect given (negative autocorrelations)
sink(file = NULL)

## 3.1.5.2 Asymmetry: Test Regression --------------------------------------------
h_arr <- c(1,5)
sign_bias_coefs <- c()
neg_size_bias_coefs <- c()
pos_sign_bias_coefs <- c()
asy_summary <- data.frame(row.names = c("XDWD"))
for (h in h_arr) {
  sign_bias <- (xtrack_returns < 0)
  neg_size_bias <- c()
  pos_sign_bias <- c()
  for (i in 1:length(xtrack_returns)) {
    neg_size_bias <- c(neg_size_bias, 100*min(xtrack_returns[i], 0))
    pos_sign_bias <- c(pos_sign_bias, 100*max(xtrack_returns[i], 0))
  }
  sign_bias <- lm(100*xtrack_returns[(h+1):length(xtrack_returns)]^2 ~ sign_bias[1:(length(xtrack_returns)-h)])
  sign_bias_coefs <- c(sign_bias_coefs, sign_bias$coef[2])
  neg_size_bias <- lm(100*xtrack_returns[(h+1):length(xtrack_returns)]^2 ~ neg_size_bias[1:(length(xtrack_returns)-h)])
  neg_size_bias_coefs <- c(neg_size_bias_coefs, neg_size_bias$coef[2])
  pos_sign_bias <- lm(100*xtrack_returns[(h+1):length(xtrack_returns)]^2 ~ pos_sign_bias[1:(length(xtrack_returns)-h)])
  pos_sign_bias_coefs <- c(pos_sign_bias_coefs, pos_sign_bias$coef[2])
}
coefs <- c(sign_bias_coefs, neg_size_bias_coefs, pos_sign_bias_coefs)
for (i in 1:6) {
  coef <- coefs[i]
  asy_summary[1,i] <- coef
}
colnames(asy_summary) <- c("SB (h=1)", "SB (h=5)", "NSB (h=1)", "NSB (h=5)", "PSB (h=1)", "PSB (h=5)")
sink(file = "Latex/Table_Leverage_Effect_regression.txt")
xtable(asy_summary) # --> INISGHT: Asymmetry given (stronger effect of negativity)
sink(file = NULL)

# 3.2 Identification of GARCH models ---------------------------------------------

# 3.2.1 Comparison of log and squared log returns --------------------------------
max_lags <- 30
pdf("Figures/Comparison_Log_Returns.pdf") 
par(mfrow=c(2,1))
acf(xtrack_returns, lag.max=max_lags, main="Xtrackers MSCI World UCITS ETF Log Returns")
acf(xtrack_returns^2, lag.max=max_lags, main="Xtrackers MSCI World UCITS ETF Squared Log Returns")
dev.off()

# 3.2.2 ACF and PACF of log and squared log returns ----------------------------
pdf("Figures/Comparison_Log_Returns.pdf") 
par(mfrow=c(2,2))
acf(xtrack_returns, lag.max=max_lags, main="ACF | Log Returns", ylim=range(-1,1))
pacf(xtrack_returns, lag.max=max_lags, main="PACF | Log Returns", ylim=range(-1,1))
acf(xtrack_returns^2, lag.max=max_lags, main="ACF | Squared Log Returns", ylim=range(-1,1))
pacf(xtrack_returns^2, lag.max=max_lags, main="PACF | Squared Log Returns", ylim=range(-1,1)) # -> (Partial) autocorrelations outside the confidence bounds are significant
dev.off()

# 3.2.3 Table with sample autocorrelations -------------------------------------
lags <- 0:30
acf_returns <- round(acf(xtrack_returns, lag.max=max_lags, pl=F)$acf, 4)
pacf_returns <- c("-", round(pacf(xtrack_returns, lag.max=max_lags, pl=F)$acf, 4))
acf_squared_returns <- round(acf(xtrack_returns^2, lag.max=max_lags, pl=F)$acf, 4)
pacf_squared_returns <- c("-", round(pacf(xtrack_returns^2, lag.max=max_lags, pl=F)$acf, 4))
acf_table <- data.frame(ACF_Returns = acf_returns, 
                        PACF_Returns = pacf_returns, 
                        ACF_Squared_Returns = acf_squared_returns, 
                        PACF_Squared_Returns = pacf_squared_returns, 
                        row.names = lags) # --> INSIGHT: We try a GARCH(1,1) model (benchmark) against other GARCH models, particularly a GARCH(1,2),
                                          #     a GARCH(2,2), as well as a GARCH(3,3) due to the partly strongly auto correlated squared log returns.

# 3.3 Model Fit of GARCH -------------------------------------------------------

# 3.3.1 Train-Test Split -------------------------------------------------------
train_size <- 0.7
xtrackers_msci$Split <- rep(x = c("Training", "Test"),
                        times = c(floor(x = train_size * nrow(x = xtrackers_msci)), ceiling(x = (1-train_size) * nrow(x = xtrackers_msci))))
train <- xtrackers_msci[xtrackers_msci$Split %in% c('Training'),]
test <- xtrackers_msci[xtrackers_msci$Split %in% c('Test'),] # --> INSIGHT: this split is interesting since we are essentially going to predict the Corona crisis.

# 3.3.2 Plot of training sample (pre-Covid) ------------------------------------
pdf("Figures/Plot_Training_Sample.pdf") 
ggplot(train, aes(x=Date, y=log_returns)) +
  geom_line() + 
  labs(title = "Training Sample - Log Returns", 
       subtitle = "Xtrackers MSCI World UCITS ETF | 08/2014 - 01/2020", 
       x = "Time",
       y = "log returns") +
  theme(panel.grid.minor = element_blank())
dev.off()
train_returns <- as.vector(train$log_returns)
test_returns <- as.vector(test$log_returns)

# 3.3.3 Model fit & Determination of model parameters --------------------------

#we start by building simple GARCH models

## ARCH(1)
m0 <- garchFit(formula = ~ garch(1,0), data=train_returns, trace=F, description="ARCH(1)") #this is the ARCH model (this is basically a pure auto-regression)
fit0 <- fitted(m0)
summary(m0)

pdf("Figures/Analysis_Residuals_0.pdf")
layout(matrix(c(1,1,2,3,4,5), nrow=3, ncol=2, byrow = TRUE))
plot(train$Date, train$log_returns, type="l", col="blue", main="ARCH(1) | Fitted Values", xlab="Time", ylab="Return")
lines(train$Date, fit0,  col=alpha("red", 0.5))
plot(train$Date, residuals(m0), type="l", main="Residuals", xlab="Time", ylab="Residual")
hist(residuals(m0), main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
qqnorm(residuals(m0), main="QQ-Plot | Residuals")
acf(residuals(m0)[-1]^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-1,1))
dev.off()

## GARCH(1,1)
m1 <- garchFit(formula = ~ garch(1,1), data=train_returns, trace=F, description="GARCH(1,1)")
fit1 <- fitted(m1)
summary(m1)

pdf("Figures/Analysis_Residuals_1.pdf") 
layout(matrix(c(1,1,2,3,4,5), nrow=3, ncol=2, byrow = TRUE))
plot(train$Date, train$log_returns, type="l", col="blue", main="GARCH(1,1) | Fitted Values", xlab="Time", ylab="Return")
lines(train$Date, fit1,  col=alpha("red", 0.5))
plot(train$Date, residuals(m1), type="l", main="Residuals", xlab="Time", ylab="Residual")
hist(residuals(m1), main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
qqnorm(residuals(m1), main="QQ-Plot | Residuals")
acf(residuals(m1)[-1]^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-1,1))
dev.off()

## GARCH(1,2)
m2 <- garchFit(formula = ~ garch(1,2), data=train_returns, trace=F, description="GARCH(1,2)")
fit2 <- fitted(m2)
summary(m2)

pdf("Figures/Analysis_Residuals_2.pdf")
layout(matrix(c(1,1,2,3,4,5), nrow=3, ncol=2, byrow = TRUE))
plot(train$Date,train$log_returns, type="l", col="blue", main="GARCH(1,2) | Fitted Values", xlab="Time", ylab="Return")
lines(train$Date, fit2,  col=alpha("red", 0.5))
plot(train$Date, residuals(m2), type="l", main="Residuals", xlab="Time", ylab="Residual")
hist(residuals(m2), main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
qqnorm(residuals(m2), main="QQ-Plot | Residuals")
acf(residuals(m2)[-1][-1]^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-1,1))
dev.off()

## GARCH(2,2)
m3 <- garchFit(formula = ~ garch(2,2), data=train_returns, trace=F, description="GARCH(2,2)")
fit3 <- fitted(m3)
summary(m3)

pdf("Figures/Analysis_Residuals_3.pdf") 
layout(matrix(c(1,1,2,3,4,5), nrow=3, ncol=2, byrow = TRUE))
plot(train$Date,train$log_returns, type="l", col="blue", main="GARCH(2,2) | Fitted Values", xlab="Time", ylab="Return")
lines(train$Date, fit3,  col=alpha("red", 0.5))
plot(train$Date, residuals(m3), type="l", main="Residuals", xlab="Time", ylab="Residual")
hist(residuals(m3), main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
qqnorm(residuals(m3), main="QQ-Plot | Residuals")
acf(residuals(m3)[-1][-1]^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-1,1))
dev.off()

## GARCH(3,3)
m4 <- garchFit(formula = ~ garch(3,3), data=train_returns, trace=F, description="GARCH(3,3)")
fit4 <- fitted(m4)
summary(m4)

pdf("Figures/Analysis_Residuals_4.pdf") 
layout(matrix(c(1,1,2,3,4,5), nrow=3, ncol=2, byrow = TRUE))
plot(train$Date,train$log_returns, type="l", col="blue", main="GARCH(3,3) | Fitted Values", xlab="Time", ylab="Return")
lines(train$Date, fit4,  col=alpha("red", 0.5))
plot(train$Date, residuals(m4), type="l", main="Residuals", xlab="Time", ylab="Residual")
hist(residuals(m4), main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
qqnorm(residuals(m4), main="QQ-Plot | Residuals")
acf(residuals(m4)[-1][-1][-1]^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-1,1))
dev.off()

## Model Summary
summary_standard <- list("GARCH(1,1)" = m1, "GARCH(1,2)" = m2, "GARCH(2,2)" = m3, "GARCH(3,3)" = m4)
stargazer(summary_standard, title = "Model summaries", dep.var.labels = "Log returns (Training Set)", single.row = T,
          nobs = F, column.sep.width = "2.5pt", model.numbers = F, column.labels = c("GARCH(1,1)","GARCH(1,2)","GARCH(2,2)", "GARCH(3,3)"))


#INSIGHT: As we saw from our analysis before, the tails are very fat and we also have an asymmetry as expected for real log returns.
#We now develop a few special GARCH models to try to account for this.

m5 <- garchFit(~aparch(1,1), data=train_returns, delta = 1, include.delta = F) #this is a T GARCH (Treshold GARCH, estimated using std. dev)
fit5 <- fitted(m5)
summary(m5)

pdf("Figures/Analysis_Residuals_5.pdf") 
layout(matrix(c(1,1,2,3,4,5), nrow=3, ncol=2, byrow = TRUE))
plot(train$Date,train$log_returns, type="l", col="blue", main="Treshold GARCH(1,1) | Fitted Values", xlab="Time", ylab="Return")
lines(train$Date, fit5,  col=alpha("red", 0.5))
plot(train$Date, residuals(m5), type="l", main="Residuals", xlab="Time", ylab="Residual")
hist(residuals(m5), main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
qqnorm(residuals(m5), main="QQ-Plot | Residuals")
acf(residuals(m5)[-1][-1][-1]^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-1,1))
dev.off()

m6 <- garchFit(~aparch(1,1), data=train_returns, delta = 2, include.delta = F) #this is a GJR GARCH
fit6 <- fitted(m6)
summary(m6)

pdf("Figures/Analysis_Residuals_6.pdf") 
layout(matrix(c(1,1,2,3,4,5), nrow=3, ncol=2, byrow = TRUE))
plot(train$Date,train$log_returns, type="l", col="blue", main="GJR GARCH(1,1) | Fitted Values", xlab="Time", ylab="Return")
lines(train$Date, fit6,  col=alpha("red", 0.5))
plot(train$Date, residuals(m6), type="l", main="Residuals", xlab="Time", ylab="Residual")
hist(residuals(m6), main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
qqnorm(residuals(m6), main="QQ-Plot | Residuals")
acf(residuals(m6)[-1][-1][-1]^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-1,1))
dev.off()

#summary special models

summary_special <- list("GARCH(1,1)" = m1, "GARCH(3,3)" = m4, "T-GARCH(1,1)" = m5, "GJR-GARCH(1,1)" = m6)
stargazer(summary_special, title = "Model summaries", dep.var.labels = "Log returns (Training Set)", single.row = T,
          nobs = F, column.sep.width = "2.5pt", model.numbers = F, column.labels = c("GARCH(1,1)","GARCH(3,3)","T-GARCH(1,1)", "GJR-GARCH(1,1)"))

#INSIGHT: The special GARCHs that account for the asymmetry unfortunately do not really produce any better results. --> what now?

# 4 Model Forecast of GARCH ------------------------------------------------------

# o Split of training (fit) and testing (forecast) samples à Use/Plot testing sample
# o Define forecast horizon/lags n
# o Calculate volatility forecasts of General GARCH and best Special GARCH

# 5 Model Forecast Comparison with Realized Volatility (RV) ---------------

# - Model Forecast Comparison with Realized Volatility (RV)





# o Choose forecast horizon/lags n as before
# o Calculate volatility forecasts by Realized Volatility (RV)
#create function to compute RV
real_vola <- sd(xtrack_returns) * sqrt(21)

# o Compare RV forecast results with General GARCH and Special GARCH forecasts

# 6 Appendix --------------------------------------------------------------

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
