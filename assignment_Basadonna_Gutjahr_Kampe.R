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
library(texreg)
library(rugarch)
library(gridExtra)
library(grid)
library(lattice)
library(Metrics)
library(aod)
library(clubSandwich)
library(tstools)
options("modelsummary_format_numeric_latex" = "plain")

# Working directory
getwd()
setwd("~/Documents/GitHub/Financial_Volatility") #NIK
#setwd("C:/Users/jan_g/OneDrive/HSG/2nd semester/Financial volatility/Assignment") #JAN

# 1 Data Import & Cleaning --------------------------------------------------------------

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
sink(file = "Latex/summary_stats.txt")
xtable(summ_stats)
sink(file = NULL)

# 3 Model Fit and Forecast Comparisons -------------------------------------------

# 3.1 Check for stylized facts of the time series --------------------------------

# 3.1.1 Test of Stylized Facts: Distribution, Moments, T-Test, JB Test -----------
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

t.test(xtrack_returns) # --> INSIGHT: Zero mean with 95% confidence
normalTest(xtrack_returns, method = "jb") # --> INSIGHT: Not normaly distirbuted
skewness(xtrackers_msci$log_returns, na.rm=TRUE) # --> INSIGHT: Slightly negatively skewed
kurtosis(xtrackers_msci$log_returns, na.rm=TRUE) # --> INSIGHT: Fat/Heavy tails (interesting, note to self)

# 3.1.1 Check for autocorrelation for different lags: Ljung-Box Test -------------
lb_lags <- c(1,3,5,7,10,15,20)
lb_test_summary <- data.frame(row.names=c("Ljung-Box Statistic", "p-value"))
for (i in 1:length(lb_lags)) {
  lb_test <- Box.test(xtrack_returns, lag=lb_lags[i], type="Ljung") 
  lb_test_summary[1,i] <- round(lb_test$statistic, digits=6)
  lb_test_summary[2,i] <- round(lb_test$p.value, digits=6)
}
colnames(lb_test_summary) <- c("Q=1","Q=3","Q=5","Q=7","Q=10","Q=15","Q=20")
sink(file = "Latex/LB_Test_Summary.txt")
xtable(lb_test_summary)
sink(file = NULL) 

# 3.1.2 Check for (G)ARCH effects: LM (Lagrange Multiplier) Test -----------------
lm_test_summary <- data.frame(row.names=c("LM Statistic", "p-value", "Result"))
lm_test <- Lm.test(xtrack_returns, lag.max = 5, alpha = 0.05) 
lm_test_summary[1,1] <- round(lm_test$statistic, digits=6)
lm_test_summary[2,1] <- round(lm_test$p.value, digits=6)
lm_test_summary[3,1] <- lm_test$Conc
colnames(lm_test_summary) <- c("LM Test")
sink(file = "Latex/LM_Test_Summary.txt")
xtable(lm_test_summary)
sink(file = NULL) 

# 3.1.3 Check for asymmetries & leverage effect: ---------------------------------

# 3.1.3.1 Leverage Effect: Autocorrelations --------------------------------------
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

## 3.1.3.2 Asymmetry: Test Regression --------------------------------------------
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

# 3.3.2 Model fit & Determination of model parameters --------------------------

## Global Model Fit Variables
p_max <- 3
q_max <- 3

## Function for rugarch summary
summary.rugarch <- function(model, modelname, include.loglike = TRUE, include.aic = TRUE, include.bic = TRUE, include.hq = TRUE) {
  coefnames <- rownames(as.data.frame(model@fit$coef))
  coefs <- model@fit$coef
  se <- as.vector(model@fit$matcoef[, c(2)])
  pvalues <-  as.vector(model@fit$matcoef[, c(4)]) 
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglike == TRUE) {
    loglike <- model@fit$LLH
    gof <- c(gof, loglike)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.aic == TRUE) {
    aic <- infocriteria(model)[c(1)]
    gof <- c(gof, aic)
    gof.names <- c(gof.names, "AIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.bic == TRUE) {
    bic <- infocriteria(model)[c(2)]
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "BIC")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.hq == TRUE) {
    bic <- infocriteria(model)[c(4)]
    gof <- c(gof, bic)
    gof.names <- c(gof.names, "HQ")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  tr <- createTexreg(
    model.name = modelname,
    coef.names = coefnames, 
    coef = coefs,
    se = se,
    pvalues = pvalues, 
    gof.names = gof.names, 
    gof = gof, 
    gof.decimal = gof.decimal
  )
  return(tr)
}

# 3.3.2.1 Standard GARCH Models ------------------------------------------------
models_standard <- c()
ics_standard <- matrix(0,9,6)
colnames(ics_standard) <- c("p","q","LLH","AIC","BIC","HQ")
summaries_standard <- c()
forecasts_standard <- c()
k <- 1
for (p in 1:p_max) {
  for(q in 1:q_max) {
    # Model Fit
    spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder=c(p,q)), distribution.model="sstd")
    model <- ugarchfit(spec=spec, data=xtrack_returns, include.mean=FALSE, out.sample=floor(0.3*nrow(xtrackers_msci)))
    forecast <- ugarchforecast(model, data=NULL, out.sample=floor(0.3*nrow(xtrackers_msci)), n.ahead = 1, n.roll = 584)
    fit <- model@fit$sigma
    res <- model@fit$residuals
    models_standard <- c(models_standard, model)
    forecasts_standard <- c(forecasts_standard, forecast)
    # Plots
    indices <- c(1,2,3,8,9,10,11,12)
    pdf(paste("Figures/Summary_Plots_GARCH_",p,"_",q,".pdf", sep="")) 
    par(mfrow = c(4,2))
    for (i in indices) {
      plot(model, which=i)
    }
    dev.off()
    # Log-Likelihood & Information Criteria
    llh <- model@fit$LLH
    aic <- infocriteria(model)[1]
    bic <- infocriteria(model)[2]
    hq <- infocriteria(model)[4]
    ics_standard[k,] <- c(p,q,llh,aic,bic,hq)
    # Plot & Residual Analysis
    pdf(paste("Figures/Analysis_Residuals_GARCH_",p,"_",q,".pdf", sep="")) 
    layout(matrix(c(1,1,2,2,3,4,5,6), nrow=4, ncol=2, byrow = TRUE))
    plot(xtrackers_msci$Date[1:(floor(0.7*nrow(xtrackers_msci))+1)], abs(xtrackers_msci$log_returns[1:(floor(0.7*nrow(xtrackers_msci))+1)]), type="l", col="blue", main="Absolute Log Return Series", xlab="Time", ylab="Abs. Return")
    plot(xtrackers_msci$Date[1:(floor(0.7*nrow(xtrackers_msci))+1)], fit, type="l", col="red", main=paste("GARCH (",p,",",q,") | Estimation of Conditional Volatility", sep=""), xlab="Time", ylab="Cond. Volatility", ylim=c(0,4.5))
    plot(xtrackers_msci$Date[1:(floor(0.7*nrow(xtrackers_msci))+1)], res, type="l", main="Residuals", xlab="Time", ylab="Residual")
    hist(res, main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
    qqnorm(res, main="QQ-Plot | Residuals")
    acf(res^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-0.5,0.5))
    dev.off()
    # Model Summary
    summaries_standard <- c(summaries_standard, summary.rugarch(model=model, modelname=paste("GARCH(",p,",",q,")", sep="")))
    k <- k+1
  }
}
# Model Summary Aggregation
sink(file="Latex/Summary_Standard_GARCH.txt")
texreg(summaries_standard, digits=4)
sink(file = NULL)

## Standard GARCHs - Best Model
ics_standard <- as.data.frame(ics_standard)
best_orders <- ics_standard[which.min(ics_standard$AIC), c("p", "q")]
spec_best <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder=c(best_orders[[1]],best_orders[[2]])), distribution.model="sstd")
best_standard_garch <- ugarchfit(spec=spec_best, data=xtrack_returns, out.sample=floor(0.3*length(xtrack_returns)))

# 3.3.2.2 Special GARCH Models ------------------------------------------------
best_orders = c(1,2)
models_special <- c()
model_types <- c("TGARCH", "eGARCH", "GJRGARCH")
model_names_summary <- c("T-GARCH", "E-GARCH", "GJR-GARCH")
ics_special <- matrix(0,9,6)
colnames(ics_special) <- c("p","q","LLH","AIC","BIC","HQ")
summaries_special <- c()
forecasts_special <- c()
k <- 1
for (model_type in model_types) {
  # Model Fit
  if (model_type == "eGARCH") {
    spec <- ugarchspec(variance.model=list(model=model_type, garchOrder=c(best_orders[[1]],best_orders[[2]]), submodel=NULL), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="sstd")
  } else {
    spec <- ugarchspec(variance.model=list(model="fGARCH", garchOrder=c(best_orders[[1]],best_orders[[2]]), submodel=model_type), mean.model=list(armaOrder=c(0,0), include.mean=FALSE), distribution.model="sstd")
  }
  model <- ugarchfit(spec=spec, data=xtrack_returns, out.sample=floor(0.3*nrow(xtrackers_msci)))
  forecast <- ugarchforecast(model, data=NULL, out.sample=floor(0.3*nrow(xtrackers_msci)), n.ahead = 1, n.roll = 584)
  fit <- model@fit$sigma
  res <- model@fit$residuals
  models_special <- c(models_special, model)
  forecasts_special <- c(forecasts_special,forecast)
  # Plots
  indices <- c(1,2,3,8,9,10,11,12)
  pdf(paste("Figures/Summary_Plots_",model_type,"_",best_orders[[1]],"_",best_orders[[2]],".pdf", sep="")) 
  par(mfrow = c(4,2))
  for (i in indices) {
    plot(model, which = i)
  }
  dev.off()
  
  # Information Criteria
  llh <- model@fit$LLH
  aic <- infocriteria(model)[1]
  bic <- infocriteria(model)[2]
  hq <- infocriteria(model)[4]
  ics_special[k,] <- c(best_orders[[1]],best_orders[[2]],llh,aic,bic,hq)
  # Plot & Residual Analysis
  
  pdf(paste("Figures/Analysis_Residuals_",model_type,"_",best_orders[[1]],"_",best_orders[[2]],".pdf", sep="")) 
  layout(matrix(c(1,1,2,2,3,4,5,6), nrow=4, ncol=2, byrow = TRUE))
  plot(xtrackers_msci$Date[1:(floor(0.7*nrow(xtrackers_msci))+1)], abs(xtrackers_msci$log_returns[1:(floor(0.7*nrow(xtrackers_msci))+1)]), type="l", col="blue", main="Absolute Log Return Series", xlab="Time", ylab="Abs. Return")
  plot(xtrackers_msci$Date[1:(floor(0.7*nrow(xtrackers_msci))+1)], fit, type="l", col="red", main=paste(model_names_summary[k]," (",best_orders[[1]],",",best_orders[[2]],") | Estimation of Conditional Volatility", sep=""), xlab="Time", ylab="Cond. Volatility", ylim=c(0,4.5))
  plot(xtrackers_msci$Date[1:(floor(0.7*nrow(xtrackers_msci))+1)], res, type="l", main="Residuals", xlab="Time", ylab="Residual")
  hist(res, main="Histogram | Residuals", breaks = 20, xlab="Residual", ylab="Count")
  qqnorm(res, main="QQ-Plot | Residuals")
  acf(res^2, lag.max=max_lags, main="ACF | Squared Residuals", ylim=range(-1,1))
  dev.off()
  # Model Summary
  summaries_special <- c(summaries_special, summary.rugarch(model=model, modelname=paste(model_names_summary[k],"(",best_orders[[1]],",",best_orders[[2]],")", sep="")))
  k <- k+1
}
# Model Summary Aggregation
sink(file="Latex/Summary_Special_GARCH.txt")
texreg(summaries_special, digits=4)
sink(file = NULL)

## Special GARCHs - Best Model
ics_special <- as.data.frame(ics_special)

# 3.3.2.3 Realized Volatility: HAR Model ---------------------------------------
## Dependent Variable
rv_d_h <- c()
for (i in 23:length(xtrackers_msci$log_returns[1:(floor(0.7*nrow(xtrackers_msci))+1)])) {
  rv_d_h <- c(rv_d_h, xtrackers_msci[i,"High"] / xtrackers_msci[i,"Low"])
}
## Independent Variables
rv_d <- rv_w <- rv_m <- c()
for (i in 22:length(xtrackers_msci$log_returns[1:(floor(0.7*nrow(xtrackers_msci))+1)])) {
  # Daily Realized Variance
  rv_d <- c(rv_d, xtrackers_msci[i,"High"] / xtrackers_msci[i,"Low"]) 
  # Weekly Average of Daily Realized Variance
  rv_w_temp <- c()
  for (j in 1:5) {
    rv_w_temp <- c(rv_w_temp, xtrackers_msci[i-j-1,"High"] / xtrackers_msci[i-j-1,"Low"])
  }
  rv_w <- c(rv_w, sum(rv_w_temp)/length(rv_w_temp))
  # Weekly Average of Daily Realized Variance
  rv_m_temp <- c()
  for (k in 1:21) {
    rv_m_temp <- c(rv_m_temp, xtrackers_msci[i-k-1,"High"] / xtrackers_msci[i-k-1,"Low"])
  }
  rv_m <- c(rv_m, sum(rv_m_temp)/length(rv_m_temp))
}

rv_d <- rv_d[-length(rv_d)]
rv_w <- rv_w[-length(rv_w)]
rv_m <- rv_m[-length(rv_m)]

## Regression with 1-lag 
har <- lm(rv_d_h ~ rv_d + rv_w + rv_m)
har_coef <- har$coefficients
# Model Summary
sink(file="Latex/Summary_HAR.txt")
texreg(summary(har, digits=4), digits=4)
sink(file = NULL)





# 4. Forecasting ----------------------------------------------------------

#showing the training and test split in action

train <- data.frame(x= xtrackers_msci$log_returns[1:(floor(0.7*nrow(xtrackers_msci))+1)])
train <- data.frame(train, D = xtrackers_msci$Date[1:1368])
test <- data.frame(y = xtrackers_msci$log_returns[(floor(0.7*nrow(xtrackers_msci))+2):1953])
test <- data.frame(test, D = xtrackers_msci$Date[1369:1953])

pdf("Figures/train_test_graph.pdf", width = 9, height = 4)
ggplot() + 
  geom_line(data = train, aes(x = D, y = x), color = "lightblue") +
  geom_line(data = test, aes(x = D, y = y), color = "red") +
  xlab('time') +
  ylab('log returns') +
  ggtitle("Training/Test Split across Time Series",
          subtitle = "Start of test sample: January, 16th, 2020")
dev.off()

# 4.1 Actual realized volatilities for the test set -----------------------


#index f on the dataframes indicates "full"
rv_d_f <- rv_w_f <- rv_m_f <- c()
for (i in 22:length(xtrackers_msci$log_returns)) {
  # Daily Realized Variance
  rv_d_f <- c(rv_d_f, xtrackers_msci[i,"High"] / xtrackers_msci[i,"Low"]) 
  # Weekly Average of Daily Realized Variance
  rv_w_temp <- c()
  for (j in 1:5) {
    rv_w_temp <- c(rv_w_temp, xtrackers_msci[i-j-1,"High"] / xtrackers_msci[i-j-1,"Low"])
  }
  rv_w_f <- c(rv_w_f, sum(rv_w_temp)/length(rv_w_temp))
  # Weekly Average of Daily Realized Variance
  rv_m_temp <- c()
  for (k in 1:21) {
    rv_m_temp <- c(rv_m_temp, xtrackers_msci[i-k-1,"High"] / xtrackers_msci[i-k-1,"Low"])
  }
  rv_m_f <- c(rv_m_f, sum(rv_m_temp)/length(rv_m_temp))
}

rv_d_f <- rv_d_f[-length(rv_d)]
rv_w_f <- rv_w_f[-length(rv_w)]
rv_m_f <- rv_m_f[-length(rv_m)]


# 4.2 Forecasts HAR -------------------------------------------------------
forecasts_har <- c()
for (x in 0:584) {
  if (x > 0){
    print(paste("HAR Model refit no: ",x,"/584",sep = ""))
    
  } else {}
  #prediction for 69 w/ data from 68
  pred_new <- predict(har, data.frame(rv_d = rv_d_f[floor(0.7*nrow(xtrackers_msci))+1+x-22],rv_w = rv_w_f[floor(0.7*nrow(xtrackers_msci))+1+x-22], rv_m = rv_m_f[floor(0.7*nrow(xtrackers_msci))+1+x-22]))
  forecasts_har <- c(forecasts_har,pred_new)
  # Refitting of the model with new true value
  ## Dependent Variable
  rv_d_h <- c()
  for (i in 23:length(xtrackers_msci$log_returns[1:(floor(0.7*nrow(xtrackers_msci))+1+x)])) {
    rv_d_h <- c(rv_d_h, xtrackers_msci[i,"High"] / xtrackers_msci[i,"Low"])
  }
  # Independent Variables
  rv_d <- rv_w <- rv_m <- c()
  for (i in 22:length(xtrackers_msci$log_returns[1:(floor(0.7*nrow(xtrackers_msci))+1+x)])) {
    # Daily Realized Variance
    rv_d <- c(rv_d, xtrackers_msci[i,"High"] / xtrackers_msci[i,"Low"])
    # Weekly Average of Daily Realized Variance
    rv_w_temp <- c()
    for (j in 1:5) {
      rv_w_temp <- c(rv_w_temp, xtrackers_msci[i-j,"High"] / xtrackers_msci[i-j,"Low"])
    }
    rv_w <- c(rv_w, sum(rv_w_temp)/length(rv_w_temp))
    # Weekly Average of Daily Realized Variance
    rv_m_temp <- c()
    for (k in 1:21) {
      rv_m_temp <- c(rv_m_temp, xtrackers_msci[i-k,"High"] / xtrackers_msci[i-k,"Low"])
    }
    rv_m <- c(rv_m, sum(rv_m_temp)/length(rv_m_temp))
  }
  rv_d <- rv_d[-length(rv_d)]
  rv_w <- rv_w[-length(rv_w)]
  rv_m <- rv_m[-length(rv_m)]
  har <- lm(rv_d_h ~ rv_d + rv_w + rv_m)
}

forecasts_har_true_values_plot <- data.frame(x = xtrackers_msci$Date[1369:1953], y= (forecasts_har-1)*100, z = abs(xtrackers_msci$log_returns[1369:1953]))

#produce one plot of the forecasted RV against the true values
pdf("Figures/har_forecast.pdf")
ggplot(forecasts_har_true_values_plot, aes(x,y,z, color = variable)) +
  geom_line(aes(y = z,  color =c("Absolute Log Returns"))) +
  geom_line(aes(y = y, color ="Forecast")) +
  ggtitle("HAR forecast vs. absolute log returns") +
  ylab("Sigma") +
  xlab("Time") +
  labs(color="Legend") +
  scale_color_manual(values = c("grey", "blue"))
dev.off()


# 4.3 forecasts GARCH -----------------------------------------------------
## GARCH forecasts we already created the predictions from our GARCH models when fitting the model (see above))
forecasts_garch <- c(forecasts_standard, forecasts_special)

#line graphs for forecast accuary for every GARCH
for (i in 1:length(forecasts_garch)) {
  
  forecasts_garch_true_values_plot <- data.frame(x = xtrackers_msci$Date[1369:1953], y= forecasts_garch[[i]]@forecast$sigmaFor[1:585], z = abs(xtrackers_msci$log_returns[1369:1953]))
  
  pdf(paste("Figures/garch_forecast_",i,".pdf", sep=""))
  print(ggplot(forecasts_garch_true_values_plot, aes(x,y,z, color = variable)) +
          geom_line(aes(y = z,  color =c("Absolute log returns"))) +
          geom_line(aes(y = y, color ="Forecast")) +
          ggtitle(paste("GARCH Model", i, "forecast vs. absolute log returns")) +
          ylab("Sigma") +
          xlab("Time") +
          labs(color="Legend") +
          scale_color_manual(values = c("grey", "red")))
  dev.off()
}






# 5. Comparison of forecasts ----------------------------------------------

# 5.1 compute performance summary and plots of GARCH models ---------------
#performance summary
performance_summary_forecasts <- data.frame(matrix(ncol = 6, nrow = 2)) #table of summary statistics

model_selection <- c(forecasts_garch[[1]], forecasts_garch[[2]], forecasts_garch[[10]],forecasts_garch[[11]],forecasts_garch[[12]]) #number of models

length(as.vector(forecasts_garch[[1]]@forecast$sigmaFor)) #number of test observations from garch models
length(forecasts_har) #number of test observations from HAR model
length(abs(xtrackers_msci$log_returns[1369:1953]))

length(as.vector(forecasts_garch[[11]]@forecast$sigmaFor))

#perf summary and model plots for all GARCH models
for (i in 1:length(model_selection)) {
  
  performance_summary_forecasts[1,i] <- mse(as.vector(model_selection[[i]]@forecast$sigmaFor)^2,xtrackers_msci$log_returns[1369:1953]^2)
  performance_summary_forecasts[2,i] <- rmse(as.vector(model_selection[[i]]@forecast$sigmaFor)^2,xtrackers_msci$log_returns[1369:1953]^2)
  
}

performance_summary_forecasts[1,6] <- mse(((forecasts_har-1)*100)^2, xtrackers_msci$log_returns[1369:1953]^2)
performance_summary_forecasts[2,6] <- rmse(((forecasts_har-1)*100)^2, xtrackers_msci$log_returns[1369:1953]^2)


#tidying the data frame
rownames(performance_summary_forecasts) <- c("MSE", "RMSE")
colnames(performance_summary_forecasts) <- c("G(1,1) (benchmark)", "G(1,2)", "T-G(1,2)", "E-G(1,2)", "GJR-G(1,2)", "HAR")

#report the performance summary in stargazer
sink(file = "Latex/performance_summary_GARCH.txt")
stargazer(round(performance_summary_forecasts,4), summary = F, column.sep.width	= "1pt",
          type = "latex", title = "Performance Statistics of all models versus the HAR estimate",
          font.size = "normalsize")
sink(file= NULL)

# 5.2 select best models and plot their forecasts against the one from RV (incl. real values) ------

forecasts_garch_har_true_values_plot <- data.frame(x = xtrackers_msci$Date[1369:1953], z = forecasts_garch[[2]]@forecast$sigmaFor[1:585], b = forecasts_garch[[10]]@forecast$sigmaFor[1:585], y = (forecasts_har-1)*100, t = abs(xtrackers_msci$log_returns[1369:1953]))

pdf("Figures/har_garch_comparison.pdf",width = 10, height = 5)
ggplot(forecasts_garch_har_true_values_plot, aes(color = variable)) +
  geom_line(aes(y = t, x = x, color =c("Absolute Log Returns"))) +
  geom_line(aes(y = z, x = x, color = c("GARCH(1,2) Forecast"))) +
  geom_line(aes(y = b, x = x, color = c("T-GARCH(1,2) Forecast"))) +
  geom_line(aes(y = y, x = x, color = c("HAR Forecast"))) +
  ggtitle("Standard and special GARCH model forecast vs HAR model forecast") +
  ylab("Sigma") +
  xlab("Time") +
  labs(color="Legend") +
  scale_color_manual(values = c("grey","green","red","black"))
dev.off()


# 5.2 HAR forecast plot ---------------------------------------------------

#create data fram wit hHAR forecasts and true values
forecasts_har_true_values_plot <- data.frame(y= (forecasts_har-1)*100, x = xtrackers_msci$Date[1369:1953], z = abs(xtrackers_msci$log_returns[1369:1953]))

### FIND THE CORRECT Series ###

#produce desired plots
ggplot(forecasts_har_true_values_plot, aes(x,y,z)) +
  geom_line(aes(y = z,  color ="True values")) +
  geom_line(aes(y = y, color ="Forecast")) +
  ggtitle("HAR forecast vs. true Sigma") +
  ylab("Sigma") +
  xlab("Time") +
  labs(color="Legend") +
  scale_color_manual(values = c("red", "grey"))


# 5.3 Performance comparison -----------------------------------------

#create the regression for multiple models

regression_data <- data.frame(t = c(xtrackers_msci$log_returns[1369:1953])^2, m1 = c(forecasts_garch[[1]]@forecast$sigmaFor)^2, m2 = c(forecasts_garch[[2]]@forecast$sigmaFor)^2, m3 = c(forecasts_garch[[11]]@forecast$sigmaFor)^2, h = ((forecasts_har-1)*100)^2)

MZ_models <- c()

for (i in 2:5) {
  
  MZ_models[[i-1]] <- lm(t ~ regression_data[,i], data = regression_data)
  
}

sink(file = "Latex/MZ_models.txt")
stargazer(MZ_models[[1]], MZ_models[[2]], MZ_models[[3]], MZ_models[[4]], covariate.labels = c("Beta","Alpha"), dep.var.labels	= c("Squared Log Returns"),
          column.labels = c("G(1,1)", "G(1,2)","T-G(1,2)","HAR"), object.names = F, omit.stat = c("ser","f","adj.rsq"))
sink(file = NULL)

#DMW Test
dmwtest_results <- data.frame(matrix(ncol=1, nrow=2))

#create the squared returns
proxy <- xtrackers_msci$log_returns[1369:1953]^2

resid_garch12 <- c(forecasts_garch[[2]]@forecast$sigmaFor)^2 - proxy
resid_tgarch12 <- c(forecasts_garch[[10]]@forecast$sigmaFor)^2 - proxy
resid_har <- c((forecasts_har-1)*100)^2 - proxy

d <- dm.test(resid_garch12,resid_tgarch12, alternative = c("two.sided"), h = 1, power = 2)

dmwtest_results[1,1] <- d$statistic
dmwtest_results[2,1] <- d$p.value
x
rownames(dmwtest_results) <- c("DMW statistic", "p value")
colnames(dmwtest_results) <- c("G(1,2) vs T-G(1,2)")

sink(file = "Latex/dmw_test.txt")
stargazer(round(dmwtest_results,4), summary = F, column.sep.width	= "1pt",
          type = "latex", title = "Performance Statistics of all models versus the HAR estimate",
          font.size = "normalsize")
sink(file = NULL)