#Group member:
#Jing Fei XU
#Michael WARNER
#Suet Wah CHU

# -----------------------------------------------------
##### Traffic at CDG Airport Time Series Analysis #####
# -----------------------------------------------------
# Project Description :
# --------------------
# 
# This file loads traffic at Paris-Charles De Gaulle Airport per month from 2000 to 2019
# We will use data to forecast the next three years of traffic per month
# ----------------------------------------------

# load the package needed
library('TSA')
library('forecast')
library('stats')
library("zoo")
library("tseries")

# Import data
CdgTraf <- read.csv("CDG.csv")

#-------------------------------
#  Explanatory Analysis
#-------------------------------

# Convert the data set from a data frame to a time series
traf <- ts(CdgTraf[,3], start=c(2000,1), frequency=12)
traf
str(traf)

# Check the start date and end date of the time series
start(traf)
end(traf)

# Show a lag/frequency of 12 months
frequency(traf)

# Show the portion of each month to a year
time(traf)

# Check if there are no missing values
sum(is.na(traf))

# Check the summary of basic statistical parameters
summary(traf)

# Create a box plot to show the descriptive figures
boxplot(CdgTraf$Total, yaxt = 'n', ylab= "Traffic", main = "Traffic of CDG airport from 2000 to 2019")
axis(2, at= pretty(CdgTraf$Total), labels = format(pretty(CdgTraf$Total), scientific = FALSE, las = 1))


#---------------------------
# step 0: stationary issue
#---------------------------

# Plot of the series:there should be no increasing variance, no trend
# acf and pacf: they should decay fast to 0
plot.ts(traf)

# Add a linear regression line to see if there is variance and 
# to check if the data is stationary
abline(lm(traf~time(traf)))

# Decompose the series into trend seasonality and random component
x <- decompose(traf)
plot(x)
x
# Data shows an increasing trend and a repeating seasonality

# increasing variance: log transformation
# trend: first order difference
# seasonality: difference with lag s=12
par(mfrow=c(2,1))
acf(ts(traf,frequency = 1))
pacf(ts(traf,frequency = 1))
# ACF does not go to zero
# PACF does not quickly go to zero.

# log transformation: to reduce the increasing variance effect
ltraf <- log(traf)
par(mfrow=c(1,1))
plot.ts(ltraf)
par(mfrow=c(2,1))
acf(ts(ltraf,frequency = 1))
pacf(ts(ltraf,frequency = 1))
# ACF of ltraf still does not go to zero
# PACF of ltraf still does not quickly go to zero

# comparison: variance
plot.ts(cbind(traf,ltraf))

# data has a strong trend
# 1st order transformation: to remove the non constant trend
dltraf <- diff(ltraf,lag=1)
plot(dltraf)  #plot of differenced data
acf(ts(dltraf,frequency = 1),main="")  # not decaying fast to 0
pacf(ts(dltraf,frequency = 1),main="")

# series appears trend-stationary, use to investigate seasonality
ggseasonplot(dltraf)

# difference of lag 12: to remove seasonality
dltraf_12 <- diff(dltraf, lag=12)

plot(dltraf_12)
# we observe some extreme values on the graph

par(mfrow=c(2,1))
acf(ts(dltraf_12,frequency = 1),main="")
pacf(ts(dltraf_12,frequency = 1),main="")
par(mfrow=c(1,1))
# ACF, PACF: mostly non significant coefficient
# except at lag1, lag11, lag 12 and lag13




#------------------------------------------------
# step 1: identification of orders p,P, q, Q, d, D
#------------------------------------------------

# fit a multiplicative ARIMA on ltraf
# d=D=1: transformation performed on ltraf

# q, Q: MA part (ACF)

# p,P: AR part (PACF)


#-------------------------------
# step 2: estimation of the model
#-------------------------------

# create some models
mod1 <- arima(ltraf,c(1,1,1),seasonal = list(order=c(1,1,1), period=12), method='ML')
mod1
# P=Q=p=q=D=d=1
# AIC= -850.82

mod2 <- arima(ltraf, c(1,1,1), seasonal = list(order=c(0,1,1),period=12), method='ML')
mod2
# P=0(to remove the sar1), Q=p=q=D=d=1
# AIC= -852.78, the second model is a little bit better

mod3 <- arima(ltraf, c(0,1,1), seasonal = list(order=c(0,1,1),period=12), method='ML')
# P=p=0(to remove the sar1 and ar1), Q=q=D=d=1
# AIC= -851.99, the third model is a little bit worse
mod3

# use the auto.arima function to find the model with the lowest AIC
auto.arima(ltraf,d=1,D=1,stepwise = FALSE, approximation = FALSE, trace = FALSE)
# p=P=0, q=3,d=D=Q=1

mod4 <- stats::arima(ltraf, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), method='ML')
mod4
# p=P=0, q=3,d=D=Q=1
# AIC= -860.1,this one is the best


# plot of the fitted value
fit1 <- fitted(mod4)
plot.ts(cbind(ltraf,fit1), plot.type = 'single',col=c('black','red'))
plot.ts(cbind(traf,2.718^fit1), plot.type = 'single',col=c('black','red'))

#-------------------------------
#step 3: validation of the model
#-------------------------------

## significance of the coefficients ##
mod4
mod4$coef
mod4$var.coef

# student test statistic
tstat <- mod4$coef/sqrt(diag(mod4$var.coef))

# pvalue of the student test
pvalue <- 2*(1-pnorm(abs(tstat)))
tstat
pvalue
#ma1, ma3 and sma1 are significant (pvalue < 5%)

## residuals analysis ##
checkresiduals(mod4)
res1 <- mod4$residuals
plot(res1)

# ACF of the residuals (white noise assumption)
acf(ts(res1, frequency = 1))
# significant coefficients at lags 11

# PACF of the residuals
pacf(ts(res1, frequency = 1))
# significant coefficients at lags 11

# Ljung-Box test to check for the significance of ACF
Box.test(res1, lag = 20, type = "Ljung-Box")
# pvalue = 63% >5% so we accept HO: white noise assumption

# Heteroscedasticity
sq.res <- (res1)^2
acf(ts(sq.res,frequency=1))

# Normality assumption
# standardized residuals (residuals divided by the standard deviation)
res_stand <- res1/sqrt(mod4$sigma2)
summary(res_stand)

# If the normality assumption is satisfied, the standardized results should lie in between -2 and 2 (with 95% of chance)
plot(res_stand)
abline(a=2, b=0, col="red")
abline(a=-2, b=0, col="red")


# we can identify some outliers in the dataset
min(res_stand)
out1 <- which(res_stand < -5.5)
out1

index(res_stand)[out1] # date of the outlier
# the date of the outlier is April 2010

res_stand[out1]# value of the outlier
traf[out1]
traf[(out1-4):(out1+4)]
traf[(out1-16):(out1-8)]
# In April 2011, the value of traffic was similar to May 2011, it is not at all the case in 2010

# 2 options: change the value to a more plausible one or add a dummy in the model


# QQ-plot
qqnorm(res1)
qqline(res1)
# the points fall along a line in the middle of the graph,
# but curve off in the extremities.
# the residuals have more extreme values than would be expected for a Normal distribution
# so the residuals are not normally distributed


# Shapiro test
shapiro.test(res1)
#pvalue <<<< 5%, we reject normality (HO: normal dist)

# Jarque.bera test 
jarque.bera.test(res1)
#pvalue <<<< 5%, we reject normality (HO: normal dist)

# constant variance
sq.res <- (res1)^2 
acf(ts(sq.res, frequency = 1))

# The first correlation is significant so there seems to be a non constant variance.
# It may be due to the outliers we have identified on the graph of the residuals
# Option 1 : remove the outliers
# Option 2 : fit a GRACH model on the variance
# Here we choose option 1


#------------------------------------------------------------
# Dummy variables - Rejected as it does not improve the model
#------------------------------------------------------------


# Dummy variables were tried to improve the model as there are outliers in our residuals  

## residuals analysis ##
checkresiduals(mod4)
plot(res1)

# If the normality assumption is satisfied, the standardized results should lie in between -2 and 2 (with 95% of chance)
plot(res_stand)
abline(a=2, b=0, col="red")
abline(a=-2, b=0, col="red")

#------------------------------
# model with one dummy variable
#------------------------------
# we can identify one outlier in the dataset
# it corresponds to the lowest value of the residuals
min(res_stand)
out1 <- which(res_stand < -5.5)
out1
index(res_stand)[out1] # date of the outlier
# the date of the outlier is April 2010

res_stand[out1] # value of the outlier
traf[out1]
traf[(out1-4):(out1+4)]
traf[(out1-16):(out1-8)]
# In April 2011, the value of traffic was similar to May 2011, it is not at all the case in 2010


# Create a dummy variable
CdgTraf$dum_1 <- 0
CdgTraf$dum_1[out1] <- 1

mod5 <- stats::arima(ltraf, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), method='ML', xreg = CdgTraf$dum_1)
mod5
# mod5 AIC = -915.89

## Prediction model is tested against the base model that has no dummy variables
# Update fit model with improved model
fit2 <- fitted(mod5)
plot.ts(cbind(ltraf,fit2), plot.type = 'single',col=c('black','red'))
plot.ts(cbind(traf,2.718^fit2), plot.type = 'single',col=c('black','red'))


# Quality of the fit
cb80_2 <- mod5$sigma2^.5*qnorm(0.9)
plot(cbind(ltraf,fit2-cb80_2,fit2+cb80_2),plot.type='single',lty=c(1,2,2))

# Proportion of points in the confidence bound
indi_2 <- (ltraf-(fit2-cb80_2))>0&(fit2+cb80_2-ltraf)>0
prop_2 <- 100*sum(indi_2)/length(indi_2)
prop_2
# Model with one dummy variable is proportion with quality of fit of 87.5%, which is 
# worse than the model of no dummy variables.


#----------------------------------------
# model with two dummy variables
#----------------------------------------

# A second model was created to test if two dummy variables in place of the 
# two largest outliers would improve the accuracy of the model.

# Check residuals with largest outlier taken out
checkresiduals(mod5)
res2 <- mod5$residuals
plot(res2)

# standardized residuals (residuals divided by the standard deviation)
res_stand_2 <- res2/sqrt(mod5$sigma2)
summary(res_stand_2)
plot(res_stand_2) 
abline(a=2, b=0, col="red")
abline(a=-2, b=0, col="red")
# Around 2014 to 2015 there is another outlier that corresponds to the new 
# lowest value of the residuals

# Out2 given the value of the current outlier
min(res_stand_2)
out2 <- which(res_stand_2 < -5.5)
out2
index(res_stand)[out2] 
# the date of the outlier is Septembre 2014

# Dummy variable created for second outlier
CdgTraf$dum_2 <- 0
CdgTraf$dum_2[out2] <- 1

# mod7 given the arima function with the two dummy variables included
mod7 <- stats::arima(ltraf, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), 
                     method='ML', xreg = cbind(CdgTraf$dum_1, CdgTraf$dum_2))
mod7
#AIC=-965.21

# Residuals for new model checked
checkresiduals(mod7)
res4 <- mod7$residuals
plot(res4)

# standardized residuals (residuals divided by the standard deviation)
res_stand_3 <- res4/sqrt(mod7$sigma2)
plot(res_stand_3)
summary(res_stand_3)

# Plot of the new fitted value with second and first largest outliers removed
fit3 <- fitted(mod7)
plot.ts(cbind(ltraf,fit3), plot.type = 'single', col=c('black','red'))

# Plot of the two models overlapped, red line is the largest outlier removed, green is the 
# second largest and first largest removed. Improvements to the model are almost trivial,
# therefore no more outliers are removed to avoid overfitting.

plot.ts(cbind(ltraf, fit2, fit3), plot.type = 'single', col=c('black','red', 'green'))


# Quality of the fit
cb80_3 <- mod7$sigma2^.5*qnorm(0.9)
plot(cbind(ltraf,fit3-cb80_3,fit3+cb80_3),plot.type='single',lty=c(1,2,2))

# Proportion of points in the confidence bound
indi_3 <- (ltraf-(fit3-cb80_3))>0&(fit3+cb80_3-ltraf)>0
prop_3 <- 100*sum(indi_3)/length(indi_3)
prop_3

# Model with two dummy variables proportion of points is worse, at 82.92%
# This shows that the model with no dummy variables should be used.



#-------------------------------------------------------------
# Step4: Prediction
#-------------------------------------------------------------

# Prediction with ARIMA (0,1,3)(0,1,1)[12] model with no dummy variables


# Quality of the fit
cb80 <- mod4$sigma2^.5*qnorm(0.9)
plot(cbind(ltraf,fit1-cb80,fit1+cb80),plot.type='single',lty=c(1,2,2))

# Proportion of points in the confidence bound
indi <- (ltraf-(fit1-cb80))>0&(fit1+cb80-ltraf)>0
prop <- 100*sum(indi)/length(indi)
prop
# If prop >80%, then the fit is considered good.(88.75%)


# Prediction 
pred <- predict(mod4, n.ahead = 36)
ts.plot(traf,2.718^pred$pred, log = "y", lty = c(1,3))
time(traf)
ts.plot(traf,xlim=c(2010,2022.500))
lines(2.718^(pred$pred), col="red")
lines(2.718^(pred$pred-1.96*pred$se),col=4,lty=2)
lines(2.718^(pred$pred+1.96*pred$se),col=4,lty=2)

# -------------------
# 4.1 In-Sample test
#--------------------


## We will use the data from 2000 to 2016 to predict the traffic from 2017 to 2019
# and then to compare our results with the real traffic data for these three years 
# in order to test the accuracy of our model
# so we just follow the same steps as before

#Import data
CdgTraf_train <- read.csv("CDG.csv")
CdgTraf_train <- CdgTraf_train[1:204,]
traf_train <- ts(CdgTraf_train[,3], start=c(2000,1), frequency=12)




#-------------------------------------------
#  4.1.1 stationary issue for the train set
#-------------------------------------------

plot.ts(traf_train)
x_train <- decompose(traf_train)
plot(x_train)
par(mfrow=c(2,1))
acf(ts(traf_train,frequency = 1))
pacf(ts(traf_train,frequency = 1)) 

# log transformation: to reduce the increasing variance effect
ltraf_train<- log(traf_train)
plot.ts(ltraf_train)

# 1st order transformation: to remove the non constant trend
dltraf_train <- diff(ltraf_train,lag=1)
plot(dltraf_train)

# difference of lag 12: to remove seasonality
dltraf_12_train <- diff(dltraf_train, lag=12)
plot(dltraf_12_train)
acf(ts(dltraf_12_train,frequency = 1),main="")
pacf(ts(dltraf_12_train,frequency = 1),main="")


#-------------------------------
# 4.1.2 estimation of the model
#-------------------------------

# use the auto.arima function to find the model with the lowest AIC
auto.arima(ltraf_train,d=1,D=1,stepwise = FALSE, approximation = FALSE, trace = FALSE)
# get the same result as before

mod6 <- stats::arima(ltraf_train, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), method='ML')
mod6 
# AIC=-696.6


# plot of the fitted value
fit1_test <- fitted(mod6)
plot.ts(cbind(ltraf_train,fit1_test), plot.type = 'single',col=c('black','red'))
plot.ts(cbind(traf_train,2.718^fit1_test), plot.type = 'single',col=c('black','red'))




#-------------------------------
# 4.1.3 validation of the model
#-------------------------------
## significance of the coefficients ##
mod6
mod6$coef
mod6$var.coef

# student test statistic
tstat_train <- mod6$coef/sqrt(diag(mod6$var.coef))

# pvalue of the student test
pvalue_train <- 2*(1-pnorm(abs(tstat_train)))
tstat_train
pvalue_train
# ma1, ma3 and sma1 are significant (pvalue < 5%)

## residuals analysis ##
checkresiduals(mod6)
res3 <- mod6$residuals
plot(res3)
acf(ts(res3, frequency = 1))
pacf(ts(res3, frequency = 1))

# Ljung-Box test 
Box.test(res3, lag = 20, type = "Ljung-Box")
# pvalue = 68.61% >5% so we accept HO: white noise assumption

# Normality assumption
# standardized residuals 
res_stand_test <- res3/sqrt(mod6$sigma2)
summary(res_stand_test)

# If the normality assumption is satisfied, the standardized results should lie in between -2 and 2 (with 95% of chance)
plot(res_stand_test)
abline(a=2, b=0, col="red")
abline(a=-2, b=0, col="red")

# QQ-plot
qqnorm(res3)
qqline(res3)

# Shapiro test
shapiro.test(res3)
# pvalue <<<< 5%, we reject normality (HO: normal dist)

# Jarque.bera test 
jarque.bera.test(res3)
# pvalue <<<< 5%, we reject normality (HO: normal dist)




#-----------------------------------------------------
# 4.1.4 Prediction of In-sample test with train set
#-----------------------------------------------------

# Quality of the fit
cb80_train <- mod6$sigma2^.5*qnorm(0.9)
plot(cbind(ltraf_train,fit1_test-cb80_train,fit1_test+cb80_train),plot.type='single',lty=c(1,2,2))

# Proportion of points in the confidence bound
indi_test <- (ltraf_train-(fit1_test-cb80_train))>0&(fit1_test+cb80_train-ltraf_train)>0
prop.test<- 100*sum(indi_test)/length(indi_test)
prop.test
# If prop >80%, then the fit is considered good.(88.24%)

pred_test <- predict(mod6, n.ahead = 36)
pred_test

# Prediction 
ts.plot(traf_train,xlim=c(2010,2019.500))
lines(2.718^(pred_test$pred), col="red")
lines(2.718^(pred_test$pred-1.96*pred_test$se),col=4,lty=2)
lines(2.718^(pred_test$pred+1.96*pred_test$se),col=4,lty=2)

ts.plot(traf,xlim=c(2010,2019.500))
lines(2.718^(pred_test$pred), col="red")
lines(2.718^(pred_test$pred-1.96*pred_test$se),col=4,lty=2)
lines(2.718^(pred_test$pred+1.96*pred_test$se),col=4,lty=2)

# Comparison between the prediction and the real data (2016~2019)
ts.plot(traf,2.718^pred_test$pred, log = "y", lty = c(1,3))
ts.plot(traf,2.718^pred_test$pred, log = "y", lty = c(1,3),xlim=c(2016,2019.500))


# -------------------------------
# 4.2 Out-of-Sample Test
# -------------------------------

# Prediction for 2020-2022
pred <- predict(mod4, n.ahead = 36)
ts.plot(traf,2.718^pred$pred, log = "y", lty = c(1,3))
time(traf)
ts.plot(traf,xlim=c(2010,2022.500))
lines(2.718^(pred$pred), col="red")
lines(2.718^(pred$pred-1.96*pred$se),col=4,lty=2)
lines(2.718^(pred$pred+1.96*pred$se),col=4,lty=2)



