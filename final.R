# Traffic at CDG Airport 
# This file loads traffic at Paris-Charles De Gaulle Airport per month from 2000 to 2019
# We will use data to forecast the next three years of traffic per month


#load the package needed
library('TSA')
library('forecast')
library(ggplot2)
library(scales)
#Import data
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


#------------------------
#step 0: stationary issue
#------------------------

# plot of the series:there should be no increasing variance, no trend
#acf and pacf: they should decay fast to 0
plot.ts(traf)

# Add a linear regression line to see if there is variance and to check if the data is stationary
abline(lm(traf~time(traf)))

# Boxplot by cycle
# Show the seasonality
# July and August have the higest no. during the year
# November, December, January and February are relatively low
boxplot(traf~cycle(traf, xlab = "Month", ylab="Traffic", main = "Traffic of CDG aiport from 2000 to 2019"))

# decompose the series into trend seasonality and random component
x <- decompose(traf, type = "multiplicative")
plot(x)
x
#Data shows an increasing trend
#data shows a repeating seasonality
#Data shows ??? in randomness


#increasing variance: log transformation
#trend: first order difference
#seasonality: difference with lag s=12
par(mfrow=c(2,1))
acf(ts(traf,frequency = 1))
pacf(ts(traf,frequency = 1))
#ACF does not go to zero
#PACF does not quickly go to zero.



#log transformation: to reduce the increasing variance effect
ltraf <- log(traf)
plot.ts(ltraf)
acf(ts(ltraf,frequency = 1))
pacf(ts(ltraf,frequency = 1))
#comparison: variance
plot.ts(cbind(traf,ltraf))

#ACF of ltraf still does not go to zero
#PACF of ltraf ??? still does not quickly go to zero


#data has a strong trend
#to remove the trend, we take the first difference
#1st order transformation: to remove the non constant trend
dltraf <- diff(ltraf,lag=1)
plot(ltraf)
plot(dltraf)  #plot of differenced data
acf(ts(dltraf,frequency = 1),main="")  # not decaying fast to 0
pacf(ts(dltraf,frequency = 1),main="")

#series appears trend-stationary, use to investigate seasonality
ggseasonplot(dltraf)

#difference of lag 12: to remove seasonality
dltraf_12 <- diff(dltraf, lag=12)

plot(dltraf_12)
#we observe some extreme values on the graph

par(mfrow=c(2,1))
acf(ts(dltraf_12,frequency = 1),main="")
pacf(ts(dltraf_12,frequency = 1),main="")
par(mfrow=c(1,1))
#ACF, PACF: mostly non significant coefficient
#except at lag1, lag11, lag 12 and lag13
test <- acf(ts(dltraf_12,frequency = 1),main="")
test

#------------------------------------------------
#step 1: identification of orders p,P, q, Q, d, D
#------------------------------------------------

#fit a multiplicative ARIMA on ltraf
#d=D=1: transformation performed on ltraf

#q, Q: MA part (ACF)

#p,P: AR part (PACF)


#-------------------------------
#step 2: estimation of the model
#-------------------------------

#TSA package: fittef; Mcleod.Li.test
library(TSA)
mod1 <- arima(ltraf,c(1,1,1),seasonal = list(order=c(1,1,1), period=12), method='ML')
mod1
#P=Q=p=q=D=d=1
#AIC= -850.82

mod2 <- arima(ltraf, c(1,1,1), seasonal = list(order=c(0,1,1),period=12), method='ML')
mod2
#P=0(to remove the sar1), Q=p=q=D=d=1
#AIC= -852.78, the second model is a little bit better

## residuals analysis ##
checkresiduals(mod1)
res1 <- mod1$residuals
plot(res1)

## residuals analysis ##
checkresiduals(mod2)
res2 <- mod2$residuals
plot(res2)


mod3 <- arima(ltraf, c(0,1,1), seasonal = list(order=c(0,1,1),period=12), method='ML')
#P=p=0(to remove the sar1 and ar1), Q=q=D=d=1
#AIC= -851.99, the third model is a little bit worse
mod3

## residuals analysis ##
checkresiduals(mod3)
res3 <- mod3$residuals
plot(res3)

#use the forecast function
auto.arima(ltraf,d=1,D=1,stepwise = FALSE, approximation = FALSE, trace = TRUE)
#p=P=0, q=3,d=D=Q=1
mod4 <- stats::arima(ltraf, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), method='ML')
mod4
#AIC= -860.1,this one is the best


#plot og the fitted value
library('stats')

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

#pvalue of the student test
pvalue <- 2*(1-pnorm(abs(tstat)))
tstat
pvalue
#ma1, ma3 and sma1 are significant (pvalue < 5%)

## residuals analysis ##
checkresiduals(mod4)
res1 <- mod4$residuals
plot(res1)

#ACF of the residuals (white noise assumption)
acf(ts(res1, frequency = 1))
# significant coefficients at lags 11

#PACF of the residuals
pacf(ts(res1, frequency = 1))
# significant coefficients at lags 11

#Ljung-Box test to check for the significance of ACF
Box.test(res1, lag = 20, type = "Ljung-Box")
# pvalue = 63% >5% so we accept HO: white noise assumption
# means the model is free of auto-correlation
# Heteroscedasticity
sq.res <- (res1)^2
acf(ts(sq.res,frequency=1))

#Normality assumption

#standardized residuals (residuals divided by the standard deviation)
res_stand <- res1/sqrt(mod4$sigma2)
summary(res_stand)

# If the normality assumption is satisfied, the standardized results should lie in between -2 and 2 (with 95% of chance)
plot(res_stand)
abline(a=2, b=0, col="red")
abline(a=-2, b=0, col="red")


# we can identify one outlier in the dataset
#it corresponds to the lowest value of the residuals
min(res_stand)
out1 <- which(res_stand < -5.5)
out1

library("zoo")
index(res_stand)[out1] #date of the outlier
# the date of the outlier is April 2010

res_stand[out1]# value of the outlier
traf[out1]
traf[(out1-4):(out1+4)]
traf[(out1-16):(out1-8)]
#In April 2011, the value of traffic was similar to May 2011, it is not at all the case in 2010

# 2 options: change the value to a more plausible one or add a dummy in the model


#QQ-plot
qqnorm(res1)
qqline(res1)

#Shapiro test
shapiro.test(res1)
#pvalue <<<< 5%, we reject normality (HO: normal dist)

#or Jarque.bera test 
library(tseries)
jarque.bera.test(res1)
#pvalue <<<< 5%, we reject normality (HO: normal dist)

#constant variance
sq.res <- (res1)^2 
acf(ts(sq.res, frequency = 1))

#The first correlation is significant so there seems to be a non constant variance.
#It may be due to the outliers we have identified on the graph of the residuals
#Other option: fit a GRACH model on the variance!


#-------------------------------------------
# Step2-bis: estimation of ARIMAX
#-------------------------------------------

# Create a dummy variable
CdgTraf$dum_1 <- 0
CdgTraf$dum_1[out1] <- 1
?arimax

mod5 <- stats::arima(ltraf, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), method='ML', xreg = CdgTraf$dum_1)
mod5
mod4

# Plot of the fitted value
fit2 <- fitted(mod5)
plot.ts(cbind(ltraf,fit2), plot.type = 'single', col=c('black','red'))


#-----------------------------------
# step3_bis: Validation of the model
#-----------------------------------

## Significance of the coefficients ##
mod5$coef
mod5$var.coef
tstat2 <- mod5$coef/sqrt(diag(mod5$var.coef))
tstat2
pvalue2 <- 2*(1-pnorm(abs(tstat2)))
pvalue2
pvalue
#except ma2, they are all significant (pvalue < 5%)

## Residuals analysis ##
res2 <- mod5$residuals
plot(res2)

# Autocorrelation of the residuals (White noise assumption)
acf(ts(res2, frequency=1))

pacf(ts(res2, frequency=1))
Box.test(res2,lag=20,type="Ljung-Box")
# pvalue = 66.5% >5% so we accept HO: white noise assumption


#Normality assumption

#standardized residuals (residuals divided by the standard deviation)
res_stand2 <- res2/sqrt(mod5$sigma2)

summary(res_stand2)

plot(res_stand2) 
abline(a=2,b=0,col="red")
abline(a=-2,b=0,col="red")

#QQ-plot
qqnorm(res2)
qqline(res2)

#Shapiro test
shapiro.test(res2)
#pvalue <<<< 5%, we reject normality (HO: normal dist)

#or Jarque.bera test 
library(tseries)
jarque.bera.test(res2)
#pvalue <<<< 5%, we reject normality (HO: normal dist)


#------------------
# Step4: Prediction
#------------------

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

#lty=line type


##We will use the data from 2000 to 2016 to predict the traffic from 2017 to 2019
#and then to compare our results with the real traffic data for these three years 
#in order to test the accuracy of our model
#so we just follow the same steps as before

#Import data
CdgTraf_test <- read.csv("CDG.csv")
CdgTraf_test <- CdgTraf_test[1:204,]
traf_test <- ts(CdgTraf_test[,3], start=c(2000,1), frequency=12)

#------------------------
#step 0: stationary issue
#------------------------
plot.ts(traf_test)
x_test <- decompose(traf_test)
plot(x_test)
par(mfrow=c(2,1))
acf(ts(traf_test,frequency = 1))
pacf(ts(traf_test,frequency = 1)) 

#log transformation: to reduce the increasing variance effect
ltraf_test <- log(traf_test)
plot.ts(ltraf_test)

#1st order transformation: to remove the non constant trend
dltraf_test <- diff(ltraf_test,lag=1)
plot(dltraf_test)

#difference of lag 12: to remove seasonality
dltraf_12_test <- diff(dltraf_test, lag=12)
plot(dltraf_12_test)
acf(ts(dltraf_12_test,frequency = 1),main="")
pacf(ts(dltraf_12_test,frequency = 1),main="")

#-------------------------------
#step 2: estimation of the model
#-------------------------------
#use the forecast function
auto.arima(ltraf_test,d=1,D=1,stepwise = FALSE, approximation = FALSE, trace = FALSE)
#get the same result as before

mod6 <- stats::arima(ltraf_test, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), method='ML')
mod6 
#AIC=-696.6


#plot og the fitted value
library('stats')

fit1_test <- fitted(mod6)
plot.ts(cbind(ltraf_test,fit1_test), plot.type = 'single',col=c('black','red'))
plot.ts(cbind(traf_test,2.718^fit1_test), plot.type = 'single',col=c('black','red'))

#-------------------------------
#step 3: validation of the model
#-------------------------------
## significance of the coefficients ##
mod6
mod6$coef
mod6$var.coef

# student test statistic
tstat_test <- mod6$coef/sqrt(diag(mod6$var.coef))

#pvalue of the student test
pvalue_test <- 2*(1-pnorm(abs(tstat_test)))
tstat_test
pvalue_test

## residuals analysis ##
checkresiduals(mod6)
res3 <- mod6$residuals
plot(res3)
acf(ts(res3, frequency = 1))
pacf(ts(res3, frequency = 1))

#Ljung-Box test to check for the significance of ACF
Box.test(res3, lag = 20, type = "Ljung-Box")
# pvalue = 68.61% >5% so we accept HO: white noise assumption

#Normality assumption

#standardized residuals (residuals divided by the standard deviation)
res_stand_test <- res3/sqrt(mod6$sigma2)
summary(res_stand_test)

# If the normality assumption is satisfied, the standardized results should lie in between -2 and 2 (with 95% of chance)
plot(res_stand_test)
abline(a=2, b=0, col="red")
abline(a=-2, b=0, col="red")

#QQ-plot
qqnorm(res3)
qqline(res3)

#Shapiro test
shapiro.test(res3)
#pvalue <<<< 5%, we reject normality (HO: normal dist)

#or Jarque.bera test 
library(tseries)
jarque.bera.test(res3)
#pvalue <<<< 5%, we reject normality (HO: normal dist)

#------------------
# Step4: Prediction
#------------------

# Quality of the fit
cb80_test <- mod6$sigma2^.5*qnorm(0.9)
plot(cbind(ltraf_test,fit1_test-cb80_test,fit1_test+cb80_test),plot.type='single',lty=c(1,2,2))

# Proportion of points in the confidence bound
indi_test <- (ltraf_test-(fit1_test-cb80_test))>0&(fit1_test+cb80_test-ltraf_test)>0
prop.test<- 100*sum(indi_test)/length(indi_test)
prop.test
# If prop >80%, then the fit is considered good.(88.24%)

pred_test <- predict(mod6, n.ahead = 36)
pred_test

# Prediction 
ts.plot(traf_test,xlim=c(2010,2019.500))
lines(2.718^(pred_test$pred), col="red")
lines(2.718^(pred_test$pred-1.96*pred_test$se),col=4,lty=2)
lines(2.718^(pred_test$pred+1.96*pred_test$se),col=4,lty=2)

ts.plot(traf,xlim=c(2010,2019.500))
lines(2.718^(pred_test$pred), col="red")
lines(2.718^(pred_test$pred-1.96*pred_test$se),col=4,lty=2)
lines(2.718^(pred_test$pred+1.96*pred_test$se),col=4,lty=2)

#comparison between the prediction and the real data (2016~2019)
ts.plot(traf,2.718^pred_test$pred, log = "y", lty = c(1,3))
ts.plot(traf,2.718^pred_test$pred, log = "y", lty = c(1,3),xlim=c(2016,2019.500))



####################################################
# Dummy variables - Rejected as it does not improve the model
####################################################


#Dummy variables were tried to improve the model as there are outliers in our residuals  

## residuals analysis ##
checkresiduals(mod4)
res1 <- mod4$residuals
plot(res1)

# If the normality assumption is satisfied, the standardized results should lie in between -2 and 2 (with 95% of chance)
plot(res_stand)
abline(a=2, b=0, col="red")
abline(a=-2, b=0, col="red")


# we can identify one outlier in the dataset
#it corresponds to the lowest value of the residuals
min(res_stand)
out1 <- which(res_stand < -5.5)
out1

res_stand[out1]# value of the outlier
traf[out1]
traf[(out1-4):(out1+4)]
traf[(out1-16):(out1-8)]
#In April 2011, the value of traffic was similar to May 2011, it is not at all the case in 2010


# Create a dummy variable
CdgTraf$dum_1 <- 0
CdgTraf$dum_1[out1] <- 1

mod5 <- stats::arima(ltraf, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), method='ML', xreg = CdgTraf$dum_1)
mod5
#mod5 AIC = -915.89

#Check residuals again with largest outlier taken out
checkresiduals(mod5)
res3 <- mod5$residuals
plot(res3)

## Prediction model is tested against the base model that has no dummy variables
#Update fit model with improved model
fit4 <- fitted(mod5)
plot.ts(cbind(ltraf,fit4), plot.type = 'single',col=c('black','red'))
plot.ts(cbind(traf,2.718^fit4), plot.type = 'single',col=c('black','red'))


# Quality of the fit

cb80 <- mod5$sigma2^.5*qnorm(0.9)
plot(cbind(ltraf,fit4-cb80,fit4+cb80),plot.type='single',lty=c(1,2,2))

# Proportion of points in the confidence bound
indi <- (ltraf-(fit4-cb80))>0&(fit4+cb80-ltraf)>0
prop <- 100*sum(indi)/length(indi)
prop
# Model with one dummy variable is proportion with quality of fit of 87.5%, which is 
# worse then the model of no dummy variables.

# A second model was created to test if two dummy variables in place of the 
# two largest outliers would improve the accuracy of the model.

#Residuals with worst outlier taken out
checkresiduals(mod5)
res3 <- mod5$residuals
plot(res3)


#standardized residuals (residuals divided by the standard deviation)
res_stand_3 <- res3/sqrt(mod5$sigma2)
summary(res_stand_3)
abline(a=2, b=0, col="red")
abline(a=-2, b=0, col="red")
# Around 2014 to 2015 there is another outlier that corresponds to the new 
# lowest value of the residuals

# Out2 given the value of the current outlier
min(res_stand_3)
out2 <- which(res_stand_3 < -5.5)
out2

# Dummy variable created for second outlier
CdgTraf$dum_2 <- 0
CdgTraf$dum_2[out2] <- 1

# mod7 given the arima function with the two dummy variables included
mod7 <- stats::arima(ltraf, c(0,1,3), seasonal = list(order=c(0,1,1),period=12), 
                     method='ML', xreg = cbind(CdgTraf$dum_1, CdgTraf$dum_2))
mod7

# Residuals for new model checked
checkresiduals(mod7)
res5 <- mod7$residuals
plot(res5)

#standardized residuals (residuals divided by the standard deviation)
res_stand_4 <- res5/sqrt(mod7$sigma2)
summary(res_stand_4)

# Plot of the new fitted value with second and first largest outliers removed
fit5 <- fitted(mod7)
plot.ts(cbind(ltraf,fit5), plot.type = 'single', col=c('black','red'))

#Plot of the two models overlapped, red line is the largest outlier removed, green is the 
#second largest and first largest removed. Improvements to the model are almost trival,
#therefore no more outliers are removed to not overfit the model to the data.

par(mfrow=c(1,1))
plot.ts(cbind(ltraf, fit4, fit5), plot.type = 'single', col=c('black','red', 'green'))


# Quality of the fit

cb80 <- mod7$sigma2^.5*qnorm(0.9)
plot(cbind(ltraf,fit5-cb80,fit5+cb80),plot.type='single',lty=c(1,2,2))

# Proportion of points in the confidence bound
indi <- (ltraf-(fit5-cb80))>0&(fit5+cb80-ltraf)>0
prop <- 100*sum(indi)/length(indi)
prop

# Model with two dummy variables proportion of points is worse, at 82.92%
# This shows that the model with no dummy variables should be used.

