load('C:\\Users\\oscar\\Desktop\\STAT 631\\Final2022.Rdata')
source('C:\\Users\\oscar\\Desktop\\STAT 631\\GARCH_RFunctions.R')
source('C:\\Users\\oscar\\Desktop\\STAT 631\\GARCH_plotFunctions.R')
library(quantmod); library(rugarch);library(forecast);library(urca)

cat("Starting at:");
head(Yt,2)

cat("Ending at:");
tail(Yt,2)

par(mfrow = c(1,2))
plot(Yt, yaxis.right = F);plot(diff(Yt), yaxis.right = F, main = "differenced Yt")
# From the plots we can observe that Yt has an upward trend and it seems to be a non stationary random walk process. After taking the difference we can observe that the series has a constant mean, but its variance is not. This suggests a log transformation.  

Xt = log(Yt);dXt = diff(Xt)[-1,]; par(mfrow = c(1,2))
plot(Xt, yaxis.right = F, main = "log Yt");plot(diff(Xt), yaxis.right = F,
                                                   main = "differenced log Yt")

# It seems to be stationary after log transforming Yt and taking the difference

# Let's check if Xt is a unit root process
# ADF unit root test on Xt

n = dim(Xt)[1]; p.max = round(12*(n/100)^.25)
adf = ur.df(Xt,type = "drift", lags = p.max, selectlags = "AIC" );
cat("ADF test with maximum lag", p.max, "and lag selected by AIC ", adf@testreg$df[1]-2,
       "\nTest Statistic = ", adf@teststat[1],"\nCritical values");
adf@cval[1,];


cat("\nReject or not Reject the null:");
adf@teststat[1] < adf@cval[1,]

# The ADF test with lag parameter p = 3 fails to reject the null of a unit root process, confirming that Xt is a unit root process.

par(mfrow = c(1,2))
acf(dXt);acf(dXt^2);


test1 = sapply(4:8,function(u) Box.test(dXt,u, type = "L")[c(1,3)])
test2 = sapply(4:8,function(u) Box.test(dXt^2,u, type = "L")[c(1,3)])
colnames(test1) = colnames(test2) = paste("df=",4:8)
cat("Ljung-Box test of dXt");
test1


cat("Ljung-Box test of squared dXt");
test2

# We conclude that Xt is of ARIMA(p, 1, q) class without GARCH effects. 
# We start with a model selected by auto.arima() in which we set d = 1.

aic.X = auto.arima(Xt, d = 1, ic = "aic")
aic.X

plot_Roots(coef(aic.X))

# Both AIC and BIC selected the same model, ARIMA(1,1,1). Both coefficients are significant, the AR and MA roots are not close numerically or graphically. The selected model does not have parameter redundancy.   

# Let's check the residuals of ARIMA(1,1,1)

res = aic.X$residuals
restest = sapply(4:8,function(u) Box.test(res,u + 2, type = "L", 2)[c(1,2,3)])
colnames(restest) = paste("lags =",(4:8) + 2);
restest

# The residuals are not correlated
# Finding a suitable distribution for the errors (white noise) in the selected model and checking the fit with a Quantile-Quantile plot.

par(mfrow = c(1,3), pty = "s");
res = resid(aic.X)
acf(res, main = "ACF of res from ARIMA(1,1,1)")
hist(res, main = "Histogram of res from ARIMA(1,1,1)")
qqnorm(res, main = "res from ARIMA(1,1,1)", xlab = "Normal quantile")
qqline(res)


dists = c("snorm", "snorm", "sstd","sged", "jsu", "nig")
AIC = BIC = Inf; aic.dist = bic.dist = c()
for(i in 1:length(dists)){
   out = fitdist(dists[i],res)
   aic = 2*tail(out$val,1) + 2*length(out$par)
   bic = 2*tail(out$val,1) + log(n)*length(out$par)
   if(aic < AIC) {AIC = aic; aic.dist = list(dist = dists[i], coef = out$par)}
   if(bic < BIC) {BIC = bic; bic.dist = list(dist = dists[i], coef = out$par)}
   }
cat("AIC selects", paste(aic.dist$dist), "\t\tBIC selects", paste(bic.dist$dist));

# Both AIC and BIC select Johnson's SU distribution for the white noise of ARIMA(1,1,1) for Xt. The QQ
# Plot with the residuals from ARIMA(1,1,1) shows the selected distribution fits well. The final model for Xt = log Yt is ARIMA(1,1,1) with i.i.d. white noise from Johnson's SU distribution.

dist = aic.dist$dist; coef = aic.dist$coef
ps = ((1:n)-1/2)/n
xs = quantile(res,ps)
ys = qdist(dist,ps,mu = coef[1],sigma = coef[2], skew = coef[3], shape = coef[4])

par(mfrow = c(1,2) )
plot(xs,ys, xlab = "sample quantiles", ylab = "jsu quantiles" )
abline(lsfit(xs,ys)$coef)

# Plot of the forecast (6 month ahead)

plot(forecast(aic.X, h = 6, level = c(68,95)), xlim = c(131,(n+6)))


# 6 month forecast with 95% prediction intervals

qs = qdist(dist,c(.025, .0975), sigma = coef[2], skew = coef[3], shape = coef[4])
pred = predict(aic.X, n.ahead = 6)
pi = sapply(1:2, function(u) pred$pred + qs[u]*pred$se)
pi = cbind(pred$pred, pi)
colnames(pi) = c("point pred", "low95", "hi95")
cat("1-6 step ahead for Xt");pi

cat("1-6 step ahead for Yt");exp(pi)
























