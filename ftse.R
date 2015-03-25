URL <- "http://ichart.finance.yahoo.com/table.csv?s=FTSEMIB.MI"
data <- read.csv(URL)
data$Date <- as.Date(strptime(data$Date, "%Y-%m-%d"))
head(data)
tail(data)
summary(data)
pt <- as.ts(data$Adj.Close, start = min(data$Date), end = max(data$Date))
rt <- as.ts(diff(log(data$Adj.Close)), start = min(data$Date), end = max(data$Date))
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = T))
plot(pt)
plot(log(pt))
plot(rt)
layout(matrix(c(1, 2), 1, byrow = T))
plot(abs(rt))
plot(rt^2)
round(c("Num of rt" = length(rt), "Min" = min(rt), "Max" = max(rt), "Mean" = mean(rt), "Median" = median(rt), "Variance" = var(rt), "Stdev" = sd(rt)), 6)

# Jarque-Bera normality test
library(tseries)
jarque.bera.test(pt) #reject H0: normality
jarque.bera.test(log(pt)) #reject H0
jarque.bera.test(rt) #reject H0


#ADF test for stationarity
library(urca)
urdrift <- function(x) {
  urtest<-ur.df(y, type = "drift", lags = x)
  c(stat = urtest@teststat[1],"5% crit. value" = urtest@cval[1, 2])
}
a <- 0:10
paste("ADF(", 1:10,")", sep = "")
names(a)<-c("DF", paste("ADF(", 1:10,")", sep = ""))
#ADF on pt
y <- pt
round(sapply(a, urdrift), 3) 
#ADF on log(pt)
y <- log(pt)
round(sapply(a, urdrift), 3)
#ADF on diff(log(pt)) = rt
y <- diff(log(pt))
round(sapply(a, urdrift), 3)
y <- rt
round(sapply(a, urdrift), 3)

#we continue the analysis using the stationary series rt
#ACF and PACF of rt and rt^2
layout(matrix(c(1, 2, 3, 4), 2, byrow = T))
acf(rt)
acf(rt^2)
pacf(rt)
pacf(rt^2)
#ARMA model fit
library(TSA)
plot(armasubsets(rt, nar = 5, nma = 5))
(ar5 <- arima(rt, order = c(5, 0, 0)))
(ar5 <- arima(rt, order = c(5, 0, 0), include.mean = F))
(ar5 <- arima(rt, order = c(5, 0, 0), fixed = c(0, 0, NA, NA, NA), include.mean = F))
(arma5_5 <- arima(rt, order = c(5, 0, 5), include.mean = F))
(arma5_3 <- arima(rt, order = c(5, 0, 3), include.mean = F))
(arma5_3 <- arima(rt, order = c(5, 0, 3), include.mean = F, fixed = c(NA, NA, NA, 0, NA, NA, NA, NA)))

#residual analysis
plot(arma5_3$res)
plot(arma5_3$res^2)
acf(arma5_3$res^2)
pacf(arma5_3$res^2)
#serial correlation test
sapply(c(1, 5, 10, 20), function(i) Box.test(rt, lag = i, type = "Ljung"))
sapply(c(1, 5, 10, 20), function(i) Box.test(rt^2, lag = i, type = "Ljung"))
sapply(c(1, 5, 10, 20), function(i) Box.test(arma5_3$res, lag = i, type = "Ljung"))
sapply(c(1, 5, 10, 20), function(i) Box.test(arma5_3$res^2, lag = i, type = "Ljung"))
library(FinTS)
sapply(c(1, 5, 10, 20), function(i) ArchTest(arma5_3$res, lags = i, demean = F)) #we reject H0: no ARCH


#GARCH fit
library(fGarch)
summary(garchFit(formula = ~ arma (5, 3) + garch(1, 1), data = rt, trace = F))
summary(garchFit(formula = ~ arma (5, 3) + garch(1, 1), data = rt, trace = F, include.mean = F))
summary(garchFit(formula = ~ arma (2, 2) + garch(1, 1), data = rt, trace = F, include.mean = F))
summary(garchFit(formula = ~ arma (2, 2) + garch(2, 1), data = rt, trace = F, include.mean = F))
summary(garchFit(formula = ~ arma (2, 2) + garch(1, 1), data = rt, trace = F, include.mean = F, cond.dist = "std"))
summary(garchFit(formula = ~ arma (2, 2) + garch(2, 1), data = rt, trace = F, include.mean = F, cond.dist = "std"))

summary(garchFit(formula = ~ garch(2, 1), data = rt, trace = F))
summary(garchFit(formula = ~ garch(2, 1), data = rt, trace = F, cond.dist = "std"))
summary(garchFit(formula = ~ garch(2, 1), data = rt, trace = F, include.mean = F))
summary(garchFit(formula = ~ garch(2, 1), data = rt, trace = F, include.mean = F, cond.dist = "std"))
summary(garchFit(formula = ~ garch(1, 1), data = rt, trace = F))
summary(garchFit(formula = ~ garch(1, 1), data = rt, trace = F, cond.dist = "std"))
summary(garchFit(formula = ~ garch(1, 1), data = rt, trace = F, include.mean = F))
summary(garchFit(formula = ~ garch(1, 1), data = rt, trace = F, include.mean = F, cond.dist = "std"))

garch1_1 <- garchFit(formula = ~ garch(1, 1), data = rt, trace = F, cond.dist = "std")
garch2_1 <- garchFit(formula = ~ garch(2, 1), data = rt, trace = F, cond.dist = "std")
arma2_2garch2_1 <- garchFit(formula = ~ arma (2, 2) + garch(2, 1), data = rt, trace = F, include.mean = F, cond.dist = "std")
arma2_2garch1_1 <- garchFit(formula = ~ arma (2, 2) + garch(1, 1), data = rt, trace = F, include.mean = F, cond.dist = "std")
plot(arma2_2garch1_1)

layout(matrix(1:16, 4, 4))
sapply(1:13, function(i) plot(garch2_1, which=i))

#GARCH simulation
spec1 <- garchSpec(model = list(mu = garch1_1@fit$coef[1], omega = garch1_1@fit$coef[2], alpha = garch1_1@fit$coef[3], beta = garch1_1@fit$coef[4])) #garch(1,1)
spec2 <- garchSpec(model = list(mu = garch2_1@fit$coef[1], omega = garch2_1@fit$coef[2], alpha = c(garch2_1@fit$coef[3], garch2_1@fit$coef[4]), beta = garch2_1@fit$coef[5])) #garch(2,1)
spec3 <- garchSpec(model = list(ar = c(arma2_2garch2_1@fit$coef[1], arma2_2garch2_1@fit$coef[2]), ma = c(arma2_2garch2_1@fit$coef[3], arma2_2garch2_1@fit$coef[4]), omega = arma2_2garch2_1@fit$coef[5], alpha = c(arma2_2garch2_1@fit$coef[6], arma2_2garch2_1@fit$coef[7]), beta = arma2_2garch2_1@fit$coef[8])) #arma(2,2) + garch(2,1)
spec4 <- garchSpec(model = list(ar = c(arma2_2garch1_1@fit$coef[1], arma2_2garch1_1@fit$coef[2]), ma = c(arma2_2garch1_1@fit$coef[3], arma2_2garch1_1@fit$coef[4]), omega = arma2_2garch1_1@fit$coef[5], alpha = arma2_2garch1_1@fit$coef[6], beta = arma2_2garch1_1@fit$coef[7])) #arma(2,2) + garch(1,1)
sts1 <- garchSim(spec1, n = 4005)
sts2 <- garchSim(spec2, n = 4005)
sts3 <- garchSim(spec3, n = 4005)
sts4 <- garchSim(spec4, n = 4005)


layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = T))
plot(sts1, ylim = c(-0.15, 0.15), main = "GARCH(1, 1)")
par(new = T)
plot(rt, ylim = c(-0.15, 0.15), col = "blue", , xaxt = "n")
plot(sts2, ylim = c(-0.15, 0.15), main = "GARCH(2, 1)")
par(new = T)
plot(rt, ylim = c(-0.15, 0.15), col = "blue", , xaxt = "n")
plot(sts3, ylim = c(-0.15, 0.15), main = "ARMA(2, 2) + GARCH(2, 1)")
par(new = T)
plot(rt, ylim = c(-0.15, 0.15), col = "blue", , xaxt = "n")
plot(sts4, ylim = c(-0.15, 0.15), main = "ARMA(2, 2) + GARCH(1, 1)")
par(new = T)
plot(rt, ylim = c(-0.15, 0.15), col = "blue", , xaxt = "n")

#forecast
predict(garch1_1, 5)
predict(garch2_1, 5)
predict(arma2_2garch2_1, 5)
predict(arma2_2garch1_1, 5)