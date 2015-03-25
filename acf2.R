acf2=function(series,max.lag=NULL,ma.test=FALSE){
  num=length(series)
  if (is.null(max.lag)) max.lag=ceiling(10+sqrt(num))
  if (max.lag > (num-1)) stop("Number of lags exceeds number of observations")
  ACF=stats::acf(series, max.lag, plot=FALSE)$acf[-1]
  PACF=pacf(series, max.lag, plot=FALSE)$acf
  LAG=1:max.lag/frequency(series)
  minA=min(ACF)
  minP=min(PACF)
  U=1.96/sqrt(num)
  L=-U
  nu<-1.96*(c(1,1+2*cumsum(ACF[-max.lag]^2))/num)^.5    # 5% bandwidths for ma.test
  minu=min(minA,minP,L,-nu)-.01    # now it includes also nu
  old.par <- par(no.readonly = TRUE)
#  par(mfrow=c(2,1), mar = c(3,3,2,0.8), oma = c(1,1.2,1,1), mgp = c(1.5,0.6,0)) replaced by layout
  par(mar = c(3,3,2,0.8), oma = c(1,1.2,1,1), mgp = c(1.5,0.6,0))
  layout(1:2)
  plot(LAG, ACF, type="h",ylim=c(minu,1), main=paste("Series: ",deparse(substitute(series))))
    abline(h=0, lty=1, col=1)
    if(ma.test==FALSE) abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,4,4))
    if(ma.test==TRUE) {
      par(new=TRUE)
      plot(LAG, nu, type="l",ylim=c(minu,1), lty=2,col=4,ylab="")
      par(new=TRUE)
      plot(LAG, -nu, type="l",ylim=c(minu,1), lty=2,col=4,ylab="")
    }
  plot(LAG, PACF, type="h",ylim=c(minu,1))
    abline(h=c(0,L,U), lty=c(1,2,2), col=c(1,4,4))
  on.exit(par(old.par))  
#  ACF<-round(ACF,2); PACF<-round(PACF,2)    
  return(cbind(ACF, PACF)) 
  }
