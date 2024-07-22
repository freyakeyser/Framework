# Here's a correlation function where we can make a nice time series for 2 correlated (or uncorrelated) terms.  This will enable 
# the inclusion of autocorrelation in the time series (only Auto-regressive terms for the moment) along with cross correlation between the time series
# and lags between the two time series.  All of this is easily obtained by looking at your input data.
# Arguments...


cor.fun <- function(ts1 = list(mn = 0.5,var=0.1),
                    ts2 = list(mn =2,var=0.1),
                    arima = list(proc = 'arima',cross.cor = 0,lag = 5,
                                 ts1.ar1=0,ts1.ar2=-0,
                                 ts2.ar1=0, ts2.ar2 =0),
                    ts = list(years=1000,start.year =200, final.ts.len = 20))
  
{
  #browser()

  # This gets the covariance strcutre for the 'error'.  Note I put it on a log scale.
  # Get the correlation matrix   
  corr.mat <- matrix(c(1,arima$cross.cor,arima$cross.cor,1),nrow = 2)
  # And convert to a covariance matrix.
  cov.mat <- corr.mat * sqrt(ts1$var) * sqrt(ts2$var)
  n <- ts$years
  # This gets us two correlated 'error' time series
  eps <- mvtnorm::rmvnorm(n = n, mean = log(c(ts1$mn,ts2$mn)), sigma = cov.mat)
  
  # So now that we have our time series generated we turn them into the ARIMA time series we want...
  if(arima$proc == 'normal') 
  {
    ts1.ts <- arima.sim(model =list(0,0,0),n = n, innov = eps[,1])
    ts2.ts <- arima.sim(n = n, model = list(0,0,0), innov = eps[,2])  
  }
  # If you wanted it to be an ARIMA model with specfied p/q/d terms
  if(arima$proc == 'arima')  
  {
    ts1.ts <- arima.sim(model =list(ar = c(arima$ts1.ar1,arima$ts1.ar2)),n = n, innov = eps[,1])
    ts2.ts <- arima.sim(model =list(ar = c(arima$ts2.ar1,arima$ts2.ar2)),n = n, innov = eps[,2]) 
  }
  # If you wanted it to be a random walk model (p=1)
  if(arima$proc == 'rand_walk')
  {
    ts1.ts <- arima.sim(model =list(1,0,0),n = n, innov = eps[,1])
    ts2.ts <- arima.sim(n = n, model = list(1,0,0), innov = eps[,2]) 
  }
  # # Save the time series (f input, ts2 output) - Mod 1
  # dat <- data.frame(year = ts$start.year:(ts$start.year+ts$final.ts.len-1),
  #                   ts1.ts = exp(ts1.ts[ts$start.year:(ts$start.year+ts$final.ts.len-1)]),
  #                   ts2.ts = exp(ts2.ts[(ts$start.year-arima$lag):(ts$start.year+ts$final.ts.len-1-arima$lag)]))
  # # Save the time series (ts2 input, f output) - Mod 2
  # Now when we add correlation to these the values decrease, so I want to rescale these time series
  # to the mean value of the timeseries....
  ts1.correct.factor <- ts1$mn/median(exp(ts1.ts[(ts$start.year+arima$lag):(ts$start.year+ts$final.ts.len-1+arima$lag)]),na.rm=T)
  ts2.correct.factor <- ts2$mn/median(exp(ts2.ts[ts$start.year:(ts$start.year+ts$final.ts.len-1)]),na.rm=T)
  dat <- data.frame(year = ts$start.year:(ts$start.year+ts$final.ts.len-1),
                    ts1.ts = ts1.correct.factor * exp(ts1.ts[(ts$start.year+arima$lag):(ts$start.year+ts$final.ts.len-1+arima$lag)]),
                    ts2.ts = ts2.correct.factor * exp(ts2.ts[ts$start.year:(ts$start.year+ts$final.ts.len-1)]))
  
  
  return(dat)
} # end function