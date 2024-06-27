# This little function is used to get the density dependence relationship between a response variable of some time and some metric of biomass
#Arguments explained

# B:  Biomass (or whatever you are comparing against) 
# dat: The data you are sampling from
# type = 'sample',
# bp = 0,
# sd:  Do you want to specify the standard deviation or select it based on the data (default behaviour picks sd from data, which is sd=NULL)
# max.B:    If the Biomass is > max.B we reduce the productivity of the variable of interest effectively to 0.
# mn.at.max:
# sd.at.max: 

bp.function <- function(B = 0, dat = prod.dat$g,ssb.ts= prod.dat$SSB, type = 'sample',sd = NULL,bp = 0,max.B =1e99,mn.at.max = 1,sd.at.max = 0.05)
{
  #browser() 
  # Split the data into low and high values, not the most computationally efficient way to do this as I'm making this happen a billion times... but makes code tidy...
  low.dat <- na.omit(dat[ssb.ts < bp])
  hi.dat <- na.omit(dat[ssb.ts >= bp])
  if(B >=bp) 
  { 
    if(type == 'dist')
    {
      if(is.null(sd)) res<-rlnorm(1,mean=log(median(hi.dat,na.rm=T)),sd=sd(log(hi.dat),na.rm=T))
      if(!is.null(sd)) res<-rlnorm(1,mean=log(median(hi.dat,na.rm=T)),sd=sd)
    }
    if(type == 'sample') res<-sample(hi.dat,1)
  } # end if(B >= bp) 
  #browser()
  # Now when biomass is smaller than the breakpoint..
  if(B < bp) 
  { 
    if(type== 'dist')
    {
      if(is.null(sd)) res<-rlnorm(1,mean=log(median(low.dat,na.rm=T)),sd(log(low.dat),na.rm=T))
      if(!is.null(sd)) res<-rlnorm(1,mean=log(median(low.dat,na.rm=T)),sd=sd)
    }
    if(type == 'sample') res<-sample(low.dat,1)
    
  } # end if(B < bp)
  
  # Finally if B is higher than anything we ever saw in the time series we can 
  # make the parameter be whatever we want with whatever variance we want.
  # Idea here to to make productivity low at super high biomasses.
  if(B > max.B) res <-rlnorm(1,log(mn.at.max),sd=sd.at.max)
#browser()
  return(res)
} # end function