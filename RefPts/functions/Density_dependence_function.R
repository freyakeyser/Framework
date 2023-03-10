# This little function is used to get the density dependence relationship between a response variable of some time and some metric of biomass
#Arguments explained
#dat = prod.dat, 
#type = 'lm',
#log.trans = T,
#response = 'RPS',
#dd.term = "SSB",
#min.cov.pred = 0,
#max.cov.pred = 20000,
# step = 100


dens.function <- function(dat = prod.dat, type = 'lm',log.trans = T,response = 'RPS',dd.term = "SSB",min.cov.pred = 0,max.cov.pred = 20000,step = 100)
{
 #browser()
 # Get a prediction object
  pred.dat <- data.frame(seq(min.cov.pred,max.cov.pred,by=step)) # Predict the RPS for each 100 tonnes of SSB
  names(pred.dat) <- 'dd'
  # Now pick the correct data from the dat object for the analysis...
  res <- dat[,which(names(dat) %in% response)]
  dd <- dat[,which(names(dat) %in% dd.term)]
  #datr <- data.frame(res = res, dd= dd)
  #names(datr) <- c("res",dd.term)
  # Get the model result you want
  if(log.trans == T)
  {
    if(type == 'glm') mod.res <- glm(res~dd, family = Gamma(link = "log"))
    if(type == 'lm')  mod.res <- lm(log(res)~dd)
    if(type == 'gam') mod.res <- gam(log(res)~s(dd))
  } else {
    if(type == 'glm') mod.res <- glm(res~dd, family = gaussian(link = "identity"))
    if(type == 'lm')  mod.res <- lm(res~dd)
    if(type == 'gam') mod.res <- gam(res~s(dd))
  } # end if/else log.trans == T
  
  # Now get the model predictions we need
  if(log.trans == T)
  {
    if(type == 'glm')  
    {
      pred.dat$y.log <- predict(mod.res,newdata = pred.dat,type = 'link')
      pred.dat$se <- predict(mod.res,newdata = pred.dat,type = 'link',se.fit=T)$se.fit
    } else { # end if(type == 'glm')  
      pred.dat$y.log <- predict(mod.res,newdata = pred.dat)                     
      pred.dat$se <- predict(mod.res,newdata = pred.dat,se.fit=T)$se.fit
    } # end the  type == glm if/else combo
    # Now get the UCI's and the y on the response scale
    pred.dat$y.res <- exp(pred.dat$y.log)
    pred.dat$uci.log <- pred.dat$y.log + 1.96*pred.dat$se
    pred.dat$lci.log <- pred.dat$y.log - 1.96*pred.dat$se
    pred.dat$uci <- exp(pred.dat$uci.log)
    pred.dat$lci <- exp(pred.dat$lci.log)
  } else { # end the if log.trans == T
    if(type == 'glm')  
    {
      pred.dat$y.res <- predict(mod.res,newdata = pred.dat,type = 'response')
      pred.dat$se <- predict(mod.res,newdata = pred.dat,type = 'response',se.fit=T)$se.fit
    } else { # end the if type = 'glm' start the else
      pred.dat$y.res <- predict(mod.res,newdata = pred.dat)                     
      pred.dat$se <- predict(mod.res,newdata = pred.dat,se.fit=T)$se.fit
    } # end the if/else glm
    pred.dat$uci <- pred.dat$y.res + 1.96*pred.dat$se
    pred.dat$lci <- pred.dat$y.res - 1.96*pred.dat$se
  } # end else log.trans == T
  
  # Now rename the pred.dat columns so we know what they are...
  names(pred.dat)[which(names(pred.dat) == 'dd')] <- dd.term
  names(pred.dat)[which(names(pred.dat) == 'y.res')] <- response
  if(log.trans == T) names(pred.dat)[which(names(pred.dat) == 'y.log')] <- paste0(response,'.log')
  
  return(list(pred.dat = pred.dat,mod.res=mod.res))
} # end dens.function
