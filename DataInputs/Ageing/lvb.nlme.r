lvb.nlme <- function(AGE.dat,random.effect="year",ran.par='both',verbose=T, GBmodel=F, k.par='estimate'){


  # Estimate LVB parameters using a
  # non-Linear mixed-effects model
browser()
  AGE.dat <- as.data.frame(AGE.dat)
  require(nlme)
  AGE.dat$raneff<-AGE.dat[,random.effect]
  ran.effects<-unique(AGE.dat$raneff)

  AGE.gdat <- groupedData(HEIGHT ~ AGE | raneff, data = AGE.dat)
  linf <- c()
  k <- c()

  if(ran.par=='both') ran.form = linf + k ~ 1
  if(ran.par=='linf') ran.form = linf ~ 1
  if(ran.par=='k') ran.form = k ~ 1

  if(k.par=='estimate') AGE.nlme <- nlme(model = HEIGHT ~ linf * (1 - exp(-k * (AGE - tzero))), data = AGE.gdat, fixed = linf + k + tzero ~ 1, random = ran.form, start = c(150, 0.3, 0), na.action = na.omit, control = list(maxIter = 50), method = "ML")
  if(k.par!='estimate'){
    AGE.gdat$k<-k.par
    AGE.nlme <- nlme(model = HEIGHT ~ linf * (1 - exp(-k * (AGE - tzero))), data = AGE.gdat, fixed = linf +  tzero ~ 1, random = ran.form, start = c(150, 0), na.action = na.omit, control = list(maxIter = 50), method = "ML")
  }
  if(is.character(AGE.dat[,random.effect]))fit <- data.frame(raneff=row.names(coef(AGE.nlme)),coef(AGE.nlme))
  if(!is.character(AGE.dat[,random.effect]))fit <- data.frame(raneff=sort(row.names(coef(AGE.nlme))),linf=coef(AGE.nlme)[order(as.numeric(row.names(coef(AGE.nlme)))),1],k=coef(AGE.nlme)[order(as.numeric(row.names(coef(AGE.nlme)))),2])
  for(i in 1:length(ran.effects)){
    linf[i] <- exp(fit[i,2])
    k[i] <- fit[i,3]
  }
  Linf <- AGE.nlme$coef$fixed[1]
  K <- AGE.nlme$coef$fixed[2]
  t0 <- AGE.nlme$coef$fixed[3]

  names(fit)[1]<-random.effect

  if(verbose) print(summary(AGE.nlme))
  AGE.dat$label<-AGE.dat[,random.effect]
  if(random.effect=="month")AGE.dat$label<-months(as.Date(paste("2009-",1:12,"-01",sep="")))[AGE.dat$raneff]
  return(list(Linf=Linf,K=K,linf=linf,k=k,t0=t0,data=AGE.dat,fit=fit))

}