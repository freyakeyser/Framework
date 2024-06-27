# We are running a model for BBn using the Bayesian SS Model. This extracts the bits of code from Freya's scripts to make things work.

funs <- c("https://raw.githubusercontent.com/Mar-Scal/Assessment_fns/master/Maps/pectinid_projector_sf.R",
          "https://raw.githubusercontent.com/Mar-Scal/Assessment_fns/master/Fishery/logs_and_fishery_data.r",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Model/projections.r",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Model/decision.r",
          "https://raw.githubusercontent.com/Mar-Scal/Assessment_fns/master/Model/post.plt.R",
          "https://raw.githubusercontent.com/Mar-Scal/Assessment_fns/master/Model/exploit.plt.R",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Model/fit.plt.R",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Model/diag.plt.R",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Model/biomass.plt.R",
          "https://raw.githubusercontent.com/Mar-Scal/Assessment_fns/master/Maps/combo_shp.R",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Model/run_model.R",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Fishery/fishery.dat.r"
           )

for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}

library(sf)
library(sp)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(R2jags)
#library(SSModel)

theme_set(theme_few(base_size = 22))


direct <- "Y:/Offshore/Assessment/"
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"

load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90RSCS_newMWSH_GBb.RData")
# Get the survey boundary file...
temp <- tempfile()

# Download this to the temp directory you created above
#if(add_layer$survey[2] == 'detailed') download.file("https://raw.githubusercontent.com/Mar-scal/GIS_layers/master/offshore_survey_strata/offshore_survey_strata.zip", temp, quiet=quiet)
download.file("https://raw.githubusercontent.com/Mar-scal/GIS_layers/master/survey_boundaries/survey_boundaries.zip", temp, quiet=T)
# Figure out what this file was saved as
temp2 <- tempfile()
# Unzip it
unzip(zipfile=temp, exdir=temp2)
offshore.strata <- combo.shp(temp2,make.sf=T,make.polys=F, quiet=T)

BBn.strata <- offshore.strata %>% dplyr::filter(ID == "BBn")

# Grab the input data for the survey.

BBn.surv <- survey.obj$BBn[[1]]

# Using teh midpoint/mixed imputation method for missing years
#if(impute=="midpoint" | impute=="mixed") {
# We didn't have a survey in 2020 so need to impute data
year2020 <- as.data.frame(lapply(X = BBn.surv[BBn.surv$year %in% 2019:2021,], MARGIN = 2, mean))
year2020$year <- 2020
names(year2020) <- names(BBn.surv)
BBn.surv <- merge(BBn.surv, year2020, all=T)


years <- min(BBn.surv$year):max(BBn.surv$year)

logs_and_fish(loc="offshore",year = 1990:2022,direct=direct)
# If you get any NA's related warnings it may be something is being treated as a Factor in one of the two files.  
# This should combine without any warnings so don't ignore warnings here.
dat.fish<-merge(new.log.dat,old.log.dat,all=T)
dat.fish$ID<-1:nrow(dat.fish)
    
# # Grab the growth data, we have ageing data from 1980's that I'm going to use to calculate growth here.
# Data is coming from ageing data in 1989, found here.... Y:\Offshore\Assessment\Data\Ageing\archive\old_ageing_from_Amy_2022\BBn height at age 1989_2.pdf 
# L.inf <- 164.4
# #to <- 1.337 # So this uses a 1 year offset that we no longer believe in, going to make this 0.337 to align more with what we now do...
# to <- -0.2
# K <- 0.2


# Run this for one or both banks
mod.dat <- NULL
cpue.dat <- NULL
proj.dat <- NULL
# Now we need to calculate the growth for the models and we also extract the fishery data for the survey year here.  

# Get to get rid of all the crap spatial data in here.
fish.tmp <- dat.fish[dat.fish$bank == "BBn"  & !is.na(dat.fish$bank) & dat.fish$lon < 0 & dat.fish$lat > 0 & !is.na(dat.fish$lon) & !is.na(dat.fish$lat),]
fish.tmp <- st_as_sf(fish.tmp, coords = c('lon','lat'))
st_crs(fish.tmp) <- 4326
fish.dat <- st_intersection(BBn.strata,fish.tmp)

# BBn survey in May, will use same set up as BBn
cpue.dat <- fishery.dat(fish.dat,bk="BBn",yr=(min(years)-1):max(years),surv='May', method='jackknife',period = "survyr", direct=direct) 
# Combine the survey and Fishery data here.
mod.dat <- merge(BBn.surv,cpue.dat,by ="year")
# Get the CV for the CPUE...
mod.dat$U.cv <- mod.dat$cpue.se/mod.dat$cpue
# Now get the data for the projection part of the year  
proj.sub <- subset(fish.dat,year %in% years & months(as.Date(fish.dat$date)) %in% c("June","July","August","September","October","November","December"))
# Now calculate the fishery statistics for the projection period
proj.dat <- fishery.dat(proj.sub,bk="BBn",yr=(min(years)-1):max(years),method='jackknife', period = "calyr", direct=direct) 	
# So first up, this condition is the weighted mean condition, this uses the GAM predicted scallop condition factor for each tow
# and the biomass from each tow to come up with an overall bank average condition factor.
# This is weight in this year, which becomes t-1 
# waa.tm1 <- mod.dat$CF*(mod.dat$l.bar/100)^3
# # Using this years average shell height we can find the exptected shell height for the scallops in the next year
# # ht = (Linf * (1-exp(-K)) + exp(-K) * height(last year))
# # laa.t is the projected size of the current years scallops into next year.
# laa.t <- L.inf*(1-exp(-K)) + exp(-K) * mod.dat$l.bar
# # The c() term in the below offsets the condition so that current year's condition slots into the previous year and repeats 
# # the condition for the final year), this effectively lines up "next year's condition" with "predictied shell height next year (laa.t)
# # This gets us the predicted weight of the current crop of scallops next year based on next years CF * laa.t^3
# # Of course we don't have next years condition thus th last condition is simply repeated
# # waa.t is using the condition from next year and the growth from next year to get next years weight
# waa.t <- c(mod.dat$CF[-1],mod.dat$CF[nrow(mod.dat)])*(laa.t/100)^3
# # Here we use the current condition factor to calculate the weight next year (since we use laa.t)
# # That's really the only difference between waa.t and waa.t2, waa.t uses next years condition to project growth
# # what waa.t2 uses the current condition to project growth.  So that's really what we are comparing here with these
# # two growth metrics isn't it, this is really just comparing impact of using current vs. future condition factor on our growth estimates.
# waa.t2 <- mod.dat$CF*(laa.t/100)^3
# # Now the growth, expected and realized.
# mod.dat$g <- waa.t/waa.tm1
# # This is using the actual condition factor and growing the scallops by laa.t
# mod.dat$g2 <- waa.t2/waa.tm1
# 
# # same thing here but for the recruits
# waa.tm1 <- mod.dat$CF*(mod.dat$l.k/100)^3
# laa.t <- L.inf*(1-exp(-K))+exp(-K)*mod.dat$l.k
# waa.t <- c(mod.dat$CF[-1],mod.dat$CF[nrow(mod.dat)])*(laa.t/100)^3
# waa.t2 <- mod.dat$CF*(laa.t/100)^3
# mod.dat$gR <- waa.t/waa.tm1
# mod.dat$gR2 <- waa.t2/waa.tm1# setwd("C:/Assessment/2014/r")
#   
### overwrite imputation for growth here using whichever method
# in 2020, the covid-19 pandemic prevented the DFO survey from occurring. An industry-lead survey of limited scope occurred, but is not suitable for inclusion in the 
# assessment models. As such, we need to fill-in the blank row for 2020 with some data. We'll try out different options for doing that here. 
# We imputed the values in the survey data earlier, but for the "mixed" imputation method, we'll handle growth separately.
# Maybe it makes more sense to use the LTM for growth but midpoint for other values. 
    
# change 2019 and 2020 values to NA    

# mod.dat$g[which(mod.dat$year %in% c(2019:2020))] <- NA
# mod.dat$g2[which(mod.dat$year %in% c(2020))] <- NA
#   
# mod.dat$gR[which(mod.dat$year %in% 2019:2020)] <- NA
# mod.dat$gR2[which(mod.dat$year %in% 2020)] <- NA
#     
# # replace the NAs with long term medians
# mod.dat$g[which(mod.dat$year %in% c(2019:2020))] <- median(mod.dat$g, na.rm=T)
# mod.dat$g2[which(mod.dat$year %in% c(2020))] <- median(mod.dat$g2, na.rm=T)
#     
# mod.dat$gR[which(mod.dat$year %in% c(2019:2020))] <- median(mod.dat$gR, na.rm=T)
# mod.dat$gR2[which(mod.dat$year %in% c(2020))] <- median(mod.dat$gR2, na.rm=T)

# Quick for Freya
#write.csv(mod.dat,"D:/Framework/SFA_25_26_2024/Model/Data/BBn_input_data_for_freya.csv")


# Changing to a new method to calculate growth...
sizes <- seq(0.025,2,by=0.05) # So I'd be using the 1.025 bin and everything bigger for the t+1 fully-recruited
surv.years <- unique(surv.dat$BBn$year)
# But not 2020...
surv.years <- surv.years[surv.years != 2020]
# The w.yst object is exactly proportional to mod.dat$I, there is an offset, but given I need proportions I think this object is perfectly fine to use.
# SO this mw.per.bin is taking the stratified biomass and dividing it by the stratified numbers in each bin, which gives us the MW in that bin. 
# There is probably a MW object out there somewhere with this in it, but it should just be the same thing as this.
mw.per.bin <- data.frame(mw.per.bin = rbind(survey.obj$BBn$shf.dat$w.yst/survey.obj$BBn$shf.dat$n.yst,rep(NA,40)),year = c(surv.years,2020))
N.per.bin <- data.frame(N.per.bin = rbind(survey.obj$BBn$shf.dat$n.yst,rep(NA,40)),year = c(surv.years,2020))
#reorder them
mw.per.bin <- mw.per.bin[order(mw.per.bin$year),]
N.per.bin <- N.per.bin[order(N.per.bin$year),]
# Get the right bins for the FRs
max.bin <- length(sizes)
bin.frs.plus <- which(sizes == 1.025):max.bin
bin.90.plus <- which(sizes == 0.925):max.bin
bin.rec <- which(sizes == 0.775):min((bin.90.plus-1))
bin.frs.minus <- min(bin.90.plus):(min(bin.90.plus)+1)

# and the right bins for the recruits

# Now make a new object
g.proper <- data.frame(year = mw.per.bin$year)
g.proper$total.abun.90 <- rowSums(N.per.bin[,bin.90.plus])
g.proper$total.abun.frs <- rowSums(N.per.bin[,bin.frs.plus])
g.proper$total.rec.abun <- rowSums(N.per.bin[,bin.rec])
g.proper$total.frs.minus <- rowSums(N.per.bin[,bin.frs.minus])
# Propotions in each bin, FRs and
N.prop.per.bin.90 <- N.per.bin[,bin.90.plus]/g.proper$total.abun.90
N.prop.per.bin.frs <- N.per.bin[,bin.frs.plus]/g.proper$total.abun.frs
# Recs
N.prop.per.bin.rec       <- N.per.bin[,bin.rec]/g.proper$total.rec.abun
N.prop.per.bin.frs.minus <- N.per.bin[,bin.frs.minus]/g.proper$total.frs.minus

# And the average mw in each of the bins of interest, first for the FRs
g.proper$mw.frs.plus <-  rowSums(mw.per.bin[,bin.frs.plus] * N.prop.per.bin.frs,na.rm=T)
g.proper$mw.90.plus <-   rowSums(mw.per.bin[,bin.90.plus] * N.prop.per.bin.90,na.rm=T)
# and for the rec
g.proper$mw.recs <-      rowSums(mw.per.bin[,bin.rec] * N.prop.per.bin.rec,na.rm=T)
g.proper$mw.frs.minus <- rowSums(mw.per.bin[,bin.frs.minus] * N.prop.per.bin.frs.minus,na.rm=T)

g.proper$g.proper <- c(g.proper$mw.frs.plus[2:length(g.proper$mw.frs.plus)]/g.proper$mw.90.plus[1:(length(g.proper$mw.90.plus)-1)],NA)
g.proper$gR.proper<- c(g.proper$mw.frs.minus[2:length(g.proper$mw.frs.minus)]/g.proper$mw.recs[1:(length(g.proper$mw.recs)-1)],NA)


g.proper[g.proper$year %in% c(1991,2020),-1] <- NA
g.proper[g.proper$year %in% c(2019),which(names(g.proper) %in% c("g.proper","gR.proper"))] <- NA

# Fill in the mean for the missing years
g.proper$g.proper[g.proper$year %in% c(1991,2019,2020,2022)] <- median(g.proper$g.proper,na.rm=T)
g.proper$gR.proper[g.proper$year %in% c(1991,2019,2020,2022)] <- median(g.proper$gR.proper,na.rm=T)
# Now use the new growth data.
mod.dat$g <- g.proper$g.proper
mod.dat$gR <- g.proper$gR.proper



strt.mod.yr <- 1994
years <-strt.mod.yr:max(mod.dat$year)
strt.retro.yr <- 2005
retro.years <- strt.retro.yr:max(years)

for(i in retro.years)
{
  yrs <- strt.mod.yr:i
  NY<- length(yrs)
  
  DD.dat <- subset(mod.dat,year %in% yrs,
                   select = c("year","n.x","I","I.cv","IR",  "IR.cv", "IPR", "IPR.cv","N","N.cv","NR","NR.cv", "NPR", "NPR.cv",
                              "w.bar","l.bar", "l.k", "w.k","CF","clappers","clappersR","CS",  "RS","catch","effort","n.y","cpue",
                              "cpue.var","cpue.se","LCI","UCI","U.cv", "g","gR"))
  
  names(DD.dat) <- c( "year","n","I","I.cv","IR",  "IR.cv", "IPR", "IPR.cv","N","N.cv","NR","NR.cv", "NPR", "NPR.cv",
                      "w.bar","l.bar", "l.k", "w.k","CF","clappers","clappersR","CS",  "RS","C","E","n.trips","U",
                      "U.var","U.se","LCI","UCI","U.cv", "g","gR") 
  # Organize the data and set up the model priors/initialization data, then run the model.
  yrs<-min(DD.dat$year):max(DD.dat$year)
  
  DD.lst<-as.list(subset(DD.dat,year %in% yrs,c("I","I.cv","IR","IR.cv","g","gR","C","U","U.cv","N","NR","clappers",
                                                "clappersR")))       
# Add a couple items to the DD.lst list...
DD.lst$NY<- length(DD.lst$C)
DD.lst$year<-min(DD.dat$year):max(DD.dat$year)
# Set up Priors.  This first bit is getting our variance correct for the CV's for Biomass, Recruit biomass, and catch rates.
# This is then added to our list of priors to get them correct.
# Biomass CV
uI=log(DD.lst$I.cv^2+1) # See Smith and Hubley 2014 for details, this is variance of the log of a CV
# DK Note:  Smith/Hubley suggest this should be 3, so why we setting it to 2???
Ip.a=2+(uI/uI)^2 # This is the alpha prior term for the prior (an inverse-gamma, i.e. gamma using 1/var); a rather funky way of setting alpha =2
Ip.b=1/(uI*((uI/uI)^2+1)) # This is the beta term for the prior, again a strangly complex way of saying 1/(2*(uI))
# Recruit biomass CV, see above comments for details.
uIR=log(DD.lst$IR.cv^2+1)
IRp.a=2+(uIR/uIR)^2
IRp.b=1/(uIR*((uIR/uIR)^2+1))
# Catch Rate CV, see above comments for details.
# uU=log(DD.lst$U.cv^2+1)
# Up.a=2+(uU/uU)^2
# Up.b=1/(uU*((uU/uU)^2+1))

DDpriors=list(
  logK=			    list(a=7,		  b=7,		d="dnorm",	l=1		),		# scaler to total biomass, a= mean  b = sd, this gives a huge range of starting values
  r=				    list(a=0, 		b=1,		d="dlnorm",	l=NY	),		# scaled recruit biomass, a= meanlog  b = sdlog
  m=				    list(a=-2,		b=2,		d="dlnorm",	l=NY	),		# natural mortality fully recruited a= meanlog  b = sdlog
  mR=				    list(a=-2,		b=2,		d="dlnorm",	l=NY	),		# natural mortality  recruits a= meanlog  b = sdlog
  S=				    list(a=1.1454e3, 		b=1e4,		d="dbeta",  l=1		),		# clapper dissolution rate a= shape1, b=shape2, 8 & 11 gives ~ normal mean of .45ish
  SR=				    list(a=2.285e3, 		b=1e4,		d="dbeta",  l=1		),		# clapper dissolution rate a= shape1, b=shape2, 8 & 11 gives ~ normal mean of .45ish
  q=				    list(a=20, 		b=40,		d="dbeta",	l=1		),		# survey catchability fully recruited a= shape1, b=shape2
  #qU=				    list(a=0,		  b=1,	  d="dunif",	l=1		),		# fishery catchability CPUE a= min, b = max
  sigma=			  list(a=0, 		b=5,		d="dunif",	l=1		),		# process error (SD) a = min, b = max
  ikappa.tau2=	list(a=3, 		b=2.2407,	d="dgamma",	l=1		),	# measurement error FR clappers  a = shape, b = scale (1/rate)
  ikappa.rho2=	list(a=3, 		b=2.2407,	d="dgamma",	l=1		),	# measurement error recruit clappers a = shape, b = scale (1/rate)
  I.precision=	list(a=Ip.a,	b=Ip.b,	d="dgamma",	l=NY	),		# measurement error variance survey FR a = shape, b = scale (1/rate)
  IR.precision=	list(a=IRp.a,	b=IRp.b,d="dgamma",	l=NY	)		# measurement error variance survey recruits a = shape, b = scale (1/rate)
  #U.precision=	list(a=Up.a,	b=Up.b,	d="dgamma",	l=NY	)		  # measurement error variance CPUE  a = shape, b = scale
)

#Prepare priors for JAGS
for(h in 1:length(DDpriors))
{
  # Get the variances for log-normal and normal converted to precisions, note that in BUGS language the precision is
  # the inverse of the squared standard deviation (which is what you specify in R).  The standard deviation is what
  # was specified in the Prior list (as it is more intuitive)
  if(DDpriors[[h]]$d%in%c("dlnorm","dnorm")) DDpriors[[h]]$b <- 1/DDpriors[[h]]$b^2
  # For a Gamma to convert to precision the precision term is  the inverse of the 'Scale" term in a typical 
  # gamma distribution parameterization, aka this is now knonwn as the rate.
  # Happily this is the same as the parameterization in R dgamma(x,shape,rate) so our b parameter is correct for posterior plots.
  if(DDpriors[[h]]$d=="dgamma")DDpriors[[h]]$b<-1/DDpriors[[h]]$b
} # end for(h in 1:length(DDpriors))
# Made a data.frame of the priors, unwrap the list and combine by row.
prior.dat<- data.frame(par=names(DDpriors),do.call("rbind",lapply(DDpriors,rbind)))
prior.lst<-list()
# Now turn this into a list
for(k in seq(1,nrow(prior.dat)*2,2))
{
  prior.lst[[k]]<-prior.dat$a[[ceiling(k/2)]]
  prior.lst[[k+1]]<-prior.dat$b[[ceiling(k/2)]]
} # end for(k in seq(1,nrow(prior.dat)*2,2))
# And give the list names
names(prior.lst)<-paste(rep(prior.dat$par,2)[order(rep(1:nrow(prior.dat),2))],rep(c('a','b'),nrow(prior.dat)),sep='.')

# Now if they haven't already been selected grab the parameters you want for the model.
parameters <- c(names(DDpriors),'K','P','B','R','mu','Imed','Ipred','Irep', 'IRmed','IRpred','IRrep',
                                                 "Cmed","Crep","CRmed","CRrep",'sIresid','sIRresid','sPresid','Iresid',
                                                 'IRresid','Presid',"Cresid","CRresid","sCresid","sCRresid")
# Run the model and see how long it takes.
# n = 400,000 and burn = 100,000, thin = 20 with 2 chains do not decrease these as retaining this much
# data is needed to stabilize the projections, it does lengthen the run time to 10-20 minutes in serial
# Running in parallel stick with that burn in but we can get away with n=200,000, burn = 100,000, thin = 20, and 6 chains
# they are longer chains than really are needed for the model to converge, but this is really being done just for the projections.
# Run the model now.
start<-Sys.time()
## Call to JAGS, do you want to run in parallel?
jags.model = "D:/Github/Framework/Model/DD_no_cpue_fix_S.bug"
out <- jags.parallel(data =  c(prior.lst,DD.lst), inits = NULL,parameters.to.save = parameters,  
                       model.file = jags.model,n.chains = 8, n.iter = 375000, n.burnin = 300000, 
                       n.thin = 20,jags.seed = 1)
print(Sys.time()-start)

# Rename the output so I retain the results 
DD.out <- list(data=c(prior.lst,DD.lst,yrs), sims.list=out$BUGSoutput$sims.list,median=out$BUGSoutput$median,
                      mean=out$BUGSoutput$mean,summary=out$BUGSoutput$summary,priors = prior.lst,parameters=parameters)

# I will also retain the MCMC object produced in case I want it for something.
mod.out <- out

save(DD.lst, DDpriors,DD.out,DD.dat,mod.out,mod.dat,cpue.dat,yrs,
     file=paste0(repo.loc,"/Results/BBn_SS_model/R_75_FR_90/Retros/BBn_SSmodel_results_",strt.mod.yr,"_",i,".RData"))
}


# Now make the figures and the bias/mr tables


# Set parameters for the run...
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
retro.years <- c(2005:2019,2021,2022)
n.retro.years <- length(retro.years)
base.year <- max(retro.years)
strt.mod.yr <- 1994
R.size <- 75
FR.size <- 90
#scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",exp(l.init.m),"_R0_",exp(l.init.R),"_",num.knots,"_knots")
#tst <- readRDS("D:/framework/SFA_25_26_2024/Model/Results/BBn/R_75_FR_90/Retros/BBn_SEAM_model_output_1994_2010_vary_m_m0_0.4_qR_0.5_4_knots.Rds")

# Load the correct base model...

trends <- NULL
for(j in retro.years)
{
  # load the data objects for the appropriate run...
  load(paste0(repo.loc,"Results/BBn_SS_model/R_75_FR_90/Retros/BBn_SSmodel_results_1994_",j,".RData"))
  
  trends[[as.character(j)]] <-  data.frame(years = strt.mod.yr:(j),
                                           B =     apply(mod.out$BUGSoutput$sims.list$B,2,median),
                                           B.LCI = apply(mod.out$BUGSoutput$sims.list$B,2,function(x){quantile(x,probs=0.025)}),
                                           B.UCI = apply(mod.out$BUGSoutput$sims.list$B,2,function(x){quantile(x,probs=0.975)}),
                                           R =     apply(mod.out$BUGSoutput$sims.list$R,2,median),
                                           R.LCI = apply(mod.out$BUGSoutput$sims.list$R,2,function(x){quantile(x,probs=0.025)}),
                                           R.UCI = apply(mod.out$BUGSoutput$sims.list$R,2,function(x){quantile(x,probs=0.975)}),
                                           m =    apply(mod.out$BUGSoutput$sims.list$m,2,median),
                                           m.LCI = apply(mod.out$BUGSoutput$sims.list$m,2,function(x){quantile(x,probs=0.025)}),
                                           m.UCI = apply(mod.out$BUGSoutput$sims.list$m,2,function(x){quantile(x,probs=0.025)}),
                                           retro.year = j)
  
  if(j == base.year) trends[[as.character(j)]]$retro.year <- "Full model"
}    
retro.base <- do.call("rbind",trends)

cols <- c(rep('#005BBB',4),rep('firebrick2',4),rep('darkgrey',4),rep('#FFD500',4),'black')
points <- c(rep(21:24,4),25) 
b.retro <- ggplot(data=retro.base %>% dplyr::filter(years < 2015),aes(x= years, y = B/1000,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + 
  geom_line(data=retro.base %>% dplyr::filter(years %in% 2016:2019),size=1)+
  geom_point(data=retro.base %>% dplyr::filter(years %in% c(1994:2014,2016:2019,2021)),size=3) + 
  scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols) +
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Fully recruited biomass (tonnes x 1000)")
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Retros/SSModel/BBn_Biomass_retro.png"),b.retro,base_width = 10,base_height =7)


r.retro <- ggplot(data=retro.base %>% dplyr::filter(years < 2015),aes(x= years, y = R/1000,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + 
  geom_line(data=retro.base %>% dplyr::filter(years %in% 2016:2019),size=1) +
  geom_point(data=retro.base %>% dplyr::filter(years %in% c(1994:2014,2016:2019,2021)),size=3) + 
  scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols)  + 
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Recruit biomass (tonnes x 1000)")
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Retros/SSModel/BBn_Rec_retro.png"),r.retro,base_width = 10,base_height =7)



m.retro <- ggplot(data=retro.base %>% dplyr::filter(years < 2015),aes(x= years, y = m,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + 
  geom_line(data=retro.base %>% dplyr::filter(years %in% 2016:2019),size=1)+
  geom_point(data=retro.base %>% dplyr::filter(years %in% c(1994:2014,2016:2019,2021)),size=3) + 
  scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols)  + 
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Natural mortality (instantaneous)")
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Retros/SSModel/BBn_mort_retro.png"),m.retro,base_width = 10,base_height =7)


# Calculate mohn's rho, it is simply the Relative bias of the estimate...
# So now bring in the final model run to get the 'true' value for the calculation

bias <- NULL

# Do a 5 year peel as recommended
#act.retro.years <- c(2010:2014,2016:2019,2021)
for(j in retro.years[-length(retro.years)])
{
  retro.dat <- retro.base %>% dplyr::filter(retro.year==j,years ==j)
  base.dat <- retro.base %>% dplyr::filter(years == j,retro.year == "Full model")
  bias[[as.character(j)]] <- data.frame(B = (retro.dat$B - base.dat$B),
                                        rel.B = (retro.dat$B - base.dat$B) / base.dat$B,
                                        R = (retro.dat$R - base.dat$R),
                                        rel.R = (retro.dat$R - base.dat$R) / base.dat$R,
                                        m = (retro.dat$m - base.dat$m),
                                        rel.m = (retro.dat$m - base.dat$m) / base.dat$m,
                                        mod = "BSSM")
}


# So now we can use this to calculate mohn's rho
bias <- do.call('rbind',bias)
# mohn's rho being fine when it is < 0.2 is found in Hutrtado-Ferro 2015 paper.
# Note we are removing 2020 from this since we don't have any survey data.
mohns.rhos <- data.frame(mr.B = sum(bias$rel.B)/n.retro.years,
                         mr.R = sum(bias$rel.R)/n.retro.years,
                         mr.m = sum(bias$rel.m)/n.retro.years,
                         mr.B5 = sum(bias$rel.B[(nrow(bias)-4):nrow(bias)])/5,
                         mr.R5 = sum(bias$rel.R[(nrow(bias)-4):nrow(bias)])/5,
                         mr.m5 = sum(bias$rel.m[(nrow(bias)-4):nrow(bias)])/5,
                         mod = "BSSM")

saveRDS(mohns.rhos,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_mohns_rose_SSModel.Rds"))
saveRDS(bias,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_bias_SSModel.Rds"))
saveRDS(retro.base,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_retro_SSModel.Rds"))

