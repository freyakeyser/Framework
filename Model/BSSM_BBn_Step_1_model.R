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
L.inf <- 164.4
#to <- 1.337 # So this uses a 1 year offset that we no longer believe in, going to make this 0.337 to align more with what we now do...
to <- -0.2
K <- 0.2


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
# ### overwrite imputation for growth here using whichever method
# # in 2020, the covid-19 pandemic prevented the DFO survey from occurring. An industry-lead survey of limited scope occurred, but is not suitable for inclusion in the 
# # assessment models. As such, we need to fill-in the blank row for 2020 with some data. We'll try out different options for doing that here. 
# # We imputed the values in the survey data earlier, but for the "mixed" imputation method, we'll handle growth separately.
# # Maybe it makes more sense to use the LTM for growth but midpoint for other values. 
#     
# # change 2019 and 2020 values to NA    
# 
# mod.dat$g[which(mod.dat$year %in% c(2020))] <- NA
# mod.dat$g2[which(mod.dat$year %in% c(2020))] <- NA
#   
# mod.dat$gR[which(mod.dat$year %in% 2020)] <- NA
# mod.dat$gR2[which(mod.dat$year %in% 2020)] <- NA
#     
# # replace the NAs with long term medians
# mod.dat$g[which(mod.dat$year %in% c(2020))] <- median(mod.dat$g, na.rm=T)
# mod.dat$g2[which(mod.dat$year %in% c(2020))] <- median(mod.dat$g2, na.rm=T)
#     
# mod.dat$gR[which(mod.dat$year %in% c(2020))] <- median(mod.dat$gR, na.rm=T)
# mod.dat$gR2[which(mod.dat$year %in% c(2020))] <- median(mod.dat$gR2, na.rm=T)

#mod.tmp <- read.csv("D:/Framework/SFA_25_26_2024/Model/Data/BBn_input_data_for_freya.csv")
#mod.tmp$g[1:29] - growth$g[1:29]

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

#mod.dat$g <- mod.dat$g -0.2
#mod.dat$gR <- mod.dat$gR -0.4

strt.mod.yr <- 1994
# Grab the data, start model at either 1986 (note that BBn data starts in 1991 so anything earlier will default to 1991)
DD.dat <- subset(mod.dat,year %in% strt.mod.yr:max(mod.dat$year),
                 select = c("year","n.x","I","I.cv","IR",  "IR.cv", "IPR", "IPR.cv","N","N.cv","NR","NR.cv", "NPR", "NPR.cv",
                            "w.bar","l.bar", "l.k", "w.k","CF","clappers","clappersR","CS",  "RS","catch","effort","n.y","cpue",
                            "cpue.var","cpue.se","LCI","UCI","U.cv", "g","gR"))

names(DD.dat) <- c( "year","n","I","I.cv","IR",  "IR.cv", "IPR", "IPR.cv","N","N.cv","NR","NR.cv", "NPR", "NPR.cv",
                    "w.bar","l.bar", "l.k", "w.k","CF","clappers","clappersR","CS",  "RS","C","E","n.trips","U",
                    "U.var","U.se","LCI","UCI","U.cv", "g","gR") 
# Organize the data and set up the model priors/initialization data, then run the model.
yrs<-min(DD.dat$year):max(DD.dat$year)
NY<- length(yrs)
DD.lst<-as.list(subset(DD.dat,year %in% yrs,c("I","I.cv","IR","IR.cv","g","gR","C","U","U.cv","N","NR","clappers",
                                                            "clappersR")))
# DK NOTE: Downweight the CV for the CPUE data. This is done to be consistent with CV used
# Previously in the model assessments. This has been flagged as an action item to investigate 
# and resolve in the next framework.
ifelse(names(DD.lst)[9] == "U.se", names(DD.lst)[9] <- "U.cv", DD.lst$U.cv <- DD.lst$U.cv*50)
# Also, if doing this we need to change the original data to represent what the model is seeing..
# So if we used the SE let's replace the U.cv data with the U.se data, if we are doing the
# 50x to match what we've done before than we need to change those data as well.
ifelse(names(DD.lst)[9] == "U.se", DD.dat$U.cv <- DD.dat$U.se, DD.dat$U.cv <- DD.dat$U.cv*50)

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


mod.out$BUGSoutput$summary

summary(DD.out$summary[588:616,5])
summary(DD.out$summary[617:645,5])
DD.out$summary[582:583,]
rownames(DD.out$summary)
#source("fn/projections.r")
# The catch since the survey for the most recent year is this, if there was no catch set this to 0.
proj.catch <- max(proj.dat$catch[proj.dat$year == max(DD.dat$year)],0)
# Get the low and upper boundaries for the decision table (this might be a silly way to do this...)

# If we are looking at one of the sub-areas we will go for 1/3 of the mean biomass estimate for the current year...
# The increment size for the decision table.  500 for GBa and 50 for BBn
step <- 10
# The URP and LRP for the bank, for the moment only GBa has been accepted so it's the only one used.
# For more info on GBa reference points: see Y:\Offshore\Assessment\Non-Github archive and documentation\Help and Documentation\GBa Reference Points Literature Review.docx 
refyears <- which(DD.out$dat$year %in% 1994:2009)
LRP <- mean(DD.out$median$B[refyears]) * 0.3 
URP <- mean(DD.out$median$B[refyears]) * 0.8 


# Get the projection scenarios of interest
proj <- seq(0,250,step) + proj.catch
# If we don't have projected catch data yet (i.e. I'm running the model before the logs have data in them..)
# The interim TAC is known for GBa and BBn,
TACi <- 200

# Now do the projections
DD.out<- projections(DD.out,C.p=proj) # C.p = potential catches in decision table

### Generate Decision Table ###
### Note that from the 2015 SSR we have these definitely set at...
#The Lower Reference Point (LRP) is 7,137 t and the Upper Stock Reference (USR) is 13,284 t.
D.tab<-decision(DD.out,"BBn", mu=0.15,refs=c(URP,LRP),post.survey.C=proj.catch, yr=2023)
write.csv(D.tab,"D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/BBn_SSModel_Decision_table.csv",row.names=F) #Write2

# For i = 1 this will just get the first bank, unfortunately if i =2 then this will pull in results for both 
#  (if running this as a loop) which is silly, but work arounds are dumber than this solution
# If you are happy and want to keep these results 


save(DD.lst, DDpriors,DD.out,DD.dat,mod.out,mod.dat,cpue.dat,proj.dat,D.tab,proj.catch,
     URP,LRP,proj,TACi,yrs,
     file="D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/BBn_SSmodel_results.RData")
saveRDS(DD.out,file = "D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/BBn_SS_mod_output.Rds")


# Here we can grab the Fully recruited and recruit biomass for the last 2 years and the median of the time series.
FR.bm <- DD.out$median$B
# We exclude the current year from the median estimate
FR.ltm <- median(DD.out$median$B[-length(DD.out$median$B)])
# Recruit biomass
rec.bm <- DD.out$median$R
# We exclude the current year from the median estimate
rec.ltm <- median(DD.out$median$R[-length(DD.out$median$R)])
# natural mortality
mort <- 1- exp(-DD.out$median$m)

# This lines up the column headers with the projected catch...
TACI<- which(DD.out$data$C.p==(TACi+proj.catch))
# This get us the predicted biomass for next year based on the projected catch
BM.proj.1yr <- DD.out$median$B.p[TACI]
# This is only useful for GBa at the moment since BBn doesn't have reference points accepted yet...
# Get the quantiles, this likely would need changed, but which quantile is > our URP (13,284 as of 2015)
B.quantiles <- quantile(DD.out$sims.list$B[,length(DD.out$sims.list$B[1,])],probs=seq(0,1,0.01))
# This is the probability (well percentage) that Biomass is below the USR
prob.below.USR <- names((which(B.quantiles > URP)[1]))


# Here we can grab the Fully recruited and recruit biomass for the last 2 years and the median of the time series.
FR.bm <- DD.out$median$B[(length(DD.out$mean$B)-1):length(DD.out$median$B)]
# We exclude the current year from the median estimate
FR.ltm <- median(DD.out$median$B[-length(DD.out$median$B)])
# Recruit biomass
rec.bm <- DD.out$median$R[(length(DD.out$median$R)-1):length(DD.out$median$R)]
# We exclude the current year from the median estimate
rec.ltm <- median(DD.out$median$R[-length(DD.out$median$R)])

# Get the percent biomass change from the projection. 0 means unchanged, + means % increase, - means % decline
percent.B.change <- (BM.proj.1yr / DD.out$median$B[length(DD.out$median$B)]) -1

####################  MODEL DIAGNOSITCS ####################  MODEL DIAGNOSITCS ####################  MODEL DIAGNOSITCS 
##### Now we can run some model diagnostics.
# Some quick diagnoistics, the maximum should be < 1.05
rhat <- summary(DD.out$summary[,8])

# Effective number of observations.  
#Not sure what our minimum should be here, but using the Rhat + looking at the chains should indicate where there are problems...
neff <- range(DD.out$summary[,9])

# extract the catchability posterior...
q <- mod.out$BUGSoutput$sims.list$q

save(mort,TACI,BM.proj.1yr,B.quantiles,percent.B.change,prob.below.USR,FR.bm,FR.ltm,rec.bm,rec.ltm,neff,rhat,q,
                           file="D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/Model_results_and_diagnostics.RData")
saveRDS(q,file = "D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/q_posterior.Rds")


# OK, so happy with model lets load up some stuff
load("D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/Model_results_and_diagnostics.RData")
load("D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/BBn_SSmodel_results.RData")


# Let's extract the catchability posterior....

# I'd like to pull out the Biomass, mortality

#####Plot model diagnostics############## 
# These plots include the posterior fits, exploitation estimate, Biomass fit to survey and CPUE, residual plot
# and the model convergence plot (which is a 700+ page pdf of the convergence of each parameter + it's ACF.)
# posterior densities for model parameters
fig <- 'png'
plotsGo <- "D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/"
language = 'en'
post.plt(DD.out,DDpriors,years=yrs, graphic=fig,multi=T,path=plotsGo)
#dev.off()
#exploitaiton time series
exploit.plt(DD.out, years=yrs, plt=c('f','m','mR'),graphic=fig,path=plotsGo)
#dev.off()

# for 2020, we have to insert NAs because of COVID non-survey, and 2015 because of new vessel non-survey.
DD.plt <- DD.out
#DD.plt$median$B[which(yrs %in% c(2015,2020))] <- NA
#DD.plt$sims.list$B[,which(yrs %in% c(2015,2020))] <- NA
#DD.plt$median$R[which(yrs %in% c(2015,2020))] <- NA
#DD.plt$sims.list$R[,which(yrs %in% c(2015,2020))] <- NA
#DD.plt$data$I[which(yrs %in% c(2015,2020))] <- NA
#DD.plt$data$IR[which(yrs %in% c(2015,2020))] <- NA

# model biomass fit to survey
#fit.plt(DD.plt, years = yrs, CI=T,graphic=fig,path=plotsGo,CV=T, language=language)
# diagnostic plot
diag.plt(DD.out, years = yrs,graphic=fig,path=plotsGo)

# Here we pull together all of the chains and look to see that they are both well mixed and that
# there is no correlation.   This is a complete crap load of plots!!
# Function to plot all the chains.
# Get the bank pulled out and figure out how many parameters we have
num.param <- length(names(DD.out$sims.list))
param.names <- names(DD.out$sims.list)
# Make the pdf, given the number of parameters in the model you don't get an option for making this plot print to screen
# if you run diagnostics you get this pdf

# Set the number of rows, this is based off the number of chains (I used 8 chains)
nchains = 8
nr <- ceiling(sqrt(nchains))
# Set the number of columns, the funky little command is used to add one to the nc if nchains is a perfect square
# As we are printing nchains + 1
ifelse(sqrt(nchains)%%1==0,  nc <- ceiling(sqrt(nchains))+1, nc <- ceiling(sqrt(nchains)))
# I always force this to make a pdf because it is a bazillion pages long...
pdf(file=paste(plotsGo,"Model_convergence.pdf",sep=""),onefile=T)
# Set up the plotting device.
par(mfrow = c(nr,nc),mar=c(2,2,3,1))
for(i in 1:num.param)
{
  # This pulls out all the plot for parameters with only one value
  if(is.vector(DD.out$sims.list[[names(DD.out$sims.list)[i]]])==T)
  {
    # Since our all of these parameters are hyperparameter (there isn't one for every year) this works, 
    #if this was a matrix we'd get an error here.
    len <- length(DD.out$sims.list[[param.names[i]]])
    # This sets us up to pull out the right "chains" from the data using the k loop below.
    bins <- seq(1,len,by = len/nchains)
    # Get the ylimits for the plot.
    ylims <- range(DD.out$sims.list[[param.names[i]]])
    colr <- rainbow(nchains) # set up a color ramp
    count <- 0 # Set up a counter
    # I need to run this loop twice, once for the first figure so we get all the lines added and a second time to 
    # add in the ACF's.  Probably a nicer way to do this, but it'll do...
    # Now run the loop and make the figure showing the mixing of the chains
    for(k in bins)
    {
      count <- count+1 # used for the color ramp
      # Get the data
      dat <-DD.out$sims.list[[param.names[i]]][k:(k+len/nchains-1)] 
      if(k ==1) plot(dat,type="l",col=colr[count], main = paste(param.names[i], "Chain"),xlab="",ylab="",ylim=ylims)
      if(k > 1) lines(dat,col=colr[count])
    } # end for(k in bins)
    
    # Now to make the ACF figures
    count  <- 0 # Reset the counter
    for(k in bins)
    {
      count <- count+1 # used to ID the chain
      # Pick it up from here.
      dat <-DD.out$sims.list[[param.names[i]]][k:(k+len/nchains-1)] 
      # And look for any signs of autocorrelation in the chains...
      acf(dat,lag.max = 10,main = paste("ACF chain",count),xlab="",ylab="",ylim=c(0,0.3))
    }# end for(k in bins)
    
  } # end if(is.vector(DD.out$sims.list[[names(DD.out$sims.list)[1]]])==T)
  
  # This pulls out all the plots for parameters with multiple values (i.e. annual estimates)
  if(is.vector(DD.out$sims.list[[names(DD.out$sims.list)[i]]])==F)
  {
    num.reps <- ncol(DD.out$sims.list[[names(DD.out$sims.list)[i]]])
    rep.names <- paste(names(DD.out$sims.list)[i],1:num.reps,sep="_")
    # Run this loop for each chain for these parameters.
    for(p in 1:num.reps)
    {
      # Get the length again (this could probably be tidied up as this number has to be the same for all parameters.)
      len <- length(DD.out$sims.list[[param.names[i]]][,p])
      # Get the bins for the loop
      bins <- seq(1,len,by = len/nchains)
      # Get the ylimits for the plot.
      ylims <- range(DD.out$sims.list[[param.names[i]]][,p])
      colr <- rainbow(nchains) # set up a color ramp
      count <- 0 # Set up a counter
      # Set up the plotting device.
      #par(mfrow = c(nr,nc),mar=c(2,2,3,1))
      # I need to run this loop twice, once for the first figure so we get all the lines added and a second time to 
      # add in the ACF's.  Probably a nicer way to do this, but it'll do...
      # Now run the loop and make the figure showing the mixing of the chains
      for(k in bins)
      {
        count <- count+1 # used for the color ramp
        # Get the data
        dat <-DD.out$sims.list[[param.names[i]]][k:(k+len/nchains-1),p]
        if(k ==1)  plot(dat,type="l",col=colr[count], main = paste(rep.names[p], "Chain"),xlab="",ylab="",ylim=ylims)
        if(k > 1) lines(dat,col=colr[count])
      } # end for(k in bins)
      count <- 0 # Reset the counter
      # And now for the ACF figures
      for(k in bins)
      {
        count <- count+1 # used for the chain index
        # Get the data
        dat <-DD.out$sims.list[[param.names[i]]][k:(k+len/nchains-1),p]
        # And look for any signs of autocorrelation in the chains...
        acf(dat,lag.max = 10,main = paste("ACF chain ",count),xlab="",ylab="",ylim=c(0,0.3))
      } # end for(k in bins)
    }# end for(p in 1:num.reps)
  } # end if(is.vector(DD.out$sims.list[[names(DD.out$sims.list)[i]]])==F)
}  # end for(i in 1:num.param)

# Next up, lets make a decision table.



# get the residuals and Process Error


# Process error
bbn.PE <- as.data.frame(t(apply(data.frame(DD.out$sims.list$Presid),2,function(x){quantile(x,probs=c(0.025,0.5,0.975),na.rm=T)})))
names(bbn.PE) <- c("LCI","median","UCI")
bbn.PE$year <- yrs
bbn.PE$fitted <- NA
bbn.PE$obs <- NA
bbn.PE$type <- "Raw process error (log)"
# standardized
bbn.stan.PE <- as.data.frame(t(apply(data.frame(DD.out$sims.list$sPresid),2,function(x){quantile(x,probs=c(0.025,0.5,0.975),na.rm=T)})))
names(bbn.stan.PE) <- c("LCI","median","UCI")
bbn.stan.PE$year <- yrs
bbn.stan.PE$fitted <- NA
bbn.stan.PE$obs <- NA
bbn.stan.PE$type <- "Standardized process error"
# I residuals...

bbn.I.resids <- as.data.frame(t(apply(data.frame(DD.out$sims.list$Iresid),2,function(x){quantile(x,probs=c(0.025,0.5,0.975),na.rm=T)})))
names(bbn.I.resids) <- c("LCI","median","UCI")
bbn.I.resids$year <- yrs
bbn.I.resids$fitted <- apply(data.frame(DD.out$sims.list$B),2,median)*median(DD.out$sims.list$q)
bbn.I.resids$obs <- DD.out$data$I
bbn.I.resids$type <- "Fully-recruited biomass residual (log)"
# Standardized version
bbn.stan.I.resids <- as.data.frame(t(apply(data.frame(DD.out$sims.list$sIresid),2,function(x){quantile(x,probs=c(0.025,0.5,0.975),na.rm=T)})))
names(bbn.stan.I.resids) <- c("LCI","median","UCI")
bbn.stan.I.resids$year <- yrs
bbn.stan.I.resids$fitted <- apply(data.frame(DD.out$sims.list$B),2,median)*median(DD.out$sims.list$q)
bbn.stan.I.resids$obs <- DD.out$data$I
bbn.stan.I.resids$type <- "Standardized fully-recruited biomass residual"
# IR resids
bbn.IR.resids <- as.data.frame(t(apply(data.frame(DD.out$sims.list$IRresid),2,function(x){quantile(x,probs=c(0.025,0.5,0.975),na.rm=T)})))
names(bbn.IR.resids) <- c("LCI","median","UCI")
bbn.IR.resids$year <- yrs
bbn.IR.resids$fitted <- apply(data.frame(DD.out$sims.list$R),2,median)*median(DD.out$sims.list$q)
bbn.IR.resids$obs <- DD.out$data$IR
bbn.IR.resids$type <- "Recruit biomass residual (log)"
# Standardized version
bbn.stan.IR.resids <- as.data.frame(t(apply(data.frame(DD.out$sims.list$sIRresid),2,function(x){quantile(x,probs=c(0.025,0.5,0.975),na.rm=T)})))
names(bbn.stan.IR.resids) <- c("LCI","median","UCI")
bbn.stan.IR.resids$year <- yrs
bbn.stan.IR.resids$fitted <- apply(data.frame(DD.out$sims.list$R),2,median)*median(DD.out$sims.list$q)
bbn.stan.IR.resids$obs <- DD.out$data$IR
bbn.stan.IR.resids$type <- "Standardized recruited biomass residual"


bbn.pe.resids <- rbind(bbn.PE,bbn.stan.PE,
                       bbn.I.resids,bbn.stan.I.resids,
                       bbn.IR.resids,bbn.stan.IR.resids)

saveRDS(bbn.pe.resids,"D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/PE_and_resids.Rds")


#ggplot(bbn.pe.resids) + geom_point(aes(x=year,y=median)) + facet_wrap(~type,scales='free_y')

median(bbn.stan.PE$median)




# Finally we can put all the key data together into one object to make some nice ts plots

res.gg <- data.frame(med = c(apply(DD.plt$sims.list$B, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T)),
                             apply(DD.plt$sims.list$R, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T)),
                             apply(DD.plt$sims.list$m, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T)),
                             apply(DD.plt$sims.list$mR, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T)),
                             apply(DD.plt$sims.list$mu, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T))),
                     uci = c(apply(DD.plt$sims.list$B, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T)),
                             apply(DD.plt$sims.list$R, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T)),
                             apply(DD.plt$sims.list$m, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T)),
                             apply(DD.plt$sims.list$mR, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T)),
                             apply(DD.plt$sims.list$mu, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T))),
                     lci = c(apply(DD.plt$sims.list$B, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T)),
                             apply(DD.plt$sims.list$R, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T)),
                             apply(DD.plt$sims.list$m, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T)),
                             apply(DD.plt$sims.list$mR, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T)),
                             apply(DD.plt$sims.list$mu, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T))),
                     term = factor(c(rep("Biomass (tonnes)",length(yrs)),
                              rep("Recruit Biomass (tonnes)",length(yrs)),
                              rep("FR Natural Mortality (90+ mm, instantaneous)",length(yrs)),
                              rep("Recruit Natural Mortality (75-90 mm, instantaneous)",length(yrs)),
                              rep("Fishing Mortality (instantaneous)",length(yrs))),
                              levels = c("Biomass (tonnes)","Recruit Biomass (tonnes)", "Fishing Mortality (instantaneous)",
                                        "FR Natural Mortality (90+ mm, instantaneous)","Recruit Natural Mortality (75-90 mm, instantaneous)")),
                     year = c(rep(yrs,5)))

# This is the super useful object to compare with other models...
saveRDS(res.gg,file = "D:/Framework/SFA_25_26_2024/Model/Results/BBn_SS_model/R_75_FR_90/B_R_M_F_summarized.Rds")

#windows(11,11)
p.ssmod.res <- ggplot(res.gg, aes(x=year,y=med)) + geom_line(linewidth=1.5) + facet_wrap(~term,scales = 'free_y') + 
                                                   geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + xlab("")
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/Key_results.png",p.ssmod.res,base_height = 7,base_width = 20)

# Just the biomass

p.ssmod.bm <- ggplot(res.gg %>% dplyr::filter(term == "Biomass (tonnes)"), aes(x=year,y=med)) + geom_line(linewidth=1.5,color='firebrick2') +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + ylab("Fully Recruited Biomass (tonnes)") + xlab("")  + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,2.5e4))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/Biomass.png",p.ssmod.bm,base_height = 8.5,base_width = 11)

# Just the recruits

p.ssmod.rec <- ggplot(res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)"), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + ylab("Recruit Biomass (tonnes)") + xlab("") + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,7e3))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/Recruits.png",p.ssmod.rec,base_height = 8.5,base_width = 11)

# Natural mortality

p.ssmod.nat.mort <- ggplot(res.gg %>% dplyr::filter(term == "FR Natural Mortality (90+ mm, instantaneous)"), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + ylab("Natural Mortality (90+ mm, instantaneous)") + xlab("") + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.9))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/FR_nm.png",p.ssmod.nat.mort,base_height = 8.5,base_width = 11)

p.ssmod.rec.nat.mort <- ggplot(res.gg %>% dplyr::filter(term == "Recruit Natural Mortality (75-90 mm, instantaneous)"), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + ylab("Natural Mortality (70-90 mm, instantaneous)") + xlab("") + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,3.3))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/Rec_nm.png",p.ssmod.rec.nat.mort,base_height = 8.5,base_width = 11)

# Exploitation Rate

p.ssmod.F <- ggplot(res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)"), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + ylab("Fishing Mortality (instantaneous)") + xlab("") + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.5))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/FR_F.png",p.ssmod.F,base_height = 8.5,base_width = 11)

# Remove the years with missing survey


# Just the biomass

p.ssmod.bm <- ggplot(res.gg %>% dplyr::filter(term == "Biomass (tonnes)" & year < 2020), aes(x=year,y=med)) + 
  geom_line(linewidth=1.5,color='firebrick2') +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + 
 # geom_line(data = res.gg %>% dplyr::filter(term == "Biomass (tonnes)" & year %in% 2016:2019),linewidth=1.5,color='firebrick2') +
 # geom_ribbon(data = res.gg %>% dplyr::filter(term == "Biomass (tonnes)" & year %in% 2016:2019),aes(x=year,ymin=lci,ymax=uci),fill='blue',color='blue',alpha = 0.2) + 
  geom_line(data = res.gg %>% dplyr::filter(term == "Biomass (tonnes)" & year %in% 2021:2022), linewidth=1.5,color='firebrick2') +
  geom_ribbon(data = res.gg %>% dplyr::filter(term == "Biomass (tonnes)" & year %in% 2021:2022), aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + 
  ylab("Fully Recruited Biomass (tonnes)") + xlab("")  + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,2.5e4))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/Biomass_no_missing_surveys.png",p.ssmod.bm,base_height = 8.5,base_width = 11)

# Just the recruits

p.ssmod.rec <- ggplot(res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)" & year < 2020), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + 
  #geom_line(data = res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)" & year %in% 2016:2019),linewidth=1.5,color='firebrick2') +
  #geom_ribbon(data = res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)" & year %in% 2016:2019),aes(x=year,ymin=lci,ymax=uci),fill='blue',color='blue',alpha = 0.2) + 
  geom_line(data = res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)" & year %in% 2021:2022), linewidth=1.5,color='firebrick2') +
  geom_ribbon(data = res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)" & year %in% 2021:2022), aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + 
  ylab("Recruit Biomass (tonnes)") + xlab("") + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,7e3))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/Recruits_no_missing_surveys.png",p.ssmod.rec,base_height = 8.5,base_width = 11)

# Natural mortality

p.ssmod.nat.mort <- ggplot(res.gg %>% dplyr::filter(term == "FR Natural Mortality (90+ mm, instantaneous)" & year < 2020), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + 
  #geom_line(data = res.gg %>% dplyr::filter(term == "FR Natural Mortality (90+ mm, instantaneous)" & year %in% 2016:2019),linewidth=1.5,color='firebrick2') +
  #geom_ribbon(data = res.gg %>% dplyr::filter(term == "FR Natural Mortality (90+ mm, instantaneous)" & year %in% 2016:2019),aes(x=year,ymin=lci,ymax=uci),fill='blue',color='blue',alpha = 0.2) + 
  geom_line(data = res.gg %>% dplyr::filter(term == "FR Natural Mortality (90+ mm, instantaneous)" & year %in% 2021:2022), linewidth=1.5,color='firebrick2') +
  geom_ribbon(data = res.gg %>% dplyr::filter(term == "FR Natural Mortality (90+ mm, instantaneous)" & year %in% 2021:2022), aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + 
  
  ylab("Natural Mortality (90+ mm, instantaneous)") + xlab("") + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.9))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/FR_nm_no_missing_surveys.png",p.ssmod.nat.mort,base_height = 8.5,base_width = 11)

# Exploitation Rate

p.ssmod.F <- ggplot(res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)" & year < 2020), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + 
  #geom_line(data = res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)" & year %in% 2016:2019),linewidth=1.5,color='firebrick2') +
  #geom_ribbon(data = res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)" & year %in% 2016:2019),aes(x=year,ymin=lci,ymax=uci),fill='blue',color='blue',alpha = 0.2) + 
  geom_line(data = res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)" & year %in% 2021:2022), linewidth=1.5,color='firebrick2') +
  geom_ribbon(data = res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)" & year %in% 2021:2022), aes(x=year,ymin=lci,ymax=uci),fill='darkblue',color='darkblue',alpha = 0.2) + 
  ylab("Fishing Mortality (instantaneous)") + xlab("") + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.5))
save_plot("D:/Framework/SFA_25_26_2024/Model/Figures/BBn/R_75_FR_90/SSModel/FR_F_no_missing_surveys.png",p.ssmod.F,base_height = 8.5,base_width = 11)
