# So here we'll try and get all the data in the correct structure for the spatial model.
# Process flow for this framework...
#Step 1: Run the SEBDAM ("Bnk"_SEBDAM.R) model, control the modelling options (initial natural mortality (init.m), recruit catchability (qR), number of knots (num.knots))
#        and the growth model being used.  I have 3 of these currently, g_original is the method currently being used.  alt_g is a simple method in which we
#        look at observed mw in the following year, exclude the effect of incoming recruits, and then figure out how much the FR scallop from last year must have grown.
#        proper_g is the final methods, it is more in depth than the other methods as it looks at the mw of scallop > 105 mm in the current year and compares that
#        to the mw of scallop 90+ in the previous year, given our von B curve the 90+ scallop in the previous year should all be 105mm+ in the following year
#        so by doing this we have an estimate of the realized growth of all the 90+ scallop from last year.
#   I am using data from 1994 to 2022 and we have decided on recruits being 75-90 mm in size and FR being 90+ mm

# Step2: Run the retrospective analyses ('Bnk'_Retro.R) for one (or more) of these SEBDAM models. The retros are looking good for this model
# Step3: Run the Prediction Evaluation Simulations ("Bnk_prediction_evaluation_simulations.R) to see how well your model predicts the following year, this looks at 8 different 
#        ways of projecting the biomass for the following year.
# Then we move to the Reference Points work.  Here we just have 2 more scripts to run
# Step4: Run the Reference Points and HCR calculations ("bnk"_Reference_Points.Rmd). This code is fairly complex.  The first half is easy enough, you just make a bunch of 
#        interesting plots.  These plots will then be used to inform the settings for your reference points simulations, this part is complicated, but effectively is 
#        simply looking at the model results and making decisions on how any density dependence or correlations could be used in the simulations.  For simplicity at the moment
#        we are mostly ignoring the correlations and the code using the correlations should get a good once over to make sure it is working properly.
# Step5: Run the Decision Tables ('bnk'_Decision_Table.R), not much point running this until you have your reference points sorted out, once you do Bob is your uncle!

library(SEBDAM)
library(tidyverse)
library(sf)
library(stringr)
library(optimx)
library(parallel)
library(INLA)
library(ggthemes)
library(cowplot)


# Download the function to go from inla to sf
funs <- c("https://raw.githubusercontent.com/Mar-Scal/Assessment_fns/master/Fishery/logs_and_fishery_data.r",
          "https://raw.githubusercontent.com/Mar-Scal/Assessment_fns/master/Maps/pectinid_projector_sf.R",
          "https://raw.githubusercontent.com/Mar-Scal/Assessment_fns/master/Maps/convert_inla_mesh_to_sf.R"
          )
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}




# Get the Sab area outline...

Sab.shape <- st_read("D:/Github/GIS_layers/survey_boundaries/Sab.shp", quiet=T)
# Bring in the survey data from the NAS
# load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/Survey_all_results.Rdata")
# surv.dat <- surv.dat$Sab
# saveRDS(surv.dat,'D:/Github/BBn_model/Results/Sab_surv.dat.RDS')
# Bring in the survey data
load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90RSCS_newMWSH_GBb.RData")
#load("D:/Framework/SFA_25_26_2024/Model/Data/testing_results_framework3.Rdata")# Need to get condition factor out of here too
#mod.dat <- survey.obj$Sab$model.dat
#saveRDS(mod.dat,'D:/Github/BBn_model/Results/Sab_model.dat.RDS')
surv.dat <- surv.dat$Sab
mod.dat <- survey.obj$Sab$model.dat
saveRDS(mod.dat,'D:/Framework/SFA_25_26_2024/Model/Data/Sab_BSSM_model.dat.RDS')

#load("F:/NAS/Offshore/Assessment/Data/Model/2022/Sab/Model_input_midpoint.RData")
# Bring in the fishery data
# logs_and_fish(loc="offshore",year = 1986:2022,direct="Y:/Offshore/Assessment/", get.marfis=F)
#  fish.dat<-merge(new.log.dat,old.log.dat,all=T)
#  fish.dat$ID<-1:nrow(fish.dat)
# 
# # Now we can clip both of these to subset it to the data that I think we need for the analysis....
# # First the fishery data
#  fish.dat <- fish.dat[!is.na(fish.dat$lon),]
#  fish.dat <- fish.dat[!is.na(fish.dat$lat),]
#  fish.dat <- fish.dat[!fish.dat$lon==0,]
#  fish.dat <- fish.dat[!fish.dat$lat==0,]
#  
#  #Sab.fish1 <- fish.dat %>% dplyr::filter(bank == "Sab")
#  Sab.fish <- fish.dat[fish.dat$bank == "Sab",]
#  #Sab.fish2 <- Sab.fish2[!is.na(Sab.fish2$lon),]
# 
# # # Now I want to put a 'survey year' on these because that's what we're gonna need for our modelling... start by porting over the year
#  Sab.fish$survey.year <- Sab.fish$year
# # # DK NOTE: Now this is going to get confusing for us and we may want to tweak SEBDAM for this, but that's a down the road job, not a playing around with model job
# # # But based on the indexing in SEBDAM, I am going to change how we index the survey year data from what we have done with offshore traditionally.
# # # Historically anything from the last half of the year goes into the following years, eg. survey.year 2002 = June 2001- May 2002.
# # # But in SEBDAM we have (B(t-1) - C(t-1)), so let's say we have year 2000 survey biomass, this says we remove the 2000 catch from that
# # # we want that catch to be the catch from June 2000 to May 2001, i.e. we remove the catch before we allow the population to grow
# # # This is what we do in our current model, but we have a different index (C(t) on our model.
# # # Basically survey year 2002 = June 2002 - May 2003 now
# # #DK note: We probably should think more about the survey year fun and how exactly we want to handle removal of catch in our models.
# Sab.fish$month <- lubridate::month(Sab.fish$date)
# Sab.fish$survey.year[Sab.fish$month %in% 1:5] <- Sab.fish$survey.year[Sab.fish$month %in% 1:5] -1
# # Add a fake 2022 data point as there were no removals in 2022
# Sab.fish[nrow(Sab.fish)+1,] <- NA
# Sab.fish$pro.repwt[nrow(Sab.fish)] <- 0
# Sab.fish$year[nrow(Sab.fish)] <- 2022
# Sab.fish$survey.year[nrow(Sab.fish)] <-2022
# Sab.fish$lon[nrow(Sab.fish)] <- -61.68767
# Sab.fish$lat[nrow(Sab.fish)] <- 43.63017
#saveRDS(Sab.fish,'D:/Framework/SFA_25_26_2024/Model/Data/Sab_fish.dat.RDS')
Sab.fish <- readRDS('D:/Framework/SFA_25_26_2024/Model/Data/Sab_fish.dat.RDS')
Sab.fish$pro.repwt <- Sab.fish$pro.repwt/1000 # It looks like what I saved is already in tonnes.


# Set up some stuff
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
mod.select <- "SEAM"
atow<-800*2.4384/10^6 # area of standard tow in km2
num.knots <- 10# 4, 8 and 12?
years <- 1994:2022
NY <- length(years)
R.size <- "75"
FR.size <- "90"
c_sys <- 32620
qR <- 0.33# log recruit catchability making this the same are FR catchability.
init.m <- 0.2 # This is for SEAM, sets first year natural mortality, going to test 2, 0.8, and 0.4
# Various explorations of the g models.
g.mod <- "g_original"
#g.mod <- 'g_1'
#g.mod <- 'alt_g'
#g.mod <- 'proper_g'
# The survey biomass index for 1995 says there were 243 tonnes of recruits that year.
#l.init.R <- log(100) # Going to test 100, 250, and 500.

# Transform Sable to 32620
Sab.shape <- Sab.shape %>% st_transform(crs = c_sys) # Sab is totally in 32620 border so think they are basically equivalent options here
# Just going to use the core area to see if that helps model and the prediction function...
#Sab.tst <- st_cast(Sab.shape, "POLYGON")
#Sab.shape <-Sab.tst[1,]
# OK, so step 1 here is getting the model input that Raphael needs for the model
# The survey data....
live.subset <- surv.dat %>% dplyr::filter(state == 'live')
dead.subset <- surv.dat %>% dplyr::filter(state== "dead")
live.input <- data.frame(I = live.subset$com.bm, IR = live.subset$rec.bm,year = live.subset$year,tow = live.subset$tow,tot.live.com = live.subset$com,lat = live.subset$lat,lon=live.subset$lon)
# Here's an option to use different size classes to see what happens...
# Here is one where we expand the recruit sizes to be 70-90 mm
#live.input <- data.frame(I = live.subset$com.bm, IR = (live.subset$`bin_70-80_bm` +live.subset$`bin_80-90_bm`) ,year = live.subset$year,tow = live.subset$tow,tot.live.com = live.subset$com,lat = live.subset$lat,lon=live.subset$lon)
# Here is one where 50-80 is recruits and 80+ is Fully recruited
# live.input <- data.frame(I = (live.subset$`bin_90-120_bm` + live.subset$`bin_120_plus_bm` + live.subset$`bin_80-90_bm`), 
#                          IR = (live.subset$`bin_50-70_bm` + live.subset$`bin_70-80_bm`),
#                          year = live.subset$year,tow = live.subset$tow,tot.live.com = live.subset$com,lat = live.subset$lat,lon=live.subset$lon)

clap.input <- data.frame(L = dead.subset$com,tow = dead.subset$tow,year = dead.subset$year)
mod.input <- left_join(live.input,clap.input,by=c('tow','year'))
mod.input$N <- round(mod.input$tot.live.com + mod.input$L)
# Looks like there are no values > 0 but < 0.5, so the low clapper numbers should all round up to 1 (which makes sense as you'd only get < 0.5 if we had tows twice as long as they should be)
mod.input$L <- round(mod.input$L) 
mod.input.sf <- st_as_sf(mod.input,coords = c('lon','lat'),remove=F,crs = 4326)
mod.input.sf <- mod.input.sf %>% st_transform(crs=c_sys)
mod.input.sf <- mod.input.sf %>% dplyr::filter(year %in% years)
mod.input.sf <- st_intersection(mod.input.sf,Sab.shape)
#mod.input.sf[nrow(mod.input.sf)+1,] <- mod.input.sf[nrow(mod.input.sf),]
#mod.input.sf$year[nrow(mod.input.sf)] <- 2015

# Now I need to get the I and IR into kg/km^2
mod.input.sf$Year <- mod.input.sf$year - (min(years)-1)
mod.input.sf$I <- mod.input.sf$I/atow
mod.input.sf$IR <- mod.input.sf$IR/atow

#survey.obj$Sab$model.dat

# Adding missing survey years.  THe L and N both need 'data', but 0's are fine, so that'll do the trick!
mod.input.sf[nrow(mod.input.sf)+1,] <- mod.input.sf[nrow(mod.input.sf),]
#mod.input.sf[nrow(mod.input.sf),] <- mod.input.sf[nrow(mod.input.sf)-1,]
mod.input.sf$year[nrow(mod.input.sf)] <- 2015
mod.input.sf$Year[nrow(mod.input.sf)] <- which(years == 2015)
mod.input.sf$I[nrow(mod.input.sf)] <- NA
mod.input.sf$IR[nrow(mod.input.sf)] <- NA
mod.input.sf$tot.live.com[nrow(mod.input.sf)] <- NA
mod.input.sf$L[nrow(mod.input.sf)] <- 0
mod.input.sf$N[nrow(mod.input.sf)] <- 0
# And repeat for 2020
mod.input.sf[nrow(mod.input.sf)+1,] <- mod.input.sf[nrow(mod.input.sf),]
#mod.input.sf[nrow(mod.input.sf),] <- mod.input.sf[nrow(mod.input.sf)-1,]
mod.input.sf$year[nrow(mod.input.sf)] <- 2020
mod.input.sf$Year[nrow(mod.input.sf)] <- which(years == 2020)
mod.input.sf$I[nrow(mod.input.sf)] <- NA
mod.input.sf$IR[nrow(mod.input.sf)] <- NA
mod.input.sf$tot.live.com[nrow(mod.input.sf)] <- NA
mod.input.sf$L[nrow(mod.input.sf)] <- 0
mod.input.sf$N[nrow(mod.input.sf)] <- 0

# Now clip mod.input.sf to the right nubmer of years...
mod.input.sf <- mod.input.sf %>% dplyr::filter(year %in% years)


# Growth!! 
# These parameters come from our von B data that was aged by Trish in 2023.
vonB.par <- data.frame(Linf = 159.2,K = 0.2, t0 = 0.2)

w.fr.current  <- mod.dat$CF*(mod.dat$l.bar/100)^3
w.rec.current <- mod.dat$CF*(mod.dat$l.k/100)^3

# Using this years average shell height we can figure out how old the scallop are on average and then we use the 
# von B and allow them to grow by 1 year, that's our average projected size next year.
len.fr.next <- NA #laa.t <- NA
len.rec.next <- NA #laa.t <- NA
for(y in 1:nrow(mod.dat))
{
  age= data.frame(age = seq(2,8,by=0.01),len = NA,l.fr = mod.dat$l.bar[y],l.rec = mod.dat$l.k[y])
  age$len <- vonB.par$Linf*(1-exp(-vonB.par$K*(age$age-vonB.par$t0)))
  age$fr.diff <- abs(age$len - age$l.fr)
  age$rec.diff <- abs(age$len - age$l.rec)
  age.fr.next <- age$age[which.min(age$fr.diff)] + 1
  age.rec.next <- age$age[which.min(age$rec.diff)] + 1
  len.fr.next[y] <-  vonB.par$Linf*(1-exp(-vonB.par$K*(age.fr.next-vonB.par$t0)))
  len.rec.next[y] <-  vonB.par$Linf*(1-exp(-vonB.par$K*(age.rec.next-vonB.par$t0)))
} # end for y in 1:NY
# The c() term in the below offsets the condition so that current year's condition slots into the previous year and repeats 
# the condition for the final year), this effectively lines up "next year's condition" with "predictied shell height next year (laa.t)
# This gets us the predicted weight of the current crop of scallops next year based on next years CF * length^3
# Get the weight of scallop next year
w.fr.next <- c(mod.dat$CF[-1],mod.dat$CF[nrow(mod.dat)])*(len.fr.next/100)^3
w.rec.next <- c(mod.dat$CF[-1],mod.dat$CF[nrow(mod.dat)])*(len.rec.next/100)^3
# Also calculate it based on 'known' condition, used for prediction evaluation figures only as we don't know it when we run the models, only
# after we get survey data.
w.fr.next.alt <- mod.dat$CF*(len.fr.next/100)^3
w.rec.next.alt <- mod.dat$CF*(len.rec.next/100)^3
# Now the growth, expected and realized.

mod.dat$g <- w.fr.next/w.fr.current
mod.dat$gR <- w.rec.next/w.rec.current
# This is using the actual condition factor and g
mod.dat$g2 <- w.fr.next.alt/w.fr.current
mod.dat$gR2 <- w.rec.next.alt/w.rec.current

# The DK growth model #1 is simply using the estimated mw of the FR scallop as the growth term...
dk.growth <- NA
for(i in 2:nrow(mod.dat)) dk.growth[i-1] <- 1+ ((mod.dat$w.bar[i] - mod.dat$w.bar[i-1])/mod.dat$w.bar[i-1])
mod.dat$g.new <- c(dk.growth,dk.growth[length(dk.growth)])

# DK growth model #2 is attempting to account for the recruit contribution to the MW estimates
head(mod.dat)
# With the recruit and FR q assumed to be approx 0.45 and natural mortality shared, there is no need to tweak the recruit and FR ratios
# at all. If we did we could multiply each by 0.9 or so and divide by whatever q is, 0.45 is the same for both currently.
alt.g <- data.frame(year = c(mod.dat$year,2015,2020),B.rec = c(unlist(mod.dat$IR),NA,NA), B.fr = c(unlist(mod.dat$I),NA,NA),wgt.fr = c(mod.dat$w.bar,NA,NA), 
                    sh.fr = c(mod.dat$l.bar,NA,NA),wgt.rec = c(mod.dat$w.k,NA,NA),sh.rec = c(mod.dat$l.k,NA,NA),
                    sh.rec.nxt = c(len.rec.next,NA,NA),wgt.rec.nxt = c(w.rec.next,NA,NA))
alt.g <- alt.g[order(alt.g$year),]

# So if we can figure out what the expected weight of the recruits will be next year using the von B we
# can then figure out how to remove that from the average weight, and I have that size above in the 
# len.rec.next vector
# So we can get the "g" for the year, which will be wgt.rec.nxt/wgt.rec
alt.g$g.rec <- 1 # I'm not going to allow the recruits to grow, just using the biomass from last year, it gets weird with that big gR term...
alt.g$B.rec.2.fr <- alt.g$B.rec * (alt.g$g.rec)
alt.g$B.rec.2.fr <- c(NA,alt.g$B.rec.2.fr[-nrow(alt.g)])
alt.g$B.FR.excluding.new.rec <- alt.g$B.fr - alt.g$B.rec.2.fr
alt.g$wgt.of.last.yrs.recs <- c(NA,alt.g$wgt.rec.nxt[-nrow(alt.g)])
alt.g$prop.rec <- alt.g$B.rec.2.fr/alt.g$B.fr
alt.g$prop.FR <- alt.g$B.FR.excluding.new.rec/alt.g$B.fr


# So using some rearranged weighted averages I should be able to get a new wgt.fr column that excludes the recruits and gets us a 'real' average size.
# Given sh Rec is between 97 and 99 we could also just calculate the change in weight of everything above 100 mm using the raw sh frequency data
# to give a possibly better version of this.
# So we can do some weighted averaging here to figure out what the wgt of the 'old' FR scallop were

alt.g$wgt.fr.excluding.new.rec <- (alt.g$wgt.fr - alt.g$prop.rec*alt.g$wgt.rec.nxt) / alt.g$prop.FR
#alt.g$wgt.fr.excluding.new.rec[find.negs] <- mean(alt.g$wgt.fr.excluding.new.rec[c(find.negs-1,find.negs+1)])
alt.g$alt.g <- c(alt.g$wgt.fr.excluding.new.rec[2:nrow(alt.g)]/alt.g$wgt.fr.excluding.new.rec[1:(nrow(alt.g)-1)],NA)
# Now make the NAs the mean
alt.g$alt.g[which(is.na(alt.g$alt.g))] <- mean(alt.g$alt.g,na.rm=T)




# The other option is to calculate the growth using w.bar of everything over 105 mm, which will be default exclude the vast majority of the
# recruits as 90 mm scallop will grow by about 17 cm, so might have a few recruits in there, but tracking the changes in that size class tells
# us what the realized growth was for the FRs that excludes the recruits
# So what we do is take the ratio of the w.bar for everything bigger than 105 mm in year 2, to the w.bar for all FR scallop in year one
# Based on the von.B the vast majority of the scallop in that ratio be the same individuals.
# So to calculate the 105 mm thing I'll need to use the shf in surv.dat...

sizes <- seq(0.025,2,by=0.05) # So I'd be using the 1.075 bin and everything bigger
# The w.yst object is exactly proportional to mod.dat$I, there is an offset, but given I need proportions I think this object is perfectly fine to use.
mw.per.bin <- data.frame(mw.per.bin = rbind(survey.obj$Sab$shf.dat$w.yst/survey.obj$Sab$shf.dat$n.yst,rep(NA,40),rep(NA,40)),year = c(mod.dat$year,2015,2020))
B.per.bin <- data.frame(B.per.bin = rbind(survey.obj$Sab$shf.dat$w.yst,rep(NA,40),rep(NA,40)),year = c(mod.dat$year,2015,2020))
#reorder them
mw.per.bin <- mw.per.bin[order(mw.per.bin$year),]
B.per.bin <- B.per.bin[order(B.per.bin$year),]
# Get the right bins for the FRs
max.bin <- length(sizes)
bin.105.plus <- which(sizes == 1.075):max.bin
bin.90.plus <- which(sizes == 0.925):max.bin
bin.rec <- which(sizes == 0.775):min((bin.90.plus-1))
bin.105.minus <- min(bin.90.plus):which(sizes == 1.025)

# and the right bins for the recruits

# Now make a new object
g.proper <- data.frame(year = mw.per.bin$year)

g.proper$total.biomass.90 <- rowSums(B.per.bin[,bin.90.plus])
g.proper$total.biomass.105 <- rowSums(B.per.bin[,bin.105.plus])
g.proper$total.rec.biomass <- rowSums(B.per.bin[,bin.rec])
g.proper$total.105.minus <- rowSums(B.per.bin[,bin.105.minus])
# Propotions in each bin, FRs and
B.prop.per.bin.90 <- B.per.bin[,bin.90.plus]/g.proper$total.biomass.90
B.prop.per.bin.105 <- B.per.bin[,bin.105.plus]/g.proper$total.biomass.105
# Recs
B.prop.per.bin.rec       <- B.per.bin[,bin.rec]/g.proper$total.rec.biomass
B.prop.per.bin.105.minus <- B.per.bin[,bin.105.minus]/g.proper$total.105.minus

# And the average mw in each of the bins of interest, first for the FRs
g.proper$mw.105.plus <-  rowSums(mw.per.bin[,bin.105.plus] * B.prop.per.bin.105,na.rm=T)
g.proper$mw.90.plus <-   rowSums(mw.per.bin[,bin.90.plus] * B.prop.per.bin.90,na.rm=T)
# and for the rec
g.proper$mw.recs <-      rowSums(mw.per.bin[,bin.rec] * B.prop.per.bin.rec,na.rm=T)
g.proper$mw.105.minus <- rowSums(mw.per.bin[,bin.105.minus] * B.prop.per.bin.105.minus,na.rm=T)

# Now calculate our g
g.proper$g.proper <- c(g.proper$mw.105.plus[2:length(g.proper$mw.105.plus)]/g.proper$mw.90.plus[1:(length(g.proper$mw.90.plus)-1)],NA)
g.proper$gR.proper<- c(g.proper$mw.105.minus[2:length(g.proper$mw.105.minus)]/g.proper$mw.recs[1:(length(g.proper$mw.recs)-1)],NA)

# And tidy up and put in mean values for the missing years (means are lower than medians mostly for these data)
g.proper[g.proper$year %in% c(1986:1991,2015,2020),-1] <- NA
g.proper[g.proper$year %in% c(2014,2019),which(names(g.proper) %in% c("g.proper","gR.proper"))] <- NA

# Fill in the mean for the missing years
g.proper$g.proper[g.proper$year %in% c(1986:1991,2014:2015,2019:2020,2022)]  <- median(g.proper$g.proper,na.rm=T)
g.proper$gR.proper[g.proper$year %in% c(1986:1991,2014:2015,2019:2020,2022)] <- median(g.proper$gR.proper,na.rm=T)


# now need to add in 2015 and 2020 to mod.dat...
mod.dat.tmp <- mod.dat
mod.dat.tmp[nrow(mod.dat.tmp)+1:2,] <- NA
mod.dat.tmp$year[nrow(mod.dat.tmp):(nrow(mod.dat.tmp)-1)] <- c(2015,2020)
mod.dat.tmp <- mod.dat.tmp[order(mod.dat.tmp$year),]

growth <- data.frame(year = mod.dat.tmp$year,g = mod.dat.tmp$g, gR = mod.dat.tmp$gR,
                     g.alt = alt.g$alt.g, gR.alt = mod.dat.tmp$gR, # I don't have a good idea how to estiamte gR growth, so using the other way
                     g.proper = g.proper$g.proper,gR.proper = g.proper$gR.proper)
# Now addin the missing growth years for g and gR
growth$g[growth$year %in% c(2015,2020)] <- median(growth$g,na.rm=T)
growth[growth$year %in% c(2015,2020),names(growth) %in% c("gR","gR.alt")] <- median(growth$gR,na.rm=T)
growth <- growth[which(!is.na(growth$g)),]
growth[nrow(growth)+1,] <- growth[nrow(growth),]
growth$year[nrow(growth)] <- max(years) + 1

growth <- growth %>% dplyr::filter(year >= min(years))

if(g.mod == 'g_original') g <- data.frame(g=growth$g,gR = growth$gR)
if(g.mod == 'alt_g') g <- data.frame(g=growth$g.alt,gR = growth$gR.alt)
if(g.mod == 'proper_g') g <- data.frame(g=growth$g.proper,gR = growth$gR.proper)
if(g.mod == 'g_1') g <- data.frame(g=growth$g/growth$g,gR = growth$gR/growth$gR)

# And now we need to get the survey year right, note that same has the months in numbers not month names, this is already taken care of...
#Sab.fish$survey.year[Sab.fish$month %in% 1:5] <- Sab.fish$survey.year[Sab.fish$month %in% 1:5]-1
#write.csv(growth,"D:/Framework/SFA_25_26_2024/Model/Data/Sab_mod_input_for_freya.csv")
# mod.tmp <- read.csv("D:/Framework/SFA_25_26_2024/Model/Data/Sab_mod_input_for_freya.csv")
# mod.tmp$g[1:29] - growth$g[1:29]

# Subset the fishery data to the correct years. We need to take on next year catch too..)
Sab.fish <- Sab.fish %>% dplyr::filter(survey.year %in% years) # c(years,(max(years)+1))
#tst <- Sab.fish.sf %>% dplyr::filter(survey.year %in% 1993:2014)
#rems <- tst %>% dplyr::group_by(year) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))
Sab.fish.sf <- st_as_sf(Sab.fish,coords = c("lon","lat"),remove =F, crs = 4326)
Sab.fish.sf <- Sab.fish.sf %>% st_transform(crs= c_sys)

# Now lets clip this to be data inside of our Sab boundary.
Sab.fish.sf <- st_intersection(Sab.fish.sf,Sab.shape)
# Check removals each fishing year calculated using this data, 1995 seems funny, but it's just how the timing of the fishery landed, it's fine
Sab.fish.by.year <- Sab.fish.sf %>% dplyr::group_by(year) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))
Sab.fish.by.survey.year <- Sab.fish.sf %>% dplyr::group_by(survey.year) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))
# So this looks reasonable in the most recent years, but I probably need to check the early years to see if we are missing any of the removals, from above check (only 12 points removed) it 
# seems like we might be fine, but need to check against our historical Removals estimates...
#tail(Sab.fish.by.year)
#tail(Sab.fish.by.survey.year)
# OK, so now let's see if we can use the catch knot thing Raph made to split this up withing the Sab domain
#We just need 3 columns for this
catch.sf <- Sab.fish.sf %>% dplyr::select(pro.repwt,survey.year)
names(catch.sf) <- c("Catch","Year","geometry")
 
# For the moment we need to have this starting at year 1.
catch.sf$Year <-  catch.sf$Year - (min(years)-1)
# Get rid of any 0s (due to survey year fun things)
#catch.sf$Catch <- catch.sf$Catch
#catch.sf$geometry <- catch.sf$geometry/1000

# A map...
# b.map <- pecjector(area= "Sab",c_sys = c_sys,add_layer = list(land = 'grey',eez = 'eez' , nafo = 'main',sfa = 'offshore',survey = c("offshore","outline")),txt.size=8,axes = "DM")
# Sab.fish.map <- b.map + geom_sf(data = Sab.fish.sf) + facet_wrap(~year) + geom_sf(data= Sab.shape,fill = NA)
# Sab.fish.map

# Set up our mesh...
#Sab.mesh <- setup_mesh(catch.sf,model_bound = Sab.shape,nknot=8,seed=20) # Seeds 20 and 66 work
#Sab.shape$geometry <-  Sab.shape$geometry/1000
#st_crs(Sab.shape) <- 32619
#mod.input.sf$geometry <-  mod.input.sf$geometry/1000
#st_crs(mod.input.sf) <- 32619

Sab.mesh <- setup_mesh(mod.input.sf,model_bound = Sab.shape,nknot=num.knots, max.edge = c(8,20),cutoff=2.5,seed=23) 
Sab.mesh.sf <- inla.mesh2sf(Sab.mesh$mesh)
Sab.mesh.sf$triangles$geometry <- Sab.mesh.sf$triangles$geometry*1000
Sab.mesh.sf$vertices$geometry <- Sab.mesh.sf$vertices$geometry*1000
st_crs(Sab.mesh.sf$triangles) <- c_sys
st_crs(Sab.mesh.sf$vertices) <- c_sys
knots.sf <- st_as_sf(as.data.frame(Sab.mesh$knots$centers), coords = c("X","Y")) 
knots.sf$geometry <- knots.sf$geometry*1000
st_crs(knots.sf) <- c_sys

# Plot the mesh
#ggplot(Sab.mesh.sf$triangles) + geom_sf() + geom_sf(data= Sab.shape,fill = NA,color = 'blue',size=2) + geom_sf(data = knots.sf,fill = NA)
# Now make the prediction grid
pred.grid<-setup_pred_grid(knots=Sab.mesh$knots,model_bound=Sab.mesh$utm_bound)
st_crs(pred.grid$grid) <- c_sys
# Plot the grid
#ggplot(pred.grid$grid) + geom_sf(aes(fill = as.factor(knotID))) + scale_fill_viridis_d()
# Get the knots on the right scale
#knots.on.right.scale <- Sab.mesh$knots
#knots.on.right.scale$centers <- knots.on.right.scale$centers*1000

#Sebdam catches
catchy <- catch_spread(catch = catch.sf,knots = Sab.mesh$knots)
catchy$sum_catches[,ncol(catchy$sum_catches)+1] <- 0 # set final year catches to 0
# For now we need to toss the first column from there
#catchy$density_catches <- catchy$density_catches[,-1]
#catchy$sum_catches <- catchy$sum_catches[,-1]
catchy

#TLM catch 
catch.tlm <- catch.sf %>% group_by(Year,.drop=F) %>% dplyr::summarise(catch = sum(Catch,na.rm=T))
catch.tlm[nrow(catch.tlm)+1,] <- catch.tlm[nrow(catch.tlm),]
catch.tlm$Year[nrow(catch.tlm)] <- max(catch.tlm$Year) + 1
catch.tlm$catch[nrow(catch.tlm)] <- 0

# ggplot(mod.input.sf) + geom_boxplot(aes(x=Year,y=I,group=Year)) + scale_y_log10()
# ggplot(mod.input.sf) + geom_boxplot(aes(x=Year,y=IR,group=Year))+ scale_y_log10()
# ggplot(mod.input.sf) + geom_boxplot(aes(x=Year,y=N,group=Year))+ scale_y_log10()
# ggplot(mod.input.sf) + geom_boxplot(aes(x=Year,y=L,group=Year))+ scale_y_log10()


# Fart around with inputs and see if that helps...
#mod.input.sf$L <- 3*mod.input.sf$L # This seemed to help when running with 10 years, but blew up with 20
#mod.input.sf$L[mod.input.sf$L > mod.input.sf$N] <- 0.5*mod.input.sf$N[mod.input.sf$L > mod.input.sf$N]
# Maybe I'm seeing too many recruits??
#mod.input.sf$IR <- mod.input.sf$IR/2.5
# What happens if I triple the number of Fully Recruited?
#mod.input.sf$I <- 3* mod.input.sf$I
# Here I try to make the IR's a fraction of the I in the following year, which may/may not be overly complicating the issue...
# for(i in 1:NY)
# {
#   if(i < NY &  years[i] != 2015) mod.input.sf$IR[mod.input.sf$year == years[i]] <- runif(mod.input.sf$IR[mod.input.sf$year == (years[i])],0.01,0.1)*mean(mod.input.sf$I[mod.input.sf$year == (years[i]+1)])
#   if(i == NY) mod.input.sf$IR[mod.input.sf$year == years[i]] <- runif(mod.input.sf$IR[mod.input.sf$year == (years[i])],0.01,0.1)*mean(mod.input.sf$I[mod.input.sf$year == (years[i])])
# }
#SEBDAM version 
if(mod.select == "SEAM")
{
  set_data<-data_setup(data=mod.input.sf,growths=data.frame(g = g$g,gR = g$gR),catch=as.data.frame(catchy$sum_catches),
                       model="SEBDAM",mesh=Sab.mesh$mesh,obs_mort=T, prior=TRUE,prior_pars=c(20,40),
                       mult_qI=T,spat_approach="spde",
                       knot_obj=Sab.mesh$knots,knot_area=pred.grid$area,separate_R_aniso =T,
                       all_se=T,weighted_mean_m = T)
  # So this will fix the mean value of m0 to be whatever the intial value is set at.  Let's see what happens!
   set_data$par$log_m0 <- log(init.m) # 0 = 1, 
  # #set_data$par$log_R0 <- l.init.R # 5.3 = 200, 5 = 148, 4 = 55, 5.9915 = 400, 4.606 = 100
  set_data$par$log_qR <- log(qR)
  # #set_data$map <-list(log_m0=factor(NA))
  set_data$map <-list(log_m0=factor(NA),log_qR = factor(NA))
  #set_data$map <-list(log_qR = factor(NA))
#set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA),log_qR = factor(NA))
}
#TLM version, Dude is this ever sensitive to the q priors! (5,12) actually looks solid in terms of results... maybe we can get so lucky with SEBDAM :-)
# Note that the Catch time series should be 1 year longer than the survey data here!!
if(mod.select == "TLM")
{
  set_data<-data_setup(data=as.data.frame(mod.input.sf),growths=data.frame(g = g$g,gR = g$gR),
                       #data=as.data.frame(mod.input.sf),growths=data.frame(g = growth$g.alt,gR = growth$gR.alt),
                       #data=as.data.frame(mod.input.sf),growths=data.frame(g = growth$g.proper,gR = growth$gR.proper),
                        catch=catch.tlm$catch, model="TLM",obs_mort=TRUE,prior=TRUE)
  set_data$par$log_q_R <- log(qR) 
  set_data$map <-list(log_q_R=factor(NA))
  }
  # Run the model
  mod.fit<-fit_model(set_data,silent=F)

# Now save the results appropriately
if(mod.select == "SEAM") 
{
  # m0.par <- signif(exp(set_data$par$log_m0),digits=2)
  # qR.par <- signif(exp(set_data$par$log_qR),digits=2)
  # 
  scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",init.m,"_qR_",qR,"_",num.knots,"_knots_",g.mod)
  
  saveRDS(mod.fit,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
  saveRDS(Sab.mesh,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
  saveRDS(pred.grid,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
  saveRDS(mod.input.sf,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,"_model_input.Rds"))
  
}

if(mod.select == "TLM") 
{
  #qR.par <- signif(exp(set_data$par$log_q_R),digits=2)
  scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR,"_",g.mod)
  
  saveRDS(mod.fit,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
}

# unique(names(mod.fit$sdrep$value))
# params<-get_parameters(mod.fit)
# str(params)
# rownames(params)
# #write.csv(params,"D:/Github/BBn_model/Results/Model/Sab_parms.csv")
# 
# pred.proc<-get_processes(mod.fit)
# str(pred.proc)
# summary(pred.proc$densities$m)
# 
# plot(c(years,max(years)+1),pred.proc$totals$totB, type = 'b',xlab="",ylab="Total Biomass (tonnes)")
# plot(c(years,max(years)+1),pred.proc$totals$totR, type = 'b',xlab="",ylab = "Recruits")
# plot(c(years,max(years)+1),pred.proc$totals$mean_m, type = 'b',xlab="",ylab = "Natural mortality")
# #TLM
# plot(c(years,max(years)+1),pred.proc$process$B, type = 'b',xlab="",ylab="Total Biomass (tonnes)")
# plot(c(years,max(years)+1),pred.proc$process$R, type = 'b',xlab="",ylab = "Recruits")
# plot(c(years,max(years)+1),pred.proc$process$m, type = 'b',xlab="",ylab = "Natural mortality")
#summary(pred.proc$process$m)

# Ideas from Raph... perhaps a knot with a lot of recruit 0's could be an issue
# If we have a knot with a 0 in the first year, that could send the model off to la-la land 
# Increasing the L/N ratio didn't help with S/m estimates
# Incresing the L/N ratio and decreasing the number of Recruits didn't help with S/m estimates
# decreasing the number of recruits by >50% seems to make the qR estimate much more sensisble
# Tripling the biomass of fully recruited also helps, but still nothing is working with the natural mortality.
# changing the length of the time series hasn't helped with much of anything, but does give unstable results.
# Using the SPDE model hasn't seemed to help
# Increasing FR biomass artificially by a factor of 3 does nothing, combining that with a more realistic 
# natural 
# Using a fixed m of 0.1 starting in 1986,1991, and 1995 all results in the qR being > 1, 1997 was false convergences with these settings.
# Moving to 1999 make qR kinda reasonable, bit high (about 0.55) but reasonable nonetheless.
# Going to 2002 and qR gets worse again... I can't get model to go with 1998 data for some reason...
# Using 2000 and qR goes above 1, start year is very sensistive.
# I tried going from 1999 with the recruits being 70-90 mm and the FR being 90+, that didn't do anything useful, qR was > 1.8.
#  I tried going from 1999 with recruits being 50-80 and FR being 80+ and..... nothing good happened.
# What does work is fixing m to be relatively high (around 0.3) and lowering the number of knots, 4 and 10 seem to work
# Other numbers end up with the 0 recruits problem.
# So here are my initial suite of models to compare.


################################################### End the initial model runs ###########################################
################################################### End the initial model runs ###########################################
################################################### End the initial model runs ###########################################

##################### Now load the model and make the figures! ##############################################

atow<-800*2.4384/10^6 # area of standard tow in km2
num.knots <- 10 # 4, 8, or 10
R.size <- 75
FR.size <- 90
#RO <- 100 # 100, 250, or 500
init.m <- 0.2 # log(0.05) # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05
qR  <- 0.33
years <- 1994:2022
NY <- length(years)
c_sys <- 32620
theme_set(theme_few(base_size = 22))
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
mod.select <- "SEAM"
g.mod <- 'g_original'
#g.mod <- 'g_1'
#g.mod <- 'alt_g'
#g.mod <- 'proper_g'
################################################### End the initial model runs ###########################################
### Make the figures for the models



if(mod.select != "TLM") scenario.select <-   scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",init.m,"_qR_",qR,"_",num.knots,"_knots_",g.mod)

if(mod.select == "TLM") scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR,"_",g.mod)


mod.fit <- readRDS(paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
if(mod.select != "TLM") catchy <- mod.fit$obj$env$data$C*mod.fit$obj$env$data$area # Get this into tonnes from catch density.
if(mod.select == "TLM") catchy <- mod.fit$obj$env$data$C

# This is only needed for SEAM.
if(mod.select != "TLM")
{
  
pred.grid <- readRDS(paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
# Now set up to run the figures
# Now set up to run the figures
matYear<-c(rep(c(years,(max(years)+1)),each=num.knots))
matYear1<-c(rep(years,each=num.knots))
knots<-rep(1:num.knots,NY+1)
knots1<-rep(1:num.knots,NY)

grid.gis <- pred.grid$grid
grid.gis$geometry <- grid.gis$geometry*1000
st_crs(grid.gis) <- 32620
# Now simplify the grid for the spatial plots...
knot.gis <- aggregate(grid.gis, list(grid.gis$knotID), function(x) x[1])


# Get the spatial data output
B<-data.frame(B=as.vector(mod.fit$report$B),Year=matYear,knotID=knots)
B.dat.plot<-left_join(knot.gis,B,by=c("knotID"))
 
# Recruits
R<-data.frame(R=as.vector(mod.fit$report$R),Year=matYear1, knotID=knots1)
R.dat.plot<-left_join(knot.gis,R,by=c("knotID"))
#Natural Mortality
m<-data.frame(m=as.vector(mod.fit$report$m),Year=matYear,knotID=knots)
m.dat.plot<-left_join(knot.gis,m,by=c("knotID"))
#Spatial q's
qI<-data.frame(qI=as.vector(mod.fit$report$qI),knotID=unique(knots))
q.dat.plot<-left_join(knot.gis,qI,by=c("knotID"))

# Explotation
# Now lets try and make a spatial exploitation rate plot. I'm not 1000% sure this is how I want to calculate this, but I think it is... sort it out :-)
F.dat<-data.frame(B=as.vector(mod.fit$report$areaB[,-ncol(mod.fit$report$areaB)]/1000),
                  C = c(as.vector(as.matrix(catchy[,-c((ncol(catchy)-1),ncol(catchy))])),rep(NA,num.knots)), Year=matYear1, knotID=knots1)
F.dat <- F.dat %>% dplyr::mutate(exploit = C/(B+C)) # Sticking with how offshore does this (C/(B+C)) C/B or some variant may be more realistic
F.dat.plot<-left_join(knot.gis,F.dat,by=c("knotID"))

F.dat.plot <- F.dat.plot %>% dplyr::filter(Year != max(years) )
# To get a weighted m time series combing the B, R, and the m data above
bmr <- data.frame(B=as.vector(mod.fit$report$areaB/1000),
                  m=as.vector(mod.fit$report$m),
                  R=c(as.vector(mod.fit$report$areaR/1000),rep(NA,num.knots)),
                  Year=matYear, knotID=knots)
# Now get an m estimate for the whole area for B, R, and both, I think easiest will be to get totB and totR by year in here...
tBR <- bmr %>% group_by(Year) %>% dplyr::summarise(totB = sum(B,na.rm=T),
                                                   totR = sum(R,na.rm=T),
                                                   totSSB = sum(B,na.rm=T) + sum(R,na.rm=T))

# Now merge that with bmr...
bmr <- left_join(bmr,tBR,by="Year")

# Now getting weighted means by year should be easy...
nat.mat <- bmr %>% group_by(Year) %>% dplyr::summarise(m.FR = sum(B*m/totB),
                                                       m.R = sum(R*m/totR),
                                                       m.all = sum((R+B)*m/totSSB))
nat.mat$Raph.m <- mod.fit$report$mean_m
# Now make them long...
nat.mat.plt <- pivot_longer(nat.mat,!Year,names_to = "method",values_to = 'm')
nat.mat.plt$method <- factor(nat.mat.plt$method,levels = c("m.all","Raph.m","m.FR","m.R"),labels = c("Weighted","Unweighted","Fully Recruited","Recruits"))



# If not already loaded, make a map
#b.map <- pecjector(area= "BBn",c_sys = c_sys,add_layer = list(land = 'grey',eez = 'eez' , nafo = 'main',sfa = 'offshore',survey = c("offshore","outline")),txt.size=8,axes = "DM")

# Final piece is that we can compare the biomass knot by knot to what would be expected based on B/R/m/g/C alone.  This is a bit messy but this does the trick...
gs <- mod.fit$obj$env$data$gI
gRs <- mod.fit$obj$env$data$gR
tmp <- NULL

for(i in 1:(NY-1))
{
  Bs.tot <- mod.fit$report$B*mod.fit$obj$env$data$area/1000
  Rs.tot <- mod.fit$report$R*mod.fit$obj$env$data$area/1000
  ms <- mod.fit$report$m
  Bst <- (exp(-ms[,i+1]))*gs[i]*(Bs.tot[,i]-catchy[,i]) 
  Rst <- (exp(-ms[,i+1]))*gRs[i]*(Rs.tot[,i]) 
  B2 <- Bst + Rst
  if(any(B2 < 0)) {B2[B2<= 0] <- 5; print("HEADS UP!! You have a negative B2 estimate (look for the 5s in B.exp).")}
  B.next <- Bs.tot[,i+1]
  B.diff <- B2 - B.next 
  B.per.diff <- 100*((B2 - B.next)/B.next)
  n.not <- length(B2)
  
  if(i == 1) tmp[[as.character(years[i])]] <- data.frame(knot = 1:n.not,Year = rep(years[i],n.not),B.exp = NA,B.mod = NA,B.fr = NA,B.rec = NA,B.diff = NA,B.per.diff = NA,
                                                         m = NA,C = NA,g = NA,gR = NA)
  tmp[[as.character(years[i+1])]] <- data.frame(knot = 1:n.not,Year = rep(years[i+1],n.not),B.exp = B2,B.mod = B.next,B.fr = Bst,B.rec = Rst,B.diff = B.diff,B.per.diff = B.per.diff,
                                                m = ms[,i+1],C = catchy[,i],g = rep(gs[i],n.not),gR = rep(gRs[i],n.not))
}

# Unpack that object
B.diff.comp <- do.call('rbind',tmp)

Bdiff.comp.dat<-data.frame(B.diff=as.vector(B.diff.comp$B.diff),B.per.diff = as.vector(B.diff.comp$B.per.diff),
                           B.exp=as.vector(B.diff.comp$B.exp),B.mod = as.vector(B.diff.comp$B.mod),
                           B.fr=as.vector(B.diff.comp$B.fr),B.rec = as.vector(B.diff.comp$B.rec),
                           Year=matYear1,knotID=knots1)
Bdiff.comp.dat.plot<-left_join(knot.gis,Bdiff.comp.dat,by=c("knotID"))
Bdiff.comp.dat.plot <- Bdiff.comp.dat.plot %>% dplyr::filter(Year != min(years))
# Save this so I can compare the 'miss' across models...

# Smaller text for spatial figures.
theme_set(theme_few(base_size = 14))

#Spatial predictions
#B
b.brk <- pretty(log(B.dat.plot$B))
b.lab <- signif(exp(b.brk),digits=2)

spatial.B.plot<- ggplot() + geom_sf(data=B.dat.plot,aes(fill=log(B)),color='grey')+
                            facet_wrap(~Year)+ 
                            scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
                            scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
                            scale_fill_viridis_c(breaks = b.brk, labels=b.lab,name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
                            theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_biomass.png"),spatial.B.plot,base_width = 10,base_height = 10)
# Remove missing survey years
spatial.B.plot<- ggplot() + geom_sf(data=B.dat.plot %>% dplyr::filter(!Year %in% c(2015,2020,2023)),aes(fill=log(B)),color='grey')+
  facet_wrap(~Year)+ 
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  scale_fill_viridis_c(breaks = b.brk, labels=b.lab,name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_biomass_no_missing_surveys.png"),spatial.B.plot,base_width = 10,base_height =7)
# Subest to 4 years
spatial.B.plot<- ggplot() + geom_sf(data=B.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(B)),color='grey')+
  facet_wrap(~Year)+ 
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  scale_fill_viridis_c(breaks = b.brk, labels=b.lab, name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_biomass_4_years.png"),spatial.B.plot,base_width = 10,base_height =7)



#R
r.brk <- pretty(log(R.dat.plot$R))
r.lab <- signif(exp(r.brk),digits=2)
spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot,aes(fill=log(R)),col='grey')+
                             facet_wrap(~Year)+ 
                             scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
                             scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
                             scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
                             theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_recruits.png"),spatial.R.plot,base_width = 10,base_height = 10)
# Remove missing survey years
spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot %>% dplyr::filter(!Year %in% c(2015,2020,2023)),aes(fill=log(R)),col='grey')+
                              facet_wrap(~Year)+ 
                              scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
                              scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
                              scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
                              theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_recruits_no_missing_surveys.png"),spatial.R.plot,base_width = 10,base_height = 10)
# 4 years
spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(R)),col='grey')+
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  facet_wrap(~Year)+ 
  scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_recruits_4_years.png"),spatial.R.plot,base_width = 10,base_height = 7)

#m
m.brk <- log(c(0.03,0.1,0.3,1))
#m.brk <- pretty(log(m.dat.plot$m))
m.lab <- signif(exp(m.brk),digits=2)

spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot,aes(fill=log(m)),color='grey')+
                              facet_wrap(~Year)+
                              scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
                              scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
                              scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality",option = "B",direction =1,begin = 0.2,end=1) + theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_mort.png"),spatial.m.plot,base_width = 10,base_height = 10)
# Remove missing survey years
spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot %>% dplyr::filter(!Year %in% c(2015,2020,2023)),aes(fill=log(m)),color='grey')+
                              facet_wrap(~Year)+
                              scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
                              scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
                              scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality",option = "B",direction =1,begin = 0.2,end=1) + theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_mort_no_missing_surveys.png"),spatial.m.plot,base_width = 10,base_height = 10)
# 4 years
spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(m)),color='grey')+
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  facet_wrap(~Year)+ 
  scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality (Inst)",option = "B",direction =1,begin = 0.2,end=1) + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_mort_4_years.png"),spatial.m.plot,base_width = 10,base_height = 7)


# q 
spatial.q.plot <- ggplot() + geom_sf(data=q.dat.plot,aes(fill=qI),col=NA)+
                             scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
                             scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
                             scale_fill_viridis_c(name="Predicted catchability (qI)",option = "C",begin = 0.2,end =0.8) + 
                             theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_catchability.png"),spatial.q.plot,base_width = 10,base_height = 10)


# OK, so lets try an make a map of the spatial exploitation rates, not sure if this is all correct yet.
F.dat.plot$exp.na <- NA
F.dat.plot$exp.na[F.dat.plot$exploit != 0] <- F.dat.plot$exploit[F.dat.plot$exploit != 0]

e.brk <- log(c(0.00015,0.001,0.005,0.02,0.08))
#e.brk <- pretty(log(F.dat.plot$exp.na))
e.lab <- signif(exp(e.brk),digits=2)

spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot,aes(fill=log(exp.na)),color='grey') +
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  facet_wrap(~Year) + scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Spatial_exploit.png"),spatial.exploit.plot,base_width = 10,base_height = 10)
# Remove missing survey years
spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot %>% dplyr::filter(!Year %in% c(2015,2020,2023)),aes(fill=log(exp.na)),color='grey') +
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  facet_wrap(~Year) + scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Spatial_exploit_no_missing_surveys.png"),spatial.exploit.plot,base_width = 10,base_height = 10)
# 4 years
spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(exp.na)),color='grey') +
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  facet_wrap(~Year) + 
  scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Spatial_exploit_4_years.png"),spatial.exploit.plot,base_width = 10,base_height = 10)


bd.brk <- pretty(Bdiff.comp.dat.plot$B.diff)
bd.lab <- bd.brk#signif(exp(b.brk),digits=2)

spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot,aes(fill=B.diff),color='grey')+
  facet_wrap(~Year)+ 
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_BDiff.png"),spatial.Bdiff.plot,base_width = 10,base_height = 10)
# Remove missing survey years
spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(!Year %in% c(2015,2020,2023)),aes(fill=B.diff),color='grey')+
  facet_wrap(~Year)+ 
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_BDiff_no_missing_surveys.png"),spatial.Bdiff.plot,base_width = 10,base_height = 10)
# 4 years
spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=B.diff),color='grey')+
  facet_wrap(~Year)+ 
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_BDiff_4_years.png"),spatial.Bdiff.plot,base_width = 10,base_height = 7)


# Same plot but the percentage miss by cell.
pb.brk <- pretty(Bdiff.comp.dat.plot$B.per.diff)
pb.lab <- pb.brk#signif(exp(b.brk),digits=2)
spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot,aes(fill=B.per.diff),color='grey')+
  facet_wrap(~Year)+ 
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_per_BDiff.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 10)
# Remove missing survey years
spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(!Year %in% c(2015,2020,2023)),aes(fill=B.per.diff),color='grey')+
  facet_wrap(~Year)+ 
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_per_BDiff_no_missing_surveys.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 10)
# 4 years
spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=B.per.diff),color='grey')+
  facet_wrap(~Year)+ 
  scale_x_continuous(breaks = c(-61.66667,-61), labels = c("61° 40'W","61° W")) +
  scale_y_continuous(breaks = c(43.33333, 43.666667, 44),labels = c("43°20'N","43°40'N","44°N")) +
  scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
  theme(axis.text.x=element_text(angle=-45,hjust=0))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Spatial_per_BDiff_4_years.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 7)

theme_set(theme_few(base_size = 22))


# And now plot them....
# m.comp <- ggplot(nat.mat.plt,aes(x=Year,y=m,group = method,color=method)) + geom_line(linewidth=1.5)  + 
#   xlab("") + ylab("Natural Mortality (instantaneous)")+
#   scale_color_manual(values = c("blue","orange","grey","black"))  + 
#   theme(legend.title=element_blank())
# save_plot(paste0(repo.loc,"Results/Figures/Sab/",mod.select,"_",scenario.select,"/natural_mortality_comparisons.png"),m.comp,base_height = 6,base_width = 10)
# # Remove missing survey years
# m.comp <- ggplot(nat.mat.plt %>% dplyr::filter(Year < 2015)) + geom_line(aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  + 
#   geom_line(data = nat.mat.plt %>% dplyr::filter(Year %in% 2016:2019), aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  + 
#   geom_line(data = nat.mat.plt %>% dplyr::filter(Year %in% 2021:2022), aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  +  
#   xlab("") + ylab("Natural Mortality (instantaneous)")+
#   scale_color_manual(values = c("blue","orange","grey","black"))  + 
#   theme(legend.title=element_blank())
# save_plot(paste0(repo.loc,"Results/Figures/Sab/",mod.select,"_",scenario.select,"/natural_mortality_comparisons_no_missing_surveys.png"),m.comp,base_height = 6,base_width = 10)
}

pred.proc <- get_processes(mod.fit)

if(mod.select == "TLM")
{
  pred.proc$log_processes$year <- c(years,max(years)+1)
  pred.proc$log_processes$totB.LCI <- exp(pred.proc$log_processes$log_B - 1.96*pred.proc$log_processes$se_log_B)
  pred.proc$log_processes$totB.UCI <- exp(pred.proc$log_processes$log_B + 1.96*pred.proc$log_processes$se_log_B)
  pred.proc$log_processes$totR.LCI <- exp(pred.proc$log_processes$log_R - 1.96*pred.proc$log_processes$se_log_R)
  pred.proc$log_processes$totR.UCI <- exp(pred.proc$log_processes$log_R + 1.96*pred.proc$log_processes$se_log_R)
  pred.proc$log_processes$m.LCI <- exp(pred.proc$log_processes$log_m - 1.96*pred.proc$log_processes$se_log_m)
  pred.proc$log_processes$m.UCI <- exp(pred.proc$log_processes$log_m + 1.96*pred.proc$log_processes$se_log_m)
  pred.proc$log_processes <- as.data.frame(pred.proc$log_processes)
}

if(mod.select != "TLM")
{
  # Get the overall estimates + the 95% CI
  pred.proc$log_processes$year <- c(years,max(years)+1)
  pred.proc$log_processes$log_B <- pred.proc$log_tot_frame$log_totB
  pred.proc$log_processes$log_R <- pred.proc$log_tot_frame$log_totR
  pred.proc$log_processes$log_m <- pred.proc$log_tot_frame$log_mean_m
  pred.proc$log_processes$totB.LCI <- exp(pred.proc$log_tot_frame$log_totB - 1.96*pred.proc$log_tot_frame$se_log_totB)
  pred.proc$log_processes$totB.UCI <- exp(pred.proc$log_tot_frame$log_totB + 1.96*pred.proc$log_tot_frame$se_log_totB)
  pred.proc$log_processes$totR.LCI <- exp(pred.proc$log_tot_frame$log_totR - 1.96*pred.proc$log_tot_frame$se_log_totR)
  pred.proc$log_processes$totR.UCI <- exp(pred.proc$log_tot_frame$log_totR + 1.96*pred.proc$log_tot_frame$se_log_totR)
  pred.proc$log_processes$m.LCI <- exp(pred.proc$log_tot_frame$log_mean_m - 1.96*pred.proc$log_tot_frame$se_log_mean_m)
  pred.proc$log_processes$m.UCI <- exp(pred.proc$log_tot_frame$log_mean_m + 1.96*pred.proc$log_tot_frame$se_log_mean_m)
  pred.proc$log_processes <- as.data.frame(pred.proc$log_processes)
}
# SEBDAM Version
# Annual explotation
if(mod.select != "TLM") catch.annual <- data.frame(totC = colSums(catchy[,-ncol(catchy)]), Year = years)
if(mod.select == "TLM") catch.annual <- data.frame(totC = catchy[-length(catchy)], Year = years)

pred.proc$log_processes <- pred.proc$log_processes %>% dplyr::filter(year < 2023)

# Not 1000% sure this is correct at this point.  It is noteworth how different this is than the map though...
if(mod.select != "TLM")
{
  ann.exploit <- data.frame(year = years,B = exp(pred.proc$log_processes$log_B), Catch = c(colSums(catchy[,-c((ncol(catchy)-1),ncol(catchy))]),NA),
                            B.LCI = pred.proc$log_processes$totB.LCI, B.UCI = pred.proc$log_processes$totB.UCI)
}
if(mod.select == "TLM")
{
  ann.exploit <- data.frame(year = years,B = exp(pred.proc$log_processes$log_B), Catch = c(NA,catchy[1:(length(catchy)-2)]),
                            B.LCI = pred.proc$log_processes$totB.LCI, B.UCI = pred.proc$log_processes$totB.UCI)
}
ann.exploit$exploit <- c(ann.exploit$Catch[1:(nrow(ann.exploit)-1)]/(ann.exploit$B[2:nrow(ann.exploit)]+ann.exploit$Catch[1:(nrow(ann.exploit)-1)]),NA)
ann.exploit$exploit.UCI <- c(ann.exploit$Catch[1:(nrow(ann.exploit)-1)]/(ann.exploit$B.UCI[2:nrow(ann.exploit)]+ann.exploit$Catch[1:(nrow(ann.exploit)-1)]),NA)
ann.exploit$exploit.LCI <- c(ann.exploit$Catch[1:(nrow(ann.exploit)-1)]/(ann.exploit$B.LCI[2:nrow(ann.exploit)]+ann.exploit$Catch[1:(nrow(ann.exploit)-1)]),NA)
ann.exploit$FM <- 1-exp(-ann.exploit$exploit)
ann.exploit$FM.LCI <- 1-exp(-ann.exploit$exploit.LCI)
ann.exploit$FM.UCI <- 1-exp(-ann.exploit$exploit.UCI)



#pred.proc$log_processes <- pred.proc$log_processes %>% dplyr::filter(year < 2023)


# Biomass time series
bm.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + ylim(c(0,1.5e4)) + 
  xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Biomass_time_series.png"),bm.ts.plot,base_width = 11,base_height = 8.5)
# Recruit time seris
rec.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Recruit Biomass (tonnes)") + ylim(c(0,1.8e3)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Recruit_time_series.png"),rec.ts.plot,base_width = 11,base_height = 8.5)
# Natural mortality time series...
mort.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Natural mortality (Instantaneous)") + ylim(c(0,1)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_nat_mort_time_series.png"),mort.ts.plot,base_width = 11,base_height = 8.5)
# Explotation Rate Time Series
exploit.plot <- ggplot(ann.exploit) + geom_line(aes(x=year,y=exploit),linewidth = 1.5) + geom_ribbon(aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') +
  xlab("") + ylab("Exploitation Rate (Proportional)") + ylim(c(0,0.1)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_exploit_time_series.png"),exploit.plot,base_width = 11,base_height = 8.5)

# Same plots but removing missing survey years....
bm.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2015))) + 
  geom_line(aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  ylim(c(0,1.5e4)) + xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Biomass_time_series_no_missing_surveys.png"),bm.ts.plot,base_width = 11,base_height = 8.5)
# Recruit time series
rec.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2015))) + 
  geom_line(aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Recruit Biomass (tonnes)") + ylim(c(0,1.8e3)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Recruit_time_series_no_missing_surveys.png"),rec.ts.plot,base_width = 11,base_height = 8.5)
# Natural mortality time series...
mort.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2015))) + 
  geom_line(aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Natural mortality (Instantaneous)") + ylim(c(0,1)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_nat_mort_time_series_no_missing_surveys.png"),mort.ts.plot,base_width = 11,base_height = 8.5)


# Explotation Rate Time Series
exploit.plot <- ggplot(ann.exploit%>% dplyr::filter(year < c(2015))) + geom_line(aes(x=year,y=exploit),linewidth = 1.5) +
  geom_ribbon(aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') +
  geom_line(data = ann.exploit%>% dplyr::filter(year %in% 2016:2019), aes(x=year,y=exploit),linewidth = 1.5) + 
  geom_ribbon(data = ann.exploit%>% dplyr::filter(year %in% 2016:2019),aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') +
  geom_line(data= ann.exploit%>% dplyr::filter(year %in% 2021:2022),aes(x=year,y=exploit),linewidth = 1.5) +
  geom_ribbon(data= ann.exploit%>% dplyr::filter(year %in% 2021:2022),aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') +
  xlab("") + ylab("Exploitation Rate (Proportional)") + ylim(c(0,0.1)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_exploit_time_series_no_missing_surveys.png"),exploit.plot,base_width = 11,base_height = 8.5)

# Save out a bunch of objects...
# The biomass difference
saveRDS(Bdiff.comp.dat.plot,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                                   mod.select,"_model_output_",scenario.select,"_B_differnce.Rds"))
# Annual exploitation
saveRDS(ann.exploit,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                                   mod.select,"_model_output_",scenario.select,"_annual.exploit.Rds"))
# The nicely summarized model object
saveRDS(pred.proc,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                           mod.select,"_model_output_",scenario.select,"_pred_proc.Rds"))
# The summarized spatial bits.
saveRDS(F.dat.plot,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                         mod.select,"_model_output_",scenario.select,"_F_spatial.Rds"))
# q spatial
saveRDS(q.dat.plot,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                          mod.select,"_model_output_",scenario.select,"_q_spatial.Rds"))
# m spatial
saveRDS(m.dat.plot,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                          mod.select,"_model_output_",scenario.select,"_m_spatial.Rds"))
# R spatial
saveRDS(R.dat.plot,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                          mod.select,"_model_output_",scenario.select,"_R_spatial.Rds"))
# B spatial
saveRDS(B.dat.plot,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                          mod.select,"_model_output_",scenario.select,"_B_spatial.Rds"))


