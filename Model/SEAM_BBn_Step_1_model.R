# So here we'll try and get all the data in the correct structure for the spatial model.
#devtools::install_github("RaphMcDo/SEBDAM@dev-branch") # to install a specific branch of SEBDAM
# Process flow for this framework...
#Step 1: Run the SEBDAM ("Bnk"_SEBDAM.R) model, control the modelling options (initial natural mortality (init.m), recruit catchability (qR), number of knots (num.knots))
#        and the growth model being used.  I have 3 of these currently, g_original is the method currently being used.  alt_g is a simple method in which we
#        look at observed mw in the following year, exclude the effect of incoming recruits, and then figure out how much the FR scallop from last year must have grown.
#        proper_g is the final methods, it is more in depth than the other methods as it looks at the mw of scallop > 105 mm in the current year and compares that
#        to the mw of scallop 90+ in the previous year, given our von B curve the 90+ scallop in the previous year should all be 105mm+ in the following year
#        so by doing this we have an estimate of the realized growth of all the 90+ scallop from last year.
#   I am using data from 1994 to 2022 and we have decided on recruits being 75-90 mm in size and FR being 90+ mm

# Step2: Run the retrospective analyses ('Bnk'_Retro.R) for one (or more) of these SEBDAM models. The retros are looking good for this model
# Step3: Run the Prediction Evaluation Simulations ("Bnk_prediction_evaluation_simulations.R) to see how well your model predicts the following year, this looks at 8 different ways of projecting the biomass for the 
#        following year.
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
library(shape)

library(ggthemes)

theme_set(theme_few(base_size = 22))

source("D:/Github/Assessment_fns/Fishery/logs_and_fishery_data.r")
source("D:/Github/Assessment_fns/Maps/pectinid_projector_sf.R")
source("D:/Github/Assessment_fns/Maps/convert_inla_mesh_to_sf.R")
source("D:/Github/SEBDAM/R/data_setup.R")
########################################################################################################
# Bring in the data and tidy it up for the analysis

bbn.shape <- st_read("D:/Github/GIS_layers/survey_boundaries/BBn.shp", quiet=T)
bbn.shape <- bbn.shape %>% st_transform(crs = 32619) # BBn is right on the 19/20 border so think they are basically equivalent options here
# Bring in the survey data
#load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90.Rdata")
load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90RSCS_newMWSH_GBb.RData")
#load("D:/Framework/SFA_25_26_2024/Model/Data/testing_results_framework3.Rdata")
#load("D:/testing_folder/testing_results_framework_75-90.Rdata")
#surv.dat <- surv.dat$BBn
#saveRDS(surv.dat,'D:/Github/BBn_model/Results/BBn_surv.dat.RDS')
#surv.dat <- readRDS('D:/Framework/SFA_25_26_2024/Model/Data/BBn_surv.dat.RDS')
surv.dat <- surv.dat$BBn
mod.dat <- survey.obj$BBn$model.dat
saveRDS(mod.dat,'D:/Framework/SFA_25_26_2024/Model/Data/BBn_BSSM_model.dat.RDS')
# For now we need to get a 2022 growth term (once we run the BBn model I can update this with the latest data), just recycling the growth for 2021 for the moment
#mod.dat <- rbind(mod.dat,mod.dat[nrow(mod.dat),])
#mod.dat$year[nrow(mod.dat)] <- 2022
# Bring in the fishery data
# logs_and_fish(loc="offshore",year = 1991:2022,un=un.ID,pw=pwd.ID,db.con=db.con,direct="Y:/Offshore/Assessment/", get.marfis=F)
# fish.dat<-merge(new.log.dat,old.log.dat,all=T)
# fish.dat$ID<-1:nrow(fish.dat)
# 
# fish.dat <- fish.dat[!is.na(fish.dat$lon),]
# fish.dat <- fish.dat[!is.na(fish.dat$lat),]
# fish.dat <- fish.dat[!fish.dat$lon==0,]
# fish.dat <- fish.dat[!fish.dat$lat==0,]
# #
# 
# # Now subset to BBn and add in the missing years of data
# #bbn.fish <- fish.dat %>% dplyr::filter(bank == "BBn")
# bbn.fish <- fish.dat[fish.dat$bank == "BBn",]
# bbn.fish <- bbn.fish[!is.na(bbn.fish$lon),]
# # There are 12 data points at 0,0 that we remove, I'm not worried about accounting for these 12 points!
# #bbn.fish <- bbn.fish %>% dplyr::filter(lat !=0 | lon != 0)
# bbn.fish$month <- lubridate::month(bbn.fish$date)
# # Now I want to put a 'survey year' on these because that's what we're gonna need for our modelling... start by porting over the year
# bbn.fish$survey.year <- bbn.fish$year
# # DK NOTE: Now this is going to get confusing for us and we may want to tweak SEBDAM for this, but that's a down the road job, not a playing around with model job
# # But based on the indexing in SEBDAM, I am going to change how we index the survey year data from what we have done with offshore traditionally.
# # Historically anything from the last half of the year goes into the following years, eg. survey.year 2002 = June 2001- May 2002.
# # But in SEBDAM we have (B(t-1) - C(t-1)), so let's say we have year 2002 survey biomass, this says we remove the 20002 catch from that
# # we want that catch to be the catch from June 2002 to May 2003, i.e. we remove the catch before we allow the population to grow
# # This is what we do in our current model, but we have a different index (C(t) on our model.
# # Basically survey year 2002 = June 2002 - May 2003 now
# #DK note: We probably should think more about the survey year fun and how exactly we want to handle removal of catch in our models.
# # We don't have removals for 2009, we need something for that, so we're adding that in here...
# # Note that in 2015 the survey was delayed until July, but there was no fishing on BBn in 2015 in June or July, so this system still works for the
# # survey year despite that.... happily!
# bbn.fish$survey.year[bbn.fish$month %in% 1:5] <- bbn.fish$survey.year[bbn.fish$month %in% 1:5] -1
# # Need to add fake data for 2009
# bbn.fish[nrow(bbn.fish)+1,] <- NA
# 
# bbn.fish$pro.repwt[nrow(bbn.fish)] <- 0
# bbn.fish$year[nrow(bbn.fish)] <- 2009
# # See DK NOte below
# bbn.fish$survey.year[nrow(bbn.fish)] <- 2009
# bbn.fish$month[nrow(bbn.fish)] <- 6
# # Getting fake lat/lon coords that are on BBn
# bbn.fish$lat[nrow(bbn.fish)] <- 42.85600
# bbn.fish$lon[nrow(bbn.fish)]  <- -65.90183
# saveRDS(bbn.fish,'D:/Framework/SFA_25_26_2024/Model/Data/BBn_fish.dat.RDS')
# 
bbn.fish <- readRDS('D:/Framework/SFA_25_26_2024/Model/Data/BBn_fish.dat.RDS')
bbn.fish$pro.repwt <- bbn.fish$pro.repwt/1000
#### Finished Data prep and clean up!
###############################################################################################

# Set parameters for the run...
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
mod.select <- "SEAM"
atow<-800*2.4384/10^6 # area of standard tow in km2
years <- 1994:2022
NY <- length(years)
R.size <- "75"
FR.size <- "90"
num.knots <- 20 # Going to test 10, 15, and 20
qR <- 0.33# This is for TMB (log recruit catchability) testing catchability of 0.5, test 0.3 and 0.1. Can be used with SEAM too in place of the m and R initiailztion.
#init.m <- 0.2 # This is for SEAM, sets first year natural mortality, going to test 0.8,0.4, 0.15, and 0.05
# Various explorations of the g models.
#g.mod <- 'g_original'
#g.mod <- 'alt_g'
#g.mod <- 'proper_g'
# I do we want to vary catchabilty spatially.
vary.q <- T
# The survey biomass index for 1994 says there were 249 tonnes of recruits that year.
#l.init.R <- log(5) # I think this is initializing the knots, i.e. each knots with this number of recruits.  Aim for 250 total
#qR.par <- log(0.9) # Trying with a fixed qR for SEAM instead

live.subset <- surv.dat %>% dplyr::filter(state == 'live' & random ==1)
dead.subset <- surv.dat %>% dplyr::filter(state== "dead" & random ==1)

live.input <- data.frame(I = live.subset$com.bm, IR = live.subset$rec.bm,year = live.subset$year,tow = live.subset$tow,tot.live.com = live.subset$com,lat = live.subset$lat,lon=live.subset$lon,tow_type = live.subset$random)
clap.input <- data.frame(L = dead.subset$com,tow = dead.subset$tow,year = dead.subset$year,tow_type = dead.subset$random)
mod.input <- left_join(live.input,clap.input,by=c('tow','year','tow_type'))
mod.input$N <- round(mod.input$tot.live.com + mod.input$L)
# Looks like there are no values > 0 but < 0.5, so the low clapper numbers should all round up to 1 (which makes sense as you'd only get < 0.5 if we had tows twice as long as they should be)
mod.input$L <- round(mod.input$L) 
mod.input.sf <- st_as_sf(mod.input,coords = c('lon','lat'),remove=F,crs = 4326)
mod.input.sf <- mod.input.sf %>% st_transform(crs=32619)
mod.input.sf <- st_intersection(mod.input.sf,bbn.shape)
mod.input.sf$Year <- mod.input.sf$year - (min(years)-1)
mod.input.sf$I <- mod.input.sf$I/atow
mod.input.sf$IR <- mod.input.sf$IR/atow
mod.input.sf <- mod.input.sf %>% dplyr::filter(year %in% years)

# We need to recalculate the growth data using the new von B curves and size bins... breaking out the ugly code to do this...
vonB.par <- data.frame(Linf = 164.4,K = 0.2, t0 = -0.2)
#vonB.par <- data.frame(Linf = 168.2,K = 0.189, t0 = -0.23)
#vonB.par <- data.frame(Linf = 148,K = 0.19, t0 = 0.11)


# Back to real code
# So first up, this condition is the weighted mean condition, this uses the GAM predicted scallop condition factor for each tow
# and the biomass from each tow to come up with an overall bank average condition factor.
# This is weight in this year, which becomes t-1 
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
# I'm not using this one, it's easy, but doesn't account for recruits properly
dk.growth <- NA
for(i in 2:nrow(mod.dat)) dk.growth[i-1] <- 1+ ((mod.dat$w.bar[i] - mod.dat$w.bar[i-1])/mod.dat$w.bar[i-1])
# Not going to use so leave out to avoid confusion
#mod.dat$g.new <- c(dk.growth,dk.growth[length(dk.growth)])

# DK growth model #2 is attempting to account for the recruit contribution to the MW estimates
head(mod.dat)

# The other option is to calculate the growth using w.bar of everything over 100 mm, which will be default exclude the vast majority of the
# recruits as 90 mm scallop will grow by about 12 cm, so might have a few recruits in there, but tracking the changes in that size class tells
# us what the realized growth was for the FRs that excludes the recruits
# So what we do is take the ratio of the w.bar for everything bigger than 100 mm in year 2, to the w.bar for all FR scallop in year one
# Based on the von.B the vast majority of the scallop in that ratio be the same individuals.
# So to calculate the 100 mm thing I'll need to use the shf in surv.dat...
# I can do the same with recruit growth can't I, everything from 90 to 100 were probably recruits last year
# so look at 75-90 last year and compare with 90 to 100 this year...

sizes <- seq(0.025,2,by=0.05) # So I'd be using the 1.075 bin and everything bigger
# The w.yst object is exactly proportional to mod.dat$I, there is an offset, but given I need proportions I think this object is perfectly fine to use.
mw.per.bin <- data.frame(mw.per.bin = rbind(survey.obj$BBn$shf.dat$w.yst/survey.obj$BBn$shf.dat$n.yst,rep(NA,40)),year = c(mod.dat$year,2020))
N.per.bin <- data.frame(N.per.bin = rbind(survey.obj$BBn$shf.dat$n.yst,rep(NA,40)),year = c(mod.dat$year,2020))
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

# now need to add in 2020 to mod.dat...
mod.dat.tmp <- mod.dat
mod.dat.tmp[nrow(mod.dat.tmp)+1,] <- NA
mod.dat.tmp$year[nrow(mod.dat.tmp)] <- 2020
mod.dat.tmp <- mod.dat.tmp[order(mod.dat.tmp$year),]

growth <- data.frame(year = mod.dat.tmp$year,g = mod.dat.tmp$g, gR = mod.dat.tmp$gR,
                     #g.alt = alt.g$alt.g, gR.alt = mod.dat.tmp$gR, # I don't have a good idea how to estiamte gR growth, so using the other way
                     g.proper = g.proper$g.proper,gR.proper = g.proper$gR.proper)
# Now addin the missing growth years for g and gR
growth$g[growth$year == 2020] <- median(growth$g,na.rm=T)
growth[growth$year == 2020,names(growth) %in% c("gR","gR.alt")] <- median(growth$gR,na.rm=T)
growth <- growth[which(!is.na(growth$g)),]
growth <- growth[which(!is.na(growth$g)),]
growth[nrow(growth)+1,] <- growth[nrow(growth),]
growth$year[nrow(growth)] <- max(years) + 1

growth <- growth %>% dplyr::filter(year >= min(years))

#if(g.mod == 'g_original') g <- data.frame(g=growth$g,gR = growth$gR)
#if(g.mod == 'alt_g') g <- data.frame(g=growth$g.alt,gR = growth$gR.alt)
g <- data.frame(g=growth$g.proper,gR = growth$gR.proper) #if(g.mod == 'proper_g') 
#if(g.mod == 'g_1') g <- data.frame(g=growth$g/growth$g,gR = growth$gR/growth$gR)

#write.csv(growth,"D:/Framework/SFA_25_26_2024/Model/Data/BBn_input_data_for_freya.csv")
#mod.tmp <- read.csv("D:/Framework/SFA_25_26_2024/Model/Data/BBn_input_data_for_freya.csv")
#mod.tmp$g[1:29] - growth$g[1:29]

# Now we can clip both of these to subset it to the data that I think we need for the analysis....
# If we run with random == 1 then we need to fill in 2020...
mod.input.sf[nrow(mod.input.sf)+1,] <- mod.input.sf[nrow(mod.input.sf),]
#mod.input.sf[nrow(mod.input.sf),] <- mod.input.sf[nrow(mod.input.sf)-1,]
mod.input.sf$year[nrow(mod.input.sf)] <- 2020
mod.input.sf$Year[nrow(mod.input.sf)] <- which(years == 2020)
mod.input.sf$I[nrow(mod.input.sf)] <- NA
mod.input.sf$IR[nrow(mod.input.sf)] <- NA
mod.input.sf$tot.live.com[nrow(mod.input.sf)] <- NA
mod.input.sf$L[nrow(mod.input.sf)] <- 0
mod.input.sf$N[nrow(mod.input.sf)] <- 0

mod.input.sf <- mod.input.sf[order(mod.input.sf$year),]

bbn.fish.sf <- st_as_sf(bbn.fish,coords = c("lon","lat"),remove =F, crs = 4326)
bbn.fish.sf <- bbn.fish.sf %>% st_transform(crs= 32619)


# Now lets clip this to be data inside of our bbn boundary.
bbn.fish.sf <- st_intersection(bbn.fish.sf,bbn.shape)

# Check removals each fishing year calculated using this data
bbn.fish.by.year <- bbn.fish.sf %>% dplyr::group_by(year) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))
bbn.fish.by.survey.year <- bbn.fish.sf %>% dplyr::group_by(survey.year,.drop=F) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))
# So this looks reasonable in the most recent years, but I probably need to check the early years to see if we are missing any of the removals, from above check (only 12 points removed) it 
# seems like we might be fine, but need to check against our historical Removals estimates...
#tail(bbn.fish.by.year)
#tail(bbn.fish.by.survey.year)

# Subset the fishery data as necessary
bbn.fish <- bbn.fish |>collapse::fsubset(survey.year %in% years)
bbn.fish.sf <- bbn.fish.sf |> collapse::fsubset(survey.year %in% years)
# OK, so now let's see if we can use the catch knot thing Raph made to split this up withing the BBn domain
#We just need 3 columns for this
catch.sf <- bbn.fish.sf %>% dplyr::select(pro.repwt,survey.year)
names(catch.sf) <- c("Catch","Year","geometry")
#catch.sf$geometry <- catch.sf$geometry/1000
#catch.sf |> data.frame() |> collapse::fgroup_by(Year) |> collapse::fsummarize(sum = sum(Catch))
# Set up the mesh

bbn.mesh <- setup_mesh(mod.input.sf,model_bound = bbn.shape,nknot=num.knots, max.edge = c(3,10),cutoff=2,seed=66) # Seeds 20 and 66 work
bbn.mesh.sf <- inla.mesh2sf(bbn.mesh$mesh)
bbn.mesh.sf$triangles$geometry <- bbn.mesh.sf$triangles$geometry*1000
bbn.mesh.sf$vertices$geometry <- bbn.mesh.sf$vertices$geometry*1000
st_crs(bbn.mesh.sf$triangles) <- 32619
st_crs(bbn.mesh.sf$vertices) <- 32619
knots.sf <- st_as_sf(as.data.frame(bbn.mesh$knots$centers), coords = c("X","Y")) 
knots.sf$geometry <- knots.sf$geometry*1000
st_crs(knots.sf) <- 32619

# Plot the mesh
#ggplot(bbn.mesh.sf$triangles) + geom_sf() + geom_sf(data= bbn.shape,fill = NA,color = 'darkblue',size=2) + geom_sf(data = knots.sf,fill = NA)
# Now make the prediction grid
pred.grid<-setup_pred_grid(knots=bbn.mesh$knots,model_bound=bbn.mesh$utm_bound)
# Plot the grid
#ggplot(pred.grid$grid) + geom_sf(aes(fill = as.factor(knotID))) + scale_fill_viridis_d()

# For the moment we need to have this starting at year 1.
catch.sf$Year <-  catch.sf$Year - (min(years)-1)
catchy <- catch_spread(catch = catch.sf,knots = bbn.mesh$knots)
catchy$sum_catches[,ncol(catchy$sum_catches)+1] <- 0
# Catch for TLM
catch.tlm <- catch.sf %>% group_by(Year,.drop=F) %>% dplyr::summarise(catch = sum(Catch,na.rm=T))
catch.tlm[nrow(catch.tlm)+1,] <- catch.tlm[nrow(catch.tlm),]
catch.tlm$Year[nrow(catch.tlm)] <- max(catch.tlm$Year) + 1
catch.tlm$catch[nrow(catch.tlm)] <- 0

# The SEBDAM set data.
if(mod.select == "SEAM")
{
  set_data<-data_setup(data=mod.input.sf,growths=data.frame(g = g$g,gR = g$gR),catch=as.data.frame(catchy$sum_catches),
                       model="SEBDAM",mesh=bbn.mesh$mesh,obs_mort=T,prior=T,prior_pars=c(20,40),
                       mult_qI=vary.q,spat_approach="spde",
                       knot_obj=bbn.mesh$knots,knot_area=pred.grid$area,separate_R_aniso = T,
                       all_se=T,weighted_mean_m = T)
  str(set_data)
  
  # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
  #set_data$par$log_m0 <- log(init.m)
  #set_data$par$log_R0 <- l.init.R 
  set_data$par$log_qR <- log(qR)
  set_data$par$log_S <- log(0.0695)
  # #set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA),log_qR = factor(NA))
  set_data$map <-list(log_qR = factor(NA),
                      log_S = factor(NA))
  #set_data$map <-list(log_qR = factor(NA))
  #set_data$map <-list(log_m0=factor(NA))
} # end if(mod.select != "TLM")



# A TLM version of the same...
if(mod.select == "TLM")
{
  set_data<-data_setup(data=as.data.frame(mod.input.sf),growths=data.frame(g = g$g,gR = g$gR),
                       catch=catch.tlm$catch, model="TLM",obs_mort=TRUE,prior=TRUE,all_se = F)
  # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
  #set_data<-fix_param(obj=set_data, pars = list(log_q_R=lqr))
  #set_data$par$log_q_R <- log(qR) # 
  #set_data$map <-list(log_q_R=factor(NA))
} # end if(mod.select == "TLM")
# Let's try setting some initial parameters 
# set_data$par$log_B0 <- 20 # If you get an NA error in the model below, set log_B0 to a super high number
#set_data$par$log_R0 <- 1 # Testing to see if I start this lower if it impacts our recruitment estimates.
# obj<-TMB::MakeADFun(data=set_data$data,parameters=set_data$par,random=set_data$random,map=set_data$map,DLL="SEBDAM",silent=F)
# Opt<-try(stats::nlminb(start=obj$par,obj=obj$fn,gr=obj$gr,control=control),T)
# Report <- obj$report()
#saveRDS(Report,"D:/Github/BBn_model/Results/Sab_woes.Rds")
#save.image("D:/Github/BBn_model/Results/BBn_input.Rdata") # To save time with troubleshooting/loading if needed. Just the years you are using so you don't get confused...

mod.fit<-fit_model(set_data,silent=F)

# Now save the results appropriately
if(mod.select == "SEAM")
{
  #m0 <- signif(exp(set_data$par$log_m0),digits=2)
  #r0 <- signif(exp(set_data$par$log_R0),digits=2)
  #qR <- signif(exp(set_data$par$log_qR),digits=2)
  # Model name
  scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR,"_",num.knots,"_knots")
  
  # And save the model
  saveRDS(mod.fit,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
  saveRDS(bbn.mesh,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
  saveRDS(pred.grid,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
  saveRDS(mod.input.sf,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_model_output_",scenario.select,"_model_input.Rds"))
  
}

if(mod.select == "TLM") 
{
  #qR <- signif(exp(set_data$par$log_q_R),digits=2)
  scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR,"_",g.mod)
  
  saveRDS(mod.fit,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
}

################################################### End the initial model runs ###########################################
################################################### End the initial model runs ###########################################
################################################### End the initial model runs ###########################################






##################### Now load the model and make the figures! ##############################################

atow<-800*2.4384/10^6 # area of standard tow in km2
num.knots <- 20 # 10, 15, and 20
#RO <- 5 # Going to test 100, 250, and 500.
init.m <- 0.2 # log(0.05) # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05
qR <- 0.33
# The different growth models.
#g.mod <- 'g_1'
#g.mod <- 'g_original'
#g.mod <- 'alt_g'
#g.mod <- 'proper_g'
#qR  <- "0_5" # This is just for TLM models, 0_5, 0_3, and 0_1
R.size <- "75"
FR.size <- "90"
years <- 1994:2022
NY <- length(years)
c_sys <- 32619
theme_set(theme_few(base_size = 22))
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
mod.select <- "SEAM"
################################################### End the initial model runs ###########################################
### Make the figures for the models


if(mod.select != "TLM") scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR,"_",num.knots,"_knots")
if(mod.select == "TLM") scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR,"_",g.mod)

#if(mod.select == "TLM")  scenario.select <- paste0(min(years),"_",max(years),"_qR_",exp(lqr),"_new_g")
# If we are going with the no extra station model this is our scenario.
#scenario.select <- "1994_2022_vary_m_m0_1_R0_150_10_knots_No_extra_stations"
mod.fit <- readRDS(paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
if(mod.select != "TLM") catchy <- mod.fit$obj$env$data$C*mod.fit$obj$env$data$area # Get this into tonnes from catch density.
if(mod.select == "TLM") catchy <- mod.fit$obj$env$data$C

#This is only needed for SEAM.
if(mod.select != "TLM")
{
  pred.grid <- readRDS(paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
  
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
  # The catch data
  # So catches from June 2021-May 2022 are called 2021 and removed from the 2021 survey biomass (this is different indexing from how we used to handle this for offshore)
  
  F.dat<-data.frame(B=as.vector(mod.fit$report$areaB[,-ncol(mod.fit$report$areaB)]/1000),
                    C = c(as.vector(as.matrix(catchy[,-c((ncol(catchy)-1),ncol(catchy))])),rep(NA,num.knots)), Year=matYear1, knotID=knots1)
  F.dat <- F.dat %>% dplyr::mutate(exploit = C/(B+C)) # Sticking with how offshore does this (C/(B+C)) C/B or some variant may be more realistic
  F.dat.plot<-left_join(knot.gis,F.dat,by=c("knotID"))
  
  F.dat.plot <- F.dat.plot %>% dplyr::filter(Year != max(years) )
  # To get a weighted m time series combing the B, R, and the m data above, no longer necessary :-)
  # bmr <- data.frame(B=as.vector(mod.fit$report$areaB/1000),
  #                   m=as.vector(mod.fit$report$m),
  #                   R=c(as.vector(mod.fit$report$areaR/1000),rep(NA,num.knots)),
  #                   Year=matYear, knotID=knots)
  # # Now get an m estimate for the whole area for B, R, and both, I think easiest will be to get totB and totR by year in here...
  # tBR <- bmr %>% group_by(Year) %>% dplyr::summarise(totB = sum(B,na.rm=T),
  #                                                    totR = sum(R,na.rm=T),
  #                                                    totSSB = sum(B,na.rm=T) + sum(R,na.rm=T))
  # 
  # # Now merge that with bmr...
  # bmr <- left_join(bmr,tBR,by="Year")
  # 
  # # Now getting weighted means by year should be easy...
  # nat.mat <- bmr %>% group_by(Year) %>% dplyr::summarise(m.FR = sum(B*m/totB),
  #                                                        m.R = sum(R*m/totR),
  #                                                        m.all = sum((R+B)*m/totSSB))
  # nat.mat$Raph.m <- mod.fit$report$mean_m
  # # Now make them long...
  # nat.mat.plt <- pivot_longer(nat.mat,!Year,names_to = "method",values_to = 'm')
  # nat.mat.plt$method <- factor(nat.mat.plt$method,levels = c("m.all","Raph.m","m.FR","m.R"),labels = c("Weighted","Unweighted","Fully Recruited","Recruits"))
  # 
  
  
  # If not already loaded, make a map
  #b.map <- pecjector(area= "BBn",c_sys = c_sys,add_layer = list(land = 'grey',eez = 'eez' , nafo = 'main',sfa = 'offshore',survey = c("offshore","outline")),txt.size=8,axes = "DM")
  
  # Final piece is that we can compare the biomass knot by knot to what would be expected based on B/R/m/g/C alone.  This is a bit messy but this does the trick...
  gs <- mod.fit$obj$env$data$gI
  gRs <- mod.fit$obj$env$data$gR
  tmp <- NULL
  Bs.tot <- mod.fit$report$B*mod.fit$obj$env$data$area/1000
  Rs.tot <- mod.fit$report$R*mod.fit$obj$env$data$area/1000
  for(i in 1:(NY-1))
  {
   
    ms <- mod.fit$report$m
    Bst <- (exp(-ms[,i+1]))*gs[i]*(Bs.tot[,i]-catchy[,i]) 
    Rst <- (exp(-ms[,i+1]))*gRs[i]*(Rs.tot[,i]) 
    B2 <- Bst + Rst
    if(any(B2 < 0)) {B2[B2<= 0] <- NA; print("HEADS UP!! You have a negative B2 estimate (look for the NAs in B.exp).")} # This feels very hacky for the moment, I think once Raph has the catch alignment fixed up this will not matter anymore, but not 100% sure!
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
  saveRDS(Bdiff.comp.dat.plot,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                                     mod.select,"_model_output_",scenario.select,"_B_differnce.Rds"))
  
  # Smaller text for spatial figures.
  theme_set(theme_few(base_size = 14))
  
  #Spatial predictions
  #B
  #B_plot <- st_transform(B_plot,crs = 4326)
  # Set up pretty breaks for the figure
  b.brk <- pretty(log(B.dat.plot$B))
  b.lab <- signif(exp(b.brk),digits=2)
  spatial.B.plot<- ggplot() + geom_sf(data=B.dat.plot%>% dplyr::filter(!Year %in% c(2023)),aes(fill=log(B)),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_viridis_c(breaks = b.brk, labels=b.lab, name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_biomass.png"),spatial.B.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.B.plot<- ggplot() + geom_sf(data=B.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=log(B)),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_viridis_c(breaks = b.brk, labels=b.lab, name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_biomass_no_missing_surveys.png"),spatial.B.plot,base_width = 10,base_height = 10)
  # Subest to 4 years
  spatial.B.plot<- ggplot() + geom_sf(data=B.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(B)),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_viridis_c(breaks = b.brk, labels=b.lab, name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_biomass_4_years.png"),spatial.B.plot,base_width = 10,base_height =7)
  
  #R
  r.brk <- pretty(log(R.dat.plot$R))
  r.lab <- signif(exp(r.brk),digits=2)
  
  spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot,aes(fill=log(R)),col='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_recruits.png"),spatial.R.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=log(R)),col='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_recruits_no_missing_surveys.png"),spatial.R.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(R)),col='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_recruits_4_years.png"),spatial.R.plot,base_width = 10,base_height = 7)
  
  #m
  m.brk <- log(c(0.003,0.007,0.02,0.05,0.15,0.4,1))
  #m.brk <- pretty(log(m.dat.plot$m))
  m.lab <- signif(exp(m.brk),digits=2)
  
  spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot %>% dplyr::filter(!Year %in% c(2023)),aes(fill=log(m)),color='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality (Inst)",option = "B",direction =1,begin = 0.2,end=1) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_mort.png"),spatial.m.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=log(m)),color='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality (Inst)",option = "B",direction =1,begin = 0.2,end=1) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_mort_no_missing_surveys.png"),spatial.m.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(m)),color='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality (Inst)",option = "B",direction =1,begin = 0.2,end=1) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_mort_4_years.png"),spatial.m.plot,base_width = 10,base_height = 7)
  
  # q 
  spatial.q.plot <- ggplot() + geom_sf(data=q.dat.plot,aes(fill=qI),col=NA)+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_viridis_c(name="Predicted catchability (qI)",option = "C",begin = 0.2,end =0.8) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_catchability.png"),spatial.q.plot,base_width = 10,base_height = 10)
  
  
  # OK, so lets try an make a map of the spatial exploitation rates, not sure if this is all correct yet.
  # Log rescalling doesn't work because we have 0's here.
  # OK, so lets try an make a map of the spatial exploitation rates, not sure if this is all correct yet.
  F.dat.plot$exp.na <- NA
  F.dat.plot$exp.na[F.dat.plot$exploit != 0] <- F.dat.plot$exploit[F.dat.plot$exploit != 0]
  
  #e.brk <- log(c(0.0015,0.01,0.08,0.4))
  e.brk <- pretty(log(F.dat.plot$exp.na))
  e.lab <- signif(exp(e.brk),digits=2)
  
  
  spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot,aes(fill=log(exp.na)),color='grey') +
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year) + 
    scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Spatial_exploit.png"),spatial.exploit.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=log(exp.na)),color='grey') +
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year) + 
    scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Spatial_exploit_no_missing_surveys.png"),spatial.exploit.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(exp.na)),color='grey') +
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    facet_wrap(~Year) + 
    scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Spatial_exploit_4_years.png"),spatial.exploit.plot,base_width = 10,base_height = 10)
  
  
  bd.brk <- pretty(Bdiff.comp.dat.plot$B.diff)
  bd.lab <- bd.brk#signif(exp(b.brk),digits=2)
  
  spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot,aes(fill=B.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_BDiff.png"),spatial.Bdiff.plot,base_width = 10,base_height = 10)
  #Remove missing survey years
  spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=B.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_BDiff_no_missing_surveys.png"),spatial.Bdiff.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=B.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_BDiff_4_years.png"),spatial.Bdiff.plot,base_width = 10,base_height = 7)
  
  # Same plot but the percentage miss by cell.
  pb.brk <- pretty(Bdiff.comp.dat.plot$B.per.diff)
  pb.lab <- pb.brk#signif(exp(b.brk),digits=2)
  spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot,aes(fill=B.per.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_per_BDiff.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=B.per.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_per_BDiff_no_missing_surveys.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=B.per.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60B018'W","60B06'W","59B054'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42B036'N","42B045'N","42B054'N")) +
    scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Spatial_per_BDiff_4_years.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 7)
  
  theme_set(theme_few(base_size = 22))
  
  
#   # And now plot them....
#   m.comp <- ggplot(nat.mat.plt,aes(x=Year,y=m,group = method,color=method)) + geom_line(linewidth=1.5)  + 
#     xlab("") + ylab("Natural Mortality (instantaneous)")+
#     scale_color_manual(values = c("blue","orange","grey","black"))  + 
#     theme(legend.title=element_blank()) + scale_x_continuous(breaks = seq(1980,2030,by=3))
#   save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/natural_mortality_comparisons.png"),m.comp,base_height = 6,base_width = 10)
#   # Remove missing survey years
#   m.comp <- ggplot(nat.mat.plt %>% dplyr::filter(Year < 2020)) + geom_line(aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  + 
#     #geom_line(data = nat.mat.plt %>% dplyr::filter(Year %in% 2016:2019), aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  + 
#     geom_line(data = nat.mat.plt %>% dplyr::filter(Year %in% 2021:2022), aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  +  
#     xlab("") + ylab("Natural Mortality (instantaneous)")+
#     scale_color_manual(values = c("blue","orange","grey","black"))  + 
#     theme(legend.title=element_blank()) + scale_x_continuous(breaks = seq(1980,2030,by=3))
#   save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/natural_mortality_comparisons_no_missing_surveys.png"),m.comp,base_height = 6,base_width = 10)
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

# The catch data
# So catches from June 2021-May 2022 are called 2021 and removed from the 2021 survey biomass (this is different indexing from how we used to handle this for offshore)
#
# SS model mu(t) <- C(t) / (B(t) + C(t)) because our model is B(t) <- B(t-1) - C(t) and C(2017) is June 2016-Aug 2017.
#  mu[2017] <- C[June 2016-Aug 2017]/(B[2017]+C[June 2016-Aug 2017]) 
# TLM and SEAM don't calculate mu, so we do it manually here, to be analogous...
# SEAM/TLM mu(t) <- C(t-1) / (B(t) + C(t-1)) because our model is B(t) <- B(t-1) - C(t-1) and C(2016) is now June 2016-Aug 2017.
# mu[2017] <- C[June 2016-Aug 2017]/(B[2017]+C[June 2016-Aug 2017]) 

if(mod.select != "TLM")
{
ann.exploit <- data.frame(year = years,B = exp(pred.proc$log_processes$log_B), Catch = c(colSums(catchy[,-c((ncol(catchy)-1),ncol(catchy))]),NA),
                          B.LCI = pred.proc$log_processes$totB.LCI, B.UCI = pred.proc$log_processes$totB.UCI)
}
if(mod.select == "TLM")
{
  ann.exploit <- data.frame(year = years,B = exp(pred.proc$log_processes$log_B), Catch = c(catchy[1:(length(catchy)-2)],NA),
                            B.LCI = pred.proc$log_processes$totB.LCI, B.UCI = pred.proc$log_processes$totB.UCI)
}

ann.exploit$exploit <- c(ann.exploit$Catch[1:(nrow(ann.exploit)-1)]/(ann.exploit$B[2:nrow(ann.exploit)]+ann.exploit$Catch[1:(nrow(ann.exploit)-1)]),NA)
ann.exploit$exploit.UCI <- c(ann.exploit$Catch[1:(nrow(ann.exploit)-1)]/(ann.exploit$B.UCI[2:nrow(ann.exploit)]+ann.exploit$Catch[1:(nrow(ann.exploit)-1)]),NA)
ann.exploit$exploit.LCI <- c(ann.exploit$Catch[1:(nrow(ann.exploit)-1)]/(ann.exploit$B.LCI[2:nrow(ann.exploit)]+ann.exploit$Catch[1:(nrow(ann.exploit)-1)]),NA)
ann.exploit$FM <- 1-exp(-ann.exploit$exploit)
ann.exploit$FM.LCI <- 1-exp(-ann.exploit$exploit.LCI)
ann.exploit$FM.UCI <- 1-exp(-ann.exploit$exploit.UCI)


# Biomass time series
bm.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') +
  xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,3.1e4))
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Biomass_time_series.png"),bm.ts.plot,base_width = 11,base_height = 8.5)
# Recruit time series
rec.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Recruit Biomass (tonnes)")  + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,7e3))
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Recruit_time_series.png"),rec.ts.plot,base_width = 11,base_height = 8.5)
# Natural mortality time series...
mort.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Natural mortality (Instantaneous)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.35))
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_nat_mort_time_series.png"),mort.ts.plot,base_width = 11,base_height = 8.5)
# Explotation Rate Time Series
exploit.plot <- ggplot(ann.exploit) + geom_line(aes(x=year,y=exploit),size=1.5) + geom_ribbon(aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') +
  xlab("") + ylab("Exploitation Rate (Proportional)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.3))
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_exploit_time_series.png"),exploit.plot,base_width = 11,base_height = 8.5)


# Same plots but remvoing missing survey years....
bm.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2020))) + 
                            geom_line(aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            # geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  
                            # geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  
                            geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            ylim(c(0,3.1e4)) + xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Biomass_time_series_no_missing_surveys.png"),bm.ts.plot,base_width = 11,base_height = 8.5)
# Recruit time seris
rec.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2020))) + 
                            geom_line(aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            # geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  
                            # geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  
                            geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            xlab("") + ylab("Recruit Biomass (tonnes)") + ylim(c(0,7e3)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_Recruit_time_series_no_missing_surveys.png"),rec.ts.plot,base_width = 11,base_height = 8.5)
# Natural mortality time series...
mort.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2020))) + 
                            geom_line(aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            #geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  
                            #geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  
                            geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
                            xlab("") + ylab("Natural mortality (Instantaneous)") + ylim(c(0,0.35)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_nat_mort_time_series_no_missing_surveys.png"),mort.ts.plot,base_width = 11,base_height = 8.5)


# Explotation Rate Time Series
exploit.plot <- ggplot(ann.exploit%>% dplyr::filter(year < c(2020))) + geom_line(aes(x=year,y=exploit),linewidth = 1.5) +
                            geom_ribbon(aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') +
                            geom_line(data= ann.exploit%>% dplyr::filter(year %in% 2021:2022),aes(x=year,y=exploit),linewidth = 1.5) +
                            geom_ribbon(data= ann.exploit%>% dplyr::filter(year %in% 2021:2022),aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.5,fill='darkblue',color='darkblue') +
                            xlab("") + ylab("Exploitation Rate (Proportional)") + ylim(c(0,0.3)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/BBn_exploit_time_series_no_missing_surveys.png"),exploit.plot,base_width = 11,base_height = 8.5)


# Save out a bunch of objects...
# The biomass difference
saveRDS(Bdiff.comp.dat.plot,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                                   mod.select,"_model_output_",scenario.select,"_B_differnce.Rds"))
# Annual exploitation
saveRDS(ann.exploit,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                           mod.select,"_model_output_",scenario.select,"_annual.exploit.Rds"))
# The nicely summarized model object
saveRDS(pred.proc,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                         mod.select,"_model_output_",scenario.select,"_pred_proc.Rds"))
# The summarized spatial bits.
saveRDS(F.dat.plot,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                          mod.select,"_model_output_",scenario.select,"_F_spatial.Rds"))
# q spatial
saveRDS(q.dat.plot,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                          mod.select,"_model_output_",scenario.select,"_q_spatial.Rds"))
# m spatial
saveRDS(m.dat.plot,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                          mod.select,"_model_output_",scenario.select,"_m_spatial.Rds"))
# R spatial
saveRDS(R.dat.plot,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                          mod.select,"_model_output_",scenario.select,"_R_spatial.Rds"))
# B spatial
saveRDS(B.dat.plot,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",
                          mod.select,"_model_output_",scenario.select,"_B_spatial.Rds"))


# We don't think OSR residuals are appropriate for a spatial model, they might make sense if you could do them 
# knot by knot, but I haven't a clue how to do that....
## One Step ahead residuals, this seems to work pretty well, but is slow so only
# calculate these for my favourite models. Also figure out what they mean!!
# osr <- oneStepPredict(mod.fit$obj,observation.name="logI", data.term.indicator = 'keep_I', method="oneStepGaussianOffMode",discrete = F)
# saveRDS(osr,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_",scenario.select,"_One_Step_residuals.Rds"))
# 
# osr <- readRDS(paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_",scenario.select,"_One_Step_residuals.Rds"))
# # So this takes about 3 seconds for every observation
# plot(osr$residual)
# # I think we can combine these with the data and look at them by year and probably spatially, go explore that tomorrow....
# osr.dat <- data.frame(osr,I = mod.fit$obj$env$data$logI,
#                       IR = mod.fit$obj$env$data$logIR)
# 
# ggplot(osr.dat) + geom_point(aes(x=I,y=residual))
# ggplot(osr.dat) + geom_point(aes(x=IR,y=residual))
# 
# 
# 
# 
# qq.plt <- gg_qq(osr$residual)
# 
# gg_qq <- function(x, distribution = "norm", ..., line.estimate = NULL, conf = 0.95,
#                   labels = names(x)){
#   q.function <- eval(parse(text = paste0("q", distribution)))
#   d.function <- eval(parse(text = paste0("d", distribution)))
#   x <- na.omit(x)
#   ord <- order(x)
#   n <- length(x)
#   P <- ppoints(length(x))
#   df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
#   
#   if(is.null(line.estimate)){
#     Q.x <- quantile(df$ord.x, c(0.25, 0.75))
#     Q.z <- q.function(c(0.25, 0.75), ...)
#     b <- diff(Q.x)/diff(Q.z)
#     coef <- c(Q.x[1] - b * Q.z[1], b)
#   } else {
#     coef <- coef(line.estimate(ord.x ~ z))
#   }
#   
#   zz <- qnorm(1 - (1 - conf)/2)
#   SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
#   fit.value <- coef[1] + coef[2] * df$z
#   df$upper <- fit.value + zz * SE
#   df$lower <- fit.value - zz * SE
#   
#   if(!is.null(labels)){ 
#     df$label <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
#   }
#   
#   p <- ggplot(df, aes(x=z, y=ord.x)) +
#     geom_point() + 32
#     geom_abline(intercept = coef[1], slope = coef[2]) +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) 
#   if(!is.null(labels)) p <- p + geom_text( aes(label = label))
#   print(p)
#   coef
# }
