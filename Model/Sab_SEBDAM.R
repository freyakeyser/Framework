# So here we'll try and get all the data in the correct structure for the spatial model.
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
load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90.Rdata")
#load("D:/Framework/SFA_25_26_2024/Model/Data/testing_results_framework3.Rdata")# Need to get condition factor out of here too
#mod.dat <- survey.obj$Sab$model.dat
#saveRDS(mod.dat,'D:/Github/BBn_model/Results/Sab_model.dat.RDS')
surv.dat <- surv.dat$Sab
mod.dat <- survey.obj$Sab$model.dat

#load("F:/NAS/Offshore/Assessment/Data/Model/2022/Sab/Model_input_midpoint.RData")
# Bring in the fishery data
# logs_and_fish(loc="offshore",year = 1986:2022,direct="Y:/Offshore/Assessment/", get.marfis=F)
# fish.dat<-merge(new.log.dat,old.log.dat,all=T)
# fish.dat$ID<-1:nrow(fish.dat)

# Now we can clip both of these to subset it to the data that I think we need for the analysis....
# First the fishery data

# Sab.fish <- fish.dat %>% dplyr::filter(bank == "Sab")
# # There are 12 data points at 0,0 that we remove, I'm not worried about accounting for these 12 points!
# Sab.fish <- Sab.fish %>% dplyr::filter(lat !=0 | lon != 0)
# # Now I want to put a 'survey year' on these because that's what we're gonna need for our modelling... start by porting over the year
# Sab.fish$survey.year <- Sab.fish$year
# # DK NOTE: Now this is going to get confusing for us and we may want to tweak SEBDAM for this, but that's a down the road job, not a playing around with model job
# # But based on the indexing in SEBDAM, I am going to change how we index the survey year data from what we have done with offshore traditionally.
# # Historically anything from the last half of the year goes into the following years, eg. survey.year 2002 = June 2001- May 2002.
# # But in SEBDAM we have (B(t-1) - C(t-1)), so let's say we have year 2000 survey biomass, this says we remove the 2000 catch from that
# # we want that catch to be the catch from June 2000 to May 2001, i.e. we remove the catch before we allow the population to grow
# # This is what we do in our current model, but we have a different index (C(t) on our model.
# # Basically survey year 2002 = June 2002 - May 2003 now
# #DK note: We probably should think more about the survey year fun and how exactly we want to handle removal of catch in our models.
# Sab.fish$month <- lubridate::month(Sab.fish$date)
# Sab.fish$survey.year[Sab.fish$month %in% c("January","February","March","April","May")] <- Sab.fish$survey.year[Sab.fish$month %in% c("January","February","March","April","May")] -1
# # Add a fake 2022 data point as there were no removals in 2022
# Sab.fish[nrow(Sab.fish)+1,] <- NA
# Sab.fish$pro.repwt[nrow(Sab.fish)] <- 0
# Sab.fish$year[nrow(Sab.fish)] <- 2022
# Sab.fish$survey.year[nrow(Sab.fish)] <-2022
# Sab.fish$lon[nrow(Sab.fish)] <- -61.68767
# Sab.fish$lat[nrow(Sab.fish)] <- 43.63017
#saveRDS(Sab.fish,'D:/Github/BBn_model/Results/Fishery_data/Sab_fish.dat.RDS')
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
lqr <- log(0.75)# This is for TMB (log recruit catchability) testing catchability of 0.5, test 0.3 and 0.1
l.init.m <- log(2) # This is for SEAM, sets first year natural mortality, going to test 2, 0.8, and 0.4
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
mod.growth.dat <- mod.dat
# These parameters come from our von B data that was aged by Trish in 2023.
L.inf <- 159.2
#to <- 1.337 # So this uses a 1 year offset that we no longer believe in, going to make this 0.337 to align more with what we now do...
to <- 0.2
K <- 0.2


# This is weight in this year, which becomes t-1 
waa.tm1 <- mod.growth.dat$CF*(mod.growth.dat$l.bar/100)^3 # Average(ish) weight of commerical sized scallop in the current year
# Using this years average shell height we can find the exptected shell height for the scallops in the next year
# ht = (Linf * (1-exp(-K)) + exp(-K) * height(last year))
# laa.t is the projected size of the current years scallops into next year.
laa.t <- L.inf*(1-exp(-K)) + exp(-K) * mod.growth.dat$l.bar
# The c() term in the below offsets the condition so that current year's condition slots into the previous year and repeats 
# the condition for the final year), this effectively lines up "next year's condition" with "predictied shell height next year (laa.t)
# This gets us the predicted weight of the current crop of scallops next year based on next years CF * laa.t^3
# Of course we don't have next years condition thus th last condition is simply repeated
# waa.t is using the condition from next year and the growth from next year to get next years weight
waa.t <- c(mod.growth.dat$CF[-1],mod.growth.dat$CF[nrow(mod.growth.dat)])*(laa.t/100)^3
# Here we use the current condition factor to calculate the weight next year (since we use laa.t)
# That's really the only difference between waa.t and waa.t2, waa.t uses next years condition to project growth
# what waa.t2 uses the current condition to project growth.  So that's really what we are comparing here with these
# two growth metrics isn't it, this is really just comparing impact of using current vs. future condition factor on our growth estimates.

waa.t2 <-mod.growth.dat$CF*(laa.t/100)^3
# Now the growth, expected and realized.
mod.growth.dat$g <- waa.t/waa.tm1
# This is using the actual condition factor and growing the scallops by laa.t
mod.growth.dat$g2 <- waa.t2/waa.tm1

# same thing here but for the recruits
waa.tm1 <-mod.growth.dat$CF*(mod.growth.dat$l.k/100)^3
laa.t <- L.inf*(1-exp(-K))+exp(-K)*mod.growth.dat$l.k
waa.t <- c(mod.growth.dat$CF[-1],mod.growth.dat$CF[nrow(mod.growth.dat)])*(laa.t/100)^3
waa.t2 <- mod.growth.dat$CF*(laa.t/100)^3
mod.growth.dat$gR <- waa.t/waa.tm1
mod.growth.dat$gR2 <- waa.t2/waa.tm1# setwd("C:/Assessment/2014/r")

# Need to replace 2019 and 2020 values because of missing 2020 survey. 2021 is NOT influenced by missing data in 2020, it's only 2019 and 2020 that need imputed
# how we have everything set up. I am going to do this the same way we have set it up for BBn and GB for consistency sake.



# Now we fill in 2014, 2015, 2019 and 2020 because of missing surveys in 2015 and 2020
mod.growth.dat[nrow(mod.growth.dat)+ 1,] <- NA
mod.growth.dat$year[nrow(mod.growth.dat)] <- 2020
mod.growth.dat[nrow(mod.growth.dat)+ 1,] <- NA
mod.growth.dat$year[nrow(mod.growth.dat)] <- 2015

mod.growth.dat <- mod.growth.dat[order(mod.growth.dat$year),]
# Fill in 2019 and 2020 with LTM growth before 2020.
mod.growth.dat$g[which(mod.growth.dat$year %in% c(2014:2015,2019:2020))] <- median(mod.growth.dat$g[mod.growth.dat$year<2020], na.rm=T)
mod.growth.dat$g2[which(mod.growth.dat$year %in% c(2015,2020))] <- median(mod.growth.dat$g2[mod.growth.dat$year<2020], na.rm=T)
mod.growth.dat$gR[which(mod.growth.dat$year %in% c(2014:2015,2019:2020))] <- median(mod.growth.dat$gR[mod.growth.dat$year<2020], na.rm=T)
mod.growth.dat$gR2[which(mod.growth.dat$year %in% c(2015,2020))] <- median(mod.growth.dat$gR2[mod.growth.dat$year<2020], na.rm=T)

# Turn this into a vector and add a value for next year.
growth <- data.frame(g = c(mod.growth.dat$g,mod.growth.dat$g2[nrow(mod.growth.dat)]), 
                     gR = c(mod.growth.dat$gR,mod.growth.dat$gR2[nrow(mod.growth.dat)]),
                     year =  c(mod.growth.dat$year,max(mod.growth.dat$year+1)))

growth <- growth %>% dplyr::filter(year %in% c(years,(max(years)+1)))

 # # We need to make up some data for 2020 since there wasn't a survey.  We could also do it this way, result is very simlar
# # First we take them mean CF between 2019 and 2021
# mod.growth.dat$CF[mod.growth.dat$year == 2020] <- median(mod.growth.dat$CF[mod.growth.dat$year %in% 2019:2021],na.rm=T)
# mod.growth.dat$l.bar[mod.growth.dat$year == 2020] <- median(mod.growth.dat$l.bar[mod.growth.dat$year %in% 2019:2021],na.rm=T)
# mod.growth.dat$l.k[mod.growth.dat$year == 2020] <- median(mod.growth.dat$l.k[mod.growth.dat$year %in% 2019:2021],na.rm=T)
#von.B <- function(L.inf,to,K) {L <- L.inf*(1-exp(-K*(age-t0)))}
# # GB
# k = 0.22
# Linf <- 149
# t0 <- 0.22
# age=3:6
# 
# Linf*(1-exp(-K*(age-t0)))
# 
# #BBN
# k = 0.19
# Linf <- 148
# t0 <- 0.11
# age=3:6
# 
# Linf*(1-exp(-K*(age-t0)))
# Here are the Sable parameters, note that I believe the to is due to differences in ageing that was done in the 
# 1980s-early 2000s which we have moved away from since 2010, this aligns us with NOAA methods
# #Sab
# k = 0.2269
# Linf <- 136
# t0 <- 1.337
# age=3:6
# 
# Linf*(1-exp(-K*(age-t0)))

#
#growth <- data.frame(g = runif(NY+1,1.1,1.2),gR = runif(NY+1,1.2,1.4))

# Subset the fishery data to the correct years. We need to take on next year catch too..)
Sab.fish <- Sab.fish %>% dplyr::filter(survey.year %in% years) # c(years,(max(years)+1))
#tst <- Sab.fish.sf %>% dplyr::filter(survey.year %in% 1993:2014)
#rems <- tst %>% dplyr::group_by(year) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))
Sab.fish.sf <- st_as_sf(Sab.fish,coords = c("lon","lat"),remove =F, crs = 4326)
Sab.fish.sf <- Sab.fish.sf %>% st_transform(crs= c_sys)

# Now lets clip this to be data inside of our Sab boundary.
Sab.fish.sf <- st_intersection(Sab.fish.sf,Sab.shape)

# Check removals each fishing year calculated using this data
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

Sab.mesh <- setup_mesh(mod.input.sf,model_bound = Sab.shape,nknot=num.knots, max.edge = c(8,20),cutoff=2.5,seed=34) 
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
set_data<-data_setup(data=mod.input.sf,growths=growth[,1:2],catch=as.data.frame(catchy$sum_catches),
                     model="SEBDAM",mesh=Sab.mesh$mesh,obs_mort=T, prior=TRUE,prior_pars=c(10,12),
                     mult_qI=T,spat_approach="spde",
                     knot_obj=Sab.mesh$knots,knot_area=pred.grid$area,separate_R_aniso =T,
                     all_se=T,weighted_mean_m = T)
# So this will fix the mean value of m0 to be whatever the intial value is set at.  Let's see what happens!
set_data$par$log_m0 <- l.init.m # 0 = 1, 
#set_data$par$log_R0 <- l.init.R # 5.3 = 200, 5 = 148, 4 = 55, 5.9915 = 400, 4.606 = 100
set_data$par$log_qR <- lqr
#set_data$map <-list(log_m0=factor(NA))
set_data$map <-list(log_m0=factor(NA),log_qR = factor(NA))
#set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA),log_qR = factor(NA))
}
#TLM version, Dude is this ever sensitive to the q priors! (5,12) actually looks solid in terms of results... maybe we can get so lucky with SEBDAM :-)
# Note that the Catch time series should be 1 year longer than the survey data here!!
if(mod.select == "TLM")
{
set_data<-data_setup(data=as.data.frame(mod.input.sf),growths=growth[,1:2],catch=catch.tlm$catch,
                      model="TLM",obs_mort=T,prior=T,prior_pars=c(10,12))
set_data$par$log_q_R <- lqr 
set_data$map <-list(log_q_R=factor(NA))
}
# Run the model
mod.fit<-fit_model(set_data,silent=F)

# Now save the results appropriately
if(mod.select != "TLM") 
{
  m0.par <- signif(exp(set_data$par$log_m0),digits=2)
  qR.par <- signif(exp(set_data$par$log_qR),digits=2)
  
  scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",m0.par,"_qR_",qR.par,"_",num.knots,"_knots")
  saveRDS(mod.fit,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
  saveRDS(Sab.mesh,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
  saveRDS(pred.grid,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
}

if(mod.select == "TLM") 
{
  qR.par <- sub("0.","0_",signif(exp(set_data$par$log_q_R),digits=2))
  scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR.par)
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
m0 <- 2 # log(0.05) # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05
qR  <- "0.75" 
years <- 1994:2022
NY <- length(years)
c_sys <- 32620
theme_set(theme_few(base_size = 22))
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
mod.select <- "SEAM"
################################################### End the initial model runs ###########################################
### Make the figures for the models



if(mod.select != "TLM") scenario.select <- scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",m0,"_qR_",qR,"_",num.knots,"_knots")
if(mod.select == "TLM") scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR)
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
                  C = c(rep(NA,num.knots),as.vector(as.matrix(catchy[,-c((ncol(catchy)-1),ncol(catchy))]))), Year=matYear1, knotID=knots1)
F.dat <- F.dat %>% dplyr::mutate(exploit = C/(B+C)) # Sticking with how offshore does this (C/(B+C)) C/B of some variant may be more realistic
F.dat.plot<-left_join(knot.gis,F.dat,by=c("knotID"))
# Remove the first year of data because it's not needed now.
F.dat.plot <- F.dat.plot %>% dplyr::filter(Year != years[1] )
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
saveRDS(Bdiff.comp.dat.plot,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",
                                   mod.select,"_model_output_",scenario.select,"_B_differnce.Rds"))

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
  ann.exploit <- data.frame(year = years,B = exp(pred.proc$log_processes$log_B), Catch = c(NA,colSums(catchy[,-c((ncol(catchy)-1),ncol(catchy))])),
                            B.LCI = pred.proc$log_processes$totB.LCI, B.UCI = pred.proc$log_processes$totB.UCI)
}
if(mod.select == "TLM")
{
  ann.exploit <- data.frame(year = years,B = exp(pred.proc$log_processes$log_B), Catch = c(NA,catchy[1:(length(catchy)-2)]),
                            B.LCI = pred.proc$log_processes$totB.LCI, B.UCI = pred.proc$log_processes$totB.UCI)
}
ann.exploit$exploit <- ann.exploit$Catch/(ann.exploit$B+ann.exploit$Catch)


#pred.proc$log_processes <- pred.proc$log_processes %>% dplyr::filter(year < 2023)


# Biomass time series
bm.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + ylim(c(0,1.6e4)) + 
  xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Biomass_time_series.png"),bm.ts.plot,base_width = 11,base_height = 8.5)
# Recruit time seris
rec.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Recruit Biomass (tonnes)") + ylim(c(0,3.5e3)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Recruit_time_series.png"),rec.ts.plot,base_width = 11,base_height = 8.5)
# Natural mortality time series...
mort.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Natural mortality (Instantaneous)") + ylim(c(0,2.5)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_nat_mort_time_series.png"),mort.ts.plot,base_width = 11,base_height = 8.5)
# Explotation Rate Time Series
exploit.plot <- ggplot(ann.exploit) + geom_line(aes(x=year,y=exploit),linewidth = 1.5) +
  xlab("") + ylab("Exploitation Rate (Proportional)") + ylim(c(0,0.2)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_exploit_time_series.png"),exploit.plot,base_width = 11,base_height = 8.5)

# Same plots but remvoing missing survey years....
bm.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2015))) + 
  geom_line(aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  ylim(c(0,1.6e4)) + xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Biomass_time_series_no_missing_surveys.png"),bm.ts.plot,base_width = 11,base_height = 8.5)
# Recruit time seris
rec.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2015))) + 
  geom_line(aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Recruit Biomass (tonnes)") + ylim(c(0,3.5e3)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_Recruit_time_series_no_missing_surveys.png"),rec.ts.plot,base_width = 11,base_height = 8.5)
# Natural mortality time series...
mort.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2015))) + 
  geom_line(aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  
  geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.2,fill='darkblue',color='darkblue') + 
  xlab("") + ylab("Natural mortality (Instantaneous)") + ylim(c(0,1.5)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_nat_mort_time_series_no_missing_surveys.png"),mort.ts.plot,base_width = 11,base_height = 8.5)


# Explotation Rate Time Series
exploit.plot <- ggplot(ann.exploit%>% dplyr::filter(year < c(2015))) + geom_line(aes(x=year,y=exploit),linewidth = 1.5) +
  geom_line(data = ann.exploit%>% dplyr::filter(year %in% 2016:2019), aes(x=year,y=exploit),linewidth = 1.5) + 
  geom_line(data= ann.exploit%>% dplyr::filter(year %in% 2021:2022),aes(x=year,y=exploit),linewidth = 1.5) +
  xlab("") + ylab("Exploitation Rate (Proportional)") + ylim(c(0,0.2)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/",mod.select,"_",scenario.select,"/Sab_exploit_time_series_no_missing_surveys.png"),exploit.plot,base_width = 11,base_height = 8.5)
