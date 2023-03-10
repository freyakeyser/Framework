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


source("D:/Github/Assessment_fns/Fishery/logs_and_fishery_data.r")
source("D:/Github/Assessment_fns/Maps/pectinid_projector_sf.R")
source("D:/Github/Assessment_fns/Maps/convert_inla_mesh_to_sf.R")

########################################################################################################
# Bring in the data and tidy it up for the analysis

bbn.shape <- st_read("D:/Github/GIS_layers/survey_boundaries/BBn.shp", quiet=T)
bbn.shape <- bbn.shape %>% st_transform(crs = 32619) # BBn is right on the 19/20 border so think they are basically equivalent options here
pred.grid <- readRDS("D:/Github/BBn_model/Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_150_10_knots_predict_grid.Rds")
# Bring in the survey data
#load("Y:/Offshore/Assessment/Data/Survey_data/2019/Survey_summary_output/Survey_all_results.Rdata")
#load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/Survey_all_results.Rdata")
#surv.dat <- surv.dat$BBn
#saveRDS(surv.dat,'D:/Github/BBn_model/Results/BBn_surv.dat.RDS')
surv.dat <- readRDS('D:/Github/BBn_model/Results/BBn_surv.dat.RDS')
# This gets our growth data, grabbing it from 2022 because we have the input data and so is way smaller
# load("Y:/Offshore/Assessment/Data/Model/2022/BBn/Model_input_mixed.RData")
# mod.dat <- mod.dat$BBn
# saveRDS(mod.dat,'D:/Github/BBn_model/Results/BBn_model.dat.RDS')
mod.dat <- readRDS('D:/Github/BBn_model/Results/BBn_model.dat.RDS')
# For now we need to get a 2022 growth term (once we run the BBn model I can update this with the latest data), just recycling the growth for 2021 for the moment
mod.dat <- rbind(mod.dat,mod.dat[nrow(mod.dat),])
mod.dat$year[nrow(mod.dat)] <- 2022
# Bring in the fishery data
# logs_and_fish(loc="offshore",year = 1991:2022,un=un.ID,pw=pwd.ID,db.con=db.con,direct="Y:/Offshore/Assessment/", get.marfis=F)
# fish.dat<-merge(new.log.dat,old.log.dat,all=T)
# fish.dat$ID<-1:nrow(fish.dat)
# # Now subset to BBn and add in the missing years of data
# bbn.fish <- fish.dat %>% dplyr::filter(bank == "BBn")
# # There are 12 data points at 0,0 that we remove, I'm not worried about accounting for these 12 points!
# bbn.fish <- bbn.fish %>% dplyr::filter(lat !=0 | lon != 0)
# # Now I want to put a 'survey year' on these because that's what we're gonna need for our modelling... start by porting over the year
# bbn.fish$survey.year <- bbn.fish$year
# # Need to add fake data for 2009
# bbn.fish[nrow(bbn.fish)+1,] <- NA
# 
# bbn.fish$pro.repwt[nrow(bbn.fish)] <- 0
# bbn.fish$year[nrow(bbn.fish)] <- 2009
# # See DK NOte below
# bbn.fish$survey.year[nrow(bbn.fish)] <- 2009
# bbn.fish$month[nrow(bbn.fish)] <- "June"
# # Getting fake lat/lon coords that are on BBn
# bbn.fish$lat[nrow(bbn.fish)] <- 42.85600
# bbn.fish$lon[nrow(bbn.fish)]  <- -65.90183
# # saveRDS(bbn.fish,'D:/Github/BBn_model/Results/BBn_fish.dat.RDS')

bbn.fish <- readRDS('D:/Github/BBn_model/Results/Fishery_data/BBn_fish.dat.RDS')
bbn.fish$pro.repwt <- bbn.fish$pro.repwt/1000
#### Finshed Data prep and clean up!
###############################################################################################


# Set parameters for the run...
repo.loc <- "D:/Github/BBn_model/"
mod.select <- "SEAM"
atow<-800*2.4384/10^6 # area of standard tow in km2
num.knots <- 10
years <- 1994:2022
NY <- length(years)


# OK, so step 1 here is getting the model input that Raphael needs for the model
# The survey data....
live.subset <- surv.dat %>% dplyr::filter(state == 'live')
dead.subset <- surv.dat %>% dplyr::filter(state== "dead")
live.input <- data.frame(I = live.subset$com.bm, IR = live.subset$rec.bm,year = live.subset$year,tow = live.subset$tow,tot.live.com = live.subset$com,lat = live.subset$lat,lon=live.subset$lon,tow_type = live.subset$random)
clap.input <- data.frame(L = dead.subset$com,tow = dead.subset$tow,year = dead.subset$year,tow_type = dead.subset$random)
mod.input <- left_join(live.input,clap.input,by=c('tow','year','tow_type'))
# Looks like there are no values > 0 but < 0.5, so the low clapper numbers should all round up to 1 (which makes sense as you'd only get < 0.5 if we had tows twice as long as they should be)
mod.input.sf <- st_as_sf(mod.input,coords = c('lon','lat'),remove=F,crs = 4326)
mod.input.sf <- mod.input.sf %>% st_transform(crs=32619)
mod.input.sf <- st_intersection(mod.input.sf,bbn.shape)
mod.input.sf$Year <- mod.input.sf$year - (min(years)-1)
mod.input.sf$I <- mod.input.sf$I/atow
mod.input.sf$IR <- mod.input.sf$IR/atow
mod.input.sf <- mod.input.sf %>% dplyr::filter(year %in% years)
#Sidebar exploration comparing extras and regular survey tows within a knot
grid.gis <- pred.grid$grid
grid.gis$geometry <- grid.gis$geometry*1000
st_crs(grid.gis) <- 32619
# Now simplify the grid for the spatial plots...
knot.gis <- aggregate(grid.gis, list(grid.gis$knotID), function(x) x[1])
ggplot(knot.gis) + geom_sf() + geom_sf_text(aes(label = knotID))
# Now intersect the model data with the knots
mod.input.knots <- st_intersection(mod.input.sf,knot.gis)

extras.by.knot.year <- mod.input.knots %>% dplyr::filter(tow_type != 1) %>% dplyr::group_by(year,knotID) %>% dplyr::summarise(mn.FR = mean(I,na.rm=T),
                                                                                                                              mn.R = mean(IR,na.rm=T),
                                                                                                                              sd.FR = sd(I,na.rm=T),
                                                                                                                              sd.R = sd(IR,na.rm=T))

reg.by.knot.year <- mod.input.knots %>% dplyr::filter(tow_type == 1) %>% dplyr::group_by(year,knotID) %>% dplyr::summarise(mn.FR.reg = mean(I,na.rm=T),
                                                                                                                           mn.R.reg = mean(IR,na.rm=T),
                                                                                                                           sd.FR.reg = sd(I,na.rm=T),
                                                                                                                           sd.R.reg = sd(IR,na.rm=T))

extras_vs_regs <- left_join(data.frame(extras.by.knot.year),data.frame(reg.by.knot.year),by = c('year',"knotID"))

extras_vs_regs$FR.diff <- extras_vs_regs$mn.FR - extras_vs_regs$mn.FR.reg
extras_vs_regs$R.diff <- extras_vs_regs$mn.R - extras_vs_regs$mn.R.reg
# Make 0's very small....
#extras_vs_regs$mn.FR.reg[is.na(extras_vs_regs$mn.FR.reg)] <- 0.1*min(extras_vs_regs$mn.FR.reg,na.rm=T)

extras_vs_regs$mn.R.reg[extras_vs_regs$mn.R.reg == 0] <- NA


extras_vs_regs$per.FR.diff <- 100*((extras_vs_regs$mn.FR - extras_vs_regs$mn.FR.reg) / extras_vs_regs$mn.FR.reg)
extras_vs_regs$per.R.diff <- 100* ((extras_vs_regs$mn.R - extras_vs_regs$mn.R.reg) /extras_vs_regs$mn.R.reg)

windows(11,11)
ggplot(extras_vs_regs) + geom_boxplot(aes(y = per.FR.diff))


ggplot(extras_vs_regs) + geom_boxplot(aes(y = per.R.diff)) + ylim(c(-100,1000))

ggplot(extras_vs_regs) + geom_text(aes(x = year, y = per.FR.diff, label= knotID)) +
  geom_hline(yintercept = c(-50,0,100),linetype = 'dashed',color='blue') +
  ylab("% Difference between Extras and Regular Tows")

ggplot(extras_vs_regs) + geom_text(aes(x = year, y = FR.diff, label= knotID)) +
  #geom_hline(yintercept = c(-50,0,100),linetype = 'dashed',color='blue') +
  ylab("Difference between Extras and Regular Tows (kg km^-2)")

ggplot(extras_vs_regs) + geom_text(aes(x = year, y = per.R.diff, label= knotID)) + # ylim(c(-100,1000)) +
  geom_hline(yintercept = c(-50,0,100),linetype = 'dashed',color='blue') +
  ylab("% Difference between Extras and Regular Tows")



summary(extras_vs_regs$per.FR.diff)
summary(extras_vs_regs$per.R.diff)

ggplot(knot.gis) + geom_sf() + geom_sf_text(aes(label = knotID))
# Here I do a special with no extra stations included...
live.subset <- surv.dat %>% dplyr::filter(state == 'live' & random == 1)
dead.subset <- surv.dat %>% dplyr::filter(state== "dead"  & random == 1)
# 