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
library(shape)


source("D:/Github/Assessment_fns/Fishery/logs_and_fishery_data.r")
source("D:/Github/Assessment_fns/Maps/pectinid_projector_sf.R")
source("D:/Github/Assessment_fns/Maps/convert_inla_mesh_to_sf.R")

########################################################################################################
# Bring in the data and tidy it up for the analysis

bbn.shape <- st_read("D:/Github/GIS_layers/survey_boundaries/BBn.shp", quiet=T)
bbn.shape <- bbn.shape %>% st_transform(crs = 32619) # BBn is right on the 19/20 border so think they are basically equivalent options here
# Bring in the survey data

surv.dat <- readRDS('D:/Github/BBn_model/Results/BBn_surv.dat.RDS')
# This gets our growth data, grabbing it from 2022 because we have the input data and so is way smaller
mod.dat <- readRDS('D:/Github/BBn_model/Results/BBn_model.dat.RDS')
# For now we need to get a 2022 growth term (once we run the BBn model I can update this with the latest data), just recycling the growth for 2021 for the moment
mod.dat <- rbind(mod.dat,mod.dat[nrow(mod.dat),])
mod.dat$year[nrow(mod.dat)] <- 2022
# Bring in the fishery data
bbn.fish <- readRDS('D:/Github/BBn_model/Results/Fishery_data/BBn_fish.dat.RDS')
bbn.fish$pro.repwt <- bbn.fish$pro.repwt/1000
#### Finshed Data prep and clean up!
###############################################################################################


# Set parameters for the run...
repo.loc <- "D:/Github/BBn_model/"
atow<-800*2.4384/10^6 # area of standard tow in km2
num.knots <- 10 # for SEAM
R0 <- 55 # for SEAM
qR <- 0.5 # for TLM
mods <- c("TLM")
n.mods <- length(mods)
retro.years <- 2010:2022
n.retro.years <- length(retro.years)

for(i in 1:n.mods)
{
  mod.select <- mods[i]
  # Bring in the data, only need to do that once...
  surv.early <- surv.dat %>% dplyr::filter(year < 2013)
  surv.late <- surv.dat %>% dplyr::filter(year %in% c(20:2019) & random ==1)
  surv.2020 <- surv.dat %>% dplyr::filter(year %in% c(2020))
  surv.20212 <- surv.dat %>% dplyr::filter(year %in% c(2021:2022) & random ==1)
  surv.data <- rbind(surv.early,surv.late,surv.2020,surv.20212)
  live.subset <- surv.data %>% dplyr::filter(state == 'live')
  dead.subset <- surv.data %>% dplyr::filter(state== "dead")
  live.input <- data.frame(I = live.subset$com.bm, IR = live.subset$rec.bm,year = live.subset$year,tow = live.subset$tow,tot.live.com = live.subset$com,lat = live.subset$lat,lon=live.subset$lon)
  clap.input <- data.frame(L = dead.subset$com,tow = dead.subset$tow,year = dead.subset$year)
  mod.input <- left_join(live.input,clap.input,by=c('tow','year'))
  mod.input$N <- round(mod.input$tot.live.com + mod.input$L)
  # Looks like there are no values > 0 but < 0.5, so the low clapper numbers should all round up to 1 (which makes sense as you'd only get < 0.5 if we had tows twice as long as they should be)
  mod.input$L <- round(mod.input$L) 
  mod.input.sfs <- st_as_sf(mod.input,coords = c('lon','lat'),remove=F,crs = 4326)
  mod.input.sfs <- mod.input.sfs %>% st_transform(crs=32619)
  mod.input.sfs <- st_intersection(mod.input.sfs,bbn.shape)
  mod.input.sfs$I <- mod.input.sfs$I/atow
  mod.input.sfs$IR <- mod.input.sfs$IR/atow
  
  # Grab the growth data
  growths <- data.frame(g = mod.dat$g, gR = mod.dat$gR,year = mod.dat$year)
  bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] <- bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] -1
  bbn.fish.sfs <- st_as_sf(bbn.fish,coords = c("lon","lat"),remove =F, crs = 4326)
  bbn.fish.sfs <- bbn.fish.sfs %>% st_transform(crs= 32619)
  # Now lets clip this to be data inside of our bbn boundary.
  bbn.fish.sfs <- st_intersection(bbn.fish.sfs,bbn.shape)
  
  for(j in 1:n.retro.years)
  {
    years <- 1994:retro.years[j]
    NY <- length(years)
    
    # if(retro.years[j] >= 2020)
    # {
    #   mod.input.sfs[nrow(mod.input.sfs)+1,] <- mod.input.sfs[nrow(mod.input.sfs),]
    #   #mod.input.sf[nrow(mod.input.sf),] <- mod.input.sf[nrow(mod.input.sf)-1,]
    #   mod.input.sfs$year[nrow(mod.input.sfs)] <- 2020
    #   mod.input.sfs$Year[nrow(mod.input.sfs)] <- which(years == 2020)
    #   mod.input.sfs$I[nrow(mod.input.sfs)] <- NA
    #   mod.input.sfs$IR[nrow(mod.input.sfs)] <- NA
    #   mod.input.sfs$tot.live.com[nrow(mod.input.sfs)] <- NA
    #   mod.input.sfs$L[nrow(mod.input.sfs)] <- 0
    #   mod.input.sfs$N[nrow(mod.input.sfs)] <- 0
    # }
    # OK, so step 1 here is getting the model input that Raphael needs for the model
    # The survey data....
    mod.input.sfs$Year <- mod.input.sfs$year - (min(years)-1)
    mod.input.sf <- mod.input.sfs %>% dplyr::filter(year %in% years)
    growth <- growths %>% dplyr::filter(year %in% c(years,(max(years)+1)))
    # Now we can clip both of these to subset it to the data that I think we need for the analysis....
    # Subset the fishery data as necessary
    bbn.fish.sf <- bbn.fish.sfs %>% dplyr::filter(survey.year %in% years)
    # OK, so now let's see if we can use the catch knot thing Raph made to split this up withing the BBn domain
    #We just need 3 columns for this
    catch.sf <- bbn.fish.sf %>% dplyr::select(pro.repwt,survey.year)
    names(catch.sf) <- c("Catch","Year","geometry")
    bbn.mesh <- setup_mesh(mod.input.sf,model_bound = bbn.shape,nknot=num.knots, max.edge = c(3,10),cutoff=2,seed=66) # Seeds 20 and 66 work
    bbn.mesh.sf <- inla.mesh2sf(bbn.mesh$mesh)
    bbn.mesh.sf$triangles$geometry <- bbn.mesh.sf$triangles$geometry*1000
    bbn.mesh.sf$vertices$geometry <- bbn.mesh.sf$vertices$geometry*1000
    st_crs(bbn.mesh.sf$triangles) <- 32619
    # Now make the prediction grid
    pred.grid<-setup_pred_grid(knots=bbn.mesh$knots,model_bound=bbn.mesh$utm_bound)
    # For the moment we need to have this starting at year 1.
    catch.sf$Year <-  catch.sf$Year - (min(years)-1)
    catchy <- catch_spread(catch = catch.sf,knots = bbn.mesh$knots)
    # A hack until Raph fixes up with model
    catchy$sum_catches <- catchy$sum_catches[,-1]
    # Catch for TLM
    catch.tlm <- catch.sf %>% group_by(Year,.drop=F) %>% dplyr::summarise(catch = sum(Catch,na.rm=T))
    
    
    # The SEBDAM set data.
    if(mod.select != "TLM")
    {
      set_data<-data_setup(data=mod.input.sf,growths=growth,catch=as.data.frame(catchy$sum_catches),
                           model="SEBDAM",mesh=bbn.mesh$mesh,obs_mort=T,prior=T,prior_pars=c(10,12),fix_m = 0.3,
                           mult_qI=T,spat_approach="spde",
                           knot_obj=bbn.mesh$knots,knot_area=pred.grid$area,separate_R_aniso = T,all_se=F)
      
      # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
      set_data$par$log_m0 <- 0
      set_data$par$log_R0 <- log(R0) # 5.3 = 200, 5.01 = 150, 4 = 55, 6 = 400
      #set_data$par$log_qR <- -1.5
      #set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA),log_qR = factor(NA))
      set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA))
      #set_data$map <-list(log_m0=factor(NA))
    } # end if(mod.select != "TLM")
    
    
    
    # A TLM version of the same...
    if(mod.select == "TLM")
    {
      names(growth) <- c("g","gR","Year")
      set_data<-data_setup(data=as.data.frame(mod.input.sf),growths=growth[,1:2],catch=catch.tlm$catch, model="TLM",obs_mort=TRUE,prior=TRUE)
      # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
      set_data$par$log_q_R <- log(qR) # Around 0.5, similar to some of the seam models.
     
      set_data$map <-list(log_q_R=factor(NA))
    } # end if(mod.select == "TLM")
    
    mod.fit<-fit_model(set_data,silent=F)
    
    # Now save the results appropriately
    if(mod.select != "TLM")
    {
      m0.par <- exp(set_data$par$log_m0)
      scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",m0.par,"_R0_",R0,"_",num.knots,"_knots")
      saveRDS(mod.fit,paste0(repo.loc,"Results/Models/BBn/Retros/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
      saveRDS(bbn.mesh,paste0(repo.loc,"Results/Models/BBn/Retros/BBn_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
      saveRDS(pred.grid,paste0(repo.loc,"Results/Models/BBn/Retros/BBn_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
      
    }
    
    if(mod.select == "TLM") 
    {
      qR.par <- sub("0.","0_",qR)
      scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR.par)
      saveRDS(mod.fit,paste0(repo.loc,"Results/Models/BBn/Retros/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
      
    }

}} # End the i and j loops


