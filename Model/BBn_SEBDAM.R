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
# Grab the growth data
growth <- data.frame(g = mod.dat$g, gR = mod.dat$gR,year = mod.dat$year)
growth <- growth %>% dplyr::filter(year %in% c(years,(max(years)+1)))
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

# DK NOTE: Now this is going to get confusing for us and we may want to tweak SEBDAM for this, but that's a down the road job, not a playing around with model job
# But based on the indexing in SEBDAM, I am going to change how we index the survey year data from what we have done with offshore traditionally.
# Historically anything from the last half of the year goes into the following years, eg. survey.year 2002 = June 2001- May 2002.  
# But in SEBDAM we have (B(t-1) - C(t-1)), so let's say we have year 2002 survey biomass, this says we remove the 20002 catch from that
# we want that catch to be the catch from June 2002 to May 2003, i.e. we remove the catch before we allow the population to grow
# This is what we do in our current model, but we have a different index (C(t) on our model.
# Basically survey year 2002 = June 2002 - May 2003 now
#DK note: We probably should think more about the survey year fun and how exactly we want to handle removal of catch in our models.
# We don't have removals for 2009, we need something for that, so we're adding that in here...
bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] <- bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] -1

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
bbn.fish <- bbn.fish %>% dplyr::filter(survey.year %in% years)
bbn.fish.sf <- bbn.fish.sf %>% dplyr::filter(survey.year %in% years)
# OK, so now let's see if we can use the catch knot thing Raph made to split this up withing the BBn domain
#We just need 3 columns for this
catch.sf <- bbn.fish.sf %>% dplyr::select(pro.repwt,survey.year)
names(catch.sf) <- c("Catch","Year","geometry")
#catch.sf$geometry <- catch.sf$geometry/1000

# A map...
 #b.map <- pecjector(area= "BBn",c_sys = 32619,add_layer = list(land = 'grey',eez = 'eez' , nafo = 'main',sfa = 'offshore',survey = c("offshore","outline")),txt.size=8,axes = "DM")
# bbn.fish.map <- b.map + geom_sf(data = bbn.fish.sf) + facet_wrap(~year) + geom_sf(data= bbn.shape,fill = NA)
# bbn.fish.map

# Set up our mesh...
#bbn.mesh <- setup_mesh(catch.sf,model_bound = bbn.shape,nknot=8,seed=20) # Seeds 20 and 66 work
#bbn.shape$geometry <-  bbn.shape$geometry/1000
#st_crs(bbn.shape) <- 32619
#mod.input.sf$geometry <-  mod.input.sf$geometry/1000
#st_crs(mod.input.sf) <- 32619

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
ggplot(bbn.mesh.sf$triangles) + geom_sf() + geom_sf(data= bbn.shape,fill = NA,color = 'blue',size=2) + geom_sf(data = knots.sf,fill = NA)
# Now make the prediction grid
pred.grid<-setup_pred_grid(knots=bbn.mesh$knots,model_bound=bbn.mesh$utm_bound)
# Plot the grid
#ggplot(pred.grid$grid) + geom_sf(aes(fill = as.factor(knotID))) + scale_fill_viridis_d()
# Get the knots on the right scale
#knots.on.right.scale <- bbn.mesh$knots
#knots.on.right.scale$centers <- knots.on.right.scale$centers*1000

# For the moment we need to have this starting at year 1.
catch.sf$Year <-  catch.sf$Year - (min(years)-1)
catchy <- catch_spread(catch = catch.sf,knots = bbn.mesh$knots)
# A hack until Raph fixes up with model
catchy$sum_catches <- catchy$sum_catches[,-1]
# For now we need to toss the first column from there
catchy
# Catch for TLM
catch.tlm <- catch.sf %>% group_by(Year,.drop=F) %>% dplyr::summarise(catch = sum(Catch,na.rm=T))



# ggplot(mod.input.sf) + geom_boxplot(aes(x=Year,y=I,group=Year)) + scale_y_log10()
# ggplot(mod.input.sf) + geom_boxplot(aes(x=Year,y=IR,group=Year))+ scale_y_log10()
# ggplot(mod.input.sf) + geom_boxplot(aes(x=Year,y=N,group=Year))+ scale_y_log10()
# ggplot(mod.input.sf) + geom_boxplot(aes(x=Year,y=L,group=Year))+ scale_y_log10()


# The SEBDAM set data.
if(mod.select != "TLM")
{
  set_data<-data_setup(data=mod.input.sf,growths=growth,catch=as.data.frame(catchy$sum_catches),
                       model="SEBDAM",mesh=bbn.mesh$mesh,obs_mort=T,prior=T,prior_pars=c(10,12),fix_m = 0.3,
                       mult_qI=T,spat_approach="spde",
                       knot_obj=bbn.mesh$knots,knot_area=pred.grid$area,separate_R_aniso = T,all_se=F)
  str(set_data)
  
  # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
  set_data$par$log_m0 <- 0
  set_data$par$log_R0 <- 4 # 5.3 = 200, 5.01 = 150, 4 = 55, 6 = 400
  #set_data$par$log_qR <- -1.5
  #set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA),log_qR = factor(NA))
  set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA))
  #set_data$map <-list(log_m0=factor(NA))
} # end if(mod.select != "TLM")



# A TLM version of the same...
if(mod.select == "TLM")
{
  set_data<-data_setup(data=as.data.frame(mod.input.sf),growths=growth[,1:2],catch=catch.tlm$catch, model="TLM",obs_mort=TRUE,prior=TRUE)
  # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
  set_data$par$log_q_R <- -0.7 # Around 0.5, similar to some of the seam models.
  #set_data$par$log_R0 <- 4 # 5.3 = 200, 5 = 148, 4 = 55, 6 = 400
  set_data$par$log_qR <- -1.5
  #set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA),log_qR = factor(NA))
  set_data$map <-list(log_q_R=factor(NA))
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
if(mod.select != "TLM")
{
  m0.par <- exp(set_data$par$log_m0)
  r0.par <- signif(exp(set_data$par$log_R0),digits=2)
  scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",m0.par,"_R0_",r0.par,"_",num.knots,"_knots")
  saveRDS(mod.fit,paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
  saveRDS(bbn.mesh,paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
  saveRDS(pred.grid,paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
}

if(mod.select == "TLM") 
{
  qR.par <- sub("0.","0_",signif(exp(set_data$par$log_q_R),digits=2))
  scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR.par)
  saveRDS(mod.fit,paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
}

# The model using only proper survey strata
# m0.par <- exp(set_data$par$log_m0)
# r0.par <- signif(exp(set_data$par$log_R0),digits=2)
# scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",m0.par,"_R0_",r0.par,"_",num.knots,"_knots_No_extra_stations")
# saveRDS(mod.fit,paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
# saveRDS(bbn.mesh,paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
# saveRDS(pred.grid,paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))



################################################### End the initial model runs ###########################################
################################################### End the initial model runs ###########################################
################################################### End the initial model runs ###########################################






##################### Now load the model and make the figures! ##############################################

atow<-800*2.4384/10^6 # area of standard tow in km2
num.knots <- 20
RO <- 150
qR  <- "0_5" # This is just for TLM models
years <- 1994:2022
NY <- length(years)
c_sys <- 32619
theme_set(theme_few(base_size = 22))
repo.loc <- "D:/Github/BBn_model/"
mod.select <- "SEAM"
################################################### End the initial model runs ###########################################
### Make the figures for the models



if(mod.select != "TLM") scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_1_R0_",RO,"_",num.knots,"_knots")
if(mod.select == "TLM") scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR)
# If we are going with the no extra station model this is our scenario.
#scenario.select <- "1994_2022_vary_m_m0_1_R0_150_10_knots_No_extra_stations"
mod.fit <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
if(mod.select != "TLM") catchy <- mod.fit$obj$env$data$C*mod.fit$obj$env$data$area # Get this into tonnes from catch density.
if(mod.select == "TLM") catchy <- mod.fit$obj$env$data$C

# This is only needed for SEAM.
if(mod.select != "TLM")
{
  pred.grid <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
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
  #
  # SS model mu(t) <- C(t) / (B(t) + C(t)) because our model is B(t) <- B(t-1) - C(t) and C(2017) is June 2016-Aug 2017.
  #  mu[2017] <- C[June 2016-Aug 2017]/(B[2017]+C[June 2016-Aug 2017]) 
  # TLM and SEAM don't calculate mu, so we do it manually here, to be analogous...
  # SEAM/TLM mu(t) <- C(t-1) / (B(t) + C(t-1)) because our model is B(t) <- B(t-1) - C(t-1) and C(2016) is now June 2016-Aug 2017.
  # mu[2017] <- C[June 2016-Aug 2017]/(B[2017]+C[June 2016-Aug 2017]) 
  F.dat<-data.frame(B=as.vector(mod.fit$report$areaB[,-ncol(mod.fit$report$areaB)]/1000),
                    C = c(rep(NA,num.knots),as.vector(as.matrix(catchy[,-ncol(catchy)]))), Year=matYear1, knotID=knots1)
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
    if(any(B2 < 0)) {B2[B2<= 0] <- 5; print("HEADS UP!! You have a negative B2 estimate (look for the 5s in B.exp).")} # This feels very hacky for the moment, I think once Raph has the catch alignment fixed up this will not matter anymore, but not 100% sure!
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
  saveRDS(Bdiff.comp.dat.plot,paste0(repo.loc,"Results/Models/BBn/BBn_",mod.select,"_model_output_",scenario.select,"_B_differnce.Rds"))
  
  
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
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_viridis_c(breaks = b.brk, labels=b.lab, name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_biomass.png"),spatial.B.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.B.plot<- ggplot() + geom_sf(data=B.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=log(B)),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_viridis_c(breaks = b.brk, labels=b.lab, name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_biomass_no_missing_surveys.png"),spatial.B.plot,base_width = 10,base_height = 10)
  # Subest to 4 years
  spatial.B.plot<- ggplot() + geom_sf(data=B.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(B)),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_viridis_c(breaks = b.brk, labels=b.lab, name="Predicted Biomass \nDensity (kg\U2022km\U207B\U00B2)",option = "A",begin=0.2) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_biomass_4_years.png"),spatial.B.plot,base_width = 10,base_height =7)
  
  #R
  r.brk <- pretty(log(R.dat.plot$R))
  r.lab <- signif(exp(r.brk),digits=2)
  
  spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot,aes(fill=log(R)),col='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_recruits.png"),spatial.R.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=log(R)),col='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_recruits_no_missing_surveys.png"),spatial.R.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.R.plot<-  ggplot() + geom_sf(data=R.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(R)),col='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = r.brk, labels = r.lab,name="Predicted Recruit \nDensity (kg\U2022km\U207B\U00B2)",end=0.8)+ 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_recruits_4_years.png"),spatial.R.plot,base_width = 10,base_height = 7)
  
  #m
  m.brk <- log(c(0.003,0.007,0.02,0.05,0.15,0.4,1))
  #m.brk <- pretty(log(m.dat.plot$m))
  m.lab <- signif(exp(m.brk),digits=2)
  
  spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot %>% dplyr::filter(!Year %in% c(2023)),aes(fill=log(m)),color='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality (Inst)",option = "B",direction =1,begin = 0.2,end=1) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_mort.png"),spatial.m.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=log(m)),color='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality (Inst)",option = "B",direction =1,begin = 0.2,end=1) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_mort_no_missing_surveys.png"),spatial.m.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.m.plot <-  ggplot() + geom_sf(data=m.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(m)),color='grey')+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year)+ 
    scale_fill_viridis_c(breaks = m.brk, labels = m.lab,name="Predicted Natural \nMortality (Inst)",option = "B",direction =1,begin = 0.2,end=1) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_mort_4_years.png"),spatial.m.plot,base_width = 10,base_height = 7)
  
  # q 
  spatial.q.plot <- ggplot() + geom_sf(data=q.dat.plot,aes(fill=qI),col=NA)+
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_viridis_c(name="Predicted catchability (qI)",option = "C",begin = 0.2,end =0.8) + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_catchability.png"),spatial.q.plot,base_width = 10,base_height = 10)
  
  
  # OK, so lets try an make a map of the spatial exploitation rates, not sure if this is all correct yet.
  # Log rescalling doesn't work because we have 0's here.
  # OK, so lets try an make a map of the spatial exploitation rates, not sure if this is all correct yet.
  F.dat.plot$exp.na <- NA
  F.dat.plot$exp.na[F.dat.plot$exploit != 0] <- F.dat.plot$exploit[F.dat.plot$exploit != 0]
  
  #e.brk <- log(c(0.0015,0.01,0.08,0.4))
  e.brk <- pretty(log(F.dat.plot$exp.na))
  e.lab <- signif(exp(e.brk),digits=2)
  
  
  spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot,aes(fill=log(exp.na)),color='grey') +
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year) + 
    scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/Spatial_exploit.png"),spatial.exploit.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=log(exp.na)),color='grey') +
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year) + 
    scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/Spatial_exploit_no_missing_surveys.png"),spatial.exploit.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.exploit.plot<- ggplot() + geom_sf(data=F.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=log(exp.na)),color='grey') +
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    facet_wrap(~Year) + 
    scale_fill_viridis_c(breaks = e.brk,labels = e.lab,name="Exploitation (Prop)",option = "D") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/Spatial_exploit_4_years.png"),spatial.exploit.plot,base_width = 10,base_height = 10)
  
  
  bd.brk <- pretty(Bdiff.comp.dat.plot$B.diff)
  bd.lab <- bd.brk#signif(exp(b.brk),digits=2)
  
  spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot,aes(fill=B.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_BDiff.png"),spatial.Bdiff.plot,base_width = 10,base_height = 10)
  #Remove missing survey years
  spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=B.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_BDiff_no_missing_surveys.png"),spatial.Bdiff.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=B.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_distiller(type = 'div',breaks = bd.brk, labels=bd.lab, name="Expected - Modeled \nBiomass (tonnes)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_BDiff_4_years.png"),spatial.Bdiff.plot,base_width = 10,base_height = 7)
  
  # Same plot but the percentage miss by cell.
  pb.brk <- pretty(Bdiff.comp.dat.plot$B.per.diff)
  pb.lab <- pb.brk#signif(exp(b.brk),digits=2)
  spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot,aes(fill=B.per.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_per_BDiff.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 10)
  # Remove missing survey years
  spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(!Year %in% c(2020,2023)),aes(fill=B.per.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_per_BDiff_no_missing_surveys.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 10)
  # 4 years
  spatial.per.Bdiff.plot<- ggplot() + geom_sf(data=Bdiff.comp.dat.plot %>% dplyr::filter(Year %in% c(2001,2009,2014,2019)),aes(fill=B.per.diff),color='grey')+
    facet_wrap(~Year)+ 
    scale_x_continuous(breaks = c(-60.3,-60.1,-59.9), labels = c("60°18'W","60°6'W","59°54'W")) +
    scale_y_continuous(breaks = c(42.6, 42.75, 42.9),labels = c("42°36'N","42°45'N","42°54'N")) +
    scale_fill_distiller(type = 'div',breaks = pb.brk, labels=pb.lab, name="Expected - Modeled \nBiomass (%)",palette = "RdBu") + 
    theme(axis.text.x=element_text(angle=-45,hjust=0))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Spatial_per_BDiff_4_years.png"),spatial.per.Bdiff.plot,base_width = 10,base_height = 7)
  
  theme_set(theme_few(base_size = 22))
  
  
  # And now plot them....
  m.comp <- ggplot(nat.mat.plt,aes(x=Year,y=m,group = method,color=method)) + geom_line(linewidth=1.5)  + 
    xlab("") + ylab("Natural Mortality (instantaneous)")+
    scale_color_manual(values = c("blue","orange","grey","black"))  + 
    theme(legend.title=element_blank()) + scale_x_continuous(breaks = seq(1980,2030,by=3))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/natural_mortality_comparisons.png"),m.comp,base_height = 6,base_width = 10)
  # Remove missing survey years
  m.comp <- ggplot(nat.mat.plt %>% dplyr::filter(Year < 2020)) + geom_line(aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  + 
    #geom_line(data = nat.mat.plt %>% dplyr::filter(Year %in% 2016:2019), aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  + 
    geom_line(data = nat.mat.plt %>% dplyr::filter(Year %in% 2021:2022), aes(x=Year,y=m,group = method,color=method),linewidth=1.5)  +  
    xlab("") + ylab("Natural Mortality (instantaneous)")+
    scale_color_manual(values = c("blue","orange","grey","black"))  + 
    theme(legend.title=element_blank()) + scale_x_continuous(breaks = seq(1980,2030,by=3))
  save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/natural_mortality_comparisons_no_missing_surveys.png"),m.comp,base_height = 6,base_width = 10)
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
if(mod.select != "TLM") catch.annual <- data.frame(totC = colSums(catchy), Year = years)
if(mod.select == "TLM") catch.annual <- data.frame(totC = catchy, Year = years)


pred.proc$log_processes <- pred.proc$log_processes %>% dplyr::filter(year < 2023)

# The catch data
# So catches from June 2021-May 2022 are called 2021 and removed from the 2021 survey biomass (this is different indexing from how we used to handle this for offshore)
#
# SS model mu(t) <- C(t) / (B(t) + C(t)) because our model is B(t) <- B(t-1) - C(t) and C(2017) is June 2016-Aug 2017.
#  mu[2017] <- C[June 2016-Aug 2017]/(B[2017]+C[June 2016-Aug 2017]) 
# TLM and SEAM don't calculate mu, so we do it manually here, to be analogous...
# SEAM/TLM mu(t) <- C(t-1) / (B(t) + C(t-1)) because our model is B(t) <- B(t-1) - C(t-1) and C(2016) is now June 2016-Aug 2017.
# mu[2017] <- C[June 2016-Aug 2017]/(B[2017]+C[June 2016-Aug 2017]) 
ann.exploit <- data.frame(year = years,B = exp(pred.proc$log_processes$log_B), Catch = c(NA,colSums(catchy[,-ncol(catchy)])),
                          B.LCI = pred.proc$log_processes$totB.LCI, B.UCI = pred.proc$log_processes$totB.UCI)
ann.exploit$exploit <- ann.exploit$Catch/(ann.exploit$B+ann.exploit$Catch)
ann.exploit$exploit.UCI <- ann.exploit$Catch/(ann.exploit$B.LCI+ann.exploit$Catch)
ann.exploit$exploit.LCI <- ann.exploit$Catch/(ann.exploit$B.UCI+ann.exploit$Catch)
ann.exploit$FM <- 1-exp(-ann.exploit$exploit)
ann.exploit$FM.LCI <- 1-exp(-ann.exploit$exploit.LCI)
ann.exploit$FM.UCI <- 1-exp(-ann.exploit$exploit.UCI)


# Biomass time series
bm.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.5,fill='blue',color='blue') +
  xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,2.35e4))
save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Biomass_time_series.png"),bm.ts.plot,base_width = 11,base_height = 8.5)
# Recruit time series
rec.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
  xlab("") + ylab("Recruit Biomass (tonnes)")  + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,1.75e4))
save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Recruit_time_series.png"),rec.ts.plot,base_width = 11,base_height = 8.5)
# Natural mortality time series...
mort.ts.plot <- ggplot(pred.proc$log_processes) + geom_line(aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) + 
  geom_ribbon(aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
  xlab("") + ylab("Natural mortality (Instantaneous)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.65))
save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_nat_mort_time_series.png"),mort.ts.plot,base_width = 11,base_height = 8.5)
# Explotation Rate Time Series
exploit.plot <- ggplot(ann.exploit) + geom_line(aes(x=year,y=exploit),size=1.5) + geom_ribbon(aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.5,fill='blue',color='blue') +
  xlab("") + ylab("Exploitation Rate (Proportional)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.35))
save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_exploit_time_series.png"),exploit.plot,base_width = 11,base_height = 8.5)


# Same plots but remvoing missing survey years....
bm.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2020))) + 
                            geom_line(aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            # geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  
                            # geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_B)),color='firebrick2',linewidth=1.5) +  
                            geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=totB.LCI,ymax=totB.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            ylim(c(0,2.35e4)) + xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Biomass_time_series_no_missing_surveys.png"),bm.ts.plot,base_width = 11,base_height = 8.5)
# Recruit time seris
rec.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2020))) + 
                            geom_line(aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            # geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  
                            # geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_R)),color='firebrick2',linewidth=1.5) +  
                            geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=totR.LCI,ymax=totR.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            xlab("") + ylab("Recruit Biomass (tonnes)") + ylim(c(0,1.75e4)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_Recruit_time_series_no_missing_surveys.png"),rec.ts.plot,base_width = 11,base_height = 8.5)
# Natural mortality time series...
mort.ts.plot <- ggplot(pred.proc$log_processes %>% dplyr::filter(year < c(2020))) + 
                            geom_line(aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  geom_ribbon(aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            #geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  
                            #geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2016:2019), aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            geom_line(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(year,exp(log_m)),color='firebrick2',linewidth=1.5) +  
                            geom_ribbon(data = pred.proc$log_processes %>% dplyr::filter(year %in% 2021:2022), aes(ymin=m.LCI,ymax=m.UCI,x=year),alpha=0.5,fill='blue',color='blue') + 
                            xlab("") + ylab("Natural mortality (Instantaneous)") + ylim(c(0,0.65)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_nat_mort_time_series_no_missing_surveys.png"),mort.ts.plot,base_width = 11,base_height = 8.5)


# Explotation Rate Time Series
exploit.plot <- ggplot(ann.exploit%>% dplyr::filter(year < c(2020))) + geom_line(aes(x=year,y=exploit),linewidth = 1.5) +
                            geom_ribbon(aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.5,fill='blue',color='blue') +
                            geom_line(data= ann.exploit%>% dplyr::filter(year %in% 2021:2022),aes(x=year,y=exploit),linewidth = 1.5) +
                            geom_ribbon(data= ann.exploit%>% dplyr::filter(year %in% 2021:2022),aes(ymin=exploit.LCI,ymax=exploit.UCI,x=year),alpha=0.5,fill='blue',color='blue') +
                            xlab("") + ylab("Exploitation Rate (Proportional)") + ylim(c(0,0.35)) + scale_x_continuous(breaks = seq(1980,2030,by=3))
save_plot(paste0(repo.loc,"Results/Figures/BBn/",mod.select,"_",scenario.select,"/BBn_exploit_time_series_no_missing_surveys.png"),exploit.plot,base_width = 11,base_height = 8.5)

