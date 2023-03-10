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
surv.dat <- readRDS('D:/Github/BBn_model/Results/Sab_surv.dat.RDS')
# Need to get condition factor out of here too

mod.dat <- readRDS('D:/Github/BBn_model/Results/Sab_model.dat.RDS')

Sab.fishs <- readRDS('D:/Github/BBn_model/Results/Fishery_data/Sab_fish.dat.RDS')
Sab.fishs$pro.repwt <- Sab.fishs$pro.repwt/1000 # It looks like what I saved is already in tonnes.



# Set parameters for the run...
repo.loc <- "D:/Github/BBn_model/"
atow<-800*2.4384/10^6 # area of standard tow in km2
num.knots <- 4 # for SEAM
R0 <- 55 # for SEAM
qR <- 0.5 # for TLM
mods <- c("TLM")
n.mods <- length(mods)
retro.years <- 2010:2022
n.retro.years <- length(retro.years)
c_sys = 32620

for(i in 1:n.mods)
{
  mod.select <- mods[i]
  # Transform Sable to 32620
  Sab.shape <- Sab.shape %>% st_transform(crs = c_sys) # Sab is totally in 32620 border so think they are basically equivalent options here
  # Just going to use the core area to see if that helps model and the prediction function...
  # The survey data....
  live.subset <- surv.dat %>% dplyr::filter(state == 'live')
  dead.subset <- surv.dat %>% dplyr::filter(state== "dead")
  live.input <- data.frame(I = live.subset$com.bm, IR = live.subset$rec.bm,year = live.subset$year,tow = live.subset$tow,tot.live.com = live.subset$com,lat = live.subset$lat,lon=live.subset$lon)
  
  clap.input <- data.frame(L = dead.subset$com,tow = dead.subset$tow,year = dead.subset$year)
  mod.input <- left_join(live.input,clap.input,by=c('tow','year'))
  mod.input$N <- round(mod.input$tot.live.com + mod.input$L)
  # Looks like there are no values > 0 but < 0.5, so the low clapper numbers should all round up to 1 (which makes sense as you'd only get < 0.5 if we had tows twice as long as they should be)
  mod.input$L <- round(mod.input$L) 
  mod.input.sfs <- st_as_sf(mod.input,coords = c('lon','lat'),remove=F,crs = 4326)
  mod.input.sfs <- mod.input.sfs %>% st_transform(crs=c_sys)
  mod.input.sfs <- st_intersection(mod.input.sfs,Sab.shape)

  # Now I need to get the I and IR into kg/km^2
  mod.input.sfs$I <- mod.input.sfs$I/atow
  mod.input.sfs$IR <- mod.input.sfs$IR/atow
  

  
  # Growth!! 
  mod.growth.dat <- mod.dat
  # # Grab the growth data, we have ageing data from 1980's that I'm going to use to calculate growth here.
  # Data is coming from ageing data in 1989, found here.... Y:\Offshore\Assessment\Data\Ageing\archive\old_ageing_from_Amy_2022\SAB height at age 1989_2.pdf 
  L.inf <- 136.628
  #to <- 1.337 # So this uses a 1 year offset that we no longer believe in, going to make this 0.337 to align more with what we now do...
  to <- 0.337
  K <- 0.2269
  
  
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
  growths <- data.frame(g = c(mod.growth.dat$g,mod.growth.dat$g2[nrow(mod.growth.dat)]), 
                       gR = c(mod.growth.dat$gR,mod.growth.dat$gR2[nrow(mod.growth.dat)]),
                       year =  c(mod.growth.dat$year,max(mod.growth.dat$year+1)))

  for(j in 1:n.retro.years)
  {
      years <- 1994:retro.years[j]
      NY <- length(years)
      #survey.obj$Sab$model.dat
      growth <- growths %>% dplyr::filter(year %in% c(years,(max(years)+1)))
      # Now clip mod.input.sfs to the right nubmer of years...
      mod.input.sfs$Year <- mod.input.sfs$year - (min(years)-1)
      # Adding missing survey years.  THe L and N both need 'data', but 0's are fine, so that'll do the trick!
      mod.input.sfs[nrow(mod.input.sfs)+1,] <- mod.input.sfs[nrow(mod.input.sfs),]
      #mod.input.sfs[nrow(mod.input.sfs),] <- mod.input.sfs[nrow(mod.input.sfs)-1,]
      mod.input.sfs$year[nrow(mod.input.sfs)] <- 2015
      if(any(years == 2015))
      {
        mod.input.sfs$Year[nrow(mod.input.sfs)] <- which(years == 2015)
        mod.input.sfs$I[nrow(mod.input.sfs)] <- NA
        mod.input.sfs$IR[nrow(mod.input.sfs)] <- NA
        mod.input.sfs$tot.live.com[nrow(mod.input.sfs)] <- NA
        mod.input.sfs$L[nrow(mod.input.sfs)] <- 0
        mod.input.sfs$N[nrow(mod.input.sfs)] <- 0
      }
      # And repeat for 2020
      if(any(years == 2020))
      {
        mod.input.sfs[nrow(mod.input.sfs)+1,] <- mod.input.sfs[nrow(mod.input.sfs),]
        #mod.input.sfs[nrow(mod.input.sfs),] <- mod.input.sfs[nrow(mod.input.sfs)-1,]
        mod.input.sfs$year[nrow(mod.input.sfs)] <- 2020
        mod.input.sfs$Year[nrow(mod.input.sfs)] <- which(years == 2020)
        mod.input.sfs$I[nrow(mod.input.sfs)] <- NA
        mod.input.sfs$IR[nrow(mod.input.sfs)] <- NA
        mod.input.sfs$tot.live.com[nrow(mod.input.sfs)] <- NA
        mod.input.sfs$L[nrow(mod.input.sfs)] <- 0
        mod.input.sfs$N[nrow(mod.input.sfs)] <- 0
      }
      # Now filter on years
      mod.input.sf <- mod.input.sfs %>% dplyr::filter(year %in% years)
      
      # Subset the fishery data to the correct years. We need to take on next year catch too..)
      Sab.fish <- Sab.fishs %>% dplyr::filter(survey.year %in% years) # c(years,(max(years)+1))
  
      Sab.fish.sf <- st_as_sf(Sab.fish,coords = c("lon","lat"),remove =F, crs = 4326)
      Sab.fish.sf <- Sab.fish.sf %>% st_transform(crs= c_sys)
      # Now lets clip this to be data inside of our Sab boundary.
      Sab.fish.sf <- st_intersection(Sab.fish.sf,Sab.shape)
      
      catch.sf <- Sab.fish.sf %>% dplyr::select(pro.repwt,survey.year)
      names(catch.sf) <- c("Catch","Year","geometry")
       
      # For the moment we need to have this starting at year 1.
      catch.sf$Year <-  catch.sf$Year - (min(years)-1)
      Sab.mesh <- setup_mesh(mod.input.sf,model_bound = Sab.shape,nknot=num.knots, max.edge = c(8,20),cutoff=2.5,seed=34) 
      Sab.mesh.sf <- inla.mesh2sf(Sab.mesh$mesh)
      Sab.mesh.sf$triangles$geometry <- Sab.mesh.sf$triangles$geometry*1000
      Sab.mesh.sf$vertices$geometry <- Sab.mesh.sf$vertices$geometry*1000
      st_crs(Sab.mesh.sf$triangles) <- c_sys
      st_crs(Sab.mesh.sf$vertices) <- c_sys
      
      # Now make the prediction grid
      pred.grid<-setup_pred_grid(knots=Sab.mesh$knots,model_bound=Sab.mesh$utm_bound)
      st_crs(pred.grid$grid) <- c_sys
      
      #Sebdam catches
      catchy <- catch_spread(catch = catch.sf,knots = Sab.mesh$knots)
      catchy$sum_catches <- catchy$sum_catches[,-1] # A hack until Raph gets new model code up and running.
      # For now we need to toss the first column from there
      #TLM catch 
      catch.tlm <- catch.sf %>% group_by(Year,.drop=F) %>% dplyr::summarise(catch = sum(Catch,na.rm=T))
      
      
      #SEBDAM version 
      if(mod.select == "SEAM")
      {
      set_data<-data_setup(data=mod.input.sf,growths=growth,catch=as.data.frame(catchy$sum_catches),
                           model="SEBDAM",mesh=Sab.mesh$mesh,obs_mort=T, prior=TRUE,prior_pars=c(10,12),
                           mult_qI=T,spat_approach="spde",
                           knot_obj=Sab.mesh$knots,knot_area=pred.grid$area,separate_R_aniso =T,all_se=FALSE)
      # So this will fix the mean value of m0 to be whatever the intial value is set at.  Let's see what happens!
      set_data$par$log_m0 <- 0 # 0 = 1
      set_data$par$log_R0 <- log(R0) # 5.3 = 200, 5 = 148, 4 = 55, 5.9915 = 400, 4.606 = 100
      set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA))
  
      }
      #TLM version, Dude is this ever sensitive to the q priors! (5,12) actually looks solid in terms of results... maybe we can get so lucky with SEBDAM :-)
      # Note that the Catch time series should be 1 year longer than the survey data here!!
      if(mod.select == "TLM")
      {
      set_data<-data_setup(data=as.data.frame(mod.input.sf),growths=growth[,1:2],catch=catch.tlm$catch,
                            model="TLM",obs_mort=T,prior=T,prior_pars=c(10,12))
      set_data$par$log_q_R <- log(qR) #-0.7 # Around 0.5, similar to some of the seam models.
  
      set_data$map <-list(log_q_R=factor(NA))
      }
  
  
      mod.fit<-fit_model(set_data,silent=F)
      
      # Now save the results appropriately
      if(mod.select != "TLM") 
      {
        m0.par <- exp(set_data$par$log_m0)
        r0.par <- signif(exp(set_data$par$log_R0),digits=2)
        scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",m0.par,"_R0_",r0.par,"_",num.knots,"_knots")
        saveRDS(mod.fit,paste0(repo.loc,"Results/Models/Sab/Retros/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
        saveRDS(Sab.mesh,paste0(repo.loc,"Results/Models/Sab/Retros/Sab_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
        saveRDS(pred.grid,paste0(repo.loc,"Results/Models/Sab/Retros/Sab_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
        
      }
      
      if(mod.select == "TLM") 
      {
        qR.par <- sub("0.","0_",signif(exp(set_data$par$log_q_R),digits=2))
        scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR.par)
        saveRDS(mod.fit,paste0(repo.loc,"Results/Models/Sab/Retros/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
        
      }
  }} # End the i and j loops