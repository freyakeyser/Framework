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
Sab.shape <- Sab.shape %>% st_transform(crs = 32620) # Sab is a solid 20
# Bring in the survey data
load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90.Rdata")
#load("D:/Framework/SFA_25_26_2024/Model/Data/testing_results_framework3.Rdata")
#surv.dat <- surv.dat$Sab
#saveRDS(surv.dat,'D:/Github/Sab_model/Results/Sab_surv.dat.RDS')
#surv.dat <- readRDS('D:/Framework/SFA_25_26_2024/Model/Data/Sab_surv.dat.RDS')
surv.dat <- surv.dat$Sab
mod.dat <- survey.obj$Sab$model.dat
Sab.fishs <- readRDS('D:/Framework/SFA_25_26_2024/Model/Data/Sab_fish.dat.RDS')
Sab.fishs$pro.repwt <- Sab.fishs$pro.repwt/1000
#### Finshed Data prep and clean up!


# Set parameters for the run...
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
atow<-800*2.4384/10^6 # area of standard tow in km2
years <- 1994:2022
num.knots <- 10 # for SEAM
#R0 <- 55 # for SEAM
lqr <- log(0.5)# This is for TMB (log recruit catchability) testing catchability of 0.5, test 0.3 and 0.1
l.init.m <- log(2) # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05R.size <- "75"
FR.size <- "90"
R.size <- "75"
mods <- c("SEAM")
n.mods <- length(mods)
retro.years <- 2010:2021
n.retro.years <- length(retro.years)
c_sys = 32620

for(i in 1:n.mods)
{
  mod.select <- mods[i]
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
  
  if(mod.select != "TLM")
  {
    # For the moment we need to have this starting at year 1.
    Sab.mesh <- setup_mesh(mod.input.sfs,model_bound = Sab.shape,nknot=num.knots, max.edge = c(8,20),cutoff=2.5,seed=34) 
    Sab.mesh.sf <- inla.mesh2sf(Sab.mesh$mesh)
    Sab.mesh.sf$triangles$geometry <- Sab.mesh.sf$triangles$geometry*1000
    Sab.mesh.sf$vertices$geometry <- Sab.mesh.sf$vertices$geometry*1000
    st_crs(Sab.mesh.sf$triangles) <- c_sys
    st_crs(Sab.mesh.sf$vertices) <- c_sys
    
    # Now make the prediction grid
    pred.grid<-setup_pred_grid(knots=Sab.mesh$knots,model_bound=Sab.mesh$utm_bound)
    st_crs(pred.grid$grid) <- c_sys
  }

  for(j in 1:n.retro.years)
  {
      years <- min(years):retro.years[j]
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
       
      # The SEBDAM set data.
      if(mod.select != "TLM")
      {
        catchy <- catch_spread(catch = catch.sf,knots = Sab.mesh$knots)
        catchy$sum_catches[,ncol(catchy$sum_catches)+1] <- 0
        set_data<-data_setup(data=mod.input.sf,growths=growth[,1:2],catch=catchy$sum_catches[],
                             model="SEBDAM",mesh=Sab.mesh$mesh,obs_mort=T,prior=T,prior_pars=c(10,12),#fix_m = 0.3,
                             mult_qI=T,spat_approach="spde",
                             knot_obj=Sab.mesh$knots,knot_area=pred.grid$area,separate_R_aniso = T,
                             all_se=T,weighted_mean_m = T)
        str(set_data)
        
        # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
        set_data$par$log_m0 <- l.init.m
        #set_data$par$log_R0 <- l.init.R 
        set_data$par$log_qR <- lqr
        #set_data$map <-list(log_m0=factor(NA),log_R0 = factor(NA),log_qR = factor(NA))
        set_data$map <-list(log_m0=factor(NA),log_qR = factor(NA))
        #set_data$map <-list(log_m0=factor(NA))
      } # end if(mod.select != "TLM")

      # A TLM version of the same...
      if(mod.select == "TLM")
      {
        catch.tlm <- catch.sf %>% group_by(Year,.drop=F) %>% dplyr::summarise(catch = sum(Catch,na.rm=T))
        catch.tlm[nrow(catch.tlm)+1,] <- catch.tlm[nrow(catch.tlm),]
        catch.tlm$Year[nrow(catch.tlm)] <- max(catch.tlm$Year) + 1
        catch.tlm$catch[nrow(catch.tlm)] <- 0
        set_data<-data_setup(data=as.data.frame(mod.input.sf),growths=growth[,1:2],
                             catch=catch.tlm$catch, model="TLM",obs_mort=TRUE,prior=TRUE)
        # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
        #set_data<-fix_param(obj=set_data, pars = list(log_q_R=lqr))
        set_data$par$log_q_R <- lqr # 
        set_data$map <-list(log_q_R=factor(NA))
      } # end if(mod.select == "TLM")
      
      mod.fit<-fit_model(set_data,silent=F)
      
      # Now save the results appropriately
      if(mod.select != "TLM") 
      {
        scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",exp(l.init.m),"_qR_",exp(lqr),"_",num.knots,"_knots")
        saveRDS(mod.fit,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
        saveRDS(Sab.mesh,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
        saveRDS(pred.grid,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
      }
      
      if(mod.select == "TLM") 
      {
        
        scenario.select <- paste0(min(years),"_",max(years),"_qR_",exp(lqr))
        saveRDS(mod.fit,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
      }
  }} # End the i and j loops



###################  SECTION 2 Make the Retro plots ##############################  SECTION 2 Make the Retro plots ###########
# Now make the retrospective plots...

# Set parameters for the run...
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
years <- 1994:2022
R.size <- "75"
FR.size <- "90"
num.knots <- 10 # Going to test 10
lqr <- log(0.5)
l.init.m <- log(2) # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05
# The survey biomass index for 1994 says there were 249 tonnes of recruits that year.
#l.init.R <- log(250) # Going to test 100, 250, and 500.

mod.select <- "SEAM"
retro.years <- 2010:2021
n.retro.years <- length(retro.years)

#scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",exp(l.init.m),"_R0_",exp(l.init.R),"_",num.knots,"_knots")
#tst <- readRDS("D:/framework/SFA_25_26_2024/Model/Results/Sab/R_75_FR_90/Retros/Sab_SEAM_model_output_1994_2010_vary_m_m0_0.4_qR_0.5_4_knots.Rds")

mod.fit <- NULL
retro.trends <- NULL
pred.proc <- NULL
for(j in retro.years)
{
  if(mod.select != "TLM") scenario.select <- paste0(min(years),"_",j,"_vary_m_m0_",exp(l.init.m),"_qR_",exp(lqr),"_",num.knots,"_knots")
  if(mod.select == "TLM") scenario.select <- paste0(min(years),"_",j,"_qR_",exp(lqr))
  mod.fit[[j]] <- readRDS(paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
  pred.proc[[j]] <- get_processes(mod.fit[[j]])
  
  if(mod.select == "TLM")
  {
    retro.trends[[as.character(j)]] <- data.frame(years = min(years):(j+1),
                                                  B =     pred.proc[[j]]$processes$B,
                                                  B.LCI = exp(pred.proc[[j]]$log_processes$log_B - 1.96*pred.proc[[j]]$log_processes$se_log_B),
                                                  B.UCI = exp(pred.proc[[j]]$log_processes$log_B + 1.96*pred.proc[[j]]$log_processes$se_log_B),
                                                  R =     pred.proc[[j]]$processes$R,
                                                  R.LCI = exp(pred.proc[[j]]$log_processes$log_R - 1.96*pred.proc[[j]]$log_processes$se_log_R),
                                                  R.UCI = exp(pred.proc[[j]]$log_processes$log_R + 1.96*pred.proc[[j]]$log_processes$se_log_R),
                                                  m =    pred.proc[[j]]$processes$m,
                                                  m.LCI = exp(pred.proc[[j]]$log_processes$log_m - 1.96*pred.proc[[j]]$log_processes$se_log_m),
                                                  m.UCI = exp(pred.proc[[j]]$log_processes$log_m + 1.96*pred.proc[[j]]$log_processes$se_log_m),
                                                  retro.year = as.character(j), mod = mod.select, scenario = scenario.select)
    retro.trends[[as.character(j)]] <- retro.trends[[as.character(j)]][-nrow(retro.trends[[as.character(j)]]),]
  }
  
  if(mod.select != "TLM")
  {
    # Get the overall estimates + the 95% CI
    retro.trends[[as.character(j)]] <- data.frame(years = min(years):(j+1),
                                                  B = exp(pred.proc[[j]]$log_tot_frame$log_totB),
                                                  R = exp(pred.proc[[j]]$log_tot_frame$log_totR),
                                                  m = pred.proc[[j]]$totals$mean_m,
                                                  B.LCI = exp(pred.proc[[j]]$log_tot_frame$log_totB - 1.96*pred.proc[[j]]$log_tot_frame$se_log_totB),
                                                  B.UCI = exp(pred.proc[[j]]$log_tot_frame$log_totB + 1.96*pred.proc[[j]]$log_tot_frame$se_log_totB),
                                                  R.LCI = exp(pred.proc[[j]]$log_tot_frame$log_totR - 1.96*pred.proc[[j]]$log_tot_frame$se_log_totR),
                                                  R.UCI = exp(pred.proc[[j]]$log_tot_frame$log_totR + 1.96*pred.proc[[j]]$log_tot_frame$se_log_totR),
                                                  m.LCI = exp(pred.proc[[j]]$log_tot_frame$log_mean_m - 1.96*pred.proc[[j]]$log_tot_frame$se_log_mean_m),
                                                  m.UCI = exp(pred.proc[[j]]$log_tot_frame$log_mean_m + 1.96*pred.proc[[j]]$log_tot_frame$se_log_mean_m),
                                                  retro.year = as.character(j), mod = mod.select, scenario = scenario.select)
  }
}

retros <- do.call("rbind",retro.trends)

# Remove the projection year for SEAM (does nothing for TLM)
retro <- retros[!is.na(retros$R),]

if(mod.select != "TLM") scenario.select <- paste0(min(years),"_",2022,"_vary_m_m0_",exp(l.init.m),"_qR_",exp(lqr),"_",num.knots,"_knots")
if(mod.select == "TLM") scenario.select <- paste0(min(years),"_",2022,"_qR_",qR)
base.mod <- readRDS(paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
pred.base <- get_processes(base.mod)

if(mod.select == "TLM")
{
  base.trend <-  data.frame(years = min(years):(j+1),
                            B =     pred.base$processes$B,
                            B.LCI = exp(pred.base$log_processes$log_B - 1.96*pred.base$log_processes$se_log_B),
                            B.UCI = exp(pred.base$log_processes$log_B + 1.96*pred.base$log_processes$se_log_B),
                            R =     pred.base$processes$R,
                            R.LCI = exp(pred.base$log_processes$log_R - 1.96*pred.base$log_processes$se_log_R),
                            R.UCI = exp(pred.base$log_processes$log_R + 1.96*pred.base$log_processes$se_log_R),
                            m =    pred.base$processes$m,
                            m.LCI = exp(pred.base$log_processes$log_m - 1.96*pred.base$log_processes$se_log_m),
                            m.UCI = exp(pred.base$log_processes$log_m + 1.96*pred.base$log_processes$se_log_m),
                            retro.year = "Base model",mod = mod.select, scenario = scenario.select)
}


if(mod.select != "TLM")
{
  # Get the overall estimates + the 95% CI
  base.trends <- data.frame(years = min(years):(2023),
                            B = exp(pred.base$log_tot_frame$log_totB),
                            R = exp(pred.base$log_tot_frame$log_totR),
                            m = pred.base$totals$mean_m,
                            B.LCI = exp(pred.base$log_tot_frame$log_totB - 1.96*pred.base$log_tot_frame$se_log_totB),
                            B.UCI = exp(pred.base$log_tot_frame$log_totB + 1.96*pred.base$log_tot_frame$se_log_totB),
                            R.LCI = exp(pred.base$log_tot_frame$log_totR - 1.96*pred.base$log_tot_frame$se_log_totR),
                            R.UCI = exp(pred.base$log_tot_frame$log_totR + 1.96*pred.base$log_tot_frame$se_log_totR),
                            m.LCI = exp(pred.base$log_tot_frame$log_mean_m - 1.96*pred.base$log_tot_frame$se_log_mean_m),
                            m.UCI = exp(pred.base$log_tot_frame$log_mean_m + 1.96*pred.base$log_tot_frame$se_log_mean_m),
                            retro.year = "Base model",mod = mod.select, scenario = scenario.select)
}

# Remove the projection year for SEAM (does nothing for TLM)
base.trend <- base.trends[!is.na(base.trends$R),]
# combine retro and base
retro.base <- rbind(retro,base.trend)


cols <- c(rep('#005BBB',4),rep('firebrick2',4),rep('darkgrey',4),rep('#FFD500',4))
points <- rep(c(21:24),4) 
b.retro <- ggplot(data=retro.base,aes(x= years, y = B,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + geom_point(size=3) + scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols) +
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Fully recruited biomass (tonnes)")
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                 "/Sab_Biomass_retro.png"),b.retro,base_width = 10,base_height =7)


r.retro <- ggplot(data=retro.base,aes(x= years, y = R,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + geom_point(size=3) + scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols)  + 
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Recruit biomass (tonnes)")
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                 "/Sab_Recruit_retro.png"),r.retro,base_width = 10,base_height =7)


m.retro <- ggplot(data=retro.base,aes(x= years, y = m,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + geom_point(size=3) + scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols)  + 
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Natural mortality (instantaneous)")
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                 "/Sab_mort_retro.png"),m.retro,base_width = 10,base_height =7)

# Calculate mohn's rho, it is simply the Relative bias of the estimate...
# So now bring in the final model run to get the 'true' value for the calculation

bias <- NULL

# Do a 5 year peel as recommend
for(j in retro.years)
{
  retro.dat <- retro %>% dplyr::filter(retro.year==j,years ==j)
  base.dat <- base.trend %>% dplyr::filter(years == j)
  bias[[as.character(j)]] <- data.frame(B = (retro.dat$B - base.dat$B),
                                        rel.B = (retro.dat$B - base.dat$B) / base.dat$B,
                                        R = (retro.dat$R - base.dat$R),
                                        rel.R = (retro.dat$R - base.dat$R) / base.dat$R,
                                        m = (retro.dat$m - base.dat$m),
                                        rel.m = (retro.dat$m - base.dat$m) / base.dat$m,
                                        mod = mod.select)
}


# So now we can use this to calculate mohn's rho
bias <- do.call('rbind',bias)
#Remove 2020 as there is no survey data....
bias <- bias[row.names(bias) != "2020",]
# mohn's rho being fine when it is < 0.2 is found in Hutrtado-Ferro 2015 paper.
# Not we are removing 2020 from this since we don't have any survey data.
mohns.rhos <- data.frame(mr.B = sum(bias$rel.B)/n.retro.years,
                         mr.R = sum(bias$rel.R)/n.retro.years,
                         mr.m = sum(bias$rel.m)/n.retro.years,
                         mr.B5 = sum(bias$rel.B[(nrow(bias)-4):nrow(bias)])/5,
                         mr.R5 = sum(bias$rel.R[(nrow(bias)-4):nrow(bias)])/5,
                         mr.m5 = sum(bias$rel.m[(nrow(bias)-4):nrow(bias)])/5,
                         mod = mod.select)

saveRDS(mohns.rhos,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_mohns_rose_",mod.select,"_",scenario.select,".Rds"))
saveRDS(bias,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_bias_",mod.select,"_",scenario.select,".Rds"))

