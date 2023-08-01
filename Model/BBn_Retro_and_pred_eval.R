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
load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90.Rdata")
#load("D:/Framework/SFA_25_26_2024/Model/Data/testing_results_framework3.Rdata")
#surv.dat <- surv.dat$BBn
#saveRDS(surv.dat,'D:/Github/BBn_model/Results/BBn_surv.dat.RDS')
#surv.dat <- readRDS('D:/Framework/SFA_25_26_2024/Model/Data/BBn_surv.dat.RDS')
surv.dat <- surv.dat$BBn
mod.dat <- survey.obj$BBn$model.dat
bbn.fish <- readRDS('D:/Framework/SFA_25_26_2024/Model/Data/BBn_fish.dat.RDS')
bbn.fish$pro.repwt <- bbn.fish$pro.repwt/1000
#### Finshed Data prep and clean up!
###############################################################################################



# Set parameters for the run...
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
atow<-800*2.4384/10^6 # area of standard tow in km2
years <- 1994:2022
NY <- length(years)
R.size <- "75"
FR.size <- "90"
num.knots <- 10 # Going to test 10, 15, and 20
lqr <- log(0.5)# This is for TMB (log recruit catchability) testing catchability of 0.5, test 0.3 and 0.1
l.init.m <- log(2) # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05
# The survey biomass index for 1994 says there were 249 tonnes of recruits that year.
#l.init.R <- log(250) # Going to test 100, 250, and 500.

mods <- c("SEAM")
n.mods <- length(mods)
retro.years <- 2010:2021
n.retro.years <- length(retro.years)

for(i in 1:n.mods)
{
  mod.select <- mods[i]
  # Bring in the data, only need to do that once...
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
  
  
  growth <- data.frame(g = mod.dat$g, gR = mod.dat$gR,year = mod.dat$year)
  # need to fudge in 2020 data.. taking average of 2019 and 2021
  growth[nrow(growth)+1,] <- growth[nrow(growth),]
  growth$year[nrow(growth)] <- 2020
  growth$g[growth$year == 2020] <- mean(growth$g[growth$year %in% c(2019,2021)])
  growth$gR[growth$year == 2020] <- mean(growth$gR[growth$year %in% c(2019,2021)])
  # Now reorder by year...
  growth <- growth[order(growth$year),]
  growths <- growth %>% dplyr::filter(year %in% c(years,(max(years)+1)))
  
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
 
  bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] <- bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] -1
  bbn.fish.sfs <- st_as_sf(bbn.fish,coords = c("lon","lat"),remove =F, crs = 4326)
  bbn.fish.sfs <- bbn.fish.sfs %>% st_transform(crs= 32619)
  # Now lets clip this to be data inside of our bbn boundary.
  bbn.fish.sfs <- st_intersection(bbn.fish.sfs,bbn.shape)
  if(mod.select != "TLM")
  {
    bbn.mesh <- setup_mesh(mod.input.sf,model_bound = bbn.shape,nknot=num.knots, max.edge = c(3,10),cutoff=2,seed=66) # Seeds 20 and 66 work
    bbn.mesh.sf <- inla.mesh2sf(bbn.mesh$mesh)
    bbn.mesh.sf$triangles$geometry <- bbn.mesh.sf$triangles$geometry*1000
    bbn.mesh.sf$vertices$geometry <- bbn.mesh.sf$vertices$geometry*1000
    st_crs(bbn.mesh.sf$triangles) <- 32619
    # Now make the prediction grid
    pred.grid<-setup_pred_grid(knots=bbn.mesh$knots,model_bound=bbn.mesh$utm_bound)
  } # end if(mod.select != "TLM")
  
  for(j in 1:n.retro.years)
  {
    years <- 1994:retro.years[j]
    NY <- length(years)

    # OK, so step 1 here is getting the model input that Raphael needs for the model
    # The survey data....
    mod.input.sfs <- mod.input.sf
    mod.input.sfs$Year <- mod.input.sfs$year - (min(years)-1)
    mod.input.sfs <- mod.input.sfs %>% dplyr::filter(year %in% years)
    growths <- growth %>% dplyr::filter(year %in% c(years,(max(years)+1)))
    # Now we can clip both of these to subset it to the data that I think we need for the analysis....
    # Subset the fishery data as necessary
    bbn.fish.sf <- bbn.fish.sfs %>% dplyr::filter(survey.year %in% years)
    # OK, so now let's see if we can use the catch knot thing Raph made to split this up withing the BBn domain
    #We just need 3 columns for this
    catch.sf <- bbn.fish.sf %>% dplyr::select(pro.repwt,survey.year)
    names(catch.sf) <- c("Catch","Year","geometry")
    # For the moment we need to have this starting at year 1.
    catch.sf$Year <-  catch.sf$Year - (min(years)-1)

    # The SEBDAM set data.
    if(mod.select != "TLM")
    {
      catchy <- catch_spread(catch = catch.sf,knots = bbn.mesh$knots)
      catchy$sum_catches[,ncol(catchy$sum_catches)+1] <- 0
      set_data<-data_setup(data=mod.input.sfs,growths=growths[,1:2],catch=catchy$sum_catches[],
                           model="SEBDAM",mesh=bbn.mesh$mesh,obs_mort=T,prior=T,prior_pars=c(10,12),#fix_m = 0.3,
                           mult_qI=T,spat_approach="spde",
                           knot_obj=bbn.mesh$knots,knot_area=pred.grid$area,separate_R_aniso = T,
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
      #m0.par <- exp(set_data$par$log_m0)
      scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",exp(l.init.m),"_qR_",exp(lqr),"_",num.knots,"_knots")
      saveRDS(mod.fit,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
      saveRDS(bbn.mesh,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
      saveRDS(pred.grid,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
      
    }
    
    if(mod.select == "TLM") 
    {
      scenario.select <- paste0(min(years),"_",max(years),"_qR_",exp(lqr))
      saveRDS(mod.fit,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
      
    }

}} # End the i and j loops



###################  SECTION 2 Make the Retro plots ##############################  SECTION 2 Make the Retro plots ###########
# Now make the retrospective plots...

# Set parameters for the run...
repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
years <- 1994:2022
R.size <- "75"
FR.size <- "90"
qR <- "0_5"# This is for TMB (log recruit catchability) testing catchability of 0.5, test 0.3 and 0.1
num.knots <- 10 # Going to test 10, 15, and 20
lqr <- log(0.5)# This is for TMB (log recruit catchability) testing catchability of 0.5, test 0.3 and 0.1
l.init.m <- log(2) 

mod.select <- "SEAM"
retro.years <- 2010:2021
n.retro.years <- length(retro.years)

#scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",exp(l.init.m),"_R0_",exp(l.init.R),"_",num.knots,"_knots")
#tst <- readRDS("D:/framework/SFA_25_26_2024/Model/Results/BBn/R_75_FR_90/Retros/BBn_SEAM_model_output_1994_2010_vary_m_m0_0.15_R0_250_15_knots.Rds")

mod.fit <- NULL
retro.trends <- NULL
pred.proc <- NULL
for(j in retro.years)
{
  if(mod.select != "TLM") scenario.select <- paste0(min(years),"_",j,"_vary_m_m0_",exp(l.init.m),"_qR_",exp(lqr),"_",num.knots,"_knots")
  if(mod.select == "TLM") scenario.select <- paste0(min(years),"_",j,"_qR_",exp(lqr))
  mod.fit[[j]] <- readRDS(paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
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

# So now bring in the final model run to get the 'true' value for the calculation
if(mod.select != "TLM") scenario.select <- paste0(min(years),"_",2022,"_vary_m_m0_",exp(l.init.m),"_qR_",exp(lqr),"_",num.knots,"_knots")
if(mod.select == "TLM") scenario.select <- paste0(min(years),"_",2022,"_qR_",qR)
base.mod <- readRDS(paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/BBn_",mod.select,"_model_output_",scenario.select,".Rds"))
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
                                #geom_line(data=base.trend,aes(x=years, y = B,color=retro.year),size=1) +
                                scale_color_manual("",values =cols) + scale_fill_manual("",values =cols) +
                                scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
                                ylab("Fully recruited biomass (tonnes)")
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                           "/BBn_Biomass_retro.png"),b.retro,base_width = 10,base_height =7)

  
r.retro <- ggplot(data=retro.base,aes(x= years, y = R,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
                                  geom_line(size=1) + geom_point(size=3) + scale_shape_manual("",values = points) + 
                                  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols)  + 
                                  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
                                  ylab("Recruit biomass (tonnes)")
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                 "/BBn_Recruit_retro.png"),r.retro,base_width = 10,base_height =7)


m.retro <- ggplot(data=retro.base,aes(x= years, y = m,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
                                  geom_line(size=1) + geom_point(size=3) + scale_shape_manual("",values = points) + 
                                  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols)  + 
                                  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
                                  ylab("Natural mortality (instantaneous)")
save_plot(paste0(repo.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                 "/BBn_mort_retro.png"),m.retro,base_width = 10,base_height =7)

# Calculate mohn's rho, it is simply the Relative bias of the estimate...

bias <- NULL

# Doing more than the recommedned 5 year peel 
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


saveRDS(mohns.rhos,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_mohns_rose_",mod.select,"_",scenario.select,".Rds"))
saveRDS(bias,paste0(repo.loc,"Results/BBn/R_",R.size,"_FR_",FR.size,"/Retros/BBn_bias_",mod.select,"_",scenario.select,".Rds"))


