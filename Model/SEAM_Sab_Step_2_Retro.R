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
#load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90.Rdata")
load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90RSCS_newMWSH_GBb.RData")

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
qR <- 0.33# This is for TMB (log recruit catchability) testing catchability of 0.5, test 0.3 and 0.1
init.m <- 0.2 # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05R.size <- "75"
FR.size <- "90"
R.size <- "75"
# The different growth models.
g.mod <- 'g_original'
#g.mod <- 'g_1'
#g.mod <- 'alt_g'
#g.mod <- 'proper_g'
# models
mods <- c("SEAM")
n.mods <- length(mods)
retro.years <- 2005:2021
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
  
  sizes <- seq(0.025,2,by=0.05) # So I'd be using the 1.025 bin and everything bigger for the t+1 fully-recruited
  surv.years <- unique(surv.dat$year)
  # The w.yst object is exactly proportional to mod.dat$I, there is an offset, but given I need proportions I think this object is perfectly fine to use.
  # SO this mw.per.bin is taking the stratified biomass and dividing it by the stratified numbers in each bin, which gives us the MW in that bin. 
  # There is probably a MW object out there somewhere with this in it, but it should just be the same thing as this.
  mw.per.bin <- data.frame(mw.per.bin = rbind(survey.obj$Sab$shf.dat$w.yst/survey.obj$Sab$shf.dat$n.yst,rep(NA,40),rep(NA,40)),year = c(surv.years,2015,2020))
  N.per.bin <- data.frame(N.per.bin = rbind(survey.obj$Sab$shf.dat$n.yst,rep(NA,40),rep(NA,40)),year = c(surv.years,2015,2020))
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
  B.prop.per.bin.90 <- N.per.bin[,bin.90.plus]/g.proper$total.abun.90
  B.prop.per.bin.frs <- N.per.bin[,bin.frs.plus]/g.proper$total.abun.frs
  # Recs
  B.prop.per.bin.rec       <- N.per.bin[,bin.rec]/g.proper$total.rec.abun
  B.prop.per.bin.frs.minus <- N.per.bin[,bin.frs.minus]/g.proper$total.frs.minus
  
  # And the average mw in each of the bins of interest, first for the FRs
  g.proper$mw.frs.plus <-  rowSums(mw.per.bin[,bin.frs.plus] * B.prop.per.bin.frs,na.rm=T)
  g.proper$mw.90.plus <-   rowSums(mw.per.bin[,bin.90.plus] * B.prop.per.bin.90,na.rm=T)
  # and for the rec
  g.proper$mw.recs <-      rowSums(mw.per.bin[,bin.rec] * B.prop.per.bin.rec,na.rm=T)
  g.proper$mw.frs.minus <- rowSums(mw.per.bin[,bin.frs.minus] * B.prop.per.bin.frs.minus,na.rm=T)
  
  g.proper$g.proper <- c(g.proper$mw.frs.plus[2:length(g.proper$mw.frs.plus)]/g.proper$mw.90.plus[1:(length(g.proper$mw.90.plus)-1)],NA)
  g.proper$gR.proper<- c(g.proper$mw.frs.minus[2:length(g.proper$mw.frs.minus)]/g.proper$mw.recs[1:(length(g.proper$mw.recs)-1)],NA)
  
  
  g.proper[g.proper$year %in% c(1991,2015,2020),-1] <- NA
  g.proper[g.proper$year %in% c(2014,2019),which(names(g.proper) %in% c("g.proper","gR.proper"))] <- NA
  
  # Fill in the mean for the missing years
  g.proper$g.proper[g.proper$year %in% c(1991,2014,2015,2019,2020,2022)] <- median(g.proper$g.proper,na.rm=T)
  g.proper$gR.proper[g.proper$year %in% c(1991,2014,2015,2019,2020,2022)] <- median(g.proper$gR.proper,na.rm=T)
  
  
  # now need to add in 2015 and 2020 to mod.dat...
  mod.dat.tmp <- mod.dat
  mod.dat.tmp[nrow(mod.dat.tmp)+1:2,] <- NA
  mod.dat.tmp$year[nrow(mod.dat.tmp):(nrow(mod.dat.tmp)-1)] <- c(2015,2020)
  mod.dat.tmp <- mod.dat.tmp[order(mod.dat.tmp$year),]
  
  growth <- data.frame(year = mod.dat.tmp$year,g = mod.dat.tmp$g, gR = mod.dat.tmp$gR,
                       #g.alt = alt.g$alt.g, gR.alt = mod.dat.tmp$gR, # I don't have a good idea how to estiamte gR growth, so using the other way
                       g.proper = g.proper$g.proper,gR.proper = g.proper$gR.proper)
  # Now addin the missing growth years for g and gR
  growth$g[growth$year %in% c(2015,2020)] <- median(growth$g,na.rm=T)
  growth[growth$year %in% c(2015,2020),names(growth) %in% c("gR","gR.alt")] <- median(growth$gR,na.rm=T)
  growth <- growth[which(!is.na(growth$g)),]
  growth[nrow(growth)+1,] <- growth[nrow(growth),]
  growth$year[nrow(growth)] <- max(years) + 1
  
  growth <- growth %>% dplyr::filter(year >= min(years))
  
  #if(g.mod == 'g_original') g <- data.frame(g=growth$g,gR = growth$gR)
  #if(g.mod == 'alt_g') g <- data.frame(g=growth$g.alt,gR = growth$gR.alt)
  g <- data.frame(g=growth$g.proper,gR = growth$gR.proper) # if(g.mod == 'proper_g') 
  #if(g.mod == 'g_1') g <- data.frame(g=growth$g/growth$g,gR = growth$gR/growth$gR)
  
  # Turn this into a vector and add a value for next year.
  #growths <- growth %>% dplyr::filter(year %in% c(years,(max(years)+1)))
  
  
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
      growths <- growth %>% dplyr::filter(year %in% c(years,(max(years)+1)))
      #if(g.mod == 'g_original') g <- data.frame(g=growths$g,gR = growths$gR)
      #if(g.mod == 'alt_g') g <- data.frame(g=growths$g.alt,gR = growths$gR.alt)
      g <- data.frame(g=growths$g.proper,gR = growths$gR.proper) #if(g.mod == 'proper_g') 
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
      
   
      # Subset the fishery data to the correct years. 
      Sab.fish <- Sab.fishs %>% dplyr::filter(survey.year %in% c(years,(max(years)+1))) # 
     
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
        
        set_data<-data_setup(data=mod.input.sf,growths=data.frame(g = g$g,gR = g$gR),catch=as.data.frame(catchy$sum_catches),
                             model="SEBDAM",mesh=Sab.mesh$mesh,obs_mort=T,prior=T,prior_pars=c(20,40),#fix_m = 0.3,
                             mult_qI=T,spat_approach="spde",
                             knot_obj=Sab.mesh$knots,knot_area=pred.grid$area,separate_R_aniso = T,
                             all_se=T,weighted_mean_m = T)

        # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
        set_data$par$log_S <- log(0.072) # I think 0.072 will be the 'right' number
        # #set_data$par$log_R0 <- l.init.R # 5.3 = 200, 5 = 148, 4 = 55, 5.9915 = 400, 4.606 = 100
        set_data$par$log_qR <- log(qR)
        # #set_data$map <-list(log_m0=factor(NA))
        set_data$map <-list(log_S=factor(NA),log_qR = factor(NA))
      } # end if(mod.select != "TLM")

      # A TLM version of the same...
      if(mod.select == "TLM")
      {
        catch.tlm <- catch.sf %>% group_by(Year,.drop=F) %>% dplyr::summarise(catch = sum(Catch,na.rm=T))
        
        set_data<-data_setup(data=as.data.frame(mod.input.sf),growths=data.frame(g = g$g,gR = g$gR),
                             catch=catch.tlm$catch, model="TLM",obs_mort=TRUE,prior=TRUE)
        # So this will fix the mean value of m0 to be whatever the initial value is set at in the data_setup step.  Let's see what happens!
        #set_data<-fix_param(obj=set_data, pars = list(log_q_R=lqr))
        set_data$par$log_q_R <- log(qR) # 
        set_data$map <-list(log_q_R=factor(NA))
      } # end if(mod.select == "TLM")
      
      mod.fit<-fit_model(set_data,silent=F)
      
      # Now save the results appropriately
      if(mod.select != "TLM") 
      {
        scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR,"_",num.knots,"_knots")
        saveRDS(mod.fit,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select,".Rds"))
        saveRDS(Sab.mesh,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select,"_mesh.Rds"))
        saveRDS(pred.grid,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select,"_predict_grid.Rds"))
      }
      
      if(mod.select == "TLM") 
      {
        
        scenario.select <- paste0(min(years),"_",max(years),"_qR_",qR,"_")
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
qR <- 0.33
#init.m <- 0.2 # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05
# g.mod <- 'g_original'
# #g.mod <- 'g_1'
# #g.mod <- 'alt_g'
# #g.mod <- 'proper_g'
# The survey biomass index for 1994 says there were 249 tonnes of recruits that year.
#l.init.R <- log(250) # Going to test 100, 250, and 500.

mod.select <- "SEAM"
retro.years <- c(2005:2014,2016:2019,2021)
n.retro.years <- length(retro.years)

#scenario.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",exp(l.init.m),"_R0_",exp(l.init.R),"_",num.knots,"_knots")
#tst <- readRDS("D:/framework/SFA_25_26_2024/Model/Results/Sab/R_75_FR_90/Retros/Sab_SEAM_model_output_1994_2010_vary_m_m0_0.4_qR_0.5_4_knots.Rds")

# Load the correct base model...

mod.fit <- NULL
retro.trends <- NULL
pred.proc <- NULL
for(j in retro.years)
{
  if(mod.select != "TLM") scenario.select.retro <- paste0(min(years),"_",j,"_qR_",qR,"_",num.knots,"_knots")
  if(mod.select == "TLM") scenario.select.retro <- paste0(min(years),"_",j,"_qR_",qR)

  mod.fit[[j]] <- readRDS(paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_",mod.select,"_model_output_",scenario.select.retro,".Rds"))
  pred.proc[[j]] <- get_processes(mod.fit[[j]])
  
  if(j == min(retro.years))  scenario.select <- paste0(min(years),"_",2022,"_qR_",qR,"_",num.knots,"_knots")
  {

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
    # Remove 2015 and 2020
    base.trend <- base.trend %>% dplyr::filter(!years %in% c(2015,2020))
  } #end if j == min(retro.years)
  
  # Now the retro data.
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

# combine retro and base
retro.base <- rbind(retro,base.trend)


cols <- c(rep('#005BBB',4),rep('firebrick2',4),rep('darkgrey',4),rep('#FFD500',3),'black')
points <- c(rep(21:24,4)) 
b.retro <- ggplot(data=retro.base %>% dplyr::filter(years < 2015),aes(x= years, y = B,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + 
  geom_line(data=retro.base %>% dplyr::filter(years %in% 2016:2019),size=1)+
  geom_point(data=retro.base %>% dplyr::filter(years %in% c(1994:2014,2016:2019,2021)),size=3) + 
  scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols) +
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Fully recruited biomass (tonnes)")
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                 "/Sab_Biomass_retro.png"),b.retro,base_width = 10,base_height =7)


r.retro <- ggplot(data=retro.base %>% dplyr::filter(years < 2015),aes(x= years, y = R,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + 
  geom_line(data=retro.base %>% dplyr::filter(years %in% 2016:2019),size=1) +
  geom_point(data=retro.base %>% dplyr::filter(years %in% c(1994:2014,2016:2019,2021)),size=3) + 
  scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols)  + 
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Recruit biomass (tonnes)")
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                 "/Sab_Recruit_retro.png"),r.retro,base_width = 10,base_height =7)


m.retro <- ggplot(data=retro.base %>% dplyr::filter(years < 2015),aes(x= years, y = m,group=retro.year,color=retro.year,shape = retro.year,fill = retro.year)) +
  geom_line(size=1) + 
  geom_line(data=retro.base %>% dplyr::filter(years %in% 2016:2019),size=1)+
  geom_point(data=retro.base %>% dplyr::filter(years %in% c(1994:2014,2016:2019,2021)),size=3) + 
  scale_shape_manual("",values = points) + 
  scale_color_manual("",values =cols) + scale_fill_manual("",values =cols)  + 
  scale_x_continuous(breaks = seq(1990,2030,by=5)) + xlab("") + 
  ylab("Natural mortality (instantaneous)")
save_plot(paste0(repo.loc,"Figures/Sab/R_",R.size,"_FR_",FR.size,"/Retros/",mod.select,"_",scenario.select,
                 "/Sab_mort_retro.png"),m.retro,base_width = 10,base_height =7)

# Calculate mohn's rho, it is simply the Relative bias of the estimate...
# So now bring in the final model run to get the 'true' value for the calculation

bias <- NULL

# Do a 5 year peel as recommended
#act.retro.years <- c(2010:2014,2016:2019,2021)
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
# mohn's rho being fine when it is < 0.2 is found in Hutrtado-Ferro 2015 paper.
# Note we are removing 2020 from this since we don't have any survey data.
mohns.rhos <- data.frame(mr.B = sum(bias$rel.B)/n.retro.years,
                         mr.R = sum(bias$rel.R)/n.retro.years,
                         mr.m = sum(bias$rel.m)/n.retro.years,
                         mr.B5 = sum(bias$rel.B[(nrow(bias)-4):nrow(bias)])/5,
                         mr.R5 = sum(bias$rel.R[(nrow(bias)-4):nrow(bias)])/5,
                         mr.m5 = sum(bias$rel.m[(nrow(bias)-4):nrow(bias)])/5,
                         mod = mod.select)

saveRDS(mohns.rhos,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_mohns_rose_",mod.select,"_",scenario.select,".Rds"))
saveRDS(bias,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_bias_",mod.select,"_",scenario.select,".Rds"))
saveRDS(retro.base,paste0(repo.loc,"Results/Sab/R_",R.size,"_FR_",FR.size,"/Retros/Sab_retro_",mod.select,"_",scenario.select,".Rds"))
#tst <- readRDS("D:/Framework/SFA_25_26_2024/Model/Results/Sab/R_75_FR_90/Retros/archive/Sab_bias_SEAM_1994_2022_vary_m_m0_2_qR_0.5_10_knots.Rds")
