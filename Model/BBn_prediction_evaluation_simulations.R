# Here are the BBn and Sable Decision Tables.
library(lubridate)
library(SEBDAM)
library(terra)
library(tidyterra)
library(ggthemes)
library(tidyverse)
library(sf)
library(cowplot)


theme_set(theme_few(base_size = 22))
u.colors <- c("#FFD500","#005BBB")
  
# Here's where I'm sticking the figures
mod.loc <- "D:/Framework/SFA_25_26_2024/Model/"
fun.loc <-"D:/Github/Framework/Model/functions/"
rp.loc <- "D:/Framework/SFA_25_26_2024/RefPts/"

# Functions we need to source for this to all work
source(paste0(fun.loc,"Decision_Table_function.R"))

# Recruit and FR sizes
R.size <- 75
FR.size <- 90
qR <- 0.33
init.m <- 0.2
g.mod <- 'g_original'
#g.mod <- "g_1"
#g.mod <- 'alt_g'
#g.mod <- 'proper_g'
num.knots <- 20

bbn.shape <- st_read("D:/Github/GIS_layers/survey_boundaries/BBn.shp", quiet=T)
bbn.shape <- bbn.shape %>% st_transform(crs = 32619) # BBn is right on the UTM boundary of 19 amd 20
# We need to get the landings for the previous year after the survey (in this case the Landings from June-December 2022)

bbn.fish <- readRDS(paste0(mod.loc,"/Data/BBn_fish.dat.RDS"))
bbn.fish$pro.repwt <- bbn.fish$pro.repwt/1000
bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] <- bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] -1
# Clip to the survey area....
bbn.fish.sf <- st_as_sf(bbn.fish,coords = c("lon","lat"),remove =F, crs = 4326)
bbn.fish.sf <- bbn.fish.sf %>% st_transform(crs= 32619)

# Now lets clip this to be data inside of our Sab boundary.
bbn.fish.sf <- st_intersection(bbn.fish.sf,bbn.shape)
bbn.fish.survey.year <- bbn.fish.sf %>% dplyr::group_by(survey.year) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))


# Inputs...
pe.years <- c(2006:2019,2021:2022) # Exclude 2020 from this (could exclude 2021 but I occasionally want the 2021 real data for plots afterwards)
n.pe.years <- length(pe.years)

# Reference Points
n.sims <- 1e6
preds <- NULL
reals <- NULL
diffs <- NULL
for(i in 1:n.pe.years)
{

# The removals for upcoming year.
rem <- bbn.fish.survey.year$tot[bbn.fish.survey.year$survey.year == pe.years[i]]
  
# You want the PE.year -1 model here as that's what we'd use to do the prediction
pred.select <- paste0("1994_",pe.years[i]-1,"_vary_m_m0_",init.m,"_qR_",qR,"_",num.knots,"_knots_",g.mod,"_vary_q=TRUE")
# Then you want the real model to compare your predictions too
realized.select <- paste0("1994_",pe.years[i],"_vary_m_m0_",init.m,"_qR_",qR,"_",num.knots,"_knots_",g.mod,"_vary_q=TRUE")
# For i =1 this is the 2010 model whose parameters we can use to get a prediction for 2011 biomass
pred.mod <- readRDS(paste0(mod.loc,"Results/BBn/R_75_FR_90/Retros/BBn_SEAM_model_output_",pred.select,".Rds"))
# This is then the 2011 model that gives us the realized value.
if(pe.years[i] <  2022) real.mod <- readRDS(paste0(mod.loc,"Results/BBn/R_75_FR_90/Retros/BBn_SEAM_model_output_",realized.select,".Rds"))
if(pe.years[i] == 2022) real.mod <- readRDS(paste0(mod.loc,"Results/BBn/R_75_FR_90/BBn_SEAM_model_output_",realized.select,".Rds"))

# Here we get the 2011 realized estiamtes 
real.res <- get_processes(real.mod)
proj.res <- get_processes(pred.mod) # and here we can get the 2010 projections for 2011
# The second last value here is the actual Biomass in the year of interest (the final one is the projection)
log.real.B <- real.res$log_tot_frame$log_totB[length(real.res$log_tot_frame$log_totB)-1]
se.log.real.B <- real.res$log_tot_frame$se_log_totB[length(real.res$log_tot_frame$se_log_totB)-1]
real.B <- exp(log.real.B)
LCI.real.B <- exp(log.real.B - 2*se.log.real.B)
UCI.real.B <- exp(log.real.B + 2*se.log.real.B)
               
# We can also extract Raph's Projected year and compare what that looks like, these have WAY more uncertainty in them.
# Note that model is run with no removals, 
log.proj.B <- proj.res$log_tot_frame$log_totB[length(proj.res$log_tot_frame$log_totB)]
se.log.proj.B <- proj.res$log_tot_frame$se_log_totB[length(proj.res$log_tot_frame$se_log_totB)]
proj.B <- exp(log.proj.B) - rem
LCI.proj.B <- exp(log.proj.B - 2*se.log.proj.B) - rem
UCI.proj.B <- exp(log.proj.B + 2*se.log.proj.B) - rem


  
# Fully productivity scenario
full     <- dec.tab(mod.select = "SEAM",data = pred.mod, catch.scenarios = rem,PSL = 0,
                    n.sims = n.sims,g.adj=1,gR.adj = 1,m.adj=1,r.adj=1)
# 0 productivity scenario
zero     <- dec.tab(mod.select = "SEAM",data = pred.mod, catch.scenarios = rem,PSL = 0,
                    n.sims = n.sims,g.adj=0,gR.adj = 0,m.adj=0,r.adj=0)
# No growth scenario
ng       <- dec.tab(mod.select = "SEAM",data = pred.mod, catch.scenarios = rem,PSL = 0,
                    n.sims = n.sims,g.adj=0,gR.adj = 0,m.adj=1,r.adj=1)
# No FR growth
ngfr     <- dec.tab(mod.select = "SEAM",data = pred.mod, catch.scenarios = rem,PSL = 0,
                    n.sims = n.sims,g.adj=0,gR.adj = 1,m.adj=1,r.adj=1)
# Average mortality scenario
med.m    <- dec.tab(mod.select = "SEAM",data = pred.mod, catch.scenarios = rem,PSL = 0,
                    n.sims = n.sims,g.adj=1,gR.adj = 1,m.adj='avg',r.adj=1)
# Average recruitment scenario
med.r    <- dec.tab(mod.select = "SEAM",data = pred.mod, catch.scenarios = rem,PSL = 0,
                    n.sims = n.sims,g.adj=1,gR.adj = 1,m.adj=1,r.adj='avg')
# Average everything scenario
med.all  <- dec.tab(mod.select = "SEAM",data = pred.mod, catch.scenarios = rem,PSL = 0,
                    n.sims = n.sims,g.adj='avg',gR.adj = 'avg',m.adj='avg',r.adj='avg')

tmp.pred <- data.frame(B=c(median(full$Biomass),median(zero$Biomass),median(ng$Biomass),
                                                     median(ngfr$Biomass),median(med.m$Biomass),median(med.r$Biomass),
                                                     median(med.all$Biomass),proj.B),
                                                 B.L95 = c(quantile(full$Biomass,probs=0.025),quantile(zero$Biomass,probs=0.025),quantile(ng$Biomass,probs=0.025),
                                                           quantile(ngfr$Biomass,probs=0.025),quantile(med.m$Biomass,probs=0.025),
                                                           quantile(med.r$Biomass,probs=0.025),quantile(med.all$Biomass,probs=0.025),LCI.proj.B),
                                                 B.U95 = c(quantile(full$Biomass,probs=0.975),quantile(zero$Biomass,probs=0.975),quantile(ng$Biomass,probs=0.975),
                                                           quantile(ngfr$Biomass,probs=0.975),quantile(med.m$Biomass,probs=0.975),
                                                           quantile(med.r$Biomass,probs=0.975),quantile(med.all$Biomass,probs=0.975),UCI.proj.B),
                                                 year = pe.years[i],
                                                 Scenario=as.factor(c("Previous year","Zero Productivity","No Growth","No FR Growth",
                                                                      "Median m","Median R","Median Productivity","Spatial Model")))
preds[[as.character(pe.years[i])]] <- tmp.pred
tmp.reals <- data.frame(B=real.B,B.L95 = LCI.real.B,B.U95 = UCI.real.B, year = pe.years[i],Scenario = "Realized")
reals[[as.character(pe.years[i])]] <- tmp.reals

diffs[[as.character(pe.years[i])]] <- data.frame(B.diff = tmp.pred$B - tmp.reals$B, B.diff.abs = abs(tmp.pred$B - tmp.reals$B),
                                                 Per.B.diff = 100*((tmp.pred$B - tmp.reals$B)/tmp.reals$B),
                                                 year = pe.years[i], Scenario= tmp.pred$Scenario)
} # end for(i in 1:n.pe.years)


B.diff <- do.call("rbind",diffs)
# We should remove 2020 and 2021 because those are real given the missing 2020 data.
B.diff <- B.diff %>% dplyr::filter(!year %in% c(2020:2021))
# The actual model
real.mods <- do.call('rbind',reals)
# All the predictions
pred <- do.call('rbind',preds)

# A summary of B.diff in a csv table will probably be useful
B.diff.summary <- data.frame(B.diff %>% dplyr::group_by(Scenario) %>% dplyr::summarise(med.B.diff = median(B.diff,na.rm=T),
                                                                                       max.B.diff = max(B.diff,na.rm=T),
                                                                                       min.B.diff = min(B.diff,na.rm=T),
                                                                                       med.abs.B.diff = median(B.diff.abs,na.rm=T),
                                                                                       max.abs.B.diff = max(B.diff.abs,na.rm=T),
                                                                                       min.abs.B.diff = min(B.diff.abs,na.rm=T),
                                                                                       per.B.diff = median(Per.B.diff,na.rm=T),
                                                                                       max.per.B.diff = max(Per.B.diff,na.rm=T),
                                                                                       min.per.B.diff = min(Per.B.diff,na.rm=T)
                                                                                       )) 
B.diff.summary[,2:ncol(B.diff.summary)] <- round(B.diff.summary[,2:ncol(B.diff.summary)],digits=0)
if(!dir.exists(paste0(mod.loc,"Results/BBn/R_75_FR_90/Pred_eval/"))) dir.create(paste0(mod.loc,"Results/BBn/R_75_FR_90/Pred_eval/"))
if(!dir.exists(paste0(mod.loc,"Results/BBn/R_75_FR_90/Pred_eval/SEAM_",realized.select))) dir.create(paste0(mod.loc,"Results/BBn/R_75_FR_90/Pred_eval/SEAM_",realized.select))

write.csv(B.diff.summary,paste0(mod.loc,"Results/BBn/R_75_FR_90/Pred_eval/SEAM_",realized.select,"/PE_table.csv"))
# Also make this a nicely formated table for both word and pdf import for later on.
tab <- kableExtra::kbl(B.diff.summary, booktabs = TRUE, escape =F, format = 'pipe')#,
# Save all these prediction evaluation results...
saveRDS(tab,paste0(mod.loc,"Results/BBn/R_75_FR_90/Pred_eval/SEAM_",realized.select,"/PE_word_table.Rds"))
tab.pdf <- kableExtra::kbl(B.diff.summary, booktabs = TRUE, escape =F, format='latex')#,
saveRDS(tab.pdf,paste0(mod.loc,"Results/BBn/R_75_FR_90/Pred_eval/SEAM_",realized.select,"/PE_pdf_table.Rds"))
saveRDS(B.diff,paste0(mod.loc,"Results/BBn/R_75_FR_90/Pred_eval/SEAM_",realized.select,"/Biomass_difference.Rds"))

if(!dir.exists(paste0(mod.loc,"Figures/BBn/R_75_FR_90/Pred_eval/"))) dir.create(paste0(mod.loc,"Figures/BBn/R_75_FR_90/Pred_eval/"))
if(!dir.exists(paste0(mod.loc,"Figures/BBn/R_75_FR_90/Pred_eval/SEAM_",realized.select))) dir.create(paste0(mod.loc,"Figures/BBn/R_75_FR_90/Pred_eval/SEAM_",realized.select))


liner <- c(1,2,1,1,2,1,2,2)
shaper <- rep(c(21,23,22,24),2)
pe.tonnes.ts <- ggplot() + geom_line(data = B.diff %>% dplyr::filter(year %in% 2006:2019), 
                                     aes(x=year,y=B.diff,group = Scenario,color=Scenario,linetype=Scenario),size=1.25) + 
  geom_point(data = B.diff,aes(x=year,y=B.diff,group = Scenario,color=Scenario,fill=Scenario,shape =Scenario),size=2) + 
  scale_color_viridis_d() + scale_fill_viridis_d() + 
  scale_shape_manual(values = shaper) + scale_linetype_manual(values=liner) + 
  scale_x_continuous(breaks = seq(1990,2024,by=2))+
  xlab("") + ylab("Predicted - Realized Biomass (tonnes)")
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/BBn_PE_biomass_ts.png"),
          pe.tonnes.ts,base_width = 11,base_height = 8.5)


pe.per.ts <- ggplot() + geom_line(data = B.diff %>% dplyr::filter(year %in% 2006:2019), 
                                  aes(x=year,y=Per.B.diff,group = Scenario,color=Scenario,linetype=Scenario),size=1.25) +
  geom_point(data = B.diff ,aes(x=year,y=Per.B.diff,group = Scenario,color=Scenario,fill=Scenario,shape=Scenario),size=2) + 
  scale_shape_manual(values = shaper) + scale_linetype_manual(values=liner) + 
  scale_color_viridis_d() + scale_fill_viridis_d() + scale_x_continuous(breaks = seq(1990,2024,by=2)) +
  xlab("") + ylab("Predicted - Realized Biomass (%)")
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/BBn_PE_per_ts.png"),
          pe.per.ts,base_width = 11,base_height = 8.5)



pe.biomass.abs <- ggplot(B.diff) + geom_violin(aes(x=Scenario,y=B.diff.abs),draw_quantiles = c(0.5)) + 
  ylab("|Predicted - Realized Biomass| (tonnes)")+ geom_hline(yintercept=0,color='#005BBB',linetype='dashed',size=1.5) +
  scale_y_continuous(breaks=c(seq(0,1e4,by=500)))
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/BBn_PE_abs_diff.png"),
          pe.biomass.abs,base_width = 17,base_height = 8.5)


pe.biomass.diff <- ggplot(B.diff) + geom_violin(aes(x=Scenario,y=B.diff),draw_quantiles = c(0.5)) + 
  ylab("Predicted - Realized Biomass (tonnes)")+ geom_hline(yintercept=0,color='#005BBB',linetype='dashed',size=1.5) +
  scale_y_continuous(breaks=c(seq(-5000,1e4,by=500)))
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/BBn_PE_diff.png"),
          pe.biomass.diff,base_width = 17,base_height = 8.5)


pe.biomass.per <- ggplot(B.diff) + geom_violin(aes(x=Scenario,y=Per.B.diff),draw_quantiles = c(0.5)) + 
  ylab("Predicted - Realized Biomass (%)") + geom_hline(yintercept=0,color='#005BBB',linetype='dashed',size=1.5) +
  scale_y_continuous(breaks=c(seq(-100,200,by=10)))
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/BBn_PE_per_diff.png"),
          pe.biomass.per,base_width = 17,base_height = 8.5)

# Plots for presentations
real.mods.bf.2020 <- real.mods %>% dplyr::filter(year < 2020)
real.mods.af.2020 <- real.mods %>% dplyr::filter(year > 2020)
# The previous year projection method
py.proj <- pred %>% dplyr::filter(Scenario == "Previous year")
py.proj.bf.2020 <- py.proj %>% dplyr::filter(year < 2020)
py.proj.af.2020 <- py.proj %>% dplyr::filter(year > 2021)
# No FR Growth
nfrg.proj <- pred %>% dplyr::filter(Scenario == "No FR Growth")
nfrg.proj.bf.2020 <- nfrg.proj %>% dplyr::filter(year < 2020)
nfrg.proj.af.2020 <- nfrg.proj %>% dplyr::filter(year > 2021)
# Surplus Production = 0
sp0.proj <- pred %>% dplyr::filter(Scenario == "Productivity Zero")
sp0.proj.bf.2020 <- sp0.proj %>% dplyr::filter(year < 2020)
sp0.proj.af.2020 <- sp0.proj %>% dplyr::filter(year > 2021)
# No Growth
ng.proj <- pred %>% dplyr::filter(Scenario == "No Growth")
ng.proj.bf.2020 <- ng.proj %>% dplyr::filter(year < 2020)
ng.proj.af.2020 <- ng.proj %>% dplyr::filter(year > 2021)


pe.real.ts <- ggplot() + geom_line(data=real.mods.bf.2020,aes(x=year,y=B),color=u.colors[2],size=1.25) + 
                            geom_ribbon(data=real.mods.bf.2020, aes(x=year,ymin=B.L95,ymax=B.U95),fill=u.colors[2],color=u.colors[2],alpha=0.3) +
                            geom_line(data=real.mods.af.2020,aes(x=year,y=B),color=u.colors[2],size=1.25) + 
                            geom_ribbon(data=real.mods.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95),fill=u.colors[2],color=u.colors[2],alpha=0.3) +
                            scale_y_continuous(name = "Fully Recruited Biomass (tonnes)",breaks =seq(0,2e4,by=1e3)) + 
                            scale_x_continuous(name='',breaks = 1990:2030) 
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/PE_real_biomass_ts.png"),
          pe.real.ts,base_width = 11,base_height = 6.5)

# Combine the data above to make a plot with a legend...
combo.py.bf.2020 <- rbind(real.mods.bf.2020,py.proj.bf.2020)
#combo.py.af.2020 <- rbind(real.mods.af.2020,py.proj.af.2020)
# Now add the projection using the decision table projection method we have been using, i.e. using the producitivy parameters from last year

pe.proj.py.ts <-  ggplot() + geom_line(data=combo.py.bf.2020,aes(x=year,y=B,color=Scenario),size=1.25) + 
                             geom_ribbon(data=combo.py.bf.2020, aes(x=year,ymin=B.L95,ymax=B.U95,fill=Scenario),alpha=0.3) +
                             geom_line(data=real.mods.af.2020,aes(x=year,y=B,color=Scenario),size=1.25) + 
                             geom_ribbon(data=real.mods.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95,fill=Scenario),alpha=0.3) +
                             geom_point(data=py.proj.af.2020,aes(x=year,y=B,color=Scenario,fill=Scenario),size=3,shape=22) + 
                             geom_errorbar(data=py.proj.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95,color=Scenario),alpha=0.3,size=4,width=0) +
                             scale_y_continuous(name = "Fully Recruited Biomass (tonnes)",breaks =seq(0,2e4,by=1e3)) + 
                             scale_color_manual(values = u.colors) + scale_fill_manual(values = u.colors) + 
                             scale_x_continuous(name='',breaks = seq(1990,2030,by=2)) 
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/PE_dt_comp_ts.png"),
          pe.proj.py.ts,base_width = 11,base_height = 6.5)

# Combine the data above to make a plot with a legend...
combo.nfrg.bf.2020 <- rbind(real.mods.bf.2020,nfrg.proj.bf.2020)
#combo.nfrg.af.2020 <- rbind(real.mods.af.2020,nfrg.proj.af.2020)
# Now add the projection using the decision table projection method we have been using, i.e. using the producitivy parameters from last year
pe.proj.nfrg.ts <-  ggplot() + geom_line(data=combo.nfrg.bf.2020,aes(x=year,y=B,color=Scenario),size=1.25) + 
                               geom_ribbon(data=combo.nfrg.bf.2020, aes(x=year,ymin=B.L95,ymax=B.U95,fill=Scenario),alpha=0.3) +
                               geom_line(data=real.mods.af.2020,aes(x=year,y=B,color=Scenario),size=1.25) + 
                               geom_ribbon(data=real.mods.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95,fill=Scenario),alpha=0.3) +
                               geom_point(data=nfrg.proj.af.2020,aes(x=year,y=B,color=Scenario,fill=Scenario),size=3,shape=22) + 
                               geom_errorbar(data=nfrg.proj.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95,color=Scenario),alpha=0.3,size=4,width=0) +
                               scale_y_continuous(name = "Fully Recruited Biomass (tonnes)",breaks =seq(0,2e4,by=1e3)) + 
                               scale_color_manual(values = u.colors) + scale_fill_manual(values = u.colors) + 
                               scale_x_continuous(name='',breaks = seq(1990,2030,by=2)) 
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/PE_nfrg_comp_ts.png"),
          pe.proj.nfrg.ts,base_width = 11,base_height = 6.5)


# Combine the data above to make a plot with a legend...
combo.sp0.bf.2020 <- rbind(real.mods.bf.2020,sp0.proj.bf.2020)
#combo.sp0.af.2020 <- rbind(real.mods.af.2020,sp0.proj.af.2020)
# Now add the projection using the decision table projection method we have been using, i.e. using the producitivy parameters from last year
pe.proj.sp0.ts <-  ggplot() + geom_line(data=combo.sp0.bf.2020,aes(x=year,y=B,color=Scenario),size=1.25) + 
                              geom_ribbon(data=combo.sp0.bf.2020, aes(x=year,ymin=B.L95,ymax=B.U95,fill=Scenario,color=Scenario),alpha=0.3) +
                              geom_line(data=real.mods.af.2020,aes(x=year,y=B,color=Scenario),size=1.25) + 
                              geom_ribbon(data=real.mods.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95,fill=Scenario),alpha=0.3) +
                              geom_point(data=sp0.proj.af.2020,aes(x=year,y=B,color=Scenario,fill=Scenario),size=3,shape=22) + 
                              geom_errorbar(data=sp0.proj.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95,color=Scenario),alpha=0.3,size=4,width=0) +
                              scale_y_continuous(name = "Fully Recruited Biomass (tonnes)",breaks =seq(0,2e4,by=1e3)) + 
                              scale_color_manual(values = u.colors) + scale_fill_manual(values = u.colors) +                   
                              scale_x_continuous(name='',breaks = seq(1990,2030,by=2)) 
save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/PE_sp0_comp_ts.png"),
          pe.proj.sp0.ts,base_width = 11,base_height = 6.5)



# Combine the data above to make a plot with a legend...
combo.ng.bf.2020 <- rbind(real.mods.bf.2020,ng.proj.bf.2020)
#combo.ng.af.2020 <- rbind(real.mods.af.2015,ng.proj.af.2020)
#combo.sp0.af.2020 <- rbind(real.mods.af.2020,sp0.proj.af.2020)
# Now add the projection using the decision table projection method we have been using, i.e. using the producitivy parameters from last year
pe.proj.ng.ts <-  ggplot()+ geom_line(data=combo.ng.bf.2020,aes(x=year,y=B,color=Scenario),size=1.25) + 
                            geom_ribbon(data=combo.ng.bf.2020, aes(x=year,ymin=B.L95,ymax=B.U95,fill=Scenario,color=Scenario),alpha=0.3) +
                            geom_line(data=real.mods.af.2020,aes(x=year,y=B,color=Scenario),size=1.25) + 
                            geom_ribbon(data=real.mods.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95,fill=Scenario),alpha=0.3) +
                            geom_point(data=ng.proj.af.2020,aes(x=year,y=B,color=Scenario,fill=Scenario),size=3,shape=22) + 
                            geom_errorbar(data=ng.proj.af.2020, aes(x=year,ymin=B.L95,ymax=B.U95,color=Scenario),alpha=0.3,size=4,width=0) +
                            scale_y_continuous(name = "Fully Recruited Biomass (tonnes)",breaks =seq(0,2e4,by=1e3)) + 
                            scale_color_manual(values = u.colors) + scale_fill_manual(values = u.colors) +                          
                            scale_x_continuous(name='',breaks = seq(1990,2030,by=2)) 

save_plot(paste0(mod.loc,"Figures/BBn/R_",R.size,"_FR_",FR.size,"/Pred_eval/SEAM_",realized.select,"/PE_ng_comp_ts.png"),
          pe.proj.ng.ts,base_width = 11,base_height = 6.5)