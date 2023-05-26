# Here we get the BBn data for use elsewhere, this is much simpler that Sable because we are just going to take the results from the 2023 update

library(sf)
library(sp)
library(tidyverse)
library(ggthemes)
library(cowplot)

theme_set(theme_few(base_size = 22))

## For BBN we will just load the model results from the 2023 assessment once available.  For the moment we'll use the 2022 results....

# OK, so happy with model lets load up some stuff
load("D:/testing_folder/data/Model_testing_results_mixed.RData")
DD.plt <- mod.out$BBn$BUGSoutput
yrs <- yrs$BBn

# DD.plt$median$B[which(yrs %in% c(2020))] <- NA
# DD.plt$sims.list$B[,which(yrs %in% c(2020))] <- NA
# DD.plt$median$R[which(yrs %in% c(2020))] <- NA
# DD.plt$sims.list$R[,which(yrs %in% c(2020))] <- NA
# DD.plt$data$I[which(yrs %in% c(2020))] <- NA
# DD.plt$data$IR[which(yrs %in% c(2020))] <- NA

# Get all the big productivity model output together

res.gg <- data.frame(med = c(apply(DD.plt$sims.list$B, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T)),
                             apply(DD.plt$sims.list$R, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T)),
                             apply(DD.plt$sims.list$m, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T)),
                             apply(DD.plt$sims.list$mR, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T)),
                             apply(DD.plt$sims.list$mu, 2, FUN = function(x) quantile(x,probs=0.5,na.rm=T))),
                     uci = c(apply(DD.plt$sims.list$B, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T)),
                             apply(DD.plt$sims.list$R, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T)),
                             apply(DD.plt$sims.list$m, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T)),
                             apply(DD.plt$sims.list$mR, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T)),
                             apply(DD.plt$sims.list$mu, 2, FUN = function(x) quantile(x,probs=0.975,na.rm=T))),
                     lci = c(apply(DD.plt$sims.list$B, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T)),
                             apply(DD.plt$sims.list$R, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T)),
                             apply(DD.plt$sims.list$m, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T)),
                             apply(DD.plt$sims.list$mR, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T)),
                             apply(DD.plt$sims.list$mu, 2, FUN = function(x) quantile(x,probs=0.025,na.rm=T))),
                     term = factor(c(rep("Biomass (tonnes)",length(yrs)),
                              rep("Recruit Biomass (tonnes)",length(yrs)),
                              rep("FR Natural Mortality (95+ mm, instantaneous)",length(yrs)),
                              rep("Recruit Natural Mortality (85-95 mm, instantaneous)",length(yrs)),
                              rep("Fishing Mortality (instantaneous)",length(yrs))),
                              levels = c("Biomass (tonnes)","Recruit Biomass (tonnes)", "Fishing Mortality (instantaneous)",
                                        "FR Natural Mortality (95+ mm, instantaneous)","Recruit Natural Mortality (85-95 mm, instantaneous)")),
                     year = c(rep(yrs,5)))

# This is the super useful object to compare with other models...
saveRDS(res.gg,file = "D:/Github/BBn_model/Results/Models/BBn_SS_model/B_R_M_F_summarized.Rds")

# Now trim this to the 1994-2022 window
res.gg <- res.gg %>% dplyr::filter(year %in% 1994:2022)

#windows(11,11)
p.ssmod.res <- ggplot(res.gg, aes(x=year,y=med)) + geom_line(linewidth=1.5) + facet_wrap(~term,scales = 'free_y') + 
                                                   geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue') 
save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/Key_results.png",p.ssmod.res,base_height = 8.5,base_width = 11)

# Just the biomass

p.ssmod.bm <- ggplot(res.gg %>% dplyr::filter(term == "Biomass (tonnes)"), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue') + ylab("Biomass (tonnes)") + xlab("") + 
  scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,2.35e4))
save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/Biomass.png",p.ssmod.bm,base_height = 8.5,base_width = 11)

# Just the recruits

p.ssmod.rec <- ggplot(res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)"), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue') + ylab("Recruit Biomass (tonnes)") + 
  xlab("") + ylab("Recruit Biomass (tonnes)")  + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,1.75e4))
  
save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/Recruits.png",p.ssmod.rec,base_height = 8.5,base_width = 11)

# Natural mortality

p.ssmod.nat.mort <- ggplot(res.gg %>% dplyr::filter(term == "FR Natural Mortality (95+ mm, instantaneous)"), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue') + xlab("") +
  ylab("Natural mortality (Instantaneous)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.65))
save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/FR_nm.png",p.ssmod.nat.mort,base_height = 8.5,base_width = 11)

# Exploitation Rate

p.ssmod.F <- ggplot(res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)"), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue') + xlab("") +
  ylab("Exploitation Rate (Proportional)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.25))
save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/FR_F.png",p.ssmod.nat.mort,base_height = 8.5,base_width = 11)


# Now all the same plots but with missing survey years removed....

# Just the biomass

p.ssmod.bm <- ggplot(res.gg %>% dplyr::filter(term == "Biomass (tonnes)" & year < 2020), aes(x=year,y=med)) + 
  geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue') +
  geom_point(data = res.gg %>% dplyr::filter(term == "Biomass (tonnes)" & year > 2020),color='firebrick2',size=1.5) +
  geom_errorbar(data = res.gg %>% dplyr::filter(term == "Biomass (tonnes)" & year > 2020),aes(x=year,ymin=lci,ymax=uci),alpha=0.5,color='blue',size=4,width=0) + 
  ylab("Biomass (tonnes)") + xlab("") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,2.35e4))
save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/Biomass_no_surveys.png",p.ssmod.bm,base_height = 8.5,base_width = 11)

# Just the recruits

p.ssmod.rec <- ggplot(res.gg.ns %>% dplyr::filter(term == "Recruit Biomass (tonnes)"& year < 2020), aes(x=year,y=med)) + 
  geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue') + 
  geom_point(data = res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)" & year > 2020),color='firebrick2',size=1.5) +
  geom_errorbar(data = res.gg %>% dplyr::filter(term == "Recruit Biomass (tonnes)" & year > 2020),aes(x=year,ymin=lci,ymax=uci),alpha=0.5,color='blue',size=4,width=0) + 
  xlab("") + ylab("Recruit Biomass (tonnes)")  + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,1.75e4))

save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/Recruits_no_surveys.png",p.ssmod.rec,base_height = 8.5,base_width = 11)

# Natural mortality

p.ssmod.nat.mort <- ggplot(res.gg.ns %>% dplyr::filter(term == "FR Natural Mortality (95+ mm, instantaneous)" & year < 2020), aes(x=year,y=med)) + 
  geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue') +
  geom_point(data = res.gg %>% dplyr::filter(term == "FR Natural Mortality (95+ mm, instantaneous)" & year > 2020),color='firebrick2',size=1.5) +
  geom_errorbar(data = res.gg %>% dplyr::filter(term == "FR Natural Mortality (95+ mm, instantaneous)" & year > 2020),aes(x=year,ymin=lci,ymax=uci),alpha=0.5,color='blue',size=4,width=0) + 
  ylab("Natural mortality (Instantaneous)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.65))
save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/FR_nm_no_surveys.png",p.ssmod.nat.mort,base_height = 8.5,base_width = 11)

# Exploitation Rate

p.ssmod.F <- ggplot(res.gg.ns %>% dplyr::filter(term == "Fishing Mortality (instantaneous)" & year < 2020), aes(x=year,y=med)) + geom_line(color='firebrick2',linewidth=1.5) +
  geom_ribbon(aes(x=year,ymin=lci,ymax=uci),alpha=0.5,fill='blue',color='blue')  +
  geom_point(data = res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)" & year > 2020),color='firebrick2',size=1.5) +
  geom_errorbar(data = res.gg %>% dplyr::filter(term == "Fishing Mortality (instantaneous)" & year > 2020),aes(x=year,ymin=lci,ymax=uci),alpha=0.5,color='blue',size=4,width=0) + 
  ylab("Exploitation Rate (Proportional)") + xlab("") + scale_x_continuous(breaks = seq(1980,2030,by=3)) + ylim(c(0,0.41))
save_plot("D:/Github/BBn_model/Results/Figures/BBn/SSModel/FR_F_no_surveys.png",p.ssmod.nat.mort,base_height = 8.5,base_width = 11)























