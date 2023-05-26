# So here we'll compare the Sable the BBn models
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





############################ Model comparisons ############################ Model comparisons ############################ Model comparisons ############################ Model comparisons ###########################
repo.loc <- "D:/Github/BBn_model/"
# Sable Bayesian model results

sab.ss.mod <- readRDS("D:/Github/BBN_model/Results/Models/Sab_SS_model/B_R_M_F_summarized.Rds")
ss.mod.B <- sab.ss.mod %>% dplyr::filter(term == "Biomass (tonnes)")
ss.mod.R <- sab.ss.mod %>% dplyr::filter(term == "Recruit Biomass (tonnes)")
ss.mod.m <- sab.ss.mod %>% dplyr::filter(term == "FR Natural Mortality (95+ mm, instantaneous)")
# Remove 1994 from results to line up with the others
ss.mod.B <- ss.mod.B[2:nrow(ss.mod.B),]
ss.mod.R <- ss.mod.R[2:nrow(ss.mod.R),]
ss.mod.m <- ss.mod.m[2:nrow(ss.mod.m),]
# And add a NA row for 2023
ss.mod.B[nrow(ss.mod.B)+1,] <- NA
ss.mod.R[nrow(ss.mod.R)+1,] <- NA
ss.mod.m[nrow(ss.mod.m)+1,] <- NA


#Loading the output from TLM BBN model
tlm.moda <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_TLM_model_output_1995_2022_qR_0_5.Rds"))
tlm.mod <- tlm.moda
tlm.mod$report$totB <- tlm.mod$report$B
tlm.mod$report$totR <- tlm.mod$report$R
tlm.mod$report$mean_m <- tlm.mod$report$m
#get_parameters(tlm.mod)
#tlm.mod$report$gI <- tlm.mod$report$g


# # If running with SEBDAM this is what ya use. My fav model here is the 100 with 4 knots, nice qR, the 'process error' isn't huge, aligns with BBn, BUT CHECK PRODUCTIVITY, OUTLOOK NOT SO GOOD I think!
seam.mod <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_55_4_knots.Rds"))
sm2 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_55_10_knots.Rds"))
sm3 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_100_4_knots.Rds"))
sm4 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_100_10_knots.Rds"))
sm4a <- sm4 # This is for the 3 model comparisons
sm5 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_400_4_knots.Rds"))
sm6 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_400_10_knots.Rds"))
pred.proc<-get_processes(seam.mod)
seam.mod$report$totR <-c(seam.mod$report$totR,NA)
sm2$report$totR <-c(sm2$report$totR,NA)
sm3$report$totR <-c(sm3$report$totR,NA)
sm4$report$totR <-c(sm4$report$totR,NA)
sm5$report$totR <-c(sm5$report$totR,NA)
sm6$report$totR <-c(sm6$report$totR,NA)


pred.proc$totals$catch <- c(tlm.mod$obj$env$data$C,0)

# If I want to compare results from a variety pack of models.
n.mods <- 8
years <- 1995:2023
NY <- length(years)
cols <- c("navyblue","orange","grey","darkgreen","black","yellow",'purple','red')[1:n.mods]
comp.mods <- data.frame(years = rep(years,n.mods),
                        mod = factor(c(rep("SS",NY),rep("TLM",NY),rep("SM_R55_K4",NY),rep("SM_55_K10",NY),rep("SM_R100_K4",NY),rep("SM_R100_K10",NY),rep("SM_R400_K4",NY),rep("SM_R400_K10",NY)),#rep("SM_R400_K4",NY),rep("SM_R100_K10",NY)),
                                     levels = c("SS","TLM","SM_R55_K4","SM_55_K10","SM_R100_K4","SM_R100_K10","SM_R400_K4","SM_R400_K10"),#,"SM_R400_K4"),
                                     labels = c("SS","TLM qR = 0.5","R0 = 55, 4 knots","R0 = 55 10 knots","R0 = 100, 4 knots","R0 = 100, 10 knots","R0 = 400, 4 knots","R0 = 400, 10 knots")),#,"R0 = 148, 40 knots" )),
                        B = c(ss.mod.B$med,tlm.mod$report$B,seam.mod$report$totB,sm2$report$totB,sm3$report$totB,sm4$report$totB,sm5$report$totB,sm6$report$totB), #
                        R = c(ss.mod.R$med,tlm.mod$report$R,seam.mod$report$totR,sm2$report$totR,sm3$report$totR,sm4$report$totR,sm5$report$totR,sm6$report$totR), #sm4$report$totR,sm5$report$totR
                        m = c(ss.mod.m$med,tlm.mod$report$mean_m,seam.mod$report$mean_m,sm2$report$mean_m,sm3$report$mean_m,sm4$report$mean_m, sm5$report$mean_m, sm6$report$mean_m),
                        g = c(tlm.mod$obj$env$data$g,tlm.mod$obj$env$data$g,seam.mod$obj$env$data$gI,sm2$obj$env$data$gI,sm3$obj$env$data$gI,sm4$obj$env$data$gI,sm5$obj$env$data$gI,sm6$obj$env$data$gI),
                        gR = c(tlm.mod$obj$env$data$gR,tlm.mod$obj$env$data$gR,seam.mod$obj$env$data$gR,sm2$obj$env$data$gR,sm3$obj$env$data$gR,sm4$obj$env$data$gR,sm5$obj$env$data$gR,sm6$obj$env$data$gR),
                        Catch = rep(c(tlm.mod$obj$env$data$C,0),n.mods),
                        B.mod = c(c(ss.mod.B$med[2:length(ss.mod.B$med)],NA),
                                  c(tlm.mod$report$B[2:length(tlm.mod$report$B)],NA),
                                  c(seam.mod$report$totB[2:length(seam.mod$report$totB)],NA),
                                  c(sm2$report$totB[2:length(sm2$report$totB)],NA),
                                  c(sm3$report$totB[2:length(sm3$report$totB)],NA),
                                  c(sm4$report$totB[2:length(sm4$report$totB)],NA),
                                  c(sm5$report$totB[2:length(sm5$report$totB)],NA),
                                  c(sm6$report$totB[2:length(sm6$report$totB)],NA)),
                        m.align =  c(c(ss.mod.m$med[2:length(ss.mod.m$med)],NA),
                                  c(tlm.mod$report$mean_m[2:length(tlm.mod$report$mean_m)],NA),
                                     c(seam.mod$report$mean_m[2:length(seam.mod$report$mean_m)],NA),
                                     c(sm2$report$mean_m[2:length(sm2$report$mean_m)],NA),
                                     c(sm3$report$mean_m[2:length(sm3$report$mean_m)],NA),
                                     c(sm4$report$mean_m[2:length(sm4$report$mean_m)],NA),
                                     c(sm5$report$mean_m[2:length(sm5$report$mean_m)],NA),
                                     c(sm6$report$mean_m[2:length(sm6$report$mean_m)],NA)))
# Now do the X different models compare....
comp.mods$B.pm <- (exp(-comp.mods$m.align))*comp.mods$g*(comp.mods$B-comp.mods$Catch) 
comp.mods$R.pm <- (exp(-comp.mods$m.align))*comp.mods$gR*(comp.mods$R) 
comp.mods$B.exp <- comp.mods$B.pm + comp.mods$R.pm
comp.mods$B.diff <- comp.mods$B.exp - comp.mods$B.mod
comp.mods$B.pd <- 100*((comp.mods$B.exp - comp.mods$B.mod)/comp.mods$B.mod)


# This does the comparisons using the knot by knot calculations.
bd1 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_55_4_knots_B_differnce.Rds"))
bd2 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_55_10_knots_B_differnce.Rds"))
bd3 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_100_4_knots_B_differnce.Rds"))
bd4 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_100_10_knots_B_differnce.Rds"))
bd5 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_400_4_knots_B_differnce.Rds"))
bd6 <- readRDS(paste0(repo.loc,"Results/Models/Sab/Sab_SEAM_model_output_1995_2022_vary_m_m0_1_R0_400_10_knots_B_differnce.Rds"))

bd1.agg <- bd1 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd1.agg$mod <- "R0 = 55, 4 knots"
bd2.agg <- bd2 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd2.agg$mod <- "R0 = 55, 10 knots"
bd3.agg <- bd3 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd3.agg$mod <- "R0 = 100, 4 knots"
bd4.agg <- bd4 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd4.agg$mod <- "R0 = 100, 10 knots"
bd5.agg <- bd5 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd5.agg$mod <- "R0 = 400, 4 knots"

bd6.agg <- bd6 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd6.agg$mod <- "R0 = 400, 10 knots"

bd.all <- rbind(bd1.agg,bd2.agg,bd3.agg,bd4.agg,bd5.agg,bd6.agg)
# Also should put TLM in here.
tlm.too <- comp.mods %>% dplyr::filter(mod == "TLM qR = 0.5") %>% dplyr::select(Year = years,B.mod = B.mod,B.exp = B.exp,B.diff=B.diff,B.pd = B.pd,mod=mod)
ss.mod.aw <- comp.mods %>% dplyr::filter(mod == "SS") %>% dplyr::select(Year = years,B.mod = B.mod,B.exp = B.exp,B.diff=B.diff,B.pd = B.pd,mod=mod)
bd.all <- rbind(bd.all,tlm.too,ss.mod.aw)
bd.all$mod <- factor(bd.all$mod,levels = c("SS","TLM qR = 0.5","R0 = 55, 4 knots","R0 = 55, 10 knots","R0 = 100, 4 knots","R0 = 100, 10 knots","R0 = 400, 4 knots","R0 = 400, 10 knots"))

# Get rid of 2023 projections....
comp.mods <- comp.mods %>% dplyr::filter(years < 2023)
bd.all <- bd.all %>% dplyr::filter(Year < 2023)



pB <- ggplot(comp.mods, aes(x=years,y=B,color=mod)) + geom_point() + 
                                                      geom_line(linewidth=1.5) + 
                                                      scale_color_manual(values =cols) + 
                                                      xlab("") + ylab("Fully Recruited Biomass (tonnes)") +
                                                      #geom_hline(yintercept = median(comp.mods$B,na.rm=T),linetype = 'dashed') +
                                                      theme(legend.title = element_blank(),legend.text = element_text(color='white'))+
                                                      guides(colour=guide_legend(override.aes=list(color='white')))
pR <- ggplot(comp.mods, aes(x=years,y=R,color=mod)) + geom_point() + 
                                                      geom_line(linewidth=1.5) + 
                                                      scale_color_manual(values = cols) + 
                                                      xlab("") + ylab("Recruit Biomass (tonnes)") +
                                                      #geom_hline(yintercept = median(comp.mods$R,na.rm=T),linetype = 'dashed') +
                                                      theme(legend.title = element_blank())
pM <- ggplot(comp.mods, aes(x=years,y=m,color=mod)) + geom_point() + 
                                                      geom_line(linewidth=1.5) + 
                                                      scale_color_manual(values = cols) + 
                                                      xlab("") + ylab("Natural Mortality (inst)") +
                                                      #geom_hline(yintercept = median(comp.mods$m,na.rm=T),linetype = 'dashed') +
                                                      theme(legend.title = element_blank(),legend.text = element_text(color='white'))+
                                                      guides(colour=guide_legend(override.aes=list(color='white')))

sab.comp <- cowplot::plot_grid(pB,pR,pM,nrow=3)

save_plot(paste0(repo.loc,"Results/Figures/Sab/Model_comparison.png"),sab.comp,base_width = 11,base_height = 11)



# The difference models
p.diff <- ggplot(comp.mods, aes(x=years,y=B.diff,color=mod)) + geom_point() + 
                                                               geom_line(linewidth=1.5) + 
                                                               scale_color_manual(values =cols) +
                                                               xlab("") + ylab ("Expected - Realized Biomass (tonnes)")+
                                                               geom_hline(yintercept = 0,linetype = 'dashed') +
                                                               theme(legend.title = element_blank())

p.per.diff <- ggplot(comp.mods, aes(x=years,y=B.pd,color=mod)) + geom_point() + 
                                                                 geom_line(linewidth=1.5) + 
                                                                 scale_color_manual(values =cols)  +  
                                                                 xlab("") + ylab ("Expected - Realized Biomass (%)")+
                                                                 geom_hline(yintercept = 0,linetype = 'dashed') +
                                                                 theme(legend.title = element_blank(),legend.text = element_text(color='white'))+
                                                                 guides(colour=guide_legend(override.aes=list(color='white')))


mod.diffs <- cowplot::plot_grid(p.diff,p.per.diff,nrow=2)

save_plot(paste0(repo.loc,"Results/Figures/Sab/Models_vs_expected.png"),mod.diffs,base_width = 11,base_height = 8.5)

##### Compare with the calculation by knot##### Compare with the calculation by knot##### Compare with the calculation by knot##### Compare with the calculation by knot##### Compare with the calculation by knot

b.diffs.knot <- ggplot(data= bd.all,aes(x=Year,y=B.diff,color=mod,group=mod)) + geom_point() + 
  geom_line(size=1.5) + 
  scale_color_manual(values = cols) + 
  scale_y_continuous(breaks=seq(-1e4,1e4,by=1000)) +
  xlab("") + ylab("Expected - Realized Biomass (tonnes)") +
  geom_hline(yintercept = 0,linetype = 'dashed') +
  theme(legend.title = element_blank())


p.diffs.knot <- ggplot(data= bd.all,aes(x=Year,y=B.pd,color=mod,group=mod)) + geom_point() + 
  geom_line(size=1.5) + 
  scale_color_manual(values = cols) + 
  scale_y_continuous(breaks=seq(-100,200,by=25)) +
  xlab("") + ylab ("Expected - Realized Biomass (%)")+
  geom_hline(yintercept = 0,linetype = 'dashed') +
  theme(legend.title = element_blank(),legend.text = element_text(color='white'))+
  guides(colour=guide_legend(override.aes=list(color='white')))

mod.diffs.knot <- cowplot::plot_grid(b.diffs.knot,p.diffs.knot,nrow=2)

save_plot(paste0(repo.loc,"Results/Figures/Sab/Models_vs_expected_knot_averaged.png"),mod.diffs.knot,base_width = 11,base_height = 8.5)




# Now we can add in uncertainty, but subset this to 3 models to compare, SS, TLM, SEAM


tlm.res <- get_processes(tlm.moda)
tlm.tmp <- NULL
#tlm biomass
tlm.tmp$B.LCI <- exp(tlm.res$log_processes$log_B - 1.96*tlm.res$log_processes$se_log_B)
tlm.tmp$B.UCI <- exp(tlm.res$log_processes$log_B + 1.96*tlm.res$log_processes$se_log_B)
tlm.tmp$B <- exp(tlm.res$log_processes$log_B)
  
tlm.tmp$R.LCI <- exp(tlm.res$log_processes$log_R - 1.96*tlm.res$log_processes$se_log_R)
tlm.tmp$R.UCI <- exp(tlm.res$log_processes$log_R + 1.96*tlm.res$log_processes$se_log_R)
tlm.tmp$R <- exp(tlm.res$log_processes$log_R)
  
tlm.tmp$m.LCI <- exp(tlm.res$log_processes$log_m - 1.96*tlm.res$log_processes$se_log_m)
tlm.tmp$m.UCI <- exp(tlm.res$log_processes$log_m + 1.96*tlm.res$log_processes$se_log_m)
tlm.tmp$m     <- exp(tlm.res$log_processes$log_m)
  
tlm.tmp <- as.data.frame(do.call('cbind',tlm.tmp))
# Now SEAM model '4'
sm4.res <- get_processes(sm4a)
sm4.tmp <- NULL
  
sm4.tmp$B <- exp(sm4.res$log_tot_frame$log_totB)
sm4.tmp$B.LCI <- exp(sm4.res$log_tot_frame$log_totB - 1.96*sm4.res$log_tot_frame$se_log_totB)
sm4.tmp$B.UCI <- exp(sm4.res$log_tot_frame$log_totB + 1.96*sm4.res$log_tot_frame$se_log_totB)

sm4.tmp$R <- exp(sm4.res$log_tot_frame$log_totR)
sm4.tmp$R.LCI <- exp(sm4.res$log_tot_frame$log_totR - 1.96*sm4.res$log_tot_frame$se_log_totR)
sm4.tmp$R.UCI <- exp(sm4.res$log_tot_frame$log_totR + 1.96*sm4.res$log_tot_frame$se_log_totR)

sm4.tmp$m <- exp(sm4.res$log_tot_frame$log_mean_m)
sm4.tmp$m.LCI <- exp(sm4.res$log_tot_frame$log_mean_m - 1.96*sm4.res$log_tot_frame$se_log_mean_m)
sm4.tmp$m.UCI <- exp(sm4.res$log_tot_frame$log_mean_m + 1.96*sm4.res$log_tot_frame$se_log_mean_m)

sm4.tmp <- as.data.frame(do.call('cbind',sm4.tmp))

biomass <- data.frame(year = rep(1995:2023,3),
                      B = c(ss.mod.B$med,tlm.tmp$B,sm4.tmp$B),
                      UCI = c(ss.mod.B$uci,tlm.tmp$B.UCI,sm4.tmp$B.UCI),
                      LCI = c(ss.mod.B$lci,tlm.tmp$B.LCI,sm4.tmp$B.LCI),
                      model = factor(c(rep("SS",29),
                                       rep("TLM",29),
                                       rep("SEAM",29)),levels = c("SS","TLM","SEAM")))
biomass <- biomass %>% dplyr::filter(year < 2023)

recruits <- data.frame(year = rep(1995:2023,3),
                       R = c(ss.mod.R$med,tlm.tmp$R,sm4.tmp$R),
                       UCI = c(ss.mod.R$uci,tlm.tmp$R.UCI,sm4.tmp$R.UCI),
                       LCI = c(ss.mod.R$lci,tlm.tmp$R.LCI,sm4.tmp$R.LCI),
                       model = factor(c(rep("SS",29),
                                        rep("TLM",29),
                                        rep("SEAM",29)),levels = c("SS","TLM","SEAM")))
recruits <- recruits %>% dplyr::filter(year < 2023)
  
nat.mort <- data.frame(year = rep(1995:2023,3),
                       m = c(ss.mod.m$med,tlm.tmp$m,sm4.tmp$m),
                       UCI = c(ss.mod.m$uci,tlm.tmp$m.UCI,sm4.tmp$m.UCI),
                       LCI = c(ss.mod.m$lci,tlm.tmp$m.LCI,sm4.tmp$m.LCI),
                       model = factor(c(rep("SS",29),
                                        rep("TLM",29),
                                        rep("SEAM",29)),levels = c("SS","TLM","SEAM")))
nat.mort <- nat.mort %>% dplyr::filter(year < 2023)
  
bm.ts.plot.w.alpha <- ggplot(biomass) + geom_line(aes(year,B, color=model),linewidth=1.5) + 
                                geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 0.5) + 
                                scale_color_manual(values = cols[1:3]) + scale_fill_manual(values = cols[c(1,2,4)]) + 
                                xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) +
                                theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/Sab/three_model_biomass_comp_transparent.png"),bm.ts.plot.w.alpha,base_width = 11,base_height = 8.5)

bm.ts.plot.no.alpha <- ggplot(biomass) + geom_line(aes(year,B, color=model),linewidth=1.5) + 
                                         geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 1) + 
                                         scale_color_manual(values = cols[1:3]) + scale_fill_manual(values = cols[c(1,2,4)]) + 
                                         xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) +
                                         theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/Sab/three_model_biomass_comp_opaque.png"),bm.ts.plot.no.alpha,base_width = 11,base_height = 8.5)



rec.ts.plot.w.alpha <- ggplot(recruits) + geom_line(aes(year,R, color=model),linewidth=1.5) + 
                                          geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 0.5) + 
                                          scale_color_manual(values = cols[1:3]) + scale_fill_manual(values = cols[c(1,2,3)]) + 
                                          xlab("") + ylab("Recruit Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))+
                                          theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/Sab/three_model_recruit_comp_transparent.png"),rec.ts.plot.w.alpha,base_width = 11,base_height = 8.5)

rec.ts.plot.no.alpha <- ggplot(recruits) + geom_line(aes(year,R, color=model),linewidth=1.5) + 
                                           geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 1) + 
                                           scale_color_manual(values = cols[1:3]) + scale_fill_manual(values = cols[c(1,2,3)]) + 
                                           xlab("") + ylab("Recruit Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))+
                                           theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/Sab/three_model_recruit_comp_opaque.png"),rec.ts.plot.no.alpha,base_width = 11,base_height = 8.5)



m.ts.plot.w.aplha <- ggplot(nat.mort) + geom_line(aes(year,m, color=model),linewidth=1.5) + 
                                        geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 0.5) + 
                                        scale_color_manual(values = cols[1:3]) + scale_fill_manual(values = cols[c(1,2,3)]) + 
                                        xlab("") + ylab("Natural mortality (instantaneous)") + scale_x_continuous(breaks = seq(1980,2030,by=3))+
                                        theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/Sab/three_model_nat_mort_comp_transparent.png"),m.ts.plot.w.aplha,base_width = 11,base_height = 8.5)



m.ts.plot.no.aplha <- ggplot(nat.mort) + geom_line(aes(year,m, color=model),linewidth=1.5) + 
                                         geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 1) + 
                                         scale_color_manual(values = cols[1:3]) + scale_fill_manual(values = cols[c(1,2,3)]) + 
                                         xlab("") + ylab("Natural mortality (instantaneous)") + scale_x_continuous(breaks = seq(1980,2030,by=3))+
                                         theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/Sab/three_model_nat_mort_comp_opaque.png"),m.ts.plot.no.aplha,base_width = 11,base_height = 8.5)


############################ BBN ############################ BBN BBN ############################ BBN BBN ############################ BBN BBN ########################## 

############################ Model comparisons ############################ Model comparisons ############################ Model comparisons ############################ Model comparisons ###########################
repo.loc <- "D:/Github/BBn_model/"
# Sable Bayesian model results

bbn.ss.mod <- readRDS("D:/Github/BBN_model/Results/Models/BBn_SS_model/B_R_M_F_summarized.Rds")
# Filter by term and start with 1994
ss.mod.B <- bbn.ss.mod %>% dplyr::filter(term == "Biomass (tonnes)" & year >=1994)
ss.mod.R <- bbn.ss.mod %>% dplyr::filter(term == "Recruit Biomass (tonnes)" & year >=1994)
ss.mod.m <- bbn.ss.mod %>% dplyr::filter(term == "FR Natural Mortality (95+ mm, instantaneous)" & year >=1994)

# And add a NA row for 2022 and  2023 We can drop one of these once we have the new model running, you'll get an error and be annoyed!
ss.mod.B[nrow(ss.mod.B)+1,] <- NA
ss.mod.R[nrow(ss.mod.R)+1,] <- NA
ss.mod.m[nrow(ss.mod.m)+1,] <- NA
# Second round of NAs
ss.mod.B[nrow(ss.mod.B)+1,] <- NA
ss.mod.R[nrow(ss.mod.R)+1,] <- NA
ss.mod.m[nrow(ss.mod.m)+1,] <- NA


#Loading the output from TLM BBN model
tlma.mod <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_TLM_model_output_1994_2022_qR_0_5.Rds"))
tlm.mod <- tlma.mod
tlm.mod$report$totB <- tlm.mod$report$B
tlm.mod$report$totR <- tlm.mod$report$R
tlm.mod$report$mean_m <- tlm.mod$report$m
#get_parameters(tlm.mod)
#tlm.mod$report$gI <- tlm.mod$report$g

# # If running with SEBDAM this is what ya use. My fav model here is the 100 with 4 knots, nice qR, the 'process error' isn't huge, aligns with BBn, BUT CHECK PRODUCTIVITY, OUTLOOK NOT SO GOOD I think!
seam.mod <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_55_10_knots.Rds"))
sm2 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_55_20_knots.Rds"))
sm3 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_150_10_knots.Rds"))
sm3a <- sm3
sm4 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_150_20_knots.Rds"))
sm5 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_150_40_knots.Rds"))

pred.proc<-get_processes(seam.mod)
seam.mod$report$totR <-c(seam.mod$report$totR,NA)
sm2$report$totR <-c(sm2$report$totR,NA)
sm3$report$totR <-c(sm3$report$totR,NA)
sm4$report$totR <-c(sm4$report$totR,NA)
sm5$report$totR <-c(sm5$report$totR,NA)
pred.proc$totals$catch <- c(tlm.mod$obj$env$data$C,0)

# If I want to compare results from a variety pack of models.
n.mods <- 7
years <- 1994:2023
NY <- length(years)
# Note I'm trying to play with the order here so the 'darkgreen' color is the one chosen for the 3 model comparison
cols <- c("navyblue","orange","grey","black","darkgreen","chartreuse2",'purple')[1:n.mods]
g.vec <- rep(c(tlm.mod$obj$env$data$g,tlm.mod$obj$env$data$g[length(tlm.mod$obj$env$data$g)]),n.mods)
gR.vec <- rep(c(tlm.mod$obj$env$data$gR,tlm.mod$obj$env$data$gR[length(tlm.mod$obj$env$data$gR)]),n.mods)
comp.mods <- data.frame(years = rep(years,n.mods),
                        mod = factor(c(rep("SS",NY),rep("TLM",NY),rep("SM_55_K10",NY),rep("SM_R55_K20",NY),rep("SM_R150_K10",NY),rep("SM_R150_K20",NY),rep("SM_R150_K40",NY)),#rep("SM_R400_K4",NY),rep("SM_R100_K10",NY)),
                                     levels = c("SS","TLM","SM_55_K10","SM_R55_K20","SM_R150_K10","SM_R150_K20","SM_R150_K40"),#,"SM_R400_K4"),
                                     labels = c("SS","TLM qR = 0.5","R0 = 55 10 knots","R0 = 55, 20 knots","R0 = 150, 10 knots","R0 = 150, 20 knots","R0 = 150, 40 knots")),#,"R0 = 148, 40 knots" )),
                        B = c(ss.mod.B$med,tlm.mod$report$B,seam.mod$report$totB,sm2$report$totB,sm3$report$totB,sm4$report$totB,sm5$report$totB), #
                        R = c(ss.mod.R$med,tlm.mod$report$R,seam.mod$report$totR,sm2$report$totR,sm3$report$totR,sm4$report$totR,sm5$report$totR), #sm4$report$totR,sm5$report$totR
                        m = c(ss.mod.m$med,tlm.mod$report$mean_m,seam.mod$report$mean_m,sm2$report$mean_m,sm3$report$mean_m,sm4$report$mean_m, sm5$report$mean_m),
                        g = c(g.vec),
                        gR = c(gR.vec),
                        Catch = rep(c(tlm.mod$obj$env$data$C,0),n.mods),
                        B.mod = c(c(ss.mod.B$med[2:length(ss.mod.B$med)],NA),
                                  c(tlm.mod$report$B[2:length(tlm.mod$report$B)],NA),
                                  c(seam.mod$report$totB[2:length(seam.mod$report$totB)],NA),
                                  c(sm2$report$totB[2:length(sm2$report$totB)],NA),
                                  c(sm3$report$totB[2:length(sm3$report$totB)],NA),
                                  c(sm4$report$totB[2:length(sm4$report$totB)],NA),
                                  c(sm5$report$totB[2:length(sm5$report$totB)],NA)),
                        m.align =  c(c(ss.mod.m$med[2:length(ss.mod.m$med)],NA),
                                     c(tlm.mod$report$mean_m[2:length(tlm.mod$report$mean_m)],NA),
                                     c(seam.mod$report$mean_m[2:length(seam.mod$report$mean_m)],NA),
                                     c(sm2$report$mean_m[2:length(sm2$report$mean_m)],NA),
                                     c(sm3$report$mean_m[2:length(sm3$report$mean_m)],NA),
                                     c(sm4$report$mean_m[2:length(sm4$report$mean_m)],NA),
                                     c(sm5$report$mean_m[2:length(sm5$report$mean_m)],NA)))
# Nowdo the X different models compare....
comp.mods$B.pm <- (exp(-comp.mods$m.align))*comp.mods$g*(comp.mods$B-comp.mods$Catch) 
comp.mods$R.pm <- (exp(-comp.mods$m.align))*comp.mods$gR*(comp.mods$R) 
comp.mods$B.exp <- comp.mods$B.pm + comp.mods$R.pm
comp.mods$B.diff <- comp.mods$B.exp - comp.mods$B.mod
comp.mods$B.pd <- 100*((comp.mods$B.exp - comp.mods$B.mod)/comp.mods$B.mod)




# This does the comparisons using the knot by knot calculations.
bd1 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_55_10_knots_B_differnce.Rds"))
bd2 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_55_20_knots_B_differnce.Rds"))
bd3 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_150_10_knots_B_differnce.Rds"))
bd4 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_150_20_knots_B_differnce.Rds"))
bd5 <- readRDS(paste0(repo.loc,"Results/Models/BBn/BBn_SEAM_model_output_1994_2022_vary_m_m0_1_R0_150_40_knots_B_differnce.Rds"))

bd1.agg <- bd1 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd1.agg$mod <- "R0 = 55, 10 knots"
bd2.agg <- bd2 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd2.agg$mod <- "R0 = 55, 20 knots"
bd3.agg <- bd3 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd3.agg$mod <- "R0 = 150, 10 knots"
bd4.agg <- bd4 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd4.agg$mod <- "R0 = 150, 20 knots"
bd5.agg <- bd5 %>% as.data.frame() %>% dplyr::group_by(Year) %>% dplyr::summarise(B.mod = sum(B.mod),
                                                                                  B.exp = sum(B.exp),
                                                                                  B.diff = sum(B.exp) - sum(B.mod),
                                                                                  B.pd = 100*(sum(B.exp) - sum(B.mod))/sum(B.mod)) 
bd5.agg$mod <- "R0 = 150, 40 knots"

bd.all <- rbind(bd1.agg,bd2.agg,bd3.agg,bd4.agg,bd5.agg)
# Also should put TLM in here.
tlm.too <- comp.mods %>% dplyr::filter(mod == "TLM qR = 0.5") %>% dplyr::select(Year = years,B.mod = B.mod,B.exp = B.exp,B.diff=B.diff,B.pd = B.pd,mod=mod)
ss.mod.aw <- comp.mods %>% dplyr::filter(mod == "SS") %>% dplyr::select(Year = years,B.mod = B.mod,B.exp = B.exp,B.diff=B.diff,B.pd = B.pd,mod=mod)
bd.all <- rbind(bd.all,tlm.too,ss.mod.aw)
bd.all$mod <- factor(bd.all$mod,levels = c("SS","TLM qR = 0.5","R0 = 55, 10 knots","R0 = 55, 20 knots","R0 = 150, 10 knots","R0 = 150, 20 knots","R0 = 150, 40 knots"))

# Get rid of 2023 projections....
comp.mods <- comp.mods %>% dplyr::filter(years < 2023)
bd.all <- bd.all %>% dplyr::filter(Year < 2023)


# Here are the aggregated method figures.
pB <- ggplot(comp.mods, aes(x=years,y=B,color=mod)) + geom_point() + 
                                                      geom_line(size=1.5) + 
                                                      scale_color_manual(values =cols) + 
                                                      scale_y_continuous(breaks=seq(0,3e4,by=1000)) +
                                                      xlab("") + ylab("Fully Recruited Biomass (tonnes)") +
                                                      theme(legend.title = element_blank(),legend.text = element_text(color='white'))+
                                                      guides(colour=guide_legend(override.aes=list(color='white')))

pR <- ggplot(comp.mods, aes(x=years,y=R,color=mod)) + geom_point() + 
                                                      geom_line(size=1.5) + 
                                                      scale_color_manual(values = cols) + 
                                                      scale_y_continuous(breaks=seq(0,3e4,by=1000)) +
                                                      xlab("") + ylab("Recruit Biomass (tonnes)") +
                                                      theme(legend.title = element_blank())

pM <- ggplot(comp.mods, aes(x=years,y=m,color=mod)) + geom_point() + 
                                                      geom_line(size=1.5) + 
                                                      scale_color_manual(values = cols) + 
                                                      xlab("") + ylab("Natural Mortality (inst)") +
                                                      theme(legend.title = element_blank(),legend.text = element_text(color='white'))+
                                                      guides(colour=guide_legend(override.aes=list(color='white')))

BBn.comp <- cowplot::plot_grid(pB,pR,pM,nrow=3)

save_plot(paste0(repo.loc,"Results/Figures/BBn/Model_comparison.png"),BBn.comp,base_width = 11,base_height = 11)

p.diff <- ggplot(comp.mods, aes(x=years,y=B.diff,color=mod)) + geom_point() + 
                                                                geom_line(size=1.5) + 
                                                                scale_color_manual(values =cols) +
                                                                scale_y_continuous(breaks=seq(-1e4,1e4,by=1000)) +
                                                                xlab("") + ylab ("Expected - Realilzed Biomass (tonnes)")+
                                                                theme(legend.title = element_blank())

p.per.diff <- ggplot(comp.mods, aes(x=years,y=B.pd,color=mod)) + geom_point() + 
                                                                  geom_line(size=1.5) + 
                                                                  scale_color_manual(values =cols)  +  
                                                                  scale_y_continuous(breaks=seq(-100,200,by=25)) +
                                                                  xlab("") + ylab ("Expected - Realized Biomass (%)")+
                                                                  theme(legend.title = element_blank(),legend.text = element_text(color='white'))+
                                                                  guides(colour=guide_legend(override.aes=list(color='white')))


mod.diffs <- cowplot::plot_grid(p.diff,p.per.diff,nrow=2)

save_plot(paste0(repo.loc,"Results/Figures/BBn/Models_vs_expected.png"),mod.diffs,base_width = 11,base_height = 8.5)


##### Compare with the calculation by knot##### Compare with the calculation by knot##### Compare with the calculation by knot##### Compare with the calculation by knot##### Compare with the calculation by knot

b.diffs.knot <- ggplot(data= bd.all,aes(x=Year,y=B.diff,color=mod,group=mod)) + geom_point() + 
                                                                                geom_line(size=1.5) + 
                                                                                scale_color_manual(values = cols) + 
                                                                                scale_y_continuous(breaks=seq(-1e4,1e4,by=1000)) +
                                                                                xlab("") + ylab("Realized - Modeled Biomass (tonnes)") +
                                                                                theme(legend.title = element_blank())


p.diffs.knot <- ggplot(data= bd.all,aes(x=Year,y=B.pd,color=mod,group=mod)) + geom_point() + 
                                                                              geom_line(size=1.5) + 
                                                                              scale_color_manual(values = cols) + 
                                                                              scale_y_continuous(breaks=seq(-100,200,by=25)) +
                                                                              xlab("") + ylab ("Realized - Modeled (%)")+
                                                                              theme(legend.title = element_blank(),legend.text = element_text(color='white'))+
                                                                              guides(colour=guide_legend(override.aes=list(color='white')))

mod.diffs.knot <- cowplot::plot_grid(b.diffs.knot,p.diffs.knot,nrow=2)

save_plot(paste0(repo.loc,"Results/Figures/BBn/Models_vs_expected_knot_averaged.png"),mod.diffs.knot,base_width = 11,base_height = 8.5)


# Now we can add in uncertainty, but subset this to 3 models to compare, SS, TLM, SEAM


tlm.res <- get_processes(tlma.mod)
tlm.tmp <- NULL
#tlm biomass
tlm.tmp$B.LCI <- exp(tlm.res$log_processes$log_B - 1.96*tlm.res$log_processes$se_log_B)
tlm.tmp$B.UCI <- exp(tlm.res$log_processes$log_B + 1.96*tlm.res$log_processes$se_log_B)
tlm.tmp$B <- exp(tlm.res$log_processes$log_B)

tlm.tmp$R.LCI <- exp(tlm.res$log_processes$log_R - 1.96*tlm.res$log_processes$se_log_R)
tlm.tmp$R.UCI <- exp(tlm.res$log_processes$log_R + 1.96*tlm.res$log_processes$se_log_R)
tlm.tmp$R <- exp(tlm.res$log_processes$log_R)

tlm.tmp$m.LCI <- exp(tlm.res$log_processes$log_m - 1.96*tlm.res$log_processes$se_log_m)
tlm.tmp$m.UCI <- exp(tlm.res$log_processes$log_m + 1.96*tlm.res$log_processes$se_log_m)
tlm.tmp$m     <- exp(tlm.res$log_processes$log_m)

tlm.tmp <- as.data.frame(do.call('cbind',tlm.tmp))
# Now SEAM model '3'
sm3.res <- get_processes(sm3a)
sm3.tmp <- NULL

sm3.tmp$B <- exp(sm3.res$log_tot_frame$log_totB)
sm3.tmp$B.LCI <- exp(sm3.res$log_tot_frame$log_totB - 1.96*sm3.res$log_tot_frame$se_log_totB)
sm3.tmp$B.UCI <- exp(sm3.res$log_tot_frame$log_totB + 1.96*sm3.res$log_tot_frame$se_log_totB)

sm3.tmp$R <- exp(sm3.res$log_tot_frame$log_totR)
sm3.tmp$R.LCI <- exp(sm3.res$log_tot_frame$log_totR - 1.96*sm3.res$log_tot_frame$se_log_totR)
sm3.tmp$R.UCI <- exp(sm3.res$log_tot_frame$log_totR + 1.96*sm3.res$log_tot_frame$se_log_totR)

sm3.tmp$m <- exp(sm3.res$log_tot_frame$log_mean_m)
sm3.tmp$m.LCI <- exp(sm3.res$log_tot_frame$log_mean_m - 1.96*sm3.res$log_tot_frame$se_log_mean_m)
sm3.tmp$m.UCI <- exp(sm3.res$log_tot_frame$log_mean_m + 1.96*sm3.res$log_tot_frame$se_log_mean_m)

sm3.tmp <- as.data.frame(do.call('cbind',sm3.tmp))



biomass <- data.frame(year = rep(1994:2023,3),
                      B = c(ss.mod.B$med,tlm.tmp$B,sm3.tmp$B),
                      UCI = c(ss.mod.B$uci,tlm.tmp$B.UCI,sm3.tmp$B.UCI),
                      LCI = c(ss.mod.B$lci,tlm.tmp$B.LCI,sm3.tmp$B.LCI),
                      model = factor(c(rep("SS",30),
                                       rep("TLM",30),
                                       rep("SEAM",30)),levels = c("TLM","SS","SEAM")))
biomass <- biomass %>% dplyr::filter(year < 2023)


recruits <- data.frame(year = rep(1994:2023,3),
                       R = c(ss.mod.R$med,tlm.tmp$R,sm3.tmp$R),
                       UCI = c(ss.mod.R$uci,tlm.tmp$R.UCI,sm3.tmp$R.UCI),
                       LCI = c(ss.mod.R$lci,tlm.tmp$R.LCI,sm3.tmp$R.LCI),
                       model = factor(c(rep("SS",30),
                                        rep("TLM",30),
                                        rep("SEAM",30)),levels = c("SEAM","SS","TLM")))
recruits <- recruits %>% dplyr::filter(year < 2023)


nat.mort <- data.frame(year = rep(1994:2023,3),
                       m = c(ss.mod.m$med,tlm.tmp$m,sm3.tmp$m),
                       UCI = c(ss.mod.m$med,tlm.tmp$m.UCI,sm3.tmp$m.UCI),
                       LCI = c(ss.mod.m$lci,tlm.tmp$m.LCI,sm3.tmp$m.LCI),
                       model = factor(c(rep("SS",30),
                                        rep("TLM",30),
                                        rep("SEAM",30)),levels = c("TLM","SEAM","SS")))
nat.mort <- nat.mort %>% dplyr::filter(year < 2023)


bm.ts.plot.w.alpha <- ggplot(biomass) + geom_line(aes(year,B, color=model),linewidth=1.5) + 
                                        geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 0.5) + 
                                        scale_color_manual(values = cols[c(2,5,1)]) + scale_fill_manual(values = cols[c(2,1,5)]) + 
                                        xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) +
                                        theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/BBn/three_model_biomass_comp_transparent.png"),bm.ts.plot.w.alpha,base_width = 11,base_height = 8.5)

bm.ts.plot.no.alpha <- ggplot(biomass) +  geom_line(aes(year,B, color=model),linewidth=1.5) + 
                                          geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 1) + 
                                          scale_color_manual(values = cols[c(2,5,1)]) + scale_fill_manual(values = cols[c(2,1,5)]) + 
                                          xlab("") + ylab("Fully Recruited Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3)) +
                                          theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/BBn/three_model_biomass_comp_opaque.png"),bm.ts.plot.no.alpha,base_width = 11,base_height = 8.5)



rec.ts.plot.w.alpha <- ggplot(recruits) + geom_line(aes(year,R, color=model),linewidth=1.5) + 
                                          geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 0.5) + 
                                          scale_color_manual(values = cols[c(2,5,1)]) + scale_fill_manual(values = cols[c(5,1,2)]) + 
                                          xlab("") + ylab("Recruit Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))+
                                          theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/BBn/three_model_recruit_comp_transparent.png"),rec.ts.plot.w.alpha,base_width = 11,base_height = 8.5)

rec.ts.plot.no.alpha <- ggplot(recruits) +  geom_line(aes(year,R, color=model),linewidth=1.5) + 
                                            geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 1) + 
                                            scale_color_manual(values = cols[c(2,5,1)]) + scale_fill_manual(values = cols[c(5,1,2)]) + 
                                            xlab("") + ylab("Recruit Biomass (tonnes)") + scale_x_continuous(breaks = seq(1980,2030,by=3))+
                                            theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/BBn/three_model_recruit_comp_opaque.png"),rec.ts.plot.no.alpha,base_width = 11,base_height = 8.5)



m.ts.plot.w.aplha <- ggplot(nat.mort) + geom_line(aes(year,m, color=model),linewidth=1.5) + 
                                        geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 0.5) + 
                                        scale_color_manual(values = cols[c(2,5,1)]) + scale_fill_manual(values = cols[c(2,5,1)]) + 
                                        xlab("") + ylab("Natural mortality (instantaneous)") + scale_x_continuous(breaks = seq(1980,2030,by=3))+
                                        theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/BBn/three_model_nat_mort_comp_transparent.png"),m.ts.plot.w.aplha,base_width = 11,base_height = 8.5)



m.ts.plot.no.aplha <- ggplot(nat.mort) +  geom_line(aes(year,m, color=model),linewidth=1.5) + 
                                          geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year,fill=model,color=model),alpha = 1) + 
                                          scale_color_manual(values = cols[c(2,5,1)]) + scale_fill_manual(values = cols[c(2,5,1)]) + 
                                          xlab("") + ylab("Natural mortality (instantaneous)") + scale_x_continuous(breaks = seq(1980,2030,by=3))+
                                          theme(legend.title = element_blank())

save_plot(paste0(repo.loc,"Results/Figures/BBn/three_model_nat_mort_comp_opaque.png"),m.ts.plot.no.aplha,base_width = 11,base_height = 8.5)
