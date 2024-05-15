# We are running a model for Sable using the Bayesian SS Model. This extracts the bits of code from Freya's scripts to make things work.
funs <- c("https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Model/projections.r",
          "https://raw.githubusercontent.com/Mar-scal/Assessment_fns/master/Model/decision.r")


for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}

library(sf)
library(sp)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(R2jags)
#library(SSModel)

theme_set(theme_few(base_size = 22))

repo.loc <- "D:/Framework/SFA_25_26_2024/Model/"
retro.years <- 2005:2021


# Load in the full model
load(paste0(repo.loc,"Results/Sab_SS_model/R_75_FR_90/Sable_SSmodel_results.RData"))
# Get the catch data for all years...
DD.base <- DD.out
Catch.dat <- data.frame(year=yrs,catch = DD.out$data$C)
# Now we can get DD objects for all the retro years....
DD.all <- list()
sims <- c("normal","no growth", "no fully-recruited growth","zp","Median mortality","Median productivity","Median recruitment")
n.sims <- length(sims)
for(i in retro.years)
{
  load(paste0(repo.loc,"Results/Sab_SS_model/R_75_FR_90/Retros/Sable_SSmodel_results_1994_",i,".RData"))
  DD.all[[as.character(i)]] <- DD.out
}

res <- list()
for(s in 1:n.sims)
{
  tmp.res <- list()
  for(y in retro.years)
  {
    mod.tmp <- DD.all[[as.character(y)]]
    if(y < 2021) 
    {
      real.B <- median(DD.all[[as.character(y+1)]]$sims.list$B[,ncol(DD.all[[as.character(y+1)]]$sims.list$B)])
      real.B.LCI <- as.numeric(quantile(DD.all[[as.character(y+1)]]$sims.list$B[,ncol(DD.all[[as.character(y+1)]]$sims.list$B)],probs=0.025))
      real.B.UCI <- as.numeric(quantile(DD.all[[as.character(y+1)]]$sims.list$B[,ncol(DD.all[[as.character(y+1)]]$sims.list$B)],probs=0.975))
    } # end y < 2021
    if(y == 2021) 
    {
      real.B <- median(DD.base$sims.list$B[,ncol(DD.base$sims.list$B)])
      real.B.LCI <- as.numeric(quantile(DD.base$sims.list$B[,ncol(DD.base$sims.list$B)],probs=0.025))
      real.B.UCI <- as.numeric(quantile(DD.base$sims.list$B[,ncol(DD.base$sims.list$B)],probs=0.975))
    } # end y = 2021
    rem <- Catch.dat$catch[Catch.dat$year == (y+1)]
    if(sims[s]=='no fully-recruited growth')  mod.tmp$data$g  <- mod.tmp$data$g/mod.tmp$data$g
    if(sims[s]=='no growth')  
    {
      mod.tmp$data$g  <- mod.tmp$data$g/mod.tmp$data$g
      mod.tmp$data$gR  <- mod.tmp$data$gR/mod.tmp$data$gR
    }
    if(sims[s]=='zp')  
    {
      mod.tmp$data$g  <- mod.tmp$data$g/mod.tmp$data$g
      mod.tmp$data$gR  <- mod.tmp$data$gR/mod.tmp$data$gR
      mod.tmp$sims.list$m[,ncol(mod.tmp$sims.list$m)] <- 0
      mod.tmp$sims.list$mR[,ncol(mod.tmp$sims.list$mR)] <- 0
      mod.tmp$sims.list$r[,ncol(mod.tmp$sims.list$r)] <- 0
    }
    
    if(sims[s]=='Median mortality')  
    {
      mod.tmp$sims.list$m[,ncol(mod.tmp$sims.list$m)] <- median(mod.tmp$median$m)
      mod.tmp$sims.list$mR[,ncol(mod.tmp$sims.list$mR)] <- median(mod.tmp$median$mR)
    }
    if(sims[s]=='Median productivity')  
    {
      mod.tmp$data$g[length(mod.tmp$data$g)]  <- median(mod.tmp$data$g)
      mod.tmp$data$gR[length(mod.tmp$data$gR)]  <- median(mod.tmp$data$gR)
      mod.tmp$sims.list$m[,ncol(mod.tmp$sims.list$m)] <- median(mod.tmp$median$m)
      mod.tmp$sims.list$mR[,ncol(mod.tmp$sims.list$mR)] <-  median(mod.tmp$median$mR)
      mod.tmp$sims.list$r[,ncol(mod.tmp$sims.list$r)] <- median(mod.tmp$median$r)
    }
    
    if(sims[s]=='Median recruitment')   mod.tmp$sims.list$r[,ncol(mod.tmp$sims.list$r)] <- median(mod.tmp$median$r)
    
    proj <- projections(mod.tmp,rem)
    proj.B <- median(proj$sims.list$Bmed.p)
    proj.B.LCI <- quantile(proj$sims.list$Bmed.p,probs = 0.025)
    proj.B.UCI <- quantile(proj$sims.list$Bmed.p,probs = 0.975)
    
    tmp.res[[as.character(y)]] <- data.frame(year=y+1,
                                             real.B = real.B,
                                             real.B.LCI = real.B.LCI,
                                             real.B.UCI = real.B.UCI,
                                             proj.B=proj.B,
                                             proj.B.LCI = proj.B.LCI,
                                             proj.B.UCI = proj.B.UCI,
                                             B.diff = proj.B-real.B,
                                             PB.diff = 100*(proj.B-real.B)/real.B,
                                             scenario = sims[s])
  } # end retro.years
  res[[sims[s]]] <- do.call('rbind',tmp.res)
}# end n.sims

res.fin <- do.call("rbind",res)

saveRDS(res.fin,paste0(repo.loc,"/Results/Sab_SS_model/R_75_FR_90/Sab_prediction_evaluation.Rds"))
res.fin <- readRDS(paste0(repo.loc,"/Results/Sab_SS_model/R_75_FR_90/Sab_prediction_evaluation.Rds"))

# Now run a decision Table with for 2023 so we can compare with what was observed...
# Given the results of the last few years we should run this with the growth terms all = 1.
# As long as we see the large negative process error term we should do it this way
# Once the process error returns to normal I'll suggest we use the 'normal' method for projection.
DD.dt <- DD.base
DD.dt$data$g[length(DD.dt$data$g)]  <- 1
#DD.dt$data$gR[length(DD.dt$data$gR)]  <- 1

# Now we can run the projection and then the decision table. Removals from June 2022 to May 2023 were 227 tonnes
proj.dt <- projections(DD.dt,seq(0,200,by=10))
# Now the dection table
dt.raw <- decision(proj.dt,bank= "Sab")

dt.clean <- data.frame(Catch = dt.raw$Catch,
                       Exploitation = round(100*dt.raw$mu,digits=1),
                       Biomass = round(apply(proj.dt$sims.list$B.p,2,median),digits=0),
                       Biomass.change.per = round(100*(apply(proj.dt$sims.list$B.p,2,median) - proj.dt$median$B[length(proj.dt$median$B)])/proj.dt$median$B[length(proj.dt$median$B)],digits=1),
                       Biomass.change.tonnes = round(apply(proj.dt$sims.list$B.p,2,median) - proj.dt$median$B[length(proj.dt$median$B)],digits=0),
                       P.decline = round(dt.raw$p.decline,digits=2))

names(dt.clean) <- c("Catch (tonnes)",
                     "Exploitation (%)",
                     "Biomass (tonnes)",
                     "Biomass change (%)",
                     "Biomass change (tonnes)",
                     "Probability of Decline")

saveRDS(dt.clean,paste0(repo.loc,"/Results/Sab_SS_model/R_75_FR_90/Sab_Decision_Table.Rds"))

tab <- kableExtra::kbl(dt.clean, booktabs = TRUE, escape =F, format = 'pipe',align = c('l','l','l','r','l'))#,
#caption = cap) %>%
#kable_styling(full_width = F) %>% row_spec(c(2:10,12,14:18,20), bold = T) %>%
#kable_styling(full_width = F) %>% row_spec(c(2,4,5,8:12,14,17,20), italic = T) %>%
#add_footnote(notation = 'number',ft.note,escape=F)
saveRDS(tab,paste0(repo.loc,"/Results/Sab_SS_model/R_75_FR_90/Sab_word_ready DT.Rds"))

# pdf version of the same
tab.pdf <- kableExtra::kbl(dt.clean, booktabs = TRUE, escape =F, format='latex',align = c('l','l','l','r','l'))#,
#caption = cap) %>%
#kable_styling(full_width = F) %>% row_spec(c(2:10,12,14:18,20), bold = T) %>%
#kable_styling(full_width = F) %>% row_spec(c(2,4,5,8:12,14,17,20), italic = T) %>%
#kable_styling(latex_options = c("hold_position","scale_down")) %>%
#add_footnote(notation = 'number',ft.note,escape=F)
saveRDS(tab.pdf,paste0(repo.loc,"/Results/Sab_SS_model/R_75_FR_90/Sab_pdf_ready DT.Rds"))

tst <- res.fin |> collapse::fgroup_by(scenario) |> collapse::fsummarise(med.pd = median(PB.diff,na.rm=T),
                                                                        mod.bd = median(B.diff,na.rm=T))

ggplot(res.fin |> collapse::fsubset(year < 2015)) + geom_violin(aes(x=scenario,y=PB.diff)) + geom_hline(yintercept=0)

res.fin


ggplot(res.fin) + geom_violin(aes(x=scenario,y=abs(B.diff/1000)),draw_quantiles = c(0.5))
ggplot(res.fin) + geom_violin(aes(x=scenario,y=B.diff/1000),draw_quantiles = c(0.5)) + geom_hline(yintercept=0)
ggplot(res.fin) + geom_violin(aes(x=scenario,y=PB.diff),draw_quantiles = c(0.5)) + geom_hline(yintercept=0)


ggplot(res.fin) + geom_point(aes(x=real.B,y=proj.B)) + 
                  facet_wrap(~scenario) + 
                  geom_abline(slope=1,intercept=0)

ggplot(res.fin) + geom_line(aes(x=year,y=proj.B)) + 
                  geom_ribbon(aes(x=year,ymin=proj.B.LCI,ymax=proj.B.UCI),fill='firebrick2',color='firebrick2',alpha=0.2)+
                  geom_line(aes(x=year,y=real.B),color='blue')+
                  geom_ribbon(aes(x=year,ymin=real.B.LCI,ymax=real.B.UCI),fill='blue',color='blue',alpha=0.2)+
                  facet_wrap(~scenario) 


