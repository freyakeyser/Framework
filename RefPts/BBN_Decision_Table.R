# Here are the BBn and Sable Decision Tables.
library(lubridate)
library(kableExtra)
library(readr)
# Here's where I'm sticking the figures
mod.loc <- "D:/Framework/SFA_25_26_2024/Model/"
fun.loc <-"D:/Github/Framework/RefPts/functions/"
rp.loc <- "D:/Framework/SFA_25_26_2024/RefPts/"

# Functions we need to source for this to all work
source(paste0(fun.loc,"Decision_Table_function.R"))

# We need to get the landings for the previous year after the survey (in this case the Landings from June-December 2022)

bbn.fish <- readRDS('D:/Framework/SFA_25_26_2024/Model/Data/BBn_fish.dat.RDS')
bbn.fish$pro.repwt <- bbn.fish$pro.repwt/1000
bbn.fish$months <- month(bbn.fish$date)
psl.months <-  6:12
bbn.psl <- bbn.fish %>% dplyr::filter(year == 2022 , months %in% psl.months)
bbn.psl.total <- sum(bbn.psl$pro.repwt,na.rm=T)


years <- 1994:2022
R.size <- "75"
FR.size <- "90"
num.knots <- 20 # Going to test 10, 15, and 20
qR <- 0.45# This is for TMB (log recruit catchability) testing catchability of 0.5, test 0.3 and 0.1
init.m <- 0.2 # This is for SEAM, sets first year natural mortality, going to test 0.4, 0.15, and 0.05
# Various explorations of the g models.
g.mod <- 'g_original'
#g.mod <- 'alt_g'
#g.mod <- 'proper_g'




#seam.select <- "1994_2022_vary_m_m0_1_qR_0.5_15_knots" # This has the old von B growth parameters...
#seam.select <- "1994_2022_vary_m_m0_1_qR_0.5_10_knots"
seam.select <- paste0(min(years),"_",max(years),"_vary_m_m0_",init.m,"_qR_",qR,"_",num.knots,"_knots_",g.mod)
bbn.mod <- readRDS(paste0(mod.loc,"Results/BBn/R_75_FR_90/BBn_SEAM_model_output_",seam.select,".Rds"))


# So we can now make our decision table
# BMSY is somewhere between 6000 and 7500, so we probably should go with around 7000 as our TRP
# So then we make the USR be 80% of the BMSY/TRP... and the USR can be 30 or 40% of that.
# Alternatively, our B0 looks to be just over 10,000, so we could do 40 and 20 % of that for USR/LRP, so LRP of 2000 and USR of 4000

catchs <- seq(0,500,by=25)
# Reference Points
TRP <- 7223 # Not sure what is appropriate here, could use BMSY ((5724), or X% of B0 (80% of B0 = 7223))
USR <- 3612 #  Two options are 80% of BMSY (4579), or 40% of B0 (3612)
LRP <- 1806 # Three options, 30% BMSY (1717), 40% BMSY (2290), or 20% B0 (1806)
RR  <-   11 # This is in %, somewhere between 12-15% seems to make sense on BBn.
g.adj <- 0
gR.adj <-0
n.sims <- 1e6

# 0 growth model.
bbn.dt <- dec.tab(mod.select = "SEAM",data = bbn.mod, catch.scenarios = catchs,PSL = bbn.psl.total,
              n.sims = n.sims,TRP = TRP,USR = USR,LRP = LRP, RR = RR,g.adj=g.adj,gR.adj=gR.adj)
if(!dir.exists(paste0(rp.loc,"Results/BBn/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/"))) dir.create(paste0(rp.loc,"Results/BBn/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/"))

saveRDS(bbn.dt,paste0(rp.loc,"Results/BBn/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,"_g_adj_",g.adj,"_gR_adj_",gR.adj,".Rds"))
write_csv(bbn.dt$table,paste0(rp.loc,"Results/BBn/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,"_g_adj_",g.adj,"_gR_adj_",gR.adj,".csv"))

# Turn the decision table into something nicely formated for Word 
tab <- kableExtra::kbl(bbn.dt$table, booktabs = TRUE, escape =F, format = 'pipe',align = c('l','l','l','r','l'))#,
                       #caption = cap) %>%
  #kable_styling(full_width = F) %>% row_spec(c(2:10,12,14:18,20), bold = T) %>%
  #kable_styling(full_width = F) %>% row_spec(c(2,4,5,8:12,14,17,20), italic = T) %>%
  #add_footnote(notation = 'number',ft.note,escape=F)
saveRDS(tab,paste0(rp.loc,"Results/BBn/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/BBn_word_ready_DT_TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,"_g_adj_",g.adj,"_gR_adj_",gR.adj,".Rds"))

# pdf version of the same
tab.pdf <- kableExtra::kbl(bbn.dt$table, booktabs = TRUE, escape =F, format='latex',align = c('l','l','l','r','l'))#,
                       #caption = cap) %>%
  #kable_styling(full_width = F) %>% row_spec(c(2:10,12,14:18,20), bold = T) %>%
  #kable_styling(full_width = F) %>% row_spec(c(2,4,5,8:12,14,17,20), italic = T) %>%
  #kable_styling(latex_options = c("hold_position","scale_down")) %>%
  #add_footnote(notation = 'number',ft.note,escape=F)
saveRDS(tab.pdf,paste0(rp.loc,"Results/BBn/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/BBn_pdf_ready_DT_TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,"_g_adj_",g.adj,"_gR_adj_",gR.adj,".Rds"))
