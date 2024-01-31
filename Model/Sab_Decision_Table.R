# Here is the Sable Decision Table.
library(lubridate)
# Here's where I'm sticking the figures
mod.loc <- "D:/Framework/SFA_25_26_2024/Model/"
fun.loc <-"D:/Github/Framework/Model/functions/"
rp.loc <- "D:/Framework/SFA_25_26_2024/RefPts/"

# Functions we need to source for this to all work
source(paste0(fun.loc,"Decision_Table_function.R"))

# We need to get the landings for the previous year after the survey (in this case the Landings from June-December 2022)


# Recruit and FR sizes
R.size <- 75
FR.size <- 90
qR <- 0.33
init.m <- 0.2
g.mod <- 'g_original'
#g.mod <- 'alt_g'
#g.mod <- 'proper_g'
num.knots <- 10

# And a Sable decision table 
seam.select <- paste0("1994_2022_vary_m_m0_",init.m,"_qR_",qR,"_",num.knots,"_knots_",g.mod)
sab.mod <- readRDS(paste0(mod.loc,"Results/Sab/R_75_FR_90/Sab_SEAM_model_output_",seam.select,".Rds"))
#rp.res <- readRDS(paste0(mod.loc,"Results/Sab/R_75_FR_90/",seam.select,"/RPs/SEAM_g_bp_8000_dist_1_0.05_r_bp_dist_3500_0.034_0.23_m_bp_dist_3000_0.54_0.1.Rds"))

sab.fish <- readRDS(paste0(mod.loc,'Data/Sab_fish.dat.RDS'))
sab.fish$pro.repwt <- sab.fish$pro.repwt/1000
sab.fish$months <- month(sab.fish$date)
psl.months <-  6:12
sab.psl <- sab.fish %>% dplyr::filter(year == 2022 , months %in% psl.months)
sab.psl.total <- sum(sab.psl$pro.repwt,na.rm=T)


# So we can now make our decision table
# BMSY is somewhere between 4380 (6%) and 5030 (5%) tonnes, B0 is a shade under 7000 tonnes
# So then we make the USR be 80% of the BMSY/TRP... and the USR can be 30 or 40% of that.


# Reference Points
TRP <- NULL # Using the B0 reference points
USR <- NULL # Typical USR's 
LRP <- NULL # So 20% of B0 
RR <- NULL # THis is in %
#RR.TRP <- 5
n.sims <- 1e6
g.adj <- gR.adj <- 0
catchs <- seq(0,200,by=10)




# Interesting, in the 'no growth' scenario, for the FR biomass we only have natural mortality acting on it, so
# taking out the biomass early means less dies and thus the biomass change from removing X tonnes is actually less than the X tonnes

sab.dt <- dec.tab(mod.select = "SEAM",data = sab.mod, catch.scenarios = catchs,PSL = sab.psl.total,
                  n.sims = n.sims,TRP = TRP,USR = USR,LRP = LRP, RR = RR,RR.TRP = NULL,g.adj=g.adj,gR.adj=gR.adj)

if(!dir.exists(paste0(mod.loc,"Results/Sab/R_75_FR_90/SEAM_",seam.select,"/"))) dir.create(paste0(mod.loc,"Results/Sab/R_75_FR_90/SEAM_",seam.select,"/"))

if(!dir.exists(paste0(mod.loc,"Results/Sab/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/"))) dir.create(paste0(mod.loc,"Results/Sab/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/"))

saveRDS(sab.dt,paste0(mod.loc,"Results/Sab/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,"_g_adj_",g.adj,"_gR_adj_",gR.adj,".Rds"))
write_csv(sab.dt$table,paste0(mod.loc,"Results/Sab/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,"_g_adj_",g.adj,"_gR_adj_",gR.adj,".csv"))

# Turn the decision table into something nicely formated for Word 
tab <- kableExtra::kbl(sab.dt$table, booktabs = TRUE, escape =F, format = 'pipe',align = c('l','l','l','r','l'))#,
#caption = cap) %>%
#kable_styling(full_width = F) %>% row_spec(c(2:10,12,14:18,20), bold = T) %>%
#kable_styling(full_width = F) %>% row_spec(c(2,4,5,8:12,14,17,20), italic = T) %>%
#add_footnote(notation = 'number',ft.note,escape=F)
saveRDS(tab,paste0(mod.loc,"Results/Sab/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/Sab_word_ready_DT_TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,"_g_adj_",g.adj,"_gR_adj_",gR.adj,".Rds"))

# pdf version of the same
tab.pdf <- kableExtra::kbl(sab.dt$table, booktabs = TRUE, escape =F, format='latex',align = c('l','l','l','r','l'))#,
#caption = cap) %>%
#kable_styling(full_width = F) %>% row_spec(c(2:10,12,14:18,20), bold = T) %>%
#kable_styling(full_width = F) %>% row_spec(c(2,4,5,8:12,14,17,20), italic = T) %>%
#kable_styling(latex_options = c("hold_position","scale_down")) %>%
#add_footnote(notation = 'number',ft.note,escape=F)
saveRDS(tab.pdf,paste0(mod.loc,"Results/Sab/R_75_FR_90/SEAM_",seam.select,"/Decision_Table/Sab_pdf_ready_DT_TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,"_g_adj_",g.adj,"_gR_adj_",gR.adj,".Rds"))
