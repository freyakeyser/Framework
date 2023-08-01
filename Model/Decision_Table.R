# Here are the BBn and Sable Decision Tables.

# Here's where I'm sticking the figures
mod.loc <- "D:/Framework/SFA_25_26_2024/Model/"
fun.loc <-"D:/Github/Framework/RefPts/functions/"
rp.loc <- "D:/Framework/SFA_25_26_2024/RefPts/"

# Functions we need to source for this to all work
source(paste0(fun.loc,"Decision_Table_function.R"))


#seam.select <- "1994_2022_vary_m_m0_1_qR_0.5_15_knots" # This has the old von B growth parameters...
#seam.select <- "1994_2022_vary_m_m0_1_qR_0.5_10_knots"
seam.select <- "1994_2022_vary_m_m0_2_qR_0.5_10_knots"
seam.mod <- readRDS(paste0(mod.loc,"Results/BBn/R_75_FR_90/BBn_SEAM_model_output_",seam.select,".Rds"))


# So now we can run a model using the B/R/m uncertainty to predict biomass next year under various catch scenarios.
catchs <- seq(0,750,by=50)
# Reference Points
TRP <- 8000
USR <- 4800
LRP <- 2400
RR <- 14 # THis is in %
n.sims <- 1e6

dt <- dec.tab(mod.select = "SEAM",data = seam.mod, catch.scenarios = catchs,n.sims = n.sims,TRP = TRP,USR = USR,LRP = LRP, RR = RR)
saveRDS(dt,paste0(rp.loc,"Results/BBn/R_75_FR_90/",seam.select,"/BBn_decision_table_TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,".Rds"))

# And a Sable decision table 
seam.select <- "1995_2022_vary_m_m0_0.4_qR_0.5_4_knots"
seam.mod <- readRDS(paste0(mod.loc,"Results/Sab/R_75_FR_90/Sab_SEAM_model_output_",seam.select,".Rds"))

# Reference Points
TRP <- 4400
USR <- 2200
LRP <- 1100
RR <- 4 # THis is in %
n.sims <- 1e6
catchs <- seq(0,150,by=10)


dt <- dec.tab(mod.select = "SEAM",data = seam.mod, catch.scenarios = catchs,n.sims = n.sims,TRP = TRP,USR = USR,LRP = LRP, RR = RR)
saveRDS(dt,paste0(rp.loc,"Results/Sab/R_75_FR_90/",seam.select,"/Sab_decision_table_TRP_",TRP,"_USR_",USR,"_LRP_",LRP,"_RR_",RR,".Rds"))
