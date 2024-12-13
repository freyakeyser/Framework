# This is going to be an overly complex function in which we can set up various productivity modules to project the population into the future.
# Will want harvest control rules in this too eventually I think

# Building the options, the suite of options are the same for each of the 'mods'. You look at the relationships betwen each productivity parameter and your biomass to determine what is most appropriate
# run:  For the moment if you put 'model_error' it grabs the error estimates from the model for Rec m, and the process error term.  If you put anything else here there is no process error 
#       in the model and the Rec/m uncertainty is coming from the data itself.  I'm still unsure how I want to handle this...
#g.mod: What is our growth model.  
      #a:  type: Basically the options are a linear model looking at RPS~SSB, or a breakpoint analysis that splits the RPS into a low and high period.  If the relationship appears to be
      #          density independent, you can do a 'bp' and set the bp to be 0.  Options are 'lm', 'glm', 'gam', 'and 'bp'.  Default is 'bp'. Also a 'cor' option to make it a correlated TS.
      #          added a new option to set the growth rate at a value you want (e.g. the median of the time series). This method keeps the gR raio I used, so if gR will be
      #          larger than gR by whatever is generally observed in the data (is prod.dat$g.ratio). You can also set g and gR, which effectively overwrites this option and just sets g and/or gR

      #b: bp.strategy: This is only used if type = 'bp'.  Options are 'sample' where we just sample growth from the historically observed growth data and 'dist' where we select growth from a 
      #                normal distribution with mean/variance observed in data  Default is 'sample'
      #c: B.rec.bp : Is there a biomass breakpoint where growth appears to be lower/higher on either side.  The model will select growth rates appropriate for the current biomass based on the breakpoint.
          # Default sets this to 0 which effectively is no breakpoint.
#rec.mod: What is our recruitment model.  This currently has a few options, pretty much infinitely more options are possible to add in here.
      #a: type: Options right now are a linear model "lm", which will effectively be a Ricker model linearized, so log(RPS) ~ SSB model using the productivity data. Also a 'cor' option to make it a correlated TS.
      #          added a new option to set the growth rate at a value you want (e.g. the median of the time series)
          # 'bp' 
#m.mod: What is our natural mortality model.  
#mod.run The model run that your projections are coming from. 
#model.yrs: The years run in the model, excluding the projection year. Defaulting to 1994:2022 which is what we did for BBn and Sable.
#base.yrs: Gives us the option to truncate the data used in the analyses to a subset of years.  Lets us explore different 'productivity regimes'.  Default is NULL which uses all the data.
#          If using this then you'd want it something like base.yrs = 2010:2022
# for type I am adding in a  option in case you want to set something at a value. This is being done for growth for now, but could be useful for anything.  You need to specify what the
# set.value is as well, default of set.value is NULL

proj.mod <- function(mods = list(tlm.mod = NA,seam.mod = seam.mod,bssm.mod=NA), n_sim = 300,exp.scenario=seq(0,0.4,0.025),n_y = 100, LRP = 1000, HCR.sim = NULL, save.results = F,
                     ci.proj = data.frame(lci=0.25,uci=0.75),run = 'model_error', base.yrs = NULL, model.yrs = 1994:2022,
                     plot_url = "D:/Github/BBn_model/Results/Figures/BBn/",
                     res_url = "D:/Github/BBn_model/Results/Models/BBn/",
                     model = "SEAM", max.SSB =NULL, 
                     mod.run = NULL,
                     # He we have our productivity models
                     g.mod =   list(type = 'bp',set.value = "off",bp.strategy = 'sample',bp= 0,mn.at.max = NULL,sd.at.max = 0.05,ar1=0,ar2=0,
                                    set.g = NULL, set.gR = NULL,cen.tend = 'mean'), 
                     rec.mod = list(type = 'bp',set.value = "off",bp.strategy = 'sample',bp= 0,mn.at.max = NULL,sd.at.max = 0.01,rec.age = 5,
                                    ar1=0,ar2=0,rps.m.cor = 0,rps.m.lag=0,cen.tend = 'mean'), # the rps.m.cor and rps.m.lag are only relevant if both rec.mod and m.mod are set as 'cor'),
                     m.mod =   list(type = 'bp',set.value = "off",bp.strategy = 'sample',bp= 0,mn.at.max = NULL,sd.at.max = 0.05,ar1=0,ar2=0,cen.tend = 'mean')
                     )
                     

{

library(mgcv)
u.colors <- c("#FFD500","#005BBB")
  
theme_set(theme_few(base_size = 18))
  # Now get the plots and results folders sorted out...
# Create the directories you want
if(!dir.exists(paste0(plot_url))) dir.create(paste0(plot_url))
if(!dir.exists(paste0(res_url))) dir.create(paste0(res_url))

# Now we can create the name of the HCR data object I want if we are HCR'ing.
if(!is.null(HCR.sim))  
{
  HCR <- paste0("T_",HCR.sim$TRP,"_",HCR.sim$TRP.exp,"_U_",HCR.sim$USR,"_",HCR.sim$USR.exp,"_L_",HCR.sim$LRP,"_",HCR.sim$LRP.exp)
  if(!dir.exists(paste0(plot_url,"HCR/",HCR,"/"))) dir.create(paste0(plot_url,"HCR/",HCR,"/"))
  if(!dir.exists(paste0(res_url,"HCR/",HCR,"/"))) dir.create(paste0(res_url,"HCR/"))
  res_sims <- paste0(res_url,"HCR/",HCR,".Rds")
  plot_sims <- paste0(plot_url,"HCR/",HCR,"/")
}

# Names and locations of save objects if we are running different combos of years, only if not doing the HCR thing.
if(is.null(HCR.sim))  
{
  if(!is.null(base.yrs)) 
  {
    res_sims <- paste0(res_url,min(base.yrs),"-",max(base.yrs),".Rds")
    plot_sims <- paste0(plot_url,min(base.yrs),"-",max(base.yrs),"/")
  }
  if(is.null(base.yrs)) 
  {
    res_sims <- paste0(res_url,min(model.yrs),"-",max(model.yrs),".Rds")
    plot_sims <- paste0(plot_url,min(model.yrs),"-",max(model.yrs),"/")
  }
} # end  if(is.null(HCR.sim))  

# Pick the appropriate data depending on your model.
if(model == "TLM") mod.fit <- mods$tlm.mod
if(model == "SEAM") mod.fit <- mods$seam.mod  
if(model == "BSSM") mod.fit <- mods$bssm.mod
  
  # If doing the process error run you need to do this, 
  #DK note: I've turned off the process error term in the model, it just results in the model producing biomass out of nothing, which I don't think is good behaviour for projections.
  if(run == 'model_error' & model == "TLM")
  {
    sigma_tau <- 0 #mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'sigma_tau'] #Biomass/process error in model
    sigma_phi <- mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'sigma_phi'] # Recruits
    sigma_m   <- mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'sigma_m'] # Natural mort
  }
  
  if(run == 'model_error' & model == "SEAM")
  {
    sigma_tau <- 0 #mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'SigmaO_B'] #Biomass
    sigma_phi <- mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'SigmaO_R'] # Recruits
    sigma_m   <- mod.fit$sdrep$value[names(mod.fit$sdrep$value) == "SigmaO_m"] # Natural mort
  }

  # This one would require some more thought if we ran with model error included.  
  if(run == 'model_error' & model == "BSSM")
  {
    sigma_tau <- 0 #mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'SigmaO_B'] #Biomass
    sigma_phi <- mean(1/mod.fit$mean$IR.precision)^0.5 # Note there is one of these every year for BSSM, while SEAM/TLM just estimate this once
    sigma_m   <- 0 # There isn't really anything for natural mortality as the distribution is on the clappers not on M.
  }


  #  If we aren't doing the process error thing then set sigma_tau to 0.  
  # For the recruit and natural mortality no error models we do something different and is handled down in the code below.
  if(run != 'model_error') sigma_tau <- 0
  
  
  # Number of fishing scenarios we are running
  if(is.null(HCR.sim)) n.f.scen <- length(exp.scenario)
  if(!is.null(HCR.sim)) n.f.scen <- 1
  # Using the recruitment inline below based on 
  #Rec<-median(mod.fit$report$totR) # Median recruitment
  #Rec<-min(mod.fit$report$totR) # Min recruitment
  # Setting up the output arrays
  Bio<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  Catch<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  mort<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  mort.r<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  g<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  gR<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  Rec <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  Surp.prod <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  r.no.fishing <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  r.realized <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  f.table <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  rps <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  m.cor.ts <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  rps.cor.ts <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  
  Exp.res <- NULL
  
  # Finally I make this dataframe to make our lives easier where I combine the main productivity parameters we need.
  # We want to scrub off the most recent year as that is a projection year.
  if(model != "BSSM")
  {
  prod.dat <- data.frame(B = mod.fit$report$totB[-length(mod.fit$report$totB)],
                          R = mod.fit$report$totR[-length(mod.fit$report$totR)], # For some reason the SEBDAM run needs a NA at the end here.
                          m = mod.fit$report$mean_m[-length(mod.fit$report$mean_m)],
                          g = as.vector(mod.fit$obj$env$data$g[-length(mod.fit$obj$env$data$g)]),
                          gR = as.vector(mod.fit$obj$env$data$gR[-length(mod.fit$obj$env$data$g)]),
                          year = model.yrs)
  }# end if(mod != "BSSM")
  
  if(model == "BSSM")
  { 
    #browser()
    # # For Sable we want the g values for this to be 1 but want to use the gR values observed.
    if(!is.null(g.mod$set.g)) mod.fit$data$g <- g.mod$set.g
    if(!is.null(g.mod$set.gR)) mod.fit$data$gR <- g.mod$set.gR
    
    prod.dat <- data.frame(year = model.yrs,
                           B = mod.fit$median$B,
                           R = mod.fit$median$R,
                           m = mod.fit$median$m,
                           mr = mod.fit$median$mR,
                           g = mod.fit$data$g, 
                           gR = mod.fit$data$gR)
  } # end if(mod == "BSSM")
  
  # Set this up so we can limit recruit biomass to not go into wild outlier land...
  max.rec.obs <- max(prod.dat$R,na.rm=T)
  # If modelling just a few years trim the data here.
  if(!is.null(base.yrs)) prod.dat <- prod.dat %>% dplyr::filter(year %in% base.yrs)
  
  # We can set the values for the parameters here...
  
  # A simple way to get gR estimate from the g estimate
  prod.dat$g.ratio <- prod.dat$gR/prod.dat$g
   # For the recruitment analysis we can get some preliminaries out of the way.
  # Using the Biomass of Recruits + Fully Recruited Scallop as our estimate for SSB.
  prod.dat$SSB <- prod.dat$B + prod.dat$R
  # Getting a vector of recruits offset so we can use that for our RPS analysis
  # First calculate Recruits per spawner based on the Recruitment Age chosen
  #browser()
  prod.dat$Recs <- c(prod.dat$R[(rec.mod$rec.age+1):nrow(prod.dat)],rep(NA,rec.mod$rec.age))
  # calculate recruits per spawner
  prod.dat$RPS <-  prod.dat$Recs/(prod.dat$SSB)
  min.SSB  <- ceiling(min(prod.dat$SSB/100,na.rm=T)) * 100 # This will round up to the 100's of tonnes, below this the RPS will be based on the observed RPS values at low biomass
  # We are going to use max.SSB to reduce productivity when the population is at a level above anything observed historically, or should we use max SSB above... 
  if(is.null(max.SSB)) max.SSB <- floor(max(prod.dat$SSB/100,na.rm=T)) * 100

  # If we are using some form a linear density dependence we run these here to the prediction object we will need to use to get a prediction of the response variable
  if(rec.mod$type %in% c('lm','glm','gam')) r.res <- dens.function(dat = prod.dat,type = rec.mod$type,min.cov.pred = min.SSB,max.cov.pred= max.SSB,step = 100,response = 'RPS',dd.term = "SSB",log.trans=T)
  if(g.mod$type %in% c('lm','glm','gam'))   g.res <- dens.function(dat = prod.dat,type = g.mod$type,min.cov.pred = min.SSB,max.cov.pred= max.SSB,step = 100,response = 'g',dd.term = "SSB",log.trans=T)
  #if(g.mod$type %in% c('lm','glm','gam'))   gr.res <- dens.function(dat = prod.dat,type = g.mod$type,min.cov.pred = min.SSB,max.cov.pred= max.SSB,step = 100,response = 'gR',dd.term = "SSB",log.trans=T)
  if(m.mod$type %in% c('lm','glm','gam'))   m.res <- dens.function(dat = prod.dat,type = m.mod$type,min.cov.pred = min.SSB,max.cov.pred= max.SSB,step = 100,response = 'm',dd.term = "SSB",log.trans=T)
 
  for(i in 1:n.f.scen) 
  {
    
    Sim.res <- NULL # reset the Sim res object, it'll just be a temporary container of the sims
    for(nn in 1:n_sim)
    {  
      # If we are running the correlation analysis for any of these, we do that outside the j loop as we can 'make' a whole time series in one shot, of course the
      # values in that time series will need to be 'replaced' if the biomass gets too high, but we'll do that below just like how we do it with the lm models.
      # Finally we can run a stand alone correlation analysis for g.mod.
      if(g.mod$type == 'cor')
      {
        # Because this is done on log scale (needed for RPS and nm really), I'm going to force the variance to be a bit artificially low here, this seems to give very reasonable growths
        if(g.mod$ar1 == 0 & g.mod$ar2 == 0)
        {
          if(g.mod$cen.tend == 'median') g.val <- median(prod.dat$g[prod.dat$g < max(g.mod$bp)],na.rm=T)
          if(g.mod$cen.tend == 'mean')  g.val <- mean(prod.dat$g[prod.dat$g < max(g.mod$bp)],na.rm=T)
          g.ts <-  cor.fun(ts1 = list(mn = g.ct,var=sd(log(prod.dat$g[prod.dat$g < max(g.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$g[prod.dat$g < max(g.mod$bp)]))),
                           ts2 = list(mn =0.1,var=0.1),
                           arima = list(proc = 'normal',cross.cor = 0,lag = 0),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
          
        } else {
          g.ts <-  cor.fun(ts1 = list(mn = g.val,var=sd(log(prod.dat$g[prod.dat$g < max(g.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$g[prod.dat$g < max(g.mod$bp)]))),
                          ts2 = list(mn =0.1,var=0.1),
                          arima = list(proc = 'arima',cross.cor = 0,lag = 0,ts1.ar1= g.mod$ar1,ts1.ar2=g.mod$ar2),
                          ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
        } # end else
      } #end g.mod$type == 'cor'
      
      # Now do the same thing for the recruits, but I also have a version where we allow rps and m to covary together, so this is the version where we only have recs doing their thing
      if(rec.mod$type == 'cor' & rec.mod$rps.m.cor ==0)
      {
        if(rec.mod$cen.tend == 'median') rps.val <- median(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)],na.rm=T)
        if(rec.mod$cen.tend == 'mean')  rps.val <- mean(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)],na.rm=T)
        # Because this is done on log scale (needed for RPS and nm really), I'm going to force the variance to be a bit artificially low here, this seems to give very reasonable growths
        if(rec.mod$ar1 == 0 & rec.mod$ar2 == 0)
        {
          rps.ts <-  cor.fun(ts1 = list(mn = rps.val,var=sd(log(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)]))),
                           ts2 = list(mn =0.1,var=0.1),
                           arima = list(proc = 'normal',cross.cor = 0,lag = 0),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts  
        } else {
          rps.ts <-  cor.fun(ts1 =list(mn = rps.val,var=sd(log(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)]))),
                           ts2 = list(mn =0.1,var=0.1),
                           arima = list(proc = 'arima',cross.cor = 0,lag = 0,ts1.ar1= rec.mod$ar1,ts1.ar2=rec.mod$ar2),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
        } # end else
      } #end rec.mod$type == 'cor'
      
      # Now do the same thing for the natural mortality, but I also have a version where we allow rps and m to covary together, this versoin is no cross-correlation
      if(m.mod$type == 'cor' & rec.mod$rps.m.cor ==0)
      {
        if(m.mod$cen.tend == 'median') m.val <- median(prod.dat$m[prod.dat$m < max(m.mod$bp)],na.rm=T)
        if(m.mod$cen.tend == 'mean')  m.val <- mean(prod.dat$m[prod.dat$m <max(m.mod$bp)],na.rm=T)
        # Because this is done on log scale (needed for RPS and nm really), I'm going to force the variance to be a bit artificially low here, this seems to give very reasonable growths
        if(m.mod$ar1 == 0 & m.mod$ar2 == 0)
        {
          m.ts <-  cor.fun(ts1 = list(mn = m.val,var=sd(log(prod.dat$m[prod.dat$m < max(m.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$m[prod.dat$m < max(m.mod$bp)]))),
                           ts2 = list(mn =0.1,var=0.1),
                           arima = list(proc = 'normal',cross.cor = 0,lag = 0),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
        } else {
          m.ts <-  cor.fun(ts1 =list(mn = m.val,var=sd(log(prod.dat$m[prod.dat$m < max(m.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$m[prod.dat$m < max(m.mod$bp)]))),
                             ts2 = list(mn =0.1,var=0.1),
                             arima = list(proc = 'arima',cross.cor = 0,lag = 0,ts1.ar1= m.mod$ar1,ts1.ar2=m.mod$ar2),
                             ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
        } # end else
      } #end m.mod$type == 'cor' 
      
      if((m.mod$type == 'cor' | rec.mod$type == 'cor') & rec.mod$rps.m.cor !=0)
      {
        if(m.mod$cen.tend == 'median') m.val <- median(prod.dat$m[prod.dat$m < max(m.mod$bp)],na.rm=T)
        if(m.mod$cen.tend == 'mean')  m.val <- mean(prod.dat$m[prod.dat$m < max(m.mod$bp)],na.rm=T)
        if(rec.mod$cen.tend == 'median') rps.val <- median(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)],na.rm=T)
        if(rec.mod$cen.tend == 'mean')  rps.val <- mean(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)],na.rm=T)
        
        # Because this is done on log scale (needed for RPS and nm really), I'm going to force the variance to be a bit artificially low here, this seems to give very reasonable growths
        if(m.mod$ar1 == 0 & m.mod$ar2 == 0 & rec.mod$ar1 == 0 & rec.mod$ar2 == 0)
        {
          tmp.res <-  cor.fun(ts1 = list(mn = m.val,var=sd(log(prod.dat$m[prod.dat$m < max(m.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$m[prod.dat$m < max(m.mod$bp)]))),
                           ts2 =  list(mn = rps.val,var=sd(log(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)]))),
                           arima = list(proc = 'normal',cross.cor = rec.mod$rps.m.cor,lag = rec.mod$rps.m.lag),
                           ts =list(years=500,start.year =300, final.ts.len = n_y+10))
          m.ts <- tmp.res$ts1.ts
          rps.ts <- tmp.res$ts2.ts
        } else {
          tmp.res <-  cor.fun(ts1 =list(mn = m.val,var=sd(log(prod.dat$m[prod.dat$m < max(m.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$m[prod.dat$m < max(m.mod$bp)]))),
                           ts2 =  list(mn = rps.val,var=sd(log(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)]),na.rm=T)/sqrt(length(prod.dat$RPS[prod.dat$RPS < max(rec.mod$bp)]))),
                           arima = list(proc = 'arima',cross.cor = rec.mod$rps.m.cor,lag = rec.mod$rps.m.lag,
                                        ts1.ar1= m.mod$ar1,ts1.ar2=m.mod$ar2,ts2.ar1= rec.mod$ar1,ts2.ar2=rec.mod$ar2),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))
          m.ts <- tmp.res$ts1.ts
          rps.ts <- tmp.res$ts2.ts
        } # end else
        #browser()
        m.cor.ts[,nn,i] <-   m.ts[1:n_y]
        rps.cor.ts[,nn,i] <- rps.ts[1:n_y]
      } #end m.mod$type == 'cor'
      
      #browser()
      for (j in 1:n_y)
      {
        #browser()
        # Get a tidy biomass and Recruit value to save on code below...
        # For the first year we take the biomass from last year as our starting point
        if (j == 1) 
        {
          # Last year of real data...
          Bio[j,nn,i] <- B.init <- prod.dat$B[nrow(prod.dat)] 
          Rec[j,nn,i] <- Rec.init <- prod.dat$R[nrow(prod.dat)] 
        } # end  if (j == 1) 
        # After the first year we have last years results.
        if (j > 1) 
        {
          B.init <- Bio[j,nn,i]
          Rec.init <- Rec[j-1,nn,i]
        } # end if (j > 1)
        
        # Get our quasi-SSB from last year, this is what we'll have to use to parameterize our density terms.
        # DK Note: Be congnizant of one year offset this causes, this is saying last year SSB influences the upcoming DD process.
        # I don't think there's much else one can do, but that is subtly different than what we are doing with the data exploration which uses growth/rec, etc from the 'same' year.
        SSB.init <- B.init + Rec.init
        ##################################################
        ### Growth Section ###
        # browser()
        # First up commercial sized growth
        # using the breakpoint method determine what growth is.  The BP method wraps 
        # all scenarios into one tidy bit of code.  Also, note that by default this includes all data because the default for bp.g.B = 0
        #browser()
        if(g.mod$type == 'bp') g[j,nn,i] <- bp.function(B = SSB.init,dat = prod.dat$g,ssb.ts=prod.dat$SSB,type = g.mod$bp.strategy, bp = g.mod$bp,cen.tend = g.mod$cen.tend)#,
                                                       #max.B = max.SSB,mn.at.max = g.mod$mn.at.max,sd.at.max = g.mod$sd.at.max)
        #max.B = max.SSB,mn.at.max = g.mod$mn.at.max,sd.at.max = g.mod$sd.at.max)
        
        #if(g.mod$type == 'sv') {g[j,nn,i] <- gR[j,nn,i] <- g.mod$set.value}
        
        # Now if we are using a linear model to predict growth then we have to do this...
        if(g.mod$type %in% c("lm",'glm','gam'))
        {
          # If below the minimum SSB ever observed assume growth is what was observed at the lowest SSB in the time series. Do not extrapolate...
          # This isn't terrible for a log normal, the huge outliers were so outlier that this doesn't really go far enough on BBn
          if(SSB.init < min.SSB) g[j,nn,i] <- rlnorm(1,g.res$pred.dat$g.log[1],g.res$pred.dat$se[1])
          # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
          # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
          #if(SSB.init > max.SSB) g[j,nn,i] <- rlnorm(1,log(g.mod$mn.at.max),g.mod$sd.at.max) 
          # Finally the most complicated one is to pick from the predictions and get a recruitment estimate based on the linear model predictions and uncertainty
          # We do this when the SSB is within observed bounds.
          if(SSB.init >= min.SSB)
          {
            pick <- which(g.res$pred.dat$SSB == round(SSB.init/100)* 100)
            g[j,nn,i] <- rlnorm(1,g.res$pred.dat$g.log[pick],g.res$pred.dat$se[pick]) 
          } # end if(SSB.init >= min.SSB & SSB.init <= max.SSB)
        } # end if(rec.mod$type %in% c("lm",'glm','gam'))
        
        if(g.mod$type == 'cor')
        {
          # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
          # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
          #if(SSB.init > max.SSB) g[j,nn,i] <- rlnorm(1,log(g.mod$mn.at.max),g.mod$sd.at.max) 
          # If above the breakpoint sample from the data above the breakpoint, this will slightly weaken the correlation
 
          
          if(SSB.init > g.mod$bp)  
          {
            if(g.mod$cen.tend == 'median') g.val <- median(prod.dat$g[prod.dat$SSB > g.mod$bp],na.rm=T)
            if(g.mod$cen.tend == 'mean')  g.val <- mean(prod.dat$g[prod.dat$SSB > g.mod$bp],na.rm=T)
            g[j,nn,i] <- rlnorm(1,log(g.val),sd(prod.dat$g[prod.dat$SSB > g.mod$bp],na.rm=T))
          }
          # For the rest grab the data from the correlation
          if(SSB.init <= g.mod$bp)  g[j,nn,i] <- g.ts[j]
        } #end if g.mod$type == 'cor'
             
        # We can see that growth differences really doesn't vary too much and that recruit growth is always larger (basically between 10 and 25% larger)
        #if(g.mod$type != 'sv') gR[j,nn,i] <-  sample(prod.dat$g.ratio,1) * g[j,nn,i] #sample(prod.dat$gR,1) # * g[j,nn,i]
        # Whenever the SSB is above the maximum observed we can trigger a clause to constrain growth
        if(SSB.init > max.SSB & !is.null(g.mod$mn.at.max)) 
        {
          g[j,nn,i] <- rlnorm(1,log(g.mod$mn.at.max),g.mod$sd.at.max) 
          #gR[j,nn,i] <- rlnorm(1,log(g.mod$mn.at.max),g.mod$sd.at.max) 
        }
        # And get a gR using the g ratio, do this if we don't have g set...
        if(is.null(g.mod$set.g)) gR[j,nn,i] <- g[j,nn,i]*sample(prod.dat$g.ratio,1)
        # If we've set g, get the g estiamtes by log-normally the g ratio rather than sampling (because that only gives us a few values for gR)
        if(!is.null(g.mod$set.g)) 
        {
          gR[j,nn,i] <- g[j,nn,i]*rnorm(1,prod.dat$g.ratio,sd(prod.dat$g.ratio,na.rm=T))
          # Make sure this stays positive (since using rnorm), if it doesn't just set it to the median value... it's like a 1 in a trillion^trillion thing...
          if(gR[j,nn,i] < 0) gR[j,nn,i] <- g[j,nn,i] * median(prod.dat$g.ratio,na.rm=T)
        }
        #browser()
        ################################# End Growth Section ##################################
        
       ############################### Start recruitment Section #######################################
        #browser()
        ### Recruit section
        # For the recruits I need to get the right biomass year to get the offset from recruits to the biomass that produced those recruits.  Easiest to do up here 
        if(j < (rec.mod$rec.age+1))  SSB.for.rec <- prod.dat$B[length(prod.dat$B)-rec.mod$rec.age+1] + prod.dat$R[length(prod.dat$R)-rec.mod$rec.age+1]
        if(j >=(rec.mod$rec.age+1)) SSB.for.rec <- Bio[j-rec.mod$rec.age,nn,i] + Rec[j-rec.mod$rec.age,nn,i]
        # Now we can get our recruit estimate...
        # Using the predictions from our linear model above...
        if(rec.mod$type %in% c("lm",'glm','gam'))
        {
          # If below the minimum SSB ever observed assume RPS is what was predicted at the lowest SSB in the time series. Do not extrapolate...
          # This isn't terrible for a log normal, the huge outliers were so outlier that this doesn't really go far enough on BBn
          if(SSB.for.rec < min.SSB) Rec[j,nn,i] <- SSB.for.rec*  rlnorm(1,r.res$pred.dat$RPS.log[1],r.res$pred.dat$se[1]) 
          # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
          # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
          #if(SSB.for.rec > max.SSB) Rec[j,nn,i] <- SSB.for.rec* rlnorm(1,log(rec.mod$mn.at.max),rec.mod$sd.at.max) 
          # Finally the most complicated one is to pick from the predictions and get a recruitment estimate based on the linear model predictions and uncertainty
          # We do this when the SSB is within observed bounds.
          if(SSB.for.rec >= min.SSB)# & SSB.for.rec <= max.SSB
          {
            pick <- which(r.res$pred.dat$SSB == round(SSB.for.rec/100)* 100)
            Rec[j,nn,i] <- rlnorm(1,r.res$pred.dat$RPS.log[pick],r.res$pred.dat$se[pick]) * SSB.for.rec
          } # end if(SSB.for.rec >= min.SSB & SSB.for.rec <= max.SSB)
        } # end if(rec.mod$type %in% c("lm",'glm','gam'))

        
        if(rec.mod$type %in% 'bp')
        {
          # DK been playing around with the best way too do this.
          if(!run %in% 'model_error') 
          {
            # Should I use the RPS or Recruitment time series for these??
             #if(SSB.for.rec >=min.SSB)  
            Rec[j,nn,i] <- SSB.for.rec * bp.function(B = SSB.for.rec,dat = prod.dat$RPS,ssb.ts=prod.dat$SSB,type = rec.mod$bp.strategy, bp = rec.mod$bp,cen.tend = rec.mod$cen.tend)#,
                                                                            #max.B = max.SSB,mn.at.max = rec.mod$mn.at.max,sd.at.max = rec.mod$sd.at.max)
            # To avoid huge outlier problem, if the recruitment estimate is higher than ever has been observed we instead run a sample
            # using the largest value and a nice normal distribution to avoid silliness
            if(Rec[j,nn,i] > max.rec.obs) Rec[j,nn,i] <- rnorm(1,max.rec.obs,0.1*max.rec.obs) # The max ends up about 50% above maximum ever seen, so reasonble.
            if(Rec[j,nn,i] < 0) Rec[j,nn,i] <- max.rec.obs # Just in case that trips and lands in negative land, very unlikely with our data (like < 1 in a billion)... but still...
          } # end if(run != 'model_error')     
          
          if(run == 'model_error') 
          {
            # Should I use the RPS or Recruitment time series for these??
            #if(SSB.for.rec >= min.SSB) 
              Rec[j,nn,i] <- SSB.for.rec * bp.function(B = SSB.for.rec,dat = prod.dat$RPS,ssb.ts=prod.dat$SSB,type = rec.mod$bp.strategy, bp = rec.mod$bp,sd = sigma_phi,cen.tend = rec.mod$cen.tend)#,
                                                                            #max.B = max.SSB,mn.at.max = rec.mod$mn.at.max,sd.at.max = rec.mod$sd.at.max)
              if(Rec[j,nn,i] > max.rec.obs) Rec[j,nn,i] <- rlnorm(1,log(max.rec.obs),0.1*max.rec.obs) # The max ends up about 50% above maximum ever seen, so reasonble.
              if(Rec[j,nn,i] < 0) Rec[j,nn,i] <- max.rec.obs # Just in case that trips and lands in negative land, very unlikely with our data (like < 1 in a billion)... but still...
              
            # If below the minimum obsreved biomass than recruitment will be based on the RPS at low biomass levels.
            #if(SSB.for.rec < min.SSB)  Rec[j,nn,i] <- SSB.for.rec * bp.function(B = SSB.for.rec,dat = prod.dat$RPS,type = rec.mod$bp.strategy, bp = rec.mod$bp,sd = sigma_phi,
            #                                                                          max.B = max.SSB,mn.at.max = rec.mod$mn.at.max,sd.at.max = rec.mod$sd.at.max)
          } # end if(run == 'model_error') 
           
         
        } # end if(rec.mod$type %in% 'bp')
     
         if(rec.mod$type == 'cor')
         {
           # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
           # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
           #if(SSB.for.rec > max.SSB) Rec[j,nn,i] <- SSB.for.rec* rlnorm(1,log(rec.mod$mn.at.max),rec.mod$sd.at.max) 
           # Do this differently for 1 or 2 breakpoints... 
           # This one breakpoint assumption is saying RPS above the breakpoint is normal, and below is the 'deviant'.
           if(length(rec.mod$bp) == 1)
           {
             
             # Above bp means we are using the correlation for that data
             if(SSB.for.rec < rec.mod$bp)    Rec[j,nn,i] <-  SSB.for.rec*rps.ts[j] #  & SSB.for.rec <= max.SSB
             # Below bp means we break the correlation.
             if(SSB.for.rec >= rec.mod$bp) 
             {
               if(rec.mod$cen.tend == 'median') rps.val <- median(prod.dat$RPS[prod.dat$SSB >= rec.mod$bp],na.rm=T)
               if(rec.mod$cen.tend == 'mean')  rps.val <- mean(prod.dat$RPS[prod.dat$SSB >= rec.mod$bp],na.rm=T)
               Rec[j,nn,i] <- SSB.for.rec* rlnorm(1,log(rps.val), sd(prod.dat$RPS[prod.dat$SSB >= rec.mod$bp],na.rm=T))
             }
           } # end if(length(rec.mod$bp) == 1)
           
          
           # If you have a couple of breakpoints
           if(length(rec.mod$bp) > 1)
           {
             if(SSB.for.rec < min(rec.mod$bp)) 
             {
               if(rec.mod$cen.tend == 'median') rps.val <- median(prod.dat$RPS[prod.dat$SSB < min(rec.mod$bp)],na.rm=T)
               if(rec.mod$cen.tend == 'mean')  rps.val <- mean(prod.dat$RPS[prod.dat$SSB < min(rec.mod$bp)],na.rm=T)
               Rec[j,nn,i] <-SSB.for.rec* rlnorm(1,log(rps.val),
                                                   sd(prod.dat$RPS[prod.dat$SSB < min(rec.mod$bp)],na.rm=T))
             } # end if(SSB.for.rec < min(rec.mod$bp)) 
             
             if(SSB.for.rec >= min(rec.mod$bp) & SSB.for.rec < max(rec.mod$bp)) # SSB.for.rec < max(rec.mod$bp) &
             {
               Rec[j,nn,i] <-SSB.for.rec*rps.ts[j]
             } # end if(SSB.for.rec >= max(rec.mod$bp)) 
             
             if(SSB.for.rec >= max(rec.mod$bp)) #  & SSB.for.rec <= max.SSB
             {
               if(rec.mod$cen.tend == 'median') rps.val <- median(prod.dat$RPS[prod.dat$SSB >= max(rec.mod$bp)],na.rm=T)
               if(rec.mod$cen.tend == 'mean')  rps.val <- mean(prod.dat$RPS[prod.dat$SSB >= max(rec.mod$bp)],na.rm=T)
               Rec[j,nn,i] <-SSB.for.rec* rlnorm(1,log(rps.val), sd(prod.dat$RPS[prod.dat$SSB >= max(rec.mod$bp)],na.rm=T))
             } # end if(SSB.for.rec >= max(rec.mod$bp) & SSB.for.rec <= max.SSB) 
           } # end if(length(rec.mod$bp) > 1)
         } # end if(rec.mod$type == 'cor')
        
        if(rec.mod$type == 'sv') Rec[j,nn,i] <- rec.mod$set.value
        
        # Whenever the SSB is above the maximum observed we can trigger a clause to constrain recruitment
        if(SSB.init > max.SSB & !is.null(rec.mod$mn.at.max)) 
        {
         Rec[j,nn,i] <- SSB.for.rec* rlnorm(1,log(rec.mod$mn.at.max),rec.mod$sd.at.max) 
        }
        
         # I'd like too be able to look at the rps time series too...
         rps[j,nn,i] <- Rec[j,nn,i]/SSB.for.rec
        
        ############################### End recruitment Section #######################################

       
        ############################### Start natural mortality Section #######################################
        # kinda a hybrid of the above two modules
        # Now if we are using a linear model to predict growth then we have to do this...
        if(m.mod$type %in% c("lm",'glm','gam'))
        {
          # If below the minimum SSB ever observed assume growth is what was observed at the lowest SSB in the time series. Do not extrapolate...
          # This isn't terrible for a log normal, the huge outliers were so outlier that this doesn't really go far enough on BBn
          if(SSB.init < min.SSB) mort[j,nn,i] <- rlnorm(1,m.res$pred.dat$m.log[1],m.res$pred.dat$se[1])
          # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
          # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
          #if(SSB.init > max.SSB) mort[j,nn,i] <- rlnorm(1,log(m.mod$mn.at.max),m.mod$sd.at.max) 
          # Finally the most complicated one is to pick from the predictions and get a recruitment estimate based on the linear model predictions and uncertainty
          # We do this when the SSB is within observed bounds.
          if(SSB.init >= min.SSB) # & SSB.init <= max.SSB
          {
            pick <- which(m.res$pred.dat$SSB == round(SSB.init/100)* 100)
            mort[j,nn,i] <- rlnorm(1,m.res$pred.dat$m.log[pick],m.res$pred.dat$se[pick]) 
          } # end if(SSB.for.rec >= min.SSB & SSB.for.rec <= max.SSB)
        } # end if(rec.mod$type %in% c("lm",'glm','gam'))
        #browser()
        if(m.mod$type %in% 'bp')
        {
          # 
          if(run != 'model_error') mort[j,nn,i] <- bp.function(B = B.init,dat = prod.dat$m, ssb.ts=prod.dat$SSB,
                                                               type = m.mod$bp.strategy, bp = m.mod$bp,cen.tend = m.mod$cen.tend)#,
                                                                              #max.B = max.SSB,mn.at.max = m.mod$mn.at.max,sd.at.max = m.mod$sd.at.max)
          if(run == 'model_error') mort[j,nn,i] <- bp.function(B = B.init,dat = prod.dat$m,ssb.ts=prod.dat$SSB,
                                                               type = m.mod$bp.strategy, bp = m.mod$bp,sd = sigma_m,cen.tend = m.mod$cen.tend)#,
                                                                              #max.B = max.SSB,mn.at.max = m.mod$mn.at.max,sd.at.max = m.mod$sd.at.max)
        } # end if(rec.mod$type %in% 'bp')    
         
         if(m.mod$type == 'cor')
         {
           # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
           # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
           #if(SSB.init > max.SSB) mort[j,nn,i] <- rlnorm(1,log(m.mod$mn.at.max),m.mod$sd.at.max) 
           # If above the breakpoint sample from the data above the breakpoint...  & SSB.init <= max.SSB
           
           if(SSB.init >= m.mod$bp) 
           {
             if(rec.mod$cen.tend == 'median') m.val <- median(prod.dat$m[prod.dat$SSB >= m.mod$bp],na.rm=T)
             if(rec.mod$cen.tend == 'mean')  m.val <- mean(prod.dat$m[prod.dat$SSB >= m.mod$bp],na.rm=T)
             mort[j,nn,i]  <- rlnorm(1,log(m.val),sd(log(prod.dat$m[prod.dat$SSB > m.mod$bp]),na.rm=T))
           }
           # For the rest grab the data from the correlation
           if(SSB.init < m.mod$bp)  mort[j,nn,i] <- m.ts[j]
         } # end if(m.mod$type == 'cor')
         
         if(m.mod$type =="sv") mort[j,nn,i] <- m.mod$set.value
        
         #browser()
         # For the recruit natural mortality we just sample from the distribution, this probably should be generalized for other scenarios, but it works
         # for what we found for BSSM.
         if(model == "BSSM") mort.r[j,nn,i] <- bp.function(B = SSB.init,ssb.ts=prod.dat$SSB,dat = prod.dat$mr,type = "dist", bp = 0,cen.tend = m.mod$cen.tend)#,
         # For the not BSSM model set the mort.r term to the same as the mort term as these are not differentiated.
         if(model != "BSSM") mort.r[j,nn,i] <- mort[j,nn,i]
        ############################### End natural mortality Section #######################################
        
         # Whenever the SSB is above the maximum observed we can trigger a clause to make natural mortality big again...
         if(SSB.init > max.SSB & !is.null(m.mod$mn.at.max)) 
         {
           mort[j,nn,i] <- mort.r[j,nn,i] <- rlnorm(1,log(m.mod$mn.at.max),m.mod$sd.at.max) 
         }
         
        ### Now if we are running the correlation on two of the above (that's the most we can do)
        #### Now we run the model and extract some metric of interest while we are at it.
        # Get our Catch based on the exp.scenario scenario and based on the previous years biomass
        # DK Note, when/how we remove the catch is slightly different from what we do in our projections, it's mostly about how we calculated F
        # In our current decision tables F would be calculated for the upcoming season using the Bio estimate, so I'm sticking in an 'F.table' for the 'F' we'd calculated in our decision tables.
         
         
        
        if(j == 1)
        {
          Bio[j,nn,i] <- B.init
          Rec[j,nn,i] <- Rec.init
        }
        if(j < n_y)
        {
          #browser()
          # Our projection forward
          if(is.null(HCR.sim)) Catch[j,nn,i]<-B.init*exp.scenario[i]
          
          if(!is.null(HCR.sim))
          {
            # LRP catch
            if(B.init < HCR.sim$LRP)  Catch[j,nn,i] <- rlnorm(1,log(HCR.sim$LRP.exp),HCR.sim$exp.sd) * B.init
            # Assume linear decline in Explotation to LRP in Cautious zone
            if(B.init > HCR.sim$LRP & B.init <= HCR.sim$USR)  Catch[j,nn,i] <- rlnorm(1,log(HCR.sim$USR.exp),HCR.sim$exp.sd)*B.init*((B.init-HCR.sim$LRP)/HCR.sim$LRP)
            # Harvest at USR explotation here
            if(B.init > HCR.sim$USR & B.init < HCR.sim$TRP)  Catch[j,nn,i] <- rlnorm(1,log(HCR.sim$USR.exp),HCR.sim$exp.sd) * B.init
            # Harvest at TRP explotation here
            if(B.init >= HCR.sim$TRP)  Catch[j,nn,i] <- rlnorm(1,log(HCR.sim$TRP.exp),HCR.sim$exp.sd)*B.init
          }
          #browser()
          # The growth term is the whole trick, if g > m, the populations grow like stink, I'm thinking we shouldn't account for adult growth in the projections
          # and just have it rely on growth of recruits.
          Bio[j+1,nn,i] <- rlnorm(1,log(exp(-mort[j,nn,i])*g[j,nn,i]*(B.init-Catch[j,nn,i])+Rec[j,nn,i]*exp(-mort.r[j,nn,i])*gR[j,nn,i]),sigma_tau)
          # The Surplus production this won't match Bio
          Surp.prod[j,nn,i]<- Bio[j+1,nn,i] - B.init + Catch[j,nn,i]
          # The realized rate of growth (r_fish)
          r.realized[j,nn,i] <- log(Bio[j+1,nn,i]/B.init)
          # The potential rate of growth (r), this method allows the catch to grow (or shrink) based on m and g.
          r.no.fishing[j,nn,i]<-   log((Bio[j+1,nn,i] + exp(-mort[j,nn,i])*g[j,nn,i]*Catch[j,nn,i])/B.init) #log(((exp(-mort[j,nn,i])*g[j,nn,i]*B.init)+Rec[j,nn,i]*exp(-mort[j,nn,i])*gR[j,nn,i]) /B.init)
          # Also calculate F using the end biomass result which is what we'd currently report as F in our decision tables
          # DK Note: Add in different ways of calculating F to see how different the F estimates might be.
          f.table[j+1,nn,i] <-  Catch[j,nn,i] / (Bio[j+1,nn,i] +  Catch[j,nn,i])
          #browser()
        } # end if(j < n_y)
       
      } # end the j loop through the years
      Sim.res[[nn]] <- data.frame(year = 1:n_y,Sim = rep(nn,n_y), B = Bio[,nn,i],Rec = Rec[,nn,i],rps = rps[,nn,i],
                                  g = g[,nn,i],gR = gR[,nn,i],M = mort[,nn,i],MR = mort.r[,nn,i],
                                 Catch = Catch[,nn,i], SP = Surp.prod[,nn,i],F.dec.table = f.table[,nn,i],
                                 r.real = r.realized[,nn,i],r.no.f = r.no.fishing[,nn,i])
      if(model == "BSSM") Sim.res[[nn]]$mr <- mort.r[,nn,i]
    } # end the n loop through the simulations
    Exp.res[[i]] <- data.frame(F.scenario = rep(exp.scenario[i],n_y*n_sim),do.call('rbind',Sim.res))
    
    print(paste("Explotation Scenario",i, "Completed"))
  } # end the i loop through the exp.scenario scenarios
  
Exp.res <-  do.call('rbind',Exp.res) 


# Save and plot the results if you want to
if(save.results == T)
{
  
  # save our results
  saveRDS(Exp.res,res_sims)
  
  # Now make the figures!
  # Set the alpha level for the lines on the figures showing the individual realizations.
  if(n_sim <= 10) alphs <- 0.5
  if(n_sim > 10) alphs <- 0.1
  if(n_sim > 100) alphs <- 1/255
  
  if(is.null(HCR.sim))
  {
    # This is summarizing the final half of the years in the simulations.
    MSY.summarized <- Exp.res %>% dplyr::filter(year >= n_y/2) %>% dplyr::group_by(F.scenario,year) %>% 
                                   dplyr::summarise(B.mn = median(B,na.rm=T),          B.U90 = quantile(B,probs=0.95,na.rm=T),           B.U50 = quantile(B,probs=0.75,na.rm=T),
                                                                                     B.L90 = quantile(B,probs=0.05,na.rm=T),           B.L50 = quantile(B,probs=0.25,na.rm=T),
                                                    Rec.mn = median(Rec,na.rm=T),      Rec.U90 = quantile(Rec,probs=0.95,na.rm=T),       Rec.U50 = quantile(Rec,probs=0.75,na.rm=T),
                                                                                     Rec.L90 = quantile(Rec,probs=0.05,na.rm=T),       Rec.L50 = quantile(Rec,probs=0.25,na.rm=T),
                                                    Catch.mn = median(Catch,na.rm=T),  Catch.U90 = quantile(Catch,probs=0.95,na.rm=T),   Catch.U50 = quantile(Catch,probs=0.75,na.rm=T),
                                                                                     Catch.L90 = quantile(Catch,probs=0.05,na.rm=T),   Catch.L50 = quantile(Catch,probs=0.25,na.rm=T),
                                                    M.mn = median(M,na.rm=T),          M.U90 = quantile(M,probs=0.95,na.rm=T),           M.U50 = quantile(M,probs=0.75,na.rm=T),
                                                                                       M.L90 = quantile(M,probs=0.05,na.rm=T),           M.L50 = quantile(M,probs=0.25,na.rm=T),
                                                    F.mn = median(F.dec.table,na.rm=T),F.U90 = quantile(F.dec.table,probs=0.95,na.rm=T), F.U50 = quantile(F.dec.table,probs=0.75,na.rm=T),
                                                                                       F.L90 = quantile(F.dec.table,probs=0.05,na.rm=T), F.L50 = quantile(F.dec.table,probs=0.25,na.rm=T),
                                                    SP.mn = median(SP,na.rm=T),        SP.U90 = quantile(SP,probs=0.95,na.rm=T),         SP.U50 = quantile(SP,probs=0.75,na.rm=T),
                                                                                     SP.L90 = quantile(SP,probs=0.05,na.rm=T),         SP.L50 = quantile(SP,probs=0.25,na.rm=T),
                                                    RPS.mn = median(rps,na.rm=T),      RPS.U90 = quantile(rps,probs=0.95,na.rm=T),       RPS.U50 = quantile(rps,probs=0.75,na.rm=T),
                                                                                      RPS.L90 = quantile(rps,probs=0.05,na.rm=T),       RPS.L50 = quantile(rps,probs=0.25,na.rm=T))
    
    Equilib.dat <-    MSY.summarized %>% dplyr::filter(year >= n_y/2) %>% dplyr::group_by(F.scenario) %>% 
                                         dplyr::summarise(B.mn = median(B.mn,na.rm=T),         B.U90 = median(B.U90,na.rm=T),        B.U50 = median(B.U50,na.rm=T),
                                                                                               B.L90 = median(B.L90,na.rm=T),        B.L50 = median(B.L50,na.rm=T),
                                                          Catch.mn = median(Catch.mn,na.rm=T), Catch.U90 = median(Catch.U90,na.rm=T),Catch.U50 = median(Catch.U50,na.rm=T),
                                                                                               Catch.L90 = median(Catch.L90,na.rm=T),Catch.L50 = median(Catch.L50,na.rm=T),
                                                          Rec.mn = median(Rec.mn,na.rm=T),     Rec.U90 = median(Rec.U90,na.rm=T),    Rec.U50 = median(Rec.U50,na.rm=T),
                                                                                               Rec.L90 = median(Rec.L90,na.rm=T),    Rec.L50 = median(Rec.L50,na.rm=T))
    
    # Get Biomass per Recruit and Yield per Recruit
    Equilib.dat$BPR <- Equilib.dat$B.mn/Equilib.dat$Rec.mn
    Equilib.dat$YPR <- Equilib.dat$Catch.mn/Equilib.dat$Rec.mn
    
    
  # Biomass
    bio_hists<-ggplot(data=Exp.res,aes(x=B/1000))+
                      geom_histogram(fill="gray",col="black")+
                      facet_wrap(~F.scenario,scales="free")+
                      ylab("Count")+xlab("Predicted Biomass (metric tonnes x 1000)") 
    
    ggsave(paste0(plot_sims,"Bio_hists.png"),plot=bio_hists,height=15,width=15)
    # Recruits
    rec_hists<-ggplot(data=Exp.res,aes(x=Rec/1000))+
                      geom_histogram(fill="gray",col="black")+
                      facet_wrap(~F.scenario,scales="free") +
                      ylab("Count")+xlab("Predicted Recruitment (metric tonnes x 1000)")
    
    ggsave(paste0(plot_sims,"Rec_hists.png"),plot=rec_hists,height=15,width=15)
    
    # Catch
    catch_hists<-ggplot(data=Exp.res,aes(x=Catch/1000))+
                        geom_histogram(fill="gray",col="black")+
                        facet_wrap(~F.scenario,scales="free")+
                        ylab("Count")+xlab("Commercial Landings (metric tonnes x 1000)") 
    
    ggsave(paste0(plot_sims,"Catch_hists.png"),plot=catch_hists,height=15,width=15)
    
    # Natural Mortality
    m_hists<-ggplot(data=Exp.res,aes(x=M))+
                    geom_histogram(fill="gray",col="black")+
                    facet_wrap(~F.scenario,scales="free")+
                    ylab("Count")+xlab("Predicted Natural Mortality") 
    
    ggsave(paste0(plot_sims,"M_hists.png"),plot=m_hists,height=15,width=15)
    
    # Fishing Mortality 
    F_hists<-ggplot(data=Exp.res,aes(x=F.dec.table))+
                    geom_histogram(fill="gray",col="black")+
                    facet_wrap(~F.scenario,scales="free")+
                    ylab("Count")+xlab("F from Decision Table") 
    
    ggsave(paste0(plot_sims,"F_hists.png"),plot=F_hists,height=15,width=15)
    
    # Surplus Production
    SP_hists<-ggplot(data=Exp.res,aes(x=SP/1000))+
                    geom_histogram(fill="gray",col="black")+
                    facet_wrap(~F.scenario,scales="free")+
                    ylab("Count")+xlab("Surplus Production (metric tonnes x 1000)") 
    
    ggsave(paste0(plot_sims,"SP_hists.png"),plot=SP_hists,height=15,width=15)
  
  
  # Now make plots of the realizations over time
  # First the Biomass
  B.real <- ggplot(Exp.res,aes(x=year,y=B/1000,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Predicted Biomass (metric tonnes x 1000)") + xlab("") + 
                  scale_color_viridis_b(alpha=alphs,end = 0.75) +
                   theme(legend.position = 'none')
  
  ggsave(paste0(plot_sims,"Biomass_realizations.png"),plot=B.real,height=15,width=15)
  # Then the Recruits
  
  R.real <- ggplot(Exp.res,aes(x=year,y=Rec/1000,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Predicted Recruitment (metric tonnes x 1000)") + xlab("") + 
                  scale_color_viridis_b(alpha=alphs,end = 0.75) +
                   theme(legend.position = 'none')
  ggsave(paste0(plot_sims,"Recruit_realizations.png"),plot=R.real,height=15,width=15)
  # Next the Catch
  C.real <- ggplot(Exp.res,aes(x=year,y=Catch,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Predicted Catch (metric tonnes)") + xlab("") + 
                  scale_color_viridis_b(alpha=alphs,end = 0.75) +
                   theme(legend.position = 'none')
  
  ggsave(paste0(plot_sims,"Catch_realizations.png"),plot=C.real,height=15,width=15)
  #Natural Mortality
  M.real <- ggplot(Exp.res,aes(x=year,y=M,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Natural Mortality (Instantaneous)") + xlab("") + 
                  scale_color_viridis_b(alpha=alphs,end = 0.75) +
                   theme(legend.position = 'none')
    
  ggsave(paste0(plot_sims,"M_realizations.png"),plot=M.real,height=15,width=15)
  
  #Fishing Mortality
  F.real <- ggplot(Exp.res,aes(x=year,y=F.dec.table,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Fishing Mortality (Proportional)") + xlab("") + 
                  scale_color_viridis_b(alpha=alphs,end = 0.75) +
                   theme(legend.position = 'none')
  
  ggsave(paste0(plot_sims,"F_realizations.png"),plot=F.real,height=15,width=15)
  
  # Surplus Production
  SP.real <- ggplot(Exp.res,aes(x=year,y=SP,group = Sim,color=Sim)) + 
                    geom_line() + facet_wrap(~F.scenario) +
                    ylab("Surplus Production (tonnes)") + xlab("") + 
                    scale_color_viridis_b(alpha=alphs,end = 0.75) +
                     theme(legend.position = 'none')
  ggsave(paste0(plot_sims,"SP_realizations.png"),plot=SP.real,height=15,width=15)
  
  
  
  
  ################# estimates with 95% CI bands
  #### Biomass
  B.ts <- ggplot(MSY.summarized,aes(x=year,y=B.mn/1000)) + geom_line() + facet_wrap(~F.scenario) +
                                                           geom_ribbon(aes(x=year,ymax = B.U90/1000,ymin=B.L90/1000),fill='firebrick2',alpha=0.2) + 
                                                           geom_ribbon(aes(x=year,ymax = B.U50/1000,ymin=B.L50/1000),fill='blue',alpha=0.2) +
                                                           xlab("") + ylab("Biomass (metric tonnes x 1000)")
  ggsave(paste0(plot_sims,"B_mn_ts.png"),plot=B.ts,height=15,width=15)
  # Recruits
  Rec.ts <- ggplot(MSY.summarized,aes(x=year,y=Rec.mn/1000)) + geom_line() + facet_wrap(~F.scenario) +
                                                               geom_ribbon(aes(x=year,ymax = Rec.U90/1000,ymin=Rec.L90/1000),fill='firebrick2',alpha=0.2) + 
                                                               geom_ribbon(aes(x=year,ymax = Rec.U50/1000,ymin=Rec.L50/1000),fill='blue',alpha=0.2) + 
                                                               xlab("") + ylab("Recruitment (metric tonnes x 1000)")
  ggsave(paste0(plot_sims,"R_mn_ts.png"),plot=Rec.ts,height=15,width=15)
  # Catch
  Catch.ts <- ggplot(MSY.summarized,aes(x=year,y=Catch.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                              geom_ribbon(aes(x=year,ymax = Catch.U90,ymin=Catch.L90),fill='firebrick2',alpha=0.2) + 
                                                              geom_ribbon(aes(x=year,ymax = Catch.U50,ymin=Catch.L50),fill='blue',alpha=0.2) + 
                                                              xlab("") + ylab("Catch (metric tonnes)")
  ggsave(paste0(plot_sims,"Catch_mn_ts.png"),plot=Catch.ts,height=15,width=15)
  # Natural mortality
  M.ts <- ggplot(MSY.summarized,aes(x=year,y=M.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                      geom_ribbon(aes(x=year,ymax = M.U90,ymin=M.L90),fill='firebrick2',alpha=0.2) + 
                                                      geom_ribbon(aes(x=year,ymax = M.U50,ymin=M.L50),fill='blue',alpha=0.2) + 
                                                       xlab("") + ylab("Natural mortality (instantaneous)")
  ggsave(paste0(plot_sims,"M_mn_ts.png"),plot=M.ts,height=15,width=15)
  # Fishing mortality
  F.ts <- ggplot(MSY.summarized,aes(x=year,y=F.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                      geom_ribbon(aes(x=year,ymax = F.U90,ymin=F.L90),fill='firebrick2',alpha=0.2) + 
                                                      geom_ribbon(aes(x=year,ymax = F.U50,ymin=F.L50),fill='blue',alpha=0.2) + 
                                                       xlab("") + ylab("Exploitation rate")
  ggsave(paste0(plot_sims,"F_mn_ts.png"),plot=F.ts,height=15,width=15)
  # Surplus Production
  SP.ts <- ggplot(MSY.summarized,aes(x=year,y=SP.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                        geom_ribbon(aes(x=year,ymax = SP.U90,ymin=SP.L90),fill='firebrick2',alpha=0.2) + 
                                                        geom_ribbon(aes(x=year,ymax = SP.U50,ymin=SP.L50),fill='blue',alpha=0.2) + 
                                                         xlab("") + ylab("Surplus Production (metric tonnes)")
  ggsave(paste0(plot_sims,"SP_mn_ts.png"),plot=SP.ts,height=15,width=15)
  
  
  ####################### Reference Points peice
  # Now for the money shot, we can make the Figures that can help inform removal reference and perhap LRP.
  #Lets take the mean of the last 25 years and use that as our steady state
  
  # Biomass Equilibrium
  B.equil.plt<- ggplot(Equilib.dat,aes(x=F.scenario,y=B.mn/1000)) + geom_point()+ geom_line()+
                                                               geom_ribbon(aes(x=F.scenario,ymin=B.L90/1000, ymax=B.U90/1000),alpha=0.05,fill='firebrick2') +
                                                               geom_ribbon(aes(x=F.scenario,ymin=B.L50/1000, ymax=B.U50/1000),alpha=0.05,fill='blue') +
                                                               ylab("Equilibrium Biomass (metric tonnes x 1000)") + xlab("Explotation Rate")
  ggsave(paste0(plot_sims,"B_equilibrum.png"),plot=B.equil.plt,height=8,width=12)
  # Catch Equilibrium
  C.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=Catch.mn)) + geom_point()+ geom_line()+
                                                                  geom_ribbon(aes(x=F.scenario,ymin=Catch.L90, ymax=Catch.U90),alpha=0.05,fill='firebrick2') +
                                                                  geom_ribbon(aes(x=F.scenario,ymin=Catch.L50, ymax=Catch.U50),alpha=0.05,fill='blue') +
                                                                  ylab("Equilibrium Catch (metric tonnes)")  + xlab("Explotation Rate")
  ggsave(paste0(plot_sims,"Catch_equilibrium.png"),plot=C.equil.plt,height=8,width=12)
  # Biomass per Recruit (Biomass)
  BPR.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=BPR)) + geom_point()+ geom_line()+
                                                               ylab("Biomass per Recruit (Biomass)")  + xlab("Explotation Rate")
  ggsave(paste0(plot_sims,"BPR_equilibrium.png"),plot=BPR.equil.plt,height=8,width=12)
  
  
  # Yield per Recruit (Biomass)
  YPR.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=YPR)) + geom_point()+ geom_line()+
                                                               ylab("Yield per Recruit (Biomass)")  + xlab("Explotation Rate (proportional)")
  ggsave(paste0(plot_sims,"YPR_equilibrium.png"),plot=YPR.equil.plt,height=8,width=12)
  
  } #end if is.null(HCR.sim)
  
  
  
  if(!is.null(HCR.sim)) 
  {
    # Summarize the data to make a nice time series plot.
    # Note that I want the mean catch here, because that's really what people are interested in for this...
    # This is summarizing the final 100 years of the simulations
    HCR.summarized <- Exp.res %>% dplyr::filter(year >= n_y/2) %>% dplyr::group_by(year) %>% 
                                  dplyr::summarise(B.mn = median(B,na.rm=T),         B.U90 = quantile(B,probs=0.95,na.rm=T),           B.U50 = quantile(B,probs=0.75,na.rm=T),
                                                                                   B.L90 = quantile(B,probs=0.05,na.rm=T),           B.L50 = quantile(B,probs=0.25,na.rm=T),
                                                  Rec.mn = median(Rec,na.rm=T),      Rec.U90 = quantile(Rec,probs=0.95,na.rm=T),       Rec.U50 = quantile(Rec,probs=0.75,na.rm=T),
                                                                                   Rec.L90 = quantile(Rec,probs=0.05,na.rm=T),       Rec.L50 = quantile(Rec,probs=0.25,na.rm=T),
                                                  Catch.mn = mean(Catch,na.rm=T),  Catch.U90 = quantile(Catch,probs=0.95,na.rm=T),   Catch.U50 = quantile(Catch,probs=0.75,na.rm=T),
                                                                                   Catch.L90 = quantile(Catch,probs=0.05,na.rm=T),   Catch.L50 = quantile(Catch,probs=0.25,na.rm=T),
                                                  M.mn = median(M,na.rm=T),          M.U90 = quantile(M,probs=0.95,na.rm=T),           M.U50 = quantile(M,probs=0.75,na.rm=T),
                                                                                   M.L90 = quantile(M,probs=0.05,na.rm=T),           M.L50 = quantile(M,probs=0.25,na.rm=T),
                                                  F.mn = median(F.dec.table,na.rm=T),F.U90 = quantile(F.dec.table,probs=0.95,na.rm=T), M.U50 = quantile(F.dec.table,probs=0.75,na.rm=T),
                                                                                   F.L90 = quantile(F.dec.table,probs=0.05,na.rm=T), M.L50 = quantile(F.dec.table,probs=0.25,na.rm=T),
                                                  SP.mn = median(SP,na.rm=T),        SP.U90 = quantile(SP,probs=0.95,na.rm=T),         SP.U50 = quantile(SP,probs=0.75,na.rm=T),
                                                                                   SP.L90 = quantile(SP,probs=0.05,na.rm=T),         SP.L50 = quantile(SP,probs=0.25,na.rm=T))
    
    
    # Biomass realizations          
    B.real <- ggplot(Exp.res,aes(x=year,y=B/1000,group = Sim)) + geom_line(alpha = alphs) + #facet_wrap(~F.scenario) +
                                                                      ylab("Predicted Biomass (metric tonnes x 1000)") + xlab("") + 
                                                                      geom_hline(yintercept = c(HCR.sim$LRP/1000,HCR.sim$USR/1000,HCR.sim$TRP/1000),
                                                                                 color=c('firebrick2',u.colors[2],u.colors[1]),linewidth = 1.5) +
                                                                      #scale_color_viridis_b(alpha=0.4,end = 0.75) +
                                                                       theme(legend.position = 'none')
    
    ggsave(paste0(plot_sims,"B_reals.png"),plot=B.real,height=15,width=15)               
    
    # Catch realizations          
    C.real <- ggplot(Exp.res,aes(x=year,y=Catch,group = Sim)) + geom_line(alpha = alphs) + #facet_wrap(~F.scenario) +
                                                                          ylab("Predicted Catch (metric tonnes)") + xlab("") + 
                                                                          #geom_hline(yintercept = c(HCR.sim$LRP/1000,HCR.sim$USR/1000,HCR.sim$TRP/1000),color=c('firebrick2','green','blue')) +
                                                                          #scale_color_viridis_b(alpha=0.4,end = 0.75) +
                                                                           theme(legend.position = 'none')
    
    ggsave(paste0(plot_sims,"C_reals.png"),plot=C.real,height=15,width=15)  
    
    
    # Biomass time series
    B.ts <- ggplot(HCR.summarized,aes(x=year,y=B.mn/1000)) + geom_line() + 
                                                        geom_ribbon(aes(x=year,ymax = B.U90/1000,ymin=B.L90/1000),fill='firebrick2',alpha=0.2) +
                                                        geom_ribbon(aes(x=year,ymax = B.U50/1000,ymin=B.L50/1000),fill='blue',alpha=0.2) + 
                                                        geom_hline(yintercept = c(HCR.sim$LRP/1000,HCR.sim$USR/1000,HCR.sim$TRP/1000),
                                                                   color=c('firebrick2',u.colors[2],u.colors[1]),linewidth = 1.5) +
                                                         xlab("") + ylab("Biomass (metric tonnes x 1000)") 
    ggsave(paste0(plot_sims,"B_mn_ts.png"),plot=B.ts,height=15,width=15)
    
    
    # Catch time series 
    C.ts <- ggplot(HCR.summarized,aes(x=year,y=Catch.mn)) + geom_line() + 
                                                            geom_ribbon(aes(x=year,ymax = Catch.U90,ymin=Catch.L90),fill='firebrick2',alpha=0.2) + 
                                                            geom_ribbon(aes(x=year,ymax = Catch.U50,ymin=Catch.L50),fill='blue',alpha=0.2) + 
                                                             #geom_hline(yintercept = c(2700,7200,9000),color=c('firebrick2','green','blue')) +
                                                             xlab("") + ylab("Catch (metric tonnes)") 
    ggsave(paste0(plot_sims,"C_mn_ts.png"),plot=C.ts,height=15,width=15)
    
  }
} # end if(save.results==T)

  
  
  
  
  
  
  return(list(Catch = Catch,Rec = Rec, rps = rps,B = Bio, mort = mort, mort.r=mort.r,g=g, gR= gR, Surp.prod = Surp.prod, 
              r.realized = r.realized, r.no.fishing=r.no.fishing,exp.scenario = exp.scenario,f.table=f.table,
              Exp.res=Exp.res,m.cor.ts = m.cor.ts,rps.cor.ts = rps.cor.ts))
} # End of the model function
