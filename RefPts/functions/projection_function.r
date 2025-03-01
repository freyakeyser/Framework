# This is going to be an overly complex function in which we can set up various productivity modules to project the population into the future.
# Will want harvest control rules in this too eventually I think

# Building the options, the suite of options are the same for each of the 'mods'. You look at the relationships betwen each productivity parameter and your biomass to determine what is most appropriate
# run:  For the moment if you put 'model_error' it grabs the error estimates from the model for Rec m, and the process error term.  If you put anything else here there is no process error 
#       in the model and the Rec/m uncertainty is coming from the data itself.  I'm still unsure how I want to handle this...
#g.mod: What is our growth model.  
      #a:  type: Basically the options are a linear model looking at RPS~SSB, or a breakpoint analysis that splits the RPS into a low and high period.  If the relationship appears to be
      #          density independent, you can do a 'bp' and set the bp to be 0.  Options are 'lm', 'glm', 'gam', and 'bp'.  Default is 'bp'

      #b: bp.strategy: This is only used if type = 'bp'.  Options are 'sample' where we just sample growth from the historically observed growth data and 'dist' where we select growth from a 
      #                normal distribution with mean/variance observed in data  Default is 'sample'
      #c: B.rec.bp : Is there a biomass breakpoint where growth appears to be lower/higher on either side.  The model will select growth rates appropriate for the current biomass based on the breakpoint.
          # Default sets this to 0 which effectively is no breakpoint.
#rec.mod: What is our recruitment model.  This currently has a few options, pretty much infinitely more options are possible to add in here.
      #a: type: Options right now are a linear model "lm", which will effectively be a Ricker model linearized, so log(RPS) ~ SSB model using the productivity data
          # 'bp' 
#m.mod: What is our natural mortality model.  
      
      

proj.mod <- function(mods = list(tlm.mod = tlm.mod,seam.mod = seam.mod), n_sim = 300,exp.scenario=seq(0,0.4,0.025),n_y = 100, LRP = 1000, HCR.sim = NULL, save.results = F,
                     ci.proj = data.frame(lci=0.25,uci=0.75),run = 'model_error',
                     plot_url = "D:/Github/BBn_model/Results/Figures/BBn/",
                     res_url = "D:/Github/BBn_model/Results/Models/BBn/",
                     model = "TLM", max.SSB =NULL,
                     # He we have our productivity models
                     g.mod =   list(type = 'bp',bp.strategy = 'sample',bp= 0,mn.at.max = 1.00,sd.at.max = 0.05,ar1=0,ar2=0), 
                     rec.mod = list(type = 'bp',bp.strategy = 'sample',bp= 0,mn.at.max = 0.01,sd.at.max = 0.01,rec.age = 5,
                                    ar1=0,ar2=0,rps.m.cor = 0,rps.m.lag=0), # the rps.m.cor and rps.m.lag are only relevant if both rec.mod and m.mod are set as 'cor'),
                     m.mod =   list(type = 'bp',bp.strategy = 'sample',bp= 0,mn.at.max = 0.30,sd.at.max = 0.05,ar1=0,ar2=0)
                     )
                     

{
  
# Now get the plots and results folders sorted out...
if(is.null(HCR.sim)) sims = paste0("Ref_points/",model)
if(!is.null(HCR.sim))  sims = paste0("HCR/",model,"_TRP_",HCR.sim$TRP,"_",HCR.sim$TRP.exp,"_USR_",HCR.sim$USR,"_",HCR.sim$USR.exp,"_LRP_",HCR.sim$LRP,"_",HCR.sim$LRP.exp,"_",HCR.sim$exp.sd)

# Names and locations of save objects
plot_sims <- paste0(plot_url,sims,"_g_",g.mod$type,"_",g.mod$bp,"_",g.mod$bp.strategy,"_",g.mod$mn.at.max,"_",g.mod$sd.at.max,"_r_",
                      rec.mod$type,"_",rec.mod$bp.strategy,"_",rec.mod$bp,"_",rec.mod$mn.at.max,"_",rec.mod$sd.at.max,rec.mod$rec.age,"_m_",
                      m.mod$type,"_",m.mod$bp.strategy,"_",m.mod$bp,"_",m.mod$mn.at.max,"_",m.mod$sd.at.max,"/")


res_sims <- paste0(res_url,sims,"_g_",g.mod$type,"_",g.mod$bp,"_",g.mod$bp.strategy,"_",g.mod$mn.at.max,"_",g.mod$sd.at.max,"_r_",
                     rec.mod$type,"_",rec.mod$bp.strategy,"_",rec.mod$bp,"_",rec.mod$mn.at.max,"_",rec.mod$sd.at.max,rec.mod$rec.age,"_m_",
                     m.mod$type,"_",m.mod$bp.strategy,"_",m.mod$bp,"_",m.mod$mn.at.max,"_",m.mod$sd.at.max,".Rds")



# Will need to do these up right when the repo goes public.  
# source("D:/Github/BBn_model/Scripts/Density_dependence_function.R")
# source("D:/Github/BBn_model/Scripts/breakpoint_function.R")

  library(mgcv)
  
  if(model == "TLM") mod.fit <- mods$tlm.mod
  if(model == "SEAM") mod.fit <- mods$seam.mod  
  
  # If doing the process error run you need to do this, 
  #DK note: I've turned off the process error term in the model, it just results in the model producing biomass out of nothing, which I don't think is good behaviour for projections.
  if(run == 'model_error' & model == "TLM")
  {
    sigma_tau <- 0 #mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'sigma_tau'] #Biomass
    sigma_phi <- mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'sigma_phi'] # Recruits
    sigma_m   <- mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'sigma_m'] # Natural mort
  }
  
  if(run == 'model_error' & model == "SEAM")
  {
    sigma_tau <- 0 #mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'SigmaO_B'] #Biomass
    sigma_phi <- mod.fit$sdrep$value[names(mod.fit$sdrep$value) == 'SigmaO_R'] # Recruits
    sigma_m   <- mod.fit$sdrep$value[names(mod.fit$sdrep$value) == "SigmaO_m"] # Natural mort
  }

  #  If we aren't doing the process error thing then set sigma_tau to 0.  
  # For the recruit and natural mortality no error models we do something different and is handled down in the code below.
  if(run != 'model_error') sigma_tau <- 0
  
  
  
  if(is.null(HCR.sim)) n.f.scen <- length(exp.scenario)
  if(!is.null(HCR.sim)) n.f.scen <- 1
  # Using the recruitment inline below based on 
  #Rec<-median(mod.fit$report$totR) # Median recruitment
  #Rec<-min(mod.fit$report$totR) # Min recruitment
  # Setting up the output arrays
  Bio<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  Catch<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  mort<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  g<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  gR<-array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  Rec <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  Surp.prod <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  r.no.fishing <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  r.realized <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  f.table <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  rps <- array(rep(NA,n_y*(n_sim)*n.f.scen),c(n_y,(n_sim),n.f.scen))
  
  Exp.res <- NULL
  #browser()
  # Finally I make this dataframe to make our lives easier where I combine the main productivity parameters we need.
  # We want to scrub off the most recent year as that is a projection year.
  prod.dat <- data.frame(B = mod.fit$report$totB[-length(mod.fit$report$totB)],
                          R = mod.fit$report$totR[-length(mod.fit$report$totR)], # For some reason the SEBDAM run needs a NA at the end here.
                          m = mod.fit$report$mean_m[-length(mod.fit$report$mean_m)],
                          g = as.vector(mod.fit$obj$env$data$g[-length(mod.fit$obj$env$data$g)]),
                          gR = as.vector(mod.fit$obj$env$data$gR[-length(mod.fit$obj$env$data$g)]))
  #browser()
  # A simple way to get gR estimate from the g estimate 
  prod.dat$g.ratio <- prod.dat$gR/prod.dat$g
   
  # For the growth data we are could divide it into high and low growth dataframes based on g.mod$bp.g.B setting
  # first step, divide the growth data if we have a bp.g.B which means we have simple density dependent growth
  # where there is evidence that growth is a function of biomass (but there isn't a tidy linear relationship)

  # For the recruitment analysis we can get some preliminaries out of the way.
  # Using the Biomass of Recruits + Fully Recruited Scallop as our estimate for SSB.
  prod.dat$SSB <- prod.dat$B + prod.dat$R
  # Getting a vector of recruits offset so we can use that for our RPS analysis
  # First calculate Recruits per spawner based on the Recruitment Age chosen
  prod.dat$Recs <- c(prod.dat$R[(rec.mod$rec.age+1):nrow(prod.dat)],rep(NA,rec.mod$rec.age))
  # calculate recruits per spawner
  prod.dat$RPS <-  prod.dat$Recs/(prod.dat$SSB)
  min.SSB  <- ceiling(min(prod.dat$SSB/100,na.rm=T)) * 100 # This will round up to the 100's of tonnes, below this the RPS will be based on the observed RPS values at low biomass
  #max.SSB.4.pred  <-  floor(max(prod.dat$SSB/100,na.rm=T)) * 100 # This will be the top end of the RPS, above this value we simply have recruitment at the lowest observed recruitment in history (with uncertainty)
  # We are going to use max.SSB to reduce productivity when the population is at a level above anything observed historically, or should we use max SSB above... hmm.... B for now
  if(is.null(max.SSB)) max.SSB <- floor(max(prod.dat$SSB/100,na.rm=T)) * 100
  
  

    
  # If we are using some form a linear density dependence we run these here to the prediction object we will need to use to get a prediction of the response variable
  if(rec.mod$type %in% c('lm','glm','gam')) r.res <- dens.function(dat = prod.dat,type = rec.mod$type,min.cov.pred = min.SSB,max.cov.pred= max.SSB,step = 100,response = 'RPS',dd.term = "SSB",log.trans=T)
  if(g.mod$type %in% c('lm','glm','gam'))   g.res <- dens.function(dat = prod.dat,type = g.mod$type,min.cov.pred = min.SSB,max.cov.pred= max.SSB,step = 100,response = 'g',dd.term = "SSB",log.trans=T)
  #if(g.mod$type %in% c('lm','glm','gam'))   gr.res <- dens.function(dat = prod.dat,type = g.mod$type,min.cov.pred = min.SSB,max.cov.pred= max.SSB,step = 100,response = 'gR',dd.term = "SSB",log.trans=T)
  if(m.mod$type %in% c('lm','glm','gam'))   m.res <- dens.function(dat = prod.dat,type = m.mod$type,min.cov.pred = min.SSB,max.cov.pred= max.SSB,step = 100,response = 'm',dd.term = "SSB",log.trans=T)

  

  for (i in 1:n.f.scen) 
  {
    #browser()
    Sim.res <- NULL # reset the Sim res object, it'll just be a temporary container of the sims
    for (n in 1:n_sim)
    {  
      # If we are running the correlation analysis for any of these, we do that outside the j loop as we can 'make' a whole time series in one shot, of course the
      # values in that time series will need to be 'replaced' if the biomass gets too high, but we'll do that below just like how we do it with the lm models.
      # Finally we can run a stand alone correlation analysis for g.mod.
      if(g.mod$type == 'cor')
      {
       #browser()
        # Because this is done on log scale (needed for RPS and nm really), I'm going to force the variance to be a bit artificially low here, this seems to give very reasonable growths
        if(g.mod$ar1 == 0 & g.mod$ar2 == 0)
        {
          g.ts <-  cor.fun(ts1 = list(mn = median(prod.dat$g,na.rm=T),var=sd(log(prod.dat$g),na.rm=T)/sqrt(length(prod.dat$g))),
                           ts2 = list(mn =0.1,var=0.1),
                           arima = list(proc = 'normal',cross.cor = 0,lag = 0),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts  
        } else {
          g.ts <-  cor.fun(ts1 = list(mn = median(prod.dat$g,na.rm=T),var=sd(log(prod.dat$g),na.rm=T)/sqrt(length(prod.dat$g))),
                          ts2 = list(mn =0.1,var=0.1),
                          arima = list(proc = 'arima',cross.cor = 0,lag = 0,ts1.ar1= g.mod$ar1,ts1.ar2=g.mod$ar2),
                          ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
        } # end else
      } #end g.mod$type == 'cor'
      
      # Now do the same thing for the recruits, but I also have a version where we allow rps and m to covary together, so this is the version where we only have recs doing their thing
      if(rec.mod$type == 'cor' & rec.mod$rps.m.cor ==0)
      {
        #browser()
        # Because this is done on log scale (needed for RPS and nm really), I'm going to force the variance to be a bit artificially low here, this seems to give very reasonable growths
        if(rec.mod$ar1 == 0 & rec.mod$ar2 == 0)
        {
          rps.ts <-  cor.fun(ts1 = list(mn = median(prod.dat$RPS,na.rm=T),var=sd(log(prod.dat$RPS),na.rm=T)/sqrt(length(prod.dat$RPS))),
                           ts2 = list(mn =0.1,var=0.1),
                           arima = list(proc = 'normal',cross.cor = 0,lag = 0),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts  
        } else {
          rps.ts <-  cor.fun(ts1 =list(mn = median(prod.dat$RPS,na.rm=T),var=sd(log(prod.dat$RPS),na.rm=T)/sqrt(length(prod.dat$RPS))),
                           ts2 = list(mn =0.1,var=0.1),
                           arima = list(proc = 'arima',cross.cor = 0,lag = 0,ts1.ar1= rec.mod$ar1,ts1.ar2=rec.mod$ar2),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
        } # end else
      } #end rec.mod$type == 'cor'
      
      # Now do the same thing for the natural mortality, but I also have a version where we allow rps and m to covary together, so this is the version where we only have recs doing their thing
      if( m.mod$type == 'cor' & rec.mod$rps.m.cor ==0)
      {
        #browser()
        # Because this is done on log scale (needed for RPS and nm really), I'm going to force the variance to be a bit artificially low here, this seems to give very reasonable growths
        if(m.mod$ar1 == 0 & m.mod$ar2 == 0)
        {
          m.ts <-  cor.fun(ts1 = list(mn = median(prod.dat$m,na.rm=T),var=sd(log(prod.dat$m),na.rm=T)/sqrt(length(prod.dat$m))),
                              ts2 = list(mn = median(prod.dat$RPS,na.rm=T),var=var(log(prod.dat$RPS),na.rm=T)),
                               arima = list(proc = 'normal',cross.cor = 0,lag = 0),
                               ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
          
 
        
        } else {
          m.ts <-  cor.fun(ts1 =list(mn = median(prod.dat$m,na.rm=T),var=sd(log(prod.dat$m),na.rm=T)/sqrt(length(prod.dat$m))),
                             ts2 = list(mn =0.1,var=0.1),
                             arima = list(proc = 'arima',cross.cor = 0,lag = 0,ts1.ar1= m.mod$ar1,ts1.ar2=m.mod$ar2),
                             ts = list(years=500,start.year =300, final.ts.len = n_y+10))$ts1.ts
        } # end else
      } #end m.mod$type == 'cor'
      
      
      if( m.mod$type == 'cor' & rec.mod$rps.m.cor !=0)
      {
        #browser()
        # Because this is done on log scale (needed for RPS and nm really), I'm going to force the variance to be a bit artificially low here, this seems to give very reasonable growths
        if(m.mod$ar1 == 0 & m.mod$ar2 == 0)
        {
          tmp.res <-  cor.fun(ts1 = list(mn = median(prod.dat$m,na.rm=T),var=sd(log(prod.dat$m),na.rm=T)/sqrt(length(prod.dat$m))),
                           ts2 =  list(mn = median(prod.dat$RPS,na.rm=T),var=sd(log(prod.dat$RPS),na.rm=T)/sqrt(length(prod.dat$RPS))),
                           arima = list(proc = 'normal',cross.cor = rec.mod$rps.m.cor,lag = rec.mod$rps.m.lag),
                           ts =list(years=500,start.year =300, final.ts.len = n_y+10))
          m.ts <- tmp.res$ts1.ts
          rps.ts <- tmp.res$ts2.ts
        } else {
          tmp.res <-  cor.fun(ts1 =list(mn = median(prod.dat$m,na.rm=T),var=sd(log(prod.dat$m),na.rm=T)/sqrt(length(prod.dat$m))),
                           ts2 =  list(mn = median(prod.dat$RPS,na.rm=T),var=sd(log(prod.dat$RPS),na.rm=T)/sqrt(length(prod.dat$RPS))),
                           arima = list(proc = 'arima',cross.cor = rec.mod$rps.m.cor,lag = rec.mod$rps.m.lag,ts1.ar1= m.mod$ar1,ts1.ar2=m.mod$ar2,ts2.ar1= rec.mod$ar1,ts2.ar2=rec.mod$ar2),
                           ts = list(years=500,start.year =300, final.ts.len = n_y+10))
          m.ts <- tmp.res$ts1.ts
          rps.ts <- tmp.res$ts2.ts
        } # end else
      } #end m.mod$type == 'cor'
      
      
      for (j in 1:n_y)
      {
        # Get a tidy biomass and Recruit value to save on code below...
        # For the first year we take the biomass from last year as our starting point
        if (j == 1) 
        {
          B.init <- mod.fit$report$totB[length(mod.fit$report$totB)]
          Rec.init <- mod.fit$report$totR[length(mod.fit$report$totR)-1]
        } # end  if (j == 1) 
        # After the first year we have last years results.
        if (j > 1) 
        {
          B.init <- Bio[j,n,i]
          Rec.init <- Rec[j-1,n,i]
        } # end if (j > 1)
        
        # Get our quasi-SSB from last year, this is what we'll have to use to parameterize our density terms.
        # DK Note: We should be thinking about the one year offset this causes, this is saying last year SSB influences the upcoming DD process.
        # I don't think there's much else one can do, but that might be subtly different than what we are doing with the data exploration which uses growth/rec, etc from the 'same' year.
        SSB.init <- B.init + Rec.init
        ##################################################
        ### Growth Section ###
        
        # First up commercial sized growth
        # using the breakpoint method determine what growth is.  The BP method wraps 
        # all scenarios into one tidy bit of code.  Also, note that by default this includes all data because the default for bp.g.B = 0

        if(g.mod$type == 'bp') g[j,n,i] <- bp.function(B = SSB.init,dat = prod.dat$g,type = g.mod$bp.strategy, bp = g.mod$bp,
                                                       max.B = max.SSB,mn.at.max = g.mod$mn.at.max,sd.at.max = g.mod$sd.at.max)
        # Now if we are using a linear model to predict growth then we have to do this...
        if(g.mod$type %in% c("lm",'glm','gam'))
        {
          # If below the minimum SSB ever observed assume growth is what was observed at the lowest SSB in the time series. Do not extrapolate...
          # This isn't terrible for a log normal, the huge outliers were so outlier that this doesn't really go far enough on BBn
          if(SSB.init < min.SSB) g[j,n,i] <- rlnorm(1,g.res$pred.dat$g.log[1],g.res$pred.dat$se[1])
          # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
          # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
          if(SSB.init > max.SSB) g[j,n,i] <- rlnorm(1,log(g.mod$mn.at.max),g.mod$sd.at.max) 
          # Finally the most complicated one is to pick from the predictions and get a recruitment estimate based on the linear model predictions and uncertainty
          # We do this when the SSB is within observed bounds.
          if(SSB.init >= min.SSB & SSB.init <= max.SSB)
          {
            pick <- which(g.res$pred.dat$SSB == round(SSB.init/100)* 100)
            g[j,n,i] <- rlnorm(1,g.res$pred.dat$g.log[pick],g.res$pred.dat$se[pick]) 
          } # end if(SSB.for.rec >= min.SSB & SSB.for.rec <= max.SSB)
        } # end if(rec.mod$type %in% c("lm",'glm','gam'))
        if(g.mod$type == 'cor')
        {
          
          # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
          # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
          if(SSB.init > max.SSB) g[j,n,i] <- rlnorm(1,log(g.mod$mn.at.max),g.mod$sd.at.max) 
          # If above the breakpoint sample from the data above the breakpoint...
          if(SSB.init > g.mod$bp & SSB.init <= max.SSB)  g[j,n,i] <- rnorm(1,median(prod.dat$g[prod.dat$SSB > g.mod$bp],na.rm=T),sd(prod.dat$g[prod.dat$SSB > g.mod$bp],na.rm=T))
          # For the rest grab the data from the correlation
          if(SSB.init <= g.mod$bp)  g[j,n,i] <- g.ts[j]
        } #end if g.mod$type == 'cor'
             # And since we know gR is related to g, let's just make gR be X% larger than g, we can just sample from the observed data to get that ratio...
        # We can see that growth differences really doesn't vary too much and that recruit growth is always larger (basically between 10 and 25% larger)
         gR[j,n,i] <- sample(prod.dat$g.ratio,1) * g[j,n,i]
         #browser()
        
        ################################# End Growth Section ##################################
        
        
        
       ############################### Start recruitment Section #######################################
        ### Recruit section
        # For the recruits I need to get the right biomass year to get the offset from recruits to the biomass that produced those recruits.  Easiest to do up here 
        if(j < (rec.mod$rec.age+1))  SSB.for.rec <- mod.fit$report$totB[length(mod.fit$report$totB)-rec.mod$rec.age] + mod.fit$report$totB[length(mod.fit$report$totB)-rec.mod$rec.age]
        if(j >=(rec.mod$rec.age+1)) SSB.for.rec <- Bio[j-rec.mod$rec.age,n,i] + Rec[j-rec.mod$rec.age,n,i]
        # Now we can get our recruit estimate...
        # Using the predictions from our linear model above...
        if(rec.mod$type %in% c("lm",'glm','gam'))
        {
          # If below the minimum SSB ever observed assume RPS is what was predicted at the lowest SSB in the time series. Do not extrapolate...
          # This isn't terrible for a log normal, the huge outliers were so outlier that this doesn't really go far enough on BBn
          if(SSB.for.rec < min.SSB) Rec[j,n,i] <- SSB.for.rec*  rlnorm(1,r.res$pred.dat$RPS.log[1],r.res$pred.dat$se[1]) 
          # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
          # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
          if(SSB.for.rec > max.SSB) Rec[j,n,i] <- SSB.for.rec* rlnorm(1,log(rec.mod$mn.at.max),rec.mod$sd.at.max) 
          # Finally the most complicated one is to pick from the predictions and get a recruitment estimate based on the linear model predictions and uncertainty
          # We do this when the SSB is within observed bounds.
          if(SSB.for.rec >= min.SSB & SSB.for.rec <= max.SSB)
          {
            pick <- which(r.res$pred.dat$SSB == round(SSB.for.rec/100)* 100)
            Rec[j,n,i] <- rlnorm(1,r.res$pred.dat$RPS.log[pick],r.res$pred.dat$se[pick]) * SSB.for.rec
          } # end if(SSB.for.rec >= min.SSB & SSB.for.rec <= max.SSB)
        } # end if(rec.mod$type %in% c("lm",'glm','gam'))

        if(rec.mod$type %in% 'bp')
        {
          # DK been playing around with the best way too do this.
          if(!run %in% 'model_error') 
          {
            #browser()
            # Should I use the RPS or Recruitment time series for these??
             #if(SSB.for.rec >=min.SSB)  
            Rec[j,n,i] <-SSB.for.rec * bp.function(B = SSB.for.rec,dat = prod.dat$RPS,type = rec.mod$bp.strategy, bp = rec.mod$bp,
                                                                            max.B = max.SSB,mn.at.max = rec.mod$mn.at.max,sd.at.max = rec.mod$sd.at.max)
            # # If below the minimum obsreved biomass than recruitment will be based on the RPS at low biomass levels.
            #if(SSB.for.rec < min.SSB)  Rec[j,n,i] <- SSB.for.rec * bp.function(B = SSB.for.rec,dat = prod.dat$RPS,type = rec.mod$bp.strategy, bp = rec.mod$bp,
            #                                                             max.B = max.SSB,mn.at.max = rec.mod$mn.at.max,sd.at.max = rec.mod$sd.at.max)
          } # end if(run != 'model_error')     
          
          if(run == 'model_error') 
          {
            
            # Should I use the RPS or Recruitment time series for these??
            #if(SSB.for.rec >= min.SSB) 
              Rec[j,n,i] <- SSB.for.rec * bp.function(B = SSB.for.rec,dat = prod.dat$RPS,type = rec.mod$bp.strategy, bp = rec.mod$bp,sd = sigma_phi,
                                                                            max.B = max.SSB,mn.at.max = rec.mod$mn.at.max,sd.at.max = rec.mod$sd.at.max)
            # If below the minimum obsreved biomass than recruitment will be based on the RPS at low biomass levels.
            #if(SSB.for.rec < min.SSB)  Rec[j,n,i] <- SSB.for.rec * bp.function(B = SSB.for.rec,dat = prod.dat$RPS,type = rec.mod$bp.strategy, bp = rec.mod$bp,sd = sigma_phi,
            #                                                                          max.B = max.SSB,mn.at.max = rec.mod$mn.at.max,sd.at.max = rec.mod$sd.at.max)
          } # end if(run == 'model_error') 
           
         
        } # end if(rec.mod$type %in% 'bp')
         
         if(rec.mod$type == 'cor')
         {
           # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
           # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
           if(SSB.for.rec > max.SSB) Rec[j,n,i] <- SSB.for.rec* rlnorm(1,log(rec.mod$mn.at.max),rec.mod$sd.at.max) 
           # If above the breakpoint sample from the data above the breakpoint...
           if(SSB.for.rec > rec.mod$bp & SSB.for.rec <= max.SSB)  Rec[j,n,i] <-SSB.for.rec* rnorm(1,median(prod.dat$RPS[prod.dat$SSB > rec.mod$bp],na.rm=T),sd(prod.dat$RPS[prod.dat$SSB > rec.mod$bp],na.rm=T))
           # For the rest grab the data from the correlation
           if(SSB.for.rec <= rec.mod$bp)  Rec[j,n,i] <- SSB.for.rec*rps.ts[j]
         }
         # I'd like too be able to look at the rps time series too...
         rps[j,n,i] <- Rec[j,n,i]/SSB.for.rec
        
        ############################### End recruitment Section #######################################

        
        ############################### Start natural mortality Section #######################################
        # kinda a hybrid of the above two modules
        # Now if we are using a linear model to predict growth then we have to do this...
        if(m.mod$type %in% c("lm",'glm','gam'))
        {
          # If below the minimum SSB ever observed assume growth is what was observed at the lowest SSB in the time series. Do not extrapolate...
          # This isn't terrible for a log normal, the huge outliers were so outlier that this doesn't really go far enough on BBn
          if(SSB.init < min.SSB) mort[j,n,i] <- rlnorm(1,m.res$pred.dat$m.log[1],m.res$pred.dat$se[1])
          # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
          # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
          if(SSB.init > max.SSB) mort[j,n,i] <- rlnorm(1,log(m.mod$mn.at.max),m.mod$sd.at.max) 
          # Finally the most complicated one is to pick from the predictions and get a recruitment estimate based on the linear model predictions and uncertainty
          # We do this when the SSB is within observed bounds.
          if(SSB.init >= min.SSB & SSB.init <= max.SSB)
          {
            pick <- which(m.res$pred.dat$SSB == round(SSB.init/100)* 100)
            mort[j,n,i] <- rlnorm(1,m.res$pred.dat$m.log[pick],m.res$pred.dat$se[pick]) 
          } # end if(SSB.for.rec >= min.SSB & SSB.for.rec <= max.SSB)
        } # end if(rec.mod$type %in% c("lm",'glm','gam'))
        
        if(m.mod$type %in% 'bp')
        {
          # 
          if(run != 'model_error') mort[j,n,i] <- bp.function(B = SSB.for.rec,dat = prod.dat$m, type = m.mod$bp.strategy, bp = m.mod$bp,
                                                                              max.B = max.SSB,mn.at.max = m.mod$mn.at.max,sd.at.max = m.mod$sd.at.max)
          
          if(run == 'model_error') mort[j,n,i] <- bp.function(B = SSB.for.rec,dat = prod.dat$m,type = m.mod$bp.strategy, bp = m.mod$bp,sd = sigma_m,
                                                                              max.B = max.SSB,mn.at.max = m.mod$mn.at.max,sd.at.max = m.mod$sd.at.max)
        } # end if(rec.mod$type %in% 'bp')    
         
         if(m.mod$type == 'cor')
         {
           # If above the maximum SSB ever observed, set recruitment to vary around the lowest recruitment numbers ever observed
           # This lets the user pick what the distribution will look like above the max.SSB, default allows for about 100 tonnes of recruitment at very high SSB.
           if(SSB.init > max.SSB) mort[j,n,i] <- rlnorm(1,log(m.mod$mn.at.max),m.mod$sd.at.max) 
           # If above the breakpoint sample from the data above the breakpoint...
           if(SSB.init > m.mod$bp & SSB.init <= max.SSB)  mort[j,n,i]  <- rlnorm(1,median(log(prod.dat$m[prod.dat$SSB > m.mod$bp]),na.rm=T),sd(log(prod.dat$m[prod.dat$SSB > m.mod$bp]),na.rm=T))
           # For the rest grab the data from the correlation
           if(SSB.init <= m.mod$bp)  mort[j,n,i] <- m.ts[j]
         }
        
        ############################### End natural mortality Section #######################################
        #browser()
        
        
        ### Now if we are running the correlation on two of the above (that's the most we can do)
        
        
        
        #### Now we run the model and extract some metric of interest while we are at it.
        # Get our Catch based on the exp.scenario scenario and based on the previous years biomass
        # DK Note, when/how we remove the catch is slightly different from what we do in our projections, it's mostly about how we calculated F
        # In our current decision tables F would be calculated for the upcoming season using the Bio estimate, so I'm sticking in an 'F.table' for the 'F' we'd calculated in our decision tables.
        
       
        if(j == 1)
        {
          Bio[j,n,i] <- B.init
          Rec[j,n,i] <- Rec.init
        }
        if(j < n_y)
        {
          # Our projection forward
          if(is.null(HCR.sim)) Catch[j,n,i]<-B.init*exp.scenario[i]
          
          if(!is.null(HCR.sim))
          {
            # LRP catch
            if(B.init < HCR.sim$LRP)  Catch[j,n,i] <- rlnorm(1,log(HCR.sim$LRP.exp),HCR.sim$exp.sd) * B.init
            # Assume linear decline in Explotation to LRP in Cautious zone
            if(B.init > HCR.sim$LRP & B.init <= HCR.sim$USR)  Catch[j,n,i] <- rlnorm(1,log(HCR.sim$USR.exp),HCR.sim$exp.sd)*B.init*((B.init-HCR.sim$LRP)/HCR.sim$LRP)
            # Harvest at USR explotation here
            if(B.init > HCR.sim$USR & B.init < HCR.sim$TRP)  Catch[j,n,i] <- rlnorm(1,log(HCR.sim$USR.exp),HCR.sim$exp.sd) * B.init
            # Harvest at TRP explotation here
            if(B.init >= HCR.sim$TRP)  Catch[j,n,i] <- rlnorm(1,log(HCR.sim$TRP.exp),HCR.sim$exp.sd)*B.init
          }
          #browser()
          Bio[j+1,n,i] <- rlnorm(1,log(exp(-mort[j,n,i])*g[j,n,i]*(B.init-Catch[j,n,i])+Rec[j,n,i]*exp(-mort[j,n,i])*gR[j,n,i]),sigma_tau)
          # The Surplus production this won't match Bio
          Surp.prod[j,n,i]<- Bio[j+1,n,i] - B.init + Catch[j,n,i]
          # The realized rate of growth (r_fish)
          r.realized[j,n,i] <- log(Bio[j+1,n,i]/B.init)
          # The potential rate of growth (r), this method allows the catch to grow (or shrink) based on m and g.
          r.no.fishing[j,n,i]<-   log((Bio[j+1,n,i] + exp(-mort[j,n,i])*g[j,n,i]*Catch[j,n,i])/B.init) #log(((exp(-mort[j,n,i])*g[j,n,i]*B.init)+Rec[j,n,i]*exp(-mort[j,n,i])*gR[j,n,i]) /B.init)
          # Also calculate F using the end biomass result which is what we'd currently report as F in our decision tables
          # DK Note: Add in different ways of calculating F to see how different the F estimates might be.
          f.table[j+1,n,i] <-  Catch[j,n,i] / (Bio[j+1,n,i] +  Catch[j,n,i])
        } # end if(j < n_y)
         #browser()
      } # end the j loop through the years
      Sim.res[[n]] <- data.frame(year = 1:n_y,Sim = rep(n,n_y), B = Bio[,n,i],Rec = Rec[,n,i],g = g[,n,i],gR = gR[,n,i],M = mort[,n,i],
                                 Catch = Catch[,n,i], SP = Surp.prod[,n,i],F.dec.table = f.table[,n,i],
                                 r.real = r.realized[,n,i],r.no.f = r.no.fishing[,n,i])
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
  if(is.null(HCR.sim))
  {
    MSY.summarized <- Exp.res %>% dplyr::group_by(F.scenario,year) %>% 
                                   dplyr::summarise(B.mn = median(B,na.rm=T),          B.U95 = quantile(B,probs=0.975,na.rm=T),           B.U50 = quantile(B,probs=0.75,na.rm=T),
                                                                                     B.L95 = quantile(B,probs=0.025,na.rm=T),           B.L50 = quantile(B,probs=0.25,na.rm=T),
                                                    Rec.mn = median(Rec,na.rm=T),      Rec.U95 = quantile(Rec,probs=0.975,na.rm=T),       Rec.U50 = quantile(Rec,probs=0.75,na.rm=T),
                                                                                     Rec.L95 = quantile(Rec,probs=0.025,na.rm=T),       Rec.L50 = quantile(Rec,probs=0.25,na.rm=T),
                                                    Catch.mn = median(Catch,na.rm=T),  Catch.U95 = quantile(Catch,probs=0.975,na.rm=T),   Catch.U50 = quantile(Catch,probs=0.75,na.rm=T),
                                                                                     Catch.L95 = quantile(Catch,probs=0.025,na.rm=T),   Catch.L50 = quantile(Catch,probs=0.25,na.rm=T),
                                                    M.mn = median(M,na.rm=T),          M.U95 = quantile(M,probs=0.975,na.rm=T),           M.U50 = quantile(M,probs=0.75,na.rm=T),
                                                                                     M.L95 = quantile(M,probs=0.025,na.rm=T),           M.L50 = quantile(M,probs=0.25,na.rm=T),
                                                    F.mn = median(F.dec.table,na.rm=T),F.U95 = quantile(F.dec.table,probs=0.975,na.rm=T), M.U50 = quantile(F.dec.table,probs=0.75,na.rm=T),
                                                                                     F.L95 = quantile(F.dec.table,probs=0.025,na.rm=T), M.L50 = quantile(F.dec.table,probs=0.25,na.rm=T),
                                                    SP.mn = median(SP,na.rm=T),        SP.U95 = quantile(SP,probs=0.975,na.rm=T),         SP.U50 = quantile(SP,probs=0.75,na.rm=T),
                                                                                     SP.L95 = quantile(SP,probs=0.025,na.rm=T),         SP.L50 = quantile(SP,probs=0.25,na.rm=T))
    
    Equilib.dat <-    MSY.summarized %>% dplyr::filter(year >= 75) %>% dplyr::group_by(F.scenario) %>% 
                                         dplyr::summarise(B.mn = median(B.mn,na.rm=T),         B.U95 = median(B.U95,na.rm=T),        B.U50 = median(B.U50,na.rm=T),
                                                                                               B.L95 = median(B.L95,na.rm=T),        B.L50 = median(B.L50,na.rm=T),
                                                          Catch.mn = median(Catch.mn,na.rm=T), Catch.U95 = median(Catch.U95,na.rm=T),Catch.U50 = median(Catch.U50,na.rm=T),
                                                                                               Catch.L95 = median(Catch.L95,na.rm=T),Catch.L50 = median(Catch.L50,na.rm=T),
                                                          Rec.mn = median(Rec.mn,na.rm=T),     Rec.U95 = median(Rec.U95,na.rm=T),    Rec.U50 = median(Rec.U50,na.rm=T),
                                                                                               Rec.L95 = median(Rec.L95,na.rm=T),    Rec.L50 = median(Rec.L50,na.rm=T))
    
    
    
    # Get Biomass per Recruit and Yield per Recruit
    Equilib.dat$BPR <- Equilib.dat$B.mn/Equilib.dat$Rec.mn
    Equilib.dat$YPR <- Equilib.dat$Catch.mn/Equilib.dat$Rec.mn
    
    
  # Biomass
    bio_hists<-ggplot(data=Exp.res,aes(x=B))+
                      geom_histogram(fill="gray",col="black")+
                      facet_wrap(~F.scenario,scales="free")+
                      ylab("Count")+xlab("Predicted Biomass (metric tonnes)") +theme_bw()
    
    ggsave(paste0(plot_sims,"Bio_hists.png"),plot=bio_hists,height=10,width=10)
    # Recruits
    rec_hists<-ggplot(data=Exp.res,aes(x=Rec))+
                      geom_histogram(fill="gray",col="black")+
                      facet_wrap(~F.scenario,scales="free") +
                      ylab("Count")+xlab("Predicted Recruitment (metric tonnes)")+ theme_bw()
    
    ggsave(paste0(plot_sims,"Rec_hists.png"),plot=rec_hists,height=10,width=10)
    
    # Catch
    catch_hists<-ggplot(data=Exp.res,aes(x=Catch))+
                        geom_histogram(fill="gray",col="black")+
                        facet_wrap(~F.scenario,scales="free")+
                        ylab("Count")+xlab("Commercial Landings (metric tonnes)") + theme_bw()
    
    ggsave(paste0(plot_sims,"Catch_hists.png"),plot=catch_hists,height=10,width=10)
    #browser()
    # Natural Mortality
    m_hists<-ggplot(data=Exp.res,aes(x=M))+
                    geom_histogram(fill="gray",col="black")+
                    facet_wrap(~F.scenario,scales="free")+
                    ylab("Count")+xlab("Predicted Natural Mortality") + theme_bw()
    
    ggsave(paste0(plot_sims,"M_hists.png"),plot=m_hists,height=10,width=10)
    
    # Fishing Mortality 
    F_hists<-ggplot(data=Exp.res,aes(x=F.dec.table))+
                    geom_histogram(fill="gray",col="black")+
                    facet_wrap(~F.scenario,scales="free")+
                    ylab("Count")+xlab("F from Decision Table") + theme_bw()
    
    ggsave(paste0(plot_sims,"F_hists.png"),plot=F_hists,height=10,width=10)
    
    # Surplus Production
    SP_hists<-ggplot(data=Exp.res,aes(x=SP))+
                    geom_histogram(fill="gray",col="black")+
                    facet_wrap(~F.scenario,scales="free")+
                    ylab("Count")+xlab("Surplus Production (metric tonnes)") + theme_bw()
    
    ggsave(paste0(plot_sims,"SP_hists.png"),plot=SP_hists,height=10,width=10)
  
  
  # Now make plots of the realizations over time
  # First the Biomass
  B.real <- ggplot(Exp.res,aes(x=year,y=B,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Predicted Biomass (metric tonnes)") + xlab("") + 
                  scale_color_viridis_b(alpha=0.4,end = 0.75) +
                  theme_bw() + theme(legend.position = 'none')
  
  ggsave(paste0(plot_sims,"Biomass_realizations.png"),plot=B.real,height=10,width=10)
  # Then the Recruits
  
  R.real <- ggplot(Exp.res,aes(x=year,y=Rec,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Predicted Recruitment (metric tonnes)") + xlab("") + 
                  scale_color_viridis_b(alpha=0.4,end = 0.75) +
                  theme_bw() + theme(legend.position = 'none')
  ggsave(paste0(plot_sims,"Recruit_realizations.png"),plot=R.real,height=10,width=10)
  # Next the Catch
  C.real <- ggplot(Exp.res,aes(x=year,y=Catch,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Predicted Catch (metric tonnes)") + xlab("") + 
                  scale_color_viridis_b(alpha=0.4,end = 0.75) +
                  theme_bw() + theme(legend.position = 'none')
  
  ggsave(paste0(plot_sims,"Catch_realizations.png"),plot=C.real,height=10,width=10)
  #Natural Mortality
  M.real <- ggplot(Exp.res,aes(x=year,y=M,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Natural Mortality (Instantaneous)") + xlab("") + 
                  scale_color_viridis_b(alpha=0.4,end = 0.75) +
                  theme_bw() + theme(legend.position = 'none')
    
  ggsave(paste0(plot_sims,"M_realizations.png"),plot=M.real,height=10,width=10)
  
  #Fishing Mortality
  F.real <- ggplot(Exp.res,aes(x=year,y=F.dec.table,group = Sim,color=Sim)) + 
                  geom_line() + facet_wrap(~F.scenario) +
                  ylab("Fishing Mortality (Proportional)") + xlab("") + 
                  scale_color_viridis_b(alpha=0.4,end = 0.75) +
                  theme_bw() + theme(legend.position = 'none')
  
  ggsave(paste0(plot_sims,"F_realizations.png"),plot=F.real,height=10,width=10)
  
  # Surplus Production
  SP.real <- ggplot(Exp.res,aes(x=year,y=SP,group = Sim,color=Sim)) + 
                    geom_line() + facet_wrap(~F.scenario) +
                    ylab("Surplus Production (tonnes)") + xlab("") + 
                    scale_color_viridis_b(alpha=0.4,end = 0.75) +
                    theme_bw() + theme(legend.position = 'none')
  ggsave(paste0(plot_sims,"SP_realizations.png"),plot=SP.real,height=10,width=10)
  
  
  
  
  ################# estimates with 95% CI bands
  #### Biomass
  B.ts <- ggplot(MSY.summarized,aes(x=year,y=B.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                      geom_ribbon(aes(x=year,ymax = B.U50,ymin=B.L50),fill='firebrick2',alpha=0.2) + 
                                                      theme_bw() + xlab("") + ylab("Biomass (metric tonnes)")
  ggsave(paste0(plot_sims,"B_mn_ts.png"),plot=B.ts,height=10,width=10)
  # Recruits
  Rec.ts <- ggplot(MSY.summarized,aes(x=year,y=Rec.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                          geom_ribbon(aes(x=year,ymax = Rec.U95,ymin=Rec.L95),fill='firebrick2',alpha=0.2) + 
                                                          theme_bw() + xlab("") + ylab("Recruitment (metric tonnes)")
  ggsave(paste0(plot_sims,"R_mn_ts.png"),plot=Rec.ts,height=10,width=10)
  # Catch
  Catch.ts <- ggplot(MSY.summarized,aes(x=year,y=Catch.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                              geom_ribbon(aes(x=year,ymax = Catch.U95,ymin=Catch.L95),fill='firebrick2',alpha=0.2) + 
                                                              theme_bw() + xlab("") + ylab("Catch (metric tonnes)")
  ggsave(paste0(plot_sims,"Catch_mn_ts.png"),plot=Catch.ts,height=10,width=10)
  # Natural mortality
  M.ts <- ggplot(MSY.summarized,aes(x=year,y=M.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                      geom_ribbon(aes(x=year,ymax = M.U95,ymin=M.L95),fill='firebrick2',alpha=0.2) + 
                                                      theme_bw() + xlab("") + ylab("Natural mortality (instantaneous)")
  ggsave(paste0(plot_sims,"M_mn_ts.png"),plot=M.ts,height=10,width=10)
  # Fishing mortality
  F.ts <- ggplot(MSY.summarized,aes(x=year,y=F.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                      geom_ribbon(aes(x=year,ymax = F.U95,ymin=F.L95),fill='firebrick2',alpha=0.2) + 
                                                      theme_bw() + xlab("") + ylab("Fishing mortality (proportional)")
  ggsave(paste0(plot_sims,"F_mn_ts.png"),plot=F.ts,height=10,width=10)
  # Surplus Production
  SP.ts <- ggplot(MSY.summarized,aes(x=year,y=SP.mn)) + geom_line() + facet_wrap(~F.scenario) +
                                                        geom_ribbon(aes(x=year,ymax = SP.U95,ymin=SP.L95),fill='firebrick2',alpha=0.2) + 
                                                        theme_bw() + xlab("") + ylab("Surplus Production (metric tonnes)")
  ggsave(paste0(plot_sims,"SP_mn_ts.png"),plot=SP.ts,height=10,width=10)
  
  
  ####################### Reference Points peice
  # Now for the money shot, we can make the Figures that can help inform removal reference and perhap LRP.
  #Lets take the mean of the last 25 years and use that as our steady state
  
  # Biomass Equilibrium
  B.equil.plt<- ggplot(Equilib.dat,aes(x=F.scenario,y=B.mn)) + geom_point()+ geom_line()+
                                                               geom_ribbon(aes(x=F.scenario,ymin=B.L95, ymax=B.U95),alpha=0.05,fill='firebrick2') +
                                                               theme_bw()+ylab("Equilibrium Biomass (metric tonnes)") + xlab("Explotation Rate (proportional)")
  ggsave(paste0(plot_sims,"B_equilibrum.png"),plot=B.equil.plt,height=6,width=6)
  # Catch Equilibrium
  C.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=Catch.mn)) + geom_point()+ geom_line()+
                                                                  geom_ribbon(aes(x=F.scenario,ymin=Catch.L95, ymax=Catch.U95),alpha=0.05,fill='firebrick2') +
                                                                  theme_bw()+ylab("Equilibrium Catch (metric tonnes)")  + xlab("Explotation Rate (proportional)")
  ggsave(paste0(plot_sims,"Catch_equilibrium.png"),plot=C.equil.plt,height=6,width=6)
  # Biomass per Recruit (Biomass)
  BPR.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=BPR)) + geom_point()+ geom_line()+
                                                               theme_bw()+ylab("Biomass per Recruit (Biomass)")  + xlab("Explotation Rate (proportional)")
  ggsave(paste0(plot_sims,"BPR_equilibrium.png"),plot=BPR.equil.plt,height=6,width=6)
  
  
  # Yield per Recruit (Biomass)
  YPR.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=YPR)) + geom_point()+ geom_line()+
                                                               theme_bw()+ylab("Yield per Recruit (Biomass)")  + xlab("Explotation Rate (proportional)")
  ggsave(paste0(plot_sims,"YPR_equilibrium.png"),plot=YPR.equil.plt,height=6,width=6)
  
  } #end if is.null(HCR.sim)
  
  
  #browser()
  if(!is.null(HCR.sim))
  {
    # Summarize the data to make a nice time series plot.
    HCR.summarized <- Exp.res %>% dplyr::group_by(year) %>% 
                                  dplyr::summarise(B.mn = median(B,na.rm=T),         B.U95 = quantile(B,probs=0.975,na.rm=T),           B.U50 = quantile(B,probs=0.75,na.rm=T),
                                                                                   B.L95 = quantile(B,probs=0.025,na.rm=T),           B.L50 = quantile(B,probs=0.25,na.rm=T),
                                                  Rec.mn = median(Rec,na.rm=T),      Rec.U95 = quantile(Rec,probs=0.975,na.rm=T),       Rec.U50 = quantile(Rec,probs=0.75,na.rm=T),
                                                                                   Rec.L95 = quantile(Rec,probs=0.025,na.rm=T),       Rec.L50 = quantile(Rec,probs=0.25,na.rm=T),
                                                  Catch.mn = median(Catch,na.rm=T),  Catch.U95 = quantile(Catch,probs=0.975,na.rm=T),   Catch.U50 = quantile(Catch,probs=0.75,na.rm=T),
                                                                                   Catch.L95 = quantile(Catch,probs=0.025,na.rm=T),   Catch.L50 = quantile(Catch,probs=0.25,na.rm=T),
                                                  M.mn = median(M,na.rm=T),          M.U95 = quantile(M,probs=0.975,na.rm=T),           M.U50 = quantile(M,probs=0.75,na.rm=T),
                                                                                   M.L95 = quantile(M,probs=0.025,na.rm=T),           M.L50 = quantile(M,probs=0.25,na.rm=T),
                                                  F.mn = median(F.dec.table,na.rm=T),F.U95 = quantile(F.dec.table,probs=0.975,na.rm=T), M.U50 = quantile(F.dec.table,probs=0.75,na.rm=T),
                                                                                   F.L95 = quantile(F.dec.table,probs=0.025,na.rm=T), M.L50 = quantile(F.dec.table,probs=0.25,na.rm=T),
                                                  SP.mn = median(SP,na.rm=T),        SP.U95 = quantile(SP,probs=0.975,na.rm=T),         SP.U50 = quantile(SP,probs=0.75,na.rm=T),
       
          
                                                  
                                                                                                                                                                      SP.L95 = quantile(SP,probs=0.025,na.rm=T),         SP.L50 = quantile(SP,probs=0.25,na.rm=T))
    # Biomass realizations          
    B.real <- ggplot(Exp.res,aes(x=year,y=B,group = Sim,color=Sim)) + geom_point() + #facet_wrap(~F.scenario) +
                                                                      ylab("Predicted Biomass (metric tonnes)") + xlab("") + 
                                                                      geom_hline(yintercept = c(hcr.strat$LRP,hcr.strat$USR,hcr.strat$TRP),color=c('firebrick2','green','blue')) +
                                                                      scale_color_viridis_b(alpha=0.4,end = 0.75) +
                                                                      theme_bw() + theme(legend.position = 'none')
    
    ggsave(paste0(plot_sims,"Biomass_realizations.png"),plot=B.real,height=10,width=10)              
    
    # Catch realizations          
    C.real <- ggplot(Exp.res,aes(x=year,y=Catch,group = Sim,color=Sim)) + geom_point() + #facet_wrap(~F.scenario) +
                                                                          ylab("Predicted Catch (metric tonnes)") + xlab("") + 
                                                                          #geom_hline(yintercept = c(hcr.strat$LRP,hcr.strat$USR,hcr.strat$TRP),color=c('firebrick2','green','blue')) +
                                                                          scale_color_viridis_b(alpha=0.4,end = 0.75) +
                                                                          theme_bw() + theme(legend.position = 'none')
    
    ggsave(paste0(plot_sims,"Catch_realizations.png"),plot=C.real,height=10,width=10)  
    
    
    # Biomass time series
    B.ts <- ggplot(HCR.summarized,aes(x=year,y=B.mn)) + geom_line() + 
                                                        geom_ribbon(aes(x=year,ymax = B.U95,ymin=B.L95),fill='firebrick2',alpha=0.2) + 
                                                        geom_hline(yintercept = c(hcr.strat$LRP,hcr.strat$USR,hcr.strat$TRP),color=c('firebrick2','green','blue')) +
                                                        theme_bw() + xlab("") + ylab("Biomass (metric tonnes)") 
    ggsave(paste0(plot_sims,"B_mn_ts.png"),plot=B.ts,height=10,width=10)
    
    
    # Catch time series 
    C.ts <- ggplot(HCR.summarized,aes(x=year,y=Catch.mn)) + geom_line() + 
                                                            geom_ribbon(aes(x=year,ymax = Catch.U95,ymin=Catch.L95),fill='firebrick2',alpha=0.2) + 
                                                            #geom_hline(yintercept = c(2700,7200,9000),color=c('firebrick2','green','blue')) +
                                                            theme_bw() + xlab("") + ylab("Catch (metric tonnes)") 
    ggsave(paste0(plot_sims,"C_mn_ts.png"),plot=C.ts,height=10,width=10)
    
    
    
    # Equilib.hcr <- HCR.summarized %>% dplyr::filter(year >= 75) %>% dplyr::summarise(B.mn = median(B.mn,na.rm=T),         B.U95 = median(B.U95,na.rm=T),        B.U50 = median(B.U50,na.rm=T),
    #                                                                                  B.L95 = median(B.L95,na.rm=T),        B.L50 = median(B.L50,na.rm=T),
    #                                                                                  Catch.mn = median(Catch.mn,na.rm=T), Catch.U95 = median(Catch.U95,na.rm=T),Catch.U50 = median(Catch.U50,na.rm=T),
    #                                                                                  Catch.L95 = median(Catch.L95,na.rm=T),Catch.L50 = median(Catch.L50,na.rm=T),
    #                                                                                  Rec.mn = median(Rec.mn,na.rm=T),     Rec.U95 = median(Rec.U95,na.rm=T),        Rec.U50 = median(Rec.U50,na.rm=T),
    #                                                                                  Rec.L95 = median(Rec.L95,na.rm=T),        Rec.L50 = median(Rec.L50,na.rm=T))
    # 
    # 
    # 
    # # Get Biomass per Recruit and Yield per Recruit
    # Equilib.hcr$BPR <- Equilib.hcr$B.mn/Equilib.hcr$Rec.mn
    # Equilib.hcr$YPR <- Equilib.hcr$Catch.mn/Equilib.hcr$Rec.mn
    # 
    # 
    # 
    # 
    # # Biomass Equilibrium
    # B.equil.plt<- ggplot(Equilib.dat,aes(x=F.scenario,y=B.mn)) + 
    #   geom_point()+ geom_line()+
    #   geom_ribbon(aes(x=F.scenario,ymin=B.L95, ymax=B.U95),alpha=0.05,fill='firebrick2') +
    #   theme_bw()+ylab("Equilibrium Biomass (metric tonnes)") + xlab("Explotation Rate (proportional)")
    # ggsave(paste0(plot_sims,"B_equilibrum.png"),plot=B.equil.plt,height=6,width=6)
    # # Catch Equilibrium
    # C.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=Catch.mn)) + 
    #   geom_point()+ geom_line()+
    #   geom_ribbon(aes(x=F.scenario,ymin=Catch.L95, ymax=Catch.U95),alpha=0.05,fill='firebrick2') +
    #   theme_bw()+ylab("Equilibrium Catch (metric tonnes)")  + xlab("Explotation Rate (proportional)")
    # ggsave(paste0(plot_sims,"Catch_equilibrium.png"),plot=C.equil.plt,height=6,width=6)
    # # Biomass per Recruit (Biomass)
    # BPR.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=BPR)) + 
    #   geom_point()+ geom_line()+
    #   #geom_ribbon(aes(x=F.scenario,ymin=Catch.L95, ymax=Catch.U95),alpha=0.05,fill='firebrick2') +
    #   theme_bw()+ylab("Biomass per Recruit (Biomass)")  + xlab("Explotation Rate (proportional)")
    # ggsave(paste0(plot_sims,"BPR_equilibrium.png"),plot=BPR.equil.plt,height=6,width=6)
    # 
    # 
    # # Yield per Recruit (Biomass)
    # YPR.equil.plt<-ggplot(Equilib.dat,aes(x=F.scenario,y=YPR)) + 
    #   geom_point()+ geom_line()+
    #   #geom_ribbon(aes(x=F.scenario,ymin=Catch.L95, ymax=Catch.U95),alpha=0.05,fill='firebrick2') +
    #   theme_bw()+ylab("Yield per Recruit (Biomass)")  + xlab("Explotation Rate (proportional)")
    # ggsave(paste0(plot_sims,"YPR_equilibrium.png"),plot=YPR.equil.plt,height=6,width=6)
  }
} # end if(save.results==T)

  
  
  
  
  
  
  return(list(Catch = Catch,Rec = Rec, rps = rps,B = Bio, mort = mort, g=g, gR= gR, Surp.prod = Surp.prod, r.realized = r.realized, r.no.fishing=r.no.fishing,exp.scenario = exp.scenario,f.table=f.table,Exp.res=Exp.res))
} # End of the model function
