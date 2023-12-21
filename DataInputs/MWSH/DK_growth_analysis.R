# ageing and growth analysis

# 	Explore doing a different model for growth (e.g. jonsen equation 2)
# •	Do it the Jonsen way
# •	Do it with a GAM
# •	Do it with spatial GAM in TMB
# •	Do it with spatial Jonsen in TMB
# •	Predict for 100mm
# •	Growth input for model

#### packages and data ########

require(ggplot2)
require(tidyverse)
require(gratia)
require(lme4)
library(bbmle)
library(arm)


# hydration data
load("C:/Users/keyserf/Documents/temp_data/testing_results_framework_75-90RSCS_oldMWSH.RData")
# used for mw data, to develop new MWSH model
#load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework_75-90_newareas_issue120.RData")

#load("C:/Users/keyserf/Documents/temp_data/testing_results_spring2022_2.Rdata")

nickname<-NULL

french <- F # set to T for french
if(french == T) {
  nickname <- paste0(nickname, "_fr")
}
nickname2 <- nickname

require(rosettafish)
funs <- c("https://raw.githubusercontent.com/freyakeyser/rosetta_shell/main/terms.csv")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs)
{
  download.file(fun,destfile = basename(fun))
  rosetta_terms <- read.csv(paste0(getwd(),"/",basename(fun)), encoding = "UTF-8")
  file.remove(paste0(getwd(),"/",basename(fun)))
}

plotsGo <- "Y:/Offshore/Assessment/Framework/SFA_25_26_2024/DataInputs/MWSH"
direct_fns <- "D:/Github/Assessment_fns/"
funs <- c("https://raw.githubusercontent.com/freyakeyser/Assessment_fns/master/Survey_and_OSAC/condFac.r")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs)
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
} # end for(un in funs)

# A qqplot function
qqplot.data <- function (dat,facet=F, ncol=NULL,...)
  # argument: vector of numbers if facet = F, or if facet =T two columns first being a vector of numbers and a second column with faceting variable
{
  # following four lines from base R's qqline()
  if(facet == F)
  {
    y <- quantile(dat[!is.na(dat)], c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]
    d <- data.frame(resids = dat)
    p <- ggplot(d, aes(sample = resids)) +
      stat_qq(size=0.5) +
      geom_abline(slope = slope, intercept = int) +
      xlab(en2fr("Theoretical", custom_terms=rosetta_terms, translate=french)) +
      ylab(en2fr("Sample", custom_terms=rosetta_terms, translate=french))
  } # end if is.null(facet)

  if(facet == T)
  {
    names(dat) <- c("var","facet")
    facets <- unique(dat[,2])
    n.facets <- length(facets)
    dat.res <- NULL
    for(i in 1:n.facets)
    {
      dt <- dat %>% dplyr::filter(facet == facets[i])
      y <- quantile(dt$var[!is.na(dt$var)], c(0.25, 0.75))
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y)/diff(x)
      int <- y[1L] - slope * x[1L]
      dat.res[[i]] <- data.frame(var = dt$var,facet = dt$facet,slope = rep(slope,nrow(dt)),int = rep(int,nrow(dt)))
    }
    dat.res <- do.call("rbind",dat.res)
    p <- ggplot(dat.res, aes(sample = var)) +
      stat_qq(size=0.5) +
      geom_abline(aes(slope = slope, intercept = int)) +
      facet_wrap(~facet, ncol=ncol) +
      xlab(en2fr("Theoretical", custom_terms=rosetta_terms, translate=french)) +
      ylab(en2fr("Sample", custom_terms=rosetta_terms, translate=french))

    }# end facet == T
  print(p)
} # end fun

# Start analysis


banker <- "Sab"
dat <- mw.dat.all[[banker]]
dat$year <- as.numeric(dat$year)
# with depth across all years (random effect is ID)
sub <- dat[complete.cases(dat),]

# let's only model sh >= 65 mm here, that gets rid of just 150 shells (going to 70 drops us by 750 which seems like too much)
sub <- sub %>% dplyr::filter(sh >= 65)
# centering this makes no difference for our predictions but if make 100 = 0, then our intercept is our condition
# which may come in handy later (or may not)
sh.cond <- 100
sub$log.sh.cen <- log(sub$sh) - log(sh.cond)
# Same idea here, if we center depth on the bank median, then our intercept is the condition at the median
# Bank depth, which should be our measure of condition
#Going to say median depth of the survey tows between 2010 and 2022 see above
# For BBn this is 75 meters, which makes loads of sense. So Center the depth at the bank mean
med.depth <- sub %>% dplyr::filter(year %in% 2010:2022) %>% dplyr::summarise(med = median(depth,na.rm=T)) %>% signif(digits=2)
sub$depth.cen <- sub$depth - med.depth$med
# Also need to update the survey data.
all.surv.dat$depth.cen <- all.surv.dat$depth - med.depth$med

# We need tow to be unique by year, ooh, so we have the 0 tow problem here, need
# to make up something that gets us a unique tow for the 0s..
sub$temp_ID <- as.factor(paste(sub$year,sub$lon,sub$lat,sep = "_"))
unique_tows <- unique(sub$temp_ID)
num.levels <- paste0("x",1:length(unique_tows))
sub$new_ID <- factor(sub$temp_ID,levels = unique_tows,labels=num.levels)
sub$new_ID <- as.character(sub$new_ID)
# Use the survey tow ID's where they exist and use the fake new_ID's where they don't.
sub$new_ID[sub$tow != "0"] <- sub$tow[sub$tow != "0"]
# 1991  only has 3 tows so is problematic for many of our models, so let's start in 1992, if we want we can
# start in 1994 as that's the first year of the model.
yrs <-sort(unique(sub$year))
yrs <- yrs[yrs > 1991]
n.yrs <- length(yrs)
sub <- sub %>% dplyr::filter(year %in% yrs)
# Color scheme...
blues <- "#0057B7"
yellows <- "#FFDD00"
# Let's do this from simple to complex First a simple MW-SH model, I also don't see any reason to log transform depth
# That just a weird add on, even centering it really isn't necessary, it's pretty well behaved.



########################################################################
# Start the analysis here. I think the logical flow for this section is
# !!!!FREYA KEYSER!!!!!! look here for what I think is a logical flow for this section.
# 1: We attempted to model the data using a single generalized mixed effects model, but unfortunately the model was unable to converge when attempting to included
#    annual variability in the intercept or slope of the MW-SH relationship or if we tried to include depth.
# 2: In lieu of this we fit several simple general linear models, the results from these models, we see the AICtab results in a table that indicate that
#    allowing SH and depth to vary by year is clearly preferred over more simple models. The residual patterns of these models weren't brilliant, but we're terrible
# 3: Given that scallop collected from a survey tow are not 'independent' and this needs to be accounted for, since we can't do this using one glmer
#    we decided to develop separate MW-SH models for each year (cite inshore for this too)
#    using a mixed effect framework with tow as a random effect (intercept only) and having both SH and a linear depth term included in the models as fixed effects.
# 4: For the condition, by centering the sh at 100 mm (i.e. 0 = 100mm) and centering the depth at the median for the bank (75m =0), the intercept
#    of the glmer is actually our condition estimate at 100 mm. The other nice thing about that is we now are able to get an estimate of uncertainty
#    right out of the model, it is dead simple solution and I think it's iron clad.
# 4: I think the AIC table for the glms should be fine to summarize the obvious need for allowing slope and depth to vary by year.
# 4: The figures I think we should include in the Res Doc have all been given a name, and I'm open to adding more figures if you think we should.

####################################################################



################## Section 1  Section 1  Section 1  ################## ################## ################## ##################
# This section was useful to show that the Gamma models have better residual structure than the Gaussian models, but we can't do much with them because
# the more complex models we'd need don't converge.
# So here I'm deviating from ZUUR who suggests making a highly complex fixed effect structure before exploring the
# random effects, but from a understanding what we are doing perspective I find this much more intuitive
# So we have million options, clearly the model we have has some issues.  We know we have a random tow effect, so let's move to lmer...
mod.lmer <- lmer(log(wmw) ~ log.sh.cen + (1|new_ID),data=sub, na.action = na.omit)
# or glmer, more statistically correct
mod.glmer <- glmer(wmw ~ log.sh.cen + (1|new_ID),data=sub, na.action = na.omit,family = gaussian(link=log))
# Get data to compare
diag.lmer <- data.frame(wmw = sub$wmw,sh = sub$sh,year = sub$year,depth=sub$depth,tow = sub$new_ID,residuals = residuals(mod.lmer))
diag.glmer <- data.frame(wmw = sub$wmw,sh = sub$sh,year = sub$year,depth=sub$depth,tow = sub$new_ID,residuals = residuals(mod.glmer))
# Normality remains a fail
qqplot.data(diag.lmer$residuals)
qqplot.data(diag.glmer$residuals)
# That year effect gets swallowed up with the tow effect, more variability later, but I think that's due to increased sampling.
ggplot(data=diag.glmer,aes(group=year,y=residuals)) + geom_boxplot() + geom_hline(yintercept=0,linetype='dashed',color= blues)
# This is an improvement, but I'd much rather not see the variance blowing up like this at the higher MWs
ggplot(data=diag.glmer,aes(x=log(wmw),y=residuals)) + geom_point() + geom_smooth(method='lm')
ggplot(data=diag.glmer,aes(x=log(sh),y=residuals)) + geom_point() + geom_smooth(method='lm')

# What about depth effects? Again it seems like the tow random term is helping here.
ggplot(data=diag.glmer,aes(depth,y=residuals)) +  geom_boxplot(aes(group = cut_width(depth, 5)))
# What about tow effects. can see that these are way better now, some tows are really highly variable.
ggplot(data=diag.glmer,aes(x=tow,y=residuals)) + geom_point()
# The wmw v residual plot is the problem left with our fits, perhaps the Gamma will help us here.
mod.gamma.glmer <- glmer(wmw ~ log.sh.cen + (1|new_ID),data=sub, na.action = na.omit,family = Gamma(link=log))
summary(mod.gamma.glmer)
diag.gamma.glmer <- data.frame(wmw = sub$wmw,sh = sub$sh,year = sub$year,depth=sub$depth,tow = sub$new_ID,residuals = residuals(mod.gamma.glmer))
# Again not very good, but it's a Gamma so I start to worry less
qqplot.data(diag.gamma.glmer$residuals)
# That year effect gets swallowed up with the tow effect, more variability later, but I think that's due to increased sampling.
ggplot(data=diag.gamma.glmer,aes(group=year,y=residuals)) + geom_boxplot() + geom_hline(yintercept=0,linetype='dashed',color= blues)
# This may be an improvement, not perfect as there is still a trend, but without going to a GAM the log-log power relationship probably isn't gonna do much better
ggplot(data=diag.gamma.glmer,aes(x=log(wmw),y=residuals)) + geom_point() + geom_smooth(method='lm')
ggplot(data=diag.gamma.glmer,aes(x=log(sh),y=residuals)) + geom_point() + geom_smooth(method='lm')
# What about depth effects? Again it seems like the tow random term is helping here.
ggplot(data=diag.gamma.glmer,aes(depth,y=residuals)) +  geom_boxplot(aes(group = cut_width(depth, 5)))
# What about tow effects. can see that these are way better now, some tows are really highly variable.
ggplot(data=diag.gamma.glmer,aes(x=tow,y=residuals)) + geom_point()
# Now we could make the random effect more complex and fit the MW-SH relationship for each tow like this.
mod.gamma.re.slope.glmer <- glmer(wmw ~ log.sh.cen + (1+ log.sh.cen|new_ID),data=sub, na.action = na.omit,family = Gamma(link=log))
summary(mod.gamma.re.slope.glmer)


diag.gamma.re.slope.glmer <- data.frame(wmw = sub$wmw,sh = sub$sh,year = sub$year,depth=sub$depth,tow = sub$new_ID,residuals = residuals(mod.gamma.re.slope.glmer))
# Again not very good, but it's a Gamma so I start to worry less
qqplot.data(diag.gamma.re.slope.glmer$residuals)
# That year effect gets swallowed up with the tow effect, more variability later, but I think that's due to increased sampling.
ggplot(data=diag.gamma.re.slope.glmer,aes(group=year,y=residuals)) + geom_boxplot() + geom_hline(yintercept=0,linetype='dashed',color= blues)
# Really no different from previous model, mostly I think because the issue is the power relationship breaks down for large SH scallop
ggplot(data=diag.gamma.re.slope.glmer,aes(x=log(wmw),y=residuals)) + geom_point() + geom_smooth(method='lm')
ggplot(data=diag.gamma.re.slope.glmer,aes(x=log(sh),y=residuals)) + geom_point() + geom_smooth(method='lm')
# What about depth effects? Again it seems like the tow random term is helping here.
ggplot(data=diag.gamma.re.slope.glmer,aes(depth,y=residuals)) +  geom_boxplot(aes(group = cut_width(depth, 5)))
# What about tow effects. can see that these are way better now, some tows are really highly variable.
ggplot(data=diag.gamma.re.slope.glmer,aes(x=tow,y=residuals)) + geom_point()

# We probably aren't going to improve much on this model, but our trouble is we can't predict on tow, so as formulated our
# model is kinda ok for fitting existing data, but not very helpful for prediction at unknown locations in a given year
# The logic of using this model is fairly simple... We want to allow each tow to have it's own MW-SH relationship
# the belief being that growth will differ across the area.
# We can easily argue the converse that this is chasing noise in the data and we can't estimate a MW-SH relationship for each tow, thus
# we assume a 'global' slope of the MW-SH relationship and just allow the intercept to change for each tow.
# The results are very similar.
# From an applied perspective I think either option is reasonable, thus I'll suggest we go for the more simple model
# it also looks like we have convergence issues with a full model.

# Now we can follow the Zuur approach, which starts with the most complex model that makes sense which I think
# is a model in which the slope and intercept of the MW-SH relationship can change each year and the effect of depth can
# change each year.

# BUT THESE MODELS aren't able to converge... so the full models are out the window!!
# mod.gamma.full.glmer <- glmer(wmw ~ log.sh.cen* as.factor(year) + depth.cen* as.factor(year)  + (1|new_ID),data=sub,
#                               na.action = na.omit,family = Gamma(link=log))
# # Now there are 4 models we can compare here, drop the depth changing by year
# mod.gamma.dep.static.glmer <- glmer(wmw ~ log.sh.cen* as.factor(year) + depth.cen  + (1|new_ID),data=sub,
#                                     na.action = na.omit,family = Gamma(link=log))
# # Or drop the slope changing by year
# mod.gamma.slope.static.glmer <- glmer(wmw ~ log.sh.cen + depth.cen* as.factor(year)  + (1|new_ID),data=sub,
#                                       na.action = na.omit,family = Gamma(link=log))
# # We could get rid of the depth term altogether but keep slope of MW-SH varying by year
# mod.gamma.no.depth.glmer <- glmer(wmw ~ log.sh.cen* as.factor(year)  + (1|new_ID),data=sub,
#                                   na.action = na.omit,family = Gamma(link=log))
# # We could get rid of the depth term altogether and just allow year to influence intercept of MW-SH
# mod.gamma.no.depth.no.slope.vary.glmer <- glmer(wmw ~ log.sh.cen+ as.factor(year)  + (1|new_ID),data=sub,
#                                                 na.action = na.omit,family = Gamma(link=log))


# So as fun as that exploration was, the answer is full models are too complex



################### This is Section 2 the simple generalized linear models ################################
# First, let's compare Gamma and Gaussian models here to see if obvious choice from above holds
mod.glm.gamma <- glm(wmw ~ log.sh.cen ,data=sub, na.action = na.omit,family = Gamma(link=log))
mod.glm.gauss <- glm(wmw ~ log.sh.cen ,data=sub, na.action = na.omit,family = gaussian(link=log))

# How about res vs potential covariates
diag.glm <- data.frame(wmw = rep(sub$wmw,2),sh = rep(sub$sh,2),year = rep(sub$year,2),
                       depth=rep(sub$depth,2),tow = paste(rep(sub$new_ID,2),rep(sub$year,2)),
                       mod = sort(rep(c("Gamma","Gaussian"),length(mod.glm.gamma$residuals))),
                       residuals = c(mod.glm.gamma$residuals,mod.glm.gauss$residuals))
# Now are either of these better than the other?
# How about the normality...
# Compare the qq plots
qqplot.data(diag.glm[,c(7,6)],facet=T) # Not in love but not terrible and shockingly similar compared to the random effects variants of the same model, weird!
# They both are not good...
ggplot(diag.glm) + geom_point(aes(x=wmw,y=residuals)) + facet_wrap(~mod)
# Very similar, no real issue here.
ggplot(diag.glm) + geom_point(aes(x=sh,y=residuals)) + facet_wrap(~mod)
# Tows remain an issue that need dealt with (i.e. random effects)
ggplot(diag.glm) + geom_point(aes(x=tow,y=residuals)) + facet_wrap(~mod)

# What about linear model log transformed... it really sucks...
mod.lm <- lm(log(wmw) ~ log.sh.cen ,data=sub, na.action = na.omit)
# Now look at model diagnostics.
# Combine the residuals into the sub data to compare against
diag.lm <- data.frame(wmw = sub$wmw,sh = sub$sh,year = sub$year,depth=sub$depth,tow = sub$new_ID,residuals = mod.lm$residuals)
# Normality is a big fail for a straight linear model, interesting how much worse it is than the Gaussian glm eh!
qqplot.data(diag.lm$residuals)
# So there is a year effect
ggplot(data=diag.lm,aes(group=year,y=residuals)) + geom_boxplot() + geom_hline(yintercept=0,linetype='dashed',color= blues)
# Clear residual trends in those residuals, so that's not good.
ggplot(data=diag.lm,aes(x=log(wmw),y=residuals)) + geom_point()
ggplot(data=diag.lm,aes(x=log(sh),y=residuals)) + geom_point() + geom_smooth(method='lm') # This is not terrible TBH

# What about depth effects? Something going on, see the shallow depths are all biased positive and deeper are mostly negatives.
ggplot(data=diag.lm,aes(depth,residuals)) + geom_boxplot(aes(group = cut_width(depth, 5)))
# What about tow effects. clearly there is vaiability with the tows, could be a year effect or depth effect, or something else.
ggplot(data=diag.lm,aes(x=tow,y=residuals)) + geom_point()




# So the takeaway from the above is we should use a glm
# Since the glms are clearly better than the lm models, lets compare the glm's with Gamma
mod.full.glm <- glm(wmw ~ log.sh.cen*as.factor(year) + depth.cen*as.factor(year),data=sub, na.action = na.omit,family=Gamma(link="log"))
mod.2.glm<- glm(wmw ~ log.sh.cen*as.factor(year) + depth.cen,data=sub, na.action = na.omit,family=Gamma(link="log"))
mod.3.glm <- glm(wmw ~ log.sh.cen+as.factor(year) + depth.cen*as.factor(year),data=sub, na.action = na.omit,family=Gamma(link="log"))
mod.4.glm <- glm(wmw ~ log.sh.cen*as.factor(year) ,data=sub, na.action = na.omit,family=Gamma(link="log"))
mod.5.glm <- glm(wmw ~ log.sh.cen+as.factor(year) ,data=sub, na.action = na.omit,family=Gamma(link="log"))
mod.6.glm <- glm(wmw ~ log.sh.cen + depth.cen,data=sub, na.action = na.omit,family=Gamma(link="log"))


#########################
# So from all of this in section 2, I think the only thing we need to show as a result is this AIC Table
# Wow is that not even close, the full model kills it, suggesting both depth and sh slope should be allowed to vary by year
# I did this with the lm models and effectively got the same result as this, which is reassuring!
AIC.comp <- as.data.frame(AICtab(mod.full.glm,mod.glm.gamma,mod.2.glm,mod.3.glm,mod.4.glm,mod.5.glm,mod.6.glm))
AIC.comp$Model <- row.names(AIC.comp)
AIC.comp$formula <- NA
for(i in 1:nrow(AIC.comp)){
  AIC.comp$formula[i] <- deparse1(formula(get(AIC.comp$Model[i])))
}

write.csv(x = AIC.comp, paste0(plotsGo, "/", banker, "/AICtable.csv"))

#########################

# Lets see what the full model diagnostics look like, might be ok...
diag.full.glm <- data.frame(wmw = sub$wmw,sh = sub$sh,year = sub$year,depth=sub$depth,tow = paste(sub$new_ID,sub$year),residuals = mod.full.glm$residuals)
# Normality is still a fail, but also not a total disaster
qqplot.data(diag.full.glm$residuals)
# This ain't awful
ggplot(data=diag.full.glm,aes(group=year,y=residuals)) + geom_boxplot(aes(group = cut_width(year, 1))) + geom_hline(yintercept=0,linetype='dashed',color= blues)
# the WMW vs residuals remain awful, but better than they've been
ggplot(data=diag.full.glm,aes(x=log(wmw),y=residuals)) + geom_point()
ggplot(data=diag.full.glm,aes(x=sh,y=residuals)) + geom_point() + geom_smooth(method='lm') # This is not terrible TBH
# Better I think, but the deep stuff seems biased high, but once we bring tow in we'll see things are ok
ggplot(data=diag.full.glm,aes(depth,residuals)) + geom_boxplot(aes(group = cut_width(depth, 5)))
# What about tow effects. clearly even our full model isn't able to deal well with these.
ggplot(data=diag.full.glm,aes(x=tow,y=residuals)) + geom_point()
### End point 2 ####


#################### Section 3, YEARLY model with random effects #############################
# But we need to include the random effects of tow as they matter, thus we'll be fitting a model for each year. The above exploration has shown that the model fits
# are much better when depth and slope of the MW-SH relationship can vary each year.  The RE part shows that the Gamma does fit the
# full dataset better than the Gaussian (he glm models are much less conclusive about this, but the Gamma is never worse)
# it seems having the intercept/slope vary for each tow is effectively equivalent, so might as well keep things simple
# So it's a glmer by year, I will simplify the model to be intercept only for the random effect (this also helps as the models don't converge in several of the early years
# if we try to estimate that slope for each tow)
mod.res <- NULL
resids <- NULL
qq.plt <- NULL
for(i in 1:n.yrs)
{
  print(yrs[i])
  dat.tmp <- sub %>% dplyr::filter(year == yrs[i])
  n.tows <- length(unique(dat.tmp$new_ID))
  # Could have had more complex depth smooth, but our residuals look fine with the linear smooth so sticking with that.
  #if(yrs[i] >= 2011) mod.res[[as.character(yrs[i])]] <- glmer(wmw ~ log.sh.cen + poly(depth.cen,3) + (1| new_ID),data = dat.tmp,family=Gamma(link=log))
  # No non-linearity if we don't have 10 tows to use, to complex...
  mod.res[[as.character(yrs[i])]] <- glmer(wmw ~ log.sh.cen + depth.cen + (1| new_ID),data = dat.tmp,family=Gamma(link=log))
  dat.tmp$residuals <- residuals(mod.res[[as.character(yrs[i])]])
  resids[[as.character(yrs[i])]] <- dat.tmp
  #qq.plt[[as.character(yrs[i])]] <- qqplot.data(dat.tmp$residuals)
}
# 1991 doesn't converge and we don't need it for the models, so I suggest we only use 1992 onwards for this.
resid <- do.call('rbind',resids)

# Here are the residual plots which I think we should show, not sure if we want the smooth on there or not?

ys <- c(1992, 2002, 2012)
depthrange <- range(resid$depth)
shrange <- range(resid$sh)
mwrange <- range(resid$wmw)

for(y in 1:length(ys)){

  if(y<length(ys)) yrange <- ys[y]:(ys[y+1]-1)
  if(y==length(ys)) yrange <- ys[y]:2022

  p.res.d <- ggplot() +
    geom_point(data=resid[resid$year %in% yrange,],aes(x=depth,y=residuals), size=0.5) +
    facet_wrap(~year, ncol=2) +
    ylim(-1,1) +
    xlim(floor(depthrange[1]),ceiling(depthrange[2]))+
    geom_hline(yintercept = 0,color=blues,linetype='dashed') +
    geom_smooth(method = 'loess',color=yellows) +
    xlab(paste0(en2fr("Depth",  custom_terms=rosetta_terms, translate=french), " (m)")) +
    ylab(en2fr("Residual",  custom_terms=rosetta_terms, translate=french)) + theme_bw()

  p.res.sh <- ggplot(resid[resid$year %in% yrange,],aes(x=sh,y=residuals)) + geom_point(size=0.5) +
    facet_wrap(~year, ncol=2) +
    ylim(-1,1) +
    xlim(floor(shrange[1]),ceiling(shrange[2]))+
    geom_hline(yintercept = 0,color=blues,linetype='dashed')+
    geom_smooth(method = 'gam',color=yellows) +
    xlab(paste0(en2fr("Shell height",  custom_terms=rosetta_terms, translate=french), " (mm)")) +
    ylab(en2fr("Residual",  custom_terms=rosetta_terms, translate=french)) + theme_bw()

  p.res.mw <- ggplot(resid[resid$year %in% yrange,],aes(x=wmw,y=residuals)) + geom_point(size=0.5) +
    facet_wrap(~year, ncol=2) +
    ylim(-1,1) +
    xlim(floor(mwrange[1]),ceiling(mwrange[2]))+
    geom_hline(yintercept = 0,color=blues,linetype='dashed') +
    geom_smooth(method = 'gam',color=yellows) +
    xlab(paste0(en2fr("Meat weight",  custom_terms=rosetta_terms, translate=french), " (g)")) +
    ylab(en2fr("Residual",  custom_terms=rosetta_terms, translate=french)) + theme_bw()#+
    #scale_x_continuous(expand=c(0.075,0.075))

  png(filename = paste0(plotsGo, "/", banker, "/MWSH_resid_depth_", y, nickname, ".png"), height=6, width=5, units="in", res=420)
  print(p.res.d)
  dev.off()

  png(filename = paste0(plotsGo, "/", banker, "/MWSH_resid_sh_", y, nickname, ".png"), height=6, width=5, units="in", res=420)
  print(p.res.sh)
  dev.off()

  png(filename = paste0(plotsGo, "/", banker, "/MWSH_resid_mw_", y, nickname, ".png"), height=6, width=5, units="in", res=420)
  print(p.res.mw)
  dev.off()
}

#landscape versions for presentation purposes

p.res.d <- ggplot() +
  geom_point(data=resid,aes(x=depth,y=residuals), size=0.5) +
  facet_wrap(~year, ncol=6) +
  geom_hline(yintercept = 0,color=blues,linetype='dashed') +
  geom_smooth(method = 'loess',color=yellows) +
  xlab(paste0(en2fr("Depth",  custom_terms=rosetta_terms, translate=french), " (m)")) +
  ylab(en2fr("Residual",  custom_terms=rosetta_terms, translate=french)) + theme_bw()

p.res.sh <- ggplot(resid,aes(x=sh,y=residuals)) + geom_point(size=0.5) +
  facet_wrap(~year, ncol=6) +
  geom_hline(yintercept = 0,color=blues,linetype='dashed')+
  geom_smooth(method = 'gam',color=yellows) +
  xlab(paste0(en2fr("Shell height",  custom_terms=rosetta_terms, translate=french), " (mm)")) +
  ylab(en2fr("Residual",  custom_terms=rosetta_terms, translate=french)) + theme_bw()

p.res.mw <- ggplot(resid,aes(x=wmw,y=residuals)) + geom_point(size=0.5) +
  facet_wrap(~year, ncol=6) +
  geom_hline(yintercept = 0,color=blues,linetype='dashed') +
  geom_smooth(method = 'gam',color=yellows) +
  xlab(paste0(en2fr("Meat weight",  custom_terms=rosetta_terms, translate=french), " (g)")) +
  ylab(en2fr("Residual",  custom_terms=rosetta_terms, translate=french)) + theme_bw()+
  scale_x_continuous(expand=c(0.075,0.075))

png(filename = paste0(plotsGo, "/", banker, "/MWSH_resid_depth_landscape", nickname, ".png"), height=6, width=10, units="in", res=420)
print(p.res.d)
dev.off()

png(filename = paste0(plotsGo, "/", banker, "/MWSH_resid_sh_landscape", nickname, ".png"), height=6, width=10, units="in", res=420)
print(p.res.sh)
dev.off()

png(filename = paste0(plotsGo, "/", banker, "/MWSH_resid_mw_landscape", nickname, ".png"), height=6, width=10, units="in", res=420)
print(p.res.mw)
dev.off()



# And here's a facet qqplot. They aren't brilliant, but I think they are reasonable
p.qqs <- qqplot.data(data.frame(var=resid$residuals,facet = resid$year),facet = T, ncol=5)
p.qqs <- p.qqs + theme_bw()


png(filename = paste0(plotsGo, "/", banker, "/MWSH_qq", nickname, ".png"), height=6, width=5, units="in", res=420)
print(p.qqs)
dev.off()

# Ok so we have our models run, now how to do the predictions, for the tows that were sampled we use the random effects
# from the model, for the other tows we use the fixed effects only since we don't have the randoms.
# Get a vector of shs to predict on
#pred.depths <- min(sub$depth),max(sub$depth))
log.sh.cen <- log(seq(2.5, 197.5, by = 5)) - log(sh.cond) #each shell height bin to predict on, making 0 being 100 mm

# Initialize some stuff
mw.res.t <- NULL
cond.pred <- NULL
all.coef <- NULL
# Takes about 1 second per year
for(i in 1:n.yrs)
{
  mw.dat <- sub %>% dplyr::filter(year== yrs[i])
  s.dat <- all.surv.dat %>% dplyr::filter(bank==banker & year == yrs[i])

  mod.r <- mod.res[[as.character(yrs[i])]]
  # Now try and do the predictions

  #get IDs for the sampled tows. Use the random effects for those.
  random.pred <- (1:nrow(s.dat))[is.element(s.dat$tow,unique(mw.dat$new_ID))]

  #get IDs for the unsampled tows. Use fixed effects for IDs that weren't sampled for meat weight shell height
  fixed.pred <- (1:nrow(s.dat))[!is.element(s.dat$tow,unique(mw.dat$new_ID))]

  #Predict using Random effects for IDs that were sampled for meat weight shell height
  temp <- matrix(NA,nrow(s.dat),40)

  for(j in random.pred)
  {
    temp[j,] <- as.vector(predict(object = mod.r,newdata=data.frame(log.sh.cen=log.sh.cen,
                                                                    depth.cen=rep(s.dat$depth.cen[j] ,40),
                                                                    new_ID=s.dat$tow[j]),
                                  re.form=NULL,type="response"))
  } # end the random loop

  #Predict using fixed effects for IDs that weren't sampled for meat weight shell height
  for(j in fixed.pred)
  {
    temp[j,] <- as.vector(predict(object=mod.r,newdata=data.frame(log.sh.cen=log.sh.cen,
                                                                  depth.cen=rep(s.dat$depth.cen[j] ,40)),
                                  re.form=~0,type="response"))
  } # end the fixed loop

  #multply temp matrix (weight) by live numbers to get weight/size bin
  s.dat[,grep("h5", colnames(s.dat))[1]:grep("h200", colnames(s.dat))] <-
    temp*s.dat[,grep("h5", colnames(s.dat))[1]:grep("h200", colnames(s.dat))]
  mw.res.t[[i]] <- s.dat


  # So finally we need to use our models to predict condition on the bank, seemingly the easiest way is to pick a depth and MW to predict at
  # for the bank and just leave it at that.  SH will be 100 mm, Going to say median depth (which is 0 as set up above) of the survey tows between 2010 and 2022 see above
  # For BBn this is 75 meters, which makes loads of sense. Downside is only way to get an SE is to try and bootstrap one, which I'm too lazy to do
  cond.est <-   as.vector(predict(object=mod.r, newdata=data.frame(log.sh.cen=0,
                                                                   depth.cen=0),
                                  re.form=~0,type="response"))
  # We can pull out the intercept now as well and see how this compares
  inter <- summary(mod.r)$coefficients[1,1]
  inter.se <- summary(mod.r)$coefficients[1,2]
  # While we are at this, let's pull out the fixed slope and the random intercepts for the mw-sh figure and the depth terms.
  slope <- summary(mod.r)$coefficients[2,1]
  slope.se <- summary(mod.r)$coefficients[2,2]
  dep.1 <- summary(mod.r)$coefficients[3,1]
  dep.1.se <- summary(mod.r)$coefficients[3,2]

  # Now extract the random terms
  rand.coef <- data.frame(rand.int = ranef(mod.r)$new_ID[[1]], rand.se = se.ranef(mod.r)$new_ID[[1]],
                          tow = attr(se.ranef(mod.r)$new_ID,'dimnames')[[1]])
  fix.coef <- data.frame(fix.int = inter, fix.int.se = inter.se, fix.slope = slope, fix.slope.se = slope.se,
                         depth1 = dep.1, depth1.se = dep.1.se,
                         tow = as.character(sort(unique(s.dat$tow))),year = yrs[i])
  all.coef[[i]] <- left_join(fix.coef,rand.coef,by='tow')

  # Object with condition prediction
  cond.pred[[i]] <- data.frame(cond = cond.est,year = yrs[i],intercept = inter,inter.se = inter.se)


} # end the i loop
#Let's rename this from whatever lousy name we use to something that makes sense, feel free to make it better.
biomass.per.sh.by.tow <- do.call("rbind",mw.res.t)
condition.ts <- do.call('rbind',cond.pred)
mw.sh.coef <- do.call('rbind',all.coef)
# Add some variables
condition.ts$cond.inter <- exp(condition.ts$intercept)
condition.ts$cond.inter.LCI <- exp(condition.ts$intercept - 1.96*condition.ts$inter.se)
condition.ts$cond.inter.UCI <- exp(condition.ts$intercept + 1.96*condition.ts$inter.se)
# Get the actual random intercepts
mw.sh.coef$ran.int.act <- mw.sh.coef$rand.int+ mw.sh.coef$fix.int

# Can see our results are identical, so there is no need to use a predict, we just pull the intercept and say we are predicting condtion at a depth of 75 meters on Browns Bank.
ggplot(condition.ts, aes(cond,cond.inter)) + geom_point() + geom_abline(intercept=0,slope=1,color=blues)

# Use the intercept which gives us a CI around the Condition, which is fun and I think it's legit the condition of a 100mm scallop at the median depth of the bank
p.cond.inter <- ggplot(condition.ts,aes(x=year,y=cond.inter)) +
  geom_ribbon(aes(ymin=cond.inter.LCI,ymax=cond.inter.UCI,x=year),color=blues,fill=blues,alpha=0.5)+
  geom_hline(yintercept = mean(condition.ts$cond,na.rm=T),color=yellows,linetype = 'dashed',size=1.5) +
  geom_line(size=1.5) +
  geom_point(size=3) +
  ylim(c(min(condition.ts$cond.inter.LCI),max(condition.ts$cond.inter.UCI)))+
  xlab("") + ylab("Meat Weight of 100 mm Scallop (grams)") +
  theme_bw()

# So now we need a new MW-SH figure and for the Res Doc we probably should show the Depth covariate, or at least
# create that figure.
slope <- mw.sh.coef %>% dplyr::filter(year == max(yrs)) %>% dplyr::pull(fix.slope)
int <- mw.sh.coef %>% dplyr::filter(year == max(yrs)) %>% dplyr::pull(fix.int)
rand.int <- mw.sh.coef[mw.sh.coef$year==max(yrs) & !is.na(mw.sh.coef$ran.int.act),] %>% dplyr::pull(ran.int.act,tow)
sh <- 65:200/100
# So now we can make the curve for every tow and our overall fixed effect.
rand.tow <- NULL
for(i in 1:length(rand.int)) rand.tow[[i]] <- data.frame(sh = 100*sh,mw = exp(rand.int)[i] * sh^slope[1],tow = names(rand.int)[i])
rand.tows <- do.call("rbind",rand.tow)
fix.mw <- data.frame(sh = 100* sh, mw = exp(int)[1] * sh^slope[1])

# Based on our current figure, I don't think this goes in RES DOC, but is basically what we want for our SS.
p.ss <- ggplot(rand.tows) +
  geom_point(data=sub %>% dplyr::filter(year == max(yrs)),aes(x=sh,y=wmw),color=blues) +
  geom_line(aes(x=sh,y=mw,group=tow),color=yellows)+
  geom_line(data = fix.mw,aes(x=sh,y=mw),color='black',size=2) +
  scale_x_continuous(limits = c(65,sub %>% dplyr::filter(year == max(yrs)) %>% dplyr::summarise(max(sh,na.rm=T)) %>% as.numeric()),
                     name = "Shell Height (mm)",
                     breaks = seq(0,200,by=10), expand = c(0.005,0.005)) +
  scale_y_continuous(limits = c(0,sub %>% dplyr::filter(year == max(yrs)) %>% dplyr::summarise(max(wmw,na.rm=T)) %>% as.numeric()),
                     name = "Meat Weight (grams)",
                     breaks = seq(0,100,by=10)) +
  theme_bw()

# I was thinking we could show the line with some uncertainty instead of showing the tows, but
# that's not so easy since we only have the parameter uncertainties and not the uncertainty of the line itself (again bootstraping may be an option, but
# I'm kinda meh on that, I think the above is fine).
# For the Res Doc we may need a MW-SH facet plot, given the number of years, I'm gonna say we make three of them, 1992-2001, 2002-2011, 2012-2022

sh <- 65:200/100
r.tows <- NULL
f.mws <- NULL
for(j in 1:n.yrs)
{
  rt <- NULL
  slope <- mw.sh.coef %>% dplyr::filter(year == yrs[j]) %>% dplyr::pull(fix.slope); slope <- slope[1]
  int <- mw.sh.coef %>% dplyr::filter(year == yrs[j]) %>% dplyr::pull(fix.int); int <- int[1]
  rand.int <- mw.sh.coef[mw.sh.coef$year==yrs[j] & !is.na(mw.sh.coef$ran.int.act),] %>% dplyr::pull(ran.int.act,tow)
  for(i in 1:length(rand.int)) rt[[i]] <- data.frame(sh = 100*sh,mw = exp(rand.int)[i] * sh^slope,tow = names(rand.int)[i])
  rts <- do.call("rbind",rt)
  r.tows[[j]] <- data.frame(rts,year = yrs[j])
  f.mws[[j]] <- data.frame(sh = 100* sh, mw = exp(int) * sh^slope,year = yrs[j])

}

r.tow <- do.call('rbind',r.tows)
f.mw <- do.call('rbind',f.mws)
y1 <- 1992:2001
y2 <- 2002:2011
y3 <- 2012:2022

# SO I THINK WE WANT THESE THREE FIGURES FOR RES DOC
# 1992-2001
p1 <- ggplot(r.tow %>% dplyr::filter(year %in% y1)) +
  geom_point(data=sub %>% dplyr::filter(year %in% y1),aes(x=sh,y=wmw),color=blues) +
  geom_line(aes(x=sh,y=mw,group=tow),color=yellows, size=0.25)+ facet_wrap(~year, ncol=2) +
  geom_line(data = f.mw %>% dplyr::filter(year %in% y1),aes(x=sh,y=mw),color='black',size=2) +
  scale_x_continuous(limits = c(65,mw.dat %>% dplyr::filter(year == max(yrs)) %>% dplyr::summarise(max(sh,na.rm=T)) %>% as.numeric()),
                     name = paste0(en2fr("Shell height",  custom_terms=rosetta_terms, translate=french), " (mm)"),
                     breaks = seq(0,200,by=20), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(0,sub %>% dplyr::summarise(max(wmw,na.rm=T)) %>% as.numeric()),
                     name = paste0(en2fr("Meat weight",  custom_terms=rosetta_terms, translate=french), " (g)"),
                     breaks = seq(0,100,by=20)) +
  theme_bw()

# 2002-2011
p2 <- ggplot(r.tow %>% dplyr::filter(year %in% y2)) +
  geom_point(data=sub %>% dplyr::filter(year %in% y2),aes(x=sh,y=wmw),color=blues) +
  geom_line(aes(x=sh,y=mw,group=tow),color=yellows, size=0.25)+ facet_wrap(~year, ncol=2) +
  geom_line(data = f.mw %>% dplyr::filter(year %in% y2),aes(x=sh,y=mw),color='black',size=2) +
  scale_x_continuous(limits = c(65,mw.dat %>% dplyr::filter(year == max(yrs)) %>% dplyr::summarise(max(sh,na.rm=T)) %>% as.numeric()),
                     name = paste0(en2fr("Shell height",  custom_terms=rosetta_terms, translate=french), " (mm)"),
                     breaks = seq(0,200,by=20), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(0,sub %>% dplyr::summarise(max(wmw,na.rm=T)) %>% as.numeric()),
                     name = paste0(en2fr("Meat weight",  custom_terms=rosetta_terms, translate=french), " (g)"),
                     breaks = seq(0,100,by=20)) +
  theme_bw()

# Now 2012-2022
p3 <- ggplot(r.tow %>% dplyr::filter(year %in% y3)) +
  geom_point(data=sub %>% dplyr::filter(year %in% y3),aes(x=sh,y=wmw),color=blues) +
  geom_line(aes(x=sh,y=mw,group=tow),color=yellows, size=0.25)+ facet_wrap(~year, ncol=2) +
  geom_line(data = f.mw %>% dplyr::filter(year %in% y3),aes(x=sh,y=mw),color='black',size=2) +
  scale_x_continuous(limits = c(65,mw.dat %>% dplyr::filter(year == max(yrs)) %>% dplyr::summarise(max(sh,na.rm=T)) %>% as.numeric()),
                     name = paste0(en2fr("Shell height",  custom_terms=rosetta_terms, translate=french), " (mm)"),
                     breaks = seq(0,200,by=20), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(0,sub %>% dplyr::summarise(max(wmw,na.rm=T)) %>% as.numeric()),
                     name = paste0(en2fr("Meat weight",  custom_terms=rosetta_terms, translate=french), " (g)"),
                     breaks = seq(0,100,by=20)) +
  theme_bw()

pall <- ggplot(r.tow %>% dplyr::filter(year %in% 1992:2022)) +
  geom_point(data=sub %>% dplyr::filter(year %in% 1992:2022),aes(x=sh,y=wmw),color=blues, size=0.5) +
  geom_line(aes(x=sh,y=mw,group=tow),color=yellows, size=0.25)+ facet_wrap(~year, ncol=6) +
  geom_line(data = f.mw %>% dplyr::filter(year %in% 1992:2022),aes(x=sh,y=mw),color='black',size=1) +
  scale_x_continuous(limits = c(65,mw.dat %>% dplyr::filter(year == max(yrs)) %>% dplyr::summarise(max(sh,na.rm=T)) %>% as.numeric()),
                     name = paste0(en2fr("Shell height",  custom_terms=rosetta_terms, translate=french), " (mm)"),
                     breaks = seq(0,200,by=25), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(0,sub %>% dplyr::summarise(max(wmw,na.rm=T)) %>% as.numeric()),
                     name = paste0(en2fr("Meat weight",  custom_terms=rosetta_terms, translate=french), " (g)"),
                     breaks = seq(0,100,by=25)) +
  theme_bw()

# Ok so final thing is to show the depth effects, I was gonna get fancy but why bother.. it's easy...

depth.effect <- data.frame(effect = NA,se=NA,year=NA)
for(i in 1:n.yrs) depth.effect[i,] <- data.frame(effect = mw.sh.coef$depth1[mw.sh.coef$year == yrs[i]][1],se = mw.sh.coef$depth1.se[mw.sh.coef$year == yrs[i]][1],year=yrs[i])
depth.effect$UCI <- depth.effect$effect + 1.96*depth.effect$se
depth.effect$LCI <- depth.effect$effect - 1.96*depth.effect$se

# Plot of the effect of depth.  In most years the slope is negative as would be expected, but it is infrequently significantly negative in any one year.
# But overall given it's always negative, there is no doubt a general decrease in mw with increasing depth.
if(french==T) depthlab <- "Effet de profondeur"
if(french==F) depthlab <- "Depth effect"
p.depth <- ggplot(depth.effect,aes(x=year,y=effect)) +
  geom_ribbon(aes(ymin=LCI,ymax=UCI,x=year),color=blues,fill=blues,alpha=0.5)+
  geom_hline(yintercept = 0,color=yellows,linetype = 'dashed',size=1.5) +
  #geom_line(size=1.5) + geom_point(size=3) +
  geom_line() + geom_point() +
  #ylim(c(min(condition.ts$cond.inter.LCI),max(condition.ts$cond.inter.UCI)))+
  xlab("") +
  ylab(depthlab) +
  theme_bw()

png(filename = paste0(plotsGo, "/", banker, "/MWSH_y1", nickname, ".png"), height=6, width=5, units="in", res=420)
print(p1)
dev.off()

png(filename = paste0(plotsGo, "/", banker, "/MWSH_y2", nickname, ".png"), height=6, width=5, units="in", res=420)
print(p2)
dev.off()

png(filename = paste0(plotsGo, "/", banker, "/MWSH_y3", nickname, ".png"), height=6, width=5, units="in", res=420)
print(p3)
dev.off()

png(filename = paste0(plotsGo, "/", banker, "/MWSH_depth", nickname, ".png"), height=4, width=5, units="in", res=420)
print(p.depth)
dev.off()

png(filename = paste0(plotsGo, "/", banker, "/MWSH_yall", nickname, ".png"), height=6, width=10, units="in", res=420)
print(pall)
dev.off()

### for the other banks after running in survey summary
load("C:/Users/keyserf/Documents/temp_data/testing_results_framework_75-90RSCS_newMWSH_GBb.RData")
nickname <- nickname2
banks <- c("Mid", "Ban", "Ger", "BBs", "GBb")
for (bank in banks){
  sh <- 65:200/100
  y3 <- rev(rev(sort(unique(cf.data[[bank]]$HtWt.fit$resid$year)))[1:10])
  y3 <- y3[!is.na(y3)]

  dat <- mw.dat.all[[bank]]
  # with depth across all years (random effect is ID)
  sub <- dat[complete.cases(dat),]

  r.tows <- NULL
  f.mws <- NULL
  n.yrs <- length(y3)
  mw.sh.coef <- cf.data[[bank]]$CF.fit$mw.sh.coef
  for(j in 1:n.yrs)
  {
    rt <- NULL
    slope <- mw.sh.coef %>% dplyr::filter(year == y3[j]) %>% dplyr::pull(fix.slope); slope <- slope[1]
    int <- mw.sh.coef %>% dplyr::filter(year == y3[j]) %>% dplyr::pull(fix.int); int <- int[1]
    rand.int <- mw.sh.coef[mw.sh.coef$year==y3[j] & !is.na(mw.sh.coef$ran.int.act),] %>% dplyr::pull(ran.int.act,tow)
    for(i in 1:length(rand.int)) rt[[i]] <- data.frame(sh = 100*sh,mw = exp(rand.int)[i] * sh^slope,tow = names(rand.int)[i])
    rts <- do.call("rbind",rt)
    r.tows[[j]] <- data.frame(rts,year = y3[j])
    f.mws[[j]] <- data.frame(sh = 100* sh, mw = exp(int) * sh^slope,year = y3[j])
  }

  r.tow <- do.call('rbind',r.tows)
  f.mw <- do.call('rbind',f.mws)

  png(filename = paste0(plotsGo, "/MWSH_y3_", bank, nickname, ".png"), height=6, width=5, units="in", res=420)
  print(
    ggplot(r.tow %>% dplyr::filter(year %in% y3)) +
      geom_point(data=cf.data[[bank]]$HtWt.fit$resid %>% dplyr::filter(year %in% y3),aes(x=sh,y=wmw),color=blues) +
      geom_line(aes(x=sh,y=mw,group=tow),color=yellows, size=0.25)+ facet_wrap(~year, ncol=2) +
      geom_line(data = f.mw %>% dplyr::filter(year %in% y3),aes(x=sh,y=mw),color='black',size=2) +
      scale_x_continuous(limits = c(65,cf.data[[bank]]$HtWt.fit$resid %>% dplyr::filter(year == max(y3)) %>% dplyr::summarise(max(sh,na.rm=T)) %>% as.numeric()),
                         name = paste0(en2fr("Shell height",  custom_terms=rosetta_terms, translate=french), " (mm)"),
                         breaks = seq(0,200,by=20), expand = c(0.01,0.01)) +
      scale_y_continuous(limits = c(0,sub %>% dplyr::summarise(max(wmw,na.rm=T)) %>% as.numeric()),
                         name = paste0(en2fr("Meat weight",  custom_terms=rosetta_terms, translate=french), " (g)"),
                         breaks = seq(0,100,by=20)) +
      theme_bw()
  )
  dev.off()
}

