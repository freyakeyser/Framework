# ageing and growth analysis

# 	Explore doing a different model for growth (e.g. jonsen equation 2)
# •	Do it the Jonsen way
# •	Do it with a GAM
# •	Do it with spatial GAM in TMB
# •	Do it with spatial Jonsen in TMB
# •	Predict for 100mm
# •	Growth input for model


require(ggplot2)
require(tidyverse)
require(gratia)

# hydration data
load("C:/Users/keyserf/Documents/temp_data/testing_results_framework_75-90.RData")

plotsGo <- "Y:/Offshore/Assessment/Framework/SFA_25_26_2024/DataInputs/MWSH"

bank <- "BBn"
dat <- mw.dat.all[[bank]]
head(dat)


### Jonsen 2009 Equation 2 method #############################################################################
# with condFac function
# (runs shwt.lme inside condFac, which uses nlme to fit model)
source("~/GitHub/Assessment_fns/Survey_and_OSAC/condFac.r")
mw.dat.mod <- dat
mw.dat.mod <- mw.dat.mod[complete.cases(mw.dat.mod),]
#gam_f for BBn, BBs, Ger, Sab. glm for Ban, Mid
if(!bank %in% c("Mid", "Ban")){
  condfac_res <- condFac(mw.dat.mod,model.type='gam_f',dirct=direct_fns)
  #condfac_res$CFyrs
}
if(bank %in% c("Mid", "Ban")){
  condfac_res <- condFac(mw.dat.mod,model.type='glm',dirct=direct_fns)
  #condfac_res$CFyrs
}
# note that we usually provide prediction data. I have not done that here yet!


### Jonsen 2009 Equation 2 with lme4

shwt.dat <- dat
shwt.dat$sh <- (shwt.dat$sh/100)^3
shwt.dat$ID <- paste0(shwt.dat$year, ".", shwt.dat$tow)
# fit MWSH model (replaces the shwt.lme step in condFac)
# use ID as unique ID for tow and year. Allow random slope (different from previous method)
require(lme4)
# this matches the nlme model in shwt.lme
lme4res <- lme4::lmer(data=shwt.dat, wmw ~ sh-1 + (sh-1|ID))
# a = sh in fit
fit <- data.frame(ID=row.names(coef(lme4res)$ID),
                  CF=coef(lme4res)$ID[[1]])
full.dat <- left_join(bank.dat[[bank]], fit)
# fit$sh is actually the CF estimate for each tow

#gam_f for BBn, BBs, Ger, Sab. glm for Ban, Mid
if(!bank %in% c("Mid", "Ban")){
  # gam_f which has location (lat/lon) and depth fit as thin plate regression splines, year is a factor, Gaussian family and identity link.
  CF.fit<-gam(CF~s(lon,lat)+s(depth)+as.factor(year),data=full.dat)
}

if(bank %in% c("Mid", "Ban")){
  # This model assumes CF varies only with depth and year, Gaussian and linear relationship, no random effects (year might be best treated as such)
  CF.fit<-glm(CF~depth+as.factor(year),data=full.dat)
}

years <- unique(full.dat[!is.na(full.dat$CF),]$year)
complete <- bank.dat[[bank]][bank.dat[[bank]]$year %in% years,]
complete <- cbind(ID=complete$ID, as.data.frame(predict(CF.fit, newdata=complete, se.fit=T)))
full.dat<- left_join(full.dat, complete)

annual.mean <- full.dat %>%
  filter(state=="live") %>%
  group_by(year) %>%
  summarize(mean.CF=mean(fit, na.rm=T))

# or predict in a specific location
pred.loc <-NULL
pred.loc[["depth"]] <- mean(subset(shwt.dat, year >=2005 & year <2015)$depth,na.rm=T)
pred.loc[["lat"]] <- mean(subset(shwt.dat, year >=2005 & year <2015)$lat,na.rm=T)
pred.loc[["lon"]] <- mean(subset(shwt.dat, year >=2005 & year <2015)$lon,na.rm=T)
# Make a new object to build predictions from, previously we were using predictions for mean data for our data
# But those change every year, this has been revised to predict on the same location every year.
CFyrs<-data.frame(year=years,depth=pred.loc[["depth"]],lon=pred.loc[["lon"]],lat=pred.loc[["lat"]])
CFyrs <- cbind(year=CFyrs$year, as.data.frame(predict(CF.fit, newdata=CFyrs, se.fit=T))) %>%
  arrange(year)


ggplot() + geom_line(data=condfac_res$CFyrs, aes(as.numeric(year), CF), group=1) +
  geom_line(data=CFyrs, aes(year, fit), colour="red", group=1)


### one step GAM? ##############################################################################################

dat$cf <- dat$wmw/dat$sh*100

# if(!bank %in% c("Mid", "Ban")){
#   # gam_f which has location (lat/lon) and depth fit as thin plate regression splines, year is a factor, Gaussian family and identity link.
#   CF.fit<-gam(CF~s(lon,lat)+s(depth)+as.factor(year),data=full.dat)
# }
#
# if(bank %in% c("Mid", "Ban")){
#   # This model assumes CF varies only with depth and year, Gaussian and linear relationship, no random effects (year might be best treated as such)
#   CF.fit<-glm(CF~depth+as.factor(year),data=full.dat)
# }


depth <- mean(dat[dat$year %in% 2005:2015,]$depth)
lon <- mean(dat[dat$year %in% 2005:2015,]$lon)
lat <- mean(dat[dat$year %in% 2005:2015,]$lat)
preds <- data.frame(year=unique(dat$year), depth=depth, lon=lon, lat=lat, sh=100)

mod_eval <- function(mod){

  if(!"lme" %in% names(mod)) {
    resid <- cbind(dat, as.data.frame(predict(mod, dat, se.fit=T, type = "link")))
    print(qqnorm(residuals(mod)))
    print(qqline(residuals(mod)))
    }
  if("lme" %in% names(mod)) {
    other <- mod[which(!names(mod)=="lme")][[1]]
    resid <- cbind(dat, as.data.frame(predict(other, dat, se.fit=T, type = "link")))
    print(qqnorm(residuals(other)))
    print(qqline(residuals(other)))
  }
  resid$resid.link <- resid$cf - resid$fit
  print(ggplot() + geom_point(data=resid, aes(sh, resid.link)))


  if(!"lme" %in% names(mod)) {
    res <- cbind(preds, as.data.frame(predict(mod, preds, se.fit=T, type = "link")))
    print(ggplot() + geom_line(data=res, aes(year, inv_link(mod)(fit), group=1)) +
            geom_point(data=res, aes(year, inv_link(mod)(fit))) +
            geom_line(data=res, aes(year,inv_link(mod)(fit-1.96*se.fit), group=1)) +
            geom_line(data=res, aes(year,inv_link(mod)(fit+1.96*se.fit), group=1)) +
            ggtitle(paste0(deparse(formula(mod)), ", ", mod$family$family, ", ", mod$family$link))
    )
  }
  if("lme" %in% names(mod)) {
    res <- cbind(preds, as.data.frame(predict(other, preds, se.fit=T, type = "link")))

    print(ggplot() + geom_line(data=res, aes(year, inv_link(other)(fit), group=1)) +
            geom_point(data=res, aes(year, inv_link(other)(fit))) +
            geom_line(data=res, aes(year,inv_link(other)(fit-1.96*se.fit), group=1)) +
            geom_line(data=res, aes(year,inv_link(other)(fit+1.96*se.fit), group=1)) +
            ggtitle(paste0(deparse(formula(other)), ", ranef: ", names(ranef(mod$lme)[3]), ",", other$family$family, ", ", other$family$link))
    )
  }


  return(list(preds=res, mod=mod))
}


# gam_cf <- gam(data=dat,
#               cf ~ s(lon,lat) + s(depth) + as.factor(year),
#               family=Gamma(link="log"))
# gam_cf <- mod_eval(gam_cf)


# glm_cf <- glm(data=dat,
#               cf ~ depth + as.factor(year),
#               family=Gamma(link="log"))
# glm_cf <- mod_eval(glm_cf)


# gamm_cf <- gamm(data=dat,
#                 cf ~ as.factor(year) + s(depth) + s(lon, lat),
#                 random=list(ID=~1),
#                 family=Gamma(link="log"))
# gamm_cf <- mod_eval(gamm_cf)

# with lon lat it goes high

# gam_cf, glm_cf, gamm_cf

# what if I use SH + depth + location + year to estimate MW for the sampled tows,
# and then expand the SHFs to get SH for each tow and predict MW from that?
# NO. This doesn't really make sense. Would have to make too many assumptions about SH (rounding the standardized counts and the 5mm bins)
# Need to skip MWSH and use condition.
# OR you could estimate MW, and then predict on 100mm individual. Let's try that.
# same model as above.

# include tow as random effect in all models (glmm/gamm)
# centre wmw and sh
# copy inshore mwsh model
# run across all years and run for individual years (individual years will need depth:year)
# include depth and not
# determine gaussian and gamma using residuals. Start with gaussian, and only do gamma if necessary (because resids suck)
# see yihao yin 2022 spa3 paper canjfish (can cite this for gaussian)
# https://github.com/Mar-scal/Inshore/blob/main/SFA29/Growth/SFA29_MeatWeightShellHeight.R

glm_mw <- glm(data=dat, wmw ~ sh + as.factor(year), family=gaussian(link="log"))
glm_mw <- mod_eval(glm_mw)

glm_mw_gamma <- glm(data=dat, wmw ~ sh + as.factor(year), family=Gamma(link="log"))
glm_mw_gamma <- mod_eval(glm_mw_gamma)

gam_mw <- gam(data=dat, wmw ~ sh + s(lon, lat) + s(depth) + as.factor(year), family=Gamma(link="log"))
# Gamma() too slow (with inverse)
# Gamma with log link is pretty slow
# gaussian(link="log") also very slow
gam_mw <- mod_eval(gam_mw)

# gaussian and log link Gamma are pretty similar

# pretty darn good

# gam_mw2 <- gam(data=dat, wmw ~ sh + s(lon, lat) + s(depth) + as.factor(year))
# gam_mw2 <- mod_eval(gam_mw2)
# gam for MW with gamma log link wins so far

# what about log transforming sh

# gam_mw3 <- gam(data=dat, wmw ~ log(sh) + s(lon, lat) + s(depth) + as.factor(year))
# gam_mw3 <- mod_eval(gam_mw3)
#super fast

# what about sh^3 maybe?

# gam_mw4 <- gam(data=dat, wmw ~ sh^3 + s(lon, lat) + s(depth) + as.factor(year))
# gam_mw4 <- mod_eval(gam_mw4)
#super fast

# what about sh^3 maybe? and log link?

gam_mw5 <- gam(data=dat, wmw ~ sh^3 + s(lon, lat) + s(depth) + as.factor(year), family=Gamma(link="log"))
gam_mw5 <- mod_eval(gam_mw5)
# IDENTICAL to gam2 (gamma log link without ^3)


#gam_mw6 <- gam(data=dat, wmw ~ sh + s(depth) + as.factor(year), family=gaussian(link="log"))
# seriously too slow.
glm_mw6 <- glm(data=dat, wmw ~ sh + as.factor(year), family=gaussian(link="log"))
glm_mw6 <- mod_eval(glm_mw6)
# pretty close to original

# gam_cf, glm_cf, gamm_cf, glm_mw, glm_mw_gamma, gam_mw, gam_mw2, gam_mw3, gam_mw4, gam_mw5, glm_mw6
# cols <- viridis::viridis(n=13)
# mods <- c(gam_cf, glm_cf, gamm_cf, glm_mw, glm_mw_gamma, gam_mw, gam_mw2, gam_mw3, gam_mw4, gam_mw5, glm_mw6)
#
# colors <- c("gam_cf" = cols[1],
#             "glm_cf" = cols[2],
#             "gamm_cf" = cols[3],
#             "glm_mw" = cols[4],
#             "glm_mw_gamma" = cols[5],
#             "gam_mw" = cols[6],
#             "gam_mw2" = cols[7],
#             "gam_mw3" = cols[8],
#             "gam_mw4" = cols[9],
#             "gam_mw5" = cols[10],
#             "glm_mw6" = cols[11],
#             "orig" = cols[12],
#             "orig2" = cols[13])

cols <- viridis::viridis(n=6)
colors <- c("glm_mw" = cols[1],
            "glm_mw_gamma" = cols[2],
            "gam_mw5" = cols[3],
            "glm_mw6" = cols[4],
            "orig" = cols[5],
            "orig2" = cols[6])

png(paste0(plotsGo, "/", bank, "/CF_model_compare.png"), height=4, width=6.5, units="in", res=420)
ggplot() + geom_line(data=condfac_res$CFyrs, aes(as.numeric(year), CF, colour="orig"),group=1) +
  geom_line(data=CFyrs, aes(year, fit, colour="orig2"), group=1) +
  #geom_line(data=gam_cf$preds, aes(as.numeric(year), inv_link(gam_cf$mod)(fit), colour="gam_cf"), group=1) +
  #geom_line(data=glm_cf$preds, aes(as.numeric(year), inv_link(glm_cf$mod)(fit), colour="glm_cf"), group=1) +
  #geom_line(data=gamm_cf$preds, aes(as.numeric(year), inv_link(gamm_cf$mod)(fit), colour="gamm_cf"), group=1) +
  geom_line(data=glm_mw$preds, aes(as.numeric(year), inv_link(glm_mw$mod)(fit), colour="glm_mw"), group=1) +
  geom_line(data=glm_mw_gamma$preds, aes(as.numeric(year), inv_link(glm_mw_gamma$mod)(fit), colour="glm_mw_gamma"), group=1) +
  #geom_line(data=gam_mw$preds, aes(as.numeric(year), inv_link(gam_mw$mod)(fit), colour="gam_mw"), group=1) +
  #geom_line(data=gam_mw2$preds, aes(as.numeric(year), inv_link(gam_mw2$mod)(fit), colour="gam_mw2"), group=1) +
  #geom_line(data=gam_mw3$preds, aes(as.numeric(year), inv_link(gam_mw3$mod)(fit), colour="gam_mw3"), group=1) +
  #geom_line(data=gam_mw4$preds, aes(as.numeric(year), inv_link(gam_mw4$mod)(fit), colour="gam_mw4"), group=1) +
  geom_line(data=gam_mw5$preds, aes(as.numeric(year), inv_link(gam_mw5$mod)(fit), colour="gam_mw5"), group=1) +
  geom_line(data=glm_mw6$preds, aes(as.numeric(year), inv_link(glm_mw6$mod)(fit), colour="glm_mw6"), group=1) +
  scale_color_manual(values = colors, name="Model name") +
  theme_bw() +
  xlab("Year")
dev.off()
# best ones so far





### Prep for TMB models ########################################################################################

require(sdmTMB)

dat$sh <- dat$sh/100

ggplot(dat, aes(lon, lat, size = wmw)) + geom_point()
ggplot(dat, aes(lon, lat, size = sh)) + geom_point()

# Adding UTMs
add_utm_columns(dat, ll_names = c("lon", "lat"))
#Note that this function guesses at an appropriate UTM projection for you and provides a URL to verify. To ensure future compatibility with prediction grids or other data, it is best to hard code the choice using the `utm_crs` argument. We will use `CRS = 3156` here to match our prediction grid:

dat <- add_utm_columns(dat, utm_crs = 32619, ll_names = c("lon", "lat"))

ggplot(dat, aes(X, Y, size = wmw)) + geom_point(shape = 21) + coord_fixed()
ggplot(dat, aes(X, Y, size = sh)) + geom_point(shape = 21) + coord_fixed()

ggplot(dat, aes(X, Y, size = wmw, colour = log(wmw + 1))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~year) + coord_fixed()

dat <- dplyr::arrange(dat, year, tow)
dat$year_factor <- as.factor(dat$year)

# Constructing a mesh

# We start by constructing an SPDE mesh with INLA. This creates some matrices that are used internally when fitting the model. We will use the shortcut function `make_mesh()` and use a cutoff (minimum triangle length of 10 km). A full analysis could explore sensitivity to this decision or other more complex meshes created with INLA.
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 3)
plot(mesh)
mesh$mesh$n # number of vertices or knots



### Fitting spatial models ########################################################################################

# fit_spatial_gam <- sdmTMB(
#   wmw ~ s(sh), # intercept only
#   data = dat,
#   family = lognormal(),#tweedie(link = "log"),
#   mesh = mesh,
#   spatial = "on",
#   control = sdmTMBcontrol(newton_loops = 1)#
#   # silent = FALSE
# )
#
# # check fit
# sanity(fit_spatial_gam)
# fit_spatial_gam
#
# #residuals
# resid <- residuals(fit_spatial_gam)
# qqnorm(resid[is.finite(resid)])
# qqline(resid[is.finite(resid)])
# dat$resid <- residuals(fit_spatial_gam)
# ggplot() +
#   geom_point(data=dat, aes(X, Y, fill=resid))+
#   facet_wrap(~year)
#
# ggplot() + geom_point(data=dat[is.finite(dat$resid),], aes(sh, resid))
# # not very starry for BBn, better for Sab
#
# #calculate residuals manually
# resid <- predict(fit_spatial_gam, dat) # where dat was the dataframe originally used for fitting model
# resid$resid.log <- resid$wmw - resid$est
# resid$resid.resp <- resid$wmw - exp(resid$est)
# ggplot() + geom_point(data=resid, aes(sh, resid.log))
# ggplot() + geom_point(data=resid, aes(sh, resid.resp))
# # not very starry for BBn or Sab
#
# # log instead?
# fit_spatial_log <- sdmTMB(
#   wmw ~ log(sh),
#   data = dat,
#   family = lognormal(),#tweedie(link = "log"),
#   mesh = mesh,
#   spatial = "on"#,
#   #control = sdmTMBcontrol(newton_loops = 1)#
#   # silent = FALSE
# )
#
# # check fit
# sanity(fit_spatial_log)
# fit_spatial_log
#
# #residuals
# resid <- residuals(fit_spatial_log)
# qqnorm(resid[is.finite(resid)])
# qqline(resid[is.finite(resid)])
# dat$resid <- residuals(fit_spatial_log)
# ggplot() +
#   geom_point(data=dat, aes(X, Y, fill=resid))+
#   facet_wrap(~year)
#
# ggplot() + geom_point(data=dat[is.finite(dat$resid),], aes(sh, resid))
# # not very starry for BBn, better for Sab
#
# #calculate residuals manually
# resid <- predict(fit_spatial_log, dat) # where dat was the dataframe originally used for fitting model
# resid$resid.log <- resid$wmw - resid$est
# resid$resid.resp <- resid$wmw - exp(resid$est)
# ggplot() + geom_point(data=resid, aes(sh, resid.log))
# ggplot() + geom_point(data=resid, aes(sh, resid.resp))
# # not very starry for either
#

### Fitting spatiotemporal model ########################################################################################

dat$year <- as.numeric(dat$year)

fit_spatial_ssh_st <- sdmTMB(
  wmw ~ log(sh), #+ 0 + as.factor(year),
  data = dat,
  family = lognormal(),#tweedie(link = "log"),
  mesh = mesh,
  spatial = "on",
  time="year",
  spatiotemporal = "iid",
  control = sdmTMBcontrol(newton_loops = 1)#,
  # silent = FALSE
)

sanity(fit_spatial_ssh_st)
fit_spatial_ssh_st

#residuals
resid <- residuals(fit_spatial_ssh_st)
png(paste0(plotsGo, "/", bank, "/sdmtmb_qqplot.png"), height=4, width=6.5, res=420, units="in")
qqnorm(resid[is.finite(resid)])
qqline(resid[is.finite(resid)])
dev.off()
dat$resid <- residuals(fit_spatial_ssh_st)
ggplot() +
  geom_point(data=dat, aes(X, Y, fill=resid))+
  facet_wrap(~year)

ggplot() + geom_point(data=dat[is.finite(dat$resid),], aes(sh, resid))

#calculate residuals manually
resid <- predict(fit_spatial_ssh_st, dat) # where dat was the dataframe originally used for fitting model
resid$resid.log <- log(resid$wmw) - resid$est
resid$resid.resp <- resid$wmw - exp(resid$est)

png(paste0(plotsGo, "/", bank, "/sdmtmb_resid.png"), height=4, width=6.5, res=420, units="in")
ggplot() + geom_point(data=resid, aes(sh, resid.log)) +theme_bw()
dev.off()
ggplot() + geom_point(data=resid, aes(sh, resid.resp))

# that looks a bit better... still not great though. Use simulations to evaluate later.

# get the survey domain
source("C:/Users/keyserf/Documents/Github/Assessment_fns/Maps/github_spatial_import.R")
survey <- github_spatial_import(subfolder = "survey_boundaries", zipname = "survey_boundaries.zip")
# convert it to the CRS used for the model
survey2 <- survey[survey$ID==bank,] %>%
  st_transform(32619)

# rasterize it
require(stars)
bbn_raster <- as.data.frame((st_rasterize(survey2, st_as_stars(st_bbox(survey2), nx = 25, ny = 25))))
bbn_raster <- bbn_raster[!is.na(bbn_raster$ID),]
bbn_raster$X <- bbn_raster$x/1000
bbn_raster$Y <- bbn_raster$y/1000

# create a dataframe for prediction with year, sh, and locations from survey domain
predict.grid <- expand.grid(year=unique(dat$year), sh = seq(0,2,0.1))
grid.years <- data.frame(X=rep(bbn_raster$X, length(unique(dat$year))),
                         Y=rep(bbn_raster$Y, length(unique(dat$year))),
                         year=rep(unique(dat$year), each=length(bbn_raster$X)))
predict.grid <- full_join(predict.grid, grid.years)

# predict on those variables
p <- predict(fit_spatial_ssh_st, newdata = predict.grid)
p$X<- p$X*1000
p$Y <- p$Y*1000

# plot the MWSH relationship for the "median" location (calculated manually)
simple.grid <- predict.grid[predict.grid$X == sort(predict.grid$X)[round(length(predict.grid$X)/2)] &
                              predict.grid$Y == sort(predict.grid$Y)[round(length(predict.grid$Y)/2)],]
overall <- predict(fit_spatial_ssh_st, newdata=simple.grid)

png(paste0(plotsGo, "/", bank, "/sdmtmb_mwsh.png"), height=4, width=6.5, res=420, units="in")
ggplot() +
  geom_point(data=dat, aes(sh*100, wmw), alpha=0.1) +
  geom_line(data=overall, aes(sh*100, exp(est), group=year), colour="red") +
  theme_bw() +
  xlab("Shell height (mm)") +
  ylab("Meat weight (g)")
dev.off()

# plot as time series (wmw where sh=100mm)
overall <-  predict(fit_spatial_ssh_st, newdata=simple.grid[simple.grid$sh==1,], se_fit=T)
png(paste0(plotsGo, "/", bank, "/sdmtmb_ts.png"), height=4, width=6.5, res=420, units="in")
ggplot() + geom_line(data=overall, aes(year, exp(est), group=1)) +
  geom_point(data=overall, aes(year, exp(est))) +
  geom_errorbar(data=overall, aes(year,
                                  ymin=exp(est-1.96*est_se),
                                  ymax=exp(est+1.96*est_se)), width=0) +
  ylab("Condition factor") +
  xlab("Year") +
  theme_bw()
dev.off()


# crop out the raster cells outside the survey domain
require(raster)
p <- dplyr::select(p, X, Y, year, sh, est, est_non_rf, est_rf, omega_s, epsilon_st)
rast <- unique(p[, c("X", "Y")])
rast <- crop(rasterFromXYZ(rast), survey2)
sp.field <- as(rast, "SpatialPolygonsDataFrame")
proj4string(sp.field) <- st_crs(32619)$proj4string # For SP need that gross full crs code, so this...
# Make it an sf object
spd <- st_as_sf(sp.field, as_points=F, merge=F)
spd <- st_intersection(spd,st_make_valid(survey2))

# join the predictions to the cropped raster polygons
p <- st_as_sf(p, coords=c("X", "Y"), crs=32619, remove=F)
# for the spatiotemporal maps, we only want to visualize wmw where sh = 1 (100mm)
spdp <- st_join(spd, p[p$sh==1,])

# remove NA years
spdp <- spdp[!is.na(spdp$year),]

# bin the response variable (wmw)
spdp$bin <- cut(exp(spdp$est), seq(0,25,2.5))

# ViSuAlIze!
# non random field
ggplot() +
  geom_sf(data=spdp, aes(fill=exp(est_non_rf))) +
  facet_wrap(~year)

# Spatial random field:
ggplot() +
  geom_sf(data=spdp, aes(fill = omega_s), colour=NA)+
  facet_wrap(~year) +
  scale_fill_gradient2()

# Spatial-temporal random field:
ggplot() +
  geom_sf(data=spdp, aes(fill = epsilon_st), colour=NA)+
  facet_wrap(~year) +
  scale_fill_gradient2()

# Overall estimate of wmw in link (log) space:
ggplot()+
  geom_sf(data=spdp, aes(fill = est), colour=NA)+
  facet_wrap(~year)

# Overall estimate of wmw: (with log-distributed colour)
ggplot() +
  geom_sf(data=spdp, aes(fill = exp(est)), colour=NA)+
  facet_wrap(~year) +
 # scale_fill_gradient2()#+
  scale_fill_viridis_c(trans="log10")

# piece de resistance:
forlims <- st_transform(spdp, 4326)
if(bank=="BBn") xdiv <- 0.3; ydiv <- 0.1; wide = 6.5*1.5; high = 5.5*1.5
if(bank=="Sab") xdiv <- 0.6; ydiv <- 0.2; wide = 6.5*1.5; high = 6*1.5
xbreaks <- seq(floor(st_bbox(forlims)$xmin[[1]]/2)*2, ceiling(st_bbox(forlims)$xmax[[1]]/2)*2, xdiv)
ybreaks <- seq(floor(st_bbox(forlims)$ymin[[1]]/2)*2, ceiling(st_bbox(forlims)$ymax[[1]]/2)*2, ydiv)
png(paste0(plotsGo, "/", bank, "/sdmtmb_cf.png"), width=wide, height=high, units="in", res=420)
ggplot() +
  geom_sf(data=spdp, aes(fill = bin), colour=NA)+
  facet_wrap(~year) +
  # scale_fill_gradient2()#+
  scale_fill_viridis_d(name="Condition factor") +
  scale_x_continuous(breaks=xbreaks) +
  scale_y_continuous(breaks=ybreaks) +
  theme_bw()
dev.off()


#### comparing all models:
colors <- c(colors, sdmtmb="red")

png(paste0(plotsGo, "/", bank, "/CF_model_compare2.png"), height=4, width=6.5, units="in", res=420)
ggplot() + geom_line(data=condfac_res$CFyrs, aes(as.numeric(year), CF, colour="orig"),group=1) +
  geom_line(data=CFyrs, aes(year, fit, colour="orig2"), group=1) +
  #geom_line(data=gam_cf$preds, aes(as.numeric(year), inv_link(gam_cf$mod)(fit), colour="gam_cf"), group=1) +
  #geom_line(data=glm_cf$preds, aes(as.numeric(year), inv_link(glm_cf$mod)(fit), colour="glm_cf"), group=1) +
  #geom_line(data=gamm_cf$preds, aes(as.numeric(year), inv_link(gamm_cf$mod)(fit), colour="gamm_cf"), group=1) +
  geom_line(data=glm_mw$preds, aes(as.numeric(year), inv_link(glm_mw$mod)(fit), colour="glm_mw"), group=1) +
  geom_line(data=glm_mw_gamma$preds, aes(as.numeric(year), inv_link(glm_mw_gamma$mod)(fit), colour="glm_mw_gamma"), group=1) +
  #geom_line(data=gam_mw$preds, aes(as.numeric(year), inv_link(gam_mw$mod)(fit), colour="gam_mw"), group=1) +
  #geom_line(data=gam_mw2$preds, aes(as.numeric(year), inv_link(gam_mw2$mod)(fit), colour="gam_mw2"), group=1) +
  #geom_line(data=gam_mw3$preds, aes(as.numeric(year), inv_link(gam_mw3$mod)(fit), colour="gam_mw3"), group=1) +
  #geom_line(data=gam_mw4$preds, aes(as.numeric(year), inv_link(gam_mw4$mod)(fit), colour="gam_mw4"), group=1) +
  geom_line(data=gam_mw5$preds, aes(as.numeric(year), inv_link(gam_mw5$mod)(fit), colour="gam_mw5"), group=1) +
  geom_line(data=glm_mw6$preds, aes(as.numeric(year), inv_link(glm_mw6$mod)(fit), colour="glm_mw6"), group=1) +
  geom_line(data=overall, aes(as.numeric(year), exp(est), colour="sdmtmb"), group=1)+
  scale_color_manual(values = colors, name="Model name") +
  theme_bw() +
  xlab("Year")
dev.off()
