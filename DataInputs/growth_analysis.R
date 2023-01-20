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

# hydration data
load("C:/Users/keyserf/Documents/temp_data/testing_results_framework.RData")

dat <- mw.dat.all$BBn
dat$sh <- dat$sh/100
require(sdmTMB)


head(dat)

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

# # ggplot alternative:
# ggplot() + inlabru::gg(mesh$mesh) + coord_fixed() +
#   geom_point(aes(X, Y), data = dat, alpha = 0.2, size = 0.5)

# Fitting a spatial model

# The most basic model we could fit would be a model with only spatial random fields. This is similar to what we were fitting yesterday. Using `silent = FALSE` lets us see what is happening, but it's awkward when running code in Rmd/Quarto chunks (with the default RStudio setting to 'Show output inline for all R Markdown document') so we will comment it out in most cases here. But it is a good idea to use it if models are running slowly or not converging to monitor progress.

# fit_spatial <- sdmTMB(
#   wmw ~ 1, # intercept only
#   data = dat,
#   family = lognormal(),#tweedie(link = "log"),
#   mesh = mesh,
#   spatial = "on"#,
#   # silent = FALSE
# )
#
# fit_spatial_sh <- sdmTMB(
#   wmw ~ sh,
#   data = dat,
#   family = lognormal(),#tweedie(link = "log"),
#   mesh = mesh,
#   spatial = "on"#,
#   # silent = FALSE
# )
#
# fit_spatial_sh_s <- sdmTMB(
#   wmw ~ s(sh),
#   data = dat,
#   family = lognormal(),#tweedie(link = "log"),
#   mesh = mesh,
#   spatial = "on",
#   control = sdmTMBcontrol(newton_loops = 1)#,
#   # silent = FALSE
# )

dat$year <- as.numeric(dat$year)
#
# fit_spatial_sh_st <- sdmTMB(
#   wmw ~ sh,
#   data = dat,
#   family = lognormal(),#tweedie(link = "log"),
#   mesh = mesh,
#   spatial = "on",
#   time="year",
#   spatiotemporal = "iid",
#   control = sdmTMBcontrol(newton_loops = 1)#,
#   # silent = FALSE
# )

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
qqnorm(resid[is.finite(resid)])
qqline(resid[is.finite(resid)])
dat$resid <- residuals(fit_spatial_ssh_st)
ggplot() +
  geom_point(data=dat, aes(X, Y, fill=resid))+
  facet_wrap(~year)

ggplot() + geom_point(data=dat[is.finite(dat$resid),], aes(sh, resid))

#calculate residuals manually
resid <- predict(fit_spatial_ssh_st, dat) # where dat was the dataframe originally used for fitting model
resid$resid.log <- log(resid$wmw) - resid$est
resid$resid.resp <- resid$wmw - exp(resid$est)

ggplot() + geom_point(data=resid, aes(sh, resid.log))
ggplot() + geom_point(data=resid, aes(sh, resid.resp))


# get the survey domain
source("C:/Users/keyserf/Documents/Github/Assessment_fns/Maps/github_spatial_import.R")
survey <- github_spatial_import(subfolder = "survey_boundaries", zipname = "survey_boundaries.zip")
# convert it to the CRS used for the model
survey2 <- survey[survey$ID=="BBn",] %>%
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
ggplot() +
  geom_point(data=dat, aes(sh, wmw), alpha=0.1) +
  geom_line(data=overall, aes(sh, exp(est), group=year), colour="red")

# plot as time series (wmw where sh=100mm)
overall <-  predict(fit_spatial_ssh_st, newdata=simple.grid[simple.grid$sh==1,], se_fit=T)
ggplot() + geom_line(data=overall, aes(year, exp(est), group=1)) +
  geom_point(data=overall, aes(year, exp(est))) +
  geom_errorbar(data=overall, aes(year,
                                  ymin=exp(est-1.96*est_se),
                                  ymax=exp(est+1.96*est_se)), width=0)


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
ggplot() +
  geom_sf(data=spdp, aes(fill = bin), colour=NA)+
  facet_wrap(~year) +
  # scale_fill_gradient2()#+
  scale_fill_viridis_d()






# Cross validation

# optional for parallel processing (just FYI)
# library(future)
# plan(multisession)

# set.seed(1) # for consistent cross validation folds
# cv_intercept <- sdmTMB_cv(
#   wmw ~ 1, # intercept only
#   data = dat,
#   family = lognormal(),#tweedie(link = "log"),
#   mesh = mesh,
#   k_folds = 8
# )
#
# cv_sh <- sdmTMB_cv(
#   wmw ~ sh, # intercept only
#   data = dat,
#   family = lognormal(),#tweedie(link = "log"),
#   mesh = mesh,
#   k_folds = 5
# )
#
# set.seed(1) # for consistent cross validation folds
# cv_sh_s <- sdmTMB_cv(
#   wmw ~ s(sh),
#   data = dat,
#   family = lognormal(), #tweedie(link = "log"),
#   mesh = mesh,
#   k_folds = 3
# )
#
# cv_sh_l <- sdmTMB_cv(
#   wmw ~ log(sh),
#   data = dat,
#   family = lognormal(), #tweedie(link = "log"),
#   mesh = mesh,
#   k_folds = 3
# )

# #These are the expected log predictive densities (ELPD) for each fold. Larger (more positive) values are better. This gives us an idea of the variance of these values across folds.
# cv_intercept$fold_elpd
# cv_sh$fold_elpd
#
# # check each fold between the two models:
# cv_sh$fold_elpd - cv_intercept$fold_elpd
#
# #These are the expected log predictive densities for all left out data combined. Larger (more positive) values are better.
# cv_intercept$elpd
# cv_sh$elpd

