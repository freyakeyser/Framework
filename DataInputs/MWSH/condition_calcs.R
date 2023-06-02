require(dplyr)
# compare cf approaches
for(b in  c("Mid", "Ban", "Sab", "BBn", "BBs")){
  original <- cf.data[[b]]$CFyrs
  original$type <- "original"
  original$year <- as.numeric(original$year)
  predicted <- data.frame(year=survey.obj[[b]]$model.dat$year, CF=survey.obj[[b]]$model.dat$CF, type="predicted")
  comb <- full_join(original, predicted) %>%
    dplyr::select(year, type, CF) %>%
    pivot_wider(id_cols = year, values_from = CF, names_from=type)

  print(ggplot() + geom_text(data=comb, aes(original, predicted, label=year)) +
          geom_abline(slope = 1, intercept=0)+
          ggtitle(b))

  comb2 <- comb %>%
    pivot_longer(!year)

  print(ggplot() + geom_point(data=comb2, aes(year, value, shape=name, colour=name)) +
          geom_line(data=comb2, aes(year, value, colour=name))+
          ggtitle(b) +
          geom_point(data=test$CFyrs, aes(as.numeric(year), CF)))

  # what happens if you convert to biomass

}


source("C:/Users/keyserf/Documents/Github/Assessment_fns/Survey_and_OSAC/condFac.r")
bnk <- "Mid"
#examples:
# CF.current$BBn[CF.current$BBn$tow==201,]
# cf.data$BBn$CF.data[cf.data$BBn$CF.data$year==2022 & cf.data$BBn$CF.data$tow==201,]

# get data ready for modelling
mw.dat.mod <- mw.dat.all$Mid[complete.cases(mw.dat.all$Mid),]

yrs<-sort(unique(mw.dat.all[["Mid"]]$year))
pred.dat <- bank.dat[["Mid"]][bank.dat[["Mid"]]$year %in% as.numeric(yrs),]

test <- condFac(mw.dat.all[["Mid"]],pred.dat,model.type='glm',dirct=direct_fns)

test$CFyrs
### SO in survey data, exclude the years w/o sampling from bank.dat for condition calculation. This trickles through
### within surv.dat object, such that those years will never get biomass or abundance estimates
### should retain abundance estimates, just not biomass. Could potentially just hide CF and Biomass estimates for these years
### but only if it has no impact on the other years.

require(ggplot2)
require(patchwork)

load("C:/Users/keyserf/Documents/temp_data/testing_results_framework_75-90.Rdata")

before <- survey.obj
beforecf <- cf.data

load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_mid.RData")
load("Y:/Offshore/Assessment/Data/Survey_data/2022/Survey_summary_output/testing_results_framework2.RData")

after <- survey.obj
aftercf <- cf.data

for (i in names(after)[!names(after) %in% c("GB", "GBa", "GBb")]){
  print(i)

  a<-ggplot() + geom_point(data=before[[i]]$model.dat, aes(year, N)) +
    geom_point(data=after[[i]]$model.dat, aes(year, N), colour="red") +
    ggtitle(i)

  b<-ggplot() + geom_boxplot(data=beforecf[[i]]$pred.dat, aes(year, CF, group=year)) +
    geom_boxplot(data=aftercf[[i]]$pred.dat, aes(year, CF, group=year), colour="red")+
    ggtitle(i)

  c<- ggplot() + geom_point(data=before[[i]]$model.dat, aes(year, I)) +
    geom_point(data=after[[i]]$model.dat, aes(year, I), colour="red")+
    ggtitle(i)

  print(a + b + c)

}

for (i in names(after)[!names(after) %in% c("GB", "GBa", "GBb")]){
  c[[i]] <- ggplot() + geom_point(data=before[[i]]$model.dat, aes(year, I)) +
    geom_point(data=after[[i]]$model.dat, aes(year, I), colour="red")+
    ggtitle(i)
}


banks <- names(after)[!names(after) %in% c("GB", "GBa", "GBb")]
c[[banks[1]]] +
  c[[banks[2]]] +
  c[[banks[3]]] +
  c[[banks[4]]] +
  c[[banks[5]]] +
  c[[banks[6]]]









# trying a gam directly on the raw data
require(mgcv)
mw.mod <- gam(log(wmw)~s(sh) + s(lon,lat)+s(depth)+as.factor(year),data=mw.dat.mod)
summary(mw.mod)
mw.dat.mod$resid.gam <- residuals(mw.mod)

ggplot(mw.dat.mod, aes(x=sh,y=resid.gam)) +geom_point() + geom_smooth(method = 'gam')
dev.off()
qqnorm(mw.dat.mod$resid.gam)
qqline(mw.dat.mod$resid.gam)

# setting up prediction data frame
pred.dat <- expand.grid(sh=50:200, depth=mean(mw.dat.mod$depth), lat=mean(mw.dat.mod$lat), lon=mean(mw.dat.mod$lon), year=unique(mw.dat.mod$year))
# predict using GAM
pred.dat$fit <- exp(predict(object=mw.mod, newdata=pred.dat))
ggplot() +geom_point(data=mw.dat.mod, aes(x=sh,y=wmw)) +
  geom_line(data=pred.dat, aes(x=sh, y=fit, group=year), colour="red")

# the current approach uses MWSH relationship with tow ID (year.tow) as random effect
CF.fit<-gam(CF~s(lon,lat)+s(depth)+as.factor(year),data=cf.data$BBn$CF.data)
summary(CF.fit)
cf.data$BBn$CF.data$resid.gam <- residuals(CF.fit)
ggplot(cf.data$BBn$CF.data, aes(x=depth,y=resid.gam)) +geom_point() + geom_smooth(method = 'gam')
dev.off()
qqnorm(cf.data$BBn$CF.data$resid.gam)
qqline(cf.data$BBn$CF.data$resid.gam)
pred.dat$fit2 <- predict(object=CF.fit, newdata=pred.dat)

mw.dat.mod$sh.mod <- (mw.dat.mod$sh/100)^3
mw.dat.mod$ID <- paste0(mw.dat.mod$year, ".", mw.dat.mod$tow)
wt.lme <- lme(fixed = wmw ~ sh.mod -1, data = mw.dat.mod, random = ~ sh.mod -1 | year/tow, method="REML",control=lmeControl(opt='optim'))
test <- lme4::lmer(data=mw.dat.mod, wmw ~ sh.mod + (1+sh.mod|year) + (1|tow))
fixef(test)
ranef(test)

pred.dat$sh.mod <- (pred.dat$sh/100)^3
A <- fixef(test)[2]
B <- fixef(test)[1]
raneff <- data.frame(year=as.factor(rownames(ranef(test)$year)), int = ranef(test)$year[1], slope=ranef(test)$year[2])
names(raneff) <- c("year", "int", "slope")
pred.dat <- left_join(pred.dat, raneff)
pred.dat$fit3 <- pred.dat$int + pred.dat$sh.mod*pred.dat$slope

cols <- c("fit"="blue", "fit2"="red", "fit3" = "green")
pred.dat$year <- as.numeric(as.character(pred.dat$year))

png("cfts_methods.png", height=6, width=8, units="in", res=420)
ggplot() +geom_point(data=cf.data$BBn$CF.data, aes(x=as.numeric(year), y=CF), alpha=0.2) +
  geom_line(data=pred.dat, aes(x=year, y=fit2, group=1, colour="fit2")) +
  geom_line(data=pred.dat[pred.dat$sh==100,], aes(x=year, y=fit, group=1, colour="fit")) +
  scale_color_manual(values = cols[c("fit", "fit2")], labels=c(formula(mw.mod), formula(CF.fit)), name="Model\n(wmw predicted at 100mm)") +
  theme_bw()+
  theme(legend.position = "top", legend.direction = "vertical")
dev.off()

png("mwsh_methods.png", height=8, width=12, units="in", res=420)
ggplot() + #geom_point(data=mw.dat.mod, aes(x=sh,y=wmw), alpha=0.2, size=0.5) +
  #geom_line(data=pred.dat, aes(x=sh, y=fit, group=year, colour="fit")) +
  #geom_line(data=pred.dat, aes(x=sh, y=fit2, group=year, colour="fit2")) +
  geom_line(data=pred.dat, aes(x=sh, y=fit3, colour="fit3", group=year))+
  #facet_wrap(~year) +
  geom_vline(xintercept=100) +
  scale_color_manual(values = cols, labels=c(formula(mw.mod), formula(CF.fit), "wmw ~ A*sh^3"), name="Model") +
  theme_bw()+
  theme(legend.position = "top", legend.direction = "vertical")
dev.off()


without <- survey.dat(shf = surv.Rand$BBn, years=unique(surv.Rand$BBn$year), RS=85, CS=95, bk="BBn", areas=survey.info[survey.info$label=="BBn"& survey.info$startyear==2021,c("Strata_ID", "towable_area", "startyear")], mw.par="CF", user.bins = bin)

without$model.dat$I
with$model.dat$I




#### from dave via teams 2022-12-22
library(lme4)
library(nlme)
library(arm)
sh <- rep(50:150,100)
a = 0.00001
b = 3
years <- sample(1996:2020,10100,replace=T)
wmw = a*sh^3 * rlnorm(10100,0,0.3)
plot(wmw~sh)
mwsh <- data.frame(wmw=wmw,sh=sh,year =years )
ggplot(mwsh,aes(x=sh,y=wmw,color = year,group=year)) + geom_point() + geom_smooth(method = 'lm') + scale_x_log10() + scale_y_log10()

tst <- lmer(log(wmw) ~ log(sh) + (1+ log(sh)|year),data = mwsh)
fixef(tst)
ranef(tst)
# Now to do it the wacko complicated offshore way...
# First just grab a year...
mwsh2 <- mwsh %>% dplyr::filter(years == 2020)
mwsh2 <- mwsh2[1:300,]
mwsh2$tow <- rep(1:30,10)
mwsh2$sh3 <- mwsh2$sh^3
wt.lme <- lme(fixed = wmw ~ sh3 -1, data = mwsh2, random = ~ sh3 -1 | tow, method="REML",control=lmeControl(opt='optim'))
ggplot(mwsh2,aes(x=sh3,y=wmw)) + geom_point()