# Port sampling

require(lubridate)
require(dplyr)
require(mgcv)
require(ggplot2)

#BBn
############### Section 3 - BBn and Sable Port sampling analysis -- Section 3###############
direct <- "Y:/Offshore/Assessment/"

# Bring in the port sampling data, which contains meat weights and meat counts for sampled fishing trips
load(paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/Port_sampling_processed_data_2022.rData"))

# survey data
load(paste0(direct,"Data/Survey_data/2022/Survey_summary_output/Survey_all_results.RData"))

# plot destination
plotsGo <- "Y:/Offshore/Assessment/Framework/SFA_25_26_2024/DataInputs/Port_sampling/Figures/"

french <- F # set to T for french
if(french==T) {fr <- "_fr"}
if(french==F) {fr <- ""}

require(rosettafish)
funs <- c("https://raw.githubusercontent.com/freyakeyser/rosetta_shell/main/terms.csv")
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs)
{
  download.file(fun,destfile = basename(fun))
  rosetta_terms <- read.csv(paste0(getwd(),"/",basename(fun)), encoding = "UTF-8")
  file.remove(paste0(getwd(),"/",basename(fun)))
}



# port.sampling$fished[1:10]
# port.sampling$ps.date[1:10]
# port.sampling$ps.fished[1:10]
# summary(port.sampling$ps.fished)

#overall summary bbn and sab
ggplot() +
  geom_point(data = port.sampling[port.sampling$bank %in% c("BBn", "Sab"),], aes(ps.fished, meat_weight)) +
  geom_smooth(data = port.sampling[port.sampling$bank %in% c("BBn", "Sab"),], aes(ps.fished, meat_weight, group=bank)) +
  facet_wrap(~bank, ncol=1) +
  theme_bw()

ggplot() +
  geom_boxplot(data = port.sampling[port.sampling$bank %in% c("BBn", "Sab"),], aes(year, meat_weight, group=year)) +
  facet_wrap(~bank, ncol=1) +
  theme_bw() +
  ylab("Meat weight (g)") +
  xlab("Year")

ggplot() +
  geom_point(data = port.sampling[port.sampling$bank %in% c("BBn", "Sab") & port.sampling$year==2022,], aes(ps.fished, meat_weight)) +
  geom_smooth(data = port.sampling[port.sampling$bank %in% c("BBn", "Sab") & port.sampling$year==2022,], aes(ps.fished, meat_weight, group=bank)) +
  facet_wrap(~bank, ncol=1) +
  scale_x_date()

ggplot() +
  geom_point(data = port.sampling[port.sampling$bank %in% c("BBn", "Sab"),], aes(month, meat_weight)) +
  geom_smooth(data = port.sampling[port.sampling$bank %in% c("BBn", "Sab"),], aes(month, meat_weight, group=bank)) +
  facet_wrap(~bank, ncol=1)

ggplot() +
  geom_boxplot(data = port.sampling[port.sampling$bank %in% c("BBn", "Sab"),], aes(month, meat_weight, group=month)) +
  facet_wrap(~bank, ncol=1)

medians <- port.sampling %>%
  group_by(bank, month, ) %>%
  summarize(monthly_med=median(meat_weight),
            monthly_mean=mean(meat_weight),
            monthly_sd=sd(meat_weight))

medians2 <- port.sampling %>%
  group_by(bank, year) %>%
  summarize(annual_med=median(meat_weight))

# port.sampling <- left_join(port.sampling, medians)
# port.sampling <- left_join(port.sampling, medians2)
# medians <- left_join(medians, medians2)

ggplot() +
  geom_line(data=medians[medians$bank %in% c("BBn", "Sab"),], aes(month, monthly_med)) +
  geom_errorbar(data=medians[medians$bank %in% c("BBn", "Sab"),], aes(month, ymin=monthly_med-(1.96*monthly_sd), ymax=monthly_med+(1.96*monthly_sd)))+
  facet_wrap(~bank, ncol=1)


# For the models to come let's center on the day we believe the ASM's arrive on the scence (May 24, 2010 according to Ginette)
# The day here represents 100 days. FK changed ps.date to ps.fished
port.sampling$mod.day <- as.numeric(port.sampling$ps.fished - as.Date("2010-05-24"))/100


##########################################################################################
########################### Bank by bank now
bank <- "BBn"

bbn.dat <- port.sampling[port.sampling$bank == bank,]

# Now we can move to a gam type of model to allow for non-linearity in the patterns....
bbn.mod <- gamm(mc ~ s(mod.day),data=bbn.dat)
summary(bbn.mod$lme)
summary(bbn.mod$gam)

# Let's take a closer look at these GAM results....
anova(bbn.mod$gam)# Careful with these p values as the smoother comes into play on these...
summary(bbn.mod) #Will be explained later

#Model validation
par(mfrow=c(3,1))
plot(bbn.mod$gam,ylim=c(-10,10))
par(mfrow=c(3,1))
plot(bbn.mod$lme)

#par(mfrow = c(2,2), mar = c(5,5,2,2))
E1 <- resid(bbn.mod$lme, type ="n")
F1 <- fitted(bbn.mod$lme)

par(mfrow = c(2,2))
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab ="Residuals",cex=0.01)
abline(h=0, lty=2)

plot(x = bbn.dat$mod.day,
     y = E1,
     xlab = "time",
     ylab = "Residuals",cex=0.01)
abline(h = 0, lty = 2)

boxplot(E1 ~ fleet, data = bbn.dat)
abline(h = 0, lty = 2)

# get ready to predict
pred.dates <- ymd(paste0(sort(rep(2006:2022,12)),"-",rep(1:12,12),"-",rep(0,144),rep(1,144)))

# center on the day we believe the ASM's arrive on the scene (May 24, 2010 according to Ginette)
pred.dat <- data.frame(date=pred.dates,
                       mod.day=as.numeric(pred.dates - as.Date("2010-05-24"))/100)

mod.pred <- predict(bbn.mod$gam, newdata = pred.dat, se = TRUE)

# Glue them all together
pred.dat$mc <- mod.pred$fit
pred.dat$ub <- mod.pred$fit + 1.96 * mod.pred$se.fit
pred.dat$lb <- mod.pred$fit - 1.96 * mod.pred$se.fit
# Also stick the acutal meat weights on here...
pred.dat$mw <- 500/pred.dat$mc
pred.dat$mw.ub <- 500/pred.dat$ub
pred.dat$mw.lb <- 500/pred.dat$lb

# # Remove all ASM predictions before June 1st 2010 since there were no active vessels before this date (well May 24th, but I'm rounding...)
# pred.dat <- pred.dat[-which(pred.dat$fleet == "ASM" & pred.dat$date < "2010-06-01"),]
# head(pred.dat)

ggplot() +
  geom_point(data = bbn.dat, aes(y = mc, x = date), shape = 16, size = 0.01, alpha=0.05)+
  geom_line(data=pred.dat, aes(date, mc))  +
  geom_ribbon(data = pred.dat, aes(x = date, ymax = ub, ymin = lb), alpha = 0.5)+
  #ylim(5,60) +
  ylab("Meat count") +
  xlab("") +
  scale_x_date(date_breaks="1 year") +
  #theme(text = element_text(size=16)) +
  theme_bw()

#ggsave(paste0(direct,"2018/Framework/Port_sampling/Fleet_port_sampling.png"),width=11,height=8.5)

### now we need a MW-SH relationship to estimate the SH of the landings

# Now what we want is to fit the MW-SH model to these data for each year back to 2006 and for both May and August survey...
# For the most recent data
years <- c(2006:2019,2021:2022)
may.mws <- NULL
a.may <- data.frame(a = rep(NA,length(years)), year = years)


# Now fit a linear model on the wmw v.s. shell height to find the intercept, note this is forced through 0.
# Now we can subset the predicted data from before and figure out what the targeted shell height is...
# pred.dat[month(pred.dat$date) %in% c(5,8),]
for(i in 1:length(years))
{

  may.mws[[as.character(years[i])]] <-  na.omit(mw[[bank]][mw[[bank]]$year == years[i],]) # get BBn and chuck the na's
  may.mws[[as.character(years[i])]]$sh_3 <- (may.mws[[as.character(years[i])]]$sh/100)^3 # cube the SH's
  may.mws[[as.character(years[i])]]$sh <- (may.mws[[as.character(years[i])]]$sh/100) # sh in decimeters
  # Run the model and save the results
  a.may$a[i] <-  lme(fixed = wmw ~ sh_3 -1, data = may.mws[[as.character(years[i])]], random = ~ sh_3 -1 | tow, method="REML")$coefficients$fixed
  # What is the shell height being targeted in may (at the time of the survey)
  # inverse of the above model using wmw value from port sampling data and sh coefficient (ignoring random effect on tow)
  pred.dat$sh[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <-
    100*(pred.dat$mw[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
           a.may$a[a.may$year == years[i]])^(1/3)
  pred.dat$sh.lb[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <-
    100*(pred.dat$mw.lb[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
           a.may$a[a.may$year == years[i]])^(1/3)
  pred.dat$sh.ub[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <-
    100*(pred.dat$mw.ub[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
           a.may$a[a.may$year == years[i]])^(1/3)
}

# Now plot up the results...
lab <- as_labeller(c('5'= "May"))
ggplot(pred.dat[month(pred.dat$date) %in% c(5),],aes(date,sh)) +
  geom_point() +
  geom_line() +
  theme(text = element_text(size=16)) + theme_bw() +
  scale_x_date(date_breaks="2 years") + xlab("") + ylab("Shell height") +
  geom_ribbon(aes(x = date, ymax = sh.ub, ymin = sh.lb), alpha = 0.5)
#ggsave(paste0(direct,"2018/Framework/Port_sampling/sh_targeted.png"),width=11,height=8.5)

# Now based on this I think we can assume that within a fishing season the fleet tends to target roughly the same sized scallop, if so
# we can actual get a feel for how the MW-SH relationship varies over the course of a season by looking at changes in the MW's over the year.
# So for each year let's get the SH targeted as the average of the May and August surveys...
SH.target <- aggregate(sh ~ year(date),pred.dat,mean)
names(SH.target) <- c("year","sh")
# Add an year month combo column to port.sampling...
port.sampling$year_month <- ymd(paste(year(port.sampling$ps.fished),month(port.sampling$ps.fished),15,sep="-"))
mw.by.month <- aggregate(meat_weight ~ year_month + year,port.sampling,mean)
names(mw.by.month) <- c("date","year","mw")

ggplot(mw.by.month,aes(x= date,mw)) + geom_point() + geom_line()

ps.mw.sh <- merge(mw.by.month,SH.target,by =c("year"))

ps.mw.sh$mc <- 500/ps.mw.sh$mw
# Condition here is kinda interesting as it standardizes the meat weight between the fleets, i.e. what is the
# meat weight the different fleets are picking up for a 100 mm scallop.  In theroy I wouldn't expect this to be any different...
# What the questions around condition are:
# 1:  How does condition vary throughout the year, when does condition appear to peak.
# 2:  For a 100 mm scallop, is there any difference between the fleets in terms of the size of the meat they capture.
# I'm not sure if my logic holds up here, but I think it does, hinges on the assumption that the targeted SH doesn't vary
# significantly throughout the year
ps.mw.sh$cond <- ps.mw.sh$mw/(ps.mw.sh$sh/100)^3
ps.mw.sh$month <- as.factor(month(ps.mw.sh$date))

ggplot() +
  geom_point(data = ps.mw.sh, aes(month, cond)) +
  geom_line(data = ps.mw.sh, aes(month, cond, group=1)) +
  facet_wrap(~year)

ggplot() +
  geom_boxplot(data = ps.mw.sh, aes(month, cond, group=month))

ggplot() +
  geom_boxplot(data = ps.mw.sh, aes(month, mw, group=month))

ggplot() +
  geom_boxplot(data = ps.mw.sh, aes(month, mc, group=month))

# OK, so I want to see what propotion of the port sampling data might have a SH of < 95 mm, based on the MW-SH data.
# This is similar to above, but I'm going to have to use the raw data, so this is gonna be slow....

years <- c(2006:2019,2021:2022)
may.mws <- NULL
bbn.dat$sh <- NA
all.may <- data.frame(a = rep(NA,length(years)), year = years)

pred.dates <- ymd(paste0(sort(rep(2006:2022,12)),"-",rep(1:12,12),"-",rep(0,144),rep(1,144)))
# For the models to come let's center on the day we believe the ASM's arrive on the scence (May 24, 2010 according to Ginette)
# This is needed to line up with mod.pred.day... I know this lumps the 31st with the day furthest away from the 31st rather than the next month
# But the idea here is everything in month x is treated the same.
bbn.dat$year_month <- date(paste0(year(bbn.dat$date),"-",bbn.dat$month,'-01'))
bbn.dat$mod.day <- as.numeric(bbn.dat$year_month - as.Date("2010-05-24"))/100
bbn.dat$year <- year(bbn.dat$date)

for(i in 1:length(years))
{
  may.mws[[as.character(years[i])]] <-  na.omit(mw[[bank]][mw[[bank]]$year == years[i],]) # get BBn and chuck the na's
  may.mws[[as.character(years[i])]]$sh_3 <- (may.mws[[as.character(years[i])]]$sh/100)^3 # cube the SH's
  may.mws[[as.character(years[i])]]$sh <- (may.mws[[as.character(years[i])]]$sh/100) # sh in decimeters
  # Run the model and save the results
  a.may$a[i] <-  lme(fixed = wmw ~ sh_3 -1, data = may.mws[[as.character(years[i])]], random = ~ sh_3 -1 | tow, method="REML")$coefficients$fixed
  # What is the shell height being targeted in may
  bbn.dat$sh[bbn.dat$year == years[i]] <- 100*(bbn.dat$meat_weight[bbn.dat$year == years[i]]/ a.may$a[i])^(1/3)
}

sh.plt <- ggplot(bbn.dat,aes(date,sh)) +
  geom_point(alpha = 0.05,size=0.2) +
  geom_smooth(method = 'gam')+
  geom_hline(yintercept = 95, linetype = 'dashed',linewidth=1.5,color="grey") +
  theme(text = element_text(size=16)) +
  theme_bw() +
  scale_x_date(date_breaks="2 years") +
  scale_y_continuous(breaks = seq(0,200,by=5)) +
  xlab("") +
  ylab("Shell height")
sh.plt


# HERE IS THE MAIN ANALYSIS FOR the BBn Res doc, shows proportion below certain sizes
n.bl.100 <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 100) %>% dplyr::summarise(n = length(sh))
n.all.bbn <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::summarise(tot = length(sh))

prop.bl.100 <- left_join(n.bl.100,n.all.bbn,'year')
prop.bl.100$prop <- prop.bl.100$n / prop.bl.100$tot

ggplot(data = prop.bl.100,aes(x=year,y=prop)) +
  geom_point() + ylab("Proportion of meats below 100 mm") +
  #theme(text = element_text(size=22)) +
  scale_x_continuous(breaks = seq(2000,2030,by=2))


n.bl.95 <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 95) %>% dplyr::summarise(n = length(sh))
prop.bl.95 <- left_join(n.bl.95,n.all.bbn,'year')
prop.bl.95$prop <- prop.bl.95$n / prop.bl.95$tot

ggplot(data = prop.bl.95,aes(x=year,y=prop)) +
  geom_point() +
  ylab("Proportion of meats below 95 mm") +
  # theme(text = element_text(size=22)) +
  scale_x_continuous(breaks = seq(2000,2030,by=2))

# Based on Ageing info, we think 75-90 mm is the best option, that would be ≈ 3-4 year olds.
n.bl.90 <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 90) %>% dplyr::summarise(n = length(sh))
prop.bl.90 <- left_join(n.bl.90,n.all.bbn,'year')
prop.bl.90$prop <- prop.bl.90$n / prop.bl.90$tot

# This is basically under 1%, which I think will be fine.
ggplot(data = prop.bl.90,aes(x=year,y=prop)) +
  geom_point() +
  ylab("Proportion of meats below 90 mm") +
  #theme(text = element_text(size=22)) +
  scale_x_continuous(breaks = seq(2000,2030,by=2))


n.bl.85 <- bbn.dat %>% dplyr::group_by(year,.drop=F) %>% dplyr::filter(sh < 85) %>% dplyr::summarise(n = length(sh))
prop.bl.85 <- left_join(n.bl.85,n.all.bbn,'year')
prop.bl.85$prop <- prop.bl.85$n / prop.bl.85$tot

# So there is basically 0 catch coming from < 85 mm, well under 1%
ggplot(data = prop.bl.85,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 90 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))


prop.bl.100$lim <- 100
prop.bl.95$lim <- 95
prop.bl.90$lim <- 90
prop.bl.85$lim <- 85

overall <- full_join(full_join(full_join(prop.bl.100, prop.bl.95), prop.bl.90), prop.bl.85)

png(paste0(plotsGo, bank, "/proportion_boxplot_", fr, ".png"), width=6.5, height=4, units = "in", res=420)
ggplot() +
  geom_boxplot(data=overall, aes(lim, prop, group=lim))+
  theme_bw() +
  ylab("Proportion of meats below fully recruited size") +
  xlab("Minimum fully recruited size (mm)")
dev.off()


# # Just look at spring-summer, going with April-August as in theory condition should be fairly stable during this period.
# # Based on Ageing info, we think 75-90 mm is the best option, that would be ≈ 3-4 year olds.
# sel.months <- 4:8
# n.bl.90.ss <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 90 & month %in% sel.months) %>% dplyr::summarise(n90 = length(sh))
# n.all.bbn.ss <-  bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(month %in% sel.months) %>% dplyr::summarise(tot = length(sh))
# prop.bl.90.ss <- left_join(n.bl.90.ss,n.all.bbn.ss,'year')
# prop.bl.90.ss$prop <- prop.bl.90.ss$n90 / prop.bl.90.ss$tot
#
# # This is basically under 1%, which I think will be fine.
# windows(11,11)
# ggplot(data = prop.bl.90.ss,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 90 mm") +
#   theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))


# Based on the ageing work, if we go with 90 mm as our end of recruitment bin, the recruits should be 75-90 mm in size.
