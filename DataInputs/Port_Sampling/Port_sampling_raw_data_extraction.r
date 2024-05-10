# This script grabs all the Port sampling data for each year, converts it to long format data, merges it with the appropriate fishery meta-data giving
# us the beautiful end result of having all the Port Sampling data tied to fishery information and being ready to analyze...
# What is found in the port sampling data is the individual meat weights of scallop measured.  The port sampling "sample" is different
# from the sample used to calculated the meat count, which if fine.  What we have here is simply a sample of the catch from each day
# and there is no reason to assume it isn't a represenative sample so we can track the size of scallop meat being landed for the last decade
# While not perfect information it certainly is useful information!!

# Frist grab the directory to work from
direct <- "Y:/Offshore/Assessment/"
direct_fns <- "C:/Users/keyserf/Documents/Github/"
# Grab the years of port sampling data.  Note that Ginette says the Port Sampling data goes all the way back to 1996
# I haven't a clue where that data is, but would certainly be interesting to compare the last 20 years of data!!
years <- 2006:2022
options(stringsAsFactors = F) # Don't make anything a Factor as it screws up random crap.

# Load in librarys we may need...
library(ggplot2)
require(readxl)
library(reshape2)
library(dplyr)
library(plyr)
library(lubridate)
library(PBSmapping)
library(ggfortify)
library(mgcv)
# We will need our fishery log functionto get the log data....
source(paste(direct_fns,"Assessment_fns/Fishery/logs_and_fishery_data.r",sep="")) #logs_and_fish is function call
source(paste(direct_fns,"Assessment_fns/Survey_and_OSAC/shwt.lme.r",sep="")) # The Meat weigth SH model we use for offshore...
#source(paste("D:/Github/Assessment_fns/Fishery/logs_and_fishery_data.r",sep="")) #logs_and_fish is function call

# And we need the Offshore Scallop fleet name file so we can link the fishery data to these names
off.fleet <- read.csv(paste0(direct,"Data/Offshore_fleet.csv"))
logs_and_fish(loc = "offshore", year = 2006:2022, get.local = T, get.marfis = F, direct=direct)

################  Section 1, processing the port sampling meat weight information  --  Section 1 ############################################################
################  Section 1, processing the port sampling meat weight information  --  Section 1 ############################################################
# Here we pull in all the port sampling data (right now 2006-2017) and tidy it up for later analysis, once happy with the results
# you can skip this Section
require(tidyverse)

# Get ready for the loop...
index <- 0
dat <- NULL
# Run this for all years we have data, this takes about 10 minutes...
for(i in 1:length(years))
{
  # # This identifies all the files that for each year of Port Sampling, this should include the vast majority of trips each year, probably would be useful
  # # to have someone check if these are exhaustive lists or not...
  # if(years[i]<2006) {
  #   files <- list.files(path = paste0(direct,"Data/Archive/PortSampling/", years[i], "/"),pattern = "\\.txt", recursive = T, ignore.case = T)
  # }
  # if(years[i]==2006) {
  #   files <- list.files(path = paste0(direct,"Data/PortSampling/PS", years[i], "/origportsamplesfiles"),pattern = "\\.xls", recursive = T, ignore.case = T)
  # }
  # if(years[i]==2007) {
  #   files <- list.files(path = paste0(direct,"Data/PortSampling/PS", years[i], "/OriginalPortSamples2007"),pattern = "\\.xls", recursive = T, ignore.case = T)
  # }
  # if(years[i]==2008) {
  #   files <- list.files(path = paste0(direct,"Data/PortSampling/PS", years[i], "/OrigianlPortSamples2008"),pattern = "\\.xls", recursive = T, ignore.case = T)
  # }
  # if(years[i]==2009) {
  #   files <- list.files(path = paste0(direct,"Data/PortSampling/PS", years[i], "/originals"),pattern = "\\.xls", recursive = T, ignore.case = T)
  # }
  # if(years[i]==2010) {
  #   files <- list.files(path = paste0(direct,"Data/PortSampling/PS", years[i], "/Originals2010"),pattern = "\\.xls", recursive = T, ignore.case = T)
  # }
  # if(years[i]>2010 & years[i]<2018) {
  #   files <- list.files(path = paste0(direct,"Data/PortSampling/PS", years[i], "/JoansOriginals"),pattern = "\\.xls", recursive = T, ignore.case = T)
  # }
  # if(years[i]>2017) {
  #   files <- list.files(path = paste0(direct,"Data/PortSampling/PS", years[i], "/JoansOriginals"),pattern = "\\.csv", recursive = T, ignore.case = T)
  # }
  # if(years[i] < 2018) files2 <- list.files(path = paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/", years[i], "/"),pattern = "\\.xls", ignore.case = T)
  # if(years[i] > 2017) files2 <- list.files(path = paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/", years[i], "/"),pattern = "\\.csv", ignore.case = T)
  # print(years[i])
  # #
  #
  # # if(years[i] <2006){
  # #
  # # }
  # if(years[i] < 2018 & years[i]>2005) {
  #   files <- gsub(x = tolower(files), pattern=".xlsx", replacement="", fixed=T)
  #   files <- gsub(x = tolower(files), pattern=".xls", replacement="", fixed=T)
  #   files2 <- gsub(x = tolower(files2), pattern=".xlsx", replacement="", fixed=T)
  #   files2 <- gsub(x = tolower(files2), pattern=".xls", replacement="", fixed=T)
  # }
  # if(years[i] > 2017) {
  #   files <- gsub(x = tolower(files), pattern=".csv", replacement="", fixed=T)
  #   files2 <- gsub(x = tolower(files2), pattern=".csv", replacement="", fixed=T)
  # }
  #
  # if(any(nchar(files)>8)) files <- str_sub(string = files,-8)
  #
  # # names that are in files2 but not in files
  # files2[which(!files2 %in% files)]
  #
  # # names that are in files but not in files2
  # files[which(!files %in% files2)]

# 2013, COME0613 # this trip doesn't exist?
  if(years[i] < 2006) files <- list.files(path = paste0(direct,"Data/Archive/PortSampling/", years[i], "/"),pattern = "\\.txt", ignore.case = T)
  if(years[i] < 2018 & years[i]>2005) files <- list.files(path = paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/", years[i], "/"),pattern = "\\.xls", ignore.case = T)
  if(years[i] > 2017) files <- list.files(path = paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/", years[i], "/"),pattern = "\\.csv", ignore.case = T)
  num.files <- length(files)

  for(j in 1:num.files)
  {
    print(index)
    index <- index + 1
    #This will pull the data from the Port Sampling file.
    if(years[i] < 2006) dat[[index]] <- read.csv(paste0(direct,"Data/Archive/PortSampling/", years[i], "/", files[j]), sep = "\t", header = F)
    if(years[i] < 2018 & years[i]>2005) dat[[index]] <- read_excel(path=paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/",years[i],"/", files[j]), sheet = 1)
    if(years[i] > 2017) dat[[index]] <- read.csv(paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/",years[i],"/", files[j]))
    # # Make all the variable names lower case as the rest of this works...
    names(dat[[index]]) <- tolower(names(dat[[index]]))
    # Remove the blank columns
    remove <- names(dat[[index]])[which(names(dat[[index]]) == "space")]
    dat[[index]] <- dat[[index]][,!names(dat[[index]]) %in% remove]

    dat[[index]] <- reshape2::melt(dat[[index]],id.vars = c("date","boat","port","id","fished"))
    # And now reorder the data by ID so the samples all stay together and the 0's come at the end so we can chuck those..
    dat[[index]] <- dat[[index]] %>% arrange(desc(value)) # First order the values from biggest to smallest
    dat[[index]] <- dat[[index]] %>% arrange(id) # Now order the ID's from smallest to largest

    # Now we can quickly add a sample ID to these in case we want it later...
    samp.ids <- unique(dat[[index]]$id)
    n.samp.ids <- length(samp.ids) # The first 1-2 digits of the ID is the day of the trip (i.e. 1 = first day, 2 = second day)
                                    # The last 2 digits decribe the location of the sample in the tub (i.e. middle/top/bottom/front/back etc.)
    for(r in 1:n.samp.ids)
    {
      dat[[index]]$sample_id[dat[[index]]$id == samp.ids[r]] <- 1:nrow(dat[[index]][dat[[index]]$id == samp.ids[r],])
    } # end for(r in 1:n.samp.ids)
    # And now we can chuck all the 0's
    dat[[index]] <- dat[[index]][which(dat[[index]]$value > 0),]

  } # end for(j in 1:num.files)
} # end the for(i in 1:length(years))

# Pull the data out from the list...
port.dat <- do.call("rbind",dat)
# Rename the value column to "meat_weight".
names(port.dat)[which(names(port.dat) == "value")] <- "meat_weight"
port.dat$meat_weight <- port.dat$meat_weight/100
# make all the vessel names All Caps to match the Offshore_Fleet csv rather than a mix of upper and lower case...
port.dat$boat <- toupper(port.dat$boat)

# Now make the dates actual dates...
port.dat$ps.date <- ymd(port.dat$date)
port.dat$date <- ymd(port.dat$date) # This will get overwritten later with the fishery log date...
port.dat$ps.fished <- ymd(port.dat$fished)
#check dates that failed to parse.

head(port.dat)

### Now we can work on merging the data to fishery log information....
port.dat$year <- year(port.dat$ps.date)
port.dat$month <- month(port.dat$ps.fished)
port.dat$day <- day(port.dat$ps.fished)
# There are some types in the fished entries for the year, so this is simple fix for those key punch mistakes...
port.dat$ps.fished <- ymd(paste(port.dat$year,port.dat$month,port.dat$day))

save.image(paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/Port_sampling_raw_data_2022.RData"))
#load(paste0("C:/Users/keyserf/Documents/temp_data/Data/PortSampling/PS_data_reorg_4_analysis/Port_sampling_raw_data_2022.RData"))

############  Section 2 - Add in teh fishery data to the port sampling...
# Now bring in the latest fishery data
logs_and_fish(loc="offshore",year = years,get.marfis = F, direct = direct, direct_fns=direct_fns)
# If you get any NA's related warnings it may be something is being treated as a Factor in one of the two files.
# This should combine without any warnings so don't ignore warnings here.
fish.dat<-merge(new.log.dat,old.log.dat,all=T)
fish.dat$ID<-1:nrow(fish.dat)

# The vessel ID's for the port sampling data...
vessel_ids <- as.character(unique(port.dat$boat))
# Loop through this for each vessel.
num.vessels <- length(vessel_ids)
# I want to add these columns to the port sampling data from the fishery data....
fish.df.to.add <- data.frame(date = NA, lon= NA,lat=NA,fleet = NA, ves=NA,vrnum=NA, vesid = NA,bank=NA,pro.repwt = NA,
                             h = NA, m=NA, hm = NA, kg.h = NA, kg.hm = NA, depth = NA, tripnum=NA)
port.dat <- data.frame(port.dat,fish.df.to.add)

head(port.dat)

# look at Port Sampling Report xxxxxx .doc files to get totals on bank for 2006

# Run this for each vessel and each day...

final.dat <- NULL
ves.port.dat <- NULL

for(k in 1:num.vessels)
{
  print(paste0("k=",k))
  # Get the vessel of interest
  vessel <- na.omit(off.fleet[off.fleet$ID_alt_Port_Sampling == vessel_ids[k], c("Pre_2008_ID","VMS_old_ID","ID")])
  # Get all the fishery data for that vessel, checking all types of identifiers
  ves.fish.dat <- fish.dat[fish.dat$vesid %in% unlist(vessel) | fish.dat$ves %in% unlist(vessel) | fish.dat$vrnum %in% unlist(vessel),]
  # Get the port sampling data for this vessel
  ves.port.dat <- port.dat[port.dat$boat == vessel_ids[k],]
  # Get the dates we port sampling data for this vessel
  dates.ps <- unique(ves.port.dat$ps.fished)
  # Now I've tried this step with using a merge, but the data isn't tidy enough for it, so I'm going for a for loop and stuffing them together
  # the slow way... m =1
  #merged.dat <- NULL
  for(m in 1:length(dates.ps))
  {
    print(paste0("m=",m))
    # Now get the fishery data for this vessel and these dates.
    fishery.data.ps <- ves.fish.dat[ves.fish.dat$date %in% dates.ps[m],names(fish.df.to.add)]

    # Now we need to take care of fishery data for when we have multiple entries for a vessel day combination
    if(nrow(fishery.data.ps) > 1 && (length(unique(fishery.data.ps$bank))==1)) # 2-4 watches
    {
      # Make a temporary object with the first row of fishery data
      tmp <- fishery.data.ps[1,]
      # Now add up catch and effort for each watch
      tmp[,c("pro.repwt","h","m","hm")] <- colSums(fishery.data.ps[,c("pro.repwt","h","m","hm")],na.rm=T)
      # Find the centroid of the lat's and lon's..  This may/may not be ideal but seems like a legit way to lump the data into a general area fished...
      tmp[,c("lat","lon")] <- colMeans(fishery.data.ps[,c("lat","lon")],na.rm=T)
      # Calculate the effort by day
      tmp[,c("kg.h","kg.hm")] <- c(tmp$pro.repwt/tmp$h, tmp$pro.repwt/tmp$hm)
      # And pop that back into the fishery data object which is now just a single line of data... hopefully
      fishery.data.ps <- tmp
    } # end if(nrow(fishery.data.ps) > 1) to deal with multiple observations each day.

    # if the trip fished multiple banks in one day, we drop the PS data because we can't reasonably assign to one bank or the other.
    if(nrow(fishery.data.ps) > 1 && (length(unique(fishery.data.ps$bank))>1)) # 2-4 watches
    {
      # Make a temporary object with the first row of fishery data
      tmp <- fishery.data.ps[1,]
      # Now add up catch and effort for each watch
      tmp[,c("pro.repwt","h","m","hm")] <- NA
      # Find the centroid of the lat's and lon's..  This may/may not be ideal but seems like a legit way to lump the data into a general area fished...
      tmp[,c("lat","lon")] <- NA
      # Calculate the effort by day
      tmp[,c("kg.h","kg.hm")] <- NA
      tmp$bank <- NA
      # And pop that back into the fishery data object which is now just a single line of data... hopefully
      fishery.data.ps <- tmp
    }

    # We then pump this data back into the port sampling data
    if(nrow(fishery.data.ps) > 0) ves.port.dat[ves.port.dat$ps.fished == dates.ps[m],names(fish.df.to.add)] <- fishery.data.ps
  } # End for(m in 1:length(dates.ps)) to sweep through each date...
  final.dat[[vessel_ids[k]]] <- ves.port.dat

} # end for(k in 1:length(num.vessels))

# And unwrap the data, make sure it is the same size as the "port.dat" object, it's the same thing with data filled in...
port.sampling <- do.call("rbind",final.dat)

table(port.sampling$year[!(port.sampling$bank %in% "GBa")])

mia <- unique(port.sampling[is.na(port.sampling$bank),c("ves", "vrnum", "vesid", "ps.date", "ps.fished")])

unique(fish.dat[fish.dat$vesid==mia$vesid[1] & fish.dat$year == year(mia$ps.fished[1]),]$date)

unique(fish.dat[fish.dat$vesid==mia$vesid[1] & fish.dat$date == mia$ps.fished[1],])

table(port.sampling$year, port.sampling$bank)


# The below list is vessel port sampling that doesn't have fishing logs associated.  Generally these are cases in which the
# port sampling dates are not matching the fishing logs, for the most part the port sampling dates are transposed by a day (i.e. log says fishing started
# on the 21, but the port sampling says 20th, I adjust the port sampling information in these cases.)
# This list should be empty once I've QA/QC'ed

#broken.ps <- aggregate(fished~boat + date,port.sampling[(which(is.na(port.sampling$fleet))),],length)
#broken.ps
#table(port.sampling$fleet) # Hopefully all FT or WF vessels
#samples.per.vessel.day <- aggregate(fished~boat+date,port.sampling,length)
#summary(samples.per.vessel.day$fished)
#quantile(samples.per.vessel.day$fished,probs=seq(0.01,1,by=0.01))

# table(port.sampling$year, port.sampling$bank)


# Convert the individual meat weights into a meat count...
port.sampling$mc <- 500/(port.sampling$meat_weight)

# Drop trips in which we don't have the lat/lon and therefor the bank.
port.sampling <- port.sampling[-which(is.na(port.sampling$bank)),]

# Now make any vessel with ASM's...ASM's
ASMs <- off.fleet[!is.na(off.fleet$ASM_date),]
for(i in 1:nrow(ASMs)) port.sampling$fleet[port.sampling$boat == ASMs$ID_alt_Port_Sampling[i] & port.sampling$date >= ASMs$ASM_date[i]] <- "ASM"

save.image(paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/Port_sampling_processed_data_2022.RData"))


################  End Section 1, processing the port sampling meat weight information  --  End Section 1 ##############################################
################  End Section 1, processing the port sampling meat weight information  --  End Section 1 ##############################################


############### Section 2 - GB Port sampling analysis -- Section 2 ###############
load(paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/Port_sampling_processed_data_2022.rData"))

# For the models to come let's center on the day we believe the ASM's arrive on the scence (May 24, 2010 according to Ginette)
# The day here respresents 100 days.
port.sampling$mod.day <- as.numeric(port.sampling$ps.date - as.Date("2010-05-24"))/100
# And let's just look at GBa
gba.dat <- port.sampling[port.sampling$bank == "GBa" & port.sampling$year>2016,]

# Plot a time series of the port sampling data by fleet
#windows(11,11)
ggplot(gba.dat,aes(fished,meat_weight,colour=fleet)) + geom_point(size=0.2) + geom_smooth()

# Again by fleet, but lump the years together
#windows(11,11)
ggplot(gba.dat,aes(year,meat_weight,colour=fleet)) + geom_point(size=0.02) + geom_smooth() + facet_wrap(~ bank)

#windows(11,11)
ggplot(gba.dat,aes(ps.fished,mc,colour=fleet)) + geom_point(size=0.02) + geom_smooth()
theme_bw() + theme(panel.grid=element_blank()) + scale_color_manual(values=c(alpha("blue",0.5),alpha("grey",0.5),alpha("red",0.5)))

mc.by.year.fleet <- aggregate(mc~year+fleet,gba.dat,median)
mc.by.year.fleet

mc.by.month.fleet <- aggregate(mc~month+fleet,gba.dat,median)
mc.by.month.fleet

# Is there a clear seasonal trend in here
#windows(11,11)
ggplot(gba.dat,aes(mc,fill=fleet)) + geom_histogram() + facet_wrap(~month)+#,scales = "free") +
  theme_bw() + theme(panel.grid=element_blank()) + scale_fill_manual(values=c(alpha("blue",0.5),alpha("grey",0.5),alpha("red",0.5)))



# First the simple model, is there a linear trend in meat counts over time...
mod.1 <- lm(mc~mod.day,gba.dat)
# quick diagnostics check, can see that overall there is a trend toward smaller scallop, think this would
# be slightly more straightforward to interpret if the date wasn't a date but was numeric...
summary(mod.1) # SO on May 24 2010 this says the mc was 28.82 and that is has declined by 0.046 every 100 days since...

# The right skew of the data is somewhat evident in the Normal Q-Q plot, given these are strictly positive
# it might make sense to go glm with a gamma, poisson, or quasipoisson.
#windows(11,11)
par(mfrow=c(2,2))
plot(mod.1)

# Here's the lm fit to the data, a small deline over time...
#windows(11,11)
ggplot(gba.dat,aes(ps.fished,mc)) + geom_point(size=0.02)+geom_smooth(method="lm")

# Next we make a slighly more sensible model, still linear but looking for a different in trends between fleets...

mod.2 <- lm(mc~mod.day*fleet,gba.dat)
mod.2.me <- lm(mc~mod.day*fleet-mod.day,gba.dat) # This gives the difference between fleets on 2010-05-24 and regression coefficients by fleet since...
summary(mod.2) # so this says the the ASM mean on May 24 2010 was 31.3 and that is has declined by 0.085 every 100 days since
# Whereas the FT fleet mean was 2.86 lower, and the Wf fleet was 4.1 lower.  The decline for the FT was not significantly different than the ASM's
# but the WF fleet decline was significantly different from the ASM (and likely the FT's), the slope being essentially twice that of the ASM fleet
summary(mod.2.me) # Says same as above but gives actual regression coefficients for the slopes rather than the differences.

# Here's the lm fit to the data,
windows(11,11)
ggplot(gba.dat,aes(mod.day,mc,colour=fleet)) + geom_point(size=0.02)+geom_smooth(method="lm")

# Now we can move to a gam type of model to allow for non-linearity in the patterns....
mod.3 <- gamm(mc ~ s(mod.day,by=as.factor(fleet)),data=gba.dat)
summary(mod.3$lme)
summary(mod.3$gam)

# Let's take a closer look at these GAM results....
anova(mod.3$gam)# Careful with these p values as the smoother comes into play on these...
summary(mod.3$lme) #Will be explained later


#Model validation
windows(11,11)
par(mfrow=c(3,1))
plot(mod.3$gam,ylim=c(-10,10))
windows(11,11)
par(mfrow=c(3,1))
plot(mod.3$lme)

#par(mfrow = c(2,2), mar = c(5,5,2,2))
E1 <- resid(mod.3$lme, type ="n")
F1 <- fitted(mod.3$lme)

windows(11,11)
par(mfrow = c(2,2))
plot(x = F1,
     y = E1,
     xlab = "Fitted values",
     ylab ="Residuals",cex=0.01)
abline(h=0, lty=2)

plot(x = gba.dat$mod.day,
     y = E1,
     xlab = "time",
     ylab = "Residuals",cex=0.01)
abline(h = 0, lty = 2)

boxplot(E1 ~ fleet, data = gba.dat)
abline(h = 0, lty = 2)

# Now we can get some predictions from the model, this is much tidier than dealing with the gigantic dataset!!
# First lets get good days to predict on, we'll go for the first of every month...
range(port.sampling$date)
pred.dates <- ymd(paste0(sort(rep(2006:2017,12)),"-",rep(1:12,12),"-",rep(0,144),rep(1,144)))
mod.pred.day <- as.numeric(pred.dates - as.Date("2010-05-24"))/100
pred.dat <- ddply(gba.dat,
                .(fleet),
                summarize,
                mod.day = seq(min(mod.pred.day),
                              max(mod.pred.day),
                              length = length(mod.pred.day)),
                date = pred.dates)


mod.pred <- predict(mod.3$gam, newdata = pred.dat, se = TRUE)

# Glue them all together
pred.dat$mc <- mod.pred$fit
pred.dat$ub <- mod.pred$fit + 1.96 * mod.pred$se.fit
pred.dat$lb <- mod.pred$fit - 1.96 * mod.pred$se.fit
# Also stick the acutal meat weights on here...
pred.dat$mw <- 500/pred.dat$mc
pred.dat$mw.ub <- 500/pred.dat$ub
pred.dat$mw.lb <- 500/pred.dat$lb

# Remove all ASM predictions before June 1st 2010 since there were no active vessels before this date (well May 24th, but I'm rounding...)
pred.dat <- pred.dat[-which(pred.dat$fleet == "ASM" & pred.dat$date < "2010-06-01"),]
head(pred.dat)

windows(11,11)
ggplot(pred.dat, aes(date,mc,colour=fleet))  + theme(text = element_text(size=16)) + theme_bw() +
           geom_point(data = gba.dat, aes(y = mc, x = date,colour = fleet),shape = 16, size = 0.01,alpha=0.05)+
           geom_line()  + ylim(5,60) + ylab("Meat count") + scale_x_date(date_breaks="1 year") + xlab("") +
           geom_ribbon(data = pred.dat, aes(x = date, ymax = ub, ymin = lb, fill = fleet, color = fleet), alpha = 0.5)
ggsave(paste0(direct,"2018/Framework/Port_sampling/Fleet_port_sampling.png"),width=11,height=8.5)

# Now let's look at the Port sampling for August and May, when we have our survey on GBa.  Using this information along with the
# MW-SH relationship for each year we can determine what size scallop the fleet is targeting, especially groovy if our May and August
# estimates line up!
# Here's my logic on this analysis...
# 1:  There is a meat-weight shell height relationship in May and August, it may differ between these two periods, whether or not the difference is
#     "signifcant" we can figure out I think.
# 2:  This meat-weight shell height relationship holds for the whole bank during this time, so using the May/August MW-SH relationship and the port sampling
#     meat weights we can determine what size (i.e. shell hieghts) that the fleets are targeting during these 2 months.
# 3:  IF the shell heights targeted by each fleet in May and August don't differ significantly WE CAN MAKE THE ASSUMPTION that the fleet, on average, is likely
#     targeting the same size scallop all year
# 4:  With this assumption in place (and given the data it does appear this is a reasonable assumption) we can make MW-SH models
#     at a finer temporal scale (e.g. calclate the average MW by fleet for each month) and estimate condition at these finer temporal scales.
#     now we have an estimate of how condition varies over the year and if condition varies between fleets (though I'm not sure if the fleet question makes sense)

# First off I'll need to load in the survey results....
load(paste0(direct,"Data/Survey_data/2017/Survey_summary_output/Survey_all_results.RData"))

# Now what we want is to fit the MW-SH model to these data for each year back to 2006 and for both May and August survey...
# For the most recent data
years <- 2006:2017
aug.mws <- NULL
may.mws <- NULL
a.may <- data.frame(a = rep(NA,length(years)), year = years)
a.aug <- data.frame(a = rep(NA,length(years)), year = years)

# Now fit a linear model on the wmw v.s. shell height to find the intercept, note this is forced through 0.
# Now we can subset the predicted data from before and figure out what the targeted shell height is...
pred.dat[month(pred.dat$date) %in% c(5,8),]


for(i in 1:length(years))
{
  aug.mws[[as.character(years[i])]] <-  na.omit(mw[["GBa"]][mw[["GBa"]]$year == years[i],]) # get GBa and chuck the na's
  may.mws[[as.character(years[i])]] <-  na.omit(mw[["GB"]][mw[["GB"]]$year == years[i],]) # get GB and chuck the na's
  aug.mws[[as.character(years[i])]]$sh_3 <- (aug.mws[[as.character(years[i])]]$sh/100)^3 # cube the SH's
  may.mws[[as.character(years[i])]]$sh_3 <- (may.mws[[as.character(years[i])]]$sh/100)^3 # cube the SH's
  aug.mws[[as.character(years[i])]]$sh <- (aug.mws[[as.character(years[i])]]$sh/100) # sh in decimeters
  may.mws[[as.character(years[i])]]$sh <- (may.mws[[as.character(years[i])]]$sh/100) # sh in decimeters
  # Run the model and save the results
  a.aug$a[i] <-  lme(fixed = wmw ~ sh_3 -1, data = aug.mws[[as.character(years[i])]], random = ~ sh_3 -1 | tow, method="REML")$coefficients$fixed
  a.may$a[i] <-  lme(fixed = wmw ~ sh_3 -1, data = may.mws[[as.character(years[i])]], random = ~ sh_3 -1 | tow, method="REML")$coefficients$fixed
  # What is the shell height being targeted in may
  pred.dat$sh[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                    a.aug$a[a.aug$year == years[i]])^(1/3)
  pred.dat$sh.lb[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw.lb[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                    a.aug$a[a.aug$year == years[i]])^(1/3)
  pred.dat$sh.ub[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw.ub[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                       a.aug$a[a.aug$year == years[i]])^(1/3)
  # What is the shell height being targeted in August
  pred.dat$sh[month(pred.dat$date) == 8 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw[month(pred.dat$date) == 8 & year(pred.dat$date) == years[i]] /
                                                                                    a.aug$a[a.aug$year == years[i]])^(1/3)
  pred.dat$sh.lb[month(pred.dat$date) == 8 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw.lb[month(pred.dat$date) == 8 & year(pred.dat$date) == years[i]] /
                                                                                    a.aug$a[a.aug$year == years[i]])^(1/3)
  pred.dat$sh.ub[month(pred.dat$date) == 8 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw.ub[month(pred.dat$date) == 8 & year(pred.dat$date) == years[i]] /
                                                                                    a.aug$a[a.aug$year == years[i]])^(1/3)

}

# Now plot up the results...
lab <- as_labeller(c('5'= "May",'8' ="August"))
ggplot(pred.dat[month(pred.dat$date) %in% c(5,8),],aes(date,sh,colour=fleet)) + geom_point() + geom_line() +
             facet_wrap(~month(date),labeller = lab) +
             theme(text = element_text(size=16)) + theme_bw() + ylim(95,115) + ylab("Shell height") + scale_x_date(date_breaks="2 years") + xlab("")  +
             geom_ribbon(aes(x = date, ymax = sh.ub, ymin = sh.lb, fill = fleet, color = fleet), alpha = 0.5)
ggsave(paste0(direct,"2018/Framework/Port_sampling/sh_targeted.png"),width=11,height=8.5)


# Now based on this I think we can assume that within a fishing season the fleet tends to target roughly the same sized scallop, if so
# we can actual get a feel for how the MW-SH relationship varies over the course of a season by looking at changes in the MW's over the year.
# So for each year let's get the SH targeted as the average of the May and August surveys...
SH.target <- aggregate(sh ~ year(date) + fleet,pred.dat,mean)
names(SH.target) <- c("year","fleet","sh")
# Add an year month combo column to port.sampling...
port.sampling$year_month <- ymd(paste(year(port.sampling$ps.fished),month(port.sampling$ps.fished),15,sep="-"))
mw.by.month <- aggregate(meat_weight ~ year_month + year+fleet,port.sampling,mean)
names(mw.by.month) <- c("month","year","fleet","mw")


windows(11,11)
ggplot(mw.by.month,aes(month,mw,colour=fleet)) + geom_point() + geom_line()

ps.mw.sh <- merge(mw.by.month,SH.target,by =c("year","fleet"))

ps.mw.sh$mc <- 500/ps.mw.sh$meat_weight
# Condition here is kinda interesting as it standardizes the meat weight between the fleets, i.e. what is the
# meat weight the different fleets are picking up for a 100 mm scallop.  In theroy I wouldn't expect this to be any different...
# What the questions around condition are:
# 1:  How does condition vary throughout the year, when does condition appear to peak.
# 2:  For a 100 mm scallop, is there any difference between the fleets in terms of the size of the meat they capture.
# I'm not sure if my logic holds up here, but I think it does, hinges on the assumption that the targeted SH doesn't vary
# significantly throughout the year
ps.mw.sh$cond <- ps.mw.sh$meat_weight/(ps.mw.sh$sh/100)^3
ps.mw.sh$month <- as.factor(month(ps.mw.sh$year_month))

windows(11,11)
ggplot(ps.mw.sh,aes(year_month,cond,colour=fleet)) + geom_point() + geom_line() + facet_wrap(~year,scales="free")

windows(11,11)
ggplot(ps.mw.sh,aes(month,cond,group=month,colour=fleet)) + geom_boxplot() +  facet_wrap(~fleet)

windows(11,11)
ggplot(ps.mw.sh,aes(month,meat_weight,group=month,colour=fleet)) + geom_boxplot()  + facet_wrap(~fleet)

windows(11,11)
ggplot(ps.mw.sh,aes(month,mc,group=month,colour=fleet)) + geom_boxplot()  + facet_wrap(~fleet)

# So we can make some models to see if month/fleet effects...
# A condition model first, does condition differ between months, yep, looks like it
# Condition on average peaks in June is about 3 higher than in January, May, July and August all better than January (August only marginally so)
mod.1 <- lm(cond~month,ps.mw.sh)
summary(mod.1) # yep

# What about if we add in fleet, doesn't look like any differences, we can't see any difference in "condition" between the fleets, this
# makes loads of sense, basically there is a size selectivity but not a quality selectivity.
mod.2 <- lm(cond~month*fleet,ps.mw.sh)
summary(mod.2)

# And confirmed that adding in fleet isn't helpful here, basically we can see a difference
AIC(mod.1,mod.2)

############### End Section 2 - GB Port sampling analysis -- End Section 2 ###############
############### End Section 2 - GB Port sampling analysis -- End Section 2 ###############



############### Section 3 - BBn and Sable Port sampling analysis -- Section 3###############


load(paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/Port_sampling_processed_data_2022.rData"))

# For the models to come let's center on the day we believe the ASM's arrive on the scence (May 24, 2010 according to Ginette)
# The day here respresents 100 days.
port.sampling$mod.day <- as.numeric(port.sampling$ps.date - as.Date("2010-05-24"))/100
# And let's just look at Sable and BBn
bbn.dat <- port.sampling[port.sampling$bank == "BBn",]

# A plot of the time series for BBn and Sab

ggplot(bbn.dat,aes(fished,meat_weight,colour=fleet)) + geom_point(size=0.2) + geom_smooth()
#ggplot(sab.dat,aes(fished,meat_weight,colour=fleet)) + geom_point(size=0.2) + geom_smooth()

mc.by.year.fleet <- aggregate(mc~year+fleet,bbn.dat,median)
mc.by.year.fleet

mc.by.month.fleet <- aggregate(mc~month+fleet,bbn.dat,median)
mc.by.month.fleet



# Now we can move to a gam type of model to allow for non-linearity in the patterns....
bbn.mod <- gamm(mc ~ s(mod.day,by=as.factor(fleet)),data=bbn.dat)
summary(bbn.mod$lme)
summary(bbn.mod$gam)

# Let's take a closer look at these GAM results....
anova(bbn.mod$gam)# Careful with these p values as the smoother comes into play on these...
summary(bbn.mod) #Will be explained later


#Model validation
windows(11,11)
par(mfrow=c(3,1))
plot(bbn.mod$gam,ylim=c(-10,10))
windows(11,11)
par(mfrow=c(3,1))
plot(bbn.mod$lme)

#par(mfrow = c(2,2), mar = c(5,5,2,2))
E1 <- resid(bbn.mod$lme, type ="n")
F1 <- fitted(bbn.mod$lme)

windows(11,11)
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

# Is there a clear seasonal trend in here
#windows(11,11)
ggplot(bbn.dat,aes(mc,fill=fleet)) + geom_histogram() + facet_wrap(~month)+#,scales = "free") +
  theme_bw() + theme(panel.grid=element_blank()) + scale_fill_manual(values=c(alpha("blue",0.5),alpha("grey",0.5),alpha("red",0.5)))

# Lets set up to look at like MW and SH targeted by fishery....

pred.dates <- ymd(paste0(sort(rep(2006:2022,12)),"-",rep(1:12,12),"-",rep(0,144),rep(1,144)))
# For the models to come let's center on the day we believe the ASM's arrive on the scence (May 24, 2010 according to Ginette)
mod.pred.day <- as.numeric(pred.dates - as.Date("2010-05-24"))/100
pred.dat <- ddply(bbn.dat,
                  .(fleet),
                  summarize,
                  mod.day = seq(min(mod.pred.day),
                                max(mod.pred.day),
                                length = length(mod.pred.day)),
                  date = pred.dates)


mod.pred <- predict(bbn.mod$gam, newdata = pred.dat, se = TRUE)

# Glue them all together
pred.dat$mc <- mod.pred$fit
pred.dat$ub <- mod.pred$fit + 1.96 * mod.pred$se.fit
pred.dat$lb <- mod.pred$fit - 1.96 * mod.pred$se.fit
# Also stick the acutal meat weights on here...
pred.dat$mw <- 500/pred.dat$mc
pred.dat$mw.ub <- 500/pred.dat$ub
pred.dat$mw.lb <- 500/pred.dat$lb

# Remove all ASM predictions before June 1st 2010 since there were no active vessels before this date (well May 24th, but I'm rounding...)
pred.dat <- pred.dat[-which(pred.dat$fleet == "ASM" & pred.dat$date < "2010-06-01"),]
head(pred.dat)

windows(11,11)
ggplot(pred.dat, aes(date,mc,colour=fleet))  + theme(text = element_text(size=16)) + theme_bw() +
  geom_point(data = bbn.dat, aes(y = mc, x = date,colour = fleet),shape = 16, size = 0.01,alpha=0.05)+
  geom_line()  + ylim(5,60) + ylab("Meat count") + scale_x_date(date_breaks="1 year") + xlab("") +
  geom_ribbon(data = pred.dat, aes(x = date, ymax = ub, ymin = lb, fill = fleet, color = fleet), alpha = 0.5)
#ggsave(paste0(direct,"2018/Framework/Port_sampling/Fleet_port_sampling.png"),width=11,height=8.5)


# First off I'll need to load in the survey results....
#load("D:/testing_folder/data/Survey_all_results.RData")
load(paste0(direct,"Data/Survey_data/2022/Survey_summary_output/Survey_all_results.RData"))

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

  may.mws[[as.character(years[i])]] <-  na.omit(mw[["BBn"]][mw[["BBn"]]$year == years[i],]) # get BBn and chuck the na's
  may.mws[[as.character(years[i])]]$sh_3 <- (may.mws[[as.character(years[i])]]$sh/100)^3 # cube the SH's
  may.mws[[as.character(years[i])]]$sh <- (may.mws[[as.character(years[i])]]$sh/100) # sh in decimeters
  # Run the model and save the results
  a.may$a[i] <-  lme(fixed = wmw ~ sh_3 -1, data = may.mws[[as.character(years[i])]], random = ~ sh_3 -1 | tow, method="REML")$coefficients$fixed
  # What is the shell height being targeted in may
  pred.dat$sh[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                     a.may$a[a.may$year == years[i]])^(1/3)
  pred.dat$sh.lb[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw.lb[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                        a.may$a[a.may$year == years[i]])^(1/3)
  pred.dat$sh.ub[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw.ub[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                        a.may$a[a.may$year == years[i]])^(1/3)
}

# Now plot up the results...
lab <- as_labeller(c('5'= "May"))
ggplot(pred.dat[month(pred.dat$date) %in% c(5),],aes(date,sh,colour=fleet)) + geom_point() + geom_line() +
                                              theme(text = element_text(size=16)) + theme_bw() + ylim(95,115) +
                                              scale_x_date(date_breaks="2 years") + xlab("") + ylab("Shell height") +
                                              geom_ribbon(aes(x = date, ymax = sh.ub, ymin = sh.lb, fill = fleet, color = fleet), alpha = 0.5)
#ggsave(paste0(direct,"2018/Framework/Port_sampling/sh_targeted.png"),width=11,height=8.5)

# Now based on this I think we can assume that within a fishing season the fleet tends to target roughly the same sized scallop, if so
# we can actual get a feel for how the MW-SH relationship varies over the course of a season by looking at changes in the MW's over the year.
# So for each year let's get the SH targeted as the average of the May and August surveys...
SH.target <- aggregate(sh ~ year(date) + fleet,pred.dat,mean)
names(SH.target) <- c("year","fleet","sh")
# Add an year month combo column to port.sampling...
port.sampling$year_month <- ymd(paste(year(port.sampling$ps.fished),month(port.sampling$ps.fished),15,sep="-"))
mw.by.month <- aggregate(meat_weight ~ year_month + year+fleet,port.sampling,mean)
names(mw.by.month) <- c("date","year","fleet","mw")


windows(11,11)
ggplot(mw.by.month,aes(x= date,mw,colour=fleet)) + geom_point() + geom_line()

ps.mw.sh <- merge(mw.by.month,SH.target,by =c("year","fleet"))

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

windows(11,11)
ggplot(ps.mw.sh,aes(month,cond,colour=fleet)) + geom_point() + geom_line() + facet_wrap(~year)

windows(11,11)
ggplot(ps.mw.sh,aes(month,cond,group=month,colour=fleet)) + geom_boxplot() +  facet_wrap(~fleet)

windows(11,11)
ggplot(ps.mw.sh,aes(month,mw,group=month,colour=fleet)) + geom_boxplot()  + facet_wrap(~fleet)

windows(11,11)
ggplot(ps.mw.sh,aes(month,mc,group=month,colour=fleet)) + geom_boxplot()  + facet_wrap(~fleet)

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
  may.mws[[as.character(years[i])]] <-  na.omit(mw[["BBn"]][mw[["BBn"]]$year == years[i],]) # get BBn and chuck the na's
  may.mws[[as.character(years[i])]]$sh_3 <- (may.mws[[as.character(years[i])]]$sh/100)^3 # cube the SH's
  may.mws[[as.character(years[i])]]$sh <- (may.mws[[as.character(years[i])]]$sh/100) # sh in decimeters
  # Run the model and save the results
  a.may$a[i] <-  lme(fixed = wmw ~ sh_3 -1, data = may.mws[[as.character(years[i])]], random = ~ sh_3 -1 | tow, method="REML")$coefficients$fixed
  # What is the shell height being targeted in may
  bbn.dat$sh[bbn.dat$year == years[i]] <- 100*(bbn.dat$meat_weight[bbn.dat$year == years[i]]/ a.may$a[i])^(1/3)
}

sh.plt <- ggplot(bbn.dat,aes(date,sh,colour=fleet)) + geom_point(alpha = 0.05,size=0.2) +
                                                      theme(text = element_text(size=16)) + theme_bw() + geom_smooth(method = 'gam')+
                                                      scale_x_date(date_breaks="2 years") + scale_y_continuous(breaks = seq(0,200,by=5)) +
                                                      xlab("") + ylab("Shell height") + geom_hline(yintercept = 95, linetype = 'dashed',linewidth=1.5,color="grey")
sh.plt


# HERE IS THE MAIN ANALYSIS FOR the BBn Res doc, shows proportion below certain sizes
n.bl.100 <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 100) %>% dplyr::summarise(n100 = length(sh))
n.all.bbn <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::summarise(tot = length(sh))

prop.bl.100 <- left_join(n.bl.100,n.all.bbn,'year')
prop.bl.100$prop <- prop.bl.100$n100 / prop.bl.100$tot

windows(11,11)
ggplot(data = prop.bl.100,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 100 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))


n.bl.95 <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 95) %>% dplyr::summarise(n95 = length(sh))
prop.bl.95 <- left_join(n.bl.95,n.all.bbn,'year')
prop.bl.95$prop <- prop.bl.95$n95 / prop.bl.95$tot

windows(11,11)
ggplot(data = prop.bl.95,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 95 mm") +
                                               theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))

# Based on Ageing info, we think 75-90 mm is the best option, that would be ≈ 3-4 year olds.
n.bl.90 <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 90) %>% dplyr::summarise(n90 = length(sh))
prop.bl.90 <- left_join(n.bl.90,n.all.bbn,'year')
prop.bl.90$prop <- prop.bl.90$n90 / prop.bl.90$tot

# This is basically under 1%, which I think will be fine.
windows(11,11)
ggplot(data = prop.bl.90,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 90 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))


n.bl.85 <- bbn.dat %>% dplyr::group_by(year,.drop=F) %>% dplyr::filter(sh < 85) %>% dplyr::summarise(n85 = length(sh))
prop.bl.85 <- left_join(n.bl.85,n.all.bbn,'year')
prop.bl.85$prop <- prop.bl.85$n85 / prop.bl.85$tot

# So there is basically 0 catch coming from < 85 mm, well under 1%
windows(11,11)
ggplot(data = prop.bl.85,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 90 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))

# Just look at spring-summer, going with April-August as in theory condition should be fairly stable during this period.
# Based on Ageing info, we think 75-90 mm is the best option, that would be ≈ 3-4 year olds.
sel.months <- 4:8
n.bl.90.ss <- bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 90 & month %in% sel.months) %>% dplyr::summarise(n90 = length(sh))
n.all.bbn.ss <-  bbn.dat %>% dplyr::group_by(year) %>% dplyr::filter(month %in% sel.months) %>% dplyr::summarise(tot = length(sh))
prop.bl.90.ss <- left_join(n.bl.90.ss,n.all.bbn.ss,'year')
prop.bl.90.ss$prop <- prop.bl.90.ss$n90 / prop.bl.90.ss$tot

# This is basically under 1%, which I think will be fine.
windows(11,11)
ggplot(data = prop.bl.90.ss,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 90 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))


# Based on the ageing work, if we go with 90 mm as our end of recruitment bin, the recruits should be 75-90 mm in size.



####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################
# Now Sable...
sab.dat <- port.sampling[port.sampling$bank == "Sab",]

years <- c(2006:2014,2016:2019,2021:2022)
may.mws <- NULL
a.may <- data.frame(a = rep(NA,length(years)), year = years)


# Now fit a linear model on the wmw v.s. shell height to find the intercept, note this is forced through 0.
# Now we can subset the predicted data from before and figure out what the targeted shell height is...
# pred.dat[month(pred.dat$date) %in% c(5,8),]


for(i in 1:length(years))
{

  may.mws[[as.character(years[i])]] <-  na.omit(mw[["Sab"]][mw[["Sab"]]$year == years[i],]) # get BBn and chuck the na's
  may.mws[[as.character(years[i])]]$sh_3 <- (may.mws[[as.character(years[i])]]$sh/100)^3 # cube the SH's
  may.mws[[as.character(years[i])]]$sh <- (may.mws[[as.character(years[i])]]$sh/100) # sh in decimeters
  # Run the model and save the results
  a.may$a[i] <-  lme(fixed = wmw ~ sh_3 -1, data = may.mws[[as.character(years[i])]], random = ~ sh_3 -1 | tow, method="REML")$coefficients$fixed
  # What is the shell height being targeted in may
  pred.dat$sh[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                     a.may$a[a.may$year == years[i]])^(1/3)
  pred.dat$sh.lb[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw.lb[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                        a.may$a[a.may$year == years[i]])^(1/3)
  pred.dat$sh.ub[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] <- 100*(pred.dat$mw.ub[month(pred.dat$date) == 5 & year(pred.dat$date) == years[i]] /
                                                                                        a.may$a[a.may$year == years[i]])^(1/3)
}

# Now plot up the results...
ggplot(pred.dat[month(pred.dat$date) %in% c(5),],aes(date,sh,colour=fleet)) + geom_point() + geom_line() +
  theme(text = element_text(size=16)) + theme_bw() + ylim(95,115) +
  scale_x_date(date_breaks="2 years") + xlab("") + ylab("Shell height") +
  geom_ribbon(aes(x = date, ymax = sh.ub, ymin = sh.lb, fill = fleet, color = fleet), alpha = 0.5)
#ggsave(paste0(direct,"2018/Framework/Port_sampling/sh_targeted.png"),width=11,height=8.5)

# Now based on this I think we can assume that within a fishing season the fleet tends to target roughly the same sized scallop, if so
# we can actual get a feel for how the MW-SH relationship varies over the course of a season by looking at changes in the MW's over the year.
# So for each year let's get the SH targeted as the average of the May and August surveys...
SH.target <- aggregate(sh ~ year(date) + fleet,pred.dat,mean)
names(SH.target) <- c("year","fleet","sh")
# Add an year month combo column to port.sampling...
port.sampling$year_month <- ymd(paste(year(port.sampling$ps.fished),month(port.sampling$ps.fished),15,sep="-"))
mw.by.month <- aggregate(meat_weight ~ year_month + year+fleet,port.sampling,mean)
names(mw.by.month) <- c("date","year","fleet","mw")


windows(11,11)
ggplot(mw.by.month,aes(date,mw,colour=fleet)) + geom_point() + geom_line()

ps.mw.sh <- merge(mw.by.month,SH.target,by =c("year","fleet"))

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

windows(11,11)
ggplot(ps.mw.sh,aes(month,cond,colour=fleet)) + geom_point() + geom_line() + facet_wrap(~year)

windows(11,11)
ggplot(ps.mw.sh,aes(month,cond,group=month,colour=fleet)) + geom_boxplot() +  facet_wrap(~fleet)

windows(11,11)
ggplot(ps.mw.sh,aes(month,mw,group=month,colour=fleet)) + geom_boxplot()  + facet_wrap(~fleet)

windows(11,11)
ggplot(ps.mw.sh,aes(month,mc,group=month,colour=fleet)) + geom_boxplot()  + facet_wrap(~fleet)

# OK, so I want to see what propotion of the port sampling data might have a SH of < 95 mm, based on the MW-SH data.
# This is similar to above, but I'm going to have to use the raw data, so this is gonna be slow....

years <- c(2006:2014,2016:2019,2021:2022)
may.mws <- NULL
sab.dat$sh <- NA
all.may <- data.frame(a = rep(NA,length(years)), year = years)

pred.dates <- ymd(paste0(sort(rep(2006:2022,12)),"-",rep(1:12,12),"-",rep(0,144),rep(1,144)))
# For the models to come let's center on the day we believe the ASM's arrive on the scence (May 24, 2010 according to Ginette)
# This is needed to line up with mod.pred.day... I know this lumps the 31st with the day furthest away from the 31st rather than the next month
# But the idea here is everything in month x is treated the same.
sab.dat$year_month <- date(paste0(year(sab.dat$date),"-",sab.dat$month,'-01'))
sab.dat$mod.day <- as.numeric(sab.dat$year_month - as.Date("2010-05-24"))/100
sab.dat$year <- year(sab.dat$date)

for(i in 1:length(years))
{
  may.mws[[as.character(years[i])]] <-  na.omit(mw[["Sab"]][mw[["Sab"]]$year == years[i],]) # get BBn and chuck the na's
  may.mws[[as.character(years[i])]]$sh_3 <- (may.mws[[as.character(years[i])]]$sh/100)^3 # cube the SH's
  may.mws[[as.character(years[i])]]$sh <- (may.mws[[as.character(years[i])]]$sh/100) # sh in decimeters
  # Run the model and save the results
  a.may$a[i] <-  lme(fixed = wmw ~ sh_3 -1, data = may.mws[[as.character(years[i])]], random = ~ sh_3 -1 | tow, method="REML")$coefficients$fixed
  # What is the shell height being targeted in may
  sab.dat$sh[sab.dat$year == years[i]] <- 100*(sab.dat$meat_weight[sab.dat$year == years[i]]/ a.may$a[i])^(1/3)
}

# I don't know what is up with these data, but the smooths are a disaster aand a lm is dumb.
sh.plt <- ggplot(sab.dat,aes(date,sh,colour=fleet)) + geom_point(alpha = 0.05,size=0.2) +
  theme(text = element_text(size=16)) + theme_bw() +# geom_smooth(method = 'lm')+
  scale_x_date(date_breaks="2 years") + scale_y_continuous(breaks = seq(0,200,by=5)) +
  xlab("") + ylab("Shell height") + geom_hline(yintercept = 90, linetype = 'dashed',linewidth=1.5,color="grey")
sh.plt


# Here are the money plots for SAB
n.bl.100 <- sab.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 100) %>% dplyr::summarise(n100 = length(sh))
n.all.sab <- sab.dat %>% dplyr::group_by(year) %>% dplyr::summarise(tot = length(sh))

prop.bl.100 <- left_join(n.bl.100,n.all.sab,'year')
prop.bl.100$prop <- prop.bl.100$n100 / prop.bl.100$tot

windows(11,11)
ggplot(data = prop.bl.100,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 100 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))


n.bl.95 <- sab.dat %>% dplyr::group_by(year,.drop=F) %>% dplyr::filter(sh < 95) %>% dplyr::summarise(n95 = length(sh))

prop.bl.95 <- left_join(n.bl.95,n.all.sab,'year')
prop.bl.95$prop <- prop.bl.95$n95 / prop.bl.95$tot

windows(11,11)
ggplot(data = prop.bl.95,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 95 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))

n.bl.90 <- sab.dat %>% dplyr::group_by(year,.drop=F) %>% dplyr::filter(sh < 90) %>% dplyr::summarise(n90 = length(sh))
prop.bl.90 <- left_join(n.bl.90,n.all.sab,'year')
prop.bl.90$prop <- prop.bl.90$n90 / prop.bl.90$tot

windows(11,11)
ggplot(data = prop.bl.90,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 90 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))


n.bl.85 <- sab.dat %>% dplyr::group_by(year,.drop=F) %>% dplyr::filter(sh < 85) %>% dplyr::summarise(n85 = length(sh))
prop.bl.85 <- left_join(n.bl.85,n.all.sab,'year')
prop.bl.85$prop <- prop.bl.85$n85 / prop.bl.85$tot

# So there is basically 0 catch coming from < 85 mm, well under 1%
windows(11,11)
ggplot(data = prop.bl.85,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 90 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))

# Just look at spring-summer, going with April-August as in theory condition should be fairly stable during this period.
# Based on Ageing info, we think 75-90 mm is the best option, that would be ≈ 3-4 year olds.
sel.months <- 4:8
n.bl.90.ss <- sab.dat %>% dplyr::group_by(year) %>% dplyr::filter(sh < 90 & month %in% sel.months) %>% dplyr::summarise(n90 = length(sh))
n.all.sab.ss <-  sab.dat %>% dplyr::group_by(year) %>% dplyr::filter(month %in% sel.months) %>% dplyr::summarise(tot = length(sh))
prop.bl.90.ss <- left_join(n.bl.90.ss,n.all.sab.ss,'year')
prop.bl.90.ss$prop <- prop.bl.90.ss$n90 / prop.bl.90.ss$tot

# This is basically under 1%, which I think will be fine.
windows(11,11)
ggplot(data = prop.bl.90.ss,aes(x=year,y=prop)) + geom_point() + ylab("Proportion of meats below 90 mm") +
  theme(text = element_text(size=22)) + scale_x_continuous(breaks = seq(2000,2030,by=2))

# So I guess the next thing is to figure out what proportion the 90-95 mm size bins is of the Fully Recruited population?