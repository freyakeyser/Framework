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
fish.dat<-merge(new.log.dat,old.log.dat,all=T)
fish.dat$ID<-1:nrow(fish.dat)
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

    #check dates
    # if(any(!year(ymd(dat[[index]]$date)) == year(ymd(dat[[index]]$fished)))){
    #   if(index==731) dat[[index]]$fished <- gsub(x=as.character(ymd(dat[[index]]$fished) + years(2)), "-", "")
    # }

    # compare to logs to assign bank
    # Get the vessel of interest
    vessel <- na.omit(off.fleet[off.fleet$ID_alt_Port_Sampling == toupper(unique(dat[[index]]$boat)), c("Pre_2008_ID","VMS_old_ID","ID")])
    # Get all the fishery data for that vessel, checking all types of identifiers
    ves.fish.dat <- fish.dat[fish.dat$date %in% ymd(dat[[index]]$fished) & (fish.dat$vesid %in% unlist(vessel) | fish.dat$ves %in% unlist(vessel) | fish.dat$vrnum %in% unlist(vessel)),]
    bank <- unique(na.omit(ves.fish.dat$bank))
    land <- unique(na.omit(ves.fish.dat$date.land))
    if(length(bank)==0) dat[[index]]$bank <- NA
    if(length(bank)==1) dat[[index]]$bank <- bank
    if(length(bank)>1) {
      # generally split trips, join specific days, and drop transition day if it gets double counted (e.g. index 2185 20140821 PRES)
      ves.fish.dat <- unique(ves.fish.dat[, c("date", "bank")])
      ves.fish.dat$fished <- as.character(gsub(x=ves.fish.dat$date, "-", ""))
      dat[[index]]$fished <- as.character(dat[[index]]$fished)
      test <- dplyr::left_join(dat[[index]], dplyr::select(ves.fish.dat, bank, fished))
      if(dim(test)[1] == dim(dat[[index]])[1]) dat[[index]] <- test
      if(!dim(test)[1] == dim(dat[[index]])[1]) {
        test2 <- test %>% group_by(fished) %>%
          dplyr::summarize(banks = length(unique(bank)))
        test <- test[test$fished %in% test2$fished[test2$banks<2],]
        dat[[index]] <- test
      }
    }
  } # end for(j in 1:num.files)
} # end the for(i in 1:length(years))

port.dat <- do.call("rbind",dat)

#save(port.dat, file = "./Port_sampling/portdat_new.RData")
load("./Port_sampling/portdat_new.RData")
port.dat.new <- port.dat
port.dat.new$year <- year(ymd(port.dat.new$date))
#table(port.dat.new$year, port.dat.new$bank)
port.dat.new$ps.date <- port.dat.new$date
port.dat.new$ps.date <- ymd(port.dat.new$ps.date)
port.dat.new$meat_weight <- port.dat.new$value
port.dat.new$meat_weight <- port.dat.new$meat_weight/100
port.dat.new <- dplyr::select(port.dat.new, -date, -value)
# load(paste0("C:/Users/keyserf/Documents/temp_data/Data/PortSampling/PS_data_reorg_4_analysis/Port_sampling_raw_data_2022.RData"))
# port.dat.old <- port.dat

direct <- "Y:/Offshore/Assessment/"
load(paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/Port_sampling_processed_data_2022.rData"))
port.dat.old <- port.sampling
port.dat.old <- dplyr::select(port.dat.old, names(port.dat.new))

dim(port.dat.old)
dim(port.dat.new)

old <- port.dat.old %>% group_by(year,bank) %>%
  dplyr::summarize(n=n(),trips=length(unique(ps.date)),
                   type="old")
new <- port.dat.new %>% group_by(year,bank) %>%
  dplyr::summarize(n=n(),trips=length(unique(ps.date)),
                   type="new")
both <- rbind(old, new)
ggplot() + geom_point(data=both, aes(year, n, colour=type)) + facet_wrap(~bank,scales="free_y")
#old has more than new

both[both$bank=="GBb" & both$year==2020,]

str(port.dat.new)
str(port.dat.old)

old2 <- port.dat.old %>% group_by(year,bank, boat, ps.date, fished) %>%
  dplyr::summarize(nold=n())
new2 <- port.dat.new %>% group_by(year,bank, boat, ps.date, fished) %>%
  dplyr::summarize(nnew=n())

test3 <- left_join(old2, new2)
any(!test3$nold==test3$nnew)

ggplot() +
  geom_point(data=test3, aes(fished, nold), colour="red") +
  geom_point(data=test3, aes(fished, nnew), colour="blue") +
  facet_wrap(~bank)
# new has more than old?
unique(test3$nold - test3$nnew)
tail(unique(test3[is.na(test3$nnew) & !test3$year==2022,c("year", "bank", "boat","ps.date"),]))
#2021 GBa   LAHA  2021-12-06
test3[test3$ps.date=="2021-12-06" & test3$boat=="LAHA",]
port.dat.new[port.dat.new$ps.date=="2021-12-06" & port.dat.new$boat=="LAHA",]
#not in port.dat.new why???
fish.dat[fish.dat$ves == "CAPE LAHAVE" & fish.dat$date %in% ymd("2021-11-17"):ymd("2021-12-05"),]
# it is in fish.dat.