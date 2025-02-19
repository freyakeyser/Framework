# Frist grab the directory to work from
direct <- "Y:/Offshore/Assessment/"
#direct <- "C:/Users/keyserf/Documents/temp_data/"
direct_fns <- "C:/Users/keyserf/Documents/Github/"
# Grab the years of port sampling data.  Note that Ginette says the Port Sampling data goes all the way back to 1996
# I haven't a clue where that data is, but would certainly be interesting to compare the last 20 years of data!!
years <- 2006:2022
options(stringsAsFactors = F) # Don't make anything a Factor as it screws up random crap.

# Load in librarys we may need...
library(ggplot2)
require(readxl)
require(xlsx)
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
logs_and_fish(loc = "offshore", year = years, get.local = T, get.marfis = F, direct=direct)
if(exists("old.log.dat") & exists("new.log.dat")) fish.dat<-merge(new.log.dat,old.log.dat,all=T)
if(exists("old.log.dat") & !exists("new.log.dat")) {
  fish.dat<-old.log.dat
  fish.dat$vrnum <- NA
  fish.dat$fished <- NA
}
if(!exists("old.log.dat") & exists("new.log.dat")) {
  fish.dat<-new.log.dat
  fish.dat$vesid <- NA
}
fish.dat$ID<-1:nrow(fish.dat)
################  Section 1, processing the port sampling meat weight information  --  Section 1 ############################################################
################  Section 1, processing the port sampling meat weight information  --  Section 1 ############################################################
# Here we pull in all the port sampling data (right now 2006-2017) and tidy it up for later analysis, once happy with the results
# you can skip this Section
require(tidyverse)

# Get ready for the loop...
index <- 0
dat <- NULL
split_days <- NULL
temp<-NULL
# Run this for all years we have data, this takes about 10 minutes...
for(i in 1:length(years))
{
  print(years[i])
  # 2013, COME0613 # this trip doesn't exist? # check 2006 files for this, I might have seen it there?
  if(years[i] < 2006) files <- list.files(path = paste0(direct,"Data/Archive/PortSampling/", years[i], "/"),pattern = "\\.txt", ignore.case = T)
  if(years[i] < 2018 & years[i]>2005) files <- list.files(path = paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/", years[i], "/"),pattern = "\\.xls", ignore.case = T)
  if(years[i] > 2017) files <- list.files(path = paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/", years[i], "/"),pattern = "\\.csv", ignore.case = T)

  # files <- files[files %in% c("Kelt0911.XLSX", "Kelt0923.XLSX", "Eepi0911.XLSX",
  # "Eepi0923.XLSX", "Ocea0911.XLSX", "Ocea0923.XLSX")]

  num.files <- length(files)

  for(j in 1:num.files)
  {
    print(index)
    index <- index + 1
    #if(index==204) browser()
    #This will pull the data from the Port Sampling file.
    if(years[i] < 2006) dat[[index]] <- read.csv(paste0(direct,"Data/Archive/PortSampling/", years[i], "/", files[j]), sep = "\t", header = F)
    if(years[i] < 2018 & years[i]>2005 & str_sub(string = direct, start=1, end=1)=="C") dat[[index]] <- read_excel(path=paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/",years[i],"/", files[j]), sheet = 1)
    if(years[i] < 2018 & years[i]>2005 & !str_sub(string = direct, start=1, end=1)=="C") dat[[index]] <- read.xlsx(file = paste0(direct,"Data/PortSampling/PS_data_reorg_4_analysis/",years[i],"/", files[j]), sheetIndex = 1)
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
    dat[[index]]$flag_datefished <- "N"

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
    temp[[index]] <- dat[[index]]

    #make boat uppercase
    dat[[index]]$boat <- toupper(dat[[index]]$boat)
    #check dates
    # fix bad dates using folder date
    if(!years[i] == unique(year(ymd(dat[[index]]$fished)))) {
      message(paste0("corrected year fished in ", files[j], " based on folder year"))
      dat[[index]]$fished <- paste0(years[i], str_sub(dat[[index]]$fished, 5,8))
      #dat[[index]]$date <- paste0(years[i], str_sub(dat[[index]]$date, 5,8))
    }
    # something might be fishy with date fished vs date sampled. They are months apart.
    if(any(!unique(month(ymd(dat[[index]]$fished))) %in% (unique(month(ymd(dat[[index]]$date)))-1):(unique(month(ymd(dat[[index]]$date)))+1))) {
      message(paste0("In ", files[j], ", month fished is ", unique(month(ymd(dat[[index]]$fished))), " but PS date is ", unique(ymd(dat[[index]]$date)), ". You better check this manually."))
      dat[[index]]$flag_datefished <- "Y"
    }

    # if(any(!year(ymd(dat[[index]]$date)) == year(ymd(dat[[index]]$fished)))){
    #   if(index==731) dat[[index]]$fished <- gsub(x=as.character(ymd(dat[[index]]$fished) + years(2)), "-", "")
    # }
#browser()
    # compare to logs to assign bank
    # Get the vessel of interest
    vessel <- na.omit(off.fleet[off.fleet$ID_alt_Port_Sampling == unique(dat[[index]]$boat), c("Pre_2008_ID","VMS_old_ID","ID")])
    # Get all the fishery data for that vessel, checking all types of identifiers
    ves.fish.dat <- fish.dat[(fish.dat$vesid %in% unlist(vessel) | fish.dat$ves %in% unlist(vessel) | fish.dat$vrnum %in% unlist(vessel)),]
    # Get the fishery data for that vessel on the appropriate dates
    ves.fish.dat <- ves.fish.dat[ves.fish.dat$date %in% ymd(dat[[index]]$fished),]
    bank <- unique(na.omit(ves.fish.dat$bank))
    # land <- unique(na.omit(ves.fish.dat$date.land))
    # if(length(land)==0){
    #   land <- ves.fish.dat %>% dplyr::group_by(trip.id) %>%
    #     dplyr::summarize(land=max(date, na.rm=T))
    #   land <- land$land
    # }
    trip.id <- unique(na.omit(ves.fish.dat$trip.id))
    if(length(bank)==0) {
      #browser()
      dat[[index]]$bank <- NA
      dat[[index]]$trip.id <- NA
    }
    if(length(bank)==1) {
      dat[[index]]$bank <- bank
      if(length(trip.id)==1) dat[[index]]$trip.id <- trip.id
      if(length(trip.id)>1) {
        dat[[index]]$trip.id <- max(trip.id)
        message("there are two trip.ids in fish.dat for this PS trip, but they're on the same bank, so I'm just using the higher one")
      }
    }
    if(length(bank)>1) {
      # generally split trips, join specific days, and drop transition day if it gets double counted (e.g. index 2185 20140821 PRES)
      ves.fish.dat <- unique(ves.fish.dat[, c("date", "bank", "trip.id")])
      ves.fish.dat$fished <- as.character(gsub(x=ves.fish.dat$date, "-", ""))
      dat[[index]]$fished <- as.character(dat[[index]]$fished)
      test <- dplyr::left_join(dat[[index]], dplyr::select(ves.fish.dat, bank, fished, trip.id))
      if(dim(test)[1] == dim(dat[[index]])[1]) dat[[index]] <- test
      if(!dim(test)[1] == dim(dat[[index]])[1]) {
        message("dropping split day")
        test2 <- test %>% dplyr::group_by(fished) %>%
          dplyr::summarize(banks = length(unique(bank)))
        split_days[[index]] <- unique(dplyr::select(test[test$fished %in% test2$fished[test2$banks>1],], -bank))
        split_days[[index]]$file <- files[j]
        test <- test[test$fished %in% test2$fished[test2$banks<2],]
        dat[[index]] <- test
      }
    }

    if(!dim(dat[[index]])[1] == dim(temp[[index]])[1]) message(paste0("File ", files[j], " - difference of ", abs(nrow(temp[[index]])-nrow(dat[[index]])), " records"))
    if(!dim(dat[[index]][complete.cases(dat[[index]]),])[1] == dim(temp[[index]][complete.cases(temp[[index]]),])[1]) message(paste0("File ", files[j], " - difference of ", abs(nrow(temp[[index]][complete.cases(temp[[index]]),])-nrow(dat[[index]][complete.cases(dat[[index]]),])), " records"))
    # dat[[index]]$file <- files[j]

    #check for duplicates
    if(any(duplicated(dat))){
      message(paste0("dropping duplicated file ", files[j]))
      dat <- dat[!which(duplicated(dat))]
    }

  } # end for(j in 1:num.files)

} # end the for(i in 1:length(years))

##### FINISH ADDING TRIP.ID SO THAT I CAN MATCH BACK TO SPATIAL LOCATIONS LATER! BUT ONLY IF IT'S IMPACTFUL!
#3h for 2006-2022, 3533 files

port.dat <- do.call("rbind",dat)
split_days <- do.call("rbind",split_days)

save(port.dat, file = "C:/Users/keyserf/Documents/temp_data/portdat_2006-2022.RData")
save(split_days, file = "C:/Users/keyserf/Documents/temp_data/splitdays_2006-2022.RData")
load("./Port_sampling/portdat_new.RData")
port.dat.new <- port.dat
port.dat.new$year <- year(ymd(port.dat.new$fished))
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

old2 <- port.dat.old %>% group_by(year,boat, ps.date, fished) %>%
  dplyr::summarize(nold=n(),
                   bankold=paste0(unique(bank), collapse=" "),
                   nbankold = length(unique(bank)))
old2[old2$nbankold>1,]
new2 <- port.dat.new %>% group_by(year,boat, ps.date, fished) %>%
  dplyr::summarize(nnew=n(),
                   banknew=paste0(unique(bank), collapse=" "),
                   nbanknew = length(unique(bank)))
new2[new2$nbanknew>1,]

test3 <- left_join(old2, new2)
any(!test3$nold==test3$nnew)

tail(test3[(!test3$banknew==test3$bankold) & !test3$banknew=="NA" &!is.na(test3$banknew),])
fish.dat$month <- month(ymd(fish.dat$date))
fish.dat[fish.dat$date=="2019-12-09" & fish.dat$ves=="ATLANTIC DESTINY",]
unique(fish.dat[fish.dat$year==2019 & fish.dat$month==12 & fish.dat$ves=="ATLANTIC DESTINY", c("tripnum", "date", "bank", "watch")]) %>% arrange(date, watch)
# definitely should be GBb instead of BBn

fish.dat[fish.dat$date=="2021-10-29" & fish.dat$ves=="ATLANTIC PROTECTOR",]
unique(fish.dat[fish.dat$year==2021 & fish.dat$month%in% 10:11 & fish.dat$ves=="ATLANTIC PROTECTOR", c("tripnum", "date", "bank", "watch")]) %>% arrange(date, watch)


# new has more than old?
unique(test3$nold - test3$nnew)
tail(unique(test3[is.na(test3$nnew) & !test3$year==2022,c("year", "bank", "boat","ps.date"),]))
#2019 Ger Prot 2019-09-04
test3[test3$ps.date=="2019-09-04" & test3$boat=="PROT",]
port.dat.new[port.dat.new$ps.date=="2019-09-04" & port.dat.new$boat=="PROT",]
port.dat.old[port.dat.old$ps.date=="2019-09-04" & port.dat.old$boat=="PROT",]
fish.dat[fish.dat$ves == "ATLANTIC PROTECTOR" & fish.dat$date %in% ymd("2019-08-20"),]
# split trip change day dropped

#2020 BBn   DEST  2020-01-10
test3[test3$ps.date=="2020-01-10" & test3$boat=="DEST",]
port.dat.new[port.dat.new$ps.date=="2020-01-10" & port.dat.new$boat=="DEST",]
port.dat.old[port.dat.old$ps.date=="2020-01-10" & port.dat.old$boat=="DEST",]
fish.dat[fish.dat$ves == "ATLANTIC DESTINY" & fish.dat$date %in% ymd("2019-08-20"),]