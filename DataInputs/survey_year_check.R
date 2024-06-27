require(tidyverse)
require(sf)

source("C:/Users/keyserf/Documents/Github/Assessment_fns/Fishery/logs_and_fishery_data.r")
logs_and_fish(loc="offshore",year = 1991:2022,un=un.ID,pw=pwd.ID,db.con=db.con,direct="Y:/Offshore/Assessment/", get.marfis=F)
fish.dat<-merge(new.log.dat,old.log.dat,all=T)
fish.dat$ID<-1:nrow(fish.dat)
# Now subset to BBn and add in the missing years of data
bbn.fish <- fish.dat %>% dplyr::filter(bank == "BBn")
#bbn.fish2 <- fish.dat[fish.dat$bank=="BBn",]
# There are 12 data points at 0,0 that we remove, I'm not worried about accounting for these 12 points!
bbn.fish <- bbn.fish %>% dplyr::filter(lat !=0 | lon != 0)
#bbn.fish2 <- bbn.fish2[!bbn.fish2$lat ==0,]
# bbn.fish2 <- bbn.fish2[!bbn.fish2$lon ==0,]
# bbn.fish2 <- bbn.fish2[!is.na(bbn.fish2$bank),]
# dim(bbn.fish2) #17971 (base)
dim(bbn.fish) # 17971 (tidy)
# Now I want to put a 'survey year' on these because that's what we're gonna need for our modelling... start by porting over the year
bbn.fish$survey.year <- bbn.fish$year
# Need to add fake data for 2009
bbn.fish[nrow(bbn.fish)+1,] <- NA

bbn.fish$pro.repwt[nrow(bbn.fish)] <- 0
bbn.fish$year[nrow(bbn.fish)] <- 2009
# See DK NOte below
bbn.fish$survey.year[nrow(bbn.fish)] <- 2009
bbn.fish$month[nrow(bbn.fish)] <- "June"
# Getting fake lat/lon coords that are on BBn
bbn.fish$lat[nrow(bbn.fish)] <- 42.85600
bbn.fish$lon[nrow(bbn.fish)]  <- -65.90183

require(lubridate)
month.num <- data.frame(month=month.name, month.num=1:12)
bbn.fish <- left_join(bbn.fish, month.num)
bbn.fish$month.num[is.na(bbn.fish$month.num)] <- month(bbn.fish$date[is.na(bbn.fish$month.num)])
month.num <-  data.frame(month.name=month.name, month.num=1:12)
bbn.fish <- left_join(bbn.fish, month.num)
bbn.fish$month[is.na(bbn.fish$month)] <- bbn.fish$month.name[is.na(bbn.fish$month)]
bbn.fish$doy <- yday(bbn.fish$date)

# DK NOTE: Now this is going to get confusing for us and we may want to tweak SEBDAM for this, but that's a down the road job, not a playing around with model job
# But based on the indexing in SEBDAM, I am going to change how we index the survey year data from what we have done with offshore traditionally.
# Historically anything from the last half of the year goes into the following years, eg. survey.year 2002 = June 2001- May 2002.
# But in SEBDAM we have (B(t-1) - C(t-1)), so let's say we have year 2002 survey biomass, this says we remove the 20002 catch from that
# we want that catch to be the catch from June 2002 to May 2003, i.e. we remove the catch before we allow the population to grow
# This is what we do in our current model, but we have a different index (C(t) on our model.
# Basically survey year 2002 = June 2002 - May 2003 now
#DK note: We probably should think more about the survey year fun and how exactly we want to handle removal of catch in our models.
# We don't have removals for 2009, we need something for that, so we're adding that in here...
bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] <- bbn.fish$survey.year[bbn.fish$month %in% c("January","February","March","April","May")] -1

bbn.fish.sf <- st_as_sf(bbn.fish,coords = c("lon","lat"),remove =F, crs = 4326)
bbn.fish.sf <- bbn.fish.sf %>% st_transform(crs= 32619)


# Now lets clip this to be data inside of our bbn boundary.
#bbn.fish.sf <- st_intersection(bbn.fish.sf,bbn.shape)

# Check removals each fishing year calculated using this data
bbn.fish.by.year <- bbn.fish.sf %>% dplyr::group_by(year) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))
bbn.fish.by.survey.year <- bbn.fish.sf %>% dplyr::group_by(survey.year,.drop=F) %>% dplyr::summarise(tot = sum(pro.repwt,na.rm=T))
# So this looks reasonable in the most recent years, but I probably need to check the early years to see if we are missing any of the removals, from above check (only 12 points removed) it
# seems like we might be fine, but need to check against our historical Removals estimates...
#tail(bbn.fish.by.year)
#tail(bbn.fish.by.survey.year)

# Subset the fishery data as necessary
years <- 1994:2022
bbn.fish <- bbn.fish %>% dplyr::filter(survey.year %in% years)
bbn.fish.sf <- bbn.fish.sf %>% dplyr::filter(survey.year %in% years)

color_mapping<- rep(c("#440154FF", "#3F4788FF", "#2D718EFF", "#29AF7FFF", "#94D840FF", "#FDE725FF"),5)

ggplot() + geom_text(data=unique(bbn.fish[, c("year", "month.num", "survey.year")]),
                     aes(year, as.factor(-month.num), label=survey.year)) +
  scale_y_discrete(labels=c(12:1)) +
  ylab("Month") +
  xlab("Fishing year") +
  ggtitle("Year labels are survey year (June to May)")

unique(bbn.fish[, c("month.num", "year", "survey.year")]) %>%
  arrange(year, month.num, survey.year) %>%
  filter(year>2019)
# survey year is june to may
