
# trying to make sense of datclass

really.old.log.dat <- read.table(paste("Y:/Offshore/Assessment/Data/Fishery_data/Logs/1955_2008/LogData1955-1980.txt",sep=""),sep="\t",
                                 header=T,stringsAsFactors = F)
less.old.log.dat <- read.table(paste("Y:/Offshore/Assessment/Data/Fishery_data/Logs/1955_2008/LogData1981-2008.txt",sep=""),sep="\t",
                               header=T,stringsAsFactors = F)

table(is.na(really.old.log.dat$lon), really.old.log.dat$datclass)
table(is.na(really.old.log.dat$lat), really.old.log.dat$datclass)
table(really.old.log.dat$lon==0, really.old.log.dat$datclass)
table(really.old.log.dat$lat==0, really.old.log.dat$datclass)
table(really.old.log.dat$pro.repwt==0, really.old.log.dat$datclass)
table(is.na(really.old.log.dat$pro.repwt), really.old.log.dat$datclass)
table(really.old.log.dat$hm==0, really.old.log.dat$datclass)

ggplot() + geom_point(data=really.old.log.dat, aes(lon, lat)) + facet_wrap(~datclass)
ggplot() + geom_histogram(data=really.old.log.dat, aes(pro.repwt)) + facet_wrap(~datclass)


ggplot() + geom_histogram(data=really.old.log.dat, aes(hm)) + facet_wrap(~datclass)

#7 means coords are 0,0
#9 might mean 0 catch recorded?

table(is.na(less.old.log.dat$lon), less.old.log.dat$datclass)
table(is.na(less.old.log.dat$lat), less.old.log.dat$datclass)
table(less.old.log.dat$lon==0, less.old.log.dat$datclass)
table(less.old.log.dat$lat==0, less.old.log.dat$datclass)
table(less.old.log.dat$pro.repwt==0, less.old.log.dat$datclass)
table(is.na(less.old.log.dat$pro.repwt), less.old.log.dat$datclass)
table(less.old.log.dat$hm==0, less.old.log.dat$datclass)

# 5 might mean coordinate issues (0s)
# 9 might mean 0 catch recorded?

logs_and_fish(loc="offshore", year=1955:2009, get.marfis = F, export = F, direct = direct)
table(old.log.dat$datclass)
table(new.log.dat$datclass)


jack.dat<-with(subset(old.log.dat,bank %in% "BBn" & year %in% 1995:2010 & datclass==1),data.frame(year=year,catch=pro.repwt,effort=hm))
#Source1 source("Y:/Offshore/Assessment/2014/r/fn/jackknife.r",local=T)
# Look at jackkniffe function along with Smith(1980) to understand what 'se' v.s. 'sd' options really mean
cpue.dat<-jackknife(jack.dat)[,-(3:4)]



