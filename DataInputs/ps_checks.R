require(tidyverse)
port.dat %>%
  dplyr::filter(year==2022) %>%
  dplyr::group_by(boat, ps.date) %>%
  dplyr::summarize(earliest = min(ymd(ps.fished)),
            latest = max(ymd(ps.fished))) %>% as.data.frame()


fish.dat %>%
  dplyr::filter(year==2022) %>%
  dplyr::select(ves, date.sail, date.land, bank) %>%
  dplyr::distinct()
# %>%
#   dplyr::group_by(ves, date.sail, date.land) %>%
#   dplyr::summarize(nbank = length(unique(bank))) %>%
#   dplyr::filter(nbank>1) %>%
#   as.data.frame()

fish.dat %>%
  dplyr::filter(year==2022 & ves=="MAUDE ADAMS" & date.sail=="2022-08-08") %>%
  dplyr::select(ves, date.sail, date.land, fished, bank)

2022-08-14 should be dropped from maude adams 2022-08-08 to 22 trip

View(port.sampling[port.sampling$ps.fished=="2022-08-14" & port.sampling$boat=="MAUD",])

75 samples removed from 2022-08-14 GBa, no effect on GBb count because it was previously assigning GBa.

View(port.sampling.1[port.sampling.1$ps.fished=="2022-08-14" & port.sampling.1$boat=="MAUD",])


fold22gbb <- list.files("Y:/Offshore/Assessment/Data/PortSampling/PS2022/JoansOriginals/")

test22gbb <- NULL
for(i in fold22gbb){
  files22gbb <- list.files(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2022/JoansOriginals/", i))
  files22gbb <- files22gbb[!grepl(x=files22gbb, pattern=".csv")]
  files22gbb <- files22gbb[!grepl(x=files22gbb, pattern=".doc")]
  files22gbb <- files22gbb[grepl(x=files22gbb, pattern="GBb")]

  for(j in files22gbb){
    files22gbb2 <- list.files(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2022/JoansOriginals/", i, "/", j))
    files22gbb2 <- files22gbb2[grepl(x=files22gbb2, pattern=".csv")]

    for(k in files22gbb2){
      psdat <- read.csv(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2022/JoansOriginals/", i, "/", j, "/", k))
      psdat$i <- i
      psdat$j <- j
      psdat$k <- k
      test22gbb <- rbind(test22gbb, psdat)
    }
  }
}
test22gbb <- test22gbb %>% pivot_longer(cols=starts_with("wt"))
dim(test22gbb)
test22gbb <- test22gbb[!test22gbb$value==0,]
dim(test22gbb) #15071 samples for GBb 2022
test22gbb <- dplyr::rename(test22gbb, variable="name", meat_weight="value", ps.date2="date")
test22gbb <- test22gbb %>% group_by(date, boat) %>%
  dplyr::mutate(sampid = 1:n()) %>%
  dplyr::mutate(unid = paste0(boat, ".", date, ".", sampid))
test22gbb$meat_weight <- test22gbb$meat_weight/100


dim(port.sampling[port.sampling$bank=="GBb" & !is.na(port.sampling$bank) & port.sampling$year==2022,])
comp22gbb <- port.sampling[port.sampling$bank=="GBb" & !is.na(port.sampling$bank) & port.sampling$year==2022,]

names(comp22gbb)
names(test22gbb)
test22gbb$ps.date2 <- ymd(test22gbb$ps.date2)
test22gbb$fished <- as.character(test22gbb$fished)
test22gbb$boat <- toupper(test22gbb$boat)

joined <- left_join(test22gbb, comp22gbb)
dim(joined)
joined[!is.na(joined$bank),] #12833 rows # 6 missing samples, the rest just weren't available yet!
unique(joined[is.na(joined$bank),c("ps.date2", "boat")]) #late fall Comeau trips - I bet we got these data too late!
names(test22gbb)[names(test22gbb) %in% names(comp22gbb)]

test22gbb2 <- unique(dplyr::select(test22gbb, boat, fished, meat_weight))
comp22gbb2 <- unique(dplyr::select(comp22gbb, boat, fished, meat_weight, bank))
joined2 <- left_join(test22gbb2, comp22gbb2)
joined2[is.na(joined2$bank),]
sort(unique(comp22gbb2$fished))
sort(unique(test22gbb2$fished))

# so do the same thing for 2021 now...
port.dat %>%
  dplyr::filter(year==2021) %>%
  dplyr::group_by(boat, ps.date) %>%
  dplyr::summarize(earliest = min(ymd(ps.fished)),
                   latest = max(ymd(ps.fished))) %>% as.data.frame()


fish.dat %>%
  dplyr::filter(year==2021) %>%
  dplyr::select(ves, date.sail, date.land, bank) %>%
  dplyr::distinct()
# %>%
#   dplyr::group_by(ves, date.sail, date.land) %>%
#   dplyr::summarize(nbank = length(unique(bank))) %>%
#   dplyr::filter(nbank>1) %>%
#   as.data.frame()

fold21gbb <- list.files("Y:/Offshore/Assessment/Data/PortSampling/PS2021/JoansOriginals/")

test21gbb <- NULL
for(i in fold21gbb){
  files21gbb <- list.files(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2021/JoansOriginals/", i))
  files21gbb <- files21gbb[!grepl(x=files21gbb, pattern=".csv")]
  files21gbb <- files21gbb[!grepl(x=files21gbb, pattern=".doc")]
  files21gbb <- files21gbb[grepl(x=files21gbb, pattern="GBb")]

  for(j in files21gbb){
    files21gbb2 <- list.files(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2021/JoansOriginals/", i, "/", j))
    files21gbb2 <- files21gbb2[grepl(x=files21gbb2, pattern=".csv")]

    for(k in files21gbb2){
      psdat <- read.csv(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2021/JoansOriginals/", i, "/", j, "/", k))
      psdat$i <- i
      psdat$j <- j
      psdat$k <- k
      test21gbb <- rbind(test21gbb, psdat)
    }
  }
}
test21gbb <- test21gbb %>% pivot_longer(cols=starts_with("wt"))
dim(test21gbb)
test21gbb <- test21gbb[!test21gbb$value==0,]
dim(test21gbb) #15071 samples for GBb 2021
test21gbb <- dplyr::rename(test21gbb, variable="name", meat_weight="value", ps.date2="date")
test21gbb <- test21gbb %>% group_by(ps.date2, boat) %>%
  dplyr::mutate(sampid = 1:n()) %>%
  dplyr::mutate(unid = paste0(boat, ".", ps.date2, ".", sampid))
test21gbb$meat_weight <- test21gbb$meat_weight/100


dim(port.sampling[port.sampling$bank=="GBb" & !is.na(port.sampling$bank) & port.sampling$year==2021,]) # 10718
comp21gbb <- port.sampling[port.sampling$bank=="GBb" & !is.na(port.sampling$bank) & port.sampling$year==2021,]

names(comp21gbb)
names(test21gbb)
test21gbb$ps.date2 <- ymd(test21gbb$ps.date2)
test21gbb$fished <- as.character(test21gbb$fished)
test21gbb$boat <- toupper(test21gbb$boat)

joined <- left_join(test21gbb, comp21gbb)
dim(joined)
joined[is.na(joined$bank),] #11052 rows # 322 missing from port.sampling df
unique(joined[is.na(joined$bank),c("ps.date2", "boat")]) # by PS date
unique(joined[is.na(joined$bank),c("fished", "boat")]) # by fishing date
names(test21gbb)[names(test21gbb) %in% names(comp21gbb)]
fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2021-07-22" & joined$boat=="FORT",]$fished)) & fish.dat$ves=="FORTUNE LADY",]$bank #
# something is fishy with the above... it was put in GBb folder, but report and log data say it is GBa only
# so this was correctly dropped from PS data

#fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2021-07-06" & joined$boat=="COME",]$fished)) & fish.dat$ves=="LADY COMEAU III",] #
#fixed

fish.dat[fish.dat$fished %in% ymd("20210427") & fish.dat$ves=="ATLANTIC PRESERVER",]$bank
# apparently port samples were collected on the 27th, but the log data has no fishing that day (sounds like they unexpectedly went in for the day)
fish.dat[fish.dat$fished=="2021-04-26" & !is.na(fish.dat$fished),]
sort(unique(fish.dat[fish.dat$tripnum=="548851",]$fished))
# port sampling/log data mismatch on one day, irreconcilable
# cannot reconcile via PS/log matching method

fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2021-08-02" & joined$boat=="FUND",]$fished)) & fish.dat$ves=="FUNDY LEADER",]$bank #
# ps date is wrong, should be 2021-09-02 I think (fishing occurred 2021-08-28 to 2021-09-01)
unique(fish.dat[fish.dat$fished=="2021-09-01" & !is.na(fish.dat$fished),]$ves)
# last watch dropped from log data for tripnum 557609
# cannot reconcile via PS/log matching method

sort(fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2021-11-03" & joined$boat=="PRES",]$fished)) & fish.dat$ves=="ATLANTIC PRESERVER",]$fished)
unique(fish.dat[fish.dat$fished=="2021-10-26" & !is.na(fish.dat$fished),]$ves)
fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2021-11-03" & joined$boat=="PRES",]$fished)) & fish.dat$ves=="ATLANTIC PRESERVER",]
# 2 watches dropped from log data for tripnum 563112 likely since they stopped fishing mid-day for a couple days due to storm
# cannot reconcile via PS/log matching method

test21gbb2 <- unique(dplyr::select(test21gbb, boat, fished, meat_weight))
comp21gbb2 <- unique(dplyr::select(comp21gbb, boat, fished, meat_weight, bank))
joined2 <- left_join(test21gbb2, comp21gbb2)
joined2[is.na(joined2$bank),]
sort(unique(comp21gbb2$fished))
sort(unique(test21gbb2$fished))


# 2020 checks
fold20gbb <- list.files("Y:/Offshore/Assessment/Data/PortSampling/PS2020/JoansOriginals/")

test20gbb <- NULL
for(i in fold20gbb){
  files20gbb <- list.files(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2020/JoansOriginals/", i))
  files20gbb <- files20gbb[!grepl(x=files20gbb, pattern=".csv")]
  files20gbb <- files20gbb[!grepl(x=files20gbb, pattern=".doc")]
  files20gbb <- files20gbb[grepl(x=files20gbb, pattern="GBb")]

  for(j in files20gbb){
    files20gbb2 <- list.files(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2020/JoansOriginals/", i, "/", j))
    files20gbb2 <- files20gbb2[grepl(x=files20gbb2, pattern=".csv")]

    for(k in files20gbb2){
      psdat <- read.csv(paste0("Y:/Offshore/Assessment/Data/PortSampling/PS2020/JoansOriginals/", i, "/", j, "/", k))
      psdat$i <- i
      psdat$j <- j
      psdat$k <- k
      test20gbb <- rbind(test20gbb, psdat)
    }
  }
}
test20gbb <- test20gbb %>% pivot_longer(cols=starts_with("wt"))
dim(test20gbb)
test20gbb <- test20gbb[!test20gbb$value==0,]
dim(test20gbb) #18149 samples for GBb 2020
test20gbb <- dplyr::rename(test20gbb, variable="name", meat_weight="value", ps.date2="date")
test20gbb <- test20gbb %>% group_by(ps.date2, boat) %>%
  dplyr::mutate(sampid = 1:n()) %>%
  dplyr::mutate(unid = paste0(boat, ".", ps.date2, ".", sampid))
test20gbb$meat_weight <- test20gbb$meat_weight/100


dim(port.sampling[port.sampling$bank=="GBb" & !is.na(port.sampling$bank) & port.sampling$year==2020,]) # 18072
comp20gbb <- port.sampling[port.sampling$bank=="GBb" & !is.na(port.sampling$bank) & port.sampling$year==2020,]

names(comp20gbb)
names(test20gbb)
test20gbb$ps.date2 <- ymd(test20gbb$ps.date2)
test20gbb$fished <- as.character(test20gbb$fished)
test20gbb$boat <- toupper(test20gbb$boat)

joined <- left_join(test20gbb, comp20gbb)
dim(joined)
joined[is.na(joined$bank),] #18181 rows # 77 missing from port.sampling df
unique(joined[is.na(joined$bank),c("ps.date2", "boat")]) # by PS date
unique(joined[is.na(joined$bank),c("fished", "boat")]) # by fishing date
names(test20gbb)[names(test20gbb) %in% names(comp20gbb)]
fish.dat[fish.dat$fished %in% "2020-11-29",] & fish.dat$ves=="FORTUNE LADY",]$bank #
# something is fishy with the above... it was put in GBb folder, but report and log data say it is GBa only
# so this was correctly dropped from PS data

#fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2020-07-06" & joined$boat=="COME",]$fished)) & fish.dat$ves=="LADY COMEAU III",] #
#fixed

fish.dat[fish.dat$fished %in% ymd("20200427") & fish.dat$ves=="ATLANTIC PRESERVER",]$bank
# apparently port samples were collected on the 27th, but the log data has no fishing that day (sounds like they unexpectedly went in for the day)
fish.dat[fish.dat$fished=="2020-04-26" & !is.na(fish.dat$fished),]
sort(unique(fish.dat[fish.dat$tripnum=="548851",]$fished))
# port sampling/log data mismatch on one day, irreconcilable
# cannot reconcile via PS/log matching method

fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2020-08-02" & joined$boat=="FUND",]$fished)) & fish.dat$ves=="FUNDY LEADER",]$bank #
# ps date is wrong, should be 2020-09-02 I think (fishing occurred 2020-08-28 to 2020-09-01)
unique(fish.dat[fish.dat$fished=="2020-09-01" & !is.na(fish.dat$fished),]$ves)
# last watch dropped from log data for tripnum 557609
# cannot reconcile via PS/log matching method

sort(fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2020-11-03" & joined$boat=="PRES",]$fished)) & fish.dat$ves=="ATLANTIC PRESERVER",]$fished)
unique(fish.dat[fish.dat$fished=="2020-10-26" & !is.na(fish.dat$fished),]$ves)
fish.dat[fish.dat$fished %in% ymd(unique(joined[joined$ps.date2=="2020-11-03" & joined$boat=="PRES",]$fished)) & fish.dat$ves=="ATLANTIC PRESERVER",]
# 2 watches dropped from log data for tripnum 563112 likely since they stopped fishing mid-day for a couple days due to storm
# cannot reconcile via PS/log matching method

test20gbb2 <- unique(dplyr::select(test20gbb, boat, fished, meat_weight))
comp20gbb2 <- unique(dplyr::select(comp20gbb, boat, fished, meat_weight, bank))
joined2 <- left_join(test20gbb2, comp20gbb2)
joined2[is.na(joined2$bank),]
sort(unique(comp20gbb2$fished))
sort(unique(test20gbb2$fished))