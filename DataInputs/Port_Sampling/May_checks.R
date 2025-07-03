# for ginette
load("./Port_sampling/portdat_2006-2022.RData")
require(lubridate)
require(tidyverse)
port.dat$year <- year(ymd(port.dat$fished))
gba <- port.dat[port.dat$bank=="GBa" & port.dat$year %in% c(2012, 2014),]
gba %>% dplyr::group_by(year) %>% dplyr::summarize(n())
gba_wide <- gba %>%
  dplyr::select(date, boat, port, id, fished, variable, bank, value) %>%
  dplyr::group_by(date, boat, port, id, fished, variable, bank) %>%
  dplyr::mutate(row = row_number()) %>%
  dplyr:: ungroup() %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  select(-row) %>%
  arrange(date, boat, fished, id) %>%
  dplyr::select(date, boat, port, id, fished, paste0("wt",1:14),bank)
str(gba_wide)

test <- gba_wide %>% pivot_longer(cols = starts_with("wt"))
test <- test[!is.na(test$value),]
test$year <- year(ymd(test$fished))
test %>% dplyr::group_by(year) %>% dplyr::summarize(n())
require(openxlsx)
write.xlsx(x = gba_wide, file = "C:/Users/keyserf/Documents/temp_data/PS_GBa_2012_2014.xlsx", na.string=0)



mid <- port.dat[port.dat$bank=="Mid" & port.dat$year %in% c(2007,2012,2016,2018),]
mid %>% dplyr::group_by(year) %>% dplyr::summarize(n())
mid_wide <- mid %>%
  dplyr::select(date, boat, port, id, fished, variable, bank, value) %>%
  dplyr::group_by(date, boat, port, id, fished, variable, bank) %>%
  dplyr::mutate(row = row_number()) %>%
  dplyr:: ungroup() %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  select(-row) %>%
  arrange(date, boat, fished, id) %>%
  dplyr::select(date, boat, port, id, fished, paste0("wt",1:14),bank)
str(mid_wide)

test <- mid_wide %>% pivot_longer(cols = starts_with("wt"))
test <- test[!is.na(test$value),]
test$year <- year(ymd(test$fished))
test %>% dplyr::group_by(year) %>% dplyr::summarize(n())
require(openxlsx)
write.xlsx(x = mid_wide, file = "C:/Users/keyserf/Documents/temp_data/PS_Mid.xlsx", na.string=0)

