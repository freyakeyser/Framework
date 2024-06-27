# for comparing strata (strata_file_review.xlsx tab 2)
sf_use_s2(TRUE)
strata <- github_spatial_import(subfolder = "offshore_survey_strata", zipname = "offshore_survey_strata.zip", quiet = T)
strata[strata$label=="GBa", c("Strt_ID", "are_km2")]
st_area(st_transform(strata[strata$label=="GBa", c("Strt_ID", "are_km2")], 32619))/1000000

ggplot() + geom_sf(data=strata[strata$label=="GBa", c("Strt_ID", "are_km2")], fill="red") + facet_wrap(~Strt_ID)


source("C:/Users/keyserf/Documents/Github/Assessment_fns/Maps/pbs_2_sf.R")
csv <- survey.detail.polys[survey.detail.polys$label=="GBa",]
out <- pbs_2_sf(pbs = csv, lon="X", lat="Y")
ggplot() + geom_sf(data=out)
st_area(st_transform(out, 32619))/1000000

out$PID <- 1:7

survey.info[survey.info$label=="GBb",]

for(i in 1:7){
  print(ggplotly(
    ggplot() +
      geom_sf(data=out[out$PID==i,], colour="red", fill="red") +
      annotate(geom="text", x=Inf, y=Inf, label=round(st_area(st_transform(out[out$PID==i,], 32619))/1000000,2), colour="red", hjust=1, vjust=1) +
      geom_sf(data=strata[strata$label=="GBa" & strata$PID==i,], colour="black", fill="black") +
      annotate(geom="text", x=Inf, y=Inf, label=round(st_area(st_transform(strata[strata$label=="GBa" & strata$PID==i,], 32619))/1000000,2), colour="black", hjust=1, vjust=2) +
      ggtitle(i)
  ))
}

st_area(st_transform(strata[strata$label=="GBa" & strata$PID==6,], 32619))/1000000 # opposite
st_area(st_transform(out[out$PID==6,], 32619))/1000000

st_area(st_transform(strata[strata$label=="GBa" & strata$PID==5,], 32619))/1000000 # because of extra lines
st_area(st_transform(out[out$PID==5,], 32619))/1000000

st_area(st_transform(strata[strata$label=="GBa" & strata$PID==3,], 32619))/1000000 # opposite
st_area(st_transform(out[out$PID==3,], 32619))/1000000

st_area(st_transform(strata[strata$label=="GBa" & strata$PID==2,], 32619))/1000000 # opposite
st_area(st_transform(out[out$PID==2,], 32619))/1000000


# 5 and 6 have an extra line up top, 3 has a hole in csv and a line at the bottom, 2 has a line at the bottom

test <- st_difference(st_transform(out[out$PID==6,], 32619),
                      st_transform(strata[strata$label=="GBa" & strata$PID==6,], 32619))
ggplot() + geom_sf(data=test) +
  annotate(geom="text", x=Inf, y=Inf, label=st_area(test)/1000000, colour="black", hjust=1, vjust=2)

