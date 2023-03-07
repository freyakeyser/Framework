load("C:/Users/keyserf/Documents/temp_data/testing_results_framework.RData")

# prep data for TS plotting
bankpertow <- NULL

# Mid, Ban
for(bank in c("Mid", "Ban")){
  df <- SS.summary[[bank]] %>%
    dplyr::select(bank, year, n, N, NR, NPR, N.cv, NR.cv, NPR.cv, I, IR, IPR, I.cv, IR.cv, IPR.cv)
  bankpertow <- rbind(bankpertow, df)
}

dim(bankpertow)

# Ger
dim(merged.survey.obj)
bankpertow <- merged.survey.obj %>%
  dplyr::mutate(bank="Ger") %>%
  dplyr::select(bank, year, n, N, NR, NPR, N.cv, NR.cv, NPR.cv, I, IR, IPR, I.cv, IR.cv, IPR.cv) %>%
  full_join(bankpertow)

dim(bankpertow)

dim(survey.obj$Ger$model.dat)
bankpertow <- survey.obj$Ger$model.dat %>%
  dplyr::mutate(bank="Ger") %>%
  dplyr::select(bank, year, n, N, NR, NPR, N.cv, NR.cv, NPR.cv, I, IR, IPR, I.cv, IR.cv, IPR.cv) %>%
  full_join(bankpertow)

dim(bankpertow)

bankpertow <- full_join(bankpertow, data.frame(bank="Ger", year=2020))

# BBn, BBs
for(b in c("BBn", "BBs", "Sab")){
  print(dim(survey.obj[[b]]$model.dat))
  bankpertow <- survey.obj[[b]]$model.dat %>%
    dplyr::mutate(bank=b) %>%
    dplyr::select(bank, year, n, N, NR, NPR, N.cv, NR.cv, NPR.cv, I, IR, IPR, I.cv, IR.cv, IPR.cv) %>%
    full_join(bankpertow)

  bankpertow <- full_join(bankpertow, data.frame(bank=b, year=2020))

  print(dim(bankpertow))
}

bankpertow <- bankpertow %>%
  pivot_longer(cols=c("N", "NR", "NPR", "I", "IPR", "IR"), values_to = "index.val", names_to="index")

cv <- bankpertow %>%
  dplyr::select(bank, year, n, N.cv, NR.cv, NPR.cv, I.cv, IPR.cv, IR.cv) %>%
  pivot_longer(cols=c("N.cv", "NR.cv", "NPR.cv", "I.cv", "IPR.cv", "IR.cv"), values_to = "cv.val", names_to="index") %>%
  dplyr::select(bank, year, n, index, cv.val)
cv$index <- gsub(x=cv$index, ".cv", "")

bankpertow <- bankpertow %>%
  dplyr::select(bank, year, n, index, index.val) %>%
  full_join(cv)

bankpertow$sfa[bankpertow$bank %in% c("Ban", "Mid", "Sab")] <- 25
bankpertow$sfa[bankpertow$bank %in% c("BBn", "BBs", "Ger")] <- 26

bankpertow$size[bankpertow$index %in% c("N", "I")] <- "FR"
bankpertow$size[bankpertow$index %in% c("NPR", "IPR")] <- "PR"
bankpertow$size[bankpertow$index %in% c("NR", "IR")] <- "R"

bankpertow$index <- gsub(x=bankpertow$index, "P", "")
bankpertow$index <- gsub(x=bankpertow$index, "R", "")


sfa <- c(25,26)
for(sfa in sfa){
  png(filename=paste0("timeseries_", sfa, "_FR.png"), width=8, height=6, units = "in", res=420)
  print(ggplot() + geom_line(data=bankpertow[bankpertow$sfa==sfa & bankpertow$size=="FR",],
                             aes(year, index.val, colour=bank)) +
          geom_point(data=bankpertow[bankpertow$sfa==sfa & bankpertow$size=="FR",],
                     aes(year, index.val, colour=bank)) +
          #geom_errorbar(data=bankpertow, aes(year, ymax=NPR+(NPR*NPR.cv), ymin=NPR-(NPR*NPR.cv), colour=bank), width=0) +
          geom_ribbon(data=bankpertow[bankpertow$sfa==sfa & bankpertow$size=="FR",],
                      aes(year, ymax=index.val+(index.val*cv.val), ymin=index.val-(index.val*cv.val), fill=bank), alpha=0.3)+
          theme_bw() +
          # facet_grid(factor(index, levels=c("N", "I"))~factor(size, levels=c("PR", "R", "FR")),
          #            scales="free_y")
          facet_wrap(~factor(size, levels=c("PR", "R", "FR"))+factor(index, levels=c("N", "I")),
                     scales="free_y", ncol=1,
                     labeller = as_labeller(
                       c(PR = "Pre-recruits", R="Recruits", FR="",
                         I="g/tow", N="n/tow")),
                     strip.position = "left") +
          theme(strip.placement="outside", strip.background=element_blank()) +
          scale_colour_brewer(type="qual", name="Management\nunit") +
          scale_fill_brewer(type="qual", name="Management\nunit") +
          #scale_y_log10() +
          ylab(NULL)+
          xlab("Year"))
  dev.off()
}

