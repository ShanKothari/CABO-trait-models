setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels/")
library(signal)
library(dplyr)
library(spectrolab)
library(tidyr)
library(stringr)
library(asdreader)

#########################################################
## Processing leaf-level spectra for all CABO projects

## Here, I write separate sections for each projects
## There's a lot of copy-pasting
## I choose to keep it that way so that I can easily go back
## and redo a single project without fussing around too much

## I delete the raw data from my computer by they can always
## be found at:
## http://data.caboscience.org/field-data/projects/

##########################################################
## defining important functions for the tidyverse bullshit
## that happens later

## this just applies a Savitzy-Golay filter with given p and n
sg_filter <- function(x, p, n) {
  x$value_sg <- sgolayfilt(x$value, p, n)
  return(x)  
}

## this function is for linear interpolation over the spectral overlap region
interpolate <- function(x) {
  wvls <- x$wavelength
  DNs <- x$calculated_value
  new_DNs <- approx(wvls, DNs, xout = inter_wvls)$y
  tmp <- data_frame(wvl = inter_wvls, DN = new_DNs)
  return(tmp)
}

## to reset to just these two functions
rm.list=setdiff(ls(), c("interpolate","sg_filter"))
rm(list=rm.list)

###########################################################
## Beauchamp-Rioux spectra

BeauchampRioux<-read.csv("UnprocessedSpectra/project_leaves_combined_BeauchampRioux.csv")
## this is important because otherwise the SG filter step acts up
BeauchampRioux<-BeauchampRioux[BeauchampRioux$wavelength>=350,]

## separate reflectance and transmittance
BeauchampRioux_ref<-BeauchampRioux[BeauchampRioux$reflectance.transmittance=="reflectance",]

## create unique ids for an individual leaf
BeauchampRioux_ref$leaf_id<-apply(BeauchampRioux_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
BeauchampRioux_ref$wvl_id<-apply(BeauchampRioux_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

## remove the spectral overlap region and interpolate over it
dup_ids_ref<-BeauchampRioux_ref$wvl_id[duplicated(BeauchampRioux_ref$wvl_id)]
BeauchampRioux_ref_no_dups<-BeauchampRioux_ref[-which(BeauchampRioux_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(BeauchampRioux_ref_no_dups$wavelength)):floor(max(BeauchampRioux_ref_no_dups$wavelength))

## apply linear interpolation step over spectral overlap
BeauchampRioux_ref_cleaned <-BeauchampRioux_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

## average across leaves for a sample
BeauchampRioux_ref_agg<- BeauchampRioux_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

## apply SG filter across groups of bands
BeauchampRioux_ref_sg_VIS <- BeauchampRioux_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
BeauchampRioux_ref_sg_NIR <- BeauchampRioux_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
BeauchampRioux_ref_sg_SWIR1 <- BeauchampRioux_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
BeauchampRioux_ref_sg_SWIR2 <- BeauchampRioux_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
BeauchampRioux_ref_sg <- bind_rows(BeauchampRioux_ref_sg_VIS,
                        BeauchampRioux_ref_sg_NIR,
                        BeauchampRioux_ref_sg_SWIR1,
                        BeauchampRioux_ref_sg_SWIR2) %>% 
  arrange(sample_id)

BeauchampRioux_ref_sg_wide<-reshape(subset(as.data.frame(BeauchampRioux_ref_sg),select= -value),
                         timevar="wvl",idvar="sample_id",direction="wide")
colnames(BeauchampRioux_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                       replacement = "",
                                       x = colnames(BeauchampRioux_ref_sg_wide))

write.csv(BeauchampRioux_ref_sg_wide,"ProcessedSpectra/BeauchampRioux_ref_processed.csv",row.names=F)

## now the same for transmittance
BeauchampRioux_trans<-BeauchampRioux[BeauchampRioux$reflectance.transmittance=="transmittance",]
BeauchampRioux_trans$leaf_id<-apply(BeauchampRioux_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
BeauchampRioux_trans$wvl_id<-apply(BeauchampRioux_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-BeauchampRioux_trans$wvl_id[duplicated(BeauchampRioux_trans$wvl_id)]
BeauchampRioux_trans_no_dups<-BeauchampRioux_trans[-which(BeauchampRioux_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(BeauchampRioux_trans_no_dups$wavelength)):floor(max(BeauchampRioux_trans_no_dups$wavelength))

BeauchampRioux_trans_cleaned <-BeauchampRioux_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

BeauchampRioux_trans_agg<- BeauchampRioux_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

BeauchampRioux_trans_sg_VIS <- BeauchampRioux_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
BeauchampRioux_trans_sg_NIR <- BeauchampRioux_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
BeauchampRioux_trans_sg_SWIR1 <- BeauchampRioux_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
BeauchampRioux_trans_sg_SWIR2 <- BeauchampRioux_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
BeauchampRioux_trans_sg <- bind_rows(BeauchampRioux_trans_sg_VIS,
                                   BeauchampRioux_trans_sg_NIR,
                                   BeauchampRioux_trans_sg_SWIR1,
                                   BeauchampRioux_trans_sg_SWIR2) %>% 
  arrange(sample_id)

BeauchampRioux_trans_sg_wide<-reshape(subset(as.data.frame(BeauchampRioux_trans_sg),select= -value),
                                    timevar="wvl",idvar="sample_id",direction="wide")
colnames(BeauchampRioux_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                           replacement = "",
                                           x = colnames(BeauchampRioux_trans_sg_wide))

write.csv(BeauchampRioux_trans_sg_wide,"ProcessedSpectra/BeauchampRioux_trans_processed.csv",row.names=F)

###########################################################
## Boucherville 2018 spectra

## now I remove comments -- 
## all annotations can be found in the Beauchamp-Rioux section

Boucherville2018<-read.csv("UnprocessedSpectra/project_leaves_combined_Boucherville2018.csv")
Boucherville2018<-Boucherville2018[Boucherville2018$wavelength>=350,]

Boucherville2018_ref<-Boucherville2018[Boucherville2018$reflectance.transmittance=="reflectance",]

Boucherville2018_ref$leaf_id<-apply(Boucherville2018_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2018_ref$wvl_id<-apply(Boucherville2018_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Boucherville2018_ref$wvl_id[duplicated(Boucherville2018_ref$wvl_id)]
Boucherville2018_ref_no_dups<-Boucherville2018_ref[-which(Boucherville2018_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Boucherville2018_ref_no_dups$wavelength)):floor(max(Boucherville2018_ref_no_dups$wavelength))

Boucherville2018_ref_cleaned <-Boucherville2018_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Boucherville2018_ref_agg<- Boucherville2018_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Boucherville2018_ref_sg_VIS <- Boucherville2018_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Boucherville2018_ref_sg_NIR <- Boucherville2018_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Boucherville2018_ref_sg_SWIR1 <- Boucherville2018_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Boucherville2018_ref_sg_SWIR2 <- Boucherville2018_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Boucherville2018_ref_sg <- bind_rows(Boucherville2018_ref_sg_VIS,
                                   Boucherville2018_ref_sg_NIR,
                                   Boucherville2018_ref_sg_SWIR1,
                                   Boucherville2018_ref_sg_SWIR2) %>% 
  arrange(sample_id)

Boucherville2018_ref_sg_wide<-reshape(subset(as.data.frame(Boucherville2018_ref_sg),select= -value),
                                    timevar="wvl",idvar="sample_id",direction="wide")
colnames(Boucherville2018_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                           replacement = "",
                                           x = colnames(Boucherville2018_ref_sg_wide))

write.csv(Boucherville2018_ref_sg_wide,"ProcessedSpectra/Boucherville2018_ref_processed.csv",row.names=F)

## now the same for transmittance
Boucherville2018_trans<-Boucherville2018[Boucherville2018$reflectance.transmittance=="transmittance",]
Boucherville2018_trans$leaf_id<-apply(Boucherville2018_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2018_trans$wvl_id<-apply(Boucherville2018_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Boucherville2018_trans$wvl_id[duplicated(Boucherville2018_trans$wvl_id)]
Boucherville2018_trans_no_dups<-Boucherville2018_trans[-which(Boucherville2018_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Boucherville2018_trans_no_dups$wavelength)):floor(max(Boucherville2018_trans_no_dups$wavelength))

Boucherville2018_trans_cleaned <-Boucherville2018_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Boucherville2018_trans_agg<- Boucherville2018_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Boucherville2018_trans_sg_VIS <- Boucherville2018_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Boucherville2018_trans_sg_NIR <- Boucherville2018_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Boucherville2018_trans_sg_SWIR1 <- Boucherville2018_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Boucherville2018_trans_sg_SWIR2 <- Boucherville2018_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Boucherville2018_trans_sg <- bind_rows(Boucherville2018_trans_sg_VIS,
                                     Boucherville2018_trans_sg_NIR,
                                     Boucherville2018_trans_sg_SWIR1,
                                     Boucherville2018_trans_sg_SWIR2) %>% 
  arrange(sample_id)

Boucherville2018_trans_sg_wide<-reshape(subset(as.data.frame(Boucherville2018_trans_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Boucherville2018_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(Boucherville2018_trans_sg_wide))

write.csv(Boucherville2018_trans_sg_wide,"ProcessedSpectra/Boucherville2018_trans_processed.csv",row.names=F)

###########################################################
## Girard spectra

Girard<-read.csv("UnprocessedSpectra/project_leaves_combined_Girard.csv")
Girard<-Girard[Girard$wavelength>=350,]

Girard_ref<-Girard[Girard$reflectance.transmittance=="reflectance",]

Girard_ref$leaf_id<-apply(Girard_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Girard_ref$wvl_id<-apply(Girard_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Girard_ref$wvl_id[duplicated(Girard_ref$wvl_id)]
Girard_ref_no_dups<-Girard_ref[-which(Girard_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Girard_ref_no_dups$wavelength)):floor(max(Girard_ref_no_dups$wavelength))

Girard_ref_cleaned <-Girard_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Girard_ref_agg<- Girard_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Girard_ref_sg_VIS <- Girard_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Girard_ref_sg_NIR <- Girard_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Girard_ref_sg_SWIR1 <- Girard_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Girard_ref_sg_SWIR2 <- Girard_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Girard_ref_sg <- bind_rows(Girard_ref_sg_VIS,
                                     Girard_ref_sg_NIR,
                                     Girard_ref_sg_SWIR1,
                                     Girard_ref_sg_SWIR2) %>% 
  arrange(sample_id)

Girard_ref_sg_wide<-reshape(subset(as.data.frame(Girard_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Girard_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(Girard_ref_sg_wide))

write.csv(Girard_ref_sg_wide,"ProcessedSpectra/Girard_ref_processed.csv",row.names=F)

## now the same for transmittance
Girard_trans<-Girard[Girard$reflectance.transmittance=="transmittance",]
Girard_trans$leaf_id<-apply(Girard_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Girard_trans$wvl_id<-apply(Girard_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Girard_trans$wvl_id[duplicated(Girard_trans$wvl_id)]
Girard_trans_no_dups<-Girard_trans[-which(Girard_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Girard_trans_no_dups$wavelength)):floor(max(Girard_trans_no_dups$wavelength))

Girard_trans_cleaned <-Girard_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Girard_trans_agg<- Girard_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Girard_trans_sg_VIS <- Girard_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Girard_trans_sg_NIR <- Girard_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Girard_trans_sg_SWIR1 <- Girard_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Girard_trans_sg_SWIR2 <- Girard_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Girard_trans_sg <- bind_rows(Girard_trans_sg_VIS,
                                       Girard_trans_sg_NIR,
                                       Girard_trans_sg_SWIR1,
                                       Girard_trans_sg_SWIR2) %>% 
  arrange(sample_id)

Girard_trans_sg_wide<-reshape(subset(as.data.frame(Girard_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Girard_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(Girard_trans_sg_wide))

write.csv(Girard_trans_sg_wide,"ProcessedSpectra/Girard_trans_processed.csv",row.names=F)

###########################################################
## Hacker spectra

Hacker2019<-read.csv("UnprocessedSpectra/project_leaves_combined_Hacker2019.csv")
Hacker2019<-Hacker2019[Hacker2019$wavelength>=350,]

Hacker2019_ref<-Hacker2019[Hacker2019$reflectance.transmittance=="reflectance",]

Hacker2019_ref$leaf_id<-apply(Hacker2019_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Hacker2019_ref$wvl_id<-apply(Hacker2019_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Hacker2019_ref$wvl_id[duplicated(Hacker2019_ref$wvl_id)]
Hacker2019_ref_no_dups<-Hacker2019_ref[-which(Hacker2019_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Hacker2019_ref_no_dups$wavelength)):floor(max(Hacker2019_ref_no_dups$wavelength))

Hacker2019_ref_cleaned <-Hacker2019_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Hacker2019_ref_agg<- Hacker2019_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Hacker2019_ref_sg_VIS <- Hacker2019_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Hacker2019_ref_sg_NIR <- Hacker2019_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Hacker2019_ref_sg_SWIR1 <- Hacker2019_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Hacker2019_ref_sg_SWIR2 <- Hacker2019_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Hacker2019_ref_sg <- bind_rows(Hacker2019_ref_sg_VIS,
                                     Hacker2019_ref_sg_NIR,
                                     Hacker2019_ref_sg_SWIR1,
                                     Hacker2019_ref_sg_SWIR2) %>% 
  arrange(sample_id)

Hacker2019_ref_sg_wide<-reshape(subset(as.data.frame(Hacker2019_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Hacker2019_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(Hacker2019_ref_sg_wide))

write.csv(Hacker2019_ref_sg_wide,"ProcessedSpectra/Hacker2019_ref_processed.csv",row.names=F)

## now the same for transmittance
Hacker2019_trans<-Hacker2019[Hacker2019$reflectance.transmittance=="transmittance",]
Hacker2019_trans$leaf_id<-apply(Hacker2019_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Hacker2019_trans$wvl_id<-apply(Hacker2019_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Hacker2019_trans$wvl_id[duplicated(Hacker2019_trans$wvl_id)]
Hacker2019_trans_no_dups<-Hacker2019_trans[-which(Hacker2019_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Hacker2019_trans_no_dups$wavelength)):floor(max(Hacker2019_trans_no_dups$wavelength))

Hacker2019_trans_cleaned <-Hacker2019_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Hacker2019_trans_agg<- Hacker2019_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Hacker2019_trans_sg_VIS <- Hacker2019_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Hacker2019_trans_sg_NIR <- Hacker2019_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Hacker2019_trans_sg_SWIR1 <- Hacker2019_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Hacker2019_trans_sg_SWIR2 <- Hacker2019_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Hacker2019_trans_sg <- bind_rows(Hacker2019_trans_sg_VIS,
                                       Hacker2019_trans_sg_NIR,
                                       Hacker2019_trans_sg_SWIR1,
                                       Hacker2019_trans_sg_SWIR2) %>% 
  arrange(sample_id)

Hacker2019_trans_sg_wide<-reshape(subset(as.data.frame(Hacker2019_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Hacker2019_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(Hacker2019_trans_sg_wide))

write.csv(Hacker2019_trans_sg_wide,"ProcessedSpectra/Hacker2019_trans_processed.csv",row.names=F)

###########################################################
## Blanchard spectra

Blanchard<-read.csv("UnprocessedSpectra/project_leaves_combined_Blanchard.csv")
Blanchard<-Blanchard[Blanchard$wavelength>=350,]

Blanchard_ref<-Blanchard[Blanchard$reflectance.transmittance=="reflectance",]

Blanchard_ref$leaf_id<-apply(Blanchard_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Blanchard_ref$wvl_id<-apply(Blanchard_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Blanchard_ref$wvl_id[duplicated(Blanchard_ref$wvl_id)]
Blanchard_ref_no_dups<-Blanchard_ref[-which(Blanchard_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Blanchard_ref_no_dups$wavelength)):floor(max(Blanchard_ref_no_dups$wavelength))

Blanchard_ref_cleaned <-Blanchard_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Blanchard_ref_agg<- Blanchard_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Blanchard_ref_sg_VIS <- Blanchard_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Blanchard_ref_sg_NIR <- Blanchard_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Blanchard_ref_sg_SWIR1 <- Blanchard_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Blanchard_ref_sg_SWIR2 <- Blanchard_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Blanchard_ref_sg <- bind_rows(Blanchard_ref_sg_VIS,
                                     Blanchard_ref_sg_NIR,
                                     Blanchard_ref_sg_SWIR1,
                                     Blanchard_ref_sg_SWIR2) %>% 
  arrange(sample_id)

Blanchard_ref_sg_wide<-reshape(subset(as.data.frame(Blanchard_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Blanchard_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(Blanchard_ref_sg_wide))

write.csv(Blanchard_ref_sg_wide,"ProcessedSpectra/Blanchard_ref_processed.csv",row.names=F)

## now the same for transmittance
Blanchard_trans<-Blanchard[Blanchard$reflectance.transmittance=="transmittance",]
Blanchard_trans$leaf_id<-apply(Blanchard_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Blanchard_trans$wvl_id<-apply(Blanchard_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Blanchard_trans$wvl_id[duplicated(Blanchard_trans$wvl_id)]
Blanchard_trans_no_dups<-Blanchard_trans[-which(Blanchard_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Blanchard_trans_no_dups$wavelength)):floor(max(Blanchard_trans_no_dups$wavelength))

Blanchard_trans_cleaned <-Blanchard_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Blanchard_trans_agg<- Blanchard_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Blanchard_trans_sg_VIS <- Blanchard_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Blanchard_trans_sg_NIR <- Blanchard_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Blanchard_trans_sg_SWIR1 <- Blanchard_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Blanchard_trans_sg_SWIR2 <- Blanchard_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Blanchard_trans_sg <- bind_rows(Blanchard_trans_sg_VIS,
                                       Blanchard_trans_sg_NIR,
                                       Blanchard_trans_sg_SWIR1,
                                       Blanchard_trans_sg_SWIR2) %>% 
  arrange(sample_id)

Blanchard_trans_sg_wide<-reshape(subset(as.data.frame(Blanchard_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Blanchard_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(Blanchard_trans_sg_wide))

write.csv(Blanchard_trans_sg_wide,"ProcessedSpectra/Blanchard_trans_processed.csv",row.names=F)

###########################################################
## Boucherville 2019 spectra

Boucherville2019<-read.csv("UnprocessedSpectra/project_leaves_combined_Boucherville2019.csv")
## temporary fix -- there are two spectra per leaf with this ID, and one of each is bad
Boucherville2019<-Boucherville2019[-which(Boucherville2019$sample_id=="44245455"),]
Boucherville2019<-Boucherville2019[Boucherville2019$wavelength>=350,]

Boucherville2019_ref<-Boucherville2019[Boucherville2019$reflectance.transmittance=="reflectance",]

Boucherville2019_ref$leaf_id<-apply(Boucherville2019_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2019_ref$wvl_id<-apply(Boucherville2019_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Boucherville2019_ref$wvl_id[duplicated(Boucherville2019_ref$wvl_id)]
Boucherville2019_ref_no_dups<-Boucherville2019_ref[-which(Boucherville2019_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Boucherville2019_ref_no_dups$wavelength)):floor(max(Boucherville2019_ref_no_dups$wavelength))

Boucherville2019_ref_cleaned <-Boucherville2019_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Boucherville2019_ref_agg<- Boucherville2019_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Boucherville2019_ref_sg_VIS <- Boucherville2019_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Boucherville2019_ref_sg_NIR <- Boucherville2019_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Boucherville2019_ref_sg_SWIR1 <- Boucherville2019_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Boucherville2019_ref_sg_SWIR2 <- Boucherville2019_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Boucherville2019_ref_sg <- bind_rows(Boucherville2019_ref_sg_VIS,
                                     Boucherville2019_ref_sg_NIR,
                                     Boucherville2019_ref_sg_SWIR1,
                                     Boucherville2019_ref_sg_SWIR2) %>% 
  arrange(sample_id)

Boucherville2019_ref_sg_wide<-reshape(subset(as.data.frame(Boucherville2019_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Boucherville2019_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(Boucherville2019_ref_sg_wide))

write.csv(Boucherville2019_ref_sg_wide,"ProcessedSpectra/Boucherville2019_ref_processed.csv",row.names=F)

## now the same for transmittance
Boucherville2019_trans<-Boucherville2019[Boucherville2019$reflectance.transmittance=="transmittance",]
Boucherville2019_trans$leaf_id<-apply(Boucherville2019_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2019_trans$wvl_id<-apply(Boucherville2019_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Boucherville2019_trans$wvl_id[duplicated(Boucherville2019_trans$wvl_id)]
Boucherville2019_trans_no_dups<-Boucherville2019_trans[-which(Boucherville2019_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Boucherville2019_trans_no_dups$wavelength)):floor(max(Boucherville2019_trans_no_dups$wavelength))

Boucherville2019_trans_cleaned <-Boucherville2019_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Boucherville2019_trans_agg<- Boucherville2019_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Boucherville2019_trans_sg_VIS <- Boucherville2019_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Boucherville2019_trans_sg_NIR <- Boucherville2019_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Boucherville2019_trans_sg_SWIR1 <- Boucherville2019_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Boucherville2019_trans_sg_SWIR2 <- Boucherville2019_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Boucherville2019_trans_sg <- bind_rows(Boucherville2019_trans_sg_VIS,
                                       Boucherville2019_trans_sg_NIR,
                                       Boucherville2019_trans_sg_SWIR1,
                                       Boucherville2019_trans_sg_SWIR2) %>% 
  arrange(sample_id)

Boucherville2019_trans_sg_wide<-reshape(subset(as.data.frame(Boucherville2019_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Boucherville2019_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(Boucherville2019_trans_sg_wide))

write.csv(Boucherville2019_trans_sg_wide,"ProcessedSpectra/Boucherville2019_trans_processed.csv",row.names=F)

###########################################################
## CABO General 2019 spectra

CABOGeneral2019<-read.csv("UnprocessedSpectra/project_leaves_combined_CABOGeneral2019.csv")
CABOGeneral2019<-CABOGeneral2019[CABOGeneral2019$wavelength>=350,]

CABOGeneral2019_ref<-CABOGeneral2019[CABOGeneral2019$reflectance.transmittance=="reflectance",]

CABOGeneral2019_ref$leaf_id<-apply(CABOGeneral2019_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneral2019_ref$wvl_id<-apply(CABOGeneral2019_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-CABOGeneral2019_ref$wvl_id[duplicated(CABOGeneral2019_ref$wvl_id)]
CABOGeneral2019_ref_no_dups<-CABOGeneral2019_ref[-which(CABOGeneral2019_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(CABOGeneral2019_ref_no_dups$wavelength)):floor(max(CABOGeneral2019_ref_no_dups$wavelength))

CABOGeneral2019_ref_cleaned <-CABOGeneral2019_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

CABOGeneral2019_ref_agg<- CABOGeneral2019_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

CABOGeneral2019_ref_sg_VIS <- CABOGeneral2019_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
CABOGeneral2019_ref_sg_NIR <- CABOGeneral2019_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
CABOGeneral2019_ref_sg_SWIR1 <- CABOGeneral2019_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
CABOGeneral2019_ref_sg_SWIR2 <- CABOGeneral2019_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
CABOGeneral2019_ref_sg <- bind_rows(CABOGeneral2019_ref_sg_VIS,
                                     CABOGeneral2019_ref_sg_NIR,
                                     CABOGeneral2019_ref_sg_SWIR1,
                                     CABOGeneral2019_ref_sg_SWIR2) %>% 
  arrange(sample_id)

CABOGeneral2019_ref_sg_wide<-reshape(subset(as.data.frame(CABOGeneral2019_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(CABOGeneral2019_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(CABOGeneral2019_ref_sg_wide))

write.csv(CABOGeneral2019_ref_sg_wide,"ProcessedSpectra/CABOGeneral2019_ref_processed.csv",row.names=F)

## now the same for transmittance
CABOGeneral2019_trans<-CABOGeneral2019[CABOGeneral2019$reflectance.transmittance=="transmittance",]
CABOGeneral2019_trans$leaf_id<-apply(CABOGeneral2019_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneral2019_trans$wvl_id<-apply(CABOGeneral2019_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-CABOGeneral2019_trans$wvl_id[duplicated(CABOGeneral2019_trans$wvl_id)]
CABOGeneral2019_trans_no_dups<-CABOGeneral2019_trans[-which(CABOGeneral2019_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(CABOGeneral2019_trans_no_dups$wavelength)):floor(max(CABOGeneral2019_trans_no_dups$wavelength))

CABOGeneral2019_trans_cleaned <-CABOGeneral2019_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

CABOGeneral2019_trans_agg<- CABOGeneral2019_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

CABOGeneral2019_trans_sg_VIS <- CABOGeneral2019_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
CABOGeneral2019_trans_sg_NIR <- CABOGeneral2019_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
CABOGeneral2019_trans_sg_SWIR1 <- CABOGeneral2019_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
CABOGeneral2019_trans_sg_SWIR2 <- CABOGeneral2019_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
CABOGeneral2019_trans_sg <- bind_rows(CABOGeneral2019_trans_sg_VIS,
                                       CABOGeneral2019_trans_sg_NIR,
                                       CABOGeneral2019_trans_sg_SWIR1,
                                       CABOGeneral2019_trans_sg_SWIR2) %>% 
  arrange(sample_id)

CABOGeneral2019_trans_sg_wide<-reshape(subset(as.data.frame(CABOGeneral2019_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(CABOGeneral2019_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(CABOGeneral2019_trans_sg_wide))

write.csv(CABOGeneral2019_trans_sg_wide,"ProcessedSpectra/CABOGeneral2019_trans_processed.csv",row.names=F)

###########################################################
## Crofts spectra

Crofts<-read.csv("UnprocessedSpectra/project_leaves_combined_Crofts.csv")
Crofts<-Crofts[Crofts$wavelength>=350,]

Crofts_ref<-Crofts[Crofts$reflectance.transmittance=="reflectance",]

Crofts_ref$leaf_id<-apply(Crofts_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Crofts_ref$wvl_id<-apply(Crofts_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Crofts_ref$wvl_id[duplicated(Crofts_ref$wvl_id)]
Crofts_ref_no_dups<-Crofts_ref[-which(Crofts_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Crofts_ref_no_dups$wavelength)):floor(max(Crofts_ref_no_dups$wavelength))

Crofts_ref_cleaned <-Crofts_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Crofts_ref_agg<- Crofts_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Crofts_ref_sg_VIS <- Crofts_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Crofts_ref_sg_NIR <- Crofts_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Crofts_ref_sg_SWIR1 <- Crofts_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Crofts_ref_sg_SWIR2 <- Crofts_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Crofts_ref_sg <- bind_rows(Crofts_ref_sg_VIS,
                                     Crofts_ref_sg_NIR,
                                     Crofts_ref_sg_SWIR1,
                                     Crofts_ref_sg_SWIR2) %>% 
  arrange(sample_id)

Crofts_ref_sg_wide<-reshape(subset(as.data.frame(Crofts_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Crofts_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(Crofts_ref_sg_wide))

write.csv(Crofts_ref_sg_wide,"ProcessedSpectra/Crofts_ref_processed.csv",row.names=F)

## now the same for transmittance
Crofts_trans<-Crofts[Crofts$reflectance.transmittance=="transmittance",]
Crofts_trans$leaf_id<-apply(Crofts_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Crofts_trans$wvl_id<-apply(Crofts_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Crofts_trans$wvl_id[duplicated(Crofts_trans$wvl_id)]
Crofts_trans_no_dups<-Crofts_trans[-which(Crofts_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Crofts_trans_no_dups$wavelength)):floor(max(Crofts_trans_no_dups$wavelength))

Crofts_trans_cleaned <-Crofts_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Crofts_trans_agg<- Crofts_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Crofts_trans_sg_VIS <- Crofts_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Crofts_trans_sg_NIR <- Crofts_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Crofts_trans_sg_SWIR1 <- Crofts_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Crofts_trans_sg_SWIR2 <- Crofts_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Crofts_trans_sg <- bind_rows(Crofts_trans_sg_VIS,
                                       Crofts_trans_sg_NIR,
                                       Crofts_trans_sg_SWIR1,
                                       Crofts_trans_sg_SWIR2) %>% 
  arrange(sample_id)

Crofts_trans_sg_wide<-reshape(subset(as.data.frame(Crofts_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Crofts_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(Crofts_trans_sg_wide))

write.csv(Crofts_trans_sg_wide,"ProcessedSpectra/Crofts_trans_processed.csv",row.names=F)

###########################################################
## Pardo spectra

Pardo<-read.csv("UnprocessedSpectra/project_leaves_combined_Pardo.csv")
Pardo<-Pardo[Pardo$wavelength>=350,]

Pardo_ref<-Pardo[Pardo$reflectance.transmittance=="reflectance",]

Pardo_ref$leaf_id<-apply(Pardo_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Pardo_ref$wvl_id<-apply(Pardo_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Pardo_ref$wvl_id[duplicated(Pardo_ref$wvl_id)]
Pardo_ref_no_dups<-Pardo_ref[-which(Pardo_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Pardo_ref_no_dups$wavelength)):floor(max(Pardo_ref_no_dups$wavelength))

Pardo_ref_cleaned <-Pardo_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Pardo_ref_agg<- Pardo_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Pardo_ref_sg_VIS <- Pardo_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Pardo_ref_sg_NIR <- Pardo_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Pardo_ref_sg_SWIR1 <- Pardo_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Pardo_ref_sg_SWIR2 <- Pardo_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Pardo_ref_sg <- bind_rows(Pardo_ref_sg_VIS,
                                     Pardo_ref_sg_NIR,
                                     Pardo_ref_sg_SWIR1,
                                     Pardo_ref_sg_SWIR2) %>% 
  arrange(sample_id)

Pardo_ref_sg_wide<-reshape(subset(as.data.frame(Pardo_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Pardo_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(Pardo_ref_sg_wide))

write.csv(Pardo_ref_sg_wide,"ProcessedSpectra/Pardo_ref_processed.csv",row.names=F)

## now the same for transmittance
Pardo_trans<-Pardo[Pardo$reflectance.transmittance=="transmittance",]
Pardo_trans$leaf_id<-apply(Pardo_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Pardo_trans$wvl_id<-apply(Pardo_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Pardo_trans$wvl_id[duplicated(Pardo_trans$wvl_id)]
Pardo_trans_no_dups<-Pardo_trans[-which(Pardo_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Pardo_trans_no_dups$wavelength)):floor(max(Pardo_trans_no_dups$wavelength))

Pardo_trans_cleaned <-Pardo_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Pardo_trans_agg<- Pardo_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Pardo_trans_sg_VIS <- Pardo_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Pardo_trans_sg_NIR <- Pardo_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Pardo_trans_sg_SWIR1 <- Pardo_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Pardo_trans_sg_SWIR2 <- Pardo_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Pardo_trans_sg <- bind_rows(Pardo_trans_sg_VIS,
                                       Pardo_trans_sg_NIR,
                                       Pardo_trans_sg_SWIR1,
                                       Pardo_trans_sg_SWIR2) %>% 
  arrange(sample_id)

Pardo_trans_sg_wide<-reshape(subset(as.data.frame(Pardo_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Pardo_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(Pardo_trans_sg_wide))

write.csv(Pardo_trans_sg_wide,"ProcessedSpectra/Pardo_trans_processed.csv",row.names=F)

###########################################################
## Phragmites temporal spectra

PhragmitesTemporal<-read.csv("UnprocessedSpectra/project_leaves_combined_PhragmitesTemporal.csv")
PhragmitesTemporal<-PhragmitesTemporal[PhragmitesTemporal$wavelength>=350,]

PhragmitesTemporal_ref<-PhragmitesTemporal[PhragmitesTemporal$reflectance.transmittance=="reflectance",]

PhragmitesTemporal_ref$leaf_id<-apply(PhragmitesTemporal_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
PhragmitesTemporal_ref$wvl_id<-apply(PhragmitesTemporal_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-PhragmitesTemporal_ref$wvl_id[duplicated(PhragmitesTemporal_ref$wvl_id)]
PhragmitesTemporal_ref_no_dups<-PhragmitesTemporal_ref[-which(PhragmitesTemporal_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(PhragmitesTemporal_ref_no_dups$wavelength)):floor(max(PhragmitesTemporal_ref_no_dups$wavelength))

PhragmitesTemporal_ref_cleaned <-PhragmitesTemporal_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

PhragmitesTemporal_ref_agg<- PhragmitesTemporal_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

PhragmitesTemporal_ref_sg_VIS <- PhragmitesTemporal_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
PhragmitesTemporal_ref_sg_NIR <- PhragmitesTemporal_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
PhragmitesTemporal_ref_sg_SWIR1 <- PhragmitesTemporal_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
PhragmitesTemporal_ref_sg_SWIR2 <- PhragmitesTemporal_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
PhragmitesTemporal_ref_sg <- bind_rows(PhragmitesTemporal_ref_sg_VIS,
                                     PhragmitesTemporal_ref_sg_NIR,
                                     PhragmitesTemporal_ref_sg_SWIR1,
                                     PhragmitesTemporal_ref_sg_SWIR2) %>% 
  arrange(sample_id)

PhragmitesTemporal_ref_sg_wide<-reshape(subset(as.data.frame(PhragmitesTemporal_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(PhragmitesTemporal_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(PhragmitesTemporal_ref_sg_wide))

write.csv(PhragmitesTemporal_ref_sg_wide,"ProcessedSpectra/PhragmitesTemporal_ref_processed.csv",row.names=F)

## now the same for transmittance
PhragmitesTemporal_trans<-PhragmitesTemporal[PhragmitesTemporal$reflectance.transmittance=="transmittance",]
PhragmitesTemporal_trans$leaf_id<-apply(PhragmitesTemporal_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
PhragmitesTemporal_trans$wvl_id<-apply(PhragmitesTemporal_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-PhragmitesTemporal_trans$wvl_id[duplicated(PhragmitesTemporal_trans$wvl_id)]
PhragmitesTemporal_trans_no_dups<-PhragmitesTemporal_trans[-which(PhragmitesTemporal_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(PhragmitesTemporal_trans_no_dups$wavelength)):floor(max(PhragmitesTemporal_trans_no_dups$wavelength))

PhragmitesTemporal_trans_cleaned <-PhragmitesTemporal_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

PhragmitesTemporal_trans_agg<- PhragmitesTemporal_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

PhragmitesTemporal_trans_sg_VIS <- PhragmitesTemporal_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
PhragmitesTemporal_trans_sg_NIR <- PhragmitesTemporal_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
PhragmitesTemporal_trans_sg_SWIR1 <- PhragmitesTemporal_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
PhragmitesTemporal_trans_sg_SWIR2 <- PhragmitesTemporal_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
PhragmitesTemporal_trans_sg <- bind_rows(PhragmitesTemporal_trans_sg_VIS,
                                       PhragmitesTemporal_trans_sg_NIR,
                                       PhragmitesTemporal_trans_sg_SWIR1,
                                       PhragmitesTemporal_trans_sg_SWIR2) %>% 
  arrange(sample_id)

PhragmitesTemporal_trans_sg_wide<-reshape(subset(as.data.frame(PhragmitesTemporal_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(PhragmitesTemporal_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(PhragmitesTemporal_trans_sg_wide))

write.csv(PhragmitesTemporal_trans_sg_wide,"ProcessedSpectra/PhragmitesTemporal_trans_processed.csv",row.names=F)

###########################################################
## CABO General Other spectra

CABOGeneralOther<-read.csv("UnprocessedSpectra/project_leaves_combined_CABOGeneralOther.csv")
CABOGeneralOther<-CABOGeneralOther[CABOGeneralOther$wavelength>=350,]

CABOGeneralOther_ref<-CABOGeneralOther[CABOGeneralOther$reflectance.transmittance=="reflectance",]

CABOGeneralOther_ref$leaf_id<-apply(CABOGeneralOther_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneralOther_ref$wvl_id<-apply(CABOGeneralOther_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-CABOGeneralOther_ref$wvl_id[duplicated(CABOGeneralOther_ref$wvl_id)]
CABOGeneralOther_ref_no_dups<-CABOGeneralOther_ref[-which(CABOGeneralOther_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(CABOGeneralOther_ref_no_dups$wavelength)):floor(max(CABOGeneralOther_ref_no_dups$wavelength))

CABOGeneralOther_ref_cleaned <-CABOGeneralOther_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

CABOGeneralOther_ref_agg<- CABOGeneralOther_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

CABOGeneralOther_ref_sg_VIS <- CABOGeneralOther_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
CABOGeneralOther_ref_sg_NIR <- CABOGeneralOther_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
CABOGeneralOther_ref_sg_SWIR1 <- CABOGeneralOther_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
CABOGeneralOther_ref_sg_SWIR2 <- CABOGeneralOther_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
CABOGeneralOther_ref_sg <- bind_rows(CABOGeneralOther_ref_sg_VIS,
                                     CABOGeneralOther_ref_sg_NIR,
                                     CABOGeneralOther_ref_sg_SWIR1,
                                     CABOGeneralOther_ref_sg_SWIR2) %>% 
  arrange(sample_id)

CABOGeneralOther_ref_sg_wide<-reshape(subset(as.data.frame(CABOGeneralOther_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(CABOGeneralOther_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(CABOGeneralOther_ref_sg_wide))

write.csv(CABOGeneralOther_ref_sg_wide,"ProcessedSpectra/CABOGeneralOther_ref_processed.csv",row.names=F)

## now the same for transmittance
CABOGeneralOther_trans<-CABOGeneralOther[CABOGeneralOther$reflectance.transmittance=="transmittance",]
CABOGeneralOther_trans$leaf_id<-apply(CABOGeneralOther_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneralOther_trans$wvl_id<-apply(CABOGeneralOther_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-CABOGeneralOther_trans$wvl_id[duplicated(CABOGeneralOther_trans$wvl_id)]
CABOGeneralOther_trans_no_dups<-CABOGeneralOther_trans[-which(CABOGeneralOther_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(CABOGeneralOther_trans_no_dups$wavelength)):floor(max(CABOGeneralOther_trans_no_dups$wavelength))

CABOGeneralOther_trans_cleaned <-CABOGeneralOther_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

CABOGeneralOther_trans_agg<- CABOGeneralOther_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

CABOGeneralOther_trans_sg_VIS <- CABOGeneralOther_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
CABOGeneralOther_trans_sg_NIR <- CABOGeneralOther_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
CABOGeneralOther_trans_sg_SWIR1 <- CABOGeneralOther_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
CABOGeneralOther_trans_sg_SWIR2 <- CABOGeneralOther_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
CABOGeneralOther_trans_sg <- bind_rows(CABOGeneralOther_trans_sg_VIS,
                                       CABOGeneralOther_trans_sg_NIR,
                                       CABOGeneralOther_trans_sg_SWIR1,
                                       CABOGeneralOther_trans_sg_SWIR2) %>% 
  arrange(sample_id)

CABOGeneralOther_trans_sg_wide<-reshape(subset(as.data.frame(CABOGeneralOther_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(CABOGeneralOther_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(CABOGeneralOther_trans_sg_wide))

write.csv(CABOGeneralOther_trans_sg_wide,"ProcessedSpectra/CABOGeneralOther_trans_processed.csv",row.names=F)

###########################################################
## Warren spectra

Warren<-read.csv("UnprocessedSpectra/project_leaves_combined_Warren.csv")
Warren<-Warren[Warren$wavelength>=350,]

Warren_ref<-Warren[Warren$reflectance.transmittance=="reflectance",]

Warren_ref$leaf_id<-apply(Warren_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Warren_ref$wvl_id<-apply(Warren_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Warren_ref$wvl_id[duplicated(Warren_ref$wvl_id)]
Warren_ref_no_dups<-Warren_ref[-which(Warren_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Warren_ref_no_dups$wavelength)):floor(max(Warren_ref_no_dups$wavelength))

Warren_ref_cleaned <-Warren_ref_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Warren_ref_agg<- Warren_ref_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Warren_ref_sg_VIS <- Warren_ref_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Warren_ref_sg_NIR <- Warren_ref_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Warren_ref_sg_SWIR1 <- Warren_ref_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Warren_ref_sg_SWIR2 <- Warren_ref_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Warren_ref_sg <- bind_rows(Warren_ref_sg_VIS,
                                     Warren_ref_sg_NIR,
                                     Warren_ref_sg_SWIR1,
                                     Warren_ref_sg_SWIR2) %>% 
  arrange(sample_id)

Warren_ref_sg_wide<-reshape(subset(as.data.frame(Warren_ref_sg),select= -value),
                                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Warren_ref_sg_wide)<-gsub(pattern = "value_sg.",
                                             replacement = "",
                                             x = colnames(Warren_ref_sg_wide))

write.csv(Warren_ref_sg_wide,"ProcessedSpectra/Warren_ref_processed.csv",row.names=F)

## now the same for transmittance
Warren_trans<-Warren[Warren$reflectance.transmittance=="transmittance",]
Warren_trans$leaf_id<-apply(Warren_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Warren_trans$wvl_id<-apply(Warren_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Warren_trans$wvl_id[duplicated(Warren_trans$wvl_id)]
Warren_trans_no_dups<-Warren_trans[-which(Warren_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Warren_trans_no_dups$wavelength)):floor(max(Warren_trans_no_dups$wavelength))

Warren_trans_cleaned <-Warren_trans_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Warren_trans_agg<- Warren_trans_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Warren_trans_sg_VIS <- Warren_trans_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Warren_trans_sg_NIR <- Warren_trans_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Warren_trans_sg_SWIR1 <- Warren_trans_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Warren_trans_sg_SWIR2 <- Warren_trans_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Warren_trans_sg <- bind_rows(Warren_trans_sg_VIS,
                                       Warren_trans_sg_NIR,
                                       Warren_trans_sg_SWIR1,
                                       Warren_trans_sg_SWIR2) %>% 
  arrange(sample_id)

Warren_trans_sg_wide<-reshape(subset(as.data.frame(Warren_trans_sg),select= -value),
                                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Warren_trans_sg_wide)<-gsub(pattern = "value_sg.",
                                               replacement = "",
                                               x = colnames(Warren_trans_sg_wide))

write.csv(Warren_trans_sg_wide,"ProcessedSpectra/Warren_trans_processed.csv",row.names=F)