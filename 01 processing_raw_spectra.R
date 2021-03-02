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
## With a few exceptions -- like the Dessain project, which
## used an ASD spectrometer -- there's a lot of copy-pasting
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

##########################################################
## Dessain
## These were done with an ASD and integrating sphere
## (the other projects with an SVC and integrating sphere)
## the spectral overlap region is already spliced out
## but I do need to do dark current correction, which isn't
## true with the other projects

zenith <- read.table('UnprocessedSpectra/DessainSpectra/c7101904_specchioformat.txt',
                     header = T)
zenith.sub <-tbl_df(dplyr::filter(zenith, wvl >= 400) )

### get paths of .asd files-----
asd_dirs <- list.dirs('UnprocessedSpectra/DessainSpectra/spectra-20201002T221442Z-001/spectra/',
                      full.names = T)
asd_paths <- list.files(asd_dirs, full.names = T, pattern = "\\.asd$")

# add filename
asd_files <- list.files(asd_dirs, full.names = F, pattern = "\\.asd$")
asd_files_short <- str_sub(asd_files, end = -5)

#### read all .asd files
dn <- get_spectra(asd_paths, type = 'radiance') # get radiances for target
wr <- get_spectra(asd_paths, type = 'white_reference') # get radiances for white reference

# only keep wvl's between 400 and 2500 nm
dn.sub <- dn[, colnames(dn) %in% 400:2500]
wr.sub <- wr[, colnames(wr) %in% 400:2500]

# convert spectra to tbl data frame
dn.df <- tbl_df(dn.sub)
wr.df <- tbl_df(wr.sub)
dn.df$fileName <- asd_files_short
wr.df$fileName <- asd_files_short

# convert to long format and order by fileName
dn.long <- gather(dn.df, wvl.char, target.rad, -fileName) %>%
  mutate(wvl = as.numeric(wvl.char)) %>%
  arrange(fileName, wvl) %>%
  dplyr::select(-wvl.char)

wr.long <- gather(wr.df, wvl.char, wr.rad, -fileName) %>%
  mutate(wvl = as.numeric(wvl.char)) %>%
  arrange(fileName, wvl) %>%
  dplyr::select(-wvl.char)

# combine the two
comb <- inner_join(dn.long, wr.long, by = c('fileName', 'wvl')) 

# check that all wavelengths in spectra and zenith calibration are identical
wvl.vec2 <- subset(dn.long, fileName == "2017-05-17-jbmcb00001")$wvl
all(zenith.sub$wvl == wvl.vec2) # TRUE?

## read in metadata
records<-read.csv("UnprocessedSpectra/DessainSpectra/spectrumRecord.csv")
measurements<-read.csv("UnprocessedSpectra/DessainSpectra/spectrumMeasurement.csv")

# merge with spectrum info
comb2 <- inner_join(comb, records, by = 'fileName')
comb3 <- inner_join(comb2, measurements, by = 'eventID') %>%
  arrange(fileName, wvl) %>% # sort by file name, specType and wavelength
  dplyr::filter(qualityFlag == 'good',
                instrumentationID == "ASD RTS-3ZC Integrating Sphere #6045-11",
                targetType != 'Reference') # only keep the good dn and drop the References

# function to correct for stray light and measure absolute reflectance for each eventID
refl_corr <- function(x) {
  stray <- zenith.sub %>% transmute(wvl, stray.rad = 0)
  if (any(x$targetType == 'Stray light') ) {
    stray <- dplyr::filter(x, targetType == 'Stray light') %>%
      group_by(wvl) %>%
      summarise(stray.rad= mean(target.rad))
  }
  x2 <- inner_join(x, stray, by = 'wvl') %>%
    inner_join(zenith.sub, by = 'wvl') %>%
    dplyr::filter(targetType != 'Stray light') %>%
    mutate(refl_corr = ( (target.rad - stray.rad ) / (wr.rad - stray.rad) ) *  c7101904)
  return(x2)
}

# run function on each folderName (i.e. combination of day / location)
refl_corr.df <- comb3 %>%
  group_by(folderName) %>%
  do(refl_corr(.))

## aggregate reps within a leaf
Dessain_spec <- refl_corr.df %>%
  group_by(parentEventID,leafNumber, wvl) %>%
  summarise(mean.refl.corr = mean(refl_corr)) %>%
  droplevels()

## for consistency with other data sets at this stage
colnames(Dessain_spec)<-c("sample_id","leaf_number","wvl","DN")

Dessain_agg<- Dessain_spec %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Dessain_sg_VIS <- Dessain_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Dessain_sg_NIR <- Dessain_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Dessain_sg_SWIR1 <- Dessain_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Dessain_sg_SWIR2 <- Dessain_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Dessain_sg <- bind_rows(Dessain_sg_VIS,
                        Dessain_sg_NIR,
                        Dessain_sg_SWIR1,
                        Dessain_sg_SWIR2) %>% 
  arrange(sample_id)

Dessain_sg_wide<-reshape(subset(as.data.frame(Dessain_sg),select= -value),
                      timevar="wvl",idvar="sample_id",direction="wide")
colnames(Dessain_sg_wide)<-gsub(pattern = "value_sg.",
                                replacement = "",
                                x = colnames(Dessain_sg_wide))

write.csv(Dessain_sg_wide,"ProcessedSpectra/Dessain_spec_processed.csv",row.names=F)

###########################################################
## Beauchamp-Rioux spectra

## like the rest, these were measured with an SVC
## and integrating sphere

BeauchampRioux_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_BeauchampRioux.csv")

## discard transmittance, keep only 400-2500 nm
BeauchampRioux_spec<-BeauchampRioux_spec[BeauchampRioux_spec$reflectance.transmittance=="reflectance",]
BeauchampRioux_spec<-BeauchampRioux_spec[BeauchampRioux_spec$wavelength>=400 & BeauchampRioux_spec$wavelength<=2500,]

## create unique ids for an individual leaf
BeauchampRioux_spec$leaf_id<-apply(BeauchampRioux_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
BeauchampRioux_spec$wvl_id<-apply(BeauchampRioux_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

## remove the spectral overlap region and interpolate over it
dup_ids<-BeauchampRioux_spec$wvl_id[duplicated(BeauchampRioux_spec$wvl_id)]
BeauchampRioux_spec_no_dups<-BeauchampRioux_spec[-which(BeauchampRioux_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(BeauchampRioux_spec_no_dups$wavelength)):floor(max(BeauchampRioux_spec_no_dups$wavelength))

## apply linear interpolation step over spectral overlap
BeauchampRioux_cleaned <-BeauchampRioux_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

## average across leaves for a sample
BeauchampRioux_agg<- BeauchampRioux_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

## apply SG filter across groups of bands
BeauchampRioux_sg_VIS <- BeauchampRioux_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
BeauchampRioux_sg_NIR <- BeauchampRioux_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
BeauchampRioux_sg_SWIR1 <- BeauchampRioux_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
BeauchampRioux_sg_SWIR2 <- BeauchampRioux_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
BeauchampRioux_sg <- bind_rows(BeauchampRioux_sg_VIS,
                        BeauchampRioux_sg_NIR,
                        BeauchampRioux_sg_SWIR1,
                        BeauchampRioux_sg_SWIR2) %>% 
  arrange(sample_id)

BeauchampRioux_sg_wide<-reshape(subset(as.data.frame(BeauchampRioux_sg),select= -value),
                         timevar="wvl",idvar="sample_id",direction="wide")
colnames(BeauchampRioux_sg_wide)<-gsub(pattern = "value_sg.",
                                       replacement = "",
                                       x = colnames(BeauchampRioux_sg_wide))

write.csv(BeauchampRioux_sg_wide,"ProcessedSpectra/BeauchampRioux_spec_processed.csv",row.names=F)

###########################################################
## Boucherville 2018 spectra

## now I remove comments -- 
## all annotations can be found in the Beauchamp-Rioux section

Boucherville2018_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Boucherville2018.csv")

Boucherville2018_spec<-Boucherville2018_spec[Boucherville2018_spec$reflectance.transmittance=="reflectance",]
Boucherville2018_spec<-Boucherville2018_spec[Boucherville2018_spec$wavelength>=400 & Boucherville2018_spec$wavelength<=2500,]

Boucherville2018_spec$leaf_id<-apply(Boucherville2018_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2018_spec$wvl_id<-apply(Boucherville2018_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-Boucherville2018_spec$wvl_id[duplicated(Boucherville2018_spec$wvl_id)]
Boucherville2018_spec_no_dups<-Boucherville2018_spec[-which(Boucherville2018_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(Boucherville2018_spec_no_dups$wavelength)):floor(max(Boucherville2018_spec_no_dups$wavelength))

Boucherville2018_cleaned <-Boucherville2018_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Boucherville2018_agg<- Boucherville2018_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Boucherville2018_sg_VIS <- Boucherville2018_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Boucherville2018_sg_NIR <- Boucherville2018_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Boucherville2018_sg_SWIR1 <- Boucherville2018_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Boucherville2018_sg_SWIR2 <- Boucherville2018_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Boucherville2018_sg <- bind_rows(Boucherville2018_sg_VIS,
                               Boucherville2018_sg_NIR,
                               Boucherville2018_sg_SWIR1,
                               Boucherville2018_sg_SWIR2) %>% 
  arrange(sample_id)

Boucherville2018_sg_wide<-reshape(subset(as.data.frame(Boucherville2018_sg),select= -value),
                                timevar="wvl",idvar="sample_id",direction="wide")
colnames(Boucherville2018_sg_wide)<-gsub(pattern = "value_sg.",
                                         replacement = "",
                                         x = colnames(Boucherville2018_sg_wide))

write.csv(Boucherville2018_sg_wide,"ProcessedSpectra/Boucherville2018_spec_processed.csv",row.names=F)

###########################################################
## Girard spectra

Girard_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Girard.csv")

Girard_spec<-Girard_spec[Girard_spec$reflectance.transmittance=="reflectance",]
Girard_spec<-Girard_spec[Girard_spec$wavelength>=400 & Girard_spec$wavelength<=2500,]

Girard_spec$leaf_id<-apply(Girard_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
Girard_spec$wvl_id<-apply(Girard_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-Girard_spec$wvl_id[duplicated(Girard_spec$wvl_id)]
Girard_spec_no_dups<-Girard_spec[-which(Girard_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(Girard_spec_no_dups$wavelength)):floor(max(Girard_spec_no_dups$wavelength))

Girard_cleaned <-Girard_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Girard_agg<- Girard_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Girard_sg_VIS <- Girard_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Girard_sg_NIR <- Girard_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Girard_sg_SWIR1 <- Girard_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Girard_sg_SWIR2 <- Girard_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Girard_sg <- bind_rows(Girard_sg_VIS,
                                 Girard_sg_NIR,
                                 Girard_sg_SWIR1,
                                 Girard_sg_SWIR2) %>% 
  arrange(sample_id)

Girard_sg_wide<-reshape(subset(as.data.frame(Girard_sg),select= -value),
                                  timevar="wvl",idvar="sample_id",direction="wide")
colnames(Girard_sg_wide)<-gsub(pattern = "value_sg.",
                                         replacement = "",
                                         x = colnames(Girard_sg_wide))

write.csv(Girard_sg_wide,"ProcessedSpectra/Girard_spec_processed.csv",row.names=F)

###########################################################
## Hacker spectra

Hacker_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Hacker.csv")

Hacker_spec<-Hacker_spec[Hacker_spec$reflectance.transmittance=="reflectance",]
Hacker_spec<-Hacker_spec[Hacker_spec$wavelength>=400 & Hacker_spec$wavelength<=2500,]

Hacker_spec$leaf_id<-apply(Hacker_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
Hacker_spec$wvl_id<-apply(Hacker_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-Hacker_spec$wvl_id[duplicated(Hacker_spec$wvl_id)]
Hacker_spec_no_dups<-Hacker_spec[-which(Hacker_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(Hacker_spec_no_dups$wavelength)):floor(max(Hacker_spec_no_dups$wavelength))

Hacker_cleaned <-Hacker_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Hacker_agg<- Hacker_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Hacker_sg_VIS <- Hacker_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Hacker_sg_NIR <- Hacker_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Hacker_sg_SWIR1 <- Hacker_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Hacker_sg_SWIR2 <- Hacker_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Hacker_sg <- bind_rows(Hacker_sg_VIS,
                       Hacker_sg_NIR,
                       Hacker_sg_SWIR1,
                       Hacker_sg_SWIR2) %>% 
  arrange(sample_id)

Hacker_sg_wide<-reshape(subset(as.data.frame(Hacker_sg),select= -value),
                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Hacker_sg_wide)<-gsub(pattern = "value_sg.",
                               replacement = "",
                               x = colnames(Hacker_sg_wide))

write.csv(Hacker_sg_wide,"ProcessedSpectra/Hacker_spec_processed.csv",row.names=F)

###########################################################
## Blanchard spectra

Blanchard_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Blanchard.csv")

Blanchard_spec<-Blanchard_spec[Blanchard_spec$reflectance.transmittance=="reflectance",]
Blanchard_spec<-Blanchard_spec[Blanchard_spec$wavelength>=400 & Blanchard_spec$wavelength<=2500,]

Blanchard_spec$leaf_id<-apply(Blanchard_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
Blanchard_spec$wvl_id<-apply(Blanchard_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-Blanchard_spec$wvl_id[duplicated(Blanchard_spec$wvl_id)]
Blanchard_spec_no_dups<-Blanchard_spec[-which(Blanchard_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(Blanchard_spec_no_dups$wavelength)):floor(max(Blanchard_spec_no_dups$wavelength))

Blanchard_cleaned <-Blanchard_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Blanchard_agg<- Blanchard_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Blanchard_sg_VIS <- Blanchard_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Blanchard_sg_NIR <- Blanchard_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Blanchard_sg_SWIR1 <- Blanchard_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Blanchard_sg_SWIR2 <- Blanchard_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Blanchard_sg <- bind_rows(Blanchard_sg_VIS,
                       Blanchard_sg_NIR,
                       Blanchard_sg_SWIR1,
                       Blanchard_sg_SWIR2) %>% 
  arrange(sample_id)

Blanchard_sg_wide<-reshape(subset(as.data.frame(Blanchard_sg),select= -value),
                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(Blanchard_sg_wide)<-gsub(pattern = "value_sg.",
                               replacement = "",
                               x = colnames(Blanchard_sg_wide))

write.csv(Blanchard_sg_wide,"ProcessedSpectra/Blanchard_spec_processed.csv",row.names=F)

###########################################################
## Boucherville 2019 spectra

Boucherville2019_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Boucherville2019.csv")
## temporary fix -- there are two spectra per leaf with this ID, and one of each is bad
Boucherville2019_spec<-Boucherville2019_spec[-which(Boucherville2019_spec$sample_id=="44245455"),]

Boucherville2019_spec<-Boucherville2019_spec[Boucherville2019_spec$reflectance.transmittance=="reflectance",]
Boucherville2019_spec<-Boucherville2019_spec[Boucherville2019_spec$wavelength>=400 & Boucherville2019_spec$wavelength<=2500,]

Boucherville2019_spec$leaf_id<-apply(Boucherville2019_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2019_spec$wvl_id<-apply(Boucherville2019_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-Boucherville2019_spec$wvl_id[duplicated(Boucherville2019_spec$wvl_id)]
Boucherville2019_spec_no_dups<-Boucherville2019_spec[-which(Boucherville2019_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(Boucherville2019_spec_no_dups$wavelength)):floor(max(Boucherville2019_spec_no_dups$wavelength))

Boucherville2019_cleaned <-Boucherville2019_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Boucherville2019_agg<- Boucherville2019_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Boucherville2019_sg_VIS <- Boucherville2019_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Boucherville2019_sg_NIR <- Boucherville2019_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Boucherville2019_sg_SWIR1 <- Boucherville2019_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Boucherville2019_sg_SWIR2 <- Boucherville2019_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Boucherville2019_sg <- bind_rows(Boucherville2019_sg_VIS,
                          Boucherville2019_sg_NIR,
                          Boucherville2019_sg_SWIR1,
                          Boucherville2019_sg_SWIR2) %>% 
  arrange(sample_id)

Boucherville2019_sg_wide<-reshape(subset(as.data.frame(Boucherville2019_sg),select= -value),
                           timevar="wvl",idvar="sample_id",direction="wide")
colnames(Boucherville2019_sg_wide)<-gsub(pattern = "value_sg.",
                                  replacement = "",
                                  x = colnames(Boucherville2019_sg_wide))

write.csv(Boucherville2019_sg_wide,"ProcessedSpectra/Boucherville2019_spec_processed.csv",row.names=F)

###########################################################
## CABO General 2019 spectra

CABOGeneral2019_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_CABOGeneral2019.csv")

CABOGeneral2019_spec<-CABOGeneral2019_spec[CABOGeneral2019_spec$reflectance.transmittance=="reflectance",]
CABOGeneral2019_spec<-CABOGeneral2019_spec[CABOGeneral2019_spec$wavelength>=400 & CABOGeneral2019_spec$wavelength<=2500,]

CABOGeneral2019_spec$leaf_id<-apply(CABOGeneral2019_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneral2019_spec$wvl_id<-apply(CABOGeneral2019_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-CABOGeneral2019_spec$wvl_id[duplicated(CABOGeneral2019_spec$wvl_id)]
CABOGeneral2019_spec_no_dups<-CABOGeneral2019_spec[-which(CABOGeneral2019_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(CABOGeneral2019_spec_no_dups$wavelength)):floor(max(CABOGeneral2019_spec_no_dups$wavelength))

CABOGeneral2019_cleaned <-CABOGeneral2019_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

CABOGeneral2019_agg<- CABOGeneral2019_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

CABOGeneral2019_sg_VIS <- CABOGeneral2019_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
CABOGeneral2019_sg_NIR <- CABOGeneral2019_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
CABOGeneral2019_sg_SWIR1 <- CABOGeneral2019_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
CABOGeneral2019_sg_SWIR2 <- CABOGeneral2019_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
CABOGeneral2019_sg <- bind_rows(CABOGeneral2019_sg_VIS,
                      CABOGeneral2019_sg_NIR,
                      CABOGeneral2019_sg_SWIR1,
                      CABOGeneral2019_sg_SWIR2) %>% 
  arrange(sample_id)

CABOGeneral2019_sg_wide<-reshape(subset(as.data.frame(CABOGeneral2019_sg),select= -value),
                       timevar="wvl",idvar="sample_id",direction="wide")
colnames(CABOGeneral2019_sg_wide)<-gsub(pattern = "value_sg.",
                              replacement = "",
                              x = colnames(CABOGeneral2019_sg_wide))

write.csv(CABOGeneral2019_sg_wide,"ProcessedSpectra/CABOGeneral2019_spec_processed.csv",row.names=F)

###########################################################
## Crofts spectra

Crofts_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Crofts.csv")

Crofts_spec<-Crofts_spec[Crofts_spec$reflectance.transmittance=="reflectance",]
Crofts_spec<-Crofts_spec[Crofts_spec$wavelength>=400 & Crofts_spec$wavelength<=2500,]

Crofts_spec$leaf_id<-apply(Crofts_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
Crofts_spec$wvl_id<-apply(Crofts_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-Crofts_spec$wvl_id[duplicated(Crofts_spec$wvl_id)]
Crofts_spec_no_dups<-Crofts_spec[-which(Crofts_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(Crofts_spec_no_dups$wavelength)):floor(max(Crofts_spec_no_dups$wavelength))

Crofts_cleaned <-Crofts_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Crofts_agg<- Crofts_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Crofts_sg_VIS <- Crofts_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Crofts_sg_NIR <- Crofts_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Crofts_sg_SWIR1 <- Crofts_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Crofts_sg_SWIR2 <- Crofts_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Crofts_sg <- bind_rows(Crofts_sg_VIS,
                                Crofts_sg_NIR,
                                Crofts_sg_SWIR1,
                                Crofts_sg_SWIR2) %>% 
  arrange(sample_id)

Crofts_sg_wide<-reshape(subset(as.data.frame(Crofts_sg),select= -value),
                                 timevar="wvl",idvar="sample_id",direction="wide")
colnames(Crofts_sg_wide)<-gsub(pattern = "value_sg.",
                                        replacement = "",
                                        x = colnames(Crofts_sg_wide))

write.csv(Crofts_sg_wide,"ProcessedSpectra/Crofts_spec_processed.csv",row.names=F)

###########################################################
## Pardo spectra

Pardo_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Pardo.csv")

Pardo_spec<-Pardo_spec[Pardo_spec$reflectance.transmittance=="reflectance",]
Pardo_spec<-Pardo_spec[Pardo_spec$wavelength>=400 & Pardo_spec$wavelength<=2500,]

Pardo_spec$leaf_id<-apply(Pardo_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
Pardo_spec$wvl_id<-apply(Pardo_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-Pardo_spec$wvl_id[duplicated(Pardo_spec$wvl_id)]
Pardo_spec_no_dups<-Pardo_spec[-which(Pardo_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(Pardo_spec_no_dups$wavelength)):floor(max(Pardo_spec_no_dups$wavelength))

Pardo_cleaned <-Pardo_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Pardo_agg<- Pardo_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Pardo_sg_VIS <- Pardo_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Pardo_sg_NIR <- Pardo_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Pardo_sg_SWIR1 <- Pardo_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Pardo_sg_SWIR2 <- Pardo_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Pardo_sg <- bind_rows(Pardo_sg_VIS,
                      Pardo_sg_NIR,
                      Pardo_sg_SWIR1,
                      Pardo_sg_SWIR2) %>% 
  arrange(sample_id)

Pardo_sg_wide<-reshape(subset(as.data.frame(Pardo_sg),select= -value),
                       timevar="wvl",idvar="sample_id",direction="wide")
colnames(Pardo_sg_wide)<-gsub(pattern = "value_sg.",
                              replacement = "",
                              x = colnames(Pardo_sg_wide))

write.csv(Pardo_sg_wide,"ProcessedSpectra/Pardo_spec_processed.csv",row.names=F)

###########################################################
## Phragmites temporal spectra

PhragmitesTemporal_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_PhragmitesTemporal.csv")

PhragmitesTemporal_spec<-PhragmitesTemporal_spec[PhragmitesTemporal_spec$reflectance.transmittance=="reflectance",]
PhragmitesTemporal_spec<-PhragmitesTemporal_spec[PhragmitesTemporal_spec$wavelength>=400 & PhragmitesTemporal_spec$wavelength<=2500,]

PhragmitesTemporal_spec$leaf_id<-apply(PhragmitesTemporal_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
PhragmitesTemporal_spec$wvl_id<-apply(PhragmitesTemporal_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-PhragmitesTemporal_spec$wvl_id[duplicated(PhragmitesTemporal_spec$wvl_id)]
PhragmitesTemporal_spec_no_dups<-PhragmitesTemporal_spec[-which(PhragmitesTemporal_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(PhragmitesTemporal_spec_no_dups$wavelength)):floor(max(PhragmitesTemporal_spec_no_dups$wavelength))

PhragmitesTemporal_cleaned <-PhragmitesTemporal_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

PhragmitesTemporal_agg<- PhragmitesTemporal_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

PhragmitesTemporal_sg_VIS <- PhragmitesTemporal_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
PhragmitesTemporal_sg_NIR <- PhragmitesTemporal_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
PhragmitesTemporal_sg_SWIR1 <- PhragmitesTemporal_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
PhragmitesTemporal_sg_SWIR2 <- PhragmitesTemporal_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
PhragmitesTemporal_sg <- bind_rows(PhragmitesTemporal_sg_VIS,
                       PhragmitesTemporal_sg_NIR,
                       PhragmitesTemporal_sg_SWIR1,
                       PhragmitesTemporal_sg_SWIR2) %>% 
  arrange(sample_id)

PhragmitesTemporal_sg_wide<-reshape(subset(as.data.frame(PhragmitesTemporal_sg),select= -value),
                        timevar="wvl",idvar="sample_id",direction="wide")
colnames(PhragmitesTemporal_sg_wide)<-gsub(pattern = "value_sg.",
                               replacement = "",
                               x = colnames(PhragmitesTemporal_sg_wide))

write.csv(PhragmitesTemporal_sg_wide,"ProcessedSpectra/PhragmitesTemporal_spec_processed.csv",row.names=F)

###########################################################
## CABO General Other spectra

CABOGeneralOther_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_CABOGeneralOther.csv")

CABOGeneralOther_spec<-CABOGeneralOther_spec[CABOGeneralOther_spec$reflectance.transmittance=="reflectance",]
CABOGeneralOther_spec<-CABOGeneralOther_spec[CABOGeneralOther_spec$wavelength>=400 & CABOGeneralOther_spec$wavelength<=2500,]

CABOGeneralOther_spec$leaf_id<-apply(CABOGeneralOther_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneralOther_spec$wvl_id<-apply(CABOGeneralOther_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-CABOGeneralOther_spec$wvl_id[duplicated(CABOGeneralOther_spec$wvl_id)]
CABOGeneralOther_spec_no_dups<-CABOGeneralOther_spec[-which(CABOGeneralOther_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(CABOGeneralOther_spec_no_dups$wavelength)):floor(max(CABOGeneralOther_spec_no_dups$wavelength))

CABOGeneralOther_cleaned <-CABOGeneralOther_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

CABOGeneralOther_agg<- CABOGeneralOther_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

CABOGeneralOther_sg_VIS <- CABOGeneralOther_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
CABOGeneralOther_sg_NIR <- CABOGeneralOther_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
CABOGeneralOther_sg_SWIR1 <- CABOGeneralOther_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
CABOGeneralOther_sg_SWIR2 <- CABOGeneralOther_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
CABOGeneralOther_sg <- bind_rows(CABOGeneralOther_sg_VIS,
                                   CABOGeneralOther_sg_NIR,
                                   CABOGeneralOther_sg_SWIR1,
                                   CABOGeneralOther_sg_SWIR2) %>% 
  arrange(sample_id)

CABOGeneralOther_sg_wide<-reshape(subset(as.data.frame(CABOGeneralOther_sg),select= -value),
                                    timevar="wvl",idvar="sample_id",direction="wide")
colnames(CABOGeneralOther_sg_wide)<-gsub(pattern = "value_sg.",
                                           replacement = "",
                                           x = colnames(CABOGeneralOther_sg_wide))

write.csv(CABOGeneralOther_sg_wide,"ProcessedSpectra/CABOGeneralOther_spec_processed.csv",row.names=F)

###########################################################
## Warren spectra

Warren_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Warren.csv")

Warren_spec<-Warren_spec[Warren_spec$reflectance.transmittance=="reflectance",]
Warren_spec<-Warren_spec[Warren_spec$wavelength>=400 & Warren_spec$wavelength<=2500,]

Warren_spec$leaf_id<-apply(Warren_spec[,c("sample_id","leaf_number")],1,paste,collapse="_")
Warren_spec$wvl_id<-apply(Warren_spec[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids<-Warren_spec$wvl_id[duplicated(Warren_spec$wvl_id)]
Warren_spec_no_dups<-Warren_spec[-which(Warren_spec$wvl_id %in% dup_ids),]
inter_wvls <- ceiling(min(Warren_spec_no_dups$wavelength)):floor(max(Warren_spec_no_dups$wavelength))

Warren_cleaned <-Warren_spec_no_dups%>%
  group_by(sample_id,leaf_number) %>%
  do(interpolate(.))

Warren_agg<- Warren_cleaned %>% 
  dplyr::group_by(sample_id, wvl) %>% 
  dplyr::summarise(value = mean(DN, na.rm = T))

Warren_sg_VIS <- Warren_agg %>%
  dplyr::filter(wvl <= 715) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 21))
Warren_sg_NIR <- Warren_agg %>%
  dplyr::filter(wvl > 715,
                wvl <= 1390 ) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 35))
Warren_sg_SWIR1 <- Warren_agg %>%
  dplyr::filter(wvl > 1390,
                wvl <= 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 3, n = 75))
Warren_sg_SWIR2 <- Warren_agg %>%
  dplyr::filter(wvl > 1880) %>% 
  group_by(sample_id) %>% 
  do(sg_filter(., p = 5, n = 175))
Warren_sg <- bind_rows(Warren_sg_VIS,
                                 Warren_sg_NIR,
                                 Warren_sg_SWIR1,
                                 Warren_sg_SWIR2) %>% 
  arrange(sample_id)

Warren_sg_wide<-reshape(subset(as.data.frame(Warren_sg),select= -value),
                                  timevar="wvl",idvar="sample_id",direction="wide")
colnames(Warren_sg_wide)<-gsub(pattern = "value_sg.",
                                         replacement = "",
                                         x = colnames(Warren_sg_wide))

write.csv(Warren_sg_wide,"ProcessedSpectra/Warren_spec_processed.csv",row.names=F)