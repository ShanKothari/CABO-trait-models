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
## With a few exceptions -- like the Dessain and Hacker 2018 projects,
## which used an ASD spectrometer -- there's a lot of copy-pasting
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
zenith.sub <-tbl_df(dplyr::filter(zenith, wvl >= 350) )

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

# convert spectra to tbl data frame
dn.df <- tbl_df(dn)
wr.df <- tbl_df(wr)
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

write.csv(Dessain_sg_wide,"ProcessedSpectra/Dessain_ref_processed.csv",row.names=F)

##########################################################
## Hacker 2018 GOP
## These were also done with an ASD and integrating sphere
## I'm not including this dataset in downstream analyses
## because the bulk samples were collected a year after the spectra
## there's also no stray current correction or calibration
## spectrum for the white reference panel
## this also does not include all the samples from
## Paul Hacker's 2018 campaign -- only the Q. garryana

# ### get paths of .asd files-----
# asd_dirs <- list.dirs('UnprocessedSpectra/Hacker2018Spectra/',
#                       full.names = T)
# asd_paths <- list.files(asd_dirs, full.names = T, pattern = "\\.asd$")
# ## IDs of reflectance spectra (as opposited to transmittance)
# refl <- grep("r0",asd_paths)
# asd_paths <-asd_paths[refl]
# 
# # add filename
# asd_files <- list.files(asd_dirs, full.names = F, pattern = "\\.asd$")
# asd_files_short <- str_sub(asd_files[refl], end = -5)
# 
# #### read all .asd files
# ref <- get_spectra(asd_paths, type = 'reflectance') # get radiances for target
# 
# # convert spectra to tbl data frame
# ref.df <- tbl_df(ref)
# ref.df$fileName <- asd_files_short
# 
# # convert to long format and order by fileName
# ref.long <- gather(ref.df, wvl.char, target.rad, -fileName) %>%
#   mutate(wvl = as.numeric(wvl.char)) %>%
#   arrange(fileName, wvl) %>%
#   dplyr::select(-wvl.char) 
# 
# ## aggregate reps within a leaf
# Hacker2018_spec <- ref.long %>%
#   group_by(substr(ref.long$fileName,1,nchar(ref.long$fileName)-6), wvl) %>%
#   summarise(mean.refl.corr = mean(target.rad)) %>%
#   droplevels()
# 
# ## for consistency with other data sets at this stage
# colnames(Hacker2018_spec)<-c("sample_id","wvl","DN")
# 
# Hacker2018_spec$leaf_number<-unlist(lapply(strsplit(Hacker2018_spec$sample_id,split="_"),function(x) x[[2]]))
# Hacker2018_spec$sample_id<-unlist(lapply(strsplit(Hacker2018_spec$sample_id,split="_"),function(x) x[[1]]))
# 
# Hacker2018_agg<- Hacker2018_spec %>% 
#   dplyr::group_by(sample_id, wvl) %>% 
#   dplyr::summarise(value = mean(DN, na.rm = T))
# 
# Hacker2018_sg_VIS <- Hacker2018_agg %>%
#   dplyr::filter(wvl <= 715) %>% 
#   group_by(sample_id) %>% 
#   do(sg_filter(., p = 3, n = 21))
# Hacker2018_sg_NIR <- Hacker2018_agg %>%
#   dplyr::filter(wvl > 715,
#                 wvl <= 1390 ) %>% 
#   group_by(sample_id) %>% 
#   do(sg_filter(., p = 3, n = 35))
# Hacker2018_sg_SWIR1 <- Hacker2018_agg %>%
#   dplyr::filter(wvl > 1390,
#                 wvl <= 1880) %>% 
#   group_by(sample_id) %>% 
#   do(sg_filter(., p = 3, n = 75))
# Hacker2018_sg_SWIR2 <- Hacker2018_agg %>%
#   dplyr::filter(wvl > 1880) %>% 
#   group_by(sample_id) %>% 
#   do(sg_filter(., p = 5, n = 175))
# Hacker2018_sg <- bind_rows(Hacker2018_sg_VIS,
#                         Hacker2018_sg_NIR,
#                         Hacker2018_sg_SWIR1,
#                         Hacker2018_sg_SWIR2) %>% 
#   arrange(sample_id)
# 
# Hacker2018_sg_wide<-reshape(subset(as.data.frame(Hacker2018_sg),select= -value),
#                          timevar="wvl",idvar="sample_id",direction="wide")
# colnames(Hacker2018_sg_wide)<-gsub(pattern = "value_sg.",
#                                 replacement = "",
#                                 x = colnames(Hacker2018_sg_wide))
# 
# write.csv(Hacker2018_sg_wide,"ProcessedSpectra/Hacker2018_spec_processed.csv",row.names=F)

###########################################################
## Beauchamp-Rioux spectra

## like the rest, these were measured with an SVC
## and integrating sphere

BeauchampRioux_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_BeauchampRioux.csv")

## separate reflectance and transmittance
BeauchampRioux_spec_ref<-BeauchampRioux_spec[BeauchampRioux_spec$reflectance.transmittance=="reflectance",]

## create unique ids for an individual leaf
BeauchampRioux_spec_ref$leaf_id<-apply(BeauchampRioux_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
BeauchampRioux_spec_ref$wvl_id<-apply(BeauchampRioux_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

## remove the spectral overlap region and interpolate over it
dup_ids_ref<-BeauchampRioux_spec_ref$wvl_id[duplicated(BeauchampRioux_spec_ref$wvl_id)]
BeauchampRioux_spec_ref_no_dups<-BeauchampRioux_spec_ref[-which(BeauchampRioux_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(BeauchampRioux_spec_ref_no_dups$wavelength)):floor(max(BeauchampRioux_spec_ref_no_dups$wavelength))

## apply linear interpolation step over spectral overlap
BeauchampRioux_ref_cleaned <-BeauchampRioux_spec_ref_no_dups%>%
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
BeauchampRioux_spec_trans<-BeauchampRioux_spec[BeauchampRioux_spec$reflectance.transmittance=="transmittance",]
BeauchampRioux_spec_trans$leaf_id<-apply(BeauchampRioux_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
BeauchampRioux_spec_trans$wvl_id<-apply(BeauchampRioux_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-BeauchampRioux_spec_trans$wvl_id[duplicated(BeauchampRioux_spec_trans$wvl_id)]
BeauchampRioux_spec_trans_no_dups<-BeauchampRioux_spec_trans[-which(BeauchampRioux_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(BeauchampRioux_spec_trans_no_dups$wavelength)):floor(max(BeauchampRioux_spec_trans_no_dups$wavelength))

BeauchampRioux_trans_cleaned <-BeauchampRioux_spec_trans_no_dups%>%
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

Boucherville2018_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Boucherville2018.csv")

Boucherville2018_spec_ref<-Boucherville2018_spec[Boucherville2018_spec$reflectance.transmittance=="reflectance",]

Boucherville2018_spec_ref$leaf_id<-apply(Boucherville2018_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2018_spec_ref$wvl_id<-apply(Boucherville2018_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Boucherville2018_spec_ref$wvl_id[duplicated(Boucherville2018_spec_ref$wvl_id)]
Boucherville2018_spec_ref_no_dups<-Boucherville2018_spec_ref[-which(Boucherville2018_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Boucherville2018_spec_ref_no_dups$wavelength)):floor(max(Boucherville2018_spec_ref_no_dups$wavelength))

Boucherville2018_ref_cleaned <-Boucherville2018_spec_ref_no_dups%>%
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
Boucherville2018_spec_trans<-Boucherville2018_spec[Boucherville2018_spec$reflectance.transmittance=="transmittance",]
Boucherville2018_spec_trans$leaf_id<-apply(Boucherville2018_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2018_spec_trans$wvl_id<-apply(Boucherville2018_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Boucherville2018_spec_trans$wvl_id[duplicated(Boucherville2018_spec_trans$wvl_id)]
Boucherville2018_spec_trans_no_dups<-Boucherville2018_spec_trans[-which(Boucherville2018_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Boucherville2018_spec_trans_no_dups$wavelength)):floor(max(Boucherville2018_spec_trans_no_dups$wavelength))

Boucherville2018_trans_cleaned <-Boucherville2018_spec_trans_no_dups%>%
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

Girard_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Girard.csv")

Girard_spec_ref<-Girard_spec[Girard_spec$reflectance.transmittance=="reflectance",]

Girard_spec_ref$leaf_id<-apply(Girard_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Girard_spec_ref$wvl_id<-apply(Girard_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Girard_spec_ref$wvl_id[duplicated(Girard_spec_ref$wvl_id)]
Girard_spec_ref_no_dups<-Girard_spec_ref[-which(Girard_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Girard_spec_ref_no_dups$wavelength)):floor(max(Girard_spec_ref_no_dups$wavelength))

Girard_ref_cleaned <-Girard_spec_ref_no_dups%>%
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
Girard_spec_trans<-Girard_spec[Girard_spec$reflectance.transmittance=="transmittance",]
Girard_spec_trans$leaf_id<-apply(Girard_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Girard_spec_trans$wvl_id<-apply(Girard_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Girard_spec_trans$wvl_id[duplicated(Girard_spec_trans$wvl_id)]
Girard_spec_trans_no_dups<-Girard_spec_trans[-which(Girard_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Girard_spec_trans_no_dups$wavelength)):floor(max(Girard_spec_trans_no_dups$wavelength))

Girard_trans_cleaned <-Girard_spec_trans_no_dups%>%
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

Hacker2019_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Hacker2019.csv")

Hacker2019_spec_ref<-Hacker2019_spec[Hacker2019_spec$reflectance.transmittance=="reflectance",]

Hacker2019_spec_ref$leaf_id<-apply(Hacker2019_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Hacker2019_spec_ref$wvl_id<-apply(Hacker2019_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Hacker2019_spec_ref$wvl_id[duplicated(Hacker2019_spec_ref$wvl_id)]
Hacker2019_spec_ref_no_dups<-Hacker2019_spec_ref[-which(Hacker2019_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Hacker2019_spec_ref_no_dups$wavelength)):floor(max(Hacker2019_spec_ref_no_dups$wavelength))

Hacker2019_ref_cleaned <-Hacker2019_spec_ref_no_dups%>%
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
Hacker2019_spec_trans<-Hacker2019_spec[Hacker2019_spec$reflectance.transmittance=="transmittance",]
Hacker2019_spec_trans$leaf_id<-apply(Hacker2019_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Hacker2019_spec_trans$wvl_id<-apply(Hacker2019_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Hacker2019_spec_trans$wvl_id[duplicated(Hacker2019_spec_trans$wvl_id)]
Hacker2019_spec_trans_no_dups<-Hacker2019_spec_trans[-which(Hacker2019_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Hacker2019_spec_trans_no_dups$wavelength)):floor(max(Hacker2019_spec_trans_no_dups$wavelength))

Hacker2019_trans_cleaned <-Hacker2019_spec_trans_no_dups%>%
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

Blanchard_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Blanchard.csv")

Blanchard_spec_ref<-Blanchard_spec[Blanchard_spec$reflectance.transmittance=="reflectance",]

Blanchard_spec_ref$leaf_id<-apply(Blanchard_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Blanchard_spec_ref$wvl_id<-apply(Blanchard_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Blanchard_spec_ref$wvl_id[duplicated(Blanchard_spec_ref$wvl_id)]
Blanchard_spec_ref_no_dups<-Blanchard_spec_ref[-which(Blanchard_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Blanchard_spec_ref_no_dups$wavelength)):floor(max(Blanchard_spec_ref_no_dups$wavelength))

Blanchard_ref_cleaned <-Blanchard_spec_ref_no_dups%>%
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
Blanchard_spec_trans<-Blanchard_spec[Blanchard_spec$reflectance.transmittance=="transmittance",]
Blanchard_spec_trans$leaf_id<-apply(Blanchard_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Blanchard_spec_trans$wvl_id<-apply(Blanchard_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Blanchard_spec_trans$wvl_id[duplicated(Blanchard_spec_trans$wvl_id)]
Blanchard_spec_trans_no_dups<-Blanchard_spec_trans[-which(Blanchard_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Blanchard_spec_trans_no_dups$wavelength)):floor(max(Blanchard_spec_trans_no_dups$wavelength))

Blanchard_trans_cleaned <-Blanchard_spec_trans_no_dups%>%
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

Boucherville2019_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Boucherville2019.csv")
## temporary fix -- there are two spectra per leaf with this ID, and one of each is bad
Boucherville2019_spec<-Boucherville2019_spec[-which(Boucherville2019_spec$sample_id=="44245455"),]

Boucherville2019_spec_ref<-Boucherville2019_spec[Boucherville2019_spec$reflectance.transmittance=="reflectance",]

Boucherville2019_spec_ref$leaf_id<-apply(Boucherville2019_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2019_spec_ref$wvl_id<-apply(Boucherville2019_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Boucherville2019_spec_ref$wvl_id[duplicated(Boucherville2019_spec_ref$wvl_id)]
Boucherville2019_spec_ref_no_dups<-Boucherville2019_spec_ref[-which(Boucherville2019_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Boucherville2019_spec_ref_no_dups$wavelength)):floor(max(Boucherville2019_spec_ref_no_dups$wavelength))

Boucherville2019_ref_cleaned <-Boucherville2019_spec_ref_no_dups%>%
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
Boucherville2019_spec_trans<-Boucherville2019_spec[Boucherville2019_spec$reflectance.transmittance=="transmittance",]
Boucherville2019_spec_trans$leaf_id<-apply(Boucherville2019_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Boucherville2019_spec_trans$wvl_id<-apply(Boucherville2019_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Boucherville2019_spec_trans$wvl_id[duplicated(Boucherville2019_spec_trans$wvl_id)]
Boucherville2019_spec_trans_no_dups<-Boucherville2019_spec_trans[-which(Boucherville2019_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Boucherville2019_spec_trans_no_dups$wavelength)):floor(max(Boucherville2019_spec_trans_no_dups$wavelength))

Boucherville2019_trans_cleaned <-Boucherville2019_spec_trans_no_dups%>%
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

CABOGeneral2019_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_CABOGeneral2019.csv")

CABOGeneral2019_spec_ref<-CABOGeneral2019_spec[CABOGeneral2019_spec$reflectance.transmittance=="reflectance",]

CABOGeneral2019_spec_ref$leaf_id<-apply(CABOGeneral2019_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneral2019_spec_ref$wvl_id<-apply(CABOGeneral2019_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-CABOGeneral2019_spec_ref$wvl_id[duplicated(CABOGeneral2019_spec_ref$wvl_id)]
CABOGeneral2019_spec_ref_no_dups<-CABOGeneral2019_spec_ref[-which(CABOGeneral2019_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(CABOGeneral2019_spec_ref_no_dups$wavelength)):floor(max(CABOGeneral2019_spec_ref_no_dups$wavelength))

CABOGeneral2019_ref_cleaned <-CABOGeneral2019_spec_ref_no_dups%>%
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
CABOGeneral2019_spec_trans<-CABOGeneral2019_spec[CABOGeneral2019_spec$reflectance.transmittance=="transmittance",]
CABOGeneral2019_spec_trans$leaf_id<-apply(CABOGeneral2019_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneral2019_spec_trans$wvl_id<-apply(CABOGeneral2019_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-CABOGeneral2019_spec_trans$wvl_id[duplicated(CABOGeneral2019_spec_trans$wvl_id)]
CABOGeneral2019_spec_trans_no_dups<-CABOGeneral2019_spec_trans[-which(CABOGeneral2019_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(CABOGeneral2019_spec_trans_no_dups$wavelength)):floor(max(CABOGeneral2019_spec_trans_no_dups$wavelength))

CABOGeneral2019_trans_cleaned <-CABOGeneral2019_spec_trans_no_dups%>%
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

Crofts_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Crofts.csv")

Crofts_spec_ref<-Crofts_spec[Crofts_spec$reflectance.transmittance=="reflectance",]

Crofts_spec_ref$leaf_id<-apply(Crofts_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Crofts_spec_ref$wvl_id<-apply(Crofts_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Crofts_spec_ref$wvl_id[duplicated(Crofts_spec_ref$wvl_id)]
Crofts_spec_ref_no_dups<-Crofts_spec_ref[-which(Crofts_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Crofts_spec_ref_no_dups$wavelength)):floor(max(Crofts_spec_ref_no_dups$wavelength))

Crofts_ref_cleaned <-Crofts_spec_ref_no_dups%>%
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
Crofts_spec_trans<-Crofts_spec[Crofts_spec$reflectance.transmittance=="transmittance",]
Crofts_spec_trans$leaf_id<-apply(Crofts_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Crofts_spec_trans$wvl_id<-apply(Crofts_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Crofts_spec_trans$wvl_id[duplicated(Crofts_spec_trans$wvl_id)]
Crofts_spec_trans_no_dups<-Crofts_spec_trans[-which(Crofts_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Crofts_spec_trans_no_dups$wavelength)):floor(max(Crofts_spec_trans_no_dups$wavelength))

Crofts_trans_cleaned <-Crofts_spec_trans_no_dups%>%
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

Pardo_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Pardo.csv")

Pardo_spec_ref<-Pardo_spec[Pardo_spec$reflectance.transmittance=="reflectance",]

Pardo_spec_ref$leaf_id<-apply(Pardo_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Pardo_spec_ref$wvl_id<-apply(Pardo_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Pardo_spec_ref$wvl_id[duplicated(Pardo_spec_ref$wvl_id)]
Pardo_spec_ref_no_dups<-Pardo_spec_ref[-which(Pardo_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Pardo_spec_ref_no_dups$wavelength)):floor(max(Pardo_spec_ref_no_dups$wavelength))

Pardo_ref_cleaned <-Pardo_spec_ref_no_dups%>%
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
Pardo_spec_trans<-Pardo_spec[Pardo_spec$reflectance.transmittance=="transmittance",]
Pardo_spec_trans$leaf_id<-apply(Pardo_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Pardo_spec_trans$wvl_id<-apply(Pardo_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Pardo_spec_trans$wvl_id[duplicated(Pardo_spec_trans$wvl_id)]
Pardo_spec_trans_no_dups<-Pardo_spec_trans[-which(Pardo_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Pardo_spec_trans_no_dups$wavelength)):floor(max(Pardo_spec_trans_no_dups$wavelength))

Pardo_trans_cleaned <-Pardo_spec_trans_no_dups%>%
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

PhragmitesTemporal_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_PhragmitesTemporal.csv")

PhragmitesTemporal_spec_ref<-PhragmitesTemporal_spec[PhragmitesTemporal_spec$reflectance.transmittance=="reflectance",]

PhragmitesTemporal_spec_ref$leaf_id<-apply(PhragmitesTemporal_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
PhragmitesTemporal_spec_ref$wvl_id<-apply(PhragmitesTemporal_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-PhragmitesTemporal_spec_ref$wvl_id[duplicated(PhragmitesTemporal_spec_ref$wvl_id)]
PhragmitesTemporal_spec_ref_no_dups<-PhragmitesTemporal_spec_ref[-which(PhragmitesTemporal_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(PhragmitesTemporal_spec_ref_no_dups$wavelength)):floor(max(PhragmitesTemporal_spec_ref_no_dups$wavelength))

PhragmitesTemporal_ref_cleaned <-PhragmitesTemporal_spec_ref_no_dups%>%
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
PhragmitesTemporal_spec_trans<-PhragmitesTemporal_spec[PhragmitesTemporal_spec$reflectance.transmittance=="transmittance",]
PhragmitesTemporal_spec_trans$leaf_id<-apply(PhragmitesTemporal_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
PhragmitesTemporal_spec_trans$wvl_id<-apply(PhragmitesTemporal_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-PhragmitesTemporal_spec_trans$wvl_id[duplicated(PhragmitesTemporal_spec_trans$wvl_id)]
PhragmitesTemporal_spec_trans_no_dups<-PhragmitesTemporal_spec_trans[-which(PhragmitesTemporal_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(PhragmitesTemporal_spec_trans_no_dups$wavelength)):floor(max(PhragmitesTemporal_spec_trans_no_dups$wavelength))

PhragmitesTemporal_trans_cleaned <-PhragmitesTemporal_spec_trans_no_dups%>%
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

CABOGeneralOther_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_CABOGeneralOther.csv")

CABOGeneralOther_spec_ref<-CABOGeneralOther_spec[CABOGeneralOther_spec$reflectance.transmittance=="reflectance",]

CABOGeneralOther_spec_ref$leaf_id<-apply(CABOGeneralOther_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneralOther_spec_ref$wvl_id<-apply(CABOGeneralOther_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-CABOGeneralOther_spec_ref$wvl_id[duplicated(CABOGeneralOther_spec_ref$wvl_id)]
CABOGeneralOther_spec_ref_no_dups<-CABOGeneralOther_spec_ref[-which(CABOGeneralOther_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(CABOGeneralOther_spec_ref_no_dups$wavelength)):floor(max(CABOGeneralOther_spec_ref_no_dups$wavelength))

CABOGeneralOther_ref_cleaned <-CABOGeneralOther_spec_ref_no_dups%>%
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
CABOGeneralOther_spec_trans<-CABOGeneralOther_spec[CABOGeneralOther_spec$reflectance.transmittance=="transmittance",]
CABOGeneralOther_spec_trans$leaf_id<-apply(CABOGeneralOther_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
CABOGeneralOther_spec_trans$wvl_id<-apply(CABOGeneralOther_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-CABOGeneralOther_spec_trans$wvl_id[duplicated(CABOGeneralOther_spec_trans$wvl_id)]
CABOGeneralOther_spec_trans_no_dups<-CABOGeneralOther_spec_trans[-which(CABOGeneralOther_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(CABOGeneralOther_spec_trans_no_dups$wavelength)):floor(max(CABOGeneralOther_spec_trans_no_dups$wavelength))

CABOGeneralOther_trans_cleaned <-CABOGeneralOther_spec_trans_no_dups%>%
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

Warren_spec<-read.csv("UnprocessedSpectra/project_leaves_combined_Warren.csv")

Warren_spec_ref<-Warren_spec[Warren_spec$reflectance.transmittance=="reflectance",]

Warren_spec_ref$leaf_id<-apply(Warren_spec_ref[,c("sample_id","leaf_number")],1,paste,collapse="_")
Warren_spec_ref$wvl_id<-apply(Warren_spec_ref[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_ref<-Warren_spec_ref$wvl_id[duplicated(Warren_spec_ref$wvl_id)]
Warren_spec_ref_no_dups<-Warren_spec_ref[-which(Warren_spec_ref$wvl_id %in% dup_ids_ref),]
inter_wvls <- ceiling(min(Warren_spec_ref_no_dups$wavelength)):floor(max(Warren_spec_ref_no_dups$wavelength))

Warren_ref_cleaned <-Warren_spec_ref_no_dups%>%
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
Warren_spec_trans<-Warren_spec[Warren_spec$reflectance.transmittance=="transmittance",]
Warren_spec_trans$leaf_id<-apply(Warren_spec_trans[,c("sample_id","leaf_number")],1,paste,collapse="_")
Warren_spec_trans$wvl_id<-apply(Warren_spec_trans[,c("leaf_id","wavelength")],1,paste,collapse="_")

dup_ids_trans<-Warren_spec_trans$wvl_id[duplicated(Warren_spec_trans$wvl_id)]
Warren_spec_trans_no_dups<-Warren_spec_trans[-which(Warren_spec_trans$wvl_id %in% dup_ids_trans),]
inter_wvls <- ceiling(min(Warren_spec_trans_no_dups$wavelength)):floor(max(Warren_spec_trans_no_dups$wavelength))

Warren_trans_cleaned <-Warren_spec_trans_no_dups%>%
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