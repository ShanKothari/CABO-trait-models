setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)

################################################
## LOPEX

LOPEX_ref<-t(read.csv("IndependentValidationData/LOPEX/lopex1993reflectance.csv",row.names = "Sample_."))
LOPEX_ref<-spectra(LOPEX_ref,
                   bands=colnames(LOPEX_ref),
                   names=rownames(LOPEX_ref))

LOPEX_traits<-read.csv("IndependentValidationData/LOPEX/lopex1993metadata.csv")
LOPEX_traits[which(LOPEX_traits== -999,arr.ind=T)]<-NA
meta(LOPEX_ref)$Nmass<-LOPEX_traits$Nitrogen....DW.
meta(LOPEX_ref)$Cmass<-LOPEX_traits$Carbon....DW.
meta(LOPEX_ref)$EWT<-LOPEX_traits$Equivalent.Water.Thickness..g.cm2.
meta(LOPEX_ref)$LMA<-(LOPEX_traits$Dry.Weight..g./1000)/(LOPEX_traits$Leaf.Area..cm2./10000)
meta(LOPEX_ref)$LDMC<-(LOPEX_traits$Dry.Weight..g./LOPEX_traits$Fresh.Weight..g.)*1000
meta(LOPEX_ref)$cellulose1<-LOPEX_traits$Cellulose_1....DW.
meta(LOPEX_ref)$cellulose2<-LOPEX_traits$Cellulose_2....DW.
meta(LOPEX_ref)$dataset<-"LOPEX"
## measures of lignin are highly inconsistent, so we ignore them here

bad_spectra_LOPEX<-c("X0176","X0177","X0178","X0179","X0180",
               "X0196","X0197","X0198","X0199","X0200",
               "X0321","X0322","X0323","X0324","X0325")
LOPEX_ref<-LOPEX_ref[-which(names(LOPEX_ref) %in% bad_spectra_LOPEX)]

LOPEX_ref<-match_sensors(LOPEX_ref,splice_at = 862,fixed_sensor = 2,interpolate_wvl = 5)
## need to do smoothing here

################################################
## ANGERS

ANGERS_ref<-t(read.csv("IndependentValidationData/ANGERS/angers2003reflectance.csv",row.names = "Sample_."))
ANGERS_ref<-spectra(ANGERS_ref,
                   bands=colnames(ANGERS_ref),
                   names=rownames(ANGERS_ref))

ANGERS_traits<-read.csv("IndependentValidationData/ANGERS/angers2003metadata.csv")
ANGERS_traits[which(ANGERS_traits== -999,arr.ind=T)]<-NA
meta(ANGERS_ref)$EWT<-ANGERS_traits$Equivalent.Water.Thickness..g.cm2.
meta(ANGERS_ref)$LMA<-ANGERS_traits$Leaf.mass.per.area..g.cm2.*10
meta(ANGERS_ref)$dataset<-"ANGERS"

bad_spectra_ANGERS<-c("X0178","X0179","X0184","X0185","X0196",
               "X0197","X0241","X0250","X0254","X0257",
               "X0258","X0269")
ANGERS_ref<-ANGERS_ref[-which(names(ANGERS_ref) %in% bad_spectra_ANGERS)]

## smooth here

#################################################
## combine

LOPEX_ref<-LOPEX_ref[,400:2400]
ANGERS_ref<-ANGERS_ref[,400:2400]
all_val_ref<-combine(LOPEX_ref,ANGERS_ref)

meta(all_val_ref)$Nmass_pred<-predict(Nmass_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_Nmass_CVmodel)[,,1]
meta(all_val_ref)$Cmass_pred<-predict(Cmass_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_Cmass_CVmodel)[,,1]
meta(all_val_ref)$EWT_pred<-predict(EWT_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_EWT_CVmodel)[,,1]
meta(all_val_ref)$LMA_pred<-predict(LMA_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_LMA_CVmodel)[,,1]
meta(all_val_ref)$LDMC_pred<-predict(LDMC_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_LDMC_CVmodel)[,,1]
meta(all_val_ref)$cellulose_pred<-predict(cellulose_mass_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_cellulose_mass_CVmodel)[,,1]

ggplot(data=meta(all_val_ref),aes(x=EWT_pred,y=EWT,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured EWT",x="Predicted EWT")+
  coord_cartesian(xlim=c(0,0.085),ylim=c(0,0.085))

ggplot(data=meta(all_val_ref),aes(x=LDMC_pred,y=LDMC*1000,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured LDMC",x="Predicted LDMC")

ggplot(data=meta(all_val_ref),aes(x=LMA_pred,y=LMA*1000,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured LMA",x="Predicted LMA")

##########################################
## Dessain

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
Dessain <- refl_corr.df %>%
  group_by(parentEventID,leafNumber, wvl) %>%
  summarise(mean.refl.corr = mean(refl_corr)) %>%
  droplevels()

## for consistency with other data sets at this stage
colnames(Dessain)<-c("sample_id","leaf_number","wvl","DN")

Dessain_agg<- Dessain %>% 
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

Dessain.df<-reshape(subset(as.data.frame(Dessain_sg),select= -value),
                         timevar="wvl",idvar="sample_id",direction="wide")
colnames(Dessain.df)<-gsub(pattern = "value_sg.",
                                replacement = "X",
                                x = colnames(Dessain.df))

Dessain.spec<-as_spectra(Dessain.df,name_idx=1)
meta(Dessain.spec)$sample_id<-gsub("spec_","",names(Dessain.spec))
Dessain.spec<-Dessain.spec[,400:2400]

Dessain.summary<-read.csv("SummaryData/DessainHerbarium.csv")
Dessain.sub<-data.frame(sample.id=Dessain.summary$parentEventID,
                        species=Dessain.summary$scientificName,
                        project="2017-Dessain-MSc")

meta(Dessain.spec)$species<-Dessain.sub$species[match(meta(Dessain.spec)$sample_id,Dessain.sub$sample.id)]
meta(Dessain.spec)$project<-Dessain.sub$project[match(meta(Dessain.spec)$sample_id,Dessain.sub$sample.id)]

Dessain.area<-read.csv("TraitData/LeafAreaWaterSamples/SLA_data_Aurelie_Dessain - Lab_data.csv")
Dessain.area$SLA_m2_kg<-as.numeric(as.character(Dessain.area$SLA_m2_kg))
Dessain.area$LDMC_mg_g<-as.numeric(as.character(Dessain.area$LDMC_mg_g))

Dessain.area.sub<-data.frame(sample_id=Dessain.area$parentEventID,
                             SLA=Dessain.area$SLA_m2_kg,
                             LDMC=Dessain.area$LDMC_mg_g)

## SLA in units m^2/kg
## LDMC in units mg/g
## EWT in units cm
meta(Dessain.spec)$SLA<-Dessain.area.sub$SLA[match(meta(Dessain.spec)$sample_id,Dessain.area.sub$sample_id)]
meta(Dessain.spec)$LDMC<-Dessain.area.sub$LDMC[match(meta(Dessain.spec)$sample_id,Dessain.area.sub$sample_id)]
meta(Dessain.spec)$EWT<-with(meta(Dessain.spec),(1/(LDMC/1000)-1)*(1/SLA*0.1))

Dessain.CN<-read.csv("TraitData/CNSamples/2017_Dessain_MSc_CN_data_total.csv")
Dessain.CN.sub<-data.frame(sample_id=Dessain.CN$Sample_id,
                           Cmass=Dessain.CN$C.....,
                           Nmass=Dessain.CN$N....)

meta(Dessain.spec)$Cmass<-Dessain.CN$Cmass[match(meta(Dessain.spec)$sample_id,Dessain.CN$sample_id)]
meta(Dessain.spec)$Nmass<-Dessain.CN$Nmass[match(meta(Dessain.spec)$sample_id,Dessain.CN$sample_id)]

Dessain.pigments<-read.csv("TraitData/Pigments/Aurelie_pigments_valeurs_brutes - Analyses_Aurelie.csv")
Dessain.pigments<-Dessain.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]
## remove records who notes contain "refaire" or "refait" since these were all redone
Dessain.pigments<-Dessain.pigments[-which(str_detect(string = Dessain.pigments$Notes,pattern = "efai")),]

meta(Dessain.spec)$chlA_mass<-as.numeric(as.character(Dessain.pigments$chlA_mg_g[match(meta(Dessain.spec)$sample_id,Dessain.pigments$sample_id)]))
meta(Dessain.spec)$chlB_mass<-as.numeric(as.character(Dessain.pigments$chlB_mg_g[match(meta(Dessain.spec)$sample_id,Dessain.pigments$sample_id)]))
meta(Dessain.spec)$car_mass<-as.numeric(as.character(Dessain.pigments$carotenoides._mg_g[match(meta(Dessain.spec)$sample_id,Dessain.pigments$sample_id)]))

##########################################################
## Hacker 2018 GOP
## These were also done with an ASD and integrating sphere
## bulk samples were collected a year after the spectra
## there's also no stray current correction or calibration
## spectrum for the white reference panel

## this also does not include all the samples from
## Paul Hacker's 2018 campaign -- only the Q. garryana
## (can I fix this?)

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
# Hacker2018 <- ref.long %>%
#   group_by(substr(ref.long$fileName,1,nchar(ref.long$fileName)-6), wvl) %>%
#   summarise(mean.refl.corr = mean(target.rad)) %>%
#   droplevels()
# 
# ## for consistency with other data sets at this stage
# colnames(Hacker2018)<-c("sample_id","wvl","DN")
# 
# Hacker2018$leaf_number<-unlist(lapply(strsplit(Hacker2018$sample_id,split="_"),function(x) x[[2]]))
# Hacker2018$sample_id<-unlist(lapply(strsplit(Hacker2018$sample_id,split="_"),function(x) x[[1]]))
# 
# Hacker2018_agg<- Hacker2018 %>% 
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
# write.csv(Hacker2018_sg_wide,"ProcessedSpectra/Hacker2018_ref_processed.csv",row.names=F)
