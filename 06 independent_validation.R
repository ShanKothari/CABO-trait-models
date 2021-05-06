setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(tidyverse)
library(spectrolab)
library(asdreader)
library(signal)

##########################################################
## defining important functions for later

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

## applying coefficients to validation spectra
apply.coefs<-function(coef.list,val.spec,intercept=T){
  if(sum(lapply(coef.list,length)==ncol(val.spec)+intercept) < length(coef.list)){
    stop("some coefficients have the wrong length")
  }
  
  coef.matrix<-matrix(unlist(coef.list),
                      nrow=length(coef.list),
                      byrow=T)
  
  if(intercept==T){
    pred.matrix<-t(t(as.matrix(val.spec) %*% t(coef.matrix[,-1]))+coef.matrix[,1])
  } else {
    pred.matrix<-as.matrix(val.spec) %*% t(coef.matrix)
  }
}

################################################
## LOPEX

LOPEX.spec<-t(read.csv("IndependentValidationData/LOPEX/lopex1993reflectance.csv",row.names = "Sample_."))
bad_spectra_LOPEX<-c("X0176","X0177","X0178","X0179","X0180",
                     "X0196","X0197","X0198","X0199","X0200",
                     "X0321","X0322","X0323","X0324","X0325")
LOPEX.spec<-LOPEX.spec[-which(rownames(LOPEX.spec) %in% bad_spectra_LOPEX),]

LOPEX.spec<-as.matrix(match_sensors(spectra(LOPEX.spec,
                                            bands=400:2500,
                                            names=rownames(LOPEX.spec)),
                                    splice_at = 862,
                                    fixed_sensor = 2,
                                    interpolate_wvl = 5))

LOPEX_sg_VIS <- t(apply(LOPEX.spec[,as.character(400:715)],1,
                        function(x) sgolayfilt(x,p=3,n=21)))
LOPEX_sg_NIR <- t(apply(LOPEX.spec[,as.character(716:1390)],1,
                        function(x) sgolayfilt(x,p=3,n=35)))
LOPEX_sg_SWIR1 <- t(apply(LOPEX.spec[,as.character(1391:1880)],1,
                          function(x) sgolayfilt(x,p=3,n=75)))
LOPEX_sg_SWIR2 <- t(apply(LOPEX.spec[,as.character(1881:2500)],1,
                          function(x) sgolayfilt(x,p=5,n=175)))
LOPEX_sg<-do.call(cbind,args = list(LOPEX_sg_VIS,
                                    LOPEX_sg_NIR,
                                    LOPEX_sg_SWIR1,
                                    LOPEX_sg_SWIR2))
LOPEX<-spectra(LOPEX_sg,
               bands = colnames(LOPEX.spec),
               names = rownames(LOPEX.spec))

LOPEX_traits<-read.csv("IndependentValidationData/LOPEX/lopex1993metadata.csv")
LOPEX_traits<-LOPEX_traits[-which(paste("X0",LOPEX_traits$Sample_.,sep="") %in% bad_spectra_LOPEX),]
LOPEX_traits[which(LOPEX_traits== -999,arr.ind=T)]<-NA
meta(LOPEX)$chlA<-LOPEX_traits$Chlorophyll_a..µg.cm2.
meta(LOPEX)$Nmass<-LOPEX_traits$Nitrogen....DW.
meta(LOPEX)$Cmass<-LOPEX_traits$Carbon....DW.
meta(LOPEX)$EWT<-LOPEX_traits$Equivalent.Water.Thickness..g.cm2.
meta(LOPEX)$LMA<-LOPEX_traits$Leaf.mass.per.area..g.cm2.*10
meta(LOPEX)$LDMC<-(LOPEX_traits$Dry.Weight..g./LOPEX_traits$Fresh.Weight..g.)*1000
meta(LOPEX)$chlA<-LOPEX_traits$Chlorophyll_a..µg.cm2./(100*meta(LOPEX)$LMA)
meta(LOPEX)$chlB<-LOPEX_traits$Chlorophyll_b..µg.cm2./(100*meta(LOPEX)$LMA)
meta(LOPEX)$car<-LOPEX_traits$Carotenoid..µg.cm2./(100*meta(LOPEX)$LMA)

## we average the two methods for estimating cellulose and lignin
meta(LOPEX)$cellulose1<-LOPEX_traits$Cellulose_1....DW.
meta(LOPEX)$cellulose2<-LOPEX_traits$Cellulose_2....DW.
meta(LOPEX)$cellulose<-rowMeans(data.frame(meta(LOPEX)$cellulose1,
                                               meta(LOPEX)$cellulose2))
meta(LOPEX)$lignin1<-LOPEX_traits$Lignin_1....DW.
meta(LOPEX)$lignin2<-LOPEX_traits$Lignin_2....DW.
meta(LOPEX)$lignin<-rowMeans(data.frame(meta(LOPEX)$lignin1,
                                                meta(LOPEX)$lignin2))

meta(LOPEX)$dataset<-"LOPEX"

################################################
## ANGERS

ANGERS.spec<-t(read.csv("IndependentValidationData/ANGERS/angers2003reflectance.csv",row.names = "Sample_."))

bad_spectra_ANGERS<-c("X0178","X0179","X0184","X0185","X0196",
                      "X0197","X0241","X0250","X0254","X0257",
                      "X0258","X0269")
ANGERS.spec<-ANGERS.spec[-which(rownames(ANGERS.spec) %in% bad_spectra_ANGERS),]

## matching sensors is probably not necessary for ANGERS
## smoothing
ANGERS_sg_VIS <- t(apply(ANGERS.spec[,as.character(400:715)],1,
                         function(x) sgolayfilt(x,p=3,n=21)))
ANGERS_sg_NIR <- t(apply(ANGERS.spec[,as.character(716:1390)],1,
                         function(x) sgolayfilt(x,p=3,n=35)))
ANGERS_sg_SWIR1 <- t(apply(ANGERS.spec[,as.character(1391:1880)],1,
                           function(x) sgolayfilt(x,p=3,n=75)))
ANGERS_sg_SWIR2 <- t(apply(ANGERS.spec[,as.character(1881:2450)],1,
                           function(x) sgolayfilt(x,p=5,n=175)))
ANGERS_sg<-do.call(cbind,args = list(ANGERS_sg_VIS,
                                     ANGERS_sg_NIR,
                                     ANGERS_sg_SWIR1,
                                     ANGERS_sg_SWIR2))
ANGERS<-spectra(ANGERS_sg,
                bands=colnames(ANGERS.spec),
                names = rownames(ANGERS.spec))

ANGERS_traits<-read.csv("IndependentValidationData/ANGERS/angers2003metadata.csv")
ANGERS_traits<-ANGERS_traits[-which(paste("X0",ANGERS_traits$Sample_.,sep="") %in% bad_spectra_ANGERS),]
ANGERS_traits[which(ANGERS_traits== -999,arr.ind=T)]<-NA
meta(ANGERS)$EWT<-ANGERS_traits$Equivalent.Water.Thickness..g.cm2.
meta(ANGERS)$LMA<-ANGERS_traits$Leaf.mass.per.area..g.cm2.*10
meta(ANGERS)$chlA<-ANGERS_traits$Chlorophyll_a..µg.cm2./(100*meta(ANGERS)$LMA)
meta(ANGERS)$chlB<-ANGERS_traits$Chlorophyll_b..µg.cm2./(100*meta(ANGERS)$LMA)
meta(ANGERS)$car<-ANGERS_traits$Carotenoid..µg.cm2./(100*meta(ANGERS)$LMA)
meta(ANGERS)$dataset<-"ANGERS"

##########################################################
## Dessain
## These were done with an ASD and integrating sphere
## (the other projects with an SVC and integrating sphere)
## the spectral overlap region is already spliced out
## but I do need to do dark current correction, which isn't
## true with the other projects

zenith <- read.table('UnprocessedSpectra/DessainSpectra/c7101904_specchioformat.txt',
                     header = T)
zenith.sub <-as_tibble(dplyr::filter(zenith, wvl >= 350) )

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
dn.df <- as_tibble(dn)
wr.df <- as_tibble(wr)
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

colnames(Dessain)<-c("sample_id","leaf_number","wvl","DN")

## now aggregate leaves within a sample
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
                                replacement = "",
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

## attach structural traits
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
meta(Dessain.spec)$LMA<-1/meta(Dessain.spec)$SLA
meta(Dessain.spec)$LDMC<-Dessain.area.sub$LDMC[match(meta(Dessain.spec)$sample_id,Dessain.area.sub$sample_id)]
meta(Dessain.spec)$EWT<-with(meta(Dessain.spec),(1/(LDMC/1000)-1)*(1/SLA*0.1))

## attach C and N
Dessain.CN<-read.csv("TraitData/CNSamples/2017_Dessain_MSc_CN_data_total.csv")
Dessain.CN.sub<-data.frame(sample_id=Dessain.CN$Sample_id,
                           Cmass=Dessain.CN$C.....,
                           Nmass=Dessain.CN$N....)

meta(Dessain.spec)$Cmass<-Dessain.CN.sub$Cmass[match(meta(Dessain.spec)$sample_id,Dessain.CN.sub$sample_id)]
meta(Dessain.spec)$Nmass<-Dessain.CN.sub$Nmass[match(meta(Dessain.spec)$sample_id,Dessain.CN.sub$sample_id)]

## attach carbon fractions
CFractions<-read.csv("TraitData/CarbonFractions/carbon_fractions_bags.csv")
meta(Dessain.spec)$solubles<-CFractions$soluble_perc[match(meta(Dessain.spec)$sample_id,CFractions$bottle_id)]
meta(Dessain.spec)$hemicellulose<-CFractions$hemicellulose_perc[match(meta(Dessain.spec)$sample_id,CFractions$bottle_id)]
meta(Dessain.spec)$cellulose<-CFractions$cellulose_perc[match(meta(Dessain.spec)$sample_id,CFractions$bottle_id)]
meta(Dessain.spec)$lignin<-CFractions$lignin_perc[match(meta(Dessain.spec)$sample_id,CFractions$bottle_id)]

## attach pigments
Dessain.pigments<-read.csv("TraitData/Pigments/Aurelie_pigments_valeurs_brutes - Analyses_Aurelie.csv")
Dessain.pigments<-Dessain.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]
## remove records who notes contain "refaire" or "refait" since these were all redone
Dessain.pigments<-Dessain.pigments[-which(str_detect(string = Dessain.pigments$Notes,pattern = "efai")),]

meta(Dessain.spec)$chlA<-as.numeric(as.character(Dessain.pigments$chlA_mg_g[match(meta(Dessain.spec)$sample_id,Dessain.pigments$sample_id)]))
meta(Dessain.spec)$chlB<-as.numeric(as.character(Dessain.pigments$chlB_mg_g[match(meta(Dessain.spec)$sample_id,Dessain.pigments$sample_id)]))
meta(Dessain.spec)$car<-as.numeric(as.character(Dessain.pigments$carotenoides._mg_g[match(meta(Dessain.spec)$sample_id,Dessain.pigments$sample_id)]))

## attach ICP data
ICP12<-read.csv("TraitData/ICP/2020-03-12 1MHCl_Etienne box1,2.csv")
ICP34<-read.csv("TraitData/ICP/2020-10-20 1MHCl_Etienne box 3,4.csv")
ICP5<-read.csv("TraitData/ICP/2020-10-21 1MHCl_Etienne box 5.csv")
ICP67<-read.csv("TraitData/ICP/2020-10-22 1MHCl_Etienne box 6,7.csv")
ICP8<-read.csv("TraitData/ICP/2020-10-23 1MHCl_Etienne box 8.csv")
ICP9<-read.csv("TraitData/ICP/2020-11-11 1MHCl_Etienne box 9.csv")
ICP1011<-read.csv("TraitData/ICP/2020-12-16 1MHCl_Etienne box 10,11.csv")
ICP_all<-do.call(rbind,args=list(ICP12,ICP34,ICP5,ICP67,
                                   ICP8,ICP9,ICP1011))
num_cols<-c("Al","B","B.1","Ca","Cu","Fe","K","Mg","Mn","Na","P","Zn")
ICP_all[,num_cols]<-data.frame(sapply(ICP_all[,num_cols],
                                        function(x) as.numeric(as.character(x))))

ICP_all$Al[ICP_all$Sample_id=="2017-08-15-jbmcb-P006"]<-NA
ICP_all$Cu[ICP_all$Sample_id=="2017-06-07-ireqa-P010"]<-NA

meta(Dessain.spec)$Al_mass<-ICP_all$Al[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$B_mass<-ICP_all$B[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$B.1_mass<-ICP_all$B.1[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$Ca_mass<-ICP_all$Ca[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$Cu_mass<-ICP_all$Cu[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$Fe_mass<-ICP_all$Fe[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$K_mass<-ICP_all$K[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$Mg_mass<-ICP_all$Mg[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$Mn_mass<-ICP_all$Mn[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$Na_mass<-ICP_all$Na[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$P_mass<-ICP_all$P[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]
meta(Dessain.spec)$Zn_mass<-ICP_all$Zn[match(meta(Dessain.spec)$sample_id,ICP_all$Sample_id)]

meta(Dessain.spec)$dataset<-"Dessain"

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

#################################################
## combine

## !!!!!!!!!
## might want to add Hacker data

LOPEX<-LOPEX[,400:2400]
ANGERS<-ANGERS[,400:2400]
Dessain.spec<-Dessain.spec[,400:2400]

LOPEX_ANGERS<-spectrolab::combine(LOPEX,ANGERS)
all_val_ref<-spectrolab::combine(LOPEX_ANGERS,Dessain.spec)

all_jack_coefs_list_ref<-readRDS("SavedResults/all_jack_coefs_list_ref.rds")
meta(all_val_ref)$solubles_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$solubles_mass,val.spec = all_val_ref))
meta(all_val_ref)$hemicellulose_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$hemicellulose_mass,val.spec = all_val_ref))
meta(all_val_ref)$cellulose_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$cellulose_mass,val.spec = all_val_ref))
meta(all_val_ref)$lignin_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$lignin_mass,val.spec = all_val_ref))
meta(all_val_ref)$Nmass_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Nmass,val.spec = all_val_ref))
meta(all_val_ref)$Cmass_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Cmass,val.spec = all_val_ref))
meta(all_val_ref)$LMA_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$LMA,val.spec = all_val_ref))
meta(all_val_ref)$LDMC_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$LDMC,val.spec = all_val_ref))
meta(all_val_ref)$EWT_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$EWT,val.spec = all_val_ref))
meta(all_val_ref)$chlA_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$chlA_mass,val.spec = all_val_ref))
meta(all_val_ref)$chlB_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$chlB_mass,val.spec = all_val_ref))
meta(all_val_ref)$car_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$car_mass,val.spec = all_val_ref))
meta(all_val_ref)$Al_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Al_mass,val.spec = all_val_ref))
meta(all_val_ref)$Ca_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Ca_mass,val.spec = all_val_ref))
meta(all_val_ref)$Cu_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Cu_mass,val.spec = all_val_ref))
meta(all_val_ref)$Fe_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Fe_mass,val.spec = all_val_ref))
meta(all_val_ref)$K_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$K_mass,val.spec = all_val_ref))
meta(all_val_ref)$Mg_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Mg_mass,val.spec = all_val_ref))
meta(all_val_ref)$Mn_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Mn_mass,val.spec = all_val_ref))
meta(all_val_ref)$Na_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Na_mass,val.spec = all_val_ref))
meta(all_val_ref)$P_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$P_mass,val.spec = all_val_ref))
meta(all_val_ref)$Zn_pred<-rowMeans(apply.coefs(all_jack_coefs_list_ref$Zn_mass,val.spec = all_val_ref))

solubles_ind_val<-ggplot(data=meta(all_val_ref),aes(x=solubles_pred,y=solubles,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured solubles",x="Predicted solubles")+
  coord_cartesian(xlim=c(20,90),ylim=c(20,90))

hemicellulose_ind_val<-ggplot(data=meta(all_val_ref),aes(x=hemicellulose_pred,y=hemicellulose,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured hemicellulose",x="Predicted hemicellulose")+
  coord_cartesian(xlim=c(-5,35),ylim=c(-5,35))

cellulose_ind_val<-ggplot(data=meta(all_val_ref),aes(x=cellulose_pred,y=cellulose,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured cellulose",x="Predicted cellulose")+
  coord_cartesian(xlim=c(0,50),ylim=c(0,50))

lignin_ind_val<-ggplot(data=meta(all_val_ref),aes(x=lignin_pred,y=lignin,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured lignin",x="Predicted lignin")+
  coord_cartesian(xlim=c(-12,25),ylim=c(-12,25))

Nmass_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Nmass_pred,y=Nmass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Nmass",x="Predicted Nmass")+
  coord_cartesian(xlim=c(0,6),ylim=c(0,6))

Cmass_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Cmass_pred,y=Cmass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Cmass",x="Predicted Cmass")+
  coord_cartesian(xlim=c(37,54),ylim=c(37,54))

LMA_ind_val<-ggplot(data=meta(all_val_ref),aes(x=LMA_pred,y=LMA,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured LMA",x="Predicted LMA")+
  coord_cartesian(xlim=c(0,0.5),ylim=c(0,0.5))

LDMC_ind_val<-ggplot(data=meta(all_val_ref),aes(x=LDMC_pred,y=LDMC,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured LDMC",x="Predicted LDMC")+
  coord_cartesian(xlim=c(-75,600),ylim=c(-75,600))

EWT_ind_val<-ggplot(data=meta(all_val_ref),aes(x=EWT_pred,y=EWT,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured EWT",x="Predicted EWT")+
  coord_cartesian(xlim=c(0,0.085),ylim=c(0,0.085))

chlA_ind_val<-ggplot(data=meta(all_val_ref),aes(x=chlA_pred,y=chlA,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured chlA",x="Predicted chlA")+
  coord_cartesian(xlim=c(0,20),ylim=c(0,20))

chlB_ind_val<-ggplot(data=meta(all_val_ref),aes(x=chlB_pred,y=chlB,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured chlB",x="Predicted chlB")+
  coord_cartesian(xlim=c(-1,7),ylim=c(-1,7))

car_ind_val<-ggplot(data=meta(all_val_ref),aes(x=car_pred,y=car,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured car",x="Predicted car")+
  coord_cartesian(xlim=c(0,5.5),ylim=c(0,5.5))

Al_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Al_pred,y=Al_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Al",x="Predicted Al")+
  coord_cartesian(xlim=c(-0.1,0.4),ylim=c(-0.1,0.4))

Ca_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Ca_pred,y=Ca_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Al",x="Predicted Al")+
  coord_cartesian(xlim=c(-5,30),ylim=c(-5,30))

Cu_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Cu_pred,y=Cu_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Cu",x="Predicted Cu")+
  coord_cartesian(xlim=c(0,0.1),ylim=c(0,0.1))

Fe_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Fe_pred,y=Fe_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Fe",x="Predicted Fe")+
  coord_cartesian(xlim=c(0,0.5),ylim=c(0,0.5))

K_ind_val<-ggplot(data=meta(all_val_ref),aes(x=K_pred,y=K_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured K",x="Predicted K")+
  coord_cartesian(xlim=c(0,40),ylim=c(0,40))

Mg_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Mg_pred,y=Mg_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Mg",x="Predicted Mg")+
  coord_cartesian(xlim=c(-1,8),ylim=c(-1,8))

Mn_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Mn_pred,y=Mn_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Mn",x="Predicted Mn")+
  coord_cartesian(xlim=c(-1,4),ylim=c(-1,4))

Na_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Na_pred,y=Na_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Na",x="Predicted Na")+
  coord_cartesian(xlim=c(-1,8),ylim=c(-1,8))

P_ind_val<-ggplot(data=meta(all_val_ref),aes(x=P_pred,y=P_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured P",x="Predicted P")+
  coord_cartesian(xlim=c(-2,7),ylim=c(-2,7))

Zn_ind_val<-ggplot(data=meta(all_val_ref),aes(x=Zn_pred,y=Zn_mass,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured Zn",x="Predicted Zn")+
  coord_cartesian(xlim=c(-0.1,0.7),ylim=c(-0.1,0.7))

pdf("Images/ind_val_plots.pdf",width = 10,height = 8)
solubles_ind_val
hemicellulose_ind_val
cellulose_ind_val
lignin_ind_val
Nmass_ind_val
Cmass_ind_val
LMA_ind_val
LDMC_ind_val
EWT_ind_val
chlA_ind_val
chlB_ind_val
car_ind_val
Al_ind_val
Ca_ind_val
Cu_ind_val
Fe_ind_val
K_ind_val
Mg_ind_val
Mn_ind_val
Na_ind_val
P_ind_val
Zn_ind_val
dev.off()
