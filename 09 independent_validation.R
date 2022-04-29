setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(tidyverse)
library(spectrolab)
library(asdreader)
library(signal)
library(ggpubr)
source("Scripts/CABO-trait-models/00 useful_functions.R")

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
meta(LOPEX)$EWT<-LOPEX_traits$Equivalent.Water.Thickness..g.cm2.*10
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
LOPEX<-LOPEX[,400:2400]

saveRDS(LOPEX,"IndependentValidationData/LOPEX/LOPEX_processed.rds")

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
meta(ANGERS)$EWT<-ANGERS_traits$Equivalent.Water.Thickness..g.cm2.*10
meta(ANGERS)$LMA<-ANGERS_traits$Leaf.mass.per.area..g.cm2.*10
meta(ANGERS)$LDMC<-with(meta(ANGERS),LMA/(LMA+EWT))*1000
meta(ANGERS)$chlA<-ANGERS_traits$Chlorophyll_a..µg.cm2./(100*meta(ANGERS)$LMA)
meta(ANGERS)$chlB<-ANGERS_traits$Chlorophyll_b..µg.cm2./(100*meta(ANGERS)$LMA)
meta(ANGERS)$car<-ANGERS_traits$Carotenoid..µg.cm2./(100*meta(ANGERS)$LMA)

meta(ANGERS)$dataset<-"ANGERS"
ANGERS<-ANGERS[,400:2400]

saveRDS(ANGERS,"IndependentValidationData/ANGERS/ANGERS_processed.rds")

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

Dessain.summary<-read.csv("SummaryData/DessainHerbarium.csv")
Dessain.sub<-data.frame(sample.id=Dessain.summary$parentEventID,
                        species=Dessain.summary$scientificName,
                        project="2017-Dessain-MSc")

meta(Dessain.spec)$species<-Dessain.sub$species[match(meta(Dessain.spec)$sample_id,Dessain.sub$sample.id)]
meta(Dessain.spec)$project<-Dessain.sub$project[match(meta(Dessain.spec)$sample_id,Dessain.sub$sample.id)]

vascan<-read.csv("SummaryData/vascan.csv")
meta(Dessain.spec)$growth.form<-vascan$growth_form[match(meta(Dessain.spec)$species,vascan$scientific_name)]
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Acer pensylvanicum Linnaeus"]<-"tree"
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Betula populifolia Marshall"]<-"tree"
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Prunus pensylvanica Linnaeus f."]<-"tree"
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Rhamnus cathartica Linnaeus"]<-"shrub"
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Prunus virginiana Linnaeus"]<-"shrub"
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Staphylea trifolia Linnaeus"]<-"shrub"
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Salix Linnaeus"]<-"shrub"
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Aesculus glabra Wildenow"]<-"tree"
meta(Dessain.spec)$growth.form[meta(Dessain.spec)$species=="Ginkgo biloba L."]<-"tree"

meta(Dessain.spec)$order<-vascan$order[match(meta(Dessain.spec)$species,vascan$scientific_name)]
meta(Dessain.spec)$family<-vascan$family[match(meta(Dessain.spec)$species,vascan$scientific_name)]
meta(Dessain.spec)$genus<-vascan$genus[match(meta(Dessain.spec)$species,vascan$scientific_name)]
which(meta(Dessain.spec)$family %in% c("Poaceae","Cyperaceae","Juncaceae","Typhaceae"))
which(meta(Dessain.spec)$family %in% c("Pinaceae","Cupressaceae"))

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
meta(Dessain.spec)$EWT<-with(meta(Dessain.spec),(1/(LDMC/1000)-1)*(1/SLA*0.1)*10)

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

NDF_bad<-c("Ulmus rubra Muhlenberg","Acer platanoides Linnaeus",
           "Ulmus americana Linnaeus","Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")
meta(Dessain.spec)$solubles[meta(Dessain.spec)$species %in% NDF_bad]<-NA
meta(Dessain.spec)$hemicellulose[meta(Dessain.spec)$species %in% NDF_bad]<-NA

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
ICP_all$Cu[ICP_all$Sample_id=="2017-08-15-jbmcb-P002"]<-NA
ICP_all$Fe[ICP_all$Sample_id=="2017-05-31-jbmtb-P011"]<-NA
ICP_all$Fe[ICP_all$Sample_id=="2017-05-26-jbmcb-P007"]<-NA
ICP_all$K[ICP_all$Sample_id=="2017-06-15-jbmcb-P003"]<-NA
ICP_all$Mg[ICP_all$Sample_id=="2017-06-03-sutod-P003"]<-NA
ICP_all$Na[ICP_all$Sample_id=="2017-06-15-jbmcb-P003"]<-NA

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
Dessain.spec<-Dessain.spec[,400:2400]

saveRDS(Dessain.spec,"IndependentValidationData/Dessain_processed.rds")

##########################################################
## Hacker 2018 GOP
## These were also done with an ASD and integrating sphere
## bulk samples were collected a year after the spectra
## there's also no stray current correction or calibration
## spectrum for the white reference panel
## Also need to reconcile identifiers between
## Paul's labels and the Fulcrum labels (using Plants app in Fulcrum)

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
# ref.df <- as_tibble(ref)
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
# Hacker2018.spec<-as_spectra(Hacker2018_sg_wide,name_idx = 1,meta_idxs = 1)
# 
# Fulcrum.summary<-read.csv("SummaryData/leaf_spectra.csv")
# Fulcrum.sub<-data.frame(sample.id=Fulcrum.summary$sample_id,
#                         species=Fulcrum.summary$scientific_name,
#                         project=Fulcrum.summary$project,
#                         site=Fulcrum.summary$site_id)
# 
# meta(Hacker2018.spec)$species<-Fulcrum.sub$species[match(meta(Hacker2018.spec)$sample_id,Fulcrum.sub$sample.id)]
# meta(Hacker2018.spec)$project<-Fulcrum.sub$project[match(meta(Hacker2018.spec)$sample_id,Fulcrum.sub$sample.id)]
# 
# ## attach structural traits
# Fulcrum.area<-read.csv("TraitData/LeafAreaWaterSamples/leaf_area_and_water_samples.csv")
# Fulcrum.area$specific_leaf_area_m2_kg<-as.numeric(as.character(Fulcrum.area$specific_leaf_area_m2_kg))
# Fulcrum.area$leaf_dry_matter_content_mg_g<-as.numeric(as.character(Fulcrum.area$leaf_dry_matter_content_mg_g))
# 
# Fulcrum.area.sub<-data.frame(sample_id=Fulcrum.area$parentEventID,
#                              SLA=Fulcrum.area$SLA_m2_kg,
#                              LDMC=Fulcrum.area$LDMC_mg_g)
# 
# ## SLA in units m^2/kg
# ## LDMC in units mg/g
# ## EWT in units cm
# 
# ## note: this doesn't work because the Hacker data
# ## have two different sets of sample IDs
# ## that need to be reconciled!
# ## so I'm not bothering to add these or other traits
# meta(Hacker2018.spec)$SLA<-Fulcrum.area.sub$SLA[match(meta(Hacker2018.spec)$sample_id,Fulcrum.area.sub$sample_id)]
# meta(Hacker2018.spec)$LMA<-1/meta(Hacker2018.spec)$SLA
# meta(Hacker2018.spec)$LDMC<-Fulcrum.area.sub$LDMC[match(meta(Hacker2018.spec)$sample_id,Fulcrum.area.sub$sample_id)]
# meta(Hacker2018.spec)$EWT<-with(meta(Hacker2018.spec),(1/(LDMC/1000)-1)*(1/SLA*0.1))
# 
# meta(Hacker2018.spec)$dataset<-"Hacker2018"

#################################################
## combine

## read in processed data if needed
LOPEX<-readRDS("IndependentValidationData/LOPEX/LOPEX_processed.rds")
ANGERS<-readRDS("IndependentValidationData/ANGERS/ANGERS_processed.rds")
Dessain.spec<-readRDS("IndependentValidationData/Dessain_processed.rds")

## !!!!!!!!!
## might want to add Hacker data

LOPEX_ANGERS<-spectrolab::combine(LOPEX,ANGERS)
all_val_ref<-spectrolab::combine(LOPEX_ANGERS,Dessain.spec)

all_jack_coefs_list_ref<-readRDS("SavedResults/all_jack_coefs_list_ref.rds")
solubles_pred<-apply.coefs(all_jack_coefs_list_ref$sol,val.spec = all_val_ref)
solubles_stat<-t(apply(solubles_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_pred_df<-data.frame(Measured=meta(all_val_ref)$solubles,
                             pred.mean=solubles_stat[,1],
                             pred.low=solubles_stat[,2],
                             pred.high=solubles_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
solubles_all<-with(solubles_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
solubles_upper<-max(solubles_all,na.rm=T)+1
solubles_lower<-min(solubles_all,na.rm=T)-1

hemicellulose_pred<-apply.coefs(all_jack_coefs_list_ref$hemi,val.spec = all_val_ref)
hemicellulose_stat<-t(apply(hemicellulose_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_pred_df<-data.frame(Measured=meta(all_val_ref)$hemicellulose,
                             pred.mean=hemicellulose_stat[,1],
                             pred.low=hemicellulose_stat[,2],
                             pred.high=hemicellulose_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
hemicellulose_all<-with(hemicellulose_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
hemicellulose_upper<-max(hemicellulose_all,na.rm=T)+1
hemicellulose_lower<-min(hemicellulose_all,na.rm=T)-1

cellulose_pred<-apply.coefs(all_jack_coefs_list_ref$cell,val.spec = all_val_ref)
cellulose_stat<-t(apply(cellulose_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_pred_df<-data.frame(Measured=meta(all_val_ref)$cellulose,
                             pred.mean=cellulose_stat[,1],
                             pred.low=cellulose_stat[,2],
                             pred.high=cellulose_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
cellulose_all<-with(cellulose_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
cellulose_upper<-max(cellulose_all,na.rm=T)+1
cellulose_lower<-min(cellulose_all,na.rm=T)-1

lignin_pred<-apply.coefs(all_jack_coefs_list_ref$lign,val.spec = all_val_ref)
lignin_stat<-t(apply(lignin_pred,1,
                     function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_pred_df<-data.frame(Measured=meta(all_val_ref)$lignin,
                           pred.mean=lignin_stat[,1],
                           pred.low=lignin_stat[,2],
                           pred.high=lignin_stat[,3],
                           dataset=meta(all_val_ref)$dataset)
lignin_all<-with(lignin_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
lignin_upper<-max(lignin_all,na.rm=T)+1
lignin_lower<-min(lignin_all,na.rm=T)-1

Nmass_pred<-apply.coefs(all_jack_coefs_list_ref$N,val.spec = all_val_ref)
Nmass_stat<-t(apply(Nmass_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Nmass_pred_df<-data.frame(Measured=meta(all_val_ref)$Nmass,
                             pred.mean=Nmass_stat[,1],
                             pred.low=Nmass_stat[,2],
                             pred.high=Nmass_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
Nmass_all<-with(Nmass_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
Nmass_upper<-max(Nmass_all,na.rm=T)+0.1
Nmass_lower<-min(Nmass_all,na.rm=T)-0.1

Cmass_pred<-apply.coefs(all_jack_coefs_list_ref$C,val.spec = all_val_ref)
Cmass_stat<-t(apply(Cmass_pred,1,
                    function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cmass_pred_df<-data.frame(Measured=meta(all_val_ref)$Cmass,
                          pred.mean=Cmass_stat[,1],
                          pred.low=Cmass_stat[,2],
                          pred.high=Cmass_stat[,3],
                          dataset=meta(all_val_ref)$dataset)
Cmass_all<-with(Cmass_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
Cmass_upper<-max(Cmass_all,na.rm=T)+1
Cmass_lower<-min(Cmass_all,na.rm=T)-1

LMA_pred<-apply.coefs(all_jack_coefs_list_ref$LMA,val.spec = all_val_ref)
LMA_stat<-t(apply(LMA_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA_pred_df<-data.frame(Measured=meta(all_val_ref)$LMA,
                             pred.mean=LMA_stat[,1],
                             pred.low=LMA_stat[,2],
                             pred.high=LMA_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
LMA_all<-with(LMA_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
LMA_upper<-max(LMA_all,na.rm=T)+0.02
LMA_lower<-min(LMA_all,na.rm=T)-0.02

LDMC_pred<-apply.coefs(all_jack_coefs_list_ref$LDMC,val.spec = all_val_ref)
LDMC_stat<-t(apply(LDMC_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LDMC_pred_df<-data.frame(Measured=meta(all_val_ref)$LDMC,
                             pred.mean=LDMC_stat[,1],
                             pred.low=LDMC_stat[,2],
                             pred.high=LDMC_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
LDMC_all<-with(LDMC_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
LDMC_upper<-max(LDMC_all,na.rm=T)+10
LDMC_lower<-min(LDMC_all,na.rm=T)-10

EWT_pred<-apply.coefs(all_jack_coefs_list_ref$EWT,val.spec = all_val_ref)
EWT_stat<-t(apply(EWT_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_pred_df<-data.frame(Measured=meta(all_val_ref)$EWT,
                             pred.mean=EWT_stat[,1],
                             pred.low=EWT_stat[,2],
                             pred.high=EWT_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
EWT_all<-with(EWT_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
EWT_upper<-max(EWT_all,na.rm=T)+0.03
EWT_lower<-min(EWT_all,na.rm=T)-0.03

chlA_pred<-apply.coefs(all_jack_coefs_list_ref$chlA,val.spec = all_val_ref)
chlA_stat<-t(apply(chlA_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_pred_df<-data.frame(Measured=meta(all_val_ref)$chlA,
                             pred.mean=chlA_stat[,1],
                             pred.low=chlA_stat[,2],
                             pred.high=chlA_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
chlA_all<-with(chlA_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
chlA_upper<-max(chlA_all,na.rm=T)+0.5
chlA_lower<-min(chlA_all,na.rm=T)-0.5

chlB_pred<-apply.coefs(all_jack_coefs_list_ref$chlB,val.spec = all_val_ref)
chlB_stat<-t(apply(chlB_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_pred_df<-data.frame(Measured=meta(all_val_ref)$chlB,
                             pred.mean=chlB_stat[,1],
                             pred.low=chlB_stat[,2],
                             pred.high=chlB_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
chlB_all<-with(chlB_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
chlB_upper<-max(chlB_all,na.rm=T)+0.3
chlB_lower<-min(chlB_all,na.rm=T)-0.3

car_pred<-apply.coefs(all_jack_coefs_list_ref$car,val.spec = all_val_ref)
car_stat<-t(apply(car_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_pred_df<-data.frame(Measured=meta(all_val_ref)$car,
                             pred.mean=car_stat[,1],
                             pred.low=car_stat[,2],
                             pred.high=car_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
car_all<-with(car_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
car_upper<-max(car_all,na.rm=T)+0.2
car_lower<-min(car_all,na.rm=T)-0.2

Al_pred<-apply.coefs(all_jack_coefs_list_ref$Al,val.spec = all_val_ref)
Al_stat<-t(apply(Al_pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_pred_df<-data.frame(Measured=meta(all_val_ref)$Al_mass,
                             pred.mean=Al_stat[,1],
                             pred.low=Al_stat[,2],
                             pred.high=Al_stat[,3],
                             dataset=meta(all_val_ref)$dataset)
Al_all<-with(Al_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
Al_upper<-max(Al_all,na.rm=T)+0.02
Al_lower<-min(Al_all,na.rm=T)-0.02

Ca_pred<-apply.coefs(all_jack_coefs_list_ref$Ca,val.spec = all_val_ref)
Ca_stat<-t(apply(Ca_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_pred_df<-data.frame(Measured=meta(all_val_ref)$Ca_mass,
                       pred.mean=Ca_stat[,1],
                       pred.low=Ca_stat[,2],
                       pred.high=Ca_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
Ca_all<-with(Ca_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
Ca_upper<-max(Ca_all,na.rm=T)+1
Ca_lower<-min(Ca_all,na.rm=T)-1

Cu_pred<-apply.coefs(all_jack_coefs_list_ref$Cu,val.spec = all_val_ref)
Cu_stat<-t(apply(Cu_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_pred_df<-data.frame(Measured=meta(all_val_ref)$Cu_mass,
                       pred.mean=Cu_stat[,1],
                       pred.low=Cu_stat[,2],
                       pred.high=Cu_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
Cu_all<-with(Cu_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
Cu_upper<-max(Cu_all,na.rm=T)+0.005
Cu_lower<-min(Cu_all,na.rm=T)-0.005

Fe_pred<-apply.coefs(all_jack_coefs_list_ref$Fe,val.spec = all_val_ref)
Fe_stat<-t(apply(Fe_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_pred_df<-data.frame(Measured=meta(all_val_ref)$Fe_mass,
                       pred.mean=Fe_stat[,1],
                       pred.low=Fe_stat[,2],
                       pred.high=Fe_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
Fe_all<-with(Fe_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
Fe_upper<-max(Fe_all,na.rm=T)+0.02
Fe_lower<-min(Fe_all,na.rm=T)-0.02

K_pred<-apply.coefs(all_jack_coefs_list_ref$K,val.spec = all_val_ref)
K_stat<-t(apply(K_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_pred_df<-data.frame(Measured=meta(all_val_ref)$K_mass,
                       pred.mean=K_stat[,1],
                       pred.low=K_stat[,2],
                       pred.high=K_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
K_all<-with(K_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
K_upper<-max(K_all,na.rm=T)+1
K_lower<-min(K_all,na.rm=T)-1

Mg_pred<-apply.coefs(all_jack_coefs_list_ref$Mg,val.spec = all_val_ref)
Mg_stat<-t(apply(Mg_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_pred_df<-data.frame(Measured=meta(all_val_ref)$Mg_mass,
                       pred.mean=Mg_stat[,1],
                       pred.low=Mg_stat[,2],
                       pred.high=Mg_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
Mg_all<-with(Mg_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
Mg_upper<-max(Mg_all,na.rm=T)+0.3
Mg_lower<-min(Mg_all,na.rm=T)-0.3

Mn_pred<-apply.coefs(all_jack_coefs_list_ref$Mn,val.spec = all_val_ref)
Mn_stat<-t(apply(Mn_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_pred_df<-data.frame(Measured=meta(all_val_ref)$Mn_mass,
                       pred.mean=Mn_stat[,1],
                       pred.low=Mn_stat[,2],
                       pred.high=Mn_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
Mn_all<-with(Mn_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
Mn_upper<-max(Mn_all,na.rm=T)+0.2
Mn_lower<-min(Mn_all,na.rm=T)-0.2

Na_pred<-apply.coefs(all_jack_coefs_list_ref$Na,val.spec = all_val_ref)
Na_stat<-t(apply(Na_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_pred_df<-data.frame(Measured=meta(all_val_ref)$Na_mass,
                       pred.mean=Na_stat[,1],
                       pred.low=Na_stat[,2],
                       pred.high=Na_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
Na_all<-with(Na_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
Na_upper<-max(Na_all,na.rm=T)+0.5
Na_lower<-min(Na_all,na.rm=T)-0.5

P_pred<-apply.coefs(all_jack_coefs_list_ref$P,val.spec = all_val_ref)
P_stat<-t(apply(P_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_pred_df<-data.frame(Measured=meta(all_val_ref)$P_mass,
                       pred.mean=P_stat[,1],
                       pred.low=P_stat[,2],
                       pred.high=P_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
P_all<-with(P_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
P_upper<-max(P_all,na.rm=T)+0.3
P_lower<-min(P_all,na.rm=T)-0.3

Zn_pred<-apply.coefs(all_jack_coefs_list_ref$Zn,val.spec = all_val_ref)
Zn_stat<-t(apply(Zn_pred,1,
                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_pred_df<-data.frame(Measured=meta(all_val_ref)$Zn_mass,
                       pred.mean=Zn_stat[,1],
                       pred.low=Zn_stat[,2],
                       pred.high=Zn_stat[,3],
                       dataset=meta(all_val_ref)$dataset)
Zn_all<-with(Zn_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
Zn_upper<-max(Zn_all,na.rm=T)+0.05
Zn_lower<-min(Zn_all,na.rm=T)-0.05

######################################
## plotting

colorBlind  <- c("#E69F00","#009E73","#F0E442","#56B4E9",
                 "#0072B2","#CC79A7","#D55E00","#999999")

solubles_ind_val<-ggplot(data=solubles_pred_df,
                         aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),
                  ylim=c(solubles_lower,solubles_upper))+
  scale_color_manual(values=colorBlind)

hemicellulose_ind_val<-ggplot(data=hemicellulose_pred_df,
                              aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),
                  ylim=c(hemicellulose_lower,hemicellulose_upper))+
  scale_color_manual(values=colorBlind)

cellulose_ind_val<-ggplot(data=cellulose_pred_df,
                          aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),
                  ylim=c(cellulose_lower,cellulose_upper))+
  scale_color_manual(values=colorBlind)

lignin_ind_val<-ggplot(data=lignin_pred_df,
                       aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),
                  ylim=c(lignin_lower,lignin_upper))+
  scale_color_manual(values=colorBlind)

Nmass_ind_val<-ggplot(data=Nmass_pred_df,
                      aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  coord_cartesian(xlim=c(Nmass_lower,Nmass_upper),
                  ylim=c(Nmass_lower,Nmass_upper))+
  scale_color_manual(values=colorBlind)

Cmass_ind_val<-ggplot(data=Cmass_pred_df,
                      aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  coord_cartesian(xlim=c(Cmass_lower,Cmass_upper),
                  ylim=c(Cmass_lower,Cmass_upper))+
  scale_color_manual(values=colorBlind)

LMA_ind_val<-ggplot(data=LMA_pred_df,
                    aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),
                  ylim=c(LMA_lower,LMA_upper))+
  scale_color_manual(values=colorBlind)

LDMC_ind_val<-ggplot(data=LDMC_pred_df,
                     aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),
                  ylim=c(LDMC_lower,LDMC_upper))+
  scale_color_manual(values=colorBlind)

EWT_ind_val<-ggplot(data=EWT_pred_df,
                    aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  scale_color_manual(values=colorBlind)

chlA_ind_val<-ggplot(data=chlA_pred_df,
                     aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),
                  ylim=c(chlA_lower,chlA_upper))+
  scale_color_manual(values=colorBlind)

chlB_ind_val<-ggplot(data=chlB_pred_df,
                     aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),
                  ylim=c(chlB_lower,chlB_upper))+
  scale_color_manual(values=colorBlind)

car_ind_val<-ggplot(data=car_pred_df,
                    aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"))+
  coord_cartesian(xlim=c(car_lower,car_upper),
                  ylim=c(car_lower,car_upper))+
  scale_color_manual(values=colorBlind)

Al_ind_val<-ggplot(data=Al_pred_df,
                   aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Al (mg g"^-1*")"),
       x=expression("Predicted Al (mg g"^-1*")"))+
  coord_cartesian(xlim=c(Al_lower,Al_upper),
                  ylim=c(Al_lower,Al_upper))+
  scale_color_manual(values=colorBlind)

Ca_ind_val<-ggplot(data=Ca_pred_df,
                   aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Ca (mg g"^-1*")"),
       x=expression("Predicted Ca (mg g"^-1*")"))+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),
                  ylim=c(Ca_lower,Ca_upper))+
  scale_color_manual(values=colorBlind)

Cu_ind_val<-ggplot(data=Cu_pred_df,
                   aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Cu (mg g"^-1*")"),
       x=expression("Predicted Cu (mg g"^-1*")"))+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),
                  ylim=c(Cu_lower,Cu_upper))+
  scale_color_manual(values=colorBlind)

Fe_ind_val<-ggplot(data=Fe_pred_df,
                   aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Fe (mg g"^-1*")"),
       x=expression("Predicted Fe (mg g"^-1*")"))+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),
                  ylim=c(Fe_lower,Fe_upper))+
  scale_color_manual(values=colorBlind)

K_ind_val<-ggplot(data=K_pred_df,
                  aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  coord_cartesian(xlim=c(K_lower,K_upper),
                  ylim=c(K_lower,K_upper))+
  scale_color_manual(values=colorBlind)

Mg_ind_val<-ggplot(data=Mg_pred_df,
                   aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Mg (mg g"^-1*")"),
       x=expression("Predicted Mg (mg g"^-1*")"))+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),
                  ylim=c(Mg_lower,Mg_upper))+
  scale_color_manual(values=colorBlind)

Mn_ind_val<-ggplot(data=Mn_pred_df,
                   aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Mn (mg g"^-1*")"),
       x=expression("Predicted Mn (mg g"^-1*")"))+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),
                  ylim=c(Mn_lower,Mn_upper))+
  scale_color_manual(values=colorBlind)

Na_ind_val<-ggplot(data=Na_pred_df,
                   aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Na (mg g"^-1*")"),
       x=expression("Predicted Na (mg g"^-1*")"))+
  coord_cartesian(xlim=c(Na_lower,Na_upper),
                  ylim=c(Na_lower,Na_upper))+
  scale_color_manual(values=colorBlind)

P_ind_val<-ggplot(data=P_pred_df,
                  aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured P (mg g"^-1*")"),
       x=expression("Predicted P (mg g"^-1*")"))+
  coord_cartesian(xlim=c(P_lower,P_upper),
                  ylim=c(P_lower,P_upper))+
  scale_color_manual(values=colorBlind)

Zn_ind_val<-ggplot(data=Zn_pred_df,
                   aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+
  theme_bw()+
  theme(text = element_text(size=25))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Zn (mg g"^-1*")"),
       x=expression("Predicted Zn (mg g"^-1*")"))+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),
                  ylim=c(Zn_lower,Zn_upper))+
  scale_color_manual(values=colorBlind)

ind_val_list<-list(solubles_mass=solubles_pred_df,
                   hemicellulose_mass=hemicellulose_pred_df,
                   cellulose_mass=cellulose_pred_df,
                   lignin_mass=lignin_pred_df,
                   Nmass=Nmass_pred_df,
                   Cmass=Cmass_pred_df,
                   LMA=LMA_pred_df,
                   LDMC=LDMC_pred_df,
                   EWT=EWT_pred_df,
                   chlA_mass=chlA_pred_df,
                   chlB_mass=chlB_pred_df,
                   car_mass=car_pred_df,
                   Al_mass=Al_pred_df,
                   Ca_mass=Ca_pred_df,
                   Cu_mass=Cu_pred_df,
                   Fe_mass=Fe_pred_df,
                   K_mass=K_pred_df,
                   Mg_mass=Mg_pred_df,
                   Mn_mass=Mn_pred_df,
                   Na_mass=Na_pred_df,
                   P_mass=P_pred_df,
                   Zn_mass=Zn_pred_df)
saveRDS(ind_val_list,"SavedResults/ind_val_list.rds")

pdf("Images/ind_val_plots1.pdf",width = 18,height = 22.5, onefile=F)
ggarrange(plotlist=list(solubles_ind_val,hemicellulose_ind_val,
                        cellulose_ind_val,lignin_ind_val,
                        Nmass_ind_val,Cmass_ind_val,
                        LMA_ind_val,LDMC_ind_val,EWT_ind_val,
                        chlA_ind_val,chlB_ind_val,car_ind_val),
          common.legend = T,legend = "bottom",
          nrow=4,ncol=3)
dev.off()

pdf("Images/ind_val_plots2.pdf",width = 18,height = 22.5, onefile=F)
ggarrange(plotlist=list(Al_ind_val,Ca_ind_val,Cu_ind_val,
                        Fe_ind_val,K_ind_val,Mg_ind_val,
                        Mn_ind_val,Na_ind_val,P_ind_val,
                        Zn_ind_val),
          common.legend = T,legend = "bottom",
          nrow=4,ncol=3)
dev.off()
