setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels/")

library(spectrolab)
library(ggplot2)
library(reshape2)
library(stringr)

all.spec<-readRDS("ProcessedSpectra/all_spectra.rds")

##############################################
## plot quantiles

# all.quantiles<-quantile(all.spec,probs=c(0.025,0.25,0.5,0.75,0.975))
# all.CV<-apply(all.spec,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
# 
# all.spec.plot<-ggplot()+
#   geom_line(aes(x=400:2400,y=as.matrix(all.quantiles)[3,]),size=1)+
#   geom_line(aes(x=400:2400,y=as.matrix(all.quantiles)[2,]),size=1,linetype="dashed")+
#   geom_line(aes(x=400:2400,y=as.matrix(all.quantiles)[4,]),size=1,linetype="dashed")+
#   geom_line(aes(x=400:2400,y=as.matrix(all.quantiles)[1,]),size=1,linetype="dotted")+
#   geom_line(aes(x=400:2400,y=as.matrix(all.quantiles)[5,]),size=1,linetype="dotted")+
#   geom_line(aes(x=400:2400,y=all.CV),size=1,color="red")+
#   theme_bw()+
#   theme(text = element_text(size=20),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   labs(x="Wavelength",y="Reflectance (or CV)")+lims(y=c(0,1))+
#   ggtitle("All CABO spectra")

##############################################
## plot projects

# meta(all.spec)$project<-factor(meta(all.spec)$project)
# all.spec.df<-as.data.frame(all.spec)
# all.spec.df$sample_name<-NULL
# 
# Girard.samples<-which(all.spec.df$project=="2018-Girard-MSc-UdeM")
# all.spec.Girard<-all.spec.df[Girard.samples,]
# Quercus.samples<-grep("Quercus",all.spec.df$species)
# all.spec.Quercus<-all.spec.df[Quercus.samples,]
# 
# all.spec.sample<-all.spec.df[sample(1:nrow(all.spec.df),200),]
# 
# all.spec.Girard.long<-melt(data=all.spec.Girard,
#                     id.vars=c("sample_id","species","project"),
#                     variable.name="wavelength",
#                     value.name="reflectance")
# all.spec.Girard.long$wavelength<-as.numeric(as.character(all.spec.Girard.long$wavelength))
# 
# all.spec.Quercus.long<-melt(data=all.spec.Quercus,
#                            id.vars=c("sample_id","species","project"),
#                            variable.name="wavelength",
#                            value.name="reflectance")
# all.spec.Quercus.long$wavelength<-as.numeric(as.character(all.spec.Quercus.long$wavelength))
# 
# all.spec.long<-melt(data=all.spec.sample,
#                            id.vars=c("sample_id","species","project"),
#                            variable.name="wavelength",
#                            value.name="reflectance")
# all.spec.long$wavelength<-as.numeric(as.character(all.spec.long$wavelength))
# 
# Quercus_plot<-ggplot()+
#   geom_line(data=all.spec.long,
#             aes(x=wavelength,y=reflectance,group=sample_id),color="gray")+
#   geom_line(data=all.spec.Quercus.long,
#             aes(x=wavelength,y=reflectance,group=sample_id),color="orange",alpha=0.2)+
#   theme_bw()+
#   theme(text = element_text(size=20),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   labs(x="Wavelength",y="Reflectance")+lims(y=c(0,1))+
#   ggtitle("Quercus spectra")
# 
# Girard_plot<-ggplot()+
#   geom_line(data=all.spec.long,
#             aes(x=wavelength,y=reflectance,group=sample_id),color="gray")+
#   geom_line(data=all.spec.Girard.long,
#             aes(x=wavelength,y=reflectance,group=sample_id),color="blue",alpha=0.2)+
#   theme_bw()+
#   theme(text = element_text(size=20),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   labs(x="Wavelength",y="Reflectance")+lims(y=c(0,1))+
#   ggtitle("Girard spectra")

#############################################
## growth form

vascan<-read.csv("SummaryData/vascan.csv")
meta(all.spec)$growth.form<-vascan$growth_form[match(meta(all.spec)$species,vascan$scientific_name)]
meta(all.spec)$family<-vascan$family[match(meta(all.spec)$species,vascan$scientific_name)]

meta(all.spec)$growth.form[meta(all.spec)$species=="Acer pensylvanicum Linnaeus"]<-"tree"
meta(all.spec)$growth.form[meta(all.spec)$species=="Betula populifolia Marshall"]<-"tree"
meta(all.spec)$growth.form[meta(all.spec)$species=="Sorbus decora (Sargent) C.K. Schneider"]<-"tree"
meta(all.spec)$growth.form[meta(all.spec)$species=="Sorbus americana Marshall"]<-"tree"
meta(all.spec)$growth.form[meta(all.spec)$species=="Prunus pensylvanica Linnaeus f."]<-"tree"
meta(all.spec)$growth.form[meta(all.spec)$species=="Frangula alnus Miller"]<-"shrub"
meta(all.spec)$growth.form[meta(all.spec)$species=="Rhamnus cathartica Linnaeus"]<-"shrub"
meta(all.spec)$growth.form[meta(all.spec)$species=="Prunus nigra Aiton"]<-"tree"
meta(all.spec)$growth.form[meta(all.spec)$species=="Salix Linnaeus"]<-"shrub"
## the next four, I'm not too sure about the tree/shrub assignment
## making a best guess based on max height and project
meta(all.spec)$growth.form[meta(all.spec)$species=="Prunus virginiana Linnaeus"]<-"shrub"
meta(all.spec)$growth.form[meta(all.spec)$species=="Staphylea trifolia Linnaeus"]<-"shrub"
meta(all.spec)$growth.form[meta(all.spec)$species=="Oemleria cerasiformis (Torrey & A. Gray ex Hooker & Arnott) J.W. Landon"]<-"shrub"
meta(all.spec)$growth.form[meta(all.spec)$species=="Crataegus monogyna Jacquin"]<-"shrub"
meta(all.spec)$growth.form[meta(all.spec)$species=="Agonis flexuosa (Willd.) Sweet"]<-"tree"

meta(all.spec)$functional.group<-as.character(meta(all.spec)$growth.form)
meta(all.spec)$functional.group[meta(all.spec)$family %in% c("Poaceae","Cyperaceae","Juncaceae")]<-"graminoid"
meta(all.spec)$functional.group[meta(all.spec)$family=="Fabaceae"]<-"legume"
meta(all.spec)$functional.group[meta(all.spec)$family %in% c("Pinaceae","Cupressaceae")]<-"conifer"

##############################################
## read leaf area/water traits

Fulcrum.area<-read.csv("TraitData/LeafAreaWaterSamples/leaf_area_and_water_samples.csv")
Dessain.area<-read.csv("TraitData/LeafAreaWaterSamples/SLA_data_Aurelie_Dessain - Lab_data.csv")

Dessain.area$SLA_m2_kg<-as.numeric(as.character(Dessain.area$SLA_m2_kg))
Dessain.area$LDMC_mg_g<-as.numeric(as.character(Dessain.area$LDMC_mg_g))

Fulcrum.area.sub<-data.frame(sample_id=Fulcrum.area$sample_id,
                             SLA=Fulcrum.area$specific_leaf_area_m2_kg,
                             LDMC=Fulcrum.area$leaf_dry_matter_content_mg_g)
Dessain.area.sub<-data.frame(sample_id=Dessain.area$parentEventID,
                             SLA=Dessain.area$SLA_m2_kg,
                             LDMC=Dessain.area$LDMC_mg_g)

all.area<-do.call(rbind,list(Fulcrum.area.sub,Dessain.area.sub))

## SLA in units m^2/kg
## LDMC in units mg/g
meta(all.spec)$SLA<-all.area$SLA[match(meta(all.spec)$sample_id,all.area$sample_id)]
meta(all.spec)$LDMC<-all.area$LDMC[match(meta(all.spec)$sample_id,all.area$sample_id)]

## remove bad values, based on notes in Fulcrum data
meta(all.spec)$LDMC[meta(all.spec)$sample_id %in% c("10290262","10966273","13404937",
                                                    "38530951","41809826","42944395",
                                                    "43711631","43713406","43717241",
                                                    "43718957","43718799","43712000","43721033",
                                                    "44060633","44142362","44148646",
                                                    "44683516","45108236")]<-NA

meta(all.spec)$SLA[meta(all.spec)$sample_id %in% c("13404937","44227362",
                                                   "44683516","45108236","44142362")]<-NA

## EWT in units cm
meta(all.spec)$EWT<-with(meta(all.spec),(1/(LDMC/1000)-1)*(1/SLA*0.1))

##############################################
## read C/N

Dessain.CN<-read.csv("TraitData/CNSamples/2017_Dessain_MSc_CN_data_total.csv")
Dessain.CN.sub<-data.frame(sample_id=Dessain.CN$Sample_id,
                           Cmass=Dessain.CN$C.....,
                           Nmass=Dessain.CN$N....)

BeauchampRioux.CN<-read.csv("TraitData/CNSamples/2018_BeauchampRioux_MSc_UdeM_CN_data_total.csv")
BeauchampRioux.CN<-BeauchampRioux.CN[-which(is.na(BeauchampRioux.CN$Sample_id)),]
BeauchampRioux.CN.sub<-data.frame(sample_id=BeauchampRioux.CN$Sample_id,
                           Cmass=BeauchampRioux.CN$C.....,
                           Nmass=BeauchampRioux.CN$N....)

## note: 20 Pinus samples missing for now
Blanchard.CN<-read.csv("TraitData/CNSamples/2019_Blanchard_CN_data_total.csv")
Blanchard.CN<-Blanchard.CN[-which(is.na(Blanchard.CN$Sample_id)),]
Blanchard.CN.sub<-data.frame(sample_id=Blanchard.CN$Sample_id,
                                  Cmass=Blanchard.CN$C.....,
                                  Nmass=Blanchard.CN$N....)

Boucherville2018.CN<-read.csv("TraitData/CNSamples/2018_Boucherville_CN_data_total.csv")
Boucherville2018.CN<-Boucherville2018.CN[-which(is.na(Boucherville2018.CN$Sample_id)),]
Boucherville2018.CN.sub<-data.frame(sample_id=Boucherville2018.CN$Sample_id,
                           Cmass=Boucherville2018.CN$C...,
                           Nmass=Boucherville2018.CN$N....)

Boucherville2019.CN<-read.csv("TraitData/CNSamples/2019-Boucherville_CN_data_total.csv")
Boucherville2019.CN.sub<-data.frame(sample_id=Boucherville2019.CN$Sample_id,
                                    Cmass=Boucherville2019.CN$C...,
                                    Nmass=Boucherville2019.CN$N....)

Crofts.CN<-read.csv("TraitData/CNSamples/2019-Crofts_CN_data_total.csv")
Crofts.CN.sub<-data.frame(sample_id=Crofts.CN$Numéro.d.échantillon,
                              Cmass=Crofts.CN$C.....,
                              Nmass=Crofts.CN$N....)

Girard.CN<-read.csv("TraitData/CNSamples/2018_Girard_CN_data_total.csv")
Girard.CN.sub<-data.frame(sample_id=Girard.CN$sample_id,
                              Cmass=Girard.CN$c_perc,
                              Nmass=Girard.CN$n_perc)

Hacker2018.CN<-read.csv("TraitData/CNSamples/2018_Hacker_CGOP2018_CN_data_total.csv")
Hacker2018.CN<-Hacker2018.CN[-which(is.na(Hacker2018.CN$Sample_id)),]
Hacker2018.CN.sub<-data.frame(sample_id=Hacker2018.CN$Sample_id,
                                    Cmass=Hacker2018.CN$C...,
                                    Nmass=Hacker2018.CN$N....)

Hacker2019.CN<-read.csv("TraitData/CNSamples/2018_Hacker_CGOP2019_CN_data_total.csv")
Hacker2019.CN<-Hacker2019.CN[-which(is.na(Hacker2019.CN$Sample_id)),]
Hacker2019.CN.sub<-data.frame(sample_id=Hacker2019.CN$Sample_id,
                              Cmass=Hacker2019.CN$C...,
                              Nmass=Hacker2019.CN$N....)

PhragmitesTemporal.CN<-read.csv("TraitData/CNSamples/2019-Phragmites-temporal_total_CN.csv")
PhragmitesTemporal.CN<-PhragmitesTemporal.CN[-which(is.na(PhragmitesTemporal.CN$bulk)),]
PhragmitesTemporal.CN.sub<-data.frame(sample_id=PhragmitesTemporal.CN$bulk,
                              Cmass=PhragmitesTemporal.CN$C...,
                              Nmass=PhragmitesTemporal.CN$N....)

CABOGeneralOther.CN<-read.csv("TraitData/CNSamples/CABO_General_CN_data_total.csv")
CABOGeneralOther.CN.sub<-data.frame(sample_id=CABOGeneralOther.CN$Sample_id,
                              Cmass=CABOGeneralOther.CN$C...,
                              Nmass=CABOGeneralOther.CN$N....)

CABOGeneral2019.CN<-read.csv("TraitData/CNSamples/2019_CABO_General_CN_data_total.csv")
CABOGeneral2019.CN.sub<-data.frame(sample_id=CABOGeneral2019.CN$Sample_id,
                                    Cmass=CABOGeneral2019.CN$C...,
                                    Nmass=CABOGeneral2019.CN$N....)

Warren.CN<-read.csv("TraitData/CNSamples/SWA_Warren_CN_data_total.csv")
Warren.CN<-Warren.CN[-which(is.na(Warren.CN$Sample_id)),]
Warren.CN.sub<-data.frame(sample_id=Warren.CN$Sample_id,
                              Cmass=Warren.CN$C...,
                              Nmass=Warren.CN$N....)

all.CN<-do.call(rbind,list(BeauchampRioux.CN.sub,Blanchard.CN.sub,Boucherville2018.CN.sub,
                           Boucherville2019.CN.sub,CABOGeneralOther.CN.sub,
                           CABOGeneral2019.CN.sub,Crofts.CN.sub,
                           Dessain.CN.sub,Girard.CN.sub,Hacker2018.CN.sub,
                           Hacker2019.CN.sub,PhragmitesTemporal.CN.sub,
                           Warren.CN.sub))

meta(all.spec)$Cmass<-all.CN$Cmass[match(meta(all.spec)$sample_id,all.CN$sample_id)]
meta(all.spec)$Nmass<-all.CN$Nmass[match(meta(all.spec)$sample_id,all.CN$sample_id)]

meta(all.spec)$Cmass[meta(all.spec)$sample_id=="38707392"]<-NA
meta(all.spec)$Nmass[meta(all.spec)$sample_id=="38707392"]<-NA

#############################################
## carbon fractions

CFractions<-read.csv("TraitData/CarbonFractions/carbon_fractions_bags.csv")
CFmatch<-match(meta(all.spec)$sample_id,CFractions$bottle_id)
meta(all.spec)$NDFmass<-CFractions$ndf_perc[CFmatch]
meta(all.spec)$ADFmass<-CFractions$adf_perc[CFmatch]
meta(all.spec)$ADLmass<-CFractions$adl_perc[CFmatch]
meta(all.spec)$solubles_mass<-CFractions$soluble_perc[CFmatch]
meta(all.spec)$hemicellulose_mass<-CFractions$hemicellulose_perc[CFmatch]
meta(all.spec)$cellulose_mass<-CFractions$cellulose_perc[CFmatch]
meta(all.spec)$lignin_mass<-CFractions$lignin_perc[CFmatch]

## removing species that have exudates that affect carbon fraction measurements
## there may also be a problem with "Polystichum munitum (Kaulfuss) C. Presl" but not confirmed
NDF_bad<-c("Ulmus rubra Muhlenberg","Acer platanoides Linnaeus",
           "Ulmus americana Linnaeus","Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")
meta(all.spec)$NDFmass[which(meta(all.spec)$species %in% NDF_bad)]<-NA
meta(all.spec)$solubles_mass[which(meta(all.spec)$species %in% NDF_bad)]<-NA
meta(all.spec)$hemicellulose_mass[which(meta(all.spec)$species %in% NDF_bad)]<-NA

meta(all.spec)$ADFmass[which(meta(all.spec)$species=="Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")]<-NA
meta(all.spec)$ADLmass[which(meta(all.spec)$species=="Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")]<-NA
meta(all.spec)$cellulose_mass[which(meta(all.spec)$species=="Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")]<-NA
meta(all.spec)$lignin_mass[which(meta(all.spec)$species=="Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")]<-NA

#############################################
## read pigments

Dessain.pigments<-read.csv("TraitData/Pigments/Aurelie_pigments_valeurs_brutes - Analyses_Aurelie.csv")
Dessain.pigments<-Dessain.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]

BeauchampRioux.pigments<-read.csv("TraitData/Pigments/RBR_pigments_valeurs_brutes - valeurs.csv")
BeauchampRioux.pigments<-BeauchampRioux.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides_mg_g","Notes")]
colnames(BeauchampRioux.pigments)[colnames(BeauchampRioux.pigments)=="carotenoides_mg_g"]<-"carotenoides._mg_g"

Blanchard.pigments<-read.csv("TraitData/Pigments/Blanchard_pigments_valeurs_brutes.csv")
Blanchard.pigments<-Blanchard.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]

Boucherville2018.pigments<-read.csv("TraitData/Pigments/Boucherville_pigments_valeurs_brutes.csv")
Boucherville2018.pigments<-Boucherville2018.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]

Boucherville2019.pigments<-read.csv("TraitData/Pigments/2019_Boucherville_pigments_valeurs_brutes.csv")
Boucherville2019.pigments<-Boucherville2019.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]
Boucherville2019.pigments<-Boucherville2019.pigments[-which(Boucherville2019.pigments$sample_id==""),]

Crofts.pigments<-read.csv("TraitData/Pigments/2019_Crofts_pigments_valeurs_brutes.csv")
Crofts.pigments<-Crofts.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]
Crofts.pigments<-Crofts.pigments[-which(Crofts.pigments$sample_id==""),]

Girard.pigments<-read.csv("TraitData/Pigments/Alizée_pigments_valeurs_brutes.csv")
Girard.pigments<-Girard.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides_mg_g","Notes")]
colnames(Girard.pigments)[colnames(Girard.pigments)=="carotenoides_mg_g"]<-"carotenoides._mg_g"
Girard.pigments<-Girard.pigments[-which(Girard.pigments$sample_id==""),]

Hacker2018.pigments<-read.csv("TraitData/Pigments/Hacker_2018_pigments_valeurs_brutes.csv")
Hacker2018.pigments<-Hacker2018.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]
Hacker2018.pigments<-Hacker2018.pigments[-which(Hacker2018.pigments$sample_id==""),]

Hacker2019.pigments<-read.csv("TraitData/Pigments/Hacker_2019_pigments_valeurs_brutes.csv")
Hacker2019.pigments<-Hacker2019.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]
Hacker2019.pigments<-Hacker2019.pigments[-which(Hacker2019.pigments$sample_id==""),]

PhragmitesTemporal.pigments<-read.csv("TraitData/Pigments/2019_Phragmites-temporal_pigments_valeurs_brutes.csv")
PhragmitesTemporal.pigments<-PhragmitesTemporal.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]
PhragmitesTemporal.pigments<-PhragmitesTemporal.pigments[-which(PhragmitesTemporal.pigments$sample_id==""),]

CABOGeneral.pigments<-read.csv("TraitData/Pigments/CABO_general_pigments_valeurs_brutes.csv")
CABOGeneral.pigments<-CABOGeneral.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]

CABOGeneral2019.pigments<-read.csv("TraitData/Pigments/2019_CABO_general_pigments_valeurs_brutes.csv")
CABOGeneral2019.pigments<-CABOGeneral2019.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]
CABOGeneral2019.pigments<-CABOGeneral2019.pigments[-which(CABOGeneral2019.pigments$sample_id==""),]

Warren.pigments<-read.csv("TraitData/Pigments/Warren_pigments_valeurs_brutes - Valeurs.csv")
Warren.pigments<-Warren.pigments[,c("sample_id","chlA_mg_g","chlB_mg_g","carotenoides._mg_g","Notes")]

all.pigments<-do.call(rbind,list(BeauchampRioux.pigments,Blanchard.pigments,Boucherville2018.pigments,
                           Boucherville2019.pigments,CABOGeneral.pigments,
                           CABOGeneral2019.pigments,Crofts.pigments,
                           Dessain.pigments,Girard.pigments,Hacker2018.pigments,
                           Hacker2019.pigments,PhragmitesTemporal.pigments,
                           Warren.pigments))

all.pigments$chlA_mg_g[all.pigments$sample_id==13404937]<-NA
all.pigments$chlB_mg_g[all.pigments$sample_id==13404937]<-NA
all.pigments$carotenoides._mg_g[all.pigments$sample_id==13404937]<-NA

## remove records who notes contain "refaire" or "refait" since these were all redone
all.pigments<-all.pigments[-which(str_detect(string = all.pigments$Notes,pattern = "efai")),]

meta(all.spec)$chlA_fresh<-as.numeric(as.character(all.pigments$chlA_mg_g[match(meta(all.spec)$sample_id,all.pigments$sample_id)]))
meta(all.spec)$chlB_fresh<-as.numeric(as.character(all.pigments$chlB_mg_g[match(meta(all.spec)$sample_id,all.pigments$sample_id)]))
meta(all.spec)$car_fresh<-as.numeric(as.character(all.pigments$carotenoides._mg_g[match(meta(all.spec)$sample_id,all.pigments$sample_id)]))

#############################################
## area-normalized chemical traits and
## normalization-independent chemical traits sensu Osnas

## LMA in kg/m^2
meta(all.spec)$LMA<-1/meta(all.spec)$SLA

## area basis, in g/cm^2
meta(all.spec)$Narea<-meta(all.spec)$Nmass*meta(all.spec)$LMA/1000
Narea_norm<-lm(log(Narea)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$Nnorm<-resid(Narea_norm)

## area basis, in g/cm^2
meta(all.spec)$Carea<-meta(all.spec)$Cmass*meta(all.spec)$LMA/1000
Carea_norm<-lm(log(Carea)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$Cnorm<-resid(Carea_norm)

## area basis, in g/cm^2
meta(all.spec)$solubles_area<-meta(all.spec)$solubles_mass*meta(all.spec)$LMA/1000
solubles_area_norm<-lm(log(solubles_area)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$solubles_norm<-resid(solubles_area_norm)

## area basis, in g/cm^2
meta(all.spec)$hemicellulose_area<-meta(all.spec)$hemicellulose_mass*meta(all.spec)$LMA/1000
hemicellulose_area_norm<-lm(log(hemicellulose_area)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$hemicellulose_norm<-resid(hemicellulose_area_norm)

## area basis, in g/cm^2
meta(all.spec)$cellulose_area<-meta(all.spec)$cellulose_mass*meta(all.spec)$LMA/1000
cellulose_area_norm<-lm(log(cellulose_area)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$cellulose_norm<-resid(cellulose_area_norm)

## area basis, in g/cm^2
meta(all.spec)$lignin_area<-meta(all.spec)$lignin_mass*meta(all.spec)$LMA/1000
lignin_area_norm<-lm(log(lignin_area)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$lignin_norm<-resid(lignin_area_norm)

## first to dry mass basis
meta(all.spec)$chlA_dry<-meta(all.spec)$chlA_fresh*1000/meta(all.spec)$LDMC
## then to area basis
meta(all.spec)$chlA_area<-meta(all.spec)$chlA_dry*meta(all.spec)$LMA/1000
chlA_area_norm<-lm(log(chlA_area)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$chlA_norm<-resid(chlA_area_norm)

## first to dry mass basis
meta(all.spec)$chlB_dry<-meta(all.spec)$chlB_fresh*1000/meta(all.spec)$LDMC
## then to area basis
meta(all.spec)$chlB_area<-meta(all.spec)$chlB_dry*meta(all.spec)$LMA/1000
chlB_area_norm<-lm(log(chlB_area)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$chlB_norm<-resid(chlB_area_norm)

## first to dry mass basis
meta(all.spec)$car_dry<-meta(all.spec)$car_fresh*1000/meta(all.spec)$LDMC
## then to area basis
meta(all.spec)$car_area<-meta(all.spec)$car_dry*meta(all.spec)$LMA/1000
car_area_norm<-lm(log(car_area)~log(0.1*LMA),data=meta(all.spec),na.action=na.exclude)
meta(all.spec)$car_norm<-resid(car_area_norm)

################
## eventually

saveRDS(all.spec,"ProcessedSpectra/all_spectra_and_traits.rds")
