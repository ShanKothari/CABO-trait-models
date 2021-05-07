setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels/")

library(spectrolab)
library(ggplot2)
library(reshape2)
library(stringr)

all.ref<-readRDS("ProcessedSpectra/all_ref.rds")
all.trans<-readRDS("ProcessedSpectra/all_trans.rds")
all.abs<-readRDS("ProcessedSpectra/all_abs.rds")

## check that IDs are the same and in the same order
if(meta(all.ref)$sample_id != meta(all.trans)$sample_id ||
   meta(all.ref)$sample_id != meta(all.abs)$sample_id){
  stop("sample ids not the same")
}

##############################################
## plot quantiles

# all.ref.quantiles<-quantile(all.ref,probs=c(0.025,0.25,0.5,0.75,0.975))
# all.ref.CV<-apply(all.ref,2,function(x) sd(x,na.rm=T)/mean(x,na.rm=T))
# 
# all.ref.plot<-ggplot()+
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

# meta(all.ref)$project<-factor(meta(all.ref)$project)
# all.ref.df<-as.data.frame(all.ref)
# all.ref.df$sample_name<-NULL
# 
# Girard.samples<-which(all.ref.df$project=="2018-Girard-MSc-UdeM")
# all.ref.Girard<-all.ref.df[Girard.samples,]
# Quercus.samples<-grep("Quercus",all.ref.df$species)
# all.ref.Quercus<-all.ref.df[Quercus.samples,]
# 
# all.ref.sample<-all.ref.df[sample(1:nrow(all.ref.df),200),]
# 
# all.ref.Girard.long<-melt(data=all.ref.Girard,
#                     id.vars=c("sample_id","species","project"),
#                     variable.name="wavelength",
#                     value.name="reflectance")
# all.ref.Girard.long$wavelength<-as.numeric(as.character(all.ref.Girard.long$wavelength))
# 
# all.ref.Quercus.long<-melt(data=all.ref.Quercus,
#                            id.vars=c("sample_id","species","project"),
#                            variable.name="wavelength",
#                            value.name="reflectance")
# all.ref.Quercus.long$wavelength<-as.numeric(as.character(all.ref.Quercus.long$wavelength))
# 
# all.ref.long<-melt(data=all.ref.sample,
#                            id.vars=c("sample_id","species","project"),
#                            variable.name="wavelength",
#                            value.name="reflectance")
# all.ref.long$wavelength<-as.numeric(as.character(all.ref.long$wavelength))
# 
# Quercus_plot<-ggplot()+
#   geom_line(data=all.ref.long,
#             aes(x=wavelength,y=reflectance,group=sample_id),color="gray")+
#   geom_line(data=all.ref.Quercus.long,
#             aes(x=wavelength,y=reflectance,group=sample_id),color="orange",alpha=0.2)+
#   theme_bw()+
#   theme(text = element_text(size=20),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   labs(x="Wavelength",y="Reflectance")+lims(y=c(0,1))+
#   ggtitle("Quercus spectra")
# 
# Girard_plot<-ggplot()+
#   geom_line(data=all.ref.long,
#             aes(x=wavelength,y=reflectance,group=sample_id),color="gray")+
#   geom_line(data=all.ref.Girard.long,
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

## attach to reflectance
meta(all.ref)$growth.form<-vascan$growth_form[match(meta(all.ref)$species,vascan$scientific_name)]
meta(all.ref)$growth.form[meta(all.ref)$species=="Acer pensylvanicum Linnaeus"]<-"tree"
meta(all.ref)$growth.form[meta(all.ref)$species=="Betula populifolia Marshall"]<-"tree"
meta(all.ref)$growth.form[meta(all.ref)$species=="Sorbus decora (Sargent) C.K. Schneider"]<-"tree"
meta(all.ref)$growth.form[meta(all.ref)$species=="Sorbus americana Marshall"]<-"tree"
meta(all.ref)$growth.form[meta(all.ref)$species=="Prunus pensylvanica Linnaeus f."]<-"tree"
meta(all.ref)$growth.form[meta(all.ref)$species=="Frangula alnus Miller"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Rhamnus cathartica Linnaeus"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Prunus nigra Aiton"]<-"tree"
meta(all.ref)$growth.form[meta(all.ref)$species=="Salix Linnaeus"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Agonis flexuosa (Willd.) Sweet"]<-"tree"
## the next four, I'm not too sure about the tree/shrub assignment
## making a best guess based on max height and project
meta(all.ref)$growth.form[meta(all.ref)$species=="Prunus virginiana Linnaeus"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Staphylea trifolia Linnaeus"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Oemleria cerasiformis (Torrey & A. Gray ex Hooker & Arnott) J.W. Landon"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Crataegus monogyna Jacquin"]<-"shrub"

meta(all.ref)$family<-vascan$family[match(meta(all.ref)$species,vascan$scientific_name)]
meta(all.ref)$genus<-vascan$genus[match(meta(all.ref)$species,vascan$scientific_name)]
levels(meta(all.ref)$family)<-c(levels(meta(all.ref)$family),"Myrtaceae")
levels(meta(all.ref)$genus)<-c(levels(meta(all.ref)$genus),"Agonis")
meta(all.ref)$family[meta(all.ref)$species=="Agonis flexuosa (Willd.) Sweet"]<-"Myrtaceae"
meta(all.ref)$genus[meta(all.ref)$species=="Agonis flexuosa (Willd.) Sweet"]<-"Agonis"

meta(all.ref)$functional.group<-as.character(meta(all.ref)$growth.form)
meta(all.ref)$functional.group[meta(all.ref)$family %in% c("Poaceae","Cyperaceae","Juncaceae")]<-"graminoid"
meta(all.ref)$functional.group[meta(all.ref)$family=="Fabaceae"]<-"legume"
meta(all.ref)$functional.group[meta(all.ref)$family %in% c("Pinaceae","Cupressaceae")]<-"conifer"

## attach to transmittance
meta(all.trans)$growth.form<-meta(all.ref)$growth.form
meta(all.trans)$family<-meta(all.ref)$family
meta(all.trans)$genus<-meta(all.ref)$genus
meta(all.trans)$functional.group<-meta(all.ref)$functional.group

## attach to absorptance
meta(all.abs)$growth.form<-meta(all.ref)$growth.form
meta(all.abs)$family<-meta(all.ref)$family
meta(all.abs)$genus<-meta(all.ref)$genus
meta(all.abs)$functional.group<-meta(all.ref)$functional.group

##############################################
## read leaf area/water traits

all.area<-read.csv("TraitData/LeafAreaWaterSamples/leaf_area_and_water_samples.csv")

all.area.sub<-data.frame(sample_id=all.area$sample_id,
                             SLA=all.area$specific_leaf_area_m2_kg,
                             LDMC=all.area$leaf_dry_matter_content_mg_g)

## remove bad values, based on notes in Fulcrum data
all.area.sub$LDMC[all.area.sub$sample_id %in% c("10290262","10966273","13404937",
                                                  "38530951","41809826","42944395",
                                                  "43711631","43713406","43717241",
                                                  "43718957","43718799","43712000","43721033",
                                                  "44060633","44142362","44148646",
                                                  "44683516","45108236")]<-NA

all.area.sub$SLA[all.area.sub$sample_id %in% c("13404937","44227362",
                                                 "44683516","45108236","44142362")]<-NA


## SLA in units m^2/kg
## LDMC in units mg/g
## EWT in units cm
meta(all.ref)$SLA<-all.area.sub$SLA[match(meta(all.ref)$sample_id,all.area.sub$sample_id)]
meta(all.ref)$LDMC<-all.area.sub$LDMC[match(meta(all.ref)$sample_id,all.area.sub$sample_id)]
meta(all.ref)$EWT<-with(meta(all.ref),(1/(LDMC/1000)-1)*(1/SLA*0.1))

meta(all.trans)$SLA<-all.area.sub$SLA[match(meta(all.trans)$sample_id,all.area.sub$sample_id)]
meta(all.trans)$LDMC<-all.area.sub$LDMC[match(meta(all.trans)$sample_id,all.area.sub$sample_id)]
meta(all.trans)$EWT<-with(meta(all.trans),(1/(LDMC/1000)-1)*(1/SLA*0.1))

meta(all.abs)$SLA<-all.area.sub$SLA[match(meta(all.abs)$sample_id,all.area.sub$sample_id)]
meta(all.abs)$LDMC<-all.area.sub$LDMC[match(meta(all.abs)$sample_id,all.area.sub$sample_id)]
meta(all.abs)$EWT<-with(meta(all.abs),(1/(LDMC/1000)-1)*(1/SLA*0.1))

##############################################
## read C/N

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
                           Girard.CN.sub,Hacker2018.CN.sub,
                           Hacker2019.CN.sub,PhragmitesTemporal.CN.sub,
                           Warren.CN.sub))

all.CN$Cmass[all.CN$sample_id=="38707392"]<-NA
all.CN$Nmass[all.CN$sample_id=="38707392"]<-NA

meta(all.ref)$Cmass<-all.CN$Cmass[match(meta(all.ref)$sample_id,all.CN$sample_id)]
meta(all.ref)$Nmass<-all.CN$Nmass[match(meta(all.ref)$sample_id,all.CN$sample_id)]

meta(all.trans)$Cmass<-all.CN$Cmass[match(meta(all.trans)$sample_id,all.CN$sample_id)]
meta(all.trans)$Nmass<-all.CN$Nmass[match(meta(all.trans)$sample_id,all.CN$sample_id)]

meta(all.abs)$Cmass<-all.CN$Cmass[match(meta(all.abs)$sample_id,all.CN$sample_id)]
meta(all.abs)$Nmass<-all.CN$Nmass[match(meta(all.abs)$sample_id,all.CN$sample_id)]

#############################################
## carbon fractions

CFractions<-read.csv("TraitData/CarbonFractions/carbon_fractions_bags.csv")
CFractions$species<-meta(all.ref)$species[match(CFractions$bottle_id,meta(all.ref)$sample_id)]

## removing species that have exudates that affect carbon fraction measurements
## there may also be a problem with "Polystichum munitum (Kaulfuss) C. Presl" but not confirmed
NDF_bad<-c("Ulmus rubra Muhlenberg","Acer platanoides Linnaeus",
           "Ulmus americana Linnaeus","Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")
CFractions$ndf_perc[which(CFractions$species %in% NDF_bad)]<-NA
CFractions$soluble_perc[which(CFractions$species %in% NDF_bad)]<-NA
CFractions$hemicellulose_perc[which(CFractions$species %in% NDF_bad)]<-NA

CFractions$adf_perc[which(CFractions$species=="Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")]<-NA
CFractions$adl_perc[which(CFractions$species=="Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")]<-NA
CFractions$cellulose_perc[which(CFractions$species=="Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")]<-NA
CFractions$lignin_perc[which(CFractions$species=="Alnus incana subsp. rugosa (Du Roi) R.T. Clausen")]<-NA

CFmatch.ref<-match(meta(all.ref)$sample_id,CFractions$bottle_id)
meta(all.ref)$NDFmass<-CFractions$ndf_perc[CFmatch.ref]
meta(all.ref)$ADFmass<-CFractions$adf_perc[CFmatch.ref]
meta(all.ref)$ADLmass<-CFractions$adl_perc[CFmatch.ref]
meta(all.ref)$solubles_mass<-CFractions$soluble_perc[CFmatch.ref]
meta(all.ref)$hemicellulose_mass<-CFractions$hemicellulose_perc[CFmatch.ref]
meta(all.ref)$cellulose_mass<-CFractions$cellulose_perc[CFmatch.ref]
meta(all.ref)$lignin_mass<-CFractions$lignin_perc[CFmatch.ref]

CFmatch.trans<-match(meta(all.trans)$sample_id,CFractions$bottle_id)
meta(all.trans)$NDFmass<-CFractions$ndf_perc[CFmatch.trans]
meta(all.trans)$ADFmass<-CFractions$adf_perc[CFmatch.trans]
meta(all.trans)$ADLmass<-CFractions$adl_perc[CFmatch.trans]
meta(all.trans)$solubles_mass<-CFractions$soluble_perc[CFmatch.trans]
meta(all.trans)$hemicellulose_mass<-CFractions$hemicellulose_perc[CFmatch.trans]
meta(all.trans)$cellulose_mass<-CFractions$cellulose_perc[CFmatch.trans]
meta(all.trans)$lignin_mass<-CFractions$lignin_perc[CFmatch.trans]

CFmatch.abs<-match(meta(all.abs)$sample_id,CFractions$bottle_id)
meta(all.abs)$NDFmass<-CFractions$ndf_perc[CFmatch.abs]
meta(all.abs)$ADFmass<-CFractions$adf_perc[CFmatch.abs]
meta(all.abs)$ADLmass<-CFractions$adl_perc[CFmatch.abs]
meta(all.abs)$solubles_mass<-CFractions$soluble_perc[CFmatch.abs]
meta(all.abs)$hemicellulose_mass<-CFractions$hemicellulose_perc[CFmatch.abs]
meta(all.abs)$cellulose_mass<-CFractions$cellulose_perc[CFmatch.abs]
meta(all.abs)$lignin_mass<-CFractions$lignin_perc[CFmatch.abs]

#############################################
## read pigments

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
                           Girard.pigments,Hacker2018.pigments,
                           Hacker2019.pigments,PhragmitesTemporal.pigments,
                           Warren.pigments))

all.pigments$chlA_mg_g[all.pigments$sample_id==13404937]<-NA
all.pigments$chlB_mg_g[all.pigments$sample_id==13404937]<-NA
all.pigments$carotenoides._mg_g[all.pigments$sample_id==13404937]<-NA

## remove records whose notes contain "refaire" or "refait" since these were all redone
all.pigments<-all.pigments[-which(str_detect(string = all.pigments$Notes,pattern = "efai")),]

meta(all.ref)$chlA_mass<-as.numeric(as.character(all.pigments$chlA_mg_g[match(meta(all.ref)$sample_id,all.pigments$sample_id)]))
meta(all.ref)$chlB_mass<-as.numeric(as.character(all.pigments$chlB_mg_g[match(meta(all.ref)$sample_id,all.pigments$sample_id)]))
meta(all.ref)$car_mass<-as.numeric(as.character(all.pigments$carotenoides._mg_g[match(meta(all.ref)$sample_id,all.pigments$sample_id)]))

meta(all.trans)$chlA_mass<-as.numeric(as.character(all.pigments$chlA_mg_g[match(meta(all.trans)$sample_id,all.pigments$sample_id)]))
meta(all.trans)$chlB_mass<-as.numeric(as.character(all.pigments$chlB_mg_g[match(meta(all.trans)$sample_id,all.pigments$sample_id)]))
meta(all.trans)$car_mass<-as.numeric(as.character(all.pigments$carotenoides._mg_g[match(meta(all.trans)$sample_id,all.pigments$sample_id)]))

meta(all.abs)$chlA_mass<-as.numeric(as.character(all.pigments$chlA_mg_g[match(meta(all.abs)$sample_id,all.pigments$sample_id)]))
meta(all.abs)$chlB_mass<-as.numeric(as.character(all.pigments$chlB_mg_g[match(meta(all.abs)$sample_id,all.pigments$sample_id)]))
meta(all.abs)$car_mass<-as.numeric(as.character(all.pigments$carotenoides._mg_g[match(meta(all.abs)$sample_id,all.pigments$sample_id)]))

#############################################
## ICP data

ICP12<-read.csv("TraitData/ICP/2020-03-12 1MHCl_Etienne box1,2.csv")
ICP34<-read.csv("TraitData/ICP/2020-10-20 1MHCl_Etienne box 3,4.csv")
ICP5<-read.csv("TraitData/ICP/2020-10-21 1MHCl_Etienne box 5.csv")
ICP67<-read.csv("TraitData/ICP/2020-10-22 1MHCl_Etienne box 6,7.csv")
ICP8<-read.csv("TraitData/ICP/2020-10-23 1MHCl_Etienne box 8.csv")
ICP9<-read.csv("TraitData/ICP/2020-11-11 1MHCl_Etienne box 9.csv")
ICP1011<-read.csv("TraitData/ICP/2020-12-16 1MHCl_Etienne box 10,11.csv")
ICP_boxes<-do.call(rbind,args=list(ICP12,ICP34,ICP5,ICP67,
                                   ICP8,ICP9,ICP1011))
num_cols<-c("Al","B","B.1","Ca","Cu","Fe","K","Mg","Mn","Na","P","Zn")
ICP_boxes[,num_cols]<-data.frame(sapply(ICP_boxes[,num_cols],
                                        function(x) as.numeric(as.character(x))))

ICP_Warren<-read.csv("TraitData/ICP/icp_leaf_element_concentrations_Warren.csv")
ICP_Warren<-data.frame(Sample_id=CFractions$bottle_id[match(ICP_Warren$leaf_chemistry_sample,CFractions$leaf_chemistry_sample)],
                       Scientific.name=NA,
                       Al=ICP_Warren$al_mg_g,
                       B=ICP_Warren$b_mg_g,
                       B.1=ICP_Warren$b_mg_g,
                       Ca=ICP_Warren$ca_mg_g,
                       Cu=ICP_Warren$cu_mg_g,
                       Fe=ICP_Warren$fe_mg_g,
                       K=ICP_Warren$k_mg_g,
                       Mg=ICP_Warren$mg_mg_g,
                       Mn=ICP_Warren$mn_mg_g,
                       Na=ICP_Warren$na_mg_g,
                       P=ICP_Warren$p_mg_g,
                       Zn=ICP_Warren$zn_mg_g)
ICP_all<-rbind(ICP_boxes,ICP_Warren)

ICP_all$Al[ICP_all$Sample_id=="9157296"]<-NA
ICP_all$Fe[ICP_all$Sample_id=="9157296"]<-NA
ICP_all$Al[ICP_all$Al<0]<-0
ICP_all$Na[ICP_all$Na<0]<-0

meta(all.ref)$Al_mass<-ICP_all$Al[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$B_mass<-ICP_all$B[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$B.1_mass<-ICP_all$B.1[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$Ca_mass<-ICP_all$Ca[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$Cu_mass<-ICP_all$Cu[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$Fe_mass<-ICP_all$Fe[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$K_mass<-ICP_all$K[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$Mg_mass<-ICP_all$Mg[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$Mn_mass<-ICP_all$Mn[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$Na_mass<-ICP_all$Na[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$P_mass<-ICP_all$P[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$Zn_mass<-ICP_all$Zn[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]

meta(all.trans)$Al_mass<-ICP_all$Al[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$B_mass<-ICP_all$B[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$B.1_mass<-ICP_all$B.1[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$Ca_mass<-ICP_all$Ca[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$Cu_mass<-ICP_all$Cu[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$Fe_mass<-ICP_all$Fe[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$K_mass<-ICP_all$K[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$Mg_mass<-ICP_all$Mg[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$Mn_mass<-ICP_all$Mn[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$Na_mass<-ICP_all$Na[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$P_mass<-ICP_all$P[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$Zn_mass<-ICP_all$Zn[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]

meta(all.abs)$Al_mass<-ICP_all$Al[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$B_mass<-ICP_all$B[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$B.1_mass<-ICP_all$B.1[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$Ca_mass<-ICP_all$Ca[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$Cu_mass<-ICP_all$Cu[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$Fe_mass<-ICP_all$Fe[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$K_mass<-ICP_all$K[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$Mg_mass<-ICP_all$Mg[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$Mn_mass<-ICP_all$Mn[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$Na_mass<-ICP_all$Na[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$P_mass<-ICP_all$P[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$Zn_mass<-ICP_all$Zn[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]

#############################################
## area-normalized chemical traits and
## normalization-independent chemical traits sensu Osnas

## LMA in kg/m^2
meta(all.ref)$LMA<-1/meta(all.ref)$SLA
meta(all.trans)$LMA<-1/meta(all.trans)$SLA
meta(all.abs)$LMA<-1/meta(all.abs)$SLA

## area basis, in g/cm^2
meta(all.ref)$Narea<-meta(all.ref)$Nmass*meta(all.ref)$LMA/1000
Narea_norm<-lm(log(Narea)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Nnorm<-resid(Narea_norm)

## area basis, in g/cm^2
meta(all.ref)$Carea<-meta(all.ref)$Cmass*meta(all.ref)$LMA/1000
Carea_norm<-lm(log(Carea)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Cnorm<-resid(Carea_norm)

## area basis, in g/cm^2
meta(all.ref)$solubles_area<-meta(all.ref)$solubles_mass*meta(all.ref)$LMA/1000
solubles_area_norm<-lm(log(solubles_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$solubles_norm<-resid(solubles_area_norm)

## area basis, in g/cm^2
meta(all.ref)$hemicellulose_area<-meta(all.ref)$hemicellulose_mass*meta(all.ref)$LMA/1000
hemicellulose_area_norm<-lm(log(hemicellulose_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$hemicellulose_norm<-resid(hemicellulose_area_norm)

## area basis, in g/cm^2
meta(all.ref)$cellulose_area<-meta(all.ref)$cellulose_mass*meta(all.ref)$LMA/1000
cellulose_area_norm<-lm(log(cellulose_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$cellulose_norm<-resid(cellulose_area_norm)

## area basis, in g/cm^2
meta(all.ref)$lignin_area<-meta(all.ref)$lignin_mass*meta(all.ref)$LMA/1000
lignin_area_norm<-lm(log(lignin_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$lignin_norm<-resid(lignin_area_norm)

## area basis, in g/cm^2
meta(all.ref)$chlA_area<-meta(all.ref)$chlA_mass*meta(all.ref)$LMA/1000
chlA_area_norm<-lm(log(chlA_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$chlA_norm<-resid(chlA_area_norm)

## area basis, in g/cm^2
meta(all.ref)$chlB_area<-meta(all.ref)$chlB_mass*meta(all.ref)$LMA/1000
chlB_area_norm<-lm(log(chlB_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$chlB_norm<-resid(chlB_area_norm)

## area basis, in g/cm^2
meta(all.ref)$car_area<-meta(all.ref)$car_mass*meta(all.ref)$LMA/1000
car_area_norm<-lm(log(car_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$car_norm<-resid(car_area_norm)

## area basis, in g/cm^2
meta(all.ref)$Al_area<-meta(all.ref)$Al_mass*meta(all.ref)$LMA/1000
meta(all.ref)$Al_area_omit<-meta(all.ref)$Al_area
meta(all.ref)$Al_area_omit[which(meta(all.ref)$Al_area_omit==0)]<-NA
Al_area_norm<-lm(log(Al_area_omit)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Al_norm<-resid(Al_area_norm)
meta(all.ref)$Al_area_omit<-NULL

## area basis, in g/cm^2
meta(all.ref)$Ca_area<-meta(all.ref)$Ca_mass*meta(all.ref)$LMA/1000
Ca_area_norm<-lm(log(Ca_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Ca_norm<-resid(Ca_area_norm)

## area basis, in g/cm^2
meta(all.ref)$Cu_area<-meta(all.ref)$Cu_mass*meta(all.ref)$LMA/1000
meta(all.ref)$Cu_area_omit<-meta(all.ref)$Cu_area
meta(all.ref)$Cu_area_omit[which(meta(all.ref)$Cu_area_omit==0)]<-NA
Cu_area_norm<-lm(log(Cu_area_omit)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Cu_norm<-resid(Cu_area_norm)
meta(all.ref)$Cu_area_omit<-NULL

## area basis, in g/cm^2
meta(all.ref)$Fe_area<-meta(all.ref)$Fe_mass*meta(all.ref)$LMA/1000
Fe_area_norm<-lm(log(Fe_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Fe_norm<-resid(Fe_area_norm)

## area basis, in g/cm^2
meta(all.ref)$K_area<-meta(all.ref)$K_mass*meta(all.ref)$LMA/1000
K_area_norm<-lm(log(K_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$K_norm<-resid(K_area_norm)

## area basis, in g/cm^2
meta(all.ref)$Mg_area<-meta(all.ref)$Mg_mass*meta(all.ref)$LMA/1000
Mg_area_norm<-lm(log(Mg_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Mg_norm<-resid(Mg_area_norm)

## area basis, in g/cm^2
meta(all.ref)$Mn_area<-meta(all.ref)$Mn_mass*meta(all.ref)$LMA/1000
meta(all.ref)$Mn_area_omit<-meta(all.ref)$Mn_area
meta(all.ref)$Mn_area_omit[which(meta(all.ref)$Mn_area_omit==0)]<-NA
Mn_area_norm<-lm(log(Mn_area_omit)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Mn_norm<-resid(Mn_area_norm)
meta(all.ref)$Mn_area_omit<-NULL

## area basis, in g/cm^2
meta(all.ref)$Na_area<-meta(all.ref)$Na_mass*meta(all.ref)$LMA/1000
meta(all.ref)$Na_area_omit<-meta(all.ref)$Na_area
meta(all.ref)$Na_area_omit[which(meta(all.ref)$Na_area_omit==0)]<-NA
Na_area_norm<-lm(log(Na_area_omit)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Na_norm<-resid(Na_area_norm)
meta(all.ref)$Na_area_omit<-NULL

## area basis, in g/cm^2
meta(all.ref)$P_area<-meta(all.ref)$P_mass*meta(all.ref)$LMA/1000
P_area_norm<-lm(log(P_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$P_norm<-resid(P_area_norm)

## area basis, in g/cm^2
meta(all.ref)$Zn_area<-meta(all.ref)$Zn_mass*meta(all.ref)$LMA/1000
Zn_area_norm<-lm(log(Zn_area)~log(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
meta(all.ref)$Zn_norm<-resid(Zn_area_norm)

################
## eventually

saveRDS(all.ref,"ProcessedSpectra/all_ref_and_traits.rds")
saveRDS(all.trans,"ProcessedSpectra/all_trans_and_traits.rds")
saveRDS(all.abs,"ProcessedSpectra/all_abs_and_traits.rds")
