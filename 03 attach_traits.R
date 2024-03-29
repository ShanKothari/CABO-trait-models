setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels/")

library(spectrolab)
library(ggplot2)
library(reshape2)
library(stringr)
library(patchwork)
library(ggfortify)

all.ref<-readRDS("ProcessedSpectra/all_ref.rds")
all.trans<-readRDS("ProcessedSpectra/all_trans.rds")
all.abs<-readRDS("ProcessedSpectra/all_abs.rds")

## check that IDs are the same and in the same order
if(any(meta(all.ref)$sample_id != meta(all.trans)$sample_id) ||
   any(meta(all.ref)$sample_id != meta(all.abs)$sample_id)){
  stop("sample ids not the same")
}

#############################################
## growth form and taxonomic information

## read data from VASCAN (Desmet & Brouillet 2013)
vascan<-read.csv("SummaryData/vascan.csv")

## attach to spectral metadata
meta(all.ref)$growth.form<-vascan$growth_form[match(meta(all.ref)$species,vascan$scientific_name)]
## manual disambiguation when VASCAN includes multiple growth forms
## in most cases, based on specific individuals sampled by CABO
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
meta(all.ref)$growth.form[meta(all.ref)$species=="Prunus virginiana Linnaeus"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Staphylea trifolia Linnaeus"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Oemleria cerasiformis (Torrey & A. Gray ex Hooker & Arnott) J.W. Landon"]<-"shrub"
meta(all.ref)$growth.form[meta(all.ref)$species=="Crataegus monogyna Jacquin"]<-"shrub"

## add taxonomic information
meta(all.ref)$order<-vascan$order[match(meta(all.ref)$species,vascan$scientific_name)]
meta(all.ref)$family<-vascan$family[match(meta(all.ref)$species,vascan$scientific_name)]
meta(all.ref)$genus<-vascan$genus[match(meta(all.ref)$species,vascan$scientific_name)]
levels(meta(all.ref)$family)<-c(levels(meta(all.ref)$family),"Myrtaceae")
levels(meta(all.ref)$genus)<-c(levels(meta(all.ref)$genus),"Agonis")
meta(all.ref)$order[meta(all.ref)$species=="Agonis flexuosa (Willd.) Sweet"]<-"Myrtales"
meta(all.ref)$family[meta(all.ref)$species=="Agonis flexuosa (Willd.) Sweet"]<-"Myrtaceae"
meta(all.ref)$genus[meta(all.ref)$species=="Agonis flexuosa (Willd.) Sweet"]<-"Agonis"

## disambiguate functional groups by families
meta(all.ref)$functional.group<-as.character(meta(all.ref)$growth.form)
meta(all.ref)$functional.group[meta(all.ref)$family %in% c("Poaceae","Cyperaceae","Juncaceae","Typhaceae")]<-"graminoid"
meta(all.ref)$functional.group[meta(all.ref)$family %in% c("Pinaceae","Cupressaceae")]<-"conifer"
meta(all.ref)$functional.group[meta(all.ref)$order %in% c("Polypodiales","Osmundales")]<-"fern"
meta(all.ref)$functional.group[meta(all.ref)$functional.group=="herb"]<-"forb"
meta(all.ref)$functional.group[meta(all.ref)$functional.group=="tree"]<-"broadleaf"

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
                         LMA=all.area$leaf_mass_per_area_g_m2/1000,
                         LDMC=all.area$leaf_dry_matter_content_mg_g,
                         EWT=all.area$equivalent_water_thickness_cm*10)

## remove bad values, based on notes in Fulcrum data
bad.LDMC<-c("10290262","10966273","13404937","38530951","41809826","42944395",
            "43711631","43713406","43717241","43718957","43718799","43712000",
            "43721033","44060633","44142362","44148646","44683516","45108236")
bad.LMA<-c("13404937","44227362","44683516","45108236","44142362")

all.area.sub$LDMC[all.area.sub$sample_id %in% bad.LDMC]<-NA
all.area.sub$LMA[all.area.sub$sample_id %in% bad.LMA]<-NA
all.area.sub$EWT[all.area.sub$sample_id %in% c(bad.LMA,bad.LDMC)]<-NA

## LMA in kg/m^2
## LDMC in units mg/g
## EWT in mm
meta(all.ref)$LMA<-all.area.sub$LMA[match(meta(all.ref)$sample_id,all.area.sub$sample_id)]
meta(all.ref)$LDMC<-all.area.sub$LDMC[match(meta(all.ref)$sample_id,all.area.sub$sample_id)]
meta(all.ref)$EWT<-all.area.sub$EWT[match(meta(all.ref)$sample_id,all.area.sub$sample_id)]

meta(all.trans)$LMA<-all.area.sub$LMA[match(meta(all.trans)$sample_id,all.area.sub$sample_id)]
meta(all.trans)$LDMC<-all.area.sub$LDMC[match(meta(all.trans)$sample_id,all.area.sub$sample_id)]
meta(all.trans)$EWT<-all.area.sub$EWT[match(meta(all.trans)$sample_id,all.area.sub$sample_id)]

meta(all.abs)$LMA<-all.area.sub$LMA[match(meta(all.abs)$sample_id,all.area.sub$sample_id)]
meta(all.abs)$LDMC<-all.area.sub$LDMC[match(meta(all.abs)$sample_id,all.area.sub$sample_id)]
meta(all.abs)$EWT<-all.area.sub$EWT[match(meta(all.abs)$sample_id,all.area.sub$sample_id)]

##############################################
## read C/N

all.CN<-read.csv("TraitData/CNSamples/c_n_leaf_concentrations.csv")

## likely error
all.CN$c_perc[all.CN$sample_id=="38707392"]<-NA
all.CN$n_perc[all.CN$sample_id=="38707392"]<-NA

meta(all.ref)$Cmass<-all.CN$c_perc[match(meta(all.ref)$sample_id,all.CN$sample_id)]
meta(all.ref)$Nmass<-all.CN$n_perc[match(meta(all.ref)$sample_id,all.CN$sample_id)]

meta(all.trans)$Cmass<-all.CN$c_perc[match(meta(all.trans)$sample_id,all.CN$sample_id)]
meta(all.trans)$Nmass<-all.CN$n_perc[match(meta(all.trans)$sample_id,all.CN$sample_id)]

meta(all.abs)$Cmass<-all.CN$c_perc[match(meta(all.abs)$sample_id,all.CN$sample_id)]
meta(all.abs)$Nmass<-all.CN$n_perc[match(meta(all.abs)$sample_id,all.CN$sample_id)]

#############################################
## carbon fractions

CFractions<-read.csv("TraitData/CarbonFractions/carbon_fractions_bags.csv")
CFractions$species<-meta(all.ref)$species[match(CFractions$bottle_id,meta(all.ref)$sample_id)]
## remove samples that were marked as to be redone
CFractions<-CFractions[-grep("redone",CFractions$sample_remarks),]

## removing species that have exudates that affect carbon fraction measurements
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
meta(all.ref)$solubles_mass<-CFractions$soluble_perc[CFmatch.ref]
meta(all.ref)$hemicellulose_mass<-CFractions$hemicellulose_perc[CFmatch.ref]
meta(all.ref)$cellulose_mass<-CFractions$cellulose_perc[CFmatch.ref]
meta(all.ref)$lignin_mass<-CFractions$lignin_perc[CFmatch.ref]

CFmatch.trans<-match(meta(all.trans)$sample_id,CFractions$bottle_id)
meta(all.trans)$solubles_mass<-CFractions$soluble_perc[CFmatch.trans]
meta(all.trans)$hemicellulose_mass<-CFractions$hemicellulose_perc[CFmatch.trans]
meta(all.trans)$cellulose_mass<-CFractions$cellulose_perc[CFmatch.trans]
meta(all.trans)$lignin_mass<-CFractions$lignin_perc[CFmatch.trans]

CFmatch.abs<-match(meta(all.abs)$sample_id,CFractions$bottle_id)
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

## likely error
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

## remove extreme outliers
ICP_all$Al[ICP_all$Sample_id=="9157296"]<-NA
ICP_all$Fe[ICP_all$Sample_id=="9157296"]<-NA
ICP_all$Cu[ICP_all$Sample_id=="12176890"]<-NA
ICP_all$Cu[ICP_all$Sample_id=="11914931"]<-NA
ICP_all$Cu[ICP_all$Sample_id=="11995935"]<-NA
ICP_all$Zn[ICP_all$Sample_id=="13384790"]<-NA
## set samples with negative values (below detection limit) to 0
ICP_all$Al[ICP_all$Al<0]<-0
ICP_all$Na[ICP_all$Na<0]<-0

## we don't build models with the boron estimates
## but they're here; details on the two methods are in
## Turner et al. 2016 https://doi.org/10.1080/00103624.2016.1228952
meta(all.ref)$Al_mass<-ICP_all$Al[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$B208.9_mass<-ICP_all$B[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
meta(all.ref)$B249.8_mass<-ICP_all$B.1[match(meta(all.ref)$sample_id,ICP_all$Sample_id)]
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
meta(all.trans)$B208.9_mass<-ICP_all$B[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
meta(all.trans)$B249.8_mass<-ICP_all$B.1[match(meta(all.trans)$sample_id,ICP_all$Sample_id)]
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
meta(all.abs)$B208.9_mass<-ICP_all$B[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
meta(all.abs)$B249.8_mass<-ICP_all$B.1[match(meta(all.abs)$sample_id,ICP_all$Sample_id)]
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

## area basis, in g/cm^2
meta(all.ref)$Narea<-meta(all.ref)$Nmass*meta(all.ref)$LMA/1000
meta(all.ref)$Carea<-meta(all.ref)$Cmass*meta(all.ref)$LMA/1000
meta(all.ref)$solubles_area<-meta(all.ref)$solubles_mass*meta(all.ref)$LMA/1000
meta(all.ref)$hemicellulose_area<-meta(all.ref)$hemicellulose_mass*meta(all.ref)$LMA/1000
meta(all.ref)$cellulose_area<-meta(all.ref)$cellulose_mass*meta(all.ref)$LMA/1000
meta(all.ref)$lignin_area<-meta(all.ref)$lignin_mass*meta(all.ref)$LMA/1000
meta(all.ref)$chlA_area<-meta(all.ref)$chlA_mass*meta(all.ref)$LMA/10000
meta(all.ref)$chlB_area<-meta(all.ref)$chlB_mass*meta(all.ref)$LMA/10000
meta(all.ref)$car_area<-meta(all.ref)$car_mass*meta(all.ref)$LMA/10000
meta(all.ref)$Al_area<-meta(all.ref)$Al_mass*meta(all.ref)$LMA/10000
meta(all.ref)$Ca_area<-meta(all.ref)$Ca_mass*meta(all.ref)$LMA/10000
meta(all.ref)$Cu_area<-meta(all.ref)$Cu_mass*meta(all.ref)$LMA/10000
meta(all.ref)$Fe_area<-meta(all.ref)$Fe_mass*meta(all.ref)$LMA/10000
meta(all.ref)$K_area<-meta(all.ref)$K_mass*meta(all.ref)$LMA/10000
meta(all.ref)$Mg_area<-meta(all.ref)$Mg_mass*meta(all.ref)$LMA/10000
meta(all.ref)$Mn_area<-meta(all.ref)$Mn_mass*meta(all.ref)$LMA/10000
meta(all.ref)$Na_area<-meta(all.ref)$Na_mass*meta(all.ref)$LMA/10000
meta(all.ref)$P_area<-meta(all.ref)$P_mass*meta(all.ref)$LMA/10000
meta(all.ref)$Zn_area<-meta(all.ref)$Zn_mass*meta(all.ref)$LMA/10000

## example line of code to calculate mass-proportionality
Narea_norm<-lm(log10(Narea)~log10(0.1*LMA),data=meta(all.ref),na.action=na.exclude)
summary(Narea_norm)

################
## save data

saveRDS(all.ref,"ProcessedSpectra/all_ref_and_traits.rds")
saveRDS(all.trans,"ProcessedSpectra/all_trans_and_traits.rds")
saveRDS(all.abs,"ProcessedSpectra/all_abs_and_traits.rds")

################
## plotting quantiles of spectra

all.ref<-all.ref[which(meta(all.ref)$project!="2019-Pardo-MSc-UdeM")]
all.trans<-all.trans[which(meta(all.trans)$project!="2019-Pardo-MSc-UdeM")]
all.abs<-all.abs[which(meta(all.abs)$project!="2019-Pardo-MSc-UdeM")]

ref_quantiles<-quantile(all.ref,probs=c(0.025,0.25,0.5,0.75,0.975))
trans_quantiles<-quantile(all.trans,probs=c(0.025,0.25,0.5,0.75,0.975))
abs_quantiles<-quantile(all.abs,probs=c(0.025,0.25,0.5,0.75,0.975))

ref_CV<-apply(all.ref,2,function(x) sd(x)/mean(x))
trans_CV<-apply(all.trans,2,function(x) sd(x)/mean(x))
abs_CV<-apply(all.abs,2,function(x) sd(x)/mean(x))
CV_df<-data.frame(wavelength=400:2400,ref_CV,trans_CV,abs_CV)

ref_spec_plot<-ggplot()+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(ref_quantiles)[1,],
                  ymax = as.matrix(ref_quantiles)[5,]),
              alpha = 0.5,fill = "light blue")+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(ref_quantiles)[2,],
                  ymax = as.matrix(ref_quantiles)[4,]),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,y=as.matrix(ref_quantiles)[3,]),size=1,color="black")+
  geom_line(aes(x=400:2400,y=ref_CV),size=1,color="red")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0,0.3,0.1,0),"in"))+
  labs(x="Wavelength (nm)",y="Reflectance (or CV)")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "A")

trans_spec_plot<-ggplot()+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(trans_quantiles)[1,],
                  ymax = as.matrix(trans_quantiles)[5,]),
              alpha = 0.5,fill = "light blue")+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(trans_quantiles)[2,],
                  ymax = as.matrix(trans_quantiles)[4,]),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,y=as.matrix(trans_quantiles)[3,]),size=1,color="black")+
  geom_line(aes(x=400:2400,y=trans_CV),size=1,color="red")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0,0.3,0.1,0),"in"))+
  labs(x="Wavelength (nm)",y="Transmittance (or CV)")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "B")

abs_spec_plot<-ggplot()+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(abs_quantiles)[1,],
                  ymax = as.matrix(abs_quantiles)[5,]),
              alpha = 0.5,fill = "light blue")+
  geom_ribbon(aes(x=400:2400,
                  ymin = as.matrix(abs_quantiles)[2,],
                  ymax = as.matrix(abs_quantiles)[4,]),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,y=as.matrix(abs_quantiles)[3,]),size=1,color="black")+
  geom_line(aes(x=400:2400,y=abs_CV),size=1,color="red")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0.3,0,0),"in"))+
  labs(x="Wavelength (nm)",y="Absorptance (or CV)")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "C")

pdf("Images/all_spec.pdf",width=7,height=11)
ref_spec_plot+trans_spec_plot+abs_spec_plot+plot_layout(ncol = 1)
dev.off()

################################
## plot spectra for each functional group

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

all.ref.df<-data.frame(sample_id=meta(all.ref)$sample_id,
                       functional.group=meta(all.ref)$functional.group,
                       as.matrix(all.ref))
## remove ferns and vines, which are represented by just one species each
all.ref.df<-all.ref.df[-which(all.ref.df$functional.group %in% c("fern","vine")),]
all.ref.long<-melt(all.ref.df,id.vars = c("sample_id","functional.group"))
all.ref.long$variable<-as.numeric(gsub("X","",all.ref.long$variable))

ref_fg_spec_plot<-ggplot(all.ref.long,aes(x=variable,y=value))+
  stat_summary(fun=median,na.rm=T,geom="line",size=1,aes(color=functional.group))+
  geom_line(data=CV_df,
            aes(x=wavelength,y=ref_CV),size=1.5,linetype="longdash")+
  theme_bw()+theme(text=element_text(size=15))+
  labs(x="Wavelength (nm)",y="Reflectance (or CV)")+
  guides(color=guide_legend("Functional group",nrow = 2))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "A")+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])

ref_fg_spec_plot_tagless<-ggplot(all.ref.long,aes(x=variable,y=value))+
  stat_summary(fun=median,na.rm=T,geom="line",size=1,aes(color=functional.group))+
  geom_line(data=CV_df,
            aes(x=wavelength,y=ref_CV),size=1.5,linetype="longdash")+
  theme_bw()+theme(text=element_text(size=15))+
  labs(x="Wavelength (nm)",y="Reflectance (or CV)")+
  guides(color=guide_legend("Functional group",nrow = 2))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])

all.trans.df<-data.frame(sample_id=meta(all.trans)$sample_id,
                         functional.group=meta(all.trans)$functional.group,
                         as.matrix(all.trans))
all.trans.df<-all.trans.df[-which(all.trans.df$functional.group %in% c("fern","vine")),]
all.trans.long<-melt(all.trans.df,id.vars = c("sample_id","functional.group"))
all.trans.long$variable<-as.numeric(gsub("X","",all.trans.long$variable))

trans_fg_spec_plot<-ggplot(all.trans.long,aes(x=variable,y=value))+
  stat_summary(fun=median,na.rm=T,geom="line",size=1,aes(color=functional.group))+
  geom_line(data=CV_df,
            aes(x=wavelength,y=trans_CV),size=1.5,linetype="longdash")+
  theme_bw()+theme(text=element_text(size=15))+
  labs(x="Wavelength (nm)",y="Transmittance (or CV)")+
  guides(color=guide_legend("Functional group",nrow = 2))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "B")+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])

all.abs.df<-data.frame(sample_id=meta(all.abs)$sample_id,
                       functional.group=meta(all.abs)$functional.group,
                       as.matrix(all.abs))
all.abs.df<-all.abs.df[-which(all.abs.df$functional.group %in% c("fern","vine")),]
all.abs.long<-melt(all.abs.df,id.vars = c("sample_id","functional.group"))
all.abs.long$variable<-as.numeric(gsub("X","",all.abs.long$variable))

abs_fg_spec_plot<-ggplot(all.abs.long,aes(x=variable,y=value))+
  stat_summary(fun=median,na.rm=T,geom="line",size=1,aes(color=functional.group))+
  geom_line(data=CV_df,
            aes(x=wavelength,y=abs_CV),size=1.5,linetype="longdash")+
  theme_bw()+theme(text=element_text(size=15))+
  labs(x="Wavelength (nm)",y="Absorptance (or CV)")+
  guides(color=guide_legend("Functional group",nrow = 2))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  labs(tag = "C")+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])

pdf("Images/all_spec_fg.pdf",width=7,height=12)
ref_fg_spec_plot+trans_fg_spec_plot+
  abs_fg_spec_plot+plot_layout(ncol = 1,guides="collect") &
  theme(legend.position = "bottom")
dev.off()

pdf("Images/ref_spec_fg.pdf",width=7,height=5)
ref_fg_spec_plot_tagless &
  theme(legend.position = "bottom")
dev.off()