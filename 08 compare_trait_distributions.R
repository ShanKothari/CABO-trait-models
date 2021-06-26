setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(ggplot2)
library(FNN)
library(lme4)
library(ggpubr)

#################################
## import TRY data

## TRY trait number reference
# SLA no petiole: 3115*
# LDMC: 47*
# Leaf N: 14*
# Leaf C: 13*
# Solubles: 85*
# Hemicellulose: 94*
# Cellulose: 92*
# Lignin: 87*
# Chlorophyll: 164*
# Chlorophyll a: 3474*
# Chlorophyll b: 3475*
# Carotenoids: 809*
# Al: 249*
# B: 250*
# Ca: 252*
# Cu: 255*
# Fe: 256*
# K: 44*
# Mg: 257*
# Mn: 258*
# Na: 260*
# P: 15*
# Zn: 268*
# Phenols: 147
# Tannins: 148

# ## leaf structure -- includes LDMC and SLA (TRY request 14553)
# TRY_struc<-read.table("C:/Users/kotha020/Downloads/TRYStructure/14553.txt",
#                       sep = "\t",fill=T,header=T,quote="")
# 
# TRY_struc_sub<-TRY_struc[which(TRY_struc$TraitID %in% c(3115,47)),]
# TRY_struc_sub$StdValue<-as.numeric(as.character(TRY_struc_sub$StdValue))
# TRY_struc_sub$ErrorRisk<-as.numeric(as.character(TRY_struc_sub$ErrorRisk))
# TRY_struc_sub<-TRY_struc_sub[!is.na(TRY_struc_sub$StdValue),]
# TRY_struc_sub<-TRY_struc_sub[which(TRY_struc_sub$ErrorRisk<3),]
# TRY_struc_sub<-TRY_struc_sub[,c("DatasetID","SpeciesName","AccSpeciesID","ObservationID",
#                             "ObsDataID","TraitID","TraitName","StdValue","UnitName","ErrorRisk")]
# 
# TRY_LDMC<-TRY_struc_sub[which(TRY_struc_sub$TraitID==47),]
# TRY_SLA<-TRY_struc_sub[which(TRY_struc_sub$TraitID==3115),]
# 
# write.csv(TRY_LDMC,"C:/Users/kotha020/Downloads/TRY_LDMC.csv")
# write.csv(TRY_SLA,"C:/Users/kotha020/Downloads/TRY_SLA.csv")
# 
# ## leaf macronutrients -- includes C, N, and P (TRY request 14554)
# TRY_CNP<-read.table("C:/Users/kotha020/Downloads/TRYNutrients/14554.txt",
#                       sep = "\t",fill=T,header=T,quote = "")
# 
# TRY_CNP_sub<-TRY_CNP[which(TRY_CNP$TraitID %in% 13:15),]
# TRY_CNP_sub$StdValue<-as.numeric(as.character(TRY_CNP_sub$StdValue))
# TRY_CNP_sub$ErrorRisk<-as.numeric(as.character(TRY_CNP_sub$ErrorRisk))
# TRY_CNP_sub<-TRY_CNP_sub[!is.na(TRY_CNP_sub$StdValue),]
# TRY_CNP_sub<-TRY_CNP_sub[which(TRY_CNP_sub$ErrorRisk<3),]
# TRY_CNP_sub<-TRY_CNP_sub[,c("DatasetID","SpeciesName","AccSpeciesID","ObservationID",
#                             "ObsDataID","TraitID","TraitName","StdValue","UnitName","ErrorRisk")]
# 
# TRY_C<-TRY_CNP_sub[which(TRY_CNP_sub$TraitID==13),]
# TRY_N<-TRY_CNP_sub[which(TRY_CNP_sub$TraitID==14),]
# TRY_P<-TRY_CNP_sub[which(TRY_CNP_sub$TraitID==15),]
# 
# write.csv(TRY_C,"C:/Users/kotha020/Downloads/TRY_C.csv")
# write.csv(TRY_N,"C:/Users/kotha020/Downloads/TRY_N.csv")
# write.csv(TRY_P,"C:/Users/kotha020/Downloads/TRY_P.csv")
# 
# ## leaf carbon fractions (TRY request 14555)
# TRY_CFrac<-read.table("C:/Users/kotha020/Downloads/TRYCFrac/14555.txt",
#                     sep = "\t",fill=T,header=T,quote = "")
# 
# TRY_CFrac_sub<-TRY_CFrac[which(TRY_CFrac$TraitID %in% c(85,87,92,94)),]
# TRY_CFrac_sub$StdValue<-as.numeric(as.character(TRY_CFrac_sub$StdValue))
# TRY_CFrac_sub$ErrorRisk<-as.numeric(as.character(TRY_CFrac_sub$ErrorRisk))
# TRY_CFrac_sub<-TRY_CFrac_sub[!is.na(TRY_CFrac_sub$StdValue),]
# TRY_CFrac_sub<-TRY_CFrac_sub[which(TRY_CFrac_sub$ErrorRisk<3),]
# TRY_CFrac_sub<-TRY_CFrac_sub[,c("DatasetID","SpeciesName","AccSpeciesID","ObservationID",
#                             "ObsDataID","TraitID","TraitName","StdValue","UnitName","ErrorRisk")]
# 
# TRY_solubles<-TRY_CFrac_sub[which(TRY_CFrac_sub$TraitID==85),]
# TRY_lignin<-TRY_CFrac_sub[which(TRY_CFrac_sub$TraitID==87),]
# TRY_cellulose<-TRY_CFrac_sub[which(TRY_CFrac_sub$TraitID==92),]
# TRY_hemicellulose<-TRY_CFrac_sub[which(TRY_CFrac_sub$TraitID==94),]
# 
# write.csv(TRY_solubles,"C:/Users/kotha020/Downloads/TRY_solubles.csv")
# write.csv(TRY_lignin,"C:/Users/kotha020/Downloads/TRY_lignin.csv")
# write.csv(TRY_cellulose,"C:/Users/kotha020/Downloads/TRY_cellulose.csv")
# write.csv(TRY_hemicellulose,"C:/Users/kotha020/Downloads/TRY_hemicellulose.csv")
# 
## leaf pigments (TRY request 14556)
# TRY_pigments<-read.table("C:/Users/kotha020/Downloads/TRYPigments/14556.txt",
#                       sep = "\t",fill=T,header=T,quote = "")
# 
# TRY_pigments_sub<-TRY_pigments[which(TRY_pigments$TraitID %in% c(164,809,3474,3475)),]
# TRY_pigments_sub$OrigValueStr<-as.numeric(as.character(TRY_pigments_sub$OrigValueStr))
# TRY_pigments_sub$StdValue<-as.numeric(as.character(TRY_pigments_sub$StdValue))
# TRY_pigments_sub$ErrorRisk<-as.numeric(as.character(TRY_pigments_sub$ErrorRisk))
## all records for ChlA, ChlB, and car have only original values
## and are missing error risk assessments,]

# TRY_Chl<-TRY_pigments_sub[which(TRY_pigments_sub$TraitID==164),]
# TRY_Chl<-TRY_Chl[which(TRY_Chl$ErrorRisk<3),]
# TRY_Chl<-TRY_Chl[,c("DatasetID","SpeciesName","AccSpeciesID","ObservationID",
#                                       "ObsDataID","TraitID","TraitName","StdValue","UnitName","ErrorRisk")]
# 
# micro_ChlA<-which(TRY_ChlA$OrigUnitStr %in% c("micro g / g","micro g/g"))
# TRY_ChlA<-TRY_pigments_sub[which(TRY_pigments_sub$TraitID==3474),]
# TRY_ChlA$StdValue<-TRY_ChlA$OrigValueStr
# TRY_ChlA$StdValue[micro_chlA]<-TRY_ChlA$OrigValueStr[micro_chlA]/1000
# TRY_ChlA$UnitName<-"mg g-1"
# TRY_ChlA<-TRY_ChlA[,c("DatasetID","SpeciesName","AccSpeciesID","ObservationID",
#                     "ObsDataID","TraitID","TraitName","StdValue","UnitName","ErrorRisk")]
# 
# 
# micro_ChlB<-which(TRY_ChlB$OrigUnitStr %in% c("micro g / g","micro g/g"))
# TRY_ChlB<-TRY_pigments_sub[which(TRY_pigments_sub$TraitID==3475),]
# TRY_ChlB$StdValue<-TRY_ChlB$OrigValueStr
# TRY_ChlB$StdValue[micro_ChlB]<-TRY_ChlB$OrigValueStr[micro_ChlB]/1000
# TRY_ChlB$UnitName<-"mg g-1"
# TRY_ChlB<-TRY_ChlB[,c("DatasetID","SpeciesName","AccSpeciesID","ObservationID",
#                     "ObsDataID","TraitID","TraitName","StdValue","UnitName","ErrorRisk")]
# 
# micro_Car<-which(TRY_Car$OrigUnitStr %in% c("micro g / g","micro g/g"))
# TRY_Car<-TRY_pigments_sub[which(TRY_pigments_sub$TraitID==809),]
# TRY_Car$StdValue<-TRY_Car$OrigValueStr
# TRY_Car$StdValue[micro_Car]<-TRY_Car$OrigValueStr[micro_Car]/1000
# TRY_Car$UnitName<-"mg g-1"
# TRY_Car<-TRY_Car[,c("DatasetID","SpeciesName","AccSpeciesID","ObservationID",
#                     "ObsDataID","TraitID","TraitName","StdValue","UnitName","ErrorRisk")]
# 
# write.csv(TRY_Chl,"C:/Users/kotha020/Downloads/TRY_Chl.csv")
# write.csv(TRY_ChlA,"C:/Users/kotha020/Downloads/TRY_ChlA.csv")
# write.csv(TRY_ChlB,"C:/Users/kotha020/Downloads/TRY_ChlB.csv")
# write.csv(TRY_Car,"C:/Users/kotha020/Downloads/TRY_Car.csv")
# 
# ## leaf micronutrients (TRY request 14557)
# TRY_Micro<-read.table("C:/Users/kotha020/Downloads/TRYMicro/14557.txt",
#                       sep = "\t",fill=T,header=T,quote = "")
# 
# TRY_Micro_sub<-TRY_Micro[which(TRY_Micro$TraitID %in% c(249,250,252,255,256,44,
#                                                         257,258,260,268)),]
# TRY_Micro_sub$StdValue<-as.numeric(as.character(TRY_Micro_sub$StdValue))
# TRY_Micro_sub$ErrorRisk<-as.numeric(as.character(TRY_Micro_sub$ErrorRisk))
# TRY_Micro_sub<-TRY_Micro_sub[!is.na(TRY_Micro_sub$StdValue),]
# TRY_Micro_sub<-TRY_Micro_sub[which(TRY_Micro_sub$ErrorRisk<3),]
# TRY_Micro_sub<-TRY_Micro_sub[,c("DatasetID","SpeciesName","AccSpeciesID","ObservationID",
#                                 "ObsDataID","TraitID","TraitName","StdValue","UnitName","ErrorRisk")]
# 
# TRY_Al<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==249),]
# TRY_B<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==250),]
# TRY_Ca<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==252),]
# TRY_Cu<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==255),]
# TRY_Fe<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==256),]
# TRY_K<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==44),]
# TRY_Mg<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==257),]
# TRY_Mn<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==258),]
# TRY_Na<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==260),]
# TRY_Zn<-TRY_Micro_sub[which(TRY_Micro_sub$TraitID==268),]
# 
# write.csv(TRY_Al,"C:/Users/kotha020/Downloads/TRY_Al.csv")
# write.csv(TRY_B,"C:/Users/kotha020/Downloads/TRY_B.csv")
# write.csv(TRY_Ca,"C:/Users/kotha020/Downloads/TRY_Ca.csv")
# write.csv(TRY_Cu,"C:/Users/kotha020/Downloads/TRY_Cu.csv")
# write.csv(TRY_Fe,"C:/Users/kotha020/Downloads/TRY_Fe.csv")
# write.csv(TRY_K,"C:/Users/kotha020/Downloads/TRY_K.csv")
# write.csv(TRY_Mn,"C:/Users/kotha020/Downloads/TRY_Mn.csv")
# write.csv(TRY_Mg,"C:/Users/kotha020/Downloads/TRY_Mg.csv")
# write.csv(TRY_Na,"C:/Users/kotha020/Downloads/TRY_Na.csv")
# write.csv(TRY_Zn,"C:/Users/kotha020/Downloads/TRY_Zn.csv")

#####################################
## read CABO trait data

ref.traits<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")
## the Pardo dataset has no trait data (yet)
ref.traits<-ref.traits[-which(meta(ref.traits)$project=="2019-Pardo-MSc-UdeM")]

trait.df<-meta(ref.traits)

TRY_SLA<-read.csv("C:/Users/kotha020/Downloads/TRY_SLA.csv")
TRY_LDMC<-read.csv("C:/Users/kotha020/Downloads/TRY_LDMC.csv")
TRY_N<-read.csv("C:/Users/kotha020/Downloads/TRY_N.csv")
TRY_C<-read.csv("C:/Users/kotha020/Downloads/TRY_C.csv")
TRY_solubles<-read.csv("C:/Users/kotha020/Downloads/TRY_solubles.csv")
TRY_hemicellulose<-read.csv("C:/Users/kotha020/Downloads/TRY_hemicellulose.csv")
TRY_cellulose<-read.csv("C:/Users/kotha020/Downloads/TRY_cellulose.csv")
TRY_lignin<-read.csv("C:/Users/kotha020/Downloads/TRY_lignin.csv")
TRY_Chl<-read.csv("C:/Users/kotha020/Downloads/TRY_Chl.csv")
TRY_ChlA<-read.csv("C:/Users/kotha020/Downloads/TRY_ChlA.csv")
TRY_ChlB<-read.csv("C:/Users/kotha020/Downloads/TRY_ChlB.csv")
TRY_Car<-read.csv("C:/Users/kotha020/Downloads/TRY_Car.csv")
TRY_Al<-read.csv("C:/Users/kotha020/Downloads/TRY_Al.csv")
TRY_Ca<-read.csv("C:/Users/kotha020/Downloads/TRY_Ca.csv")
TRY_Cu<-read.csv("C:/Users/kotha020/Downloads/TRY_Cu.csv")
TRY_Fe<-read.csv("C:/Users/kotha020/Downloads/TRY_Fe.csv")
TRY_K<-read.csv("C:/Users/kotha020/Downloads/TRY_K.csv")
TRY_Mg<-read.csv("C:/Users/kotha020/Downloads/TRY_Mg.csv")
TRY_Mn<-read.csv("C:/Users/kotha020/Downloads/TRY_Mn.csv")
TRY_Na<-read.csv("C:/Users/kotha020/Downloads/TRY_Na.csv")
TRY_P<-read.csv("C:/Users/kotha020/Downloads/TRY_P.csv")
TRY_Zn<-read.csv("C:/Users/kotha020/Downloads/TRY_Zn.csv")

TRY_SLA$LMA<-1/TRY_SLA$StdValue

####################################
## plot trait distributions

## find projects with more than 30 samples (all but CABO-General)
common.proj<-names(table(trait.df$project))[table(trait.df$project)>30]

ggplot(data=trait.df[trait.df$project %in% common.proj,],
       aes(x=chlA_mass,color=functional.group))+
  geom_density(size=1.25)+
  theme_bw()+
  scale_color_brewer(palette="Set3")

LMA_KL<-list()
LDMC_KL<-list()
Nmass_KL<-list()
hemicellulose_mass_KL<-list()
lignin_mass_KL<-list()
chlA_mass_KL<-list()
EWT_KL<-list()

for(i in common.proj){
  print(i)
  trait.df.proj<-trait.df[trait.df$project==i,]
  trait.df.other<-trait.df[trait.df$project!=i,]
  
  LMA_KL[[i]]<-KL.divergence(na.omit(trait.df.other$LMA),
                             na.omit(trait.df.proj$LMA))
  
  LDMC_KL[[i]]<-KL.divergence(na.omit(trait.df.other$LDMC),
                              na.omit(trait.df.proj$LDMC))

  Nmass_KL[[i]]<-KL.divergence(na.omit(trait.df.other$Nmass),
                               na.omit(trait.df.proj$Nmass))

  hemicellulose_mass_KL[[i]]<-KL.divergence(na.omit(trait.df.other$hemicellulose_mass),
                                     na.omit(trait.df.proj$hemicellulose_mass))
  
  lignin_mass_KL[[i]]<-KL.divergence(na.omit(trait.df.other$lignin_mass),
                                     na.omit(trait.df.proj$lignin_mass))
  
  chlA_mass_KL[[i]]<-KL.divergence(na.omit(trait.df.other$chlA_mass),
                                   na.omit(trait.df.proj$chlA_mass))
  
  EWT_KL[[i]]<-KL.divergence(na.omit(trait.df.other$EWT),
                             na.omit(trait.df.proj$EWT))
  
}

## LMA density plot
LMA_diff<-common.proj[which.max(unlist(lapply(LMA_KL,function(x) x[[5]])))]
trait.df$LMA_cat<-"CABO"
# trait.df$LMA_cat[trait.df$project==LMA_diff]<-LMA_diff
LMA_sub<-trait.df[trait.df$project %in% common.proj,
                  c("sample_id","LMA","LMA_cat")]

TRY_SLA_sub<-data.frame(sample_id=TRY_SLA$ObsDataID,
                        LMA=TRY_SLA$LMA,
                        LMA_cat="TRY")
all.LMA<-rbind(LMA_sub,TRY_SLA_sub)

LMA_density<-ggplot(data=all.LMA,
                    aes(x=LMA,color=LMA_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression("LMA (kg m"^-2*")"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$LMA)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_SLA$LMA)),")",sep="")))

## EWT density plot
EWT_diff<-common.proj[which.max(unlist(lapply(EWT_KL,function(x) x[[5]])))]
trait.df$EWT_cat<-"CABO"
# trait.df$EWT_cat[trait.df$project==EWT_diff]<-EWT_diff
EWT_sub<-trait.df[trait.df$project %in% common.proj,
                   c("sample_id","EWT","EWT_cat")]

EWT_density<-ggplot(data=EWT_sub,
                     aes(x=EWT,color=EWT_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x="EWT (cm)")+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$EWT)),")",sep="")))

## LDMC density plot
LDMC_diff<-common.proj[which.max(unlist(lapply(LDMC_KL,function(x) x[[5]])))]
trait.df$LDMC_cat<-"CABO"
# trait.df$LDMC_cat[trait.df$project==LDMC_diff]<-LDMC_diff
LDMC_sub<-trait.df[trait.df$project %in% common.proj,
                  c("sample_id","LDMC","LDMC_cat")]

TRY_LDMC_sub<-data.frame(sample_id=TRY_LDMC$ObsDataID,
                        LDMC=TRY_LDMC$StdValue*1000,
                        LDMC_cat="TRY")
all.LDMC<-rbind(LDMC_sub,TRY_LDMC_sub)

LDMC_density<-ggplot(data=all.LDMC,
                     aes(x=LDMC,color=LDMC_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.75,0.7))+
  labs(y="Density",x=expression("LDMC (mg g"^-1*")"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$LDMC)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_LDMC$StdValue)),")",sep="")))

## Nmass density plot
Nmass_diff<-common.proj[which.max(unlist(lapply(Nmass_KL,function(x) x[[5]])))]
trait.df$Nmass_cat<-"CABO"
# trait.df$Nmass_cat[trait.df$project==Nmass_diff]<-Nmass_diff
Nmass_sub<-trait.df[trait.df$project %in% common.proj,
                   c("sample_id","Nmass","Nmass_cat")]

TRY_N_sub<-data.frame(sample_id=TRY_N$ObsDataID,
                         Nmass=TRY_N$StdValue/10,
                         Nmass_cat="TRY")
all.Nmass<-rbind(Nmass_sub,TRY_N_sub)

Nmass_density<-ggplot(data=all.Nmass,
                     aes(x=Nmass,color=Nmass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression("N"[mass]~" (%)"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Nmass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_N$StdValue)),")",sep="")))

## Cmass density plot
# Cmass_diff<-common.proj[which.max(unlist(lapply(Cmass_KL,function(x) x[[5]])))]
trait.df$Cmass_cat<-"CABO"
# trait.df$Cmass_cat[trait.df$project==Cmass_diff]<-Cmass_diff
Cmass_sub<-trait.df[trait.df$project %in% common.proj,
                    c("sample_id","Cmass","Cmass_cat")]

TRY_C_sub<-data.frame(sample_id=TRY_C$ObsDataID,
                      Cmass=TRY_C$StdValue/10,
                      Cmass_cat="TRY")
all.Cmass<-rbind(Cmass_sub,TRY_C_sub)

Cmass_density<-ggplot(data=all.Cmass,
                      aes(x=Cmass,color=Cmass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.25,0.85))+
  labs(y="Density",x=expression("C"[mass]~" (%)"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Cmass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_C$StdValue)),")",sep="")))

## solubles density plot
# solubles_mass_diff<-common.proj[which.max(unlist(lapply(solubles_mass_KL,function(x) x[[5]])))]
trait.df$solubles_mass_cat<-"CABO"
# trait.df$solubles_mass_cat[trait.df$project==solubles_mass_diff]<-solubles_mass_diff
solubles_mass_sub<-trait.df[trait.df$project %in% common.proj,
                                 c("sample_id","solubles_mass","solubles_mass_cat")]

TRY_solubles_sub<-data.frame(sample_id=TRY_solubles$ObsDataID,
                                  solubles_mass=TRY_solubles$StdValue/10,
                                  solubles_mass_cat="TRY")
all.solubles_mass<-rbind(solubles_mass_sub,TRY_solubles_sub)

solubles_mass_density<-ggplot(data=all.solubles_mass,
                                   aes(x=solubles_mass,color=solubles_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.25,0.85))+
  labs(y="Density",x=expression("Solubles (%)"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$solubles_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_solubles$StdValue)),")",sep="")))

## hemicellulose density plot
hemicellulose_mass_diff<-common.proj[which.max(unlist(lapply(hemicellulose_mass_KL,function(x) x[[5]])))]
trait.df$hemicellulose_mass_cat<-"CABO"
# trait.df$hemicellulose_mass_cat[trait.df$project==hemicellulose_mass_diff]<-hemicellulose_mass_diff
hemicellulose_mass_sub<-trait.df[trait.df$project %in% common.proj,
                    c("sample_id","hemicellulose_mass","hemicellulose_mass_cat")]

TRY_hemicellulose_sub<-data.frame(sample_id=TRY_hemicellulose$ObsDataID,
                      hemicellulose_mass=TRY_hemicellulose$StdValue/10,
                      hemicellulose_mass_cat="TRY")
all.hemicellulose_mass<-rbind(hemicellulose_mass_sub,TRY_hemicellulose_sub)

hemicellulose_mass_density<-ggplot(data=all.hemicellulose_mass,
                      aes(x=hemicellulose_mass,color=hemicellulose_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression("Hemicellulose (%)"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$hemicellulose_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_hemicellulose$StdValue)),")",sep="")))

## cellulose density plot
# cellulose_mass_diff<-common.proj[which.max(unlist(lapply(cellulose_mass_KL,function(x) x[[5]])))]
trait.df$cellulose_mass_cat<-"CABO"
# trait.df$cellulose_mass_cat[trait.df$project==cellulose_mass_diff]<-cellulose_mass_diff
cellulose_mass_sub<-trait.df[trait.df$project %in% common.proj,
                                 c("sample_id","cellulose_mass","cellulose_mass_cat")]

TRY_cellulose_sub<-data.frame(sample_id=TRY_cellulose$ObsDataID,
                                  cellulose_mass=TRY_cellulose$StdValue/10,
                                  cellulose_mass_cat="TRY")
all.cellulose_mass<-rbind(cellulose_mass_sub,TRY_cellulose_sub)

cellulose_mass_density<-ggplot(data=all.cellulose_mass,
                                   aes(x=cellulose_mass,color=cellulose_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression("Cellulose (%)"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$cellulose_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_cellulose$StdValue)),")",sep="")))

## lignin density plot
lignin_mass_diff<-common.proj[which.max(unlist(lapply(lignin_mass_KL,function(x) x[[5]])))]
trait.df$lignin_mass_cat<-"CABO"
# trait.df$lignin_mass_cat[trait.df$project==lignin_mass_diff]<-lignin_mass_diff
lignin_mass_sub<-trait.df[trait.df$project %in% common.proj,
                                 c("sample_id","lignin_mass","lignin_mass_cat")]

TRY_lignin_sub<-data.frame(sample_id=TRY_lignin$ObsDataID,
                                  lignin_mass=TRY_lignin$StdValue/10,
                                  lignin_mass_cat="TRY")
all.lignin_mass<-rbind(lignin_mass_sub,TRY_lignin_sub)

lignin_mass_density<-ggplot(data=all.lignin_mass,
                                   aes(x=lignin_mass,color=lignin_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression("Lignin (%)"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$lignin_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_lignin$StdValue)),")",sep="")))

## chlA density plot
# chlA_mass_diff<-common.proj[which.max(unlist(lapply(chlA_mass_KL,function(x) x[[5]])))]
trait.df$chlA_mass_cat<-"CABO"
# trait.df$chlA_mass_cat[trait.df$project==chlA_mass_diff]<-chlA_mass_diff
chlA_mass_sub<-trait.df[trait.df$project %in% common.proj,
                                 c("sample_id","chlA_mass","chlA_mass_cat")]

TRY_ChlA_sub<-data.frame(sample_id=TRY_ChlA$ObsDataID,
                                  chlA_mass=TRY_ChlA$StdValue,
                                  chlA_mass_cat="TRY")
all.chlA_mass<-rbind(chlA_mass_sub,TRY_ChlA_sub)

chlA_mass_density<-ggplot(data=all.chlA_mass,
                                   aes(x=chlA_mass,color=chlA_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Chl ",italic("a")," (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$chlA_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_ChlA$StdValue)),")",sep="")))

## chlB density plot
# chlB_mass_diff<-common.proj[which.max(unlist(lapply(chlB_mass_KL,function(x) x[[5]])))]
trait.df$chlB_mass_cat<-"CABO"
# trait.df$chlB_mass_cat[trait.df$project==chlB_mass_diff]<-chlB_mass_diff
chlB_mass_sub<-trait.df[trait.df$project %in% common.proj,
                        c("sample_id","chlB_mass","chlB_mass_cat")]

TRY_ChlB_sub<-data.frame(sample_id=TRY_ChlB$ObsDataID,
                         chlB_mass=TRY_ChlB$StdValue,
                         chlB_mass_cat="TRY")
all.chlB_mass<-rbind(chlB_mass_sub,TRY_ChlB_sub)

chlB_mass_density<-ggplot(data=all.chlB_mass,
                          aes(x=chlB_mass,color=chlB_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Chl ",italic("b")," (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$chlB_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_ChlB$StdValue)),")",sep="")))

## car density plot
# car_mass_diff<-common.proj[which.max(unlist(lapply(car_mass_KL,function(x) x[[5]])))]
trait.df$car_mass_cat<-"CABO"
# trait.df$car_mass_cat[trait.df$project==car_mass_diff]<-car_mass_diff
car_mass_sub<-trait.df[trait.df$project %in% common.proj,
                        c("sample_id","car_mass","car_mass_cat")]

TRY_Car_sub<-data.frame(sample_id=TRY_Car$ObsDataID,
                         car_mass=TRY_Car$StdValue,
                         car_mass_cat="TRY")
all.car_mass<-rbind(car_mass_sub,TRY_Car_sub)

car_mass_density<-ggplot(data=all.car_mass,
                          aes(x=car_mass,color=car_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression("Carotenoids (mg g"^-1*")"))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$car_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Car$StdValue)),")",sep="")))

## Al density plot
# Al_mass_diff<-common.proj[which.max(unlist(lapply(Al_mass_KL,function(x) x[[5]])))]
trait.df$Al_mass_cat<-"CABO"
# trait.df$Al_mass_cat[trait.df$project==Al_mass_diff]<-Al_mass_diff
Al_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","Al_mass","Al_mass_cat")]

TRY_Al_sub<-data.frame(sample_id=TRY_Al$ObsDataID,
                        Al_mass=TRY_Al$StdValue,
                        Al_mass_cat="TRY")
all.Al_mass<-rbind(Al_mass_sub,TRY_Al_sub)

Al_mass_density<-ggplot(data=all.Al_mass,
                         aes(x=Al_mass,color=Al_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Al (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Al_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Al$StdValue)),")",sep="")))+
  coord_cartesian(xlim=c(0,1))

## Ca density plot
# Ca_mass_diff<-common.proj[which.max(unlist(lapply(Ca_mass_KL,function(x) x[[5]])))]
trait.df$Ca_mass_cat<-"CABO"
# trait.df$Ca_mass_cat[trait.df$project==Ca_mass_diff]<-Ca_mass_diff
Ca_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","Ca_mass","Ca_mass_cat")]

TRY_Ca_sub<-data.frame(sample_id=TRY_Ca$ObsDataID,
                        Ca_mass=TRY_Ca$StdValue,
                        Ca_mass_cat="TRY")
all.Ca_mass<-rbind(Ca_mass_sub,TRY_Ca_sub)

Ca_mass_density<-ggplot(data=all.Ca_mass,
                         aes(x=Ca_mass,color=Ca_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Ca (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Ca_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Ca$StdValue)),")",sep="")))

## Cu density plot
# Cu_mass_diff<-common.proj[which.max(unlist(lapply(Cu_mass_KL,function(x) x[[5]])))]
trait.df$Cu_mass_cat<-"CABO"
# trait.df$Cu_mass_cat[trait.df$project==Cu_mass_diff]<-Cu_mass_diff
Cu_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","Cu_mass","Cu_mass_cat")]

TRY_Cu_sub<-data.frame(sample_id=TRY_Cu$ObsDataID,
                        Cu_mass=TRY_Cu$StdValue,
                        Cu_mass_cat="TRY")
all.Cu_mass<-rbind(Cu_mass_sub,TRY_Cu_sub)

Cu_mass_density<-ggplot(data=all.Cu_mass,
                         aes(x=Cu_mass,color=Cu_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Cu (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Cu_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Cu$StdValue)),")",sep="")))

## Fe density plot
# Fe_mass_diff<-common.proj[which.max(unlist(lapply(Fe_mass_KL,function(x) x[[5]])))]
trait.df$Fe_mass_cat<-"CABO"
# trait.df$Fe_mass_cat[trait.df$project==Fe_mass_diff]<-Fe_mass_diff
Fe_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","Fe_mass","Fe_mass_cat")]

TRY_Fe_sub<-data.frame(sample_id=TRY_Fe$ObsDataID,
                        Fe_mass=TRY_Fe$StdValue,
                        Fe_mass_cat="TRY")
all.Fe_mass<-rbind(Fe_mass_sub,TRY_Fe_sub)

Fe_mass_density<-ggplot(data=all.Fe_mass,
                         aes(x=Fe_mass,color=Fe_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Fe (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Fe_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Fe$StdValue)),")",sep="")))

## K density plot
# K_mass_diff<-common.proj[which.max(unlist(lapply(K_mass_KL,function(x) x[[5]])))]
trait.df$K_mass_cat<-"CABO"
# trait.df$K_mass_cat[trait.df$project==K_mass_diff]<-K_mass_diff
K_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","K_mass","K_mass_cat")]

TRY_K_sub<-data.frame(sample_id=TRY_K$ObsDataID,
                        K_mass=TRY_K$StdValue,
                        K_mass_cat="TRY")
all.K_mass<-rbind(K_mass_sub,TRY_K_sub)

K_mass_density<-ggplot(data=all.K_mass,
                         aes(x=K_mass,color=K_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("K (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$K_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_K$StdValue)),")",sep="")))

## Mg density plot
# Mg_mass_diff<-common.proj[which.max(unlist(lapply(Mg_mass_KL,function(x) x[[5]])))]
trait.df$Mg_mass_cat<-"CABO"
# trait.df$Mg_mass_cat[trait.df$project==Mg_mass_diff]<-Mg_mass_diff
Mg_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","Mg_mass","Mg_mass_cat")]

TRY_Mg_sub<-data.frame(sample_id=TRY_Mg$ObsDataID,
                        Mg_mass=TRY_Mg$StdValue,
                        Mg_mass_cat="TRY")
all.Mg_mass<-rbind(Mg_mass_sub,TRY_Mg_sub)

Mg_mass_density<-ggplot(data=all.Mg_mass,
                         aes(x=Mg_mass,color=Mg_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Mg (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Mg_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Mg$StdValue)),")",sep="")))

## Mn density plot
# Mn_mass_diff<-common.proj[which.max(unlist(lapply(Mn_mass_KL,function(x) x[[5]])))]
trait.df$Mn_mass_cat<-"CABO"
# trait.df$Mn_mass_cat[trait.df$project==Mn_mass_diff]<-Mn_mass_diff
Mn_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","Mn_mass","Mn_mass_cat")]

TRY_Mn_sub<-data.frame(sample_id=TRY_Mn$ObsDataID,
                        Mn_mass=TRY_Mn$StdValue,
                        Mn_mass_cat="TRY")
all.Mn_mass<-rbind(Mn_mass_sub,TRY_Mn_sub)

Mn_mass_density<-ggplot(data=all.Mn_mass,
                         aes(x=Mn_mass,color=Mn_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Mn (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Mn_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Mn$StdValue)),")",sep="")))

## Na density plot
# Na_mass_diff<-common.proj[which.max(unlist(lapply(Na_mass_KL,function(x) x[[5]])))]
trait.df$Na_mass_cat<-"CABO"
# trait.df$Na_mass_cat[trait.df$project==Na_mass_diff]<-Na_mass_diff
Na_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","Na_mass","Na_mass_cat")]

TRY_Na_sub<-data.frame(sample_id=TRY_Na$ObsDataID,
                        Na_mass=TRY_Na$StdValue,
                        Na_mass_cat="TRY")
all.Na_mass<-rbind(Na_mass_sub,TRY_Na_sub)

Na_mass_density<-ggplot(data=all.Na_mass,
                         aes(x=Na_mass,color=Na_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Na (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Na_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Na$StdValue)),")",sep="")))+
  coord_cartesian(xlim=c(0,15))

## P density plot
# P_mass_diff<-common.proj[which.max(unlist(lapply(P_mass_KL,function(x) x[[5]])))]
trait.df$P_mass_cat<-"CABO"
# trait.df$P_mass_cat[trait.df$project==P_mass_diff]<-P_mass_diff
P_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","P_mass","P_mass_cat")]

TRY_P_sub<-data.frame(sample_id=TRY_P$ObsDataID,
                        P_mass=TRY_P$StdValue,
                        P_mass_cat="TRY")
all.P_mass<-rbind(P_mass_sub,TRY_P_sub)

P_mass_density<-ggplot(data=all.P_mass,
                         aes(x=P_mass,color=P_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("P (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$P_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_P$StdValue)),")",sep="")))

## Zn density plot
# Zn_mass_diff<-common.proj[which.max(unlist(lapply(Zn_mass_KL,function(x) x[[5]])))]
trait.df$Zn_mass_cat<-"CABO"
# trait.df$Zn_mass_cat[trait.df$project==Zn_mass_diff]<-Zn_mass_diff
Zn_mass_sub<-trait.df[trait.df$project %in% common.proj,
                       c("sample_id","Zn_mass","Zn_mass_cat")]

TRY_Zn_sub<-data.frame(sample_id=TRY_Zn$ObsDataID,
                        Zn_mass=TRY_Zn$StdValue,
                        Zn_mass_cat="TRY")
all.Zn_mass<-rbind(Zn_mass_sub,TRY_Zn_sub)

Zn_mass_density<-ggplot(data=all.Zn_mass,
                         aes(x=Zn_mass,color=Zn_mass_cat))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7))+
  labs(y="Density",x=expression(paste("Zn (mg g"^-1*")")))+
  scale_color_discrete(labels=c(paste("CABO (n=",sum(!is.na(trait.df$Zn_mass)),")",sep=""),
                                paste("TRY (n=",sum(!is.na(TRY_Zn$StdValue)),")",sep="")))

pdf("Images/density_plot1.pdf",width=24,height=30)
ggarrange(plotlist = list(LMA_density,EWT_density,LDMC_density,
                          Nmass_density,Cmass_density,
                          solubles_mass_density,hemicellulose_mass_density,
                          cellulose_mass_density,lignin_mass_density,
                          chlA_mass_density,chlB_mass_density,
                          car_mass_density),
          ncol = 3,nrow=4)
dev.off()

pdf("Images/density_plot2.pdf",width=24,height=30)
ggarrange(plotlist = list(Al_mass_density,
                          Ca_mass_density,Cu_mass_density,
                          Fe_mass_density,K_mass_density,
                          Mg_mass_density,Mn_mass_density,
                          Na_mass_density,P_mass_density,
                          Zn_mass_density),
          ncol = 3,nrow=4)
dev.off()

############################
## variance partitioning

variances<-data.frame(effect=c("species","genus","site","family","residual"))

varpart_LMA<-lmer(LMA~1+(1|family/genus/species)+(1|site),data=trait.df)
variances$LMA<-as.data.frame(VarCorr(varpart_LMA))$vcov

varpart_LDMC<-lmer(LDMC~1+(1|family/genus/species)+(1|site),data=trait.df)
variances$LDMC<-as.data.frame(VarCorr(varpart_LDMC))$vcov

varpart_Nmass<-lmer(Nmass~1+(1|family/genus/species)+(1|site),data=trait.df)
variances$Nmass<-as.data.frame(VarCorr(varpart_Nmass))$vcov

varpart_EWT<-lmer(EWT~1+(1|family/genus/species)+(1|site),data=trait.df)
variances$EWT<-as.data.frame(VarCorr(varpart_EWT))$vcov

varpart_chlA_mass<-lmer(chlA_mass~1+(1|family/genus/species)+(1|site),data=trait.df)
variances$chlA_mass<-as.data.frame(VarCorr(varpart_chlA_mass))$vcov

varpart_hemicellulose_mass<-lmer(hemicellulose_mass~1+(1|family/genus/species)+(1|site),data=trait.df)
variances$hemicellulose_mass<-as.data.frame(VarCorr(varpart_hemicellulose_mass))$vcov

varpart_lignin_mass<-lmer(lignin_mass~1+(1|family/genus/species)+(1|site),data=trait.df)
variances$lignin_mass<-as.data.frame(VarCorr(varpart_lignin_mass))$vcov

variances_long<-gather(variances, trait, variance, LMA:lignin_mass, factor_key=TRUE)

ggplot(variances_long, aes(x = trait, y = variance,fill=effect)) + 
  geom_bar(position = "fill", stat="identity")+
  scale_y_continuous(labels = scales::percent_format())
  