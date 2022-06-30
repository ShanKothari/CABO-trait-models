setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(caret)
library(pls)
library(reshape2)
library(statip)
library(lme4)
library(FNN)
library(ggplot2)
library(patchwork)
source("Scripts/CABO-trait-models/00 useful_functions.R")

spec.traits<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")
## the Pardo dataset has no trait data (yet)
spec.traits<-spec.traits[-which(meta(spec.traits)$project=="2019-Pardo-MSc-UdeM")]

## third color (for ferns) removed
colorBlind  <- c("#E69F00","#009E73","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

##########################################
## model transferral among functional groups

LMA_pred_list<-list()
EWT_pred_list<-list()
LDMC_pred_list<-list()
Nmass_pred_list<-list()
Cmass_pred_list<-list()
solubles_mass_pred_list<-list()
hemicellulose_mass_pred_list<-list()
cellulose_mass_pred_list<-list()
lignin_mass_pred_list<-list()
chlA_mass_pred_list<-list()
chlB_mass_pred_list<-list()
car_mass_pred_list<-list()

LMA_dist<-list()
EWT_dist<-list()
LDMC_dist<-list()
Nmass_dist<-list()
Cmass_dist<-list()
solubles_mass_dist<-list()
hemicellulose_mass_dist<-list()
cellulose_mass_dist<-list()
lignin_mass_dist<-list()
chlA_mass_dist<-list()
chlB_mass_dist<-list()
car_mass_dist<-list()

LMA_KL<-list()
EWT_KL<-list()
LDMC_KL<-list()
Nmass_KL<-list()
Cmass_KL<-list()
solubles_mass_KL<-list()
hemicellulose_mass_KL<-list()
cellulose_mass_KL<-list()
lignin_mass_KL<-list()
chlA_mass_KL<-list()
chlB_mass_KL<-list()
car_mass_KL<-list()

common.fg<-names(table(meta(spec.traits)$functional.group))[table(meta(spec.traits)$functional.group)>30]

for(i in common.fg){
  print(i)
  spec.traits.fg<-spec.traits[meta(spec.traits)$functional.group==i]
  spec.traits.other<-spec.traits[meta(spec.traits)$functional.group!=i]
  
  LMA_cal<-plsr(meta(spec.traits.other)$LMA~as.matrix(spec.traits.other),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_LMA<-selectNcomp(LMA_cal, method = "onesigma", plot = F)
  LMA_val<-predict(LMA_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_LMA)
  LMA_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$LMA,
                                 val_pred=as.vector(LMA_val))
  LMA_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$LMA),
                           na.omit(meta(spec.traits.other)$LMA),
                           lower=-Inf,
                           upper=Inf)
  LMA_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$LMA),
                             na.omit(meta(spec.traits.fg)$LMA))
  
  EWT_cal<-plsr(meta(spec.traits.other)$EWT~as.matrix(spec.traits.other),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_EWT<-selectNcomp(EWT_cal, method = "onesigma", plot = F)
  EWT_val<-predict(EWT_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_EWT)
  EWT_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$EWT,
                                 val_pred=as.vector(EWT_val))
  EWT_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$EWT),
                           na.omit(meta(spec.traits.other)$EWT),
                           lower=-Inf,
                           upper=Inf)
  EWT_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$EWT),
                             na.omit(meta(spec.traits.fg)$EWT))
  
  LDMC_cal<-plsr(meta(spec.traits.other)$LDMC~as.matrix(spec.traits.other),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_LDMC<-selectNcomp(LDMC_cal, method = "onesigma", plot = F)
  LDMC_val<-predict(LDMC_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_LDMC)
  LDMC_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$LDMC,
                                  val_pred=as.vector(LDMC_val))
  LDMC_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$LDMC),
                            na.omit(meta(spec.traits.other)$LDMC),
                            lower=-Inf,
                            upper=Inf)
  LDMC_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$LDMC),
                             na.omit(meta(spec.traits.fg)$LDMC))
  
  Nmass_cal<-plsr(meta(spec.traits.other)$Nmass~as.matrix(spec.traits.other),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_Nmass<-selectNcomp(Nmass_cal, method = "onesigma", plot = F)
  Nmass_val<-predict(Nmass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_Nmass)
  Nmass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$Nmass,
                                   val_pred=as.vector(Nmass_val))
  Nmass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$Nmass),
                             na.omit(meta(spec.traits.other)$Nmass),
                             lower=-Inf,
                             upper=Inf)
  Nmass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$Nmass),
                             na.omit(meta(spec.traits.fg)$Nmass))
  
  Cmass_cal<-plsr(meta(spec.traits.other)$Cmass~as.matrix(spec.traits.other),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_Cmass<-selectNcomp(Cmass_cal, method = "onesigma", plot = F)
  Cmass_val<-predict(Cmass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_Cmass)
  Cmass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$Cmass,
                                   val_pred=as.vector(Cmass_val))
  Cmass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$Cmass),
                             na.omit(meta(spec.traits.other)$Cmass),
                             lower=0,
                             upper=100)
  Cmass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$Cmass),
                             na.omit(meta(spec.traits.fg)$Cmass))
  
  solubles_mass_cal<-plsr(meta(spec.traits.other)$solubles_mass~as.matrix(spec.traits.other),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_solubles_mass<-selectNcomp(solubles_mass_cal, method = "onesigma", plot = F)
  solubles_mass_val<-predict(solubles_mass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_solubles_mass)
  solubles_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$solubles_mass,
                                         val_pred=as.vector(solubles_mass_val))
  solubles_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$solubles_mass),
                                   na.omit(meta(spec.traits.other)$solubles_mass),
                                   lower=-Inf,
                                   upper=Inf)
  solubles_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$solubles_mass),
                             na.omit(meta(spec.traits.fg)$solubles_mass))
  
  hemicellulose_mass_cal<-plsr(meta(spec.traits.other)$hemicellulose_mass~as.matrix(spec.traits.other),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_hemicellulose_mass<-selectNcomp(hemicellulose_mass_cal, method = "onesigma", plot = F)
  hemicellulose_mass_val<-predict(hemicellulose_mass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_hemicellulose_mass)
  hemicellulose_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$hemicellulose_mass,
                                         val_pred=as.vector(hemicellulose_mass_val))
  hemicellulose_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$hemicellulose_mass),
                                   na.omit(meta(spec.traits.other)$hemicellulose_mass),
                                   lower=-Inf,
                                   upper=Inf)
  hemicellulose_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$hemicellulose_mass),
                             na.omit(meta(spec.traits.fg)$hemicellulose_mass))
  
  cellulose_mass_cal<-plsr(meta(spec.traits.other)$cellulose_mass~as.matrix(spec.traits.other),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_cellulose_mass<-selectNcomp(cellulose_mass_cal, method = "onesigma", plot = F)
  cellulose_mass_val<-predict(cellulose_mass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_cellulose_mass)
  cellulose_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$cellulose_mass,
                                         val_pred=as.vector(cellulose_mass_val))
  cellulose_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$cellulose_mass),
                                   na.omit(meta(spec.traits.other)$cellulose_mass),
                                   lower=-Inf,
                                   upper=Inf)
  cellulose_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$cellulose_mass),
                             na.omit(meta(spec.traits.fg)$cellulose_mass))
  
  lignin_mass_cal<-plsr(meta(spec.traits.other)$lignin_mass~as.matrix(spec.traits.other),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_lignin_mass<-selectNcomp(lignin_mass_cal, method = "onesigma", plot = F)
  lignin_mass_val<-predict(lignin_mass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_lignin_mass)
  lignin_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$lignin_mass,
                                         val_pred=as.vector(lignin_mass_val))
  lignin_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$lignin_mass),
                                   na.omit(meta(spec.traits.other)$lignin_mass),
                                   lower=Inf,
                                   upper=Inf)
  lignin_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$lignin_mass),
                             na.omit(meta(spec.traits.fg)$lignin_mass))
  
  chlA_mass_cal<-plsr(meta(spec.traits.other)$chlA_mass~as.matrix(spec.traits.other),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_chlA_mass<-selectNcomp(chlA_mass_cal, method = "onesigma", plot = F)
  chlA_mass_val<-predict(chlA_mass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_chlA_mass)
  chlA_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$chlA_mass,
                                        val_pred=as.vector(chlA_mass_val))
  chlA_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$chlA_mass),
                                  na.omit(meta(spec.traits.other)$chlA_mass),
                                  lower=-Inf,
                                  upper=Inf)
  chlA_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$chlA_mass),
                             na.omit(meta(spec.traits.fg)$chlA_mass))
  
  chlB_mass_cal<-plsr(meta(spec.traits.other)$chlB_mass~as.matrix(spec.traits.other),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_chlB_mass<-selectNcomp(chlB_mass_cal, method = "onesigma", plot = F)
  chlB_mass_val<-predict(chlB_mass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_chlB_mass)
  chlB_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$chlB_mass,
                                       val_pred=as.vector(chlB_mass_val))
  chlB_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$chlB_mass),
                                 na.omit(meta(spec.traits.other)$chlB_mass),
                                 lower=-Inf,
                                 upper=Inf)
  chlB_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$chlB_mass),
                             na.omit(meta(spec.traits.fg)$chlB_mass))
  
  car_mass_cal<-plsr(meta(spec.traits.other)$car_mass~as.matrix(spec.traits.other),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_car_mass<-selectNcomp(car_mass_cal, method = "onesigma", plot = F)
  car_mass_val<-predict(car_mass_cal,newdata=as.matrix(spec.traits.fg),ncomp=ncomp_car_mass)
  car_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.fg)$car_mass,
                                       val_pred=as.vector(car_mass_val))
  car_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.fg)$car_mass),
                                 na.omit(meta(spec.traits.other)$car_mass),
                                 lower=-Inf,
                                 upper=Inf)
  car_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$car_mass),
                             na.omit(meta(spec.traits.fg)$car_mass))
  
}

###############################
## measured vs predicted transfer plots

LMA_named_pred_list<-mapply(`[<-`, LMA_pred_list, 'functional.group', value = names(LMA_pred_list), SIMPLIFY = FALSE)
LMA_pred_df<-do.call(rbind,LMA_named_pred_list)
LMA_lower<-min(c(LMA_pred_df$measured,LMA_pred_df$val_pred),na.rm=T)-0.01
LMA_upper<-max(c(LMA_pred_df$measured,LMA_pred_df$val_pred),na.rm=T)+0.01
LMA_transfer_plot<-ggplot(LMA_pred_df,
                         aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

LDMC_named_pred_list<-mapply(`[<-`, LDMC_pred_list, 'functional.group', value = names(LDMC_pred_list), SIMPLIFY = FALSE)
LDMC_pred_df<-do.call(rbind,LDMC_named_pred_list)
LDMC_lower<-min(c(LDMC_pred_df$measured,LDMC_pred_df$val_pred),na.rm=T)-0.01
LDMC_upper<-max(c(LDMC_pred_df$measured,LDMC_pred_df$val_pred),na.rm=T)+0.01
LDMC_transfer_plot<-ggplot(LDMC_pred_df,
                          aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

EWT_named_pred_list<-mapply(`[<-`, EWT_pred_list, 'functional.group', value = names(EWT_pred_list), SIMPLIFY = FALSE)
EWT_pred_df<-do.call(rbind,EWT_named_pred_list)
EWT_lower<-min(c(EWT_pred_df$measured,EWT_pred_df$val_pred),na.rm=T)-0.01
EWT_upper<-max(c(EWT_pred_df$measured,EWT_pred_df$val_pred),na.rm=T)+0.01
EWT_transfer_plot<-ggplot(EWT_pred_df,
                          aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=20))+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Nmass_named_pred_list<-mapply(`[<-`, Nmass_pred_list, 'functional.group', value = names(Nmass_pred_list), SIMPLIFY = FALSE)
Nmass_pred_df<-do.call(rbind,Nmass_named_pred_list)
Nmass_lower<-min(c(Nmass_pred_df$measured,Nmass_pred_df$val_pred),na.rm=T)-0.01
Nmass_upper<-max(c(Nmass_pred_df$measured,Nmass_pred_df$val_pred),na.rm=T)+0.01
Nmass_transfer_plot<-ggplot(Nmass_pred_df,
                          aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Nmass_lower,Nmass_upper),ylim=c(Nmass_lower,Nmass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cmass_named_pred_list<-mapply(`[<-`, Cmass_pred_list, 'functional.group', value = names(Cmass_pred_list), SIMPLIFY = FALSE)
Cmass_pred_df<-do.call(rbind,Cmass_named_pred_list)
Cmass_lower<-min(c(Cmass_pred_df$measured,Cmass_pred_df$val_pred),na.rm=T)-0.01
Cmass_upper<-max(c(Cmass_pred_df$measured,Cmass_pred_df$val_pred),na.rm=T)+0.01
Cmass_transfer_plot<-ggplot(Cmass_pred_df,
                            aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cmass_lower,Cmass_upper),ylim=c(Cmass_lower,Cmass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

solubles_mass_named_pred_list<-mapply(`[<-`, solubles_mass_pred_list, 'functional.group', value = names(solubles_mass_pred_list), SIMPLIFY = FALSE)
solubles_mass_pred_df<-do.call(rbind,solubles_mass_named_pred_list)
solubles_mass_lower<-min(c(solubles_mass_pred_df$measured,solubles_mass_pred_df$val_pred),na.rm=T)-0.01
solubles_mass_upper<-max(c(solubles_mass_pred_df$measured,solubles_mass_pred_df$val_pred),na.rm=T)+0.01
solubles_mass_transfer_plot<-ggplot(solubles_mass_pred_df,
                            aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_mass_lower,solubles_mass_upper),
                  ylim=c(solubles_mass_lower,solubles_mass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured solubles (%)"),
       x=expression("Predicted solubles (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

hemicellulose_mass_named_pred_list<-mapply(`[<-`, hemicellulose_mass_pred_list, 'functional.group', value = names(hemicellulose_mass_pred_list), SIMPLIFY = FALSE)
hemicellulose_mass_pred_df<-do.call(rbind,hemicellulose_mass_named_pred_list)
hemicellulose_mass_lower<-min(c(hemicellulose_mass_pred_df$measured,hemicellulose_mass_pred_df$val_pred),na.rm=T)-0.01
hemicellulose_mass_upper<-max(c(hemicellulose_mass_pred_df$measured,hemicellulose_mass_pred_df$val_pred),na.rm=T)+0.01
hemicellulose_mass_transfer_plot<-ggplot(hemicellulose_mass_pred_df,
                            aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_mass_lower,hemicellulose_mass_upper),
                  ylim=c(hemicellulose_mass_lower,hemicellulose_mass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured hemicellulose (%)"),
       x=expression("Predicted hemicellulose (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

cellulose_mass_named_pred_list<-mapply(`[<-`, cellulose_mass_pred_list, 'functional.group', value = names(cellulose_mass_pred_list), SIMPLIFY = FALSE)
cellulose_mass_pred_df<-do.call(rbind,cellulose_mass_named_pred_list)
cellulose_mass_lower<-min(c(cellulose_mass_pred_df$measured,cellulose_mass_pred_df$val_pred),na.rm=T)-0.01
cellulose_mass_upper<-max(c(cellulose_mass_pred_df$measured,cellulose_mass_pred_df$val_pred),na.rm=T)+0.01
cellulose_mass_transfer_plot<-ggplot(cellulose_mass_pred_df,
                            aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_mass_lower,cellulose_mass_upper),
                  ylim=c(cellulose_mass_lower,cellulose_mass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured cellulose (%)"),
       x=expression("Predicted cellulose (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

lignin_mass_named_pred_list<-mapply(`[<-`, lignin_mass_pred_list, 'functional.group', value = names(lignin_mass_pred_list), SIMPLIFY = FALSE)
lignin_mass_pred_df<-do.call(rbind,lignin_mass_named_pred_list)
lignin_mass_lower<-min(c(lignin_mass_pred_df$measured,lignin_mass_pred_df$val_pred),na.rm=T)-0.01
lignin_mass_upper<-max(c(lignin_mass_pred_df$measured,lignin_mass_pred_df$val_pred),na.rm=T)+0.01
lignin_mass_transfer_plot<-ggplot(lignin_mass_pred_df,
                            aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_mass_lower,lignin_mass_upper),
                  ylim=c(lignin_mass_lower,lignin_mass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured lignin (%)"),
       x=expression("Predicted lignin (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_mass_named_pred_list<-mapply(`[<-`, chlA_mass_pred_list, 'functional.group', value = names(chlA_mass_pred_list), SIMPLIFY = FALSE)
chlA_mass_pred_df<-do.call(rbind,chlA_mass_named_pred_list)
chlA_mass_lower<-min(c(chlA_mass_pred_df$measured,chlA_mass_pred_df$val_pred),na.rm=T)-0.01
chlA_mass_upper<-max(c(chlA_mass_pred_df$measured,chlA_mass_pred_df$val_pred),na.rm=T)+0.01
chlA_mass_transfer_plot<-ggplot(chlA_mass_pred_df,
                            aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_mass_lower,chlA_mass_upper),ylim=c(chlA_mass_lower,chlA_mass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlB_mass_named_pred_list<-mapply(`[<-`, chlB_mass_pred_list, 'functional.group', value = names(chlB_mass_pred_list), SIMPLIFY = FALSE)
chlB_mass_pred_df<-do.call(rbind,chlB_mass_named_pred_list)
chlB_mass_lower<-min(c(chlB_mass_pred_df$measured,chlB_mass_pred_df$val_pred),na.rm=T)-0.01
chlB_mass_upper<-max(c(chlB_mass_pred_df$measured,chlB_mass_pred_df$val_pred),na.rm=T)+0.01
chlB_mass_transfer_plot<-ggplot(chlB_mass_pred_df,
                                aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_mass_lower,chlB_mass_upper),ylim=c(chlB_mass_lower,chlB_mass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

car_mass_named_pred_list<-mapply(`[<-`, car_mass_pred_list, 'functional.group', value = names(car_mass_pred_list), SIMPLIFY = FALSE)
car_mass_pred_df<-do.call(rbind,car_mass_named_pred_list)
car_mass_lower<-min(c(car_mass_pred_df$measured,car_mass_pred_df$val_pred),na.rm=T)-0.01
car_mass_upper<-max(c(car_mass_pred_df$measured,car_mass_pred_df$val_pred),na.rm=T)+0.01
car_mass_transfer_plot<-ggplot(car_mass_pred_df,
                                aes(y=measured,x=val_pred,color=functional.group))+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_mass_lower,car_mass_upper),ylim=c(car_mass_lower,car_mass_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

pdf("Images/model_transfer.pdf",width=16,height=19)
(LMA_transfer_plot+LDMC_transfer_plot+EWT_transfer_plot)/
  (Nmass_transfer_plot+Cmass_transfer_plot+solubles_mass_transfer_plot)/
  (hemicellulose_mass_transfer_plot+cellulose_mass_transfer_plot+lignin_mass_transfer_plot)/
  (chlA_mass_transfer_plot+chlB_mass_transfer_plot+car_mass_transfer_plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

model_transfer_list<-list(LMA=LMA_pred_df,
                          LDMC=LDMC_pred_df,
                          EWT=EWT_pred_df,
                          N=Nmass_pred_df,
                          C=Cmass_pred_df,
                          sol=solubles_mass_pred_df,
                          hemi=hemicellulose_mass_pred_df,
                          cell=cellulose_mass_pred_df,
                          lign=lignin_mass_pred_df,
                          chlA=chlA_mass_pred_df,
                          chlB=chlB_mass_pred_df,
                          car=car_mass_pred_df)

##################################
## summary statistic plots

LMA_sum_df<-data.frame(functional.group=names(LMA_pred_list),
                       RMSD=unlist(lapply(LMA_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                       perRMSD=unlist(lapply(LMA_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                       SD=unlist(lapply(LMA_pred_list,function(x) sd(x$measured,na.rm=T))),
                       R2=unlist(lapply(LMA_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                       nsamp=unlist(lapply(LMA_pred_list,nrow)),
                       hellinger=unlist(LMA_dist),
                       KL=unlist(lapply(LMA_KL,function(x) x[[5]])),
                       trait="LMA")

EWT_sum_df<-data.frame(functional.group=names(EWT_pred_list),
                       RMSD=unlist(lapply(EWT_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                       perRMSD=unlist(lapply(EWT_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                       SD=unlist(lapply(EWT_pred_list,function(x) sd(x$measured,na.rm=T))),
                       R2=unlist(lapply(EWT_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                       nsamp=unlist(lapply(EWT_pred_list,nrow)),
                       hellinger=unlist(EWT_dist),
                       KL=unlist(lapply(EWT_KL,function(x) x[[5]])),
                       trait="EWT")

LDMC_sum_df<-data.frame(functional.group=names(LDMC_pred_list),
                        RMSD=unlist(lapply(LDMC_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                        perRMSD=unlist(lapply(LDMC_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                        SD=unlist(lapply(LDMC_pred_list,function(x) sd(x$measured,na.rm=T))),
                        R2=unlist(lapply(LDMC_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                        nsamp=unlist(lapply(LDMC_pred_list,nrow)),
                        hellinger=unlist(LDMC_dist),
                        KL=unlist(lapply(LDMC_KL,function(x) x[[5]])),
                        trait="LDMC")

Nmass_sum_df<-data.frame(functional.group=names(Nmass_pred_list),
                         RMSD=unlist(lapply(Nmass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                         perRMSD=unlist(lapply(Nmass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                         SD=unlist(lapply(Nmass_pred_list,function(x) sd(x$measured,na.rm=T))),
                         R2=unlist(lapply(Nmass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                         nsamp=unlist(lapply(Nmass_pred_list,nrow)),
                         hellinger=unlist(Nmass_dist),
                         KL=unlist(lapply(Nmass_KL,function(x) x[[5]])),
                         trait="N")

## not sure why KL divergence for C mass comes up infinite
## for K<8 or so
Cmass_sum_df<-data.frame(functional.group=names(Cmass_pred_list),
                         RMSD=unlist(lapply(Cmass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                         perRMSD=unlist(lapply(Cmass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                         SD=unlist(lapply(Cmass_pred_list,function(x) sd(x$measured,na.rm=T))),
                         R2=unlist(lapply(Cmass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                         nsamp=unlist(lapply(Cmass_pred_list,nrow)),
                         hellinger=unlist(Cmass_dist),
                         KL=unlist(lapply(Cmass_KL,function(x) x[[10]])),
                         trait="C")

solubles_mass_sum_df<-data.frame(functional.group=names(solubles_mass_pred_list),
                                 RMSD=unlist(lapply(solubles_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                                 perRMSD=unlist(lapply(solubles_mass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                                 SD=unlist(lapply(solubles_mass_pred_list,function(x) sd(x$measured,na.rm=T))),
                                 R2=unlist(lapply(solubles_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                                 nsamp=unlist(lapply(solubles_mass_pred_list,nrow)),
                                 hellinger=unlist(solubles_mass_dist),
                                 KL=unlist(lapply(solubles_mass_KL,function(x) x[[5]])),
                                 trait="solubles")

hemicellulose_mass_sum_df<-data.frame(functional.group=names(hemicellulose_mass_pred_list),
                                      RMSD=unlist(lapply(hemicellulose_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                                      perRMSD=unlist(lapply(hemicellulose_mass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                                      SD=unlist(lapply(hemicellulose_mass_pred_list,function(x) sd(x$measured,na.rm=T))),
                                      R2=unlist(lapply(hemicellulose_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                                      nsamp=unlist(lapply(hemicellulose_mass_pred_list,nrow)),
                                      hellinger=unlist(hemicellulose_mass_dist),
                                      KL=unlist(lapply(hemicellulose_mass_KL,function(x) x[[5]])),
                                      trait="hemicellulose")

cellulose_mass_sum_df<-data.frame(functional.group=names(cellulose_mass_pred_list),
                                  RMSD=unlist(lapply(cellulose_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                                  perRMSD=unlist(lapply(cellulose_mass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                                  SD=unlist(lapply(cellulose_mass_pred_list,function(x) sd(x$measured,na.rm=T))),
                                  R2=unlist(lapply(cellulose_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                                  nsamp=unlist(lapply(cellulose_mass_pred_list,nrow)),
                                  hellinger=unlist(cellulose_mass_dist),
                                  KL=unlist(lapply(cellulose_mass_KL,function(x) x[[5]])),
                                  trait="cellulose")

lignin_mass_sum_df<-data.frame(functional.group=names(lignin_mass_pred_list),
                               RMSD=unlist(lapply(lignin_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                               perRMSD=unlist(lapply(lignin_mass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                               SD=unlist(lapply(lignin_mass_pred_list,function(x) sd(x$measured,na.rm=T))),
                               R2=unlist(lapply(lignin_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                               nsamp=unlist(lapply(lignin_mass_pred_list,nrow)),
                               hellinger=unlist(lignin_mass_dist),
                               KL=unlist(lapply(lignin_mass_KL,function(x) x[[5]])),
                               trait="lignin")

chlA_mass_sum_df<-data.frame(functional.group=names(chlA_mass_pred_list),
                             RMSD=unlist(lapply(chlA_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                             perRMSD=unlist(lapply(chlA_mass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                             SD=unlist(lapply(chlA_mass_pred_list,function(x) sd(x$measured,na.rm=T))),
                             R2=unlist(lapply(chlA_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                             nsamp=unlist(lapply(chlA_mass_pred_list,nrow)),
                             hellinger=unlist(chlA_mass_dist),
                             KL=unlist(lapply(chlA_mass_KL,function(x) x[[5]])),
                             trait="Chl a")

chlB_mass_sum_df<-data.frame(functional.group=names(chlB_mass_pred_list),
                             RMSD=unlist(lapply(chlB_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                             perRMSD=unlist(lapply(chlB_mass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                             SD=unlist(lapply(chlB_mass_pred_list,function(x) sd(x$measured,na.rm=T))),
                             R2=unlist(lapply(chlB_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                             nsamp=unlist(lapply(chlB_mass_pred_list,nrow)),
                             hellinger=unlist(chlB_mass_dist),
                             KL=unlist(lapply(chlB_mass_KL,function(x) x[[5]])),
                             trait="Chl b")

car_mass_sum_df<-data.frame(functional.group=names(car_mass_pred_list),
                            RMSD=unlist(lapply(car_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                            perRMSD=unlist(lapply(car_mass_pred_list,function(x) percentRMSD(x$measured,x$val_pred,0.025,0.975,na.rm=T))),
                            SD=unlist(lapply(car_mass_pred_list,function(x) sd(x$measured,na.rm=T))),
                            R2=unlist(lapply(car_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                            nsamp=unlist(lapply(car_mass_pred_list,nrow)),
                            hellinger=unlist(car_mass_dist),
                            KL=unlist(lapply(car_mass_KL,function(x) x[[5]])),
                            trait="carotenoids")

summary_df<-do.call(rbind,list(LMA_sum_df,LDMC_sum_df,EWT_sum_df,
                               Nmass_sum_df,Cmass_sum_df,solubles_mass_sum_df,
                               hemicellulose_mass_sum_df,cellulose_mass_sum_df,
                               lignin_mass_sum_df,chlA_mass_sum_df,
                               chlB_mass_sum_df,car_mass_sum_df))

## 2.5% trimmed range across the dataset
summary_df$tr95<-NA
summary_df$tr95[summary_df$trait=="LMA"]<-quantile(meta(spec.traits)$LMA,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$LMA,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="EWT"]<-quantile(meta(spec.traits)$EWT,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$EWT,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="LDMC"]<-quantile(meta(spec.traits)$LDMC,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$LDMC,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="N"]<-quantile(meta(spec.traits)$Nmass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$Nmass,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="C"]<-quantile(meta(spec.traits)$Cmass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$Cmass,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="solubles"]<-quantile(meta(spec.traits)$solubles_mass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$solubles_mass,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="hemicellulose"]<-quantile(meta(spec.traits)$hemicellulose_mass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$hemicellulose_mass,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="cellulose"]<-quantile(meta(spec.traits)$cellulose_mass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$cellulose_mass,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="lignin"]<-quantile(meta(spec.traits)$lignin_mass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$lignin_mass,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="Chl a"]<-quantile(meta(spec.traits)$chlA_mass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$chlA_mass,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="Chl b"]<-quantile(meta(spec.traits)$chlB_mass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$chlB_mass,probs=0.025,na.rm=T)
summary_df$tr95[summary_df$trait=="carotenoids"]<-quantile(meta(spec.traits)$car_mass,probs=0.975,na.rm=T)-quantile(meta(spec.traits)$car_mass,probs=0.025,na.rm=T)
summary_df$perRMSDfull<-summary_df$RMSD/summary_df$tr95*100

## output summary statistics for tables
write.csv(summary_df[summary_df$functional.group=="shrub",],file = "SavedResults/plsr_summ_mt.csv")

pdf("Images/HellingerPerRMSE.pdf",height=6,width=8,onefile=T)
ggplot(data=summary_df,aes(x=hellinger,y=RMSD/tr95*100,
                           color=trait))+
  geom_point(size=2,aes(shape=functional.group))+
  geom_smooth(method="lm",se=F)+theme_bw()+
  theme(text=element_text(size=20))+
  labs(x="Hellinger distance",y=expression("%RMSE"[full]),tag="B")+
  guides(color=guide_legend("Trait"),
         shape="none")+
  scale_color_discrete(labels=expression("LMA","LDMC","EWT","N",
                                   "C","solubles","hemicellulose",
                                   "cellulose","lignin",
                                   "Chl "~italic("a"),
                                   "Chl "~italic("b"),"carotenoids"))
dev.off()

## bringing in randomized validation data
val_models<-readRDS("SavedResults/all_jack_df_list_ref.rds")
summary_df$R2_val<-NA
summary_df$R2_val[summary_df$trait=="LMA"]<-summary(lm(Measured~pred.mean,val_models[["LMA"]]))$r.squared
summary_df$R2_val[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,data=val_models[["EWT"]]))$r.squared
summary_df$R2_val[summary_df$trait=="LDMC"]<-summary(lm(Measured~pred.mean,val_models[["LDMC"]]))$r.squared
summary_df$R2_val[summary_df$trait=="N"]<-summary(lm(Measured~pred.mean,val_models[["N"]]))$r.squared
summary_df$R2_val[summary_df$trait=="C"]<-summary(lm(Measured~pred.mean,val_models[["C"]]))$r.squared
summary_df$R2_val[summary_df$trait=="solubles"]<-summary(lm(Measured~pred.mean,val_models[["sol"]]))$r.squared
summary_df$R2_val[summary_df$trait=="hemicellulose"]<-summary(lm(Measured~pred.mean,val_models[["hemi"]]))$r.squared
summary_df$R2_val[summary_df$trait=="cellulose"]<-summary(lm(Measured~pred.mean,val_models[["cell"]]))$r.squared
summary_df$R2_val[summary_df$trait=="lignin"]<-summary(lm(Measured~pred.mean,val_models[["lign"]]))$r.squared
summary_df$R2_val[summary_df$trait=="Chl a"]<-summary(lm(Measured~pred.mean,val_models[["chlA"]]))$r.squared
summary_df$R2_val[summary_df$trait=="Chl b"]<-summary(lm(Measured~pred.mean,val_models[["chlB"]]))$r.squared
summary_df$R2_val[summary_df$trait=="carotenoids"]<-summary(lm(Measured~pred.mean,val_models[["car"]]))$r.squared

ind_val_list<-readRDS("SavedResults/ind_val_list.rds")
summary_df$R2_Dessain<-NA
summary_df$R2_LOPEX<-NA
summary_df$R2_ANGERS<-NA

LMA_pred_df<-ind_val_list[["LMA"]]
summary_df$R2_Dessain[summary_df$trait=="LMA"]<-summary(lm(Measured~pred.mean,
                                                           LMA_pred_df[which(LMA_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="LMA"]<-summary(lm(Measured~pred.mean,
                                                         LMA_pred_df[which(LMA_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="LMA"]<-summary(lm(Measured~pred.mean,
                                                          LMA_pred_df[which(LMA_pred_df$dataset=="ANGERS"),]))$r.squared

LDMC_pred_df<-ind_val_list[["LDMC"]]
summary_df$R2_Dessain[summary_df$trait=="LDMC"]<-summary(lm(Measured~pred.mean,
                                                           LDMC_pred_df[which(LDMC_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="LDMC"]<-summary(lm(Measured~pred.mean,
                                                         LDMC_pred_df[which(LDMC_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="LDMC"]<-summary(lm(Measured~pred.mean,
                                                          LDMC_pred_df[which(LDMC_pred_df$dataset=="ANGERS"),]))$r.squared

EWT_pred_df<-ind_val_list[["EWT"]]
summary_df$R2_Dessain[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,
                                                           EWT_pred_df[which(EWT_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,
                                                         EWT_pred_df[which(EWT_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,
                                                          EWT_pred_df[which(EWT_pred_df$dataset=="ANGERS"),]))$r.squared

solubles_mass_pred_df<-ind_val_list[["solubles_mass"]]
summary_df$R2_Dessain[summary_df$trait=="solubles"]<-summary(lm(Measured~pred.mean,
                                                           solubles_mass_pred_df[which(solubles_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="solubles"]<-summary(lm(Measured~pred.mean,
                                                         solubles_mass_pred_df[which(solubles_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="solubles"]<-summary(lm(Measured~pred.mean,
                                                          solubles_mass_pred_df[which(solubles_mass_pred_df$dataset=="ANGERS"),]))$r.squared

hemicellulose_mass_pred_df<-ind_val_list[["hemicellulose_mass"]]
summary_df$R2_Dessain[summary_df$trait=="hemicellulose"]<-summary(lm(Measured~pred.mean,
                                                                     hemicellulose_mass_pred_df[which(hemicellulose_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="hemicellulose"]<-summary(lm(Measured~pred.mean,
                                                                   hemicellulose_mass_pred_df[which(hemicellulose_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="hemicellulose"]<-summary(lm(Measured~pred.mean,
                                                                    hemicellulose_mass_pred_df[which(hemicellulose_mass_pred_df$dataset=="ANGERS"),]))$r.squared

cellulose_mass_pred_df<-ind_val_list[["cellulose_mass"]]
summary_df$R2_Dessain[summary_df$trait=="cellulose"]<-summary(lm(Measured~pred.mean,
                                                           cellulose_mass_pred_df[which(cellulose_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="cellulose"]<-summary(lm(Measured~pred.mean,
                                                         cellulose_mass_pred_df[which(cellulose_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="cellulose"]<-summary(lm(Measured~pred.mean,
                                                          cellulose_mass_pred_df[which(cellulose_mass_pred_df$dataset=="ANGERS"),]))$r.squared

lignin_mass_pred_df<-ind_val_list[["lignin_mass"]]
summary_df$R2_Dessain[summary_df$trait=="lignin"]<-summary(lm(Measured~pred.mean,
                                                           lignin_mass_pred_df[which(lignin_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="lignin"]<-summary(lm(Measured~pred.mean,
                                                         lignin_mass_pred_df[which(lignin_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="lignin"]<-summary(lm(Measured~pred.mean,
                                                          lignin_mass_pred_df[which(lignin_mass_pred_df$dataset=="ANGERS"),]))$r.squared

Nmass_pred_df<-ind_val_list[["Nmass"]]
summary_df$R2_Dessain[summary_df$trait=="N"]<-summary(lm(Measured~pred.mean,
                                                           Nmass_pred_df[which(Nmass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="N"]<-summary(lm(Measured~pred.mean,
                                                         Nmass_pred_df[which(Nmass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="N"]<-summary(lm(Measured~pred.mean,
                                                          Nmass_pred_df[which(Nmass_pred_df$dataset=="ANGERS"),]))$r.squared

Cmass_pred_df<-ind_val_list[["Cmass"]]
summary_df$R2_Dessain[summary_df$trait=="C"]<-summary(lm(Measured~pred.mean,
                                                           Cmass_pred_df[which(Cmass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="C"]<-summary(lm(Measured~pred.mean,
                                                         Cmass_pred_df[which(Cmass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="C"]<-summary(lm(Measured~pred.mean,
                                                          Cmass_pred_df[which(Cmass_pred_df$dataset=="ANGERS"),]))$r.squared

chlA_mass_pred_df<-ind_val_list[["chlA_mass"]]
summary_df$R2_Dessain[summary_df$trait=="Chl a"]<-summary(lm(Measured~pred.mean,
                                                           chlA_mass_pred_df[which(chlA_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="Chl a"]<-summary(lm(Measured~pred.mean,
                                                         chlA_mass_pred_df[which(chlA_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="Chl a"]<-summary(lm(Measured~pred.mean,
                                                          chlA_mass_pred_df[which(chlA_mass_pred_df$dataset=="ANGERS"),]))$r.squared

chlB_mass_pred_df<-ind_val_list[["chlB_mass"]]
summary_df$R2_Dessain[summary_df$trait=="Chl b"]<-summary(lm(Measured~pred.mean,
                                                           chlB_mass_pred_df[which(chlB_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="Chl b"]<-summary(lm(Measured~pred.mean,
                                                         chlB_mass_pred_df[which(chlB_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="Chl b"]<-summary(lm(Measured~pred.mean,
                                                          chlB_mass_pred_df[which(chlB_mass_pred_df$dataset=="ANGERS"),]))$r.squared

car_mass_pred_df<-ind_val_list[["car_mass"]]
summary_df$R2_Dessain[summary_df$trait=="carotenoids"]<-summary(lm(Measured~pred.mean,
                                                           car_mass_pred_df[which(car_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="carotenoids"]<-summary(lm(Measured~pred.mean,
                                                         car_mass_pred_df[which(car_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="carotenoids"]<-summary(lm(Measured~pred.mean,
                                                          car_mass_pred_df[which(car_mass_pred_df$dataset=="ANGERS"),]))$r.squared


R2_Asner<-data.frame(trait=c("LMA","LDMC","EWT","N","C","solubles",
                             "hemicellulose","cellulose","lignin",
                             "Chl a","Chl b","carotenoids"),
                     R2=c(0.86,NA,0.88,0.77,0.71,0.63,0.60,0.77,
                          0.62,0.8,0.79,0.76))
summary_df$R2_Asner<-R2_Asner$R2[match(summary_df$trait,R2_Asner$trait)]

ind_val_long<-melt(summary_df[summary_df$functional.group=="broadleaf",
                              c("trait","R2_Dessain","R2_LOPEX","R2_ANGERS")],
                   id.vars="trait")
ind_val_long$variable<-gsub(pattern="R2_",replacement="",x = ind_val_long$variable)

pdf("Images/R2_summary.pdf",height = 6,width=8,onefile=F)
ggplot(data=summary_df,aes(x=trait,y=R2))+
  geom_point(color=colorBlind[7],aes(shape=functional.group),size=2)+
  geom_point(data=ind_val_long,aes(x=trait,y=value,color=variable),size=2)+
  # geom_point(data=summary_df,aes(x=trait,y=R2_Dessain),color="red",shape=1,size=2)+
  # geom_point(data=summary_df,aes(x=trait,y=R2_LOPEX),color="red",shape=2,size=2)+
  # geom_point(data=summary_df,aes(x=trait,y=R2_ANGERS),color="red",shape=3,size=2)+
  geom_point(data=summary_df,aes(x=trait,y=R2_val),color="black",size=3)+
  # geom_point(data=summary_df,aes(x=trait,y=R2_Asner),color="black",size=2,shape=1)+
  theme_bw()+
  scale_x_discrete(labels=expression("LMA","LDMC","EWT","N",
                            "C","sol","hemi","cell",
                            "lign","Chl "~italic("a"),
                            "Chl "~italic("b"),"car"))+
  scale_color_manual(values=colorBlind[c(1,2,4)])+
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0))+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Trait",y=expression(italic("R"^2)),
       shape="Functional group",color="Dataset",tag="A")
dev.off()

write.csv(summary_df,"SavedResults/transfer_summary.csv")
