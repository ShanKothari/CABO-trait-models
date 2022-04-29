setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(pls)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
source("Scripts/VIP.R")
source("Scripts/CABO-trait-models/00 useful_functions.R")

ref.train<-readRDS("ProcessedSpectra/ref_train.rds")
ref.test<-readRDS("ProcessedSpectra/ref_test.rds")

##########################################
## to dos

## try Type II regression?
## add together Chl a and b?

###########################################
## start building first-pass models using
## K-fold cross validation to get the
## optimal number of components

solubles_area_CVmodel<-plsr(meta(ref.train)$solubles_area~as.matrix(ref.train),
                            ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_area_CVmodel <- selectNcomp(solubles_area_CVmodel, method = "onesigma", plot = FALSE)
solubles_area_valid <- which(!is.na(meta(ref.train)$solubles_area))
solubles_area_pred<-data.frame(ID=meta(ref.train)$sample_id[solubles_area_valid],
                               Species=meta(ref.train)$species[solubles_area_valid],
                               Project=meta(ref.train)$project[solubles_area_valid],
                               measured=meta(ref.train)$solubles_area[solubles_area_valid],
                               val_pred=solubles_area_CVmodel$validation$pred[,,ncomp_solubles_area_CVmodel])
ggplot(solubles_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting % solubles from fresh-leaf spectra")

hemicellulose_area_CVmodel<-plsr(meta(ref.train)$hemicellulose_area~as.matrix(ref.train),
                                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_area_CVmodel <- selectNcomp(hemicellulose_area_CVmodel, method = "onesigma", plot = FALSE)
hemicellulose_area_valid <- which(!is.na(meta(ref.train)$hemicellulose_area))
hemicellulose_area_pred<-data.frame(ID=meta(ref.train)$sample_id[hemicellulose_area_valid],
                                    Species=meta(ref.train)$species[hemicellulose_area_valid],
                                    Project=meta(ref.train)$project[hemicellulose_area_valid],
                                    measured=meta(ref.train)$hemicellulose_area[hemicellulose_area_valid],
                                    val_pred=hemicellulose_area_CVmodel$validation$pred[,,ncomp_hemicellulose_area_CVmodel])
ggplot(hemicellulose_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting % hemicellulose from fresh-leaf spectra")

cellulose_area_CVmodel<-plsr(meta(ref.train)$cellulose_area~as.matrix(ref.train),
                             ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_area_CVmodel <- selectNcomp(cellulose_area_CVmodel, method = "onesigma", plot = FALSE)
cellulose_area_valid <- which(!is.na(meta(ref.train)$cellulose_area))
cellulose_area_pred<-data.frame(ID=meta(ref.train)$sample_id[cellulose_area_valid],
                                Species=meta(ref.train)$species[cellulose_area_valid],
                                Project=meta(ref.train)$project[cellulose_area_valid],
                                measured=meta(ref.train)$cellulose_area[cellulose_area_valid],
                                val_pred=cellulose_area_CVmodel$validation$pred[,,ncomp_cellulose_area_CVmodel])
ggplot(cellulose_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting % cellulose from fresh-leaf spectra")

lignin_area_CVmodel<-plsr(meta(ref.train)$lignin_area~as.matrix(ref.train),
                          ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_area_CVmodel <- selectNcomp(lignin_area_CVmodel, method = "onesigma", plot = FALSE)
lignin_area_valid <- which(!is.na(meta(ref.train)$lignin_area))
lignin_area_pred<-data.frame(ID=meta(ref.train)$sample_id[lignin_area_valid],
                             Species=meta(ref.train)$species[lignin_area_valid],
                             Project=meta(ref.train)$project[lignin_area_valid],
                             measured=meta(ref.train)$lignin_area[lignin_area_valid],
                             val_pred=lignin_area_CVmodel$validation$pred[,,ncomp_lignin_area_CVmodel])
ggplot(lignin_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting % lignin from fresh-leaf spectra")

Carea_CVmodel<-plsr(meta(ref.train)$Carea~as.matrix(ref.train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Carea_CVmodel <- selectNcomp(Carea_CVmodel, method = "onesigma", plot = FALSE)
Carea_valid <- which(!is.na(meta(ref.train)$Carea))
Carea_pred<-data.frame(ID=meta(ref.train)$sample_id[Carea_valid],
                       Species=meta(ref.train)$species[Carea_valid],
                       Project=meta(ref.train)$project[Carea_valid],
                       measured=meta(ref.train)$Carea[Carea_valid],
                       val_pred=Carea_CVmodel$validation$pred[,,ncomp_Carea_CVmodel])
ggplot(Carea_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting %C from fresh-leaf spectra")

Narea_CVmodel<-plsr(meta(ref.train)$Narea~as.matrix(ref.train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Narea_CVmodel <- selectNcomp(Narea_CVmodel, method = "onesigma", plot = FALSE)
Narea_valid <- which(!is.na(meta(ref.train)$Narea))
Narea_pred<-data.frame(ID=meta(ref.train)$sample_id[Narea_valid],
                       Species=meta(ref.train)$species[Narea_valid],
                       Project=meta(ref.train)$project[Narea_valid],
                       measured=meta(ref.train)$Narea[Narea_valid],
                       val_pred=Narea_CVmodel$validation$pred[,,ncomp_Narea_CVmodel])
ggplot(Narea_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting %N from fresh-leaf spectra")

chlA_area_CVmodel<-plsr(meta(ref.train)$chlA_area~as.matrix(ref.train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_area_CVmodel <- selectNcomp(chlA_area_CVmodel, method = "onesigma", plot = FALSE)
chlA_area_valid <- which(!is.na(meta(ref.train)$chlA_area))
chlA_area_pred<-data.frame(ID=meta(ref.train)$sample_id[chlA_area_valid],
                           Species=meta(ref.train)$species[chlA_area_valid],
                           Project=meta(ref.train)$project[chlA_area_valid],
                           measured=meta(ref.train)$chlA_area[chlA_area_valid],
                           val_pred=chlA_area_CVmodel$validation$pred[,,ncomp_chlA_area_CVmodel])
ggplot(chlA_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting chlA from fresh-leaf spectra")

chlB_area_CVmodel<-plsr(meta(ref.train)$chlB_area~as.matrix(ref.train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_area_CVmodel <- selectNcomp(chlB_area_CVmodel, method = "onesigma", plot = FALSE)
chlB_area_valid <- which(!is.na(meta(ref.train)$chlB_area))
chlB_area_pred<-data.frame(ID=meta(ref.train)$sample_id[chlB_area_valid],
                           Species=meta(ref.train)$species[chlB_area_valid],
                           Project=meta(ref.train)$project[chlB_area_valid],
                           measured=meta(ref.train)$chlB_area[chlB_area_valid],
                           val_pred=chlB_area_CVmodel$validation$pred[,,ncomp_chlB_area_CVmodel])
ggplot(chlB_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting chlB from fresh-leaf spectra")

car_area_CVmodel<-plsr(meta(ref.train)$car_area~as.matrix(ref.train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_area_CVmodel <- selectNcomp(car_area_CVmodel, method = "onesigma", plot = FALSE)
car_area_valid <- which(!is.na(meta(ref.train)$car_area))
car_area_pred<-data.frame(ID=meta(ref.train)$sample_id[car_area_valid],
                          Species=meta(ref.train)$species[car_area_valid],
                          Project=meta(ref.train)$project[car_area_valid],
                          measured=meta(ref.train)$car_area[car_area_valid],
                          val_pred=car_area_CVmodel$validation$pred[,,ncomp_car_area_CVmodel])
ggplot(car_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting carotenoids from fresh-leaf spectra")

Al_area_CVmodel<-plsr(meta(ref.train)$Al_area~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_area_CVmodel <- selectNcomp(Al_area_CVmodel, method = "onesigma", plot = FALSE)
Al_area_valid <- which(!is.na(meta(ref.train)$Al_area))
Al_area_pred<-data.frame(ID=meta(ref.train)$sample_id[Al_area_valid],
                         Species=meta(ref.train)$species[Al_area_valid],
                         Project=meta(ref.train)$project[Al_area_valid],
                         measured=meta(ref.train)$Al_area[Al_area_valid],
                         val_pred=Al_area_CVmodel$validation$pred[,,ncomp_Al_area_CVmodel])
ggplot(Al_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Al from fresh-leaf spectra")

Ca_area_CVmodel<-plsr(meta(ref.train)$Ca_area~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_area_CVmodel <- selectNcomp(Ca_area_CVmodel, method = "onesigma", plot = FALSE)
Ca_area_valid <- which(!is.na(meta(ref.train)$Ca_area))
Ca_area_pred<-data.frame(ID=meta(ref.train)$sample_id[Ca_area_valid],
                         Species=meta(ref.train)$species[Ca_area_valid],
                         Project=meta(ref.train)$project[Ca_area_valid],
                         measured=meta(ref.train)$Ca_area[Ca_area_valid],
                         val_pred=Ca_area_CVmodel$validation$pred[,,ncomp_Ca_area_CVmodel])
ggplot(Ca_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Ca from fresh-leaf spectra")

Cu_area_CVmodel<-plsr(meta(ref.train)$Cu_area~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_area_CVmodel <- max(c(selectNcomp(Cu_area_CVmodel, method = "onesigma", plot = FALSE),1))
Cu_area_valid <- which(!is.na(meta(ref.train)$Cu_area))
Cu_area_pred<-data.frame(ID=meta(ref.train)$sample_id[Cu_area_valid],
                         Species=meta(ref.train)$species[Cu_area_valid],
                         Project=meta(ref.train)$project[Cu_area_valid],
                         measured=meta(ref.train)$Cu_area[Cu_area_valid],
                         val_pred=Cu_area_CVmodel$validation$pred[,,ncomp_Cu_area_CVmodel])
ggplot(Cu_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Cu from fresh-leaf spectra")

Fe_area_CVmodel<-plsr(meta(ref.train)$Fe_area~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_area_CVmodel <- max(c(selectNcomp(Fe_area_CVmodel, method = "onesigma", plot = FALSE),1))
Fe_area_valid <- which(!is.na(meta(ref.train)$Fe_area))
Fe_area_pred<-data.frame(ID=meta(ref.train)$sample_id[Fe_area_valid],
                         Species=meta(ref.train)$species[Fe_area_valid],
                         Project=meta(ref.train)$project[Fe_area_valid],
                         measured=meta(ref.train)$Fe_area[Fe_area_valid],
                         val_pred=Fe_area_CVmodel$validation$pred[,,ncomp_Fe_area_CVmodel])
ggplot(Fe_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Fe from fresh-leaf spectra")

K_area_CVmodel<-plsr(meta(ref.train)$K_area~as.matrix(ref.train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_area_CVmodel <- selectNcomp(K_area_CVmodel, method = "onesigma", plot = FALSE)
K_area_valid <- which(!is.na(meta(ref.train)$K_area))
K_area_pred<-data.frame(ID=meta(ref.train)$sample_id[K_area_valid],
                        Species=meta(ref.train)$species[K_area_valid],
                        Project=meta(ref.train)$project[K_area_valid],
                        measured=meta(ref.train)$K_area[K_area_valid],
                        val_pred=K_area_CVmodel$validation$pred[,,ncomp_K_area_CVmodel])
ggplot(K_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting K from fresh-leaf spectra")

Mg_area_CVmodel<-plsr(meta(ref.train)$Mg_area~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_area_CVmodel <- selectNcomp(Mg_area_CVmodel, method = "onesigma", plot = FALSE)
Mg_area_valid <- which(!is.na(meta(ref.train)$Mg_area))
Mg_area_pred<-data.frame(ID=meta(ref.train)$sample_id[Mg_area_valid],
                         Species=meta(ref.train)$species[Mg_area_valid],
                         Project=meta(ref.train)$project[Mg_area_valid],
                         measured=meta(ref.train)$Mg_area[Mg_area_valid],
                         val_pred=Mg_area_CVmodel$validation$pred[,,ncomp_Mg_area_CVmodel])
ggplot(Mg_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Mg from fresh-leaf spectra")

Mn_area_CVmodel<-plsr(meta(ref.train)$Mn_area~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_area_CVmodel <- selectNcomp(Mn_area_CVmodel, method = "onesigma", plot = FALSE)
Mn_area_valid <- which(!is.na(meta(ref.train)$Mn_area))
Mn_area_pred<-data.frame(ID=meta(ref.train)$sample_id[Mn_area_valid],
                         Species=meta(ref.train)$species[Mn_area_valid],
                         Project=meta(ref.train)$project[Mn_area_valid],
                         measured=meta(ref.train)$Mn_area[Mn_area_valid],
                         val_pred=Mn_area_CVmodel$validation$pred[,,ncomp_Mn_area_CVmodel])
ggplot(Mn_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Mn from fresh-leaf spectra")

Na_area_CVmodel<-plsr(meta(ref.train)$Na_area~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_area_CVmodel <- selectNcomp(Na_area_CVmodel, method = "onesigma", plot = FALSE)
Na_area_valid <- which(!is.na(meta(ref.train)$Na_area))
Na_area_pred<-data.frame(ID=meta(ref.train)$sample_id[Na_area_valid],
                         Species=meta(ref.train)$species[Na_area_valid],
                         Project=meta(ref.train)$project[Na_area_valid],
                         measured=meta(ref.train)$Na_area[Na_area_valid],
                         val_pred=Na_area_CVmodel$validation$pred[,,ncomp_Na_area_CVmodel])
ggplot(Na_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Na from fresh-leaf spectra")

P_area_CVmodel<-plsr(meta(ref.train)$P_area~as.matrix(ref.train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_area_CVmodel <- selectNcomp(P_area_CVmodel, method = "onesigma", plot = FALSE)
P_area_valid <- which(!is.na(meta(ref.train)$P_area))
P_area_pred<-data.frame(ID=meta(ref.train)$sample_id[P_area_valid],
                        Species=meta(ref.train)$species[P_area_valid],
                        Project=meta(ref.train)$project[P_area_valid],
                        measured=meta(ref.train)$P_area[P_area_valid],
                        val_pred=P_area_CVmodel$validation$pred[,,ncomp_P_area_CVmodel])
ggplot(P_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting P from fresh-leaf spectra")

Zn_area_CVmodel<-plsr(meta(ref.train)$Zn_area~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_area_CVmodel <- selectNcomp(Zn_area_CVmodel, method = "onesigma", plot = FALSE)
Zn_area_valid <- which(!is.na(meta(ref.train)$Zn_area))
Zn_area_pred<-data.frame(ID=meta(ref.train)$sample_id[Zn_area_valid],
                         Species=meta(ref.train)$species[Zn_area_valid],
                         Project=meta(ref.train)$project[Zn_area_valid],
                         measured=meta(ref.train)$Zn_area[Zn_area_valid],
                         val_pred=Zn_area_CVmodel$validation$pred[,,ncomp_Zn_area_CVmodel])
ggplot(Zn_area_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Zn from fresh-leaf spectra")

#######################################################
## VIP from calibration models

VIP1.df<-data.frame(Carea=VIP(Carea_CVmodel)[ncomp_Carea_CVmodel,],
                    Narea=VIP(Narea_CVmodel)[ncomp_Narea_CVmodel,],
                    wavelength=400:2400)

VIP2.df<-data.frame(solubles_area=VIP(solubles_area_CVmodel)[ncomp_solubles_area_CVmodel,],
                    hemicellulose_area=VIP(hemicellulose_area_CVmodel)[ncomp_hemicellulose_area_CVmodel,],
                    cellulose_area=VIP(cellulose_area_CVmodel)[ncomp_cellulose_area_CVmodel,],
                    lignin_area=VIP(lignin_area_CVmodel)[ncomp_lignin_area_CVmodel,],
                    wavelength=400:2400)

VIP3.df<-data.frame(chlA_area=VIP(chlA_area_CVmodel)[ncomp_chlA_area_CVmodel,],
                    chlB_area=VIP(chlB_area_CVmodel)[ncomp_chlB_area_CVmodel,],
                    car_area=VIP(car_area_CVmodel)[ncomp_car_area_CVmodel,],
                    wavelength=400:2400)

VIP4.df<-data.frame(Al_area=VIP(Al_area_CVmodel)[ncomp_Al_area_CVmodel,],
                    Ca_area=VIP(Ca_area_CVmodel)[ncomp_Ca_area_CVmodel,],
                    Cu_area=VIP(Cu_area_CVmodel)[ncomp_Cu_area_CVmodel,],
                    Fe_area=VIP(Fe_area_CVmodel)[ncomp_Fe_area_CVmodel,],
                    K_area=VIP(K_area_CVmodel)[ncomp_K_area_CVmodel,],
                    wavelength=400:2400)

VIP5.df<-data.frame(Mg_area=VIP(Mg_area_CVmodel)[ncomp_Mg_area_CVmodel,],
                    Mn_area=VIP(Mn_area_CVmodel)[ncomp_Mn_area_CVmodel,],
                    Na_area=VIP(Na_area_CVmodel)[ncomp_Na_area_CVmodel,],
                    P_area=VIP(P_area_CVmodel)[ncomp_P_area_CVmodel,],
                    Zn_area=VIP(Zn_area_CVmodel)[ncomp_Zn_area_CVmodel,],
                    wavelength=400:2400)

VIP1.long<-melt(VIP1.df,id.vars = "wavelength")
VIP2.long<-melt(VIP2.df,id.vars = "wavelength")
VIP3.long<-melt(VIP3.df,id.vars = "wavelength")
VIP4.long<-melt(VIP4.df,id.vars = "wavelength")
VIP5.long<-melt(VIP5.df,id.vars = "wavelength")

focal_palette=palette(brewer.pal(8,name="Set2")[c(1,2,4,7,8)])

VIP1.plot<-ggplot(VIP1.long,aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1)+theme_bw()+
  theme(text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength",color="Trait")+
  coord_cartesian(ylim=c(0,4))+
  scale_color_manual(labels=c("C","N","LMA","LDMC","EWT"),
                     values=focal_palette)+
  geom_abline(slope=0,intercept=0.8,linetype="dashed",size=1)

VIP2.plot<-ggplot(VIP2.long,aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1)+theme_bw()+
  theme(text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength",color="Trait")+
  coord_cartesian(ylim=c(0,4))+
  scale_color_manual(labels=c("sol","hemi","cell","lign"),
                     values=focal_palette)+
  geom_abline(slope=0,intercept=0.8,linetype="dashed",size=1)

VIP3.plot<-ggplot(VIP3.long,aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1)+theme_bw()+
  theme(text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength",color="Trait")+
  coord_cartesian(ylim=c(0,4))+
  scale_color_manual(labels=c(expression("Chl"~italic("a")),
                              expression("Chl"~italic("b")),
                              expression("Car")),
                     values=focal_palette)+
  geom_abline(slope=0,intercept=0.8,linetype="dashed",size=1)

VIP4.plot<-ggplot(VIP4.long,aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1)+theme_bw()+
  theme(text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y="VIP",x="Wavelength",color="Trait")+
  coord_cartesian(ylim=c(0,4))+
  scale_color_manual(labels=c("Al","Ca","Cu","Fe","K"),
                     values=focal_palette)+
  geom_abline(slope=0,intercept=0.8,linetype="dashed",size=1)

VIP5.plot<-ggplot(VIP5.long,aes(x=wavelength,y=value,color=variable))+
  geom_line(size=1)+theme_bw()+
  theme(text=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(y="VIP",x="Wavelength",color="Trait")+
  coord_cartesian(ylim=c(0,4))+
  scale_color_manual(labels=c("Mg","Mn","Na","P","Zn"),
                     values=focal_palette)+
  geom_abline(slope=0,intercept=0.8,linetype="dashed",size=1)

pdf("Images/VIP_ref_area.pdf",height=12,width=6,onefile = F)
VIP1.plot / VIP2.plot / VIP3.plot / VIP4.plot / VIP5.plot
dev.off()

#######################################################
## jackknife analyses

solubles_area.jack.coefs<-list()
hemicellulose_area.jack.coefs<-list()
cellulose_area.jack.coefs<-list()
lignin_area.jack.coefs<-list()
chlA_area.jack.coefs<-list()
chlB_area.jack.coefs<-list()
car_area.jack.coefs<-list()
Carea.jack.coefs<-list()
Narea.jack.coefs<-list()
Al_area.jack.coefs<-list()
Ca_area.jack.coefs<-list()
Cu_area.jack.coefs<-list()
Fe_area.jack.coefs<-list()
K_area.jack.coefs<-list()
Mg_area.jack.coefs<-list()
Mn_area.jack.coefs<-list()
Na_area.jack.coefs<-list()
P_area.jack.coefs<-list()
Zn_area.jack.coefs<-list()

solubles_area.jack.stats<-list()
hemicellulose_area.jack.stats<-list()
cellulose_area.jack.stats<-list()
lignin_area.jack.stats<-list()
chlA_area.jack.stats<-list()
chlB_area.jack.stats<-list()
car_area.jack.stats<-list()
Carea.jack.stats<-list()
Narea.jack.stats<-list()
Al_area.jack.stats<-list()
Ca_area.jack.stats<-list()
Cu_area.jack.stats<-list()
Fe_area.jack.stats<-list()
K_area.jack.stats<-list()
Mg_area.jack.stats<-list()
Mn_area.jack.stats<-list()
Na_area.jack.stats<-list()
P_area.jack.stats<-list()
Zn_area.jack.stats<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n.cal.spec<-nrow(ref.train)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-ref.train[train.jack]
  val.jack<-ref.train[test.jack]
  
  solubles_area.jack<-plsr(meta(calib.jack)$solubles_area~as.matrix(calib.jack),
                           ncomp=30,method = "oscorespls",validation="none")
  hemicellulose_area.jack<-plsr(meta(calib.jack)$hemicellulose_area~as.matrix(calib.jack),
                                ncomp=30,method = "oscorespls",validation="none")
  cellulose_area.jack<-plsr(meta(calib.jack)$cellulose_area~as.matrix(calib.jack),
                            ncomp=30,method = "oscorespls",validation="none")
  lignin_area.jack<-plsr(meta(calib.jack)$lignin_area~as.matrix(calib.jack),
                         ncomp=30,method = "oscorespls",validation="none")
  chlA_area.jack<-plsr(meta(calib.jack)$chlA_area~as.matrix(calib.jack),
                       ncomp=30,method = "oscorespls",validation="none")
  chlB_area.jack<-plsr(meta(calib.jack)$chlB_area~as.matrix(calib.jack),
                       ncomp=30,method = "oscorespls",validation="none")
  car_area.jack<-plsr(meta(calib.jack)$car_area~as.matrix(calib.jack),
                      ncomp=30,method = "oscorespls",validation="none")
  Carea.jack<-plsr(meta(calib.jack)$Carea~as.matrix(calib.jack),
                   ncomp=30,method = "oscorespls",validation="none")
  Narea.jack<-plsr(meta(calib.jack)$Narea~as.matrix(calib.jack),
                   ncomp=30,method = "oscorespls",validation="none")
  Al_area.jack<-plsr(meta(calib.jack)$Al_area~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Ca_area.jack<-plsr(meta(calib.jack)$Ca_area~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Cu_area.jack<-plsr(meta(calib.jack)$Cu_area~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Fe_area.jack<-plsr(meta(calib.jack)$Fe_area~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  K_area.jack<-plsr(meta(calib.jack)$K_area~as.matrix(calib.jack),
                    ncomp=30,method = "oscorespls",validation="none")
  Mg_area.jack<-plsr(meta(calib.jack)$Mg_area~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Mn_area.jack<-plsr(meta(calib.jack)$Mn_area~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Na_area.jack<-plsr(meta(calib.jack)$Na_area~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  P_area.jack<-plsr(meta(calib.jack)$P_area~as.matrix(calib.jack),
                    ncomp=30,method = "oscorespls",validation="none")
  Zn_area.jack<-plsr(meta(calib.jack)$Zn_area~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  
  solubles_area.jack.val.pred<-as.vector(predict(solubles_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_solubles_area_CVmodel)[,,1])
  solubles_area.jack.val.fit<-lm(solubles_area.jack.val.pred~meta(val.jack)$solubles_area)
  solubles_area.jack.stats[[i]]<-c(R2=summary(solubles_area.jack.val.fit)$r.squared,
                                   RMSE=RMSD(meta(val.jack)$solubles_area,solubles_area.jack.val.pred),
                                   percentRMSE=percentRMSD(meta(val.jack)$solubles_area,solubles_area.jack.val.pred,0.025,0.975),
                                   bias=mean(solubles_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$solubles_area,na.rm=T))
  
  hemicellulose_area.jack.val.pred<-as.vector(predict(hemicellulose_area.jack,newdata=as.matrix(val.jack),
                                                      ncomp=ncomp_hemicellulose_area_CVmodel)[,,1])
  hemicellulose_area.jack.val.fit<-lm(hemicellulose_area.jack.val.pred~meta(val.jack)$hemicellulose_area)
  hemicellulose_area.jack.stats[[i]]<-c(R2=summary(hemicellulose_area.jack.val.fit)$r.squared,
                                        RMSE=RMSD(meta(val.jack)$hemicellulose_area,hemicellulose_area.jack.val.pred),
                                        percentRMSE=percentRMSD(meta(val.jack)$hemicellulose_area,hemicellulose_area.jack.val.pred,0.025,0.975),
                                        bias=mean(hemicellulose_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$hemicellulose_area,na.rm=T))
  
  cellulose_area.jack.val.pred<-as.vector(predict(cellulose_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_cellulose_area_CVmodel)[,,1])
  cellulose_area.jack.val.fit<-lm(cellulose_area.jack.val.pred~meta(val.jack)$cellulose_area)
  cellulose_area.jack.stats[[i]]<-c(R2=summary(cellulose_area.jack.val.fit)$r.squared,
                                    RMSE=RMSD(meta(val.jack)$cellulose_area,cellulose_area.jack.val.pred),
                                    percentRMSE=percentRMSD(meta(val.jack)$cellulose_area,cellulose_area.jack.val.pred,0.025,0.975),
                                    bias=mean(cellulose_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$cellulose_area,na.rm=T))
  
  lignin_area.jack.val.pred<-as.vector(predict(lignin_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_lignin_area_CVmodel)[,,1])
  lignin_area.jack.val.fit<-lm(lignin_area.jack.val.pred~meta(val.jack)$lignin_area)
  lignin_area.jack.stats[[i]]<-c(R2=summary(lignin_area.jack.val.fit)$r.squared,
                                 RMSE=RMSD(meta(val.jack)$lignin_area,lignin_area.jack.val.pred),
                                 percentRMSE=percentRMSD(meta(val.jack)$lignin_area,lignin_area.jack.val.pred,0.025,0.975),
                                 bias=mean(lignin_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$lignin_area,na.rm=T))
  
  chlA_area.jack.val.pred<-as.vector(predict(chlA_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlA_area_CVmodel)[,,1])
  chlA_area.jack.val.fit<-lm(chlA_area.jack.val.pred~meta(val.jack)$chlA_area)
  chlA_area.jack.stats[[i]]<-c(R2=summary(chlA_area.jack.val.fit)$r.squared,
                               RMSE=RMSD(meta(val.jack)$chlA_area,chlA_area.jack.val.pred),
                               percentRMSE=percentRMSD(meta(val.jack)$chlA_area,chlA_area.jack.val.pred,0.025,0.975),
                               bias=mean(chlA_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlA_area,na.rm=T))
  
  chlB_area.jack.val.pred<-as.vector(predict(chlB_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlB_area_CVmodel)[,,1])
  chlB_area.jack.val.fit<-lm(chlB_area.jack.val.pred~meta(val.jack)$chlB_area)
  chlB_area.jack.stats[[i]]<-c(R2=summary(chlB_area.jack.val.fit)$r.squared,
                               RMSE=RMSD(meta(val.jack)$chlB_area,chlB_area.jack.val.pred),
                               percentRMSE=percentRMSD(meta(val.jack)$chlB_area,chlB_area.jack.val.pred,0.025,0.975),
                               bias=mean(chlB_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlB_area,na.rm=T))
  
  car_area.jack.val.pred<-as.vector(predict(car_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_car_area_CVmodel)[,,1])
  car_area.jack.val.fit<-lm(car_area.jack.val.pred~meta(val.jack)$car_area)
  car_area.jack.stats[[i]]<-c(R2=summary(car_area.jack.val.fit)$r.squared,
                              RMSE=RMSD(meta(val.jack)$car_area,car_area.jack.val.pred),
                              percentRMSE=percentRMSD(meta(val.jack)$car_area,car_area.jack.val.pred,0.025,0.975),
                              bias=mean(car_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$car_area,na.rm=T))
  
  Carea.jack.val.pred<-as.vector(predict(Carea.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Carea_CVmodel)[,,1])
  Carea.jack.val.fit<-lm(Carea.jack.val.pred~meta(val.jack)$Carea)
  Carea.jack.stats[[i]]<-c(R2=summary(Carea.jack.val.fit)$r.squared,
                           RMSE=RMSD(meta(val.jack)$Carea,Carea.jack.val.pred),
                           percentRMSE=percentRMSD(meta(val.jack)$Carea,Carea.jack.val.pred,0.025,0.975),
                           bias=mean(Carea.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Carea,na.rm=T))
  
  Narea.jack.val.pred<-as.vector(predict(Narea.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Narea_CVmodel)[,,1])
  Narea.jack.val.fit<-lm(Narea.jack.val.pred~meta(val.jack)$Narea)
  Narea.jack.stats[[i]]<-c(R2=summary(Narea.jack.val.fit)$r.squared,
                           RMSE=RMSD(meta(val.jack)$Narea,Narea.jack.val.pred),
                           percentRMSE=percentRMSD(meta(val.jack)$Narea,Narea.jack.val.pred,0.025,0.975),
                           bias=mean(Narea.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Narea,na.rm=T))
  
  Al_area.jack.val.pred<-as.vector(predict(Al_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Al_area_CVmodel)[,,1])
  Al_area.jack.val.fit<-lm(Al_area.jack.val.pred~meta(val.jack)$Al_area)
  Al_area.jack.stats[[i]]<-c(R2=summary(Al_area.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Al_area,Al_area.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Al_area,Al_area.jack.val.pred,0.025,0.975),
                             bias=mean(Al_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Al_area,na.rm=T))
  
  Ca_area.jack.val.pred<-as.vector(predict(Ca_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Ca_area_CVmodel)[,,1])
  Ca_area.jack.val.fit<-lm(Ca_area.jack.val.pred~meta(val.jack)$Ca_area)
  Ca_area.jack.stats[[i]]<-c(R2=summary(Ca_area.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Ca_area,Ca_area.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Ca_area,Ca_area.jack.val.pred,0.025,0.975),
                             bias=mean(Ca_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Ca_area,na.rm=T))
  
  Cu_area.jack.val.pred<-as.vector(predict(Cu_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Cu_area_CVmodel)[,,1])
  Cu_area.jack.val.fit<-lm(Cu_area.jack.val.pred~meta(val.jack)$Cu_area)
  Cu_area.jack.stats[[i]]<-c(R2=summary(Cu_area.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Cu_area,Cu_area.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Cu_area,Cu_area.jack.val.pred,0.025,0.975),
                             bias=mean(Cu_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Cu_area,na.rm=T))
  
  Fe_area.jack.val.pred<-as.vector(predict(Fe_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Fe_area_CVmodel)[,,1])
  Fe_area.jack.val.fit<-lm(Fe_area.jack.val.pred~meta(val.jack)$Fe_area)
  Fe_area.jack.stats[[i]]<-c(R2=summary(Fe_area.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Fe_area,Fe_area.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Fe_area,Fe_area.jack.val.pred,0.025,0.975),
                             bias=mean(Fe_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Fe_area,na.rm=T))
  
  K_area.jack.val.pred<-as.vector(predict(K_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_K_area_CVmodel)[,,1])
  K_area.jack.val.fit<-lm(K_area.jack.val.pred~meta(val.jack)$K_area)
  K_area.jack.stats[[i]]<-c(R2=summary(K_area.jack.val.fit)$r.squared,
                            RMSE=RMSD(meta(val.jack)$K_area,K_area.jack.val.pred),
                            percentRMSE=percentRMSD(meta(val.jack)$K_area,K_area.jack.val.pred,0.025,0.975),
                            bias=mean(K_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$K_area,na.rm=T))
  
  Mg_area.jack.val.pred<-as.vector(predict(Mg_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Mg_area_CVmodel)[,,1])
  Mg_area.jack.val.fit<-lm(Mg_area.jack.val.pred~meta(val.jack)$Mg_area)
  Mg_area.jack.stats[[i]]<-c(R2=summary(Mg_area.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Mg_area,Mg_area.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Mg_area,Mg_area.jack.val.pred,0.025,0.975),
                             bias=mean(Mg_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Mg_area,na.rm=T))
  
  Mn_area.jack.val.pred<-as.vector(predict(Mn_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Mn_area_CVmodel)[,,1])
  Mn_area.jack.val.fit<-lm(Mn_area.jack.val.pred~meta(val.jack)$Mn_area)
  Mn_area.jack.stats[[i]]<-c(R2=summary(Mn_area.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Mn_area,Mn_area.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Mn_area,Mn_area.jack.val.pred,0.025,0.975),
                             bias=mean(Mn_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Mn_area,na.rm=T))
  
  Na_area.jack.val.pred<-as.vector(predict(Na_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Na_area_CVmodel)[,,1])
  Na_area.jack.val.fit<-lm(Na_area.jack.val.pred~meta(val.jack)$Na_area)
  Na_area.jack.stats[[i]]<-c(R2=summary(Na_area.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Na_area,Na_area.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Na_area,Na_area.jack.val.pred,0.025,0.975),
                             bias=mean(Na_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Na_area,na.rm=T))
  
  P_area.jack.val.pred<-as.vector(predict(P_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_P_area_CVmodel)[,,1])
  P_area.jack.val.fit<-lm(P_area.jack.val.pred~meta(val.jack)$P_area)
  P_area.jack.stats[[i]]<-c(R2=summary(P_area.jack.val.fit)$r.squared,
                            RMSE=RMSD(meta(val.jack)$P_area,P_area.jack.val.pred),
                            percentRMSE=percentRMSD(meta(val.jack)$P_area,P_area.jack.val.pred,0.025,0.975),
                            bias=mean(P_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$P_area,na.rm=T))
  
  Zn_area.jack.val.pred<-as.vector(predict(Zn_area.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Zn_area_CVmodel)[,,1])
  Zn_area.jack.val.fit<-lm(Zn_area.jack.val.pred~meta(val.jack)$Zn_area)
  Zn_area.jack.stats[[i]]<-c(R2=summary(Zn_area.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Zn_area,Zn_area.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Zn_area,Zn_area.jack.val.pred,0.025,0.975),
                             bias=mean(Zn_area.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Zn_area,na.rm=T))
  
  solubles_area.jack.coefs[[i]]<-as.vector(coef(solubles_area.jack,ncomp=ncomp_solubles_area_CVmodel,intercept=TRUE))
  hemicellulose_area.jack.coefs[[i]]<-as.vector(coef(hemicellulose_area.jack,ncomp=ncomp_hemicellulose_area_CVmodel,intercept=TRUE))
  cellulose_area.jack.coefs[[i]]<-as.vector(coef(cellulose_area.jack,ncomp=ncomp_cellulose_area_CVmodel,intercept=TRUE))
  lignin_area.jack.coefs[[i]]<-as.vector(coef(lignin_area.jack,ncomp=ncomp_lignin_area_CVmodel,intercept=TRUE))
  chlA_area.jack.coefs[[i]]<-as.vector(coef(chlA_area.jack,ncomp=ncomp_chlA_area_CVmodel,intercept=TRUE))
  chlB_area.jack.coefs[[i]]<-as.vector(coef(chlB_area.jack,ncomp=ncomp_chlB_area_CVmodel,intercept=TRUE))
  car_area.jack.coefs[[i]]<-as.vector(coef(car_area.jack,ncomp=ncomp_car_area_CVmodel,intercept=TRUE))
  Carea.jack.coefs[[i]]<-as.vector(coef(Carea.jack,ncomp=ncomp_Carea_CVmodel,intercept=TRUE))
  Narea.jack.coefs[[i]]<-as.vector(coef(Narea.jack,ncomp=ncomp_Narea_CVmodel,intercept=TRUE))
  Al_area.jack.coefs[[i]]<-as.vector(coef(Al_area.jack,ncomp=ncomp_Al_area_CVmodel,intercept=TRUE))
  Ca_area.jack.coefs[[i]]<-as.vector(coef(Ca_area.jack,ncomp=ncomp_Ca_area_CVmodel,intercept=TRUE))
  Cu_area.jack.coefs[[i]]<-as.vector(coef(Cu_area.jack,ncomp=ncomp_Cu_area_CVmodel,intercept=TRUE))
  Fe_area.jack.coefs[[i]]<-as.vector(coef(Fe_area.jack,ncomp=ncomp_Fe_area_CVmodel,intercept=TRUE))
  K_area.jack.coefs[[i]]<-as.vector(coef(K_area.jack,ncomp=ncomp_K_area_CVmodel,intercept=TRUE))
  Mg_area.jack.coefs[[i]]<-as.vector(coef(Mg_area.jack,ncomp=ncomp_Mg_area_CVmodel,intercept=TRUE))
  Mn_area.jack.coefs[[i]]<-as.vector(coef(Mn_area.jack,ncomp=ncomp_Mn_area_CVmodel,intercept=TRUE))
  Na_area.jack.coefs[[i]]<-as.vector(coef(Na_area.jack,ncomp=ncomp_Na_area_CVmodel,intercept=TRUE))
  P_area.jack.coefs[[i]]<-as.vector(coef(P_area.jack,ncomp=ncomp_P_area_CVmodel,intercept=TRUE))
  Zn_area.jack.coefs[[i]]<-as.vector(coef(Zn_area.jack,ncomp=ncomp_Zn_area_CVmodel,intercept=TRUE))
  
}

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

solubles_area.jack.pred<-apply.coefs(solubles_area.jack.coefs,as.matrix(ref.test))
solubles_area.jack.stat<-t(apply(solubles_area.jack.pred,1,
                                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_area.jack.df<-data.frame(pred.mean=solubles_area.jack.stat[,1],
                                  pred.low=solubles_area.jack.stat[,2],
                                  pred.high=solubles_area.jack.stat[,3],
                                  Measured=meta(ref.test)$solubles_area,
                                  ncomp=ncomp_solubles_area_CVmodel,
                                  Species=meta(ref.test)$species,
                                  Project=meta(ref.test)$project,
                                  functional.group=meta(ref.test)$functional.group,
                                  ID=meta(ref.test)$sample_id)

solubles_area.val.plot<-ggplot(solubles_area.jack.df,
                               aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-3,25),ylim=c(-3,25))+
  theme(text = element_text(size=20),
        legend.position = c(0, 0.25))+
  labs(y=expression("Measured solubles (mg cm"^-2*")"),
       x=expression("Measured solubles (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

hemicellulose_area.jack.pred<-apply.coefs(hemicellulose_area.jack.coefs,as.matrix(ref.test))
hemicellulose_area.jack.stat<-t(apply(hemicellulose_area.jack.pred,1,
                                      function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_area.jack.df<-data.frame(pred.mean=hemicellulose_area.jack.stat[,1],
                                       pred.low=hemicellulose_area.jack.stat[,2],
                                       pred.high=hemicellulose_area.jack.stat[,3],
                                       Measured=meta(ref.test)$hemicellulose_area,
                                       ncomp=ncomp_hemicellulose_area_CVmodel,
                                       Species=meta(ref.test)$species,
                                       Project=meta(ref.test)$project,
                                       functional.group=meta(ref.test)$functional.group,
                                       ID=meta(ref.test)$sample_id)

hemicellulose_area.val.plot<-ggplot(hemicellulose_area.jack.df,
                                    aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-1,5.5),ylim=c(-1,5.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured hemicellulose (mg cm"^-2*")"),
       x=expression("Measured hemicellulose (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

cellulose_area.jack.pred<-apply.coefs(cellulose_area.jack.coefs,as.matrix(ref.test))
cellulose_area.jack.stat<-t(apply(cellulose_area.jack.pred,1,
                                  function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_area.jack.df<-data.frame(pred.mean=cellulose_area.jack.stat[,1],
                                   pred.low=cellulose_area.jack.stat[,2],
                                   pred.high=cellulose_area.jack.stat[,3],
                                   Measured=meta(ref.test)$cellulose_area,
                                   ncomp=ncomp_cellulose_area_CVmodel,
                                   Species=meta(ref.test)$species,
                                   Project=meta(ref.test)$project,
                                   functional.group=meta(ref.test)$functional.group,
                                   ID=meta(ref.test)$sample_id)

cellulose_area.val.plot<-ggplot(cellulose_area.jack.df,
                                aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-1,7),ylim=c(-1,7))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured cellulose (mg cm"^-2*")"),
       x=expression("Measured cellulose (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

lignin_area.jack.pred<-apply.coefs(lignin_area.jack.coefs,as.matrix(ref.test))
lignin_area.jack.stat<-t(apply(lignin_area.jack.pred,1,
                               function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_area.jack.df<-data.frame(pred.mean=lignin_area.jack.stat[,1],
                                pred.low=lignin_area.jack.stat[,2],
                                pred.high=lignin_area.jack.stat[,3],
                                Measured=meta(ref.test)$lignin_area,
                                ncomp=ncomp_lignin_area_CVmodel,
                                Species=meta(ref.test)$species,
                                Project=meta(ref.test)$project,
                                functional.group=meta(ref.test)$functional.group,
                                ID=meta(ref.test)$sample_id)

lignin_area.val.plot<-ggplot(lignin_area.jack.df,
                             aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.5,4.5),ylim=c(-0.5,4.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured lignin (mg cm"^-2*")"),
       x=expression("Measured lignin (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_area.jack.pred<-apply.coefs(chlA_area.jack.coefs,as.matrix(ref.test))
chlA_area.jack.stat<-t(apply(chlA_area.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_area.jack.df<-data.frame(pred.mean=chlA_area.jack.stat[,1],
                              pred.low=chlA_area.jack.stat[,2],
                              pred.high=chlA_area.jack.stat[,3],
                              Measured=meta(ref.test)$chlA_area,
                              ncomp=ncomp_chlA_area_CVmodel,
                              Species=meta(ref.test)$species,
                              Project=meta(ref.test)$project,
                              functional.group=meta(ref.test)$functional.group,
                              ID=meta(ref.test)$sample_id)

chlA_area.val.plot<-ggplot(chlA_area.jack.df,
                           aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.25),ylim=c(0,0.25))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg cm"^-2*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlB_area.jack.pred<-apply.coefs(chlB_area.jack.coefs,as.matrix(ref.test))
chlB_area.jack.stat<-t(apply(chlB_area.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_area.jack.df<-data.frame(pred.mean=chlB_area.jack.stat[,1],
                              pred.low=chlB_area.jack.stat[,2],
                              pred.high=chlB_area.jack.stat[,3],
                              Measured=meta(ref.test)$chlB_area,
                              ncomp=ncomp_chlB_area_CVmodel,
                              Species=meta(ref.test)$species,
                              Project=meta(ref.test)$project,
                              functional.group=meta(ref.test)$functional.group,
                              ID=meta(ref.test)$sample_id)

chlB_area.val.plot<-ggplot(chlB_area.jack.df,
                           aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.075),ylim=c(0,0.075))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg cm"^-2*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

car_area.jack.pred<-apply.coefs(car_area.jack.coefs,as.matrix(ref.test))
car_area.jack.stat<-t(apply(car_area.jack.pred,1,
                            function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_area.jack.df<-data.frame(pred.mean=car_area.jack.stat[,1],
                             pred.low=car_area.jack.stat[,2],
                             pred.high=car_area.jack.stat[,3],
                             Measured=meta(ref.test)$car_area,
                             ncomp=ncomp_car_area_CVmodel,
                             Species=meta(ref.test)$species,
                             Project=meta(ref.test)$project,
                             functional.group=meta(ref.test)$functional.group,
                             ID=meta(ref.test)$sample_id)

car_area.val.plot<-ggplot(car_area.jack.df,
                          aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.05),ylim=c(0,0.05))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured carotenoids (mg cm"^-2*")"),
       x=expression("Predicted carotenoids (mg cm"^-2*")"))+
  scale_color_manual(values=colorBlind)

Carea.jack.pred<-apply.coefs(Carea.jack.coefs,as.matrix(ref.test))
Carea.jack.stat<-t(apply(Carea.jack.pred,1,
                         function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Carea.jack.df<-data.frame(pred.mean=Carea.jack.stat[,1],
                          pred.low=Carea.jack.stat[,2],
                          pred.high=Carea.jack.stat[,3],
                          Measured=meta(ref.test)$Carea,
                          ncomp=ncomp_Carea_CVmodel,
                          Species=meta(ref.test)$species,
                          Project=meta(ref.test)$project,
                          functional.group=meta(ref.test)$functional.group,
                          ID=meta(ref.test)$sample_id)

Carea.val.plot<-ggplot(Carea.jack.df,
                       aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,19),ylim=c(0,19))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured C (mg cm"^-2*")"),
       x=expression("Predicted C (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Narea.jack.pred<-apply.coefs(Narea.jack.coefs,as.matrix(ref.test))
Narea.jack.stat<-t(apply(Narea.jack.pred,1,
                         function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Narea.jack.df<-data.frame(pred.mean=Narea.jack.stat[,1],
                          pred.low=Narea.jack.stat[,2],
                          pred.high=Narea.jack.stat[,3],
                          Measured=meta(ref.test)$Narea,
                          ncomp=ncomp_Narea_CVmodel,
                          Species=meta(ref.test)$species,
                          Project=meta(ref.test)$project,
                          functional.group=meta(ref.test)$functional.group,
                          ID=meta(ref.test)$sample_id)

Narea.val.plot<-ggplot(Narea.jack.df,
                       aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.42),ylim=c(0,0.42))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured N (mg cm"^-2*")"),
       x=expression("Predicted N (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Al_area.jack.pred<-apply.coefs(Al_area.jack.coefs,as.matrix(ref.test))
Al_area.jack.stat<-t(apply(Al_area.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_area.jack.df<-data.frame(pred.mean=Al_area.jack.stat[,1],
                            pred.low=Al_area.jack.stat[,2],
                            pred.high=Al_area.jack.stat[,3],
                            Measured=meta(ref.test)$Al_area,
                            ncomp=ncomp_Al_area_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Al_area.val.plot<-ggplot(Al_area.jack.df,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,1.8),ylim=c(-0.1,1.8))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Al (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Al (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Ca_area.jack.pred<-apply.coefs(Ca_area.jack.coefs,as.matrix(ref.test))
Ca_area.jack.stat<-t(apply(Ca_area.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_area.jack.df<-data.frame(pred.mean=Ca_area.jack.stat[,1],
                            pred.low=Ca_area.jack.stat[,2],
                            pred.high=Ca_area.jack.stat[,3],
                            Measured=meta(ref.test)$Ca_area,
                            ncomp=ncomp_Ca_area_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Ca_area.val.plot<-ggplot(Ca_area.jack.df,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-50,300),ylim=c(-50,300))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Ca (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Ca (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cu_area.jack.pred<-apply.coefs(Cu_area.jack.coefs,as.matrix(ref.test))
Cu_area.jack.stat<-t(apply(Cu_area.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_area.jack.df<-data.frame(pred.mean=Cu_area.jack.stat[,1],
                            pred.low=Cu_area.jack.stat[,2],
                            pred.high=Cu_area.jack.stat[,3],
                            Measured=meta(ref.test)$Cu_area,
                            ncomp=ncomp_Cu_area_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Cu_area.val.plot<-ggplot(Cu_area.jack.df,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.025,0.25),ylim=c(-0.025,0.25))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Cu (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Cu (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Fe_area.jack.pred<-apply.coefs(Fe_area.jack.coefs,as.matrix(ref.test))
Fe_area.jack.stat<-t(apply(Fe_area.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_area.jack.df<-data.frame(pred.mean=Fe_area.jack.stat[,1],
                            pred.low=Fe_area.jack.stat[,2],
                            pred.high=Fe_area.jack.stat[,3],
                            Measured=meta(ref.test)$Fe_area,
                            ncomp=ncomp_Fe_area_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Fe_area.val.plot<-ggplot(Fe_area.jack.df,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2),ylim=c(0,2))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Fe (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Fe (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

K_area.jack.pred<-apply.coefs(K_area.jack.coefs,as.matrix(ref.test))
K_area.jack.stat<-t(apply(K_area.jack.pred,1,
                          function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_area.jack.df<-data.frame(pred.mean=K_area.jack.stat[,1],
                           pred.low=K_area.jack.stat[,2],
                           pred.high=K_area.jack.stat[,3],
                           Measured=meta(ref.test)$K_area,
                           ncomp=ncomp_K_area_CVmodel,
                           Species=meta(ref.test)$species,
                           Project=meta(ref.test)$project,
                           functional.group=meta(ref.test)$functional.group,
                           ID=meta(ref.test)$sample_id)

K_area.val.plot<-ggplot(K_area.jack.df,
                        aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,150),ylim=c(0,150))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured K (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted K (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mg_area.jack.pred<-apply.coefs(Mg_area.jack.coefs,as.matrix(ref.test))
Mg_area.jack.stat<-t(apply(Mg_area.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_area.jack.df<-data.frame(pred.mean=Mg_area.jack.stat[,1],
                            pred.low=Mg_area.jack.stat[,2],
                            pred.high=Mg_area.jack.stat[,3],
                            Measured=meta(ref.test)$Mg_area,
                            ncomp=ncomp_Mg_area_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Mg_area.val.plot<-ggplot(Mg_area.jack.df,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,60),ylim=c(0,60))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Mg (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Mg (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mn_area.jack.pred<-apply.coefs(Mn_area.jack.coefs,as.matrix(ref.test))
Mn_area.jack.stat<-t(apply(Mn_area.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_area.jack.df<-data.frame(pred.mean=Mn_area.jack.stat[,1],
                            pred.low=Mn_area.jack.stat[,2],
                            pred.high=Mn_area.jack.stat[,3],
                            Measured=meta(ref.test)$Mn_area,
                            ncomp=ncomp_Mn_area_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Mn_area.val.plot<-ggplot(Mn_area.jack.df,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-5,12),ylim=c(-5,12))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Mn (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Mn (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Na_area.jack.pred<-apply.coefs(Na_area.jack.coefs,as.matrix(ref.test))
Na_area.jack.stat<-t(apply(Na_area.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_area.jack.df<-data.frame(pred.mean=Na_area.jack.stat[,1],
                            pred.low=Na_area.jack.stat[,2],
                            pred.high=Na_area.jack.stat[,3],
                            Measured=meta(ref.test)$Na_area,
                            ncomp=ncomp_Na_area_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Na_area.val.plot<-ggplot(Na_area.jack.df,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-10,60),ylim=c(-10,60))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Na (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Na (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

P_area.jack.pred<-apply.coefs(P_area.jack.coefs,as.matrix(ref.test))
P_area.jack.stat<-t(apply(P_area.jack.pred,1,
                          function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_area.jack.df<-data.frame(pred.mean=P_area.jack.stat[,1],
                           pred.low=P_area.jack.stat[,2],
                           pred.high=P_area.jack.stat[,3],
                           Measured=meta(ref.test)$P_area,
                           ncomp=ncomp_P_area_CVmodel,
                           Species=meta(ref.test)$species,
                           Project=meta(ref.test)$project,
                           functional.group=meta(ref.test)$functional.group,
                           ID=meta(ref.test)$sample_id)

P_area.val.plot<-ggplot(P_area.jack.df,
                        aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,45),ylim=c(0,45))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured P (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted P (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Zn_area.jack.pred<-apply.coefs(Zn_area.jack.coefs,as.matrix(ref.test))
Zn_area.jack.stat<-t(apply(Zn_area.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_area.jack.df<-data.frame(pred.mean=Zn_area.jack.stat[,1],
                            pred.low=Zn_area.jack.stat[,2],
                            pred.high=Zn_area.jack.stat[,3],
                            Measured=meta(ref.test)$Zn_area,
                            ncomp=ncomp_Zn_area_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Zn_area.val.plot<-ggplot(Zn_area.jack.df,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-1,2.5),ylim=c(-1,2.5))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Zn (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Zn (",mu,"g cm"^-2*")")))+
  scale_color_manual(values=colorBlind)

all.jack.coef.list<-list(solubles_area=solubles_area.jack.coefs,
                         hemicellulose_area=hemicellulose_area.jack.coefs,
                         cellulose_area=cellulose_area.jack.coefs,
                         lignin_area=lignin_area.jack.coefs,
                         chlA_area=chlA_area.jack.coefs,
                         chlB_area=chlB_area.jack.coefs,
                         car_area=car_area.jack.coefs,
                         Narea=Narea.jack.coefs,
                         Carea=Carea.jack.coefs,
                         Al_area=Al_area.jack.coefs,
                         Ca_area=Ca_area.jack.coefs,
                         Cu_area=Cu_area.jack.coefs,
                         Fe_area=Fe_area.jack.coefs,
                         K_area=K_area.jack.coefs,
                         Mg_area=Mg_area.jack.coefs,
                         Mn_area=Mn_area.jack.coefs,
                         Na_area=Na_area.jack.coefs,
                         P_area=P_area.jack.coefs,
                         Zn_area=Zn_area.jack.coefs)
saveRDS(all.jack.coef.list,"SavedResults/all_jack_coefs_list_ref_area.rds")

all.jack.df.list<-list(solubles_area=solubles_area.jack.df,
                       hemicellulose_area=hemicellulose_area.jack.df,
                       cellulose_area=cellulose_area.jack.df,
                       lignin_area=lignin_area.jack.df,
                       chlA_area=chlA_area.jack.df,
                       chlB_area=chlB_area.jack.df,
                       car_area=car_area.jack.df,
                       Narea=Narea.jack.df,
                       Carea=Carea.jack.df,
                       Al_area=Al_area.jack.df,
                       Ca_area=Ca_area.jack.df,
                       Cu_area=Cu_area.jack.df,
                       Fe_area=Fe_area.jack.df,
                       K_area=K_area.jack.df,
                       Mg_area=Mg_area.jack.df,
                       Mn_area=Mn_area.jack.df,
                       Na_area=Na_area.jack.df,
                       P_area=P_area.jack.df,
                       Zn_area=Zn_area.jack.df)
saveRDS(all.jack.df.list,"SavedResults/all_jack_df_list_ref_area.rds")

all.jack.stats.list<-list(sol=solubles_area.jack.stats,
                          hemi=hemicellulose_area.jack.stats,
                          cell=cellulose_area.jack.stats,
                          lign=lignin_area.jack.stats,
                          chlA=chlA_area.jack.stats,
                          chlB=chlB_area.jack.stats,
                          car=car_area.jack.stats,
                          N=Narea.jack.stats,
                          C=Carea.jack.stats,
                          Al=Al_area.jack.stats,
                          Ca=Ca_area.jack.stats,
                          Cu=Cu_area.jack.stats,
                          Fe=Fe_area.jack.stats,
                          K=K_area.jack.stats,
                          Mg=Mg_area.jack.stats,
                          Mn=Mn_area.jack.stats,
                          Na=Na_area.jack.stats,
                          P=P_area.jack.stats,
                          Zn=Zn_area.jack.stats)
saveRDS(all.jack.stats.list,"SavedResults/all_jack_stats_list_ref_area.rds")

pdf("Images/val_plots_ref_area1.pdf",width = 16,height = 15)
(Narea.val.plot+Carea.val.plot+solubles_area.val.plot)/
  (hemicellulose_area.val.plot+cellulose_area.val.plot+lignin_area.val.plot)/
  (chlA_area.val.plot+chlB_area.val.plot+car_area.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Images/val_plots_ref_area2.pdf",width = 16,height = 19)
(Al_area.val.plot + Ca_area.val.plot + Cu_area.val.plot) / 
  (Fe_area.val.plot + K_area.val.plot + Mg_area.val.plot) / 
  (Mn_area.val.plot + Na_area.val.plot + P_area.val.plot) / 
  (Zn_area.val.plot + guide_area() + guide_area()) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()
