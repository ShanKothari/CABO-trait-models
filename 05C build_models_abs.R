setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(pls)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
source("Scripts/VIP.R")

abs.train<-readRDS("ProcessedSpectra/abs_train.rds")
abs.test<-readRDS("ProcessedSpectra/abs_test.rds")

##########################################
## to dos

## try Type II regression?
## add together Chl a and b?

#########################################
## define functions

## RMSD between predicted and observed values
RMSD<-function(measured,predicted){
  not.na<-which(!is.na(measured) & !is.na(predicted))
  return(sqrt(sum((measured-predicted)^2,na.rm=T)/(length(not.na)-1)))
}

## percent RMSD (based on data quantiles)
## set min and max to 0 and 1 for range as denominator
## or to 0.25 and 0.75 for IQR as denominator
percentRMSD<-function(measured,predicted,min,max,na.rm=T){
  RMSD_data<-RMSD(measured,predicted)
  range<-unname(quantile(measured,probs=max,na.rm=na.rm)-quantile(measured,probs=min,na.rm=na.rm))
  return(RMSD_data/range)
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

###########################################
## start building first-pass models using
## K-fold cross validation to get the
## optimal number of components

solubles_mass_CVmodel<-plsr(meta(abs.train)$solubles_mass~as.matrix(abs.train),
                            ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_mass_CVmodel <- selectNcomp(solubles_mass_CVmodel, method = "onesigma", plot = FALSE)
solubles_mass_valid <- which(!is.na(meta(abs.train)$solubles_mass))
solubles_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[solubles_mass_valid],
                               Species=meta(abs.train)$species[solubles_mass_valid],
                               Project=meta(abs.train)$project[solubles_mass_valid],
                               measured=meta(abs.train)$solubles_mass[solubles_mass_valid],
                               val_pred=solubles_mass_CVmodel$validation$pred[,,ncomp_solubles_mass_CVmodel])
ggplot(solubles_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting % solubles from fresh-leaf spectra")

hemicellulose_mass_CVmodel<-plsr(meta(abs.train)$hemicellulose_mass~as.matrix(abs.train),
                                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_mass_CVmodel <- selectNcomp(hemicellulose_mass_CVmodel, method = "onesigma", plot = FALSE)
hemicellulose_mass_valid <- which(!is.na(meta(abs.train)$hemicellulose_mass))
hemicellulose_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[hemicellulose_mass_valid],
                                    Species=meta(abs.train)$species[hemicellulose_mass_valid],
                                    Project=meta(abs.train)$project[hemicellulose_mass_valid],
                                    measured=meta(abs.train)$hemicellulose_mass[hemicellulose_mass_valid],
                                    val_pred=hemicellulose_mass_CVmodel$validation$pred[,,ncomp_hemicellulose_mass_CVmodel])
ggplot(hemicellulose_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting % hemicellulose from fresh-leaf spectra")

cellulose_mass_CVmodel<-plsr(meta(abs.train)$cellulose_mass~as.matrix(abs.train),
                             ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_mass_CVmodel <- selectNcomp(cellulose_mass_CVmodel, method = "onesigma", plot = FALSE)
cellulose_mass_valid <- which(!is.na(meta(abs.train)$cellulose_mass))
cellulose_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[cellulose_mass_valid],
                                Species=meta(abs.train)$species[cellulose_mass_valid],
                                Project=meta(abs.train)$project[cellulose_mass_valid],
                                measured=meta(abs.train)$cellulose_mass[cellulose_mass_valid],
                                val_pred=cellulose_mass_CVmodel$validation$pred[,,ncomp_cellulose_mass_CVmodel])
ggplot(cellulose_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting % cellulose from fresh-leaf spectra")

lignin_mass_CVmodel<-plsr(meta(abs.train)$lignin_mass~as.matrix(abs.train),
                          ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_mass_CVmodel <- selectNcomp(lignin_mass_CVmodel, method = "onesigma", plot = FALSE)
lignin_mass_valid <- which(!is.na(meta(abs.train)$lignin_mass))
lignin_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[lignin_mass_valid],
                             Species=meta(abs.train)$species[lignin_mass_valid],
                             Project=meta(abs.train)$project[lignin_mass_valid],
                             measured=meta(abs.train)$lignin_mass[lignin_mass_valid],
                             val_pred=lignin_mass_CVmodel$validation$pred[,,ncomp_lignin_mass_CVmodel])
ggplot(lignin_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting % lignin from fresh-leaf spectra")

Cmass_CVmodel<-plsr(meta(abs.train)$Cmass~as.matrix(abs.train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cmass_CVmodel <- selectNcomp(Cmass_CVmodel, method = "onesigma", plot = FALSE)
Cmass_valid <- which(!is.na(meta(abs.train)$Cmass))
Cmass_pred<-data.frame(ID=meta(abs.train)$sample_id[Cmass_valid],
                       Species=meta(abs.train)$species[Cmass_valid],
                       Project=meta(abs.train)$project[Cmass_valid],
                       measured=meta(abs.train)$Cmass[Cmass_valid],
                       val_pred=Cmass_CVmodel$validation$pred[,,ncomp_Cmass_CVmodel])
ggplot(Cmass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting %C from fresh-leaf spectra")

Nmass_CVmodel<-plsr(meta(abs.train)$Nmass~as.matrix(abs.train),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Nmass_CVmodel <- selectNcomp(Nmass_CVmodel, method = "onesigma", plot = FALSE)
Nmass_valid <- which(!is.na(meta(abs.train)$Nmass))
Nmass_pred<-data.frame(ID=meta(abs.train)$sample_id[Nmass_valid],
                       Species=meta(abs.train)$species[Nmass_valid],
                       Project=meta(abs.train)$project[Nmass_valid],
                       measured=meta(abs.train)$Nmass[Nmass_valid],
                       val_pred=Nmass_CVmodel$validation$pred[,,ncomp_Nmass_CVmodel])
ggplot(Nmass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting %N from fresh-leaf spectra")

# Narea_CVmodel<-plsr(meta(abs.train)$Narea~as.matrix(abs.train),
#                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
# ncomp_Narea_CVmodel <- selectNcomp(Narea_CVmodel, method = "onesigma", plot = FALSE)
# Narea_valid <- which(!is.na(meta(abs.train)$Narea))
# Narea_pred<-data.frame(ID=meta(abs.train)$sample_id[Narea_valid],
#                        Species=meta(abs.train)$species[Narea_valid],
#                        Project=meta(abs.train)$project[Narea_valid],
#                        measured=meta(abs.train)$Narea[Narea_valid],
#                        val_pred=Narea_CVmodel$validation$pred[,,ncomp_Narea_CVmodel])
# ggplot(Narea_pred,aes(y=measured,x=val_pred,color=Project))+
#   geom_point(size=2)+geom_smooth(method="lm",se=F)+
#   theme_bw()+
#   geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
#   ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
#   theme(text = element_text(size=20),
#         legend.position = c(0.8, 0.2))+
#   labs(y="Measured",x="Predicted")+
#   ggtitle("Predicting Narea from fresh-leaf spectra")
# 
# Nnorm_CVmodel<-plsr(meta(abs.train)$Nnorm~as.matrix(abs.train),
#                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
# ncomp_Nnorm_CVmodel <- selectNcomp(Nnorm_CVmodel, method = "onesigma", plot = FALSE)
# Nnorm_valid <- which(!is.na(meta(abs.train)$Nnorm))
# Nnorm_pred<-data.frame(ID=meta(abs.train)$sample_id[Nnorm_valid],
#                        Species=meta(abs.train)$species[Nnorm_valid],
#                        Project=meta(abs.train)$project[Nnorm_valid],
#                        measured=meta(abs.train)$Nnorm[Nnorm_valid],
#                        val_pred=Nnorm_CVmodel$validation$pred[,,ncomp_Nnorm_CVmodel])
# ggplot(Nnorm_pred,aes(y=measured,x=val_pred,color=Project))+
#   geom_point(size=2)+geom_smooth(method="lm",se=F)+
#   theme_bw()+
#   geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
#   ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
#   theme(text = element_text(size=20),
#         legend.position = c(0.8, 0.2))+
#   labs(y="Measured",x="Predicted")+
#   ggtitle("Predicting N (norm-ind) from fresh-leaf spectra")

EWT_CVmodel<-plsr(meta(abs.train)$EWT~as.matrix(abs.train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_CVmodel <- selectNcomp(EWT_CVmodel, method = "onesigma", plot = FALSE)
EWT_valid <- which(!is.na(meta(abs.train)$EWT))
EWT_pred<-data.frame(ID=meta(abs.train)$sample_id[EWT_valid],
                     Species=meta(abs.train)$species[EWT_valid],
                     Project=meta(abs.train)$project[EWT_valid],
                     measured=meta(abs.train)$EWT[EWT_valid],
                     val_pred=EWT_CVmodel$validation$pred[,,ncomp_EWT_CVmodel])
ggplot(EWT_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting %EWT from fresh-leaf spectra")

LDMC_CVmodel<-plsr(meta(abs.train)$LDMC~as.matrix(abs.train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LDMC_CVmodel <- selectNcomp(LDMC_CVmodel, method = "onesigma", plot = FALSE)
LDMC_valid <- which(!is.na(meta(abs.train)$LDMC))
LDMC_pred<-data.frame(ID=meta(abs.train)$sample_id[LDMC_valid],
                      Species=meta(abs.train)$species[LDMC_valid],
                      Project=meta(abs.train)$project[LDMC_valid],
                      measured=meta(abs.train)$LDMC[LDMC_valid],
                      val_pred=LDMC_CVmodel$validation$pred[,,ncomp_LDMC_CVmodel])
ggplot(LDMC_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting LDMC from fresh-leaf spectra")

LMA_CVmodel<-plsr(meta(abs.train)$LMA~as.matrix(abs.train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_CVmodel <- selectNcomp(LMA_CVmodel, method = "onesigma", plot = FALSE)
LMA_valid <- which(!is.na(meta(abs.train)$LMA))
LMA_pred<-data.frame(ID=meta(abs.train)$sample_id[LMA_valid],
                     Species=meta(abs.train)$species[LMA_valid],
                     Project=meta(abs.train)$project[LMA_valid],
                     measured=meta(abs.train)$LMA[LMA_valid],
                     val_pred=LMA_CVmodel$validation$pred[,,ncomp_LMA_CVmodel])
ggplot(LMA_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting LMA from fresh-leaf spectra")

chlA_mass_CVmodel<-plsr(meta(abs.train)$chlA_mass~as.matrix(abs.train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_mass_CVmodel <- selectNcomp(chlA_mass_CVmodel, method = "onesigma", plot = FALSE)
chlA_mass_valid <- which(!is.na(meta(abs.train)$chlA_mass))
chlA_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[chlA_mass_valid],
                           Species=meta(abs.train)$species[chlA_mass_valid],
                           Project=meta(abs.train)$project[chlA_mass_valid],
                           measured=meta(abs.train)$chlA_mass[chlA_mass_valid],
                           val_pred=chlA_mass_CVmodel$validation$pred[,,ncomp_chlA_mass_CVmodel])
ggplot(chlA_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting chlA from fresh-leaf spectra")

chlB_mass_CVmodel<-plsr(meta(abs.train)$chlB_mass~as.matrix(abs.train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_mass_CVmodel <- selectNcomp(chlB_mass_CVmodel, method = "onesigma", plot = FALSE)
chlB_mass_valid <- which(!is.na(meta(abs.train)$chlB_mass))
chlB_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[chlB_mass_valid],
                           Species=meta(abs.train)$species[chlB_mass_valid],
                           Project=meta(abs.train)$project[chlB_mass_valid],
                           measured=meta(abs.train)$chlB_mass[chlB_mass_valid],
                           val_pred=chlB_mass_CVmodel$validation$pred[,,ncomp_chlB_mass_CVmodel])
ggplot(chlB_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting chlB from fresh-leaf spectra")

car_mass_CVmodel<-plsr(meta(abs.train)$car_mass~as.matrix(abs.train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_mass_CVmodel <- selectNcomp(car_mass_CVmodel, method = "onesigma", plot = FALSE)
car_mass_valid <- which(!is.na(meta(abs.train)$car_mass))
car_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[car_mass_valid],
                          Species=meta(abs.train)$species[car_mass_valid],
                          Project=meta(abs.train)$project[car_mass_valid],
                          measured=meta(abs.train)$car_mass[car_mass_valid],
                          val_pred=car_mass_CVmodel$validation$pred[,,ncomp_car_mass_CVmodel])
ggplot(car_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting carotenoids from fresh-leaf spectra")

Al_mass_CVmodel<-plsr(meta(abs.train)$Al_mass~as.matrix(abs.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_mass_CVmodel <- selectNcomp(Al_mass_CVmodel, method = "onesigma", plot = FALSE)
Al_mass_valid <- which(!is.na(meta(abs.train)$Al_mass))
Al_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[Al_mass_valid],
                         Species=meta(abs.train)$species[Al_mass_valid],
                         Project=meta(abs.train)$project[Al_mass_valid],
                         measured=meta(abs.train)$Al_mass[Al_mass_valid],
                         val_pred=Al_mass_CVmodel$validation$pred[,,ncomp_Al_mass_CVmodel])
ggplot(Al_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Al from fresh-leaf spectra")

Ca_mass_CVmodel<-plsr(meta(abs.train)$Ca_mass~as.matrix(abs.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_mass_CVmodel <- selectNcomp(Ca_mass_CVmodel, method = "onesigma", plot = FALSE)
Ca_mass_valid <- which(!is.na(meta(abs.train)$Ca_mass))
Ca_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[Ca_mass_valid],
                         Species=meta(abs.train)$species[Ca_mass_valid],
                         Project=meta(abs.train)$project[Ca_mass_valid],
                         measured=meta(abs.train)$Ca_mass[Ca_mass_valid],
                         val_pred=Ca_mass_CVmodel$validation$pred[,,ncomp_Ca_mass_CVmodel])
ggplot(Ca_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Ca from fresh-leaf spectra")

Cu_mass_CVmodel<-plsr(meta(abs.train)$Cu_mass~as.matrix(abs.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_mass_CVmodel <- selectNcomp(Cu_mass_CVmodel, method = "onesigma", plot = FALSE)
Cu_mass_valid <- which(!is.na(meta(abs.train)$Cu_mass))
Cu_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[Cu_mass_valid],
                         Species=meta(abs.train)$species[Cu_mass_valid],
                         Project=meta(abs.train)$project[Cu_mass_valid],
                         measured=meta(abs.train)$Cu_mass[Cu_mass_valid],
                         val_pred=Cu_mass_CVmodel$validation$pred[,,ncomp_Cu_mass_CVmodel])
ggplot(Cu_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Cu from fresh-leaf spectra")

Fe_mass_CVmodel<-plsr(meta(abs.train)$Fe_mass~as.matrix(abs.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_mass_CVmodel <- selectNcomp(Fe_mass_CVmodel, method = "onesigma", plot = FALSE)
Fe_mass_valid <- which(!is.na(meta(abs.train)$Fe_mass))
Fe_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[Fe_mass_valid],
                         Species=meta(abs.train)$species[Fe_mass_valid],
                         Project=meta(abs.train)$project[Fe_mass_valid],
                         measured=meta(abs.train)$Fe_mass[Fe_mass_valid],
                         val_pred=Fe_mass_CVmodel$validation$pred[,,ncomp_Fe_mass_CVmodel])
ggplot(Fe_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Fe from fresh-leaf spectra")

K_mass_CVmodel<-plsr(meta(abs.train)$K_mass~as.matrix(abs.train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_mass_CVmodel <- selectNcomp(K_mass_CVmodel, method = "onesigma", plot = FALSE)
K_mass_valid <- which(!is.na(meta(abs.train)$K_mass))
K_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[K_mass_valid],
                        Species=meta(abs.train)$species[K_mass_valid],
                        Project=meta(abs.train)$project[K_mass_valid],
                        measured=meta(abs.train)$K_mass[K_mass_valid],
                        val_pred=K_mass_CVmodel$validation$pred[,,ncomp_K_mass_CVmodel])
ggplot(K_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting K from fresh-leaf spectra")

Mg_mass_CVmodel<-plsr(meta(abs.train)$Mg_mass~as.matrix(abs.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_mass_CVmodel <- selectNcomp(Mg_mass_CVmodel, method = "onesigma", plot = FALSE)
Mg_mass_valid <- which(!is.na(meta(abs.train)$Mg_mass))
Mg_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[Mg_mass_valid],
                         Species=meta(abs.train)$species[Mg_mass_valid],
                         Project=meta(abs.train)$project[Mg_mass_valid],
                         measured=meta(abs.train)$Mg_mass[Mg_mass_valid],
                         val_pred=Mg_mass_CVmodel$validation$pred[,,ncomp_Mg_mass_CVmodel])
ggplot(Mg_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Mg from fresh-leaf spectra")

Mn_mass_CVmodel<-plsr(meta(abs.train)$Mn_mass~as.matrix(abs.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_mass_CVmodel <- selectNcomp(Mn_mass_CVmodel, method = "onesigma", plot = FALSE)
Mn_mass_valid <- which(!is.na(meta(abs.train)$Mn_mass))
Mn_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[Mn_mass_valid],
                         Species=meta(abs.train)$species[Mn_mass_valid],
                         Project=meta(abs.train)$project[Mn_mass_valid],
                         measured=meta(abs.train)$Mn_mass[Mn_mass_valid],
                         val_pred=Mn_mass_CVmodel$validation$pred[,,ncomp_Mn_mass_CVmodel])
ggplot(Mn_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Mn from fresh-leaf spectra")

Na_mass_CVmodel<-plsr(meta(abs.train)$Na_mass~as.matrix(abs.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_mass_CVmodel <- selectNcomp(Na_mass_CVmodel, method = "onesigma", plot = FALSE)
Na_mass_valid <- which(!is.na(meta(abs.train)$Na_mass))
Na_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[Na_mass_valid],
                         Species=meta(abs.train)$species[Na_mass_valid],
                         Project=meta(abs.train)$project[Na_mass_valid],
                         measured=meta(abs.train)$Na_mass[Na_mass_valid],
                         val_pred=Na_mass_CVmodel$validation$pred[,,ncomp_Na_mass_CVmodel])
ggplot(Na_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Na from fresh-leaf spectra")

P_mass_CVmodel<-plsr(meta(abs.train)$P_mass~as.matrix(abs.train),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_mass_CVmodel <- selectNcomp(P_mass_CVmodel, method = "onesigma", plot = FALSE)
P_mass_valid <- which(!is.na(meta(abs.train)$P_mass))
P_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[P_mass_valid],
                        Species=meta(abs.train)$species[P_mass_valid],
                        Project=meta(abs.train)$project[P_mass_valid],
                        measured=meta(abs.train)$P_mass[P_mass_valid],
                        val_pred=P_mass_CVmodel$validation$pred[,,ncomp_P_mass_CVmodel])
ggplot(P_mass_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting P from fresh-leaf spectra")

Zn_mass_CVmodel<-plsr(meta(abs.train)$Zn_mass~as.matrix(abs.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_mass_CVmodel <- selectNcomp(Zn_mass_CVmodel, method = "onesigma", plot = FALSE)
Zn_mass_valid <- which(!is.na(meta(abs.train)$Zn_mass))
Zn_mass_pred<-data.frame(ID=meta(abs.train)$sample_id[Zn_mass_valid],
                         Species=meta(abs.train)$species[Zn_mass_valid],
                         Project=meta(abs.train)$project[Zn_mass_valid],
                         measured=meta(abs.train)$Zn_mass[Zn_mass_valid],
                         val_pred=Zn_mass_CVmodel$validation$pred[,,ncomp_Zn_mass_CVmodel])
ggplot(Zn_mass_pred,aes(y=measured,x=val_pred,color=Project))+
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

VIP1.df<-data.frame(Cmass=VIP(Cmass_CVmodel)[ncomp_Cmass_CVmodel,],
                    Nmass=VIP(Nmass_CVmodel)[ncomp_Nmass_CVmodel,],
                    LMA=VIP(LMA_CVmodel)[ncomp_LMA_CVmodel,],
                    LDMC=VIP(LDMC_CVmodel)[ncomp_LDMC_CVmodel,],
                    EWT=VIP(EWT_CVmodel)[ncomp_EWT_CVmodel,],
                    wavelength=400:2400)

VIP2.df<-data.frame(NDFmass=VIP(NDFmass_CVmodel)[ncomp_NDFmass_CVmodel,],
                    ADFmass=VIP(ADFmass_CVmodel)[ncomp_ADFmass_CVmodel,],
                    ADLmass=VIP(ADLmass_CVmodel)[ncomp_ADLmass_CVmodel,],
                    wavelength=400:2400)

VIPN.df<-data.frame(Nmass=VIP(Nmass_CVmodel)[ncomp_Nmass_CVmodel,],
                    Narea=VIP(Narea_CVmodel)[ncomp_Narea_CVmodel,],
                    Nnorm=VIP(Nnorm_CVmodel)[ncomp_Nnorm_CVmodel,],
                    LMA=VIP(LMA_CVmodel)[ncomp_LMA_CVmodel,],
                    wavelength=400:2400)

VIP1.long<-melt(VIP1.df,id.vars = "wavelength")
VIP2.long<-melt(VIP2.df,id.vars = "wavelength")
VIPN.long<-melt(VIPN.df,id.vars = "wavelength")

VIP1.plot<-ggplot(VIP1.long,aes(x=wavelength,y=value,color=variable))+
  geom_line()+theme_bw()+
  theme(text=element_text(size=20))+
  labs(y="PLSR VIP",x="Wavelength")

VIP2.plot<-ggplot(VIP2.long,aes(x=wavelength,y=value,color=variable))+
  geom_line()+theme_bw()+
  theme(text=element_text(size=20))+
  labs(y="PLSR VIP",x="Wavelength")

VIPN.plot<-ggplot(VIPN.long,aes(x=wavelength,y=value,color=variable))+
  geom_line()+theme_bw()+
  theme(text=element_text(size=20))+
  labs(y="PLSR VIP",x="Wavelength")

#######################################################
## jackknife analyses

solubles_mass.jack.coefs<-list()
hemicellulose_mass.jack.coefs<-list()
cellulose_mass.jack.coefs<-list()
lignin_mass.jack.coefs<-list()
chlA_mass.jack.coefs<-list()
chlB_mass.jack.coefs<-list()
car_mass.jack.coefs<-list()
Cmass.jack.coefs<-list()
Nmass.jack.coefs<-list()
EWT.jack.coefs<-list()
LMA.jack.coefs<-list()
LDMC.jack.coefs<-list()
Al_mass.jack.coefs<-list()
Ca_mass.jack.coefs<-list()
Cu_mass.jack.coefs<-list()
Fe_mass.jack.coefs<-list()
K_mass.jack.coefs<-list()
Mg_mass.jack.coefs<-list()
Mn_mass.jack.coefs<-list()
Na_mass.jack.coefs<-list()
P_mass.jack.coefs<-list()
Zn_mass.jack.coefs<-list()

solubles_mass.jack.stats<-list()
hemicellulose_mass.jack.stats<-list()
cellulose_mass.jack.stats<-list()
lignin_mass.jack.stats<-list()
chlA_mass.jack.stats<-list()
chlB_mass.jack.stats<-list()
car_mass.jack.stats<-list()
Cmass.jack.stats<-list()
Nmass.jack.stats<-list()
EWT.jack.stats<-list()
LMA.jack.stats<-list()
LDMC.jack.stats<-list()
Al_mass.jack.stats<-list()
Ca_mass.jack.stats<-list()
Cu_mass.jack.stats<-list()
Fe_mass.jack.stats<-list()
K_mass.jack.stats<-list()
Mg_mass.jack.stats<-list()
Mn_mass.jack.stats<-list()
Na_mass.jack.stats<-list()
P_mass.jack.stats<-list()
Zn_mass.jack.stats<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n.cal.spec<-nrow(abs.train)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-abs.train[train.jack]
  val.jack<-abs.train[test.jack]
  
  solubles_mass.jack<-plsr(meta(calib.jack)$solubles_mass~as.matrix(calib.jack),
                           ncomp=30,method = "oscorespls",validation="none")
  hemicellulose_mass.jack<-plsr(meta(calib.jack)$hemicellulose_mass~as.matrix(calib.jack),
                                ncomp=30,method = "oscorespls",validation="none")
  cellulose_mass.jack<-plsr(meta(calib.jack)$cellulose_mass~as.matrix(calib.jack),
                            ncomp=30,method = "oscorespls",validation="none")
  lignin_mass.jack<-plsr(meta(calib.jack)$lignin_mass~as.matrix(calib.jack),
                         ncomp=30,method = "oscorespls",validation="none")
  chlA_mass.jack<-plsr(meta(calib.jack)$chlA_mass~as.matrix(calib.jack),
                       ncomp=30,method = "oscorespls",validation="none")
  chlB_mass.jack<-plsr(meta(calib.jack)$chlB_mass~as.matrix(calib.jack),
                       ncomp=30,method = "oscorespls",validation="none")
  car_mass.jack<-plsr(meta(calib.jack)$car_mass~as.matrix(calib.jack),
                      ncomp=30,method = "oscorespls",validation="none")
  Cmass.jack<-plsr(meta(calib.jack)$Cmass~as.matrix(calib.jack),
                   ncomp=30,method = "oscorespls",validation="none")
  Nmass.jack<-plsr(meta(calib.jack)$Nmass~as.matrix(calib.jack),
                   ncomp=30,method = "oscorespls",validation="none")
  EWT.jack<-plsr(meta(calib.jack)$EWT~as.matrix(calib.jack),
                 ncomp=30,method = "oscorespls",validation="none")
  LMA.jack<-plsr(meta(calib.jack)$LMA~as.matrix(calib.jack),
                 ncomp=30,method = "oscorespls",validation="none")
  LDMC.jack<-plsr(meta(calib.jack)$LDMC~as.matrix(calib.jack),
                  ncomp=30,method = "oscorespls",validation="none")
  Al_mass.jack<-plsr(meta(calib.jack)$Al_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Ca_mass.jack<-plsr(meta(calib.jack)$Ca_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Cu_mass.jack<-plsr(meta(calib.jack)$Cu_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Fe_mass.jack<-plsr(meta(calib.jack)$Fe_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  K_mass.jack<-plsr(meta(calib.jack)$K_mass~as.matrix(calib.jack),
                    ncomp=30,method = "oscorespls",validation="none")
  Mg_mass.jack<-plsr(meta(calib.jack)$Mg_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Mn_mass.jack<-plsr(meta(calib.jack)$Mn_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  Na_mass.jack<-plsr(meta(calib.jack)$Na_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  P_mass.jack<-plsr(meta(calib.jack)$P_mass~as.matrix(calib.jack),
                    ncomp=30,method = "oscorespls",validation="none")
  Zn_mass.jack<-plsr(meta(calib.jack)$Zn_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  
  solubles_mass.jack.val.pred<-as.vector(predict(solubles_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_solubles_mass_CVmodel)[,,1])
  solubles_mass.jack.val.fit<-lm(solubles_mass.jack.val.pred~meta(val.jack)$solubles_mass)
  solubles_mass.jack.stats[[i]]<-c(R2=summary(solubles_mass.jack.val.fit)$r.squared,
                                   RMSE=RMSD(meta(val.jack)$solubles_mass,solubles_mass.jack.val.pred),
                                   perRMSE=percentRMSD(meta(val.jack)$solubles_mass,solubles_mass.jack.val.pred,0.025,0.975),
                                   bias=mean(solubles_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$solubles_mass,na.rm=T))
  
  hemicellulose_mass.jack.val.pred<-as.vector(predict(hemicellulose_mass.jack,newdata=as.matrix(val.jack),
                                                      ncomp=ncomp_hemicellulose_mass_CVmodel)[,,1])
  hemicellulose_mass.jack.val.fit<-lm(hemicellulose_mass.jack.val.pred~meta(val.jack)$hemicellulose_mass)
  hemicellulose_mass.jack.stats[[i]]<-c(R2=summary(hemicellulose_mass.jack.val.fit)$r.squared,
                                        RMSE=RMSD(meta(val.jack)$hemicellulose_mass,hemicellulose_mass.jack.val.pred),
                                        perRMSE=percentRMSD(meta(val.jack)$hemicellulose_mass,
                                                                hemicellulose_mass.jack.val.pred,0.025,0.975),
                                        bias=mean(hemicellulose_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$hemicellulose_mass,na.rm=T))
  
  cellulose_mass.jack.val.pred<-as.vector(predict(cellulose_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_cellulose_mass_CVmodel)[,,1])
  cellulose_mass.jack.val.fit<-lm(cellulose_mass.jack.val.pred~meta(val.jack)$cellulose_mass)
  cellulose_mass.jack.stats[[i]]<-c(R2=summary(cellulose_mass.jack.val.fit)$r.squared,
                                    RMSE=RMSD(meta(val.jack)$cellulose_mass,cellulose_mass.jack.val.pred),
                                    perRMSE=percentRMSD(meta(val.jack)$cellulose_mass,cellulose_mass.jack.val.pred,0.025,0.975),
                                    bias=mean(cellulose_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$cellulose_mass,na.rm=T))
  
  lignin_mass.jack.val.pred<-as.vector(predict(lignin_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_lignin_mass_CVmodel)[,,1])
  lignin_mass.jack.val.fit<-lm(lignin_mass.jack.val.pred~meta(val.jack)$lignin_mass)
  lignin_mass.jack.stats[[i]]<-c(R2=summary(lignin_mass.jack.val.fit)$r.squared,
                                 RMSE=RMSD(meta(val.jack)$lignin_mass,lignin_mass.jack.val.pred),
                                 perRMSE=percentRMSD(meta(val.jack)$lignin_mass,lignin_mass.jack.val.pred,0.025,0.975),
                                 bias=mean(lignin_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$lignin_mass,na.rm=T))
  
  chlA_mass.jack.val.pred<-as.vector(predict(chlA_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlA_mass_CVmodel)[,,1])
  chlA_mass.jack.val.fit<-lm(chlA_mass.jack.val.pred~meta(val.jack)$chlA_mass)
  chlA_mass.jack.stats[[i]]<-c(R2=summary(chlA_mass.jack.val.fit)$r.squared,
                               RMSE=RMSD(meta(val.jack)$chlA_mass,chlA_mass.jack.val.pred),
                               perRMSE=percentRMSD(meta(val.jack)$chlA_mass,chlA_mass.jack.val.pred,0.025,0.975),
                               bias=mean(chlA_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlA_mass,na.rm=T))
  
  chlB_mass.jack.val.pred<-as.vector(predict(chlB_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlB_mass_CVmodel)[,,1])
  chlB_mass.jack.val.fit<-lm(chlB_mass.jack.val.pred~meta(val.jack)$chlB_mass)
  chlB_mass.jack.stats[[i]]<-c(R2=summary(chlB_mass.jack.val.fit)$r.squared,
                               RMSE=RMSD(meta(val.jack)$chlB_mass,chlB_mass.jack.val.pred),
                               perRMSE=percentRMSD(meta(val.jack)$chlB_mass,chlB_mass.jack.val.pred,0.025,0.975),
                               bias=mean(chlB_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlB_mass,na.rm=T))
  
  car_mass.jack.val.pred<-as.vector(predict(car_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_car_mass_CVmodel)[,,1])
  car_mass.jack.val.fit<-lm(car_mass.jack.val.pred~meta(val.jack)$car_mass)
  car_mass.jack.stats[[i]]<-c(R2=summary(car_mass.jack.val.fit)$r.squared,
                              RMSE=RMSD(meta(val.jack)$car_mass,car_mass.jack.val.pred),
                              perRMSE=percentRMSD(meta(val.jack)$car_mass,car_mass.jack.val.pred,0.025,0.975),
                              bias=mean(car_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$car_mass,na.rm=T))
  
  Cmass.jack.val.pred<-as.vector(predict(Cmass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Cmass_CVmodel)[,,1])
  Cmass.jack.val.fit<-lm(Cmass.jack.val.pred~meta(val.jack)$Cmass)
  Cmass.jack.stats[[i]]<-c(R2=summary(Cmass.jack.val.fit)$r.squared,
                           RMSE=RMSD(meta(val.jack)$Cmass,Cmass.jack.val.pred),
                           perRMSE=percentRMSD(meta(val.jack)$Cmass,Cmass.jack.val.pred,0.025,0.975),
                           bias=mean(Cmass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Cmass,na.rm=T))
  
  Nmass.jack.val.pred<-as.vector(predict(Nmass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Nmass_CVmodel)[,,1])
  Nmass.jack.val.fit<-lm(Nmass.jack.val.pred~meta(val.jack)$Nmass)
  Nmass.jack.stats[[i]]<-c(R2=summary(Nmass.jack.val.fit)$r.squared,
                           RMSE=RMSD(meta(val.jack)$Nmass,Nmass.jack.val.pred),
                           perRMSE=percentRMSD(meta(val.jack)$Nmass,Nmass.jack.val.pred,0.025,0.975),
                           bias=mean(Nmass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Nmass,na.rm=T))
  
  EWT.jack.val.pred<-as.vector(predict(EWT.jack,newdata=as.matrix(val.jack),ncomp=ncomp_EWT_CVmodel)[,,1])
  EWT.jack.val.fit<-lm(EWT.jack.val.pred~meta(val.jack)$EWT)
  EWT.jack.stats[[i]]<-c(R2=summary(EWT.jack.val.fit)$r.squared,
                         RMSE=RMSD(meta(val.jack)$EWT,EWT.jack.val.pred),
                         perRMSE=percentRMSD(meta(val.jack)$EWT,EWT.jack.val.pred,0.025,0.975),
                         bias=mean(EWT.jack.val.pred,na.rm=T)-mean(meta(val.jack)$EWT,na.rm=T))
  
  LMA.jack.val.pred<-as.vector(predict(LMA.jack,newdata=as.matrix(val.jack),ncomp=ncomp_LMA_CVmodel)[,,1])
  LMA.jack.val.fit<-lm(LMA.jack.val.pred~meta(val.jack)$LMA)
  LMA.jack.stats[[i]]<-c(R2=summary(LMA.jack.val.fit)$r.squared,
                         RMSE=RMSD(meta(val.jack)$LMA,LMA.jack.val.pred),
                         perRMSE=percentRMSD(meta(val.jack)$LMA,LMA.jack.val.pred,0.025,0.975),
                         bias=mean(LMA.jack.val.pred,na.rm=T)-mean(meta(val.jack)$LMA,na.rm=T))
  
  LDMC.jack.val.pred<-as.vector(predict(LDMC.jack,newdata=as.matrix(val.jack),ncomp=ncomp_LDMC_CVmodel)[,,1])
  LDMC.jack.val.fit<-lm(LDMC.jack.val.pred~meta(val.jack)$LDMC)
  LDMC.jack.stats[[i]]<-c(R2=summary(LDMC.jack.val.fit)$r.squared,
                          RMSE=RMSD(meta(val.jack)$LDMC,LDMC.jack.val.pred),
                          perRMSE=percentRMSD(meta(val.jack)$LDMC,LDMC.jack.val.pred,0.025,0.975),
                          bias=mean(LDMC.jack.val.pred,na.rm=T)-mean(meta(val.jack)$LDMC,na.rm=T))
  
  Al_mass.jack.val.pred<-as.vector(predict(Al_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Al_mass_CVmodel)[,,1])
  Al_mass.jack.val.fit<-lm(Al_mass.jack.val.pred~meta(val.jack)$Al_mass)
  Al_mass.jack.stats[[i]]<-c(R2=summary(Al_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Al_mass,Al_mass.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$Al_mass,Al_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Al_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Al_mass,na.rm=T))
  
  Ca_mass.jack.val.pred<-as.vector(predict(Ca_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Ca_mass_CVmodel)[,,1])
  Ca_mass.jack.val.fit<-lm(Ca_mass.jack.val.pred~meta(val.jack)$Ca_mass)
  Ca_mass.jack.stats[[i]]<-c(R2=summary(Ca_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Ca_mass,Ca_mass.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$Ca_mass,Ca_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Ca_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Ca_mass,na.rm=T))
  
  Cu_mass.jack.val.pred<-as.vector(predict(Cu_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Cu_mass_CVmodel)[,,1])
  Cu_mass.jack.val.fit<-lm(Cu_mass.jack.val.pred~meta(val.jack)$Cu_mass)
  Cu_mass.jack.stats[[i]]<-c(R2=summary(Cu_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Cu_mass,Cu_mass.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$Cu_mass,Cu_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Cu_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Cu_mass,na.rm=T))
  
  Fe_mass.jack.val.pred<-as.vector(predict(Fe_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Fe_mass_CVmodel)[,,1])
  Fe_mass.jack.val.fit<-lm(Fe_mass.jack.val.pred~meta(val.jack)$Fe_mass)
  Fe_mass.jack.stats[[i]]<-c(R2=summary(Fe_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Fe_mass,Fe_mass.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$Fe_mass,Fe_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Fe_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Fe_mass,na.rm=T))
  
  K_mass.jack.val.pred<-as.vector(predict(K_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_K_mass_CVmodel)[,,1])
  K_mass.jack.val.fit<-lm(K_mass.jack.val.pred~meta(val.jack)$K_mass)
  K_mass.jack.stats[[i]]<-c(R2=summary(K_mass.jack.val.fit)$r.squared,
                            RMSE=RMSD(meta(val.jack)$K_mass,K_mass.jack.val.pred),
                            perRMSE=percentRMSD(meta(val.jack)$K_mass,K_mass.jack.val.pred,0.025,0.975),
                            bias=mean(K_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$K_mass,na.rm=T))
  
  Mg_mass.jack.val.pred<-as.vector(predict(Mg_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Mg_mass_CVmodel)[,,1])
  Mg_mass.jack.val.fit<-lm(Mg_mass.jack.val.pred~meta(val.jack)$Mg_mass)
  Mg_mass.jack.stats[[i]]<-c(R2=summary(Mg_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Mg_mass,Mg_mass.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$Mg_mass,Mg_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Mg_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Mg_mass,na.rm=T))
  
  Mn_mass.jack.val.pred<-as.vector(predict(Mn_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Mn_mass_CVmodel)[,,1])
  Mn_mass.jack.val.fit<-lm(Mn_mass.jack.val.pred~meta(val.jack)$Mn_mass)
  Mn_mass.jack.stats[[i]]<-c(R2=summary(Mn_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Mn_mass,Mn_mass.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$Mn_mass,Mn_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Mn_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Mn_mass,na.rm=T))
  
  Na_mass.jack.val.pred<-as.vector(predict(Na_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Na_mass_CVmodel)[,,1])
  Na_mass.jack.val.fit<-lm(Na_mass.jack.val.pred~meta(val.jack)$Na_mass)
  Na_mass.jack.stats[[i]]<-c(R2=summary(Na_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Na_mass,Na_mass.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$Na_mass,Na_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Na_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Na_mass,na.rm=T))
  
  P_mass.jack.val.pred<-as.vector(predict(P_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_P_mass_CVmodel)[,,1])
  P_mass.jack.val.fit<-lm(P_mass.jack.val.pred~meta(val.jack)$P_mass)
  P_mass.jack.stats[[i]]<-c(R2=summary(P_mass.jack.val.fit)$r.squared,
                            RMSE=RMSD(meta(val.jack)$P_mass,P_mass.jack.val.pred),
                            perRMSE=percentRMSD(meta(val.jack)$P_mass,P_mass.jack.val.pred,0.025,0.975),
                            bias=mean(P_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$P_mass,na.rm=T))
  
  Zn_mass.jack.val.pred<-as.vector(predict(Zn_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Zn_mass_CVmodel)[,,1])
  Zn_mass.jack.val.fit<-lm(Zn_mass.jack.val.pred~meta(val.jack)$Zn_mass)
  Zn_mass.jack.stats[[i]]<-c(R2=summary(Zn_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Zn_mass,Zn_mass.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$Zn_mass,Zn_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Zn_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Zn_mass,na.rm=T))
  
  solubles_mass.jack.coefs[[i]]<-as.vector(coef(solubles_mass.jack,ncomp=ncomp_solubles_mass_CVmodel,intercept=TRUE))
  hemicellulose_mass.jack.coefs[[i]]<-as.vector(coef(hemicellulose_mass.jack,ncomp=ncomp_hemicellulose_mass_CVmodel,intercept=TRUE))
  cellulose_mass.jack.coefs[[i]]<-as.vector(coef(cellulose_mass.jack,ncomp=ncomp_cellulose_mass_CVmodel,intercept=TRUE))
  lignin_mass.jack.coefs[[i]]<-as.vector(coef(lignin_mass.jack,ncomp=ncomp_lignin_mass_CVmodel,intercept=TRUE))
  chlA_mass.jack.coefs[[i]]<-as.vector(coef(chlA_mass.jack,ncomp=ncomp_chlA_mass_CVmodel,intercept=TRUE))
  chlB_mass.jack.coefs[[i]]<-as.vector(coef(chlB_mass.jack,ncomp=ncomp_chlB_mass_CVmodel,intercept=TRUE))
  car_mass.jack.coefs[[i]]<-as.vector(coef(car_mass.jack,ncomp=ncomp_car_mass_CVmodel,intercept=TRUE))
  Cmass.jack.coefs[[i]]<-as.vector(coef(Cmass.jack,ncomp=ncomp_Cmass_CVmodel,intercept=TRUE))
  Nmass.jack.coefs[[i]]<-as.vector(coef(Nmass.jack,ncomp=ncomp_Nmass_CVmodel,intercept=TRUE))
  EWT.jack.coefs[[i]]<-as.vector(coef(EWT.jack,ncomp=ncomp_EWT_CVmodel,intercept=TRUE))
  LMA.jack.coefs[[i]]<-as.vector(coef(LMA.jack,ncomp=ncomp_LMA_CVmodel,intercept=TRUE))
  LDMC.jack.coefs[[i]]<-as.vector(coef(LDMC.jack,ncomp=ncomp_LDMC_CVmodel,intercept=TRUE))
  Al_mass.jack.coefs[[i]]<-as.vector(coef(Al_mass.jack,ncomp=ncomp_Al_mass_CVmodel,intercept=TRUE))
  Ca_mass.jack.coefs[[i]]<-as.vector(coef(Ca_mass.jack,ncomp=ncomp_Ca_mass_CVmodel,intercept=TRUE))
  Cu_mass.jack.coefs[[i]]<-as.vector(coef(Cu_mass.jack,ncomp=ncomp_Cu_mass_CVmodel,intercept=TRUE))
  Fe_mass.jack.coefs[[i]]<-as.vector(coef(Fe_mass.jack,ncomp=ncomp_Fe_mass_CVmodel,intercept=TRUE))
  K_mass.jack.coefs[[i]]<-as.vector(coef(K_mass.jack,ncomp=ncomp_K_mass_CVmodel,intercept=TRUE))
  Mg_mass.jack.coefs[[i]]<-as.vector(coef(Mg_mass.jack,ncomp=ncomp_Mg_mass_CVmodel,intercept=TRUE))
  Mn_mass.jack.coefs[[i]]<-as.vector(coef(Mn_mass.jack,ncomp=ncomp_Mn_mass_CVmodel,intercept=TRUE))
  Na_mass.jack.coefs[[i]]<-as.vector(coef(Na_mass.jack,ncomp=ncomp_Na_mass_CVmodel,intercept=TRUE))
  P_mass.jack.coefs[[i]]<-as.vector(coef(P_mass.jack,ncomp=ncomp_P_mass_CVmodel,intercept=TRUE))
  Zn_mass.jack.coefs[[i]]<-as.vector(coef(Zn_mass.jack,ncomp=ncomp_Zn_mass_CVmodel,intercept=TRUE))
  
}

solubles_mass.jack.pred<-apply.coefs(solubles_mass.jack.coefs,as.matrix(abs.test))
solubles_mass.jack.stat<-t(apply(solubles_mass.jack.pred,1,
                                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_mass.jack.df<-data.frame(pred.mean=solubles_mass.jack.stat[,1],
                                  pred.low=solubles_mass.jack.stat[,2],
                                  pred.high=solubles_mass.jack.stat[,3],
                                  Measured=meta(abs.test)$solubles_mass,
                                  ncomp=ncomp_solubles_mass_CVmodel,
                                  Species=meta(abs.test)$species,
                                  Project=meta(abs.test)$project,
                                  functional.group=meta(abs.test)$functional.group,
                                  ID=meta(abs.test)$sample_id)

solubles_mass.val.plot<-ggplot(solubles_mass.jack.df,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,95),ylim=c(25,95))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)

hemicellulose_mass.jack.pred<-apply.coefs(hemicellulose_mass.jack.coefs,as.matrix(abs.test))
hemicellulose_mass.jack.stat<-t(apply(hemicellulose_mass.jack.pred,1,
                                      function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_mass.jack.df<-data.frame(pred.mean=hemicellulose_mass.jack.stat[,1],
                                       pred.low=hemicellulose_mass.jack.stat[,2],
                                       pred.high=hemicellulose_mass.jack.stat[,3],
                                       Measured=meta(abs.test)$hemicellulose_mass,
                                       ncomp=ncomp_hemicellulose_mass_CVmodel,
                                       Species=meta(abs.test)$species,
                                       Project=meta(abs.test)$project,
                                       functional.group=meta(abs.test)$functional.group,
                                       ID=meta(abs.test)$sample_id)

hemicellulose_mass.val.plot<-ggplot(hemicellulose_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,44),ylim=c(0,44))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)

cellulose_mass.jack.pred<-apply.coefs(cellulose_mass.jack.coefs,as.matrix(abs.test))
cellulose_mass.jack.stat<-t(apply(cellulose_mass.jack.pred,1,
                                  function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_mass.jack.df<-data.frame(pred.mean=cellulose_mass.jack.stat[,1],
                                   pred.low=cellulose_mass.jack.stat[,2],
                                   pred.high=cellulose_mass.jack.stat[,3],
                                   Measured=meta(abs.test)$cellulose_mass,
                                   ncomp=ncomp_cellulose_mass_CVmodel,
                                   Species=meta(abs.test)$species,
                                   Project=meta(abs.test)$project,
                                   functional.group=meta(abs.test)$functional.group,
                                   ID=meta(abs.test)$sample_id)

cellulose_mass.val.plot<-ggplot(cellulose_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,34),ylim=c(0,34))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)

lignin_mass.jack.pred<-apply.coefs(lignin_mass.jack.coefs,as.matrix(abs.test))
lignin_mass.jack.stat<-t(apply(lignin_mass.jack.pred,1,
                               function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_mass.jack.df<-data.frame(pred.mean=lignin_mass.jack.stat[,1],
                                pred.low=lignin_mass.jack.stat[,2],
                                pred.high=lignin_mass.jack.stat[,3],
                                Measured=meta(abs.test)$lignin_mass,
                                ncomp=ncomp_lignin_mass_CVmodel,
                                Species=meta(abs.test)$species,
                                Project=meta(abs.test)$project,
                                functional.group=meta(abs.test)$functional.group,
                                ID=meta(abs.test)$sample_id)

lignin_mass.val.plot<-ggplot(lignin_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-2.5,23.5),ylim=c(-2.5,23.5))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)

chlA_mass.jack.pred<-apply.coefs(chlA_mass.jack.coefs,as.matrix(abs.test))
chlA_mass.jack.stat<-t(apply(chlA_mass.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_mass.jack.df<-data.frame(pred.mean=chlA_mass.jack.stat[,1],
                              pred.low=chlA_mass.jack.stat[,2],
                              pred.high=chlA_mass.jack.stat[,3],
                              Measured=meta(abs.test)$chlA_mass,
                              ncomp=ncomp_chlA_mass_CVmodel,
                              Species=meta(abs.test)$species,
                              Project=meta(abs.test)$project,
                              functional.group=meta(abs.test)$functional.group,
                              ID=meta(abs.test)$sample_id)

chlA_mass.val.plot<-ggplot(chlA_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-2,25),ylim=c(-2,25))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)

chlB_mass.jack.pred<-apply.coefs(chlB_mass.jack.coefs,as.matrix(abs.test))
chlB_mass.jack.stat<-t(apply(chlB_mass.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_mass.jack.df<-data.frame(pred.mean=chlB_mass.jack.stat[,1],
                              pred.low=chlB_mass.jack.stat[,2],
                              pred.high=chlB_mass.jack.stat[,3],
                              Measured=meta(abs.test)$chlB_mass,
                              ncomp=ncomp_chlB_mass_CVmodel,
                              Species=meta(abs.test)$species,
                              Project=meta(abs.test)$project,
                              functional.group=meta(abs.test)$functional.group,
                              ID=meta(abs.test)$sample_id)

chlB_mass.val.plot<-ggplot(chlB_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-1.5,9),ylim=c(-1.5,9))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)

car_mass.jack.pred<-apply.coefs(car_mass.jack.coefs,as.matrix(abs.test))
car_mass.jack.stat<-t(apply(car_mass.jack.pred,1,
                            function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_mass.jack.df<-data.frame(pred.mean=car_mass.jack.stat[,1],
                             pred.low=car_mass.jack.stat[,2],
                             pred.high=car_mass.jack.stat[,3],
                             Measured=meta(abs.test)$car_mass,
                             ncomp=ncomp_car_mass_CVmodel,
                             Species=meta(abs.test)$species,
                             Project=meta(abs.test)$project,
                             functional.group=meta(abs.test)$functional.group,
                             ID=meta(abs.test)$sample_id)

car_mass.val.plot<-ggplot(car_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.5,5.5),ylim=c(-0.5,5.5))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"),
       color="Functional group")

Cmass.jack.pred<-apply.coefs(Cmass.jack.coefs,as.matrix(abs.test))
Cmass.jack.stat<-t(apply(Cmass.jack.pred,1,
                         function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cmass.jack.df<-data.frame(pred.mean=Cmass.jack.stat[,1],
                          pred.low=Cmass.jack.stat[,2],
                          pred.high=Cmass.jack.stat[,3],
                          Measured=meta(abs.test)$Cmass,
                          ncomp=ncomp_Cmass_CVmodel,
                          Species=meta(abs.test)$species,
                          Project=meta(abs.test)$project,
                          functional.group=meta(abs.test)$functional.group,
                          ID=meta(abs.test)$sample_id)

Cmass.val.plot<-ggplot(Cmass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  guides(color=F)

Nmass.jack.pred<-apply.coefs(Nmass.jack.coefs,as.matrix(abs.test))
Nmass.jack.stat<-t(apply(Nmass.jack.pred,1,
                         function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Nmass.jack.df<-data.frame(pred.mean=Nmass.jack.stat[,1],
                          pred.low=Nmass.jack.stat[,2],
                          pred.high=Nmass.jack.stat[,3],
                          Measured=meta(abs.test)$Nmass,
                          ncomp=ncomp_Nmass_CVmodel,
                          Species=meta(abs.test)$species,
                          Project=meta(abs.test)$project,
                          functional.group=meta(abs.test)$functional.group,
                          ID=meta(abs.test)$sample_id)

Nmass.val.plot<-ggplot(Nmass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,5.5),ylim=c(0,5.5))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)

EWT.jack.pred<-apply.coefs(EWT.jack.coefs,as.matrix(abs.test))
EWT.jack.stat<-t(apply(EWT.jack.pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT.jack.df<-data.frame(pred.mean=EWT.jack.stat[,1],
                        pred.low=EWT.jack.stat[,2],
                        pred.high=EWT.jack.stat[,3],
                        Measured=meta(abs.test)$EWT,
                        ncomp=ncomp_EWT_CVmodel,
                        Species=meta(abs.test)$species,
                        Project=meta(abs.test)$project,
                        functional.group=meta(abs.test)$functional.group,
                        ID=meta(abs.test)$sample_id)

EWT.val.plot<-ggplot(EWT.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.68),ylim=c(0,0.68))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)

LMA.jack.pred<-apply.coefs(LMA.jack.coefs,as.matrix(abs.test))
LMA.jack.stat<-t(apply(LMA.jack.pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA.jack.df<-data.frame(pred.mean=LMA.jack.stat[,1],
                        pred.low=LMA.jack.stat[,2],
                        pred.high=LMA.jack.stat[,3],
                        Measured=meta(abs.test)$LMA,
                        ncomp=ncomp_LMA_CVmodel,
                        Species=meta(abs.test)$species,
                        Project=meta(abs.test)$project,
                        functional.group=meta(abs.test)$functional.group,
                        ID=meta(abs.test)$sample_id)

LMA.val.plot<-ggplot(LMA.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.05,0.4),ylim=c(-0.05,0.4))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)

LDMC.jack.pred<-apply.coefs(LDMC.jack.coefs,as.matrix(abs.test))
LDMC.jack.stat<-t(apply(LDMC.jack.pred,1,
                        function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LDMC.jack.df<-data.frame(pred.mean=LDMC.jack.stat[,1],
                         pred.low=LDMC.jack.stat[,2],
                         pred.high=LDMC.jack.stat[,3],
                         Measured=meta(abs.test)$LDMC,
                         ncomp=ncomp_LDMC_CVmodel,
                         Species=meta(abs.test)$species,
                         Project=meta(abs.test)$project,
                         functional.group=meta(abs.test)$functional.group,
                         ID=meta(abs.test)$sample_id)

LDMC.val.plot<-ggplot(LDMC.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,600),ylim=c(0,600))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)

Al_mass.jack.pred<-apply.coefs(Al_mass.jack.coefs,as.matrix(abs.test))
Al_mass.jack.stat<-t(apply(Al_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_mass.jack.df<-data.frame(pred.mean=Al_mass.jack.stat[,1],
                            pred.low=Al_mass.jack.stat[,2],
                            pred.high=Al_mass.jack.stat[,3],
                            Measured=meta(abs.test)$Al_mass,
                            ncomp=ncomp_Al_mass_CVmodel,
                            Species=meta(abs.test)$species,
                            Project=meta(abs.test)$project,
                            functional.group=meta(abs.test)$functional.group,
                            ID=meta(abs.test)$sample_id)

Al_mass.val.plot<-ggplot(Al_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.025,0.25),ylim=c(-0.025,0.25))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Al (mg g"^-1*")"),x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)

Ca_mass.jack.pred<-apply.coefs(Ca_mass.jack.coefs,as.matrix(abs.test))
Ca_mass.jack.stat<-t(apply(Ca_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_mass.jack.df<-data.frame(pred.mean=Ca_mass.jack.stat[,1],
                            pred.low=Ca_mass.jack.stat[,2],
                            pred.high=Ca_mass.jack.stat[,3],
                            Measured=meta(abs.test)$Ca_mass,
                            ncomp=ncomp_Ca_mass_CVmodel,
                            Species=meta(abs.test)$species,
                            Project=meta(abs.test)$project,
                            functional.group=meta(abs.test)$functional.group,
                            ID=meta(abs.test)$sample_id)

Ca_mass.val.plot<-ggplot(Ca_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-5,35),ylim=c(-5,35))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Ca (mg g"^-1*")"),
       x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)

Cu_mass.jack.pred<-apply.coefs(Cu_mass.jack.coefs,as.matrix(abs.test))
Cu_mass.jack.stat<-t(apply(Cu_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_mass.jack.df<-data.frame(pred.mean=Cu_mass.jack.stat[,1],
                            pred.low=Cu_mass.jack.stat[,2],
                            pred.high=Cu_mass.jack.stat[,3],
                            Measured=meta(abs.test)$Cu_mass,
                            ncomp=ncomp_Cu_mass_CVmodel,
                            Species=meta(abs.test)$species,
                            Project=meta(abs.test)$project,
                            functional.group=meta(abs.test)$functional.group,
                            ID=meta(abs.test)$sample_id)

Cu_mass.val.plot<-ggplot(Cu_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.002,0.02),ylim=c(-0.002,0.02))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Cu (mg g"^-1*")"),
       x=expression("Predicted Cu (mg g"^-1*")"),
       color="Functional group")

Fe_mass.jack.pred<-apply.coefs(Fe_mass.jack.coefs,as.matrix(abs.test))
Fe_mass.jack.stat<-t(apply(Fe_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_mass.jack.df<-data.frame(pred.mean=Fe_mass.jack.stat[,1],
                            pred.low=Fe_mass.jack.stat[,2],
                            pred.high=Fe_mass.jack.stat[,3],
                            Measured=meta(abs.test)$Fe_mass,
                            ncomp=ncomp_Fe_mass_CVmodel,
                            Species=meta(abs.test)$species,
                            Project=meta(abs.test)$project,
                            functional.group=meta(abs.test)$functional.group,
                            ID=meta(abs.test)$sample_id)

Fe_mass.val.plot<-ggplot(Fe_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.17),ylim=c(0,0.17))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Fe (mg g"^-1*")"),
       x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)

K_mass.jack.pred<-apply.coefs(K_mass.jack.coefs,as.matrix(abs.test))
K_mass.jack.stat<-t(apply(K_mass.jack.pred,1,
                          function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_mass.jack.df<-data.frame(pred.mean=K_mass.jack.stat[,1],
                           pred.low=K_mass.jack.stat[,2],
                           pred.high=K_mass.jack.stat[,3],
                           Measured=meta(abs.test)$K_mass,
                           ncomp=ncomp_K_mass_CVmodel,
                           Species=meta(abs.test)$species,
                           Project=meta(abs.test)$project,
                           functional.group=meta(abs.test)$functional.group,
                           ID=meta(abs.test)$sample_id)

K_mass.val.plot<-ggplot(K_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,18),ylim=c(0,18))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)

Mg_mass.jack.pred<-apply.coefs(Mg_mass.jack.coefs,as.matrix(abs.test))
Mg_mass.jack.stat<-t(apply(Mg_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_mass.jack.df<-data.frame(pred.mean=Mg_mass.jack.stat[,1],
                            pred.low=Mg_mass.jack.stat[,2],
                            pred.high=Mg_mass.jack.stat[,3],
                            Measured=meta(abs.test)$Mg_mass,
                            ncomp=ncomp_Mg_mass_CVmodel,
                            Species=meta(abs.test)$species,
                            Project=meta(abs.test)$project,
                            functional.group=meta(abs.test)$functional.group,
                            ID=meta(abs.test)$sample_id)

Mg_mass.val.plot<-ggplot(Mg_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.5,7),ylim=c(-0.5,7))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Mg (mg g"^-1*")"),
       x=expression("Predicted Mg (mg g"^-1*")"))+
  guides(color=F)

Mn_mass.jack.pred<-apply.coefs(Mn_mass.jack.coefs,as.matrix(abs.test))
Mn_mass.jack.stat<-t(apply(Mn_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_mass.jack.df<-data.frame(pred.mean=Mn_mass.jack.stat[,1],
                            pred.low=Mn_mass.jack.stat[,2],
                            pred.high=Mn_mass.jack.stat[,3],
                            Measured=meta(abs.test)$Mn_mass,
                            ncomp=ncomp_Mn_mass_CVmodel,
                            Species=meta(abs.test)$species,
                            Project=meta(abs.test)$project,
                            functional.group=meta(abs.test)$functional.group,
                            ID=meta(abs.test)$sample_id)

Mn_mass.val.plot<-ggplot(Mn_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,1.1),ylim=c(-0.1,1.1))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Mn (mg g"^-1*")"),
       x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)

Na_mass.jack.pred<-apply.coefs(Na_mass.jack.coefs,as.matrix(abs.test))
Na_mass.jack.stat<-t(apply(Na_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_mass.jack.df<-data.frame(pred.mean=Na_mass.jack.stat[,1],
                            pred.low=Na_mass.jack.stat[,2],
                            pred.high=Na_mass.jack.stat[,3],
                            Measured=meta(abs.test)$Na_mass,
                            ncomp=ncomp_Na_mass_CVmodel,
                            Species=meta(abs.test)$species,
                            Project=meta(abs.test)$project,
                            functional.group=meta(abs.test)$functional.group,
                            ID=meta(abs.test)$sample_id)

Na_mass.val.plot<-ggplot(Na_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.3,4.2),ylim=c(-0.3,4.2))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Na (mg g"^-1*")"),
       x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)

P_mass.jack.pred<-apply.coefs(P_mass.jack.coefs,as.matrix(abs.test))
P_mass.jack.stat<-t(apply(P_mass.jack.pred,1,
                          function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_mass.jack.df<-data.frame(pred.mean=P_mass.jack.stat[,1],
                           pred.low=P_mass.jack.stat[,2],
                           pred.high=P_mass.jack.stat[,3],
                           Measured=meta(abs.test)$P_mass,
                           ncomp=ncomp_P_mass_CVmodel,
                           Species=meta(abs.test)$species,
                           Project=meta(abs.test)$project,
                           functional.group=meta(abs.test)$functional.group,
                           ID=meta(abs.test)$sample_id)

P_mass.val.plot<-ggplot(P_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.5,7.5),ylim=c(-0.5,7.5))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured P (mg g"^-1*")"),
       x=expression("Predicted P (mg g"^-1*")"))+
  guides(color=F)

Zn_mass.jack.pred<-apply.coefs(Zn_mass.jack.coefs,as.matrix(abs.test))
Zn_mass.jack.stat<-t(apply(Zn_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_mass.jack.df<-data.frame(pred.mean=Zn_mass.jack.stat[,1],
                            pred.low=Zn_mass.jack.stat[,2],
                            pred.high=Zn_mass.jack.stat[,3],
                            Measured=meta(abs.test)$Zn_mass,
                            ncomp=ncomp_Zn_mass_CVmodel,
                            Species=meta(abs.test)$species,
                            Project=meta(abs.test)$project,
                            functional.group=meta(abs.test)$functional.group,
                            ID=meta(abs.test)$sample_id)

Zn_mass.val.plot<-ggplot(Zn_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.05,0.35),ylim=c(-0.05,0.35))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Zn (mg g"^-1*")"),
       x=expression("Predicted Zn (mg g"^-1*")"))+
  guides(color=F)

all.jack.coef.list<-list(sol=solubles_mass.jack.coefs,
                         hemi=hemicellulose_mass.jack.coefs,
                         cell=cellulose_mass.jack.coefs,
                         lign=lignin_mass.jack.coefs,
                         chlA=chlA_mass.jack.coefs,
                         chlB=chlB_mass.jack.coefs,
                         car=car_mass.jack.coefs,
                         N=Nmass.jack.coefs,
                         C=Cmass.jack.coefs,
                         EWT=EWT.jack.coefs,
                         LDMC=LDMC.jack.coefs,
                         LMA=LMA.jack.coefs,
                         Al=Al_mass.jack.coefs,
                         Ca=Ca_mass.jack.coefs,
                         Cu=Cu_mass.jack.coefs,
                         Fe=Fe_mass.jack.coefs,
                         K=K_mass.jack.coefs,
                         Mg=Mg_mass.jack.coefs,
                         Mn=Mn_mass.jack.coefs,
                         Na=Na_mass.jack.coefs,
                         P=P_mass.jack.coefs,
                         Zn=Zn_mass.jack.coefs)
saveRDS(all.jack.coef.list,"SavedResults/all_jack_coefs_list_abs.rds")

all.jack.df.list<-list(sol=solubles_mass.jack.df,
                       hemi=hemicellulose_mass.jack.df,
                       cell=cellulose_mass.jack.df,
                       lign=lignin_mass.jack.df,
                       chlA=chlA_mass.jack.df,
                       chlB=chlB_mass.jack.df,
                       car=car_mass.jack.df,
                       N=Nmass.jack.df,
                       C=Cmass.jack.df,
                       EWT=EWT.jack.df,
                       LDMC=LDMC.jack.df,
                       LMA=LMA.jack.df,
                       Al=Al_mass.jack.df,
                       Ca=Ca_mass.jack.df,
                       Cu=Cu_mass.jack.df,
                       Fe=Fe_mass.jack.df,
                       K=K_mass.jack.df,
                       Mg=Mg_mass.jack.df,
                       Mn=Mn_mass.jack.df,
                       Na=Na_mass.jack.df,
                       P=P_mass.jack.df,
                       Zn=Zn_mass.jack.df)
saveRDS(all.jack.df.list,"SavedResults/all_jack_df_list_abs.rds")

all.jack.stats.list<-list(sol=solubles_mass.jack.stats,
                          hemi=hemicellulose_mass.jack.stats,
                          cell=cellulose_mass.jack.stats,
                          lign=lignin_mass.jack.stats,
                          chlA=chlA_mass.jack.stats,
                          chlB=chlB_mass.jack.stats,
                          car=car_mass.jack.stats,
                          N=Nmass.jack.stats,
                          C=Cmass.jack.stats,
                          EWT=EWT.jack.stats,
                          LDMC=LDMC.jack.stats,
                          LMA=LMA.jack.stats,
                          Al=Al_mass.jack.stats,
                          Ca=Ca_mass.jack.stats,
                          Cu=Cu_mass.jack.stats,
                          Fe=Fe_mass.jack.stats,
                          K=K_mass.jack.stats,
                          Mg=Mg_mass.jack.stats,
                          Mn=Mn_mass.jack.stats,
                          Na=Na_mass.jack.stats,
                          P=P_mass.jack.stats,
                          Zn=Zn_mass.jack.stats)
saveRDS(all.jack.stats.list,"SavedResults/all_jack_stats_list_abs.rds")

pdf("Images/val_plots_abs1.pdf",width = 16,height = 20)
(EWT.val.plot + LDMC.val.plot + LMA.val.plot) / 
  (cellulose_mass.val.plot + solubles_mass.val.plot + Cmass.val.plot) / 
  (hemicellulose_mass.val.plot + Nmass.val.plot + chlB_mass.val.plot) / 
  (chlA_mass.val.plot + lignin_mass.val.plot + car_mass.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Images/val_plots_abs2.pdf",width = 16,height = 20)
(Ca_mass.val.plot + K_mass.val.plot + Zn_mass.val.plot) / 
  (P_mass.val.plot + Mg_mass.val.plot + Na_mass.val.plot) / 
  (Fe_mass.val.plot + Mn_mass.val.plot + Al_mass.val.plot) / 
  (Cu_mass.val.plot + guide_area() + guide_area()) &
  plot_layout(guides="collect")
dev.off()

######################
## violin plots

abs_val_summary<-data.frame(variable=names(all.jack.df.list),
                            perRMSE=unlist(lapply(all.jack.df.list,
                                                  function(x) percentRMSD(x$Measured,x$pred.mean,0.025,0.975))),
                            R2=unlist(lapply(all.jack.df.list,
                                             function(x) summary(lm(Measured~pred.mean,data=x))$r.squared)))

R2.df<-data.frame(LMA=unlist(lapply(LMA.jack.stats,function(x) x[["R2"]])),
                  LDMC=unlist(lapply(LDMC.jack.stats,function(x) x[["R2"]])),
                  EWT=unlist(lapply(EWT.jack.stats,function(x) x[["R2"]])),
                  sol=unlist(lapply(solubles_mass.jack.stats,function(x) x[["R2"]])),
                  hemi=unlist(lapply(hemicellulose_mass.jack.stats,function(x) x[["R2"]])),
                  cell=unlist(lapply(cellulose_mass.jack.stats,function(x) x[["R2"]])),
                  lign=unlist(lapply(lignin_mass.jack.stats,function(x) x[["R2"]])),
                  C=unlist(lapply(Cmass.jack.stats,function(x) x[["R2"]])),
                  N=unlist(lapply(Nmass.jack.stats,function(x) x[["R2"]])),
                  chlA=unlist(lapply(chlA_mass.jack.stats,function(x) x[["R2"]])),
                  chlB=unlist(lapply(chlB_mass.jack.stats,function(x) x[["R2"]])),
                  car=unlist(lapply(car_mass.jack.stats,function(x) x[["R2"]])),
                  Al=unlist(lapply(Al_mass.jack.stats,function(x) x[["R2"]])),
                  Ca=unlist(lapply(Ca_mass.jack.stats,function(x) x[["R2"]])),
                  Cu=unlist(lapply(Cu_mass.jack.stats,function(x) x[["R2"]])),
                  Fe=unlist(lapply(Fe_mass.jack.stats,function(x) x[["R2"]])),
                  K=unlist(lapply(K_mass.jack.stats,function(x) x[["R2"]])),
                  Mg=unlist(lapply(Mg_mass.jack.stats,function(x) x[["R2"]])),
                  Mn=unlist(lapply(Mn_mass.jack.stats,function(x) x[["R2"]])),
                  Na=unlist(lapply(Na_mass.jack.stats,function(x) x[["R2"]])),
                  P=unlist(lapply(P_mass.jack.stats,function(x) x[["R2"]])),
                  Zn=unlist(lapply(Zn_mass.jack.stats,function(x) x[["R2"]])))

R2.long<-melt(R2.df)
val_R2<-ggplot(R2.long,aes(y=value,x=variable))+
  geom_violin()+
  geom_point(data=abs_val_summary,
             aes(x=variable,y=R2),color="red",size=2)+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(y=expression(italic("R"^2)))+
  ggtitle("Fresh-leaf spectra")+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1))

perRMSE.df<-data.frame(LMA=unlist(lapply(LMA.jack.stats,function(x) 100*x[["perRMSE"]])),
                       LDMC=unlist(lapply(LDMC.jack.stats,function(x) 100*x[["perRMSE"]])),
                       EWT=unlist(lapply(EWT.jack.stats,function(x) 100*x[["perRMSE"]])),
                       sol=unlist(lapply(solubles_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       hemi=unlist(lapply(hemicellulose_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       cell=unlist(lapply(cellulose_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       lign=unlist(lapply(lignin_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       C=unlist(lapply(Cmass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       N=unlist(lapply(Nmass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       chlA=unlist(lapply(chlA_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       chlB=unlist(lapply(chlB_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       car=unlist(lapply(car_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       Al=unlist(lapply(Al_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       Ca=unlist(lapply(Ca_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       Cu=unlist(lapply(Cu_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       Fe=unlist(lapply(Fe_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       K=unlist(lapply(K_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       Mg=unlist(lapply(Mg_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       Mn=unlist(lapply(Mn_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       Na=unlist(lapply(Na_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       P=unlist(lapply(P_mass.jack.stats,function(x) 100*x[["perRMSE"]])),
                       Zn=unlist(lapply(Zn_mass.jack.stats,function(x) 100*x[["perRMSE"]])))

perRMSE.long<-melt(perRMSE.df)
val_perRMSE<-ggplot(perRMSE.long,aes(y=value,x=variable))+
  geom_violin()+
  geom_point(data=abs_val_summary,
             aes(x=variable,y=perRMSE*100),color="red",size=2)+
  theme_bw()+
  theme(text = element_text(size=20),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
  labs(y="%RMSE")+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,max(perRMSE.long$value)*1.1))

pdf("Images/violin_plots_abs.pdf",width=8,height=6,onefile=F)
egg::ggarrange(plots = list(val_R2,val_perRMSE),
               nrow=2,ncol=1)
dev.off()