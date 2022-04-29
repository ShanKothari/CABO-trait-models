setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(pls)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
source("Scripts/VIP.R")

ref.train<-readRDS("ProcessedSpectra/ref_train.rds")
ref.test<-readRDS("ProcessedSpectra/ref_test.rds")

ref.train.bn<-t(apply(as.matrix(ref.train),1,function(x) x/sqrt(sum(x^2))))
ref.train.bn<-spectra(ref.train.bn,bands=bands(ref.train),names=names(ref.train))
meta(ref.train.bn)<-meta(ref.train)

ref.test.bn<-t(apply(as.matrix(ref.test),1,function(x) x/sqrt(sum(x^2))))
ref.test.bn<-spectra(ref.test.bn,bands=bands(ref.test),names=names(ref.test))
meta(ref.test.bn)<-meta(ref.test)

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

solubles_mass_CVmodel<-plsr(meta(ref.train.bn)$solubles_mass~as.matrix(ref.train.bn),
                            ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_mass_CVmodel <- selectNcomp(solubles_mass_CVmodel, method = "onesigma", plot = FALSE)
solubles_mass_valid <- which(!is.na(meta(ref.train.bn)$solubles_mass))
solubles_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[solubles_mass_valid],
                               Species=meta(ref.train.bn)$species[solubles_mass_valid],
                               Project=meta(ref.train.bn)$project[solubles_mass_valid],
                               measured=meta(ref.train.bn)$solubles_mass[solubles_mass_valid],
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

hemicellulose_mass_CVmodel<-plsr(meta(ref.train.bn)$hemicellulose_mass~as.matrix(ref.train.bn),
                                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_mass_CVmodel <- selectNcomp(hemicellulose_mass_CVmodel, method = "onesigma", plot = FALSE)
hemicellulose_mass_valid <- which(!is.na(meta(ref.train.bn)$hemicellulose_mass))
hemicellulose_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[hemicellulose_mass_valid],
                                    Species=meta(ref.train.bn)$species[hemicellulose_mass_valid],
                                    Project=meta(ref.train.bn)$project[hemicellulose_mass_valid],
                                    measured=meta(ref.train.bn)$hemicellulose_mass[hemicellulose_mass_valid],
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

cellulose_mass_CVmodel<-plsr(meta(ref.train.bn)$cellulose_mass~as.matrix(ref.train.bn),
                             ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_mass_CVmodel <- selectNcomp(cellulose_mass_CVmodel, method = "onesigma", plot = FALSE)
cellulose_mass_valid <- which(!is.na(meta(ref.train.bn)$cellulose_mass))
cellulose_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[cellulose_mass_valid],
                                Species=meta(ref.train.bn)$species[cellulose_mass_valid],
                                Project=meta(ref.train.bn)$project[cellulose_mass_valid],
                                measured=meta(ref.train.bn)$cellulose_mass[cellulose_mass_valid],
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

lignin_mass_CVmodel<-plsr(meta(ref.train.bn)$lignin_mass~as.matrix(ref.train.bn),
                          ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_mass_CVmodel <- selectNcomp(lignin_mass_CVmodel, method = "onesigma", plot = FALSE)
lignin_mass_valid <- which(!is.na(meta(ref.train.bn)$lignin_mass))
lignin_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[lignin_mass_valid],
                             Species=meta(ref.train.bn)$species[lignin_mass_valid],
                             Project=meta(ref.train.bn)$project[lignin_mass_valid],
                             measured=meta(ref.train.bn)$lignin_mass[lignin_mass_valid],
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

Cmass_CVmodel<-plsr(meta(ref.train.bn)$Cmass~as.matrix(ref.train.bn),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cmass_CVmodel <- selectNcomp(Cmass_CVmodel, method = "onesigma", plot = FALSE)
Cmass_valid <- which(!is.na(meta(ref.train.bn)$Cmass))
Cmass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Cmass_valid],
                       Species=meta(ref.train.bn)$species[Cmass_valid],
                       Project=meta(ref.train.bn)$project[Cmass_valid],
                       measured=meta(ref.train.bn)$Cmass[Cmass_valid],
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

Nmass_CVmodel<-plsr(meta(ref.train.bn)$Nmass~as.matrix(ref.train.bn),
                    ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Nmass_CVmodel <- selectNcomp(Nmass_CVmodel, method = "onesigma", plot = FALSE)
Nmass_valid <- which(!is.na(meta(ref.train.bn)$Nmass))
Nmass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Nmass_valid],
                       Species=meta(ref.train.bn)$species[Nmass_valid],
                       Project=meta(ref.train.bn)$project[Nmass_valid],
                       measured=meta(ref.train.bn)$Nmass[Nmass_valid],
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

EWT_CVmodel<-plsr(meta(ref.train.bn)$EWT~as.matrix(ref.train.bn),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_CVmodel <- selectNcomp(EWT_CVmodel, method = "onesigma", plot = FALSE)
EWT_valid <- which(!is.na(meta(ref.train.bn)$EWT))
EWT_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[EWT_valid],
                     Species=meta(ref.train.bn)$species[EWT_valid],
                     Project=meta(ref.train.bn)$project[EWT_valid],
                     measured=meta(ref.train.bn)$EWT[EWT_valid],
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

LDMC_CVmodel<-plsr(meta(ref.train.bn)$LDMC~as.matrix(ref.train.bn),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LDMC_CVmodel <- selectNcomp(LDMC_CVmodel, method = "onesigma", plot = FALSE)
LDMC_valid <- which(!is.na(meta(ref.train.bn)$LDMC))
LDMC_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[LDMC_valid],
                      Species=meta(ref.train.bn)$species[LDMC_valid],
                      Project=meta(ref.train.bn)$project[LDMC_valid],
                      measured=meta(ref.train.bn)$LDMC[LDMC_valid],
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

LMA_CVmodel<-plsr(meta(ref.train.bn)$LMA~as.matrix(ref.train.bn),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_CVmodel <- selectNcomp(LMA_CVmodel, method = "onesigma", plot = FALSE)
LMA_valid <- which(!is.na(meta(ref.train.bn)$LMA))
LMA_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[LMA_valid],
                     Species=meta(ref.train.bn)$species[LMA_valid],
                     Project=meta(ref.train.bn)$project[LMA_valid],
                     measured=meta(ref.train.bn)$LMA[LMA_valid],
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

chlA_mass_CVmodel<-plsr(meta(ref.train.bn)$chlA_mass~as.matrix(ref.train.bn),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_mass_CVmodel <- selectNcomp(chlA_mass_CVmodel, method = "onesigma", plot = FALSE)
chlA_mass_valid <- which(!is.na(meta(ref.train.bn)$chlA_mass))
chlA_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[chlA_mass_valid],
                           Species=meta(ref.train.bn)$species[chlA_mass_valid],
                           Project=meta(ref.train.bn)$project[chlA_mass_valid],
                           measured=meta(ref.train.bn)$chlA_mass[chlA_mass_valid],
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

chlB_mass_CVmodel<-plsr(meta(ref.train.bn)$chlB_mass~as.matrix(ref.train.bn),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_mass_CVmodel <- selectNcomp(chlB_mass_CVmodel, method = "onesigma", plot = FALSE)
chlB_mass_valid <- which(!is.na(meta(ref.train.bn)$chlB_mass))
chlB_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[chlB_mass_valid],
                           Species=meta(ref.train.bn)$species[chlB_mass_valid],
                           Project=meta(ref.train.bn)$project[chlB_mass_valid],
                           measured=meta(ref.train.bn)$chlB_mass[chlB_mass_valid],
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

car_mass_CVmodel<-plsr(meta(ref.train.bn)$car_mass~as.matrix(ref.train.bn),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_mass_CVmodel <- selectNcomp(car_mass_CVmodel, method = "onesigma", plot = FALSE)
car_mass_valid <- which(!is.na(meta(ref.train.bn)$car_mass))
car_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[car_mass_valid],
                          Species=meta(ref.train.bn)$species[car_mass_valid],
                          Project=meta(ref.train.bn)$project[car_mass_valid],
                          measured=meta(ref.train.bn)$car_mass[car_mass_valid],
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

Al_mass_CVmodel<-plsr(meta(ref.train.bn)$Al_mass~as.matrix(ref.train.bn),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_mass_CVmodel <- selectNcomp(Al_mass_CVmodel, method = "onesigma", plot = FALSE)
Al_mass_valid <- which(!is.na(meta(ref.train.bn)$Al_mass))
Al_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Al_mass_valid],
                         Species=meta(ref.train.bn)$species[Al_mass_valid],
                         Project=meta(ref.train.bn)$project[Al_mass_valid],
                         measured=meta(ref.train.bn)$Al_mass[Al_mass_valid],
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

Ca_mass_CVmodel<-plsr(meta(ref.train.bn)$Ca_mass~as.matrix(ref.train.bn),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_mass_CVmodel <- selectNcomp(Ca_mass_CVmodel, method = "onesigma", plot = FALSE)
Ca_mass_valid <- which(!is.na(meta(ref.train.bn)$Ca_mass))
Ca_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Ca_mass_valid],
                         Species=meta(ref.train.bn)$species[Ca_mass_valid],
                         Project=meta(ref.train.bn)$project[Ca_mass_valid],
                         measured=meta(ref.train.bn)$Ca_mass[Ca_mass_valid],
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

Cu_mass_CVmodel<-plsr(meta(ref.train.bn)$Cu_mass~as.matrix(ref.train.bn),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_mass_CVmodel <- selectNcomp(Cu_mass_CVmodel, method = "onesigma", plot = FALSE)
Cu_mass_valid <- which(!is.na(meta(ref.train.bn)$Cu_mass))
Cu_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Cu_mass_valid],
                         Species=meta(ref.train.bn)$species[Cu_mass_valid],
                         Project=meta(ref.train.bn)$project[Cu_mass_valid],
                         measured=meta(ref.train.bn)$Cu_mass[Cu_mass_valid],
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

Fe_mass_CVmodel<-plsr(meta(ref.train.bn)$Fe_mass~as.matrix(ref.train.bn),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_mass_CVmodel <- selectNcomp(Fe_mass_CVmodel, method = "onesigma", plot = FALSE)
Fe_mass_valid <- which(!is.na(meta(ref.train.bn)$Fe_mass))
Fe_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Fe_mass_valid],
                         Species=meta(ref.train.bn)$species[Fe_mass_valid],
                         Project=meta(ref.train.bn)$project[Fe_mass_valid],
                         measured=meta(ref.train.bn)$Fe_mass[Fe_mass_valid],
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

K_mass_CVmodel<-plsr(meta(ref.train.bn)$K_mass~as.matrix(ref.train.bn),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_mass_CVmodel <- selectNcomp(K_mass_CVmodel, method = "onesigma", plot = FALSE)
K_mass_valid <- which(!is.na(meta(ref.train.bn)$K_mass))
K_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[K_mass_valid],
                        Species=meta(ref.train.bn)$species[K_mass_valid],
                        Project=meta(ref.train.bn)$project[K_mass_valid],
                        measured=meta(ref.train.bn)$K_mass[K_mass_valid],
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

Mg_mass_CVmodel<-plsr(meta(ref.train.bn)$Mg_mass~as.matrix(ref.train.bn),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_mass_CVmodel <- selectNcomp(Mg_mass_CVmodel, method = "onesigma", plot = FALSE)
Mg_mass_valid <- which(!is.na(meta(ref.train.bn)$Mg_mass))
Mg_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Mg_mass_valid],
                         Species=meta(ref.train.bn)$species[Mg_mass_valid],
                         Project=meta(ref.train.bn)$project[Mg_mass_valid],
                         measured=meta(ref.train.bn)$Mg_mass[Mg_mass_valid],
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

Mn_mass_CVmodel<-plsr(meta(ref.train.bn)$Mn_mass~as.matrix(ref.train.bn),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_mass_CVmodel <- selectNcomp(Mn_mass_CVmodel, method = "onesigma", plot = FALSE)
Mn_mass_valid <- which(!is.na(meta(ref.train.bn)$Mn_mass))
Mn_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Mn_mass_valid],
                         Species=meta(ref.train.bn)$species[Mn_mass_valid],
                         Project=meta(ref.train.bn)$project[Mn_mass_valid],
                         measured=meta(ref.train.bn)$Mn_mass[Mn_mass_valid],
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

Na_mass_CVmodel<-plsr(meta(ref.train.bn)$Na_mass~as.matrix(ref.train.bn),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_mass_CVmodel <- selectNcomp(Na_mass_CVmodel, method = "onesigma", plot = FALSE)
Na_mass_valid <- which(!is.na(meta(ref.train.bn)$Na_mass))
Na_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Na_mass_valid],
                         Species=meta(ref.train.bn)$species[Na_mass_valid],
                         Project=meta(ref.train.bn)$project[Na_mass_valid],
                         measured=meta(ref.train.bn)$Na_mass[Na_mass_valid],
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

P_mass_CVmodel<-plsr(meta(ref.train.bn)$P_mass~as.matrix(ref.train.bn),
                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_mass_CVmodel <- selectNcomp(P_mass_CVmodel, method = "onesigma", plot = FALSE)
P_mass_valid <- which(!is.na(meta(ref.train.bn)$P_mass))
P_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[P_mass_valid],
                        Species=meta(ref.train.bn)$species[P_mass_valid],
                        Project=meta(ref.train.bn)$project[P_mass_valid],
                        measured=meta(ref.train.bn)$P_mass[P_mass_valid],
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

Zn_mass_CVmodel<-plsr(meta(ref.train.bn)$Zn_mass~as.matrix(ref.train.bn),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_mass_CVmodel <- selectNcomp(Zn_mass_CVmodel, method = "onesigma", plot = FALSE)
Zn_mass_valid <- which(!is.na(meta(ref.train.bn)$Zn_mass))
Zn_mass_pred<-data.frame(ID=meta(ref.train.bn)$sample_id[Zn_mass_valid],
                         Species=meta(ref.train.bn)$species[Zn_mass_valid],
                         Project=meta(ref.train.bn)$project[Zn_mass_valid],
                         measured=meta(ref.train.bn)$Zn_mass[Zn_mass_valid],
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
  
  n.cal.spec<-nrow(ref.train.bn)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-ref.train.bn[train.jack]
  val.jack<-ref.train.bn[test.jack]
  
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
                                        perRMSE=percentRMSD(meta(val.jack)$hemicellulose_mass,hemicellulose_mass.jack.val.pred,0.025,0.975),
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

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

solubles_mass.jack.pred<-apply.coefs(solubles_mass.jack.coefs,as.matrix(ref.test.bn))
solubles_mass.jack.stat<-t(apply(solubles_mass.jack.pred,1,
                                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_mass.jack.df<-data.frame(pred.mean=solubles_mass.jack.stat[,1],
                                  pred.low=solubles_mass.jack.stat[,2],
                                  pred.high=solubles_mass.jack.stat[,3],
                                  Measured=meta(ref.test.bn)$solubles_mass,
                                  ncomp=ncomp_solubles_mass_CVmodel,
                                  Species=meta(ref.test.bn)$species,
                                  Project=meta(ref.test.bn)$project,
                                  functional.group=meta(ref.test.bn)$functional.group,
                                  ID=meta(ref.test.bn)$sample_id)
solubles_all<-with(solubles_mass.jack.df,c(pred.low,pred.high,Measured))
solubles_upper<-max(solubles_all,na.rm=T)+1
solubles_lower<-min(solubles_all,na.rm=T)-1

solubles_mass.val.plot<-ggplot(solubles_mass.jack.df,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),
                  ylim=c(solubles_lower,solubles_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

hemicellulose_mass.jack.pred<-apply.coefs(hemicellulose_mass.jack.coefs,as.matrix(ref.test.bn))
hemicellulose_mass.jack.stat<-t(apply(hemicellulose_mass.jack.pred,1,
                                      function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_mass.jack.df<-data.frame(pred.mean=hemicellulose_mass.jack.stat[,1],
                                       pred.low=hemicellulose_mass.jack.stat[,2],
                                       pred.high=hemicellulose_mass.jack.stat[,3],
                                       Measured=meta(ref.test.bn)$hemicellulose_mass,
                                       ncomp=ncomp_hemicellulose_mass_CVmodel,
                                       Species=meta(ref.test.bn)$species,
                                       Project=meta(ref.test.bn)$project,
                                       functional.group=meta(ref.test.bn)$functional.group,
                                       ID=meta(ref.test.bn)$sample_id)
hemicellulose_all<-with(hemicellulose_mass.jack.df,c(pred.low,pred.high,Measured))
hemicellulose_upper<-max(hemicellulose_all,na.rm=T)+1
hemicellulose_lower<-min(hemicellulose_all,na.rm=T)-1

hemicellulose_mass.val.plot<-ggplot(hemicellulose_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),
                  ylim=c(hemicellulose_lower,hemicellulose_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured hemicellulose (%)",
       x="Predicted hemicellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

cellulose_mass.jack.pred<-apply.coefs(cellulose_mass.jack.coefs,as.matrix(ref.test.bn))
cellulose_mass.jack.stat<-t(apply(cellulose_mass.jack.pred,1,
                                  function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_mass.jack.df<-data.frame(pred.mean=cellulose_mass.jack.stat[,1],
                                   pred.low=cellulose_mass.jack.stat[,2],
                                   pred.high=cellulose_mass.jack.stat[,3],
                                   Measured=meta(ref.test.bn)$cellulose_mass,
                                   ncomp=ncomp_cellulose_mass_CVmodel,
                                   Species=meta(ref.test.bn)$species,
                                   Project=meta(ref.test.bn)$project,
                                   functional.group=meta(ref.test.bn)$functional.group,
                                   ID=meta(ref.test.bn)$sample_id)
cellulose_all<-with(cellulose_mass.jack.df,c(pred.low,pred.high,Measured))
cellulose_upper<-max(cellulose_all,na.rm=T)+1
cellulose_lower<-min(cellulose_all,na.rm=T)-1

cellulose_mass.val.plot<-ggplot(cellulose_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),
                  ylim=c(cellulose_lower,cellulose_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

lignin_mass.jack.pred<-apply.coefs(lignin_mass.jack.coefs,as.matrix(ref.test.bn))
lignin_mass.jack.stat<-t(apply(lignin_mass.jack.pred,1,
                               function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_mass.jack.df<-data.frame(pred.mean=lignin_mass.jack.stat[,1],
                                pred.low=lignin_mass.jack.stat[,2],
                                pred.high=lignin_mass.jack.stat[,3],
                                Measured=meta(ref.test.bn)$lignin_mass,
                                ncomp=ncomp_lignin_mass_CVmodel,
                                Species=meta(ref.test.bn)$species,
                                Project=meta(ref.test.bn)$project,
                                functional.group=meta(ref.test.bn)$functional.group,
                                ID=meta(ref.test.bn)$sample_id)
lignin_all<-with(lignin_mass.jack.df,c(pred.low,pred.high,Measured))
lignin_upper<-max(lignin_all,na.rm=T)+0.5
lignin_lower<-min(lignin_all,na.rm=T)-0.5

lignin_mass.val.plot<-ggplot(lignin_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),
                  ylim=c(lignin_lower,lignin_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_mass.jack.pred<-apply.coefs(chlA_mass.jack.coefs,as.matrix(ref.test.bn))
chlA_mass.jack.stat<-t(apply(chlA_mass.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_mass.jack.df<-data.frame(pred.mean=chlA_mass.jack.stat[,1],
                              pred.low=chlA_mass.jack.stat[,2],
                              pred.high=chlA_mass.jack.stat[,3],
                              Measured=meta(ref.test.bn)$chlA_mass,
                              ncomp=ncomp_chlA_mass_CVmodel,
                              Species=meta(ref.test.bn)$species,
                              Project=meta(ref.test.bn)$project,
                              functional.group=meta(ref.test.bn)$functional.group,
                              ID=meta(ref.test.bn)$sample_id)
chlA_all<-with(chlA_mass.jack.df,c(pred.low,pred.high,Measured))
chlA_upper<-max(chlA_all,na.rm=T)+1
chlA_lower<-min(chlA_all,na.rm=T)-1

chlA_mass.val.plot<-ggplot(chlA_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),
                  ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlB_mass.jack.pred<-apply.coefs(chlB_mass.jack.coefs,as.matrix(ref.test.bn))
chlB_mass.jack.stat<-t(apply(chlB_mass.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_mass.jack.df<-data.frame(pred.mean=chlB_mass.jack.stat[,1],
                              pred.low=chlB_mass.jack.stat[,2],
                              pred.high=chlB_mass.jack.stat[,3],
                              Measured=meta(ref.test.bn)$chlB_mass,
                              ncomp=ncomp_chlB_mass_CVmodel,
                              Species=meta(ref.test.bn)$species,
                              Project=meta(ref.test.bn)$project,
                              functional.group=meta(ref.test.bn)$functional.group,
                              ID=meta(ref.test.bn)$sample_id)
chlB_all<-with(chlB_mass.jack.df,c(pred.low,pred.high,Measured))
chlB_upper<-max(chlB_all,na.rm=T)+0.5
chlB_lower<-min(chlB_all,na.rm=T)-0.5

chlB_mass.val.plot<-ggplot(chlB_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),
                  ylim=c(chlB_lower,chlB_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

car_mass.jack.pred<-apply.coefs(car_mass.jack.coefs,as.matrix(ref.test.bn))
car_mass.jack.stat<-t(apply(car_mass.jack.pred,1,
                            function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_mass.jack.df<-data.frame(pred.mean=car_mass.jack.stat[,1],
                             pred.low=car_mass.jack.stat[,2],
                             pred.high=car_mass.jack.stat[,3],
                             Measured=meta(ref.test.bn)$car_mass,
                             ncomp=ncomp_car_mass_CVmodel,
                             Species=meta(ref.test.bn)$species,
                             Project=meta(ref.test.bn)$project,
                             functional.group=meta(ref.test.bn)$functional.group,
                             ID=meta(ref.test.bn)$sample_id)
car_all<-with(car_mass.jack.df,c(pred.low,pred.high,Measured))
car_upper<-max(car_all,na.rm=T)+0.3
car_lower<-min(car_all,na.rm=T)-0.3

car_mass.val.plot<-ggplot(car_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_lower,car_upper),
                  ylim=c(car_lower,car_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

Cmass.jack.pred<-apply.coefs(Cmass.jack.coefs,as.matrix(ref.test.bn))
Cmass.jack.stat<-t(apply(Cmass.jack.pred,1,
                         function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cmass.jack.df<-data.frame(pred.mean=Cmass.jack.stat[,1],
                          pred.low=Cmass.jack.stat[,2],
                          pred.high=Cmass.jack.stat[,3],
                          Measured=meta(ref.test.bn)$Cmass,
                          ncomp=ncomp_Cmass_CVmodel,
                          Species=meta(ref.test.bn)$species,
                          Project=meta(ref.test.bn)$project,
                          functional.group=meta(ref.test.bn)$functional.group,
                          ID=meta(ref.test.bn)$sample_id)
C_all<-with(Cmass.jack.df,c(pred.low,pred.high,Measured))
C_upper<-max(C_all,na.rm=T)+1
C_lower<-min(C_all,na.rm=T)-1

Cmass.val.plot<-ggplot(Cmass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(C_lower,C_upper),
                  ylim=c(C_lower,C_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Nmass.jack.pred<-apply.coefs(Nmass.jack.coefs,as.matrix(ref.test.bn))
Nmass.jack.stat<-t(apply(Nmass.jack.pred,1,
                         function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Nmass.jack.df<-data.frame(pred.mean=Nmass.jack.stat[,1],
                          pred.low=Nmass.jack.stat[,2],
                          pred.high=Nmass.jack.stat[,3],
                          Measured=meta(ref.test.bn)$Nmass,
                          ncomp=ncomp_Nmass_CVmodel,
                          Species=meta(ref.test.bn)$species,
                          Project=meta(ref.test.bn)$project,
                          functional.group=meta(ref.test.bn)$functional.group,
                          ID=meta(ref.test.bn)$sample_id)
N_all<-with(Nmass.jack.df,c(pred.low,pred.high,Measured))
N_upper<-max(N_all,na.rm=T)+0.1
N_lower<-min(N_all,na.rm=T)-0.1

Nmass.val.plot<-ggplot(Nmass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(N_lower,N_upper),
                  ylim=c(N_lower,N_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

EWT.jack.pred<-apply.coefs(EWT.jack.coefs,as.matrix(ref.test.bn))
EWT.jack.stat<-t(apply(EWT.jack.pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT.jack.df<-data.frame(pred.mean=EWT.jack.stat[,1],
                        pred.low=EWT.jack.stat[,2],
                        pred.high=EWT.jack.stat[,3],
                        Measured=meta(ref.test.bn)$EWT,
                        ncomp=ncomp_EWT_CVmodel,
                        Species=meta(ref.test.bn)$species,
                        Project=meta(ref.test.bn)$project,
                        functional.group=meta(ref.test.bn)$functional.group,
                        ID=meta(ref.test.bn)$sample_id)
EWT_all<-with(EWT.jack.df,c(pred.low,pred.high,Measured))
EWT_upper<-max(EWT_all,na.rm=T)+0.03
EWT_lower<-min(EWT_all,na.rm=T)-0.03

EWT.val.plot<-ggplot(EWT.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

LMA.jack.pred<-apply.coefs(LMA.jack.coefs,as.matrix(ref.test.bn))
LMA.jack.stat<-t(apply(LMA.jack.pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA.jack.df<-data.frame(pred.mean=LMA.jack.stat[,1],
                        pred.low=LMA.jack.stat[,2],
                        pred.high=LMA.jack.stat[,3],
                        Measured=meta(ref.test.bn)$LMA,
                        ncomp=ncomp_LMA_CVmodel,
                        Species=meta(ref.test.bn)$species,
                        Project=meta(ref.test.bn)$project,
                        functional.group=meta(ref.test.bn)$functional.group,
                        ID=meta(ref.test.bn)$sample_id)
LMA_all<-with(LMA.jack.df,c(pred.low,pred.high,Measured))
LMA_upper<-max(LMA_all,na.rm=T)+0.02
LMA_lower<-min(LMA_all,na.rm=T)-0.02

LMA.val.plot<-ggplot(LMA.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),
                  ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

LDMC.jack.pred<-apply.coefs(LDMC.jack.coefs,as.matrix(ref.test.bn))
LDMC.jack.stat<-t(apply(LDMC.jack.pred,1,
                        function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LDMC.jack.df<-data.frame(pred.mean=LDMC.jack.stat[,1],
                         pred.low=LDMC.jack.stat[,2],
                         pred.high=LDMC.jack.stat[,3],
                         Measured=meta(ref.test.bn)$LDMC,
                         ncomp=ncomp_LDMC_CVmodel,
                         Species=meta(ref.test.bn)$species,
                         Project=meta(ref.test.bn)$project,
                         functional.group=meta(ref.test.bn)$functional.group,
                         ID=meta(ref.test.bn)$sample_id)
LDMC_all<-with(LDMC.jack.df,c(pred.low,pred.high,Measured))
LDMC_upper<-max(LDMC_all,na.rm=T)+10
LDMC_lower<-min(LDMC_all,na.rm=T)-10

LDMC.val.plot<-ggplot(LDMC.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),
                  ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Al_mass.jack.pred<-apply.coefs(Al_mass.jack.coefs,as.matrix(ref.test.bn))
Al_mass.jack.stat<-t(apply(Al_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_mass.jack.df<-data.frame(pred.mean=Al_mass.jack.stat[,1],
                            pred.low=Al_mass.jack.stat[,2],
                            pred.high=Al_mass.jack.stat[,3],
                            Measured=meta(ref.test.bn)$Al_mass,
                            ncomp=ncomp_Al_mass_CVmodel,
                            Species=meta(ref.test.bn)$species,
                            Project=meta(ref.test.bn)$project,
                            functional.group=meta(ref.test.bn)$functional.group,
                            ID=meta(ref.test.bn)$sample_id)
Al_all<-with(Al_mass.jack.df,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Al_upper<-max(Al_all,na.rm=T)+0.02
Al_lower<-min(Al_all,na.rm=T)-0.02

Al_mass.val.plot<-ggplot(Al_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Al_lower,Al_upper),
                  ylim=c(Al_lower,Al_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Al (mg g"^-1*")"),x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Ca_mass.jack.pred<-apply.coefs(Ca_mass.jack.coefs,as.matrix(ref.test.bn))
Ca_mass.jack.stat<-t(apply(Ca_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_mass.jack.df<-data.frame(pred.mean=Ca_mass.jack.stat[,1],
                            pred.low=Ca_mass.jack.stat[,2],
                            pred.high=Ca_mass.jack.stat[,3],
                            Measured=meta(ref.test.bn)$Ca_mass,
                            ncomp=ncomp_Ca_mass_CVmodel,
                            Species=meta(ref.test.bn)$species,
                            Project=meta(ref.test.bn)$project,
                            functional.group=meta(ref.test.bn)$functional.group,
                            ID=meta(ref.test.bn)$sample_id)
Ca_all<-with(Ca_mass.jack.df,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Ca_upper<-max(Ca_all,na.rm=T)+1
Ca_lower<-min(Ca_all,na.rm=T)-1

Ca_mass.val.plot<-ggplot(Ca_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),
                  ylim=c(Ca_lower,Ca_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Ca (mg g"^-1*")"),
       x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cu_mass.jack.pred<-apply.coefs(Cu_mass.jack.coefs,as.matrix(ref.test.bn))
Cu_mass.jack.stat<-t(apply(Cu_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_mass.jack.df<-data.frame(pred.mean=Cu_mass.jack.stat[,1],
                            pred.low=Cu_mass.jack.stat[,2],
                            pred.high=Cu_mass.jack.stat[,3],
                            Measured=meta(ref.test.bn)$Cu_mass,
                            ncomp=ncomp_Cu_mass_CVmodel,
                            Species=meta(ref.test.bn)$species,
                            Project=meta(ref.test.bn)$project,
                            functional.group=meta(ref.test.bn)$functional.group,
                            ID=meta(ref.test.bn)$sample_id)
Cu_all<-with(Cu_mass.jack.df,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Cu_upper<-max(Cu_all,na.rm=T)+0.002
Cu_lower<-min(Cu_all,na.rm=T)-0.002

Cu_mass.val.plot<-ggplot(Cu_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),
                  ylim=c(Cu_lower,Cu_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Cu (mg g"^-1*")"),
       x=expression("Predicted Cu (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

Fe_mass.jack.pred<-apply.coefs(Fe_mass.jack.coefs,as.matrix(ref.test.bn))
Fe_mass.jack.stat<-t(apply(Fe_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_mass.jack.df<-data.frame(pred.mean=Fe_mass.jack.stat[,1],
                            pred.low=Fe_mass.jack.stat[,2],
                            pred.high=Fe_mass.jack.stat[,3],
                            Measured=meta(ref.test.bn)$Fe_mass,
                            ncomp=ncomp_Fe_mass_CVmodel,
                            Species=meta(ref.test.bn)$species,
                            Project=meta(ref.test.bn)$project,
                            functional.group=meta(ref.test.bn)$functional.group,
                            ID=meta(ref.test.bn)$sample_id)
Fe_all<-with(Fe_mass.jack.df,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Fe_upper<-max(Fe_all,na.rm=T)+0.01
Fe_lower<-min(Fe_all,na.rm=T)-0.01

Fe_mass.val.plot<-ggplot(Fe_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),
                  ylim=c(Fe_lower,Fe_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Fe (mg g"^-1*")"),
       x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

K_mass.jack.pred<-apply.coefs(K_mass.jack.coefs,as.matrix(ref.test.bn))
K_mass.jack.stat<-t(apply(K_mass.jack.pred,1,
                          function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_mass.jack.df<-data.frame(pred.mean=K_mass.jack.stat[,1],
                           pred.low=K_mass.jack.stat[,2],
                           pred.high=K_mass.jack.stat[,3],
                           Measured=meta(ref.test.bn)$K_mass,
                           ncomp=ncomp_K_mass_CVmodel,
                           Species=meta(ref.test.bn)$species,
                           Project=meta(ref.test.bn)$project,
                           functional.group=meta(ref.test.bn)$functional.group,
                           ID=meta(ref.test.bn)$sample_id)
K_all<-with(K_mass.jack.df,c(pred.low[!is.na(Measured)],
                             pred.high[!is.na(Measured)],
                             Measured))
K_upper<-max(K_all,na.rm=T)+1
K_lower<-min(K_all,na.rm=T)-1

K_mass.val.plot<-ggplot(K_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),
                  ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mg_mass.jack.pred<-apply.coefs(Mg_mass.jack.coefs,as.matrix(ref.test.bn))
Mg_mass.jack.stat<-t(apply(Mg_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_mass.jack.df<-data.frame(pred.mean=Mg_mass.jack.stat[,1],
                            pred.low=Mg_mass.jack.stat[,2],
                            pred.high=Mg_mass.jack.stat[,3],
                            Measured=meta(ref.test.bn)$Mg_mass,
                            ncomp=ncomp_Mg_mass_CVmodel,
                            Species=meta(ref.test.bn)$species,
                            Project=meta(ref.test.bn)$project,
                            functional.group=meta(ref.test.bn)$functional.group,
                            ID=meta(ref.test.bn)$sample_id)
Mg_all<-with(Mg_mass.jack.df,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Mg_upper<-max(Mg_all,na.rm=T)+0.5
Mg_lower<-min(Mg_all,na.rm=T)-0.5

Mg_mass.val.plot<-ggplot(Mg_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),
                  ylim=c(Mg_lower,Mg_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Mg (mg g"^-1*")"),
       x=expression("Predicted Mg (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mn_mass.jack.pred<-apply.coefs(Mn_mass.jack.coefs,as.matrix(ref.test.bn))
Mn_mass.jack.stat<-t(apply(Mn_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_mass.jack.df<-data.frame(pred.mean=Mn_mass.jack.stat[,1],
                            pred.low=Mn_mass.jack.stat[,2],
                            pred.high=Mn_mass.jack.stat[,3],
                            Measured=meta(ref.test.bn)$Mn_mass,
                            ncomp=ncomp_Mn_mass_CVmodel,
                            Species=meta(ref.test.bn)$species,
                            Project=meta(ref.test.bn)$project,
                            functional.group=meta(ref.test.bn)$functional.group,
                            ID=meta(ref.test.bn)$sample_id)
Mn_all<-with(Mn_mass.jack.df,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Mn_upper<-max(Mn_all,na.rm=T)+0.05
Mn_lower<-min(Mn_all,na.rm=T)-0.05

Mn_mass.val.plot<-ggplot(Mn_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),
                  ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Mn (mg g"^-1*")"),
       x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Na_mass.jack.pred<-apply.coefs(Na_mass.jack.coefs,as.matrix(ref.test.bn))
Na_mass.jack.stat<-t(apply(Na_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_mass.jack.df<-data.frame(pred.mean=Na_mass.jack.stat[,1],
                            pred.low=Na_mass.jack.stat[,2],
                            pred.high=Na_mass.jack.stat[,3],
                            Measured=meta(ref.test.bn)$Na_mass,
                            ncomp=ncomp_Na_mass_CVmodel,
                            Species=meta(ref.test.bn)$species,
                            Project=meta(ref.test.bn)$project,
                            functional.group=meta(ref.test.bn)$functional.group,
                            ID=meta(ref.test.bn)$sample_id)
Na_all<-with(Na_mass.jack.df,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Na_upper<-max(Na_all,na.rm=T)+0.2
Na_lower<-min(Na_all,na.rm=T)-0.2

Na_mass.val.plot<-ggplot(Na_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Na_lower,Na_upper),
                  ylim=c(Na_lower,Na_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Na (mg g"^-1*")"),
       x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

P_mass.jack.pred<-apply.coefs(P_mass.jack.coefs,as.matrix(ref.test.bn))
P_mass.jack.stat<-t(apply(P_mass.jack.pred,1,
                          function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_mass.jack.df<-data.frame(pred.mean=P_mass.jack.stat[,1],
                           pred.low=P_mass.jack.stat[,2],
                           pred.high=P_mass.jack.stat[,3],
                           Measured=meta(ref.test.bn)$P_mass,
                           ncomp=ncomp_P_mass_CVmodel,
                           Species=meta(ref.test.bn)$species,
                           Project=meta(ref.test.bn)$project,
                           functional.group=meta(ref.test.bn)$functional.group,
                           ID=meta(ref.test.bn)$sample_id)
P_all<-with(P_mass.jack.df,c(pred.low[!is.na(Measured)],
                             pred.high[!is.na(Measured)],
                             Measured))
P_upper<-max(P_all,na.rm=T)+0.5
P_lower<-min(P_all,na.rm=T)-0.5

P_mass.val.plot<-ggplot(P_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(P_lower,P_upper),
                  ylim=c(P_lower,P_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured P (mg g"^-1*")"),
       x=expression("Predicted P (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Zn_mass.jack.pred<-apply.coefs(Zn_mass.jack.coefs,as.matrix(ref.test.bn))
Zn_mass.jack.stat<-t(apply(Zn_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_mass.jack.df<-data.frame(pred.mean=Zn_mass.jack.stat[,1],
                            pred.low=Zn_mass.jack.stat[,2],
                            pred.high=Zn_mass.jack.stat[,3],
                            Measured=meta(ref.test.bn)$Zn_mass,
                            ncomp=ncomp_Zn_mass_CVmodel,
                            Species=meta(ref.test.bn)$species,
                            Project=meta(ref.test.bn)$project,
                            functional.group=meta(ref.test.bn)$functional.group,
                            ID=meta(ref.test.bn)$sample_id)
Zn_all<-with(Zn_mass.jack.df,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Zn_upper<-max(Zn_all,na.rm=T)+0.03
Zn_lower<-min(Zn_all,na.rm=T)-0.03

Zn_mass.val.plot<-ggplot(Zn_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),
                  ylim=c(Zn_lower,Zn_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Zn (mg g"^-1*")"),
       x=expression("Predicted Zn (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

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
saveRDS(all.jack.coef.list,"SavedResults/all_jack_coefs_list_ref_bn.rds")

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
saveRDS(all.jack.df.list,"SavedResults/all_jack_df_list_ref_bn.rds")

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
saveRDS(all.jack.stats.list,"SavedResults/all_jack_stats_list_ref_bn.rds")