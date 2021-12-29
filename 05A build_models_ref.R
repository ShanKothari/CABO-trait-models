setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(pls)
library(reshape2)
library(RColorBrewer)
library(patchwork)
source("Scripts/VIP.R")

ref.train<-readRDS("ProcessedSpectra/ref_train.rds")
ref.test<-readRDS("ProcessedSpectra/ref_test.rds")

##########################################
## to dos

## try Type II regression?
## eliminate outliers for Cu and Fe
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
  range<-quantile(measured,probs=max,na.rm=na.rm)-quantile(measured,probs=min,na.rm=na.rm)
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

solubles_mass_CVmodel<-plsr(meta(ref.train)$solubles_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_mass_CVmodel <- selectNcomp(solubles_mass_CVmodel, method = "onesigma", plot = FALSE)
solubles_mass_valid <- which(!is.na(meta(ref.train)$solubles_mass))
solubles_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[solubles_mass_valid],
                         Species=meta(ref.train)$species[solubles_mass_valid],
                         Project=meta(ref.train)$project[solubles_mass_valid],
                         measured=meta(ref.train)$solubles_mass[solubles_mass_valid],
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

hemicellulose_mass_CVmodel<-plsr(meta(ref.train)$hemicellulose_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_mass_CVmodel <- selectNcomp(hemicellulose_mass_CVmodel, method = "onesigma", plot = FALSE)
hemicellulose_mass_valid <- which(!is.na(meta(ref.train)$hemicellulose_mass))
hemicellulose_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[hemicellulose_mass_valid],
                         Species=meta(ref.train)$species[hemicellulose_mass_valid],
                         Project=meta(ref.train)$project[hemicellulose_mass_valid],
                         measured=meta(ref.train)$hemicellulose_mass[hemicellulose_mass_valid],
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

cellulose_mass_CVmodel<-plsr(meta(ref.train)$cellulose_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_mass_CVmodel <- selectNcomp(cellulose_mass_CVmodel, method = "onesigma", plot = FALSE)
cellulose_mass_valid <- which(!is.na(meta(ref.train)$cellulose_mass))
cellulose_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[cellulose_mass_valid],
                         Species=meta(ref.train)$species[cellulose_mass_valid],
                         Project=meta(ref.train)$project[cellulose_mass_valid],
                         measured=meta(ref.train)$cellulose_mass[cellulose_mass_valid],
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

lignin_mass_CVmodel<-plsr(meta(ref.train)$lignin_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_mass_CVmodel <- selectNcomp(lignin_mass_CVmodel, method = "onesigma", plot = FALSE)
lignin_mass_valid <- which(!is.na(meta(ref.train)$lignin_mass))
lignin_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[lignin_mass_valid],
                         Species=meta(ref.train)$species[lignin_mass_valid],
                         Project=meta(ref.train)$project[lignin_mass_valid],
                         measured=meta(ref.train)$lignin_mass[lignin_mass_valid],
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

Cmass_CVmodel<-plsr(meta(ref.train)$Cmass~as.matrix(ref.train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cmass_CVmodel <- selectNcomp(Cmass_CVmodel, method = "onesigma", plot = FALSE)
Cmass_valid <- which(!is.na(meta(ref.train)$Cmass))
Cmass_pred<-data.frame(ID=meta(ref.train)$sample_id[Cmass_valid],
                            Species=meta(ref.train)$species[Cmass_valid],
                            Project=meta(ref.train)$project[Cmass_valid],
                            measured=meta(ref.train)$Cmass[Cmass_valid],
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

Nmass_CVmodel<-plsr(meta(ref.train)$Nmass~as.matrix(ref.train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Nmass_CVmodel <- selectNcomp(Nmass_CVmodel, method = "onesigma", plot = FALSE)
Nmass_valid <- which(!is.na(meta(ref.train)$Nmass))
Nmass_pred<-data.frame(ID=meta(ref.train)$sample_id[Nmass_valid],
                      Species=meta(ref.train)$species[Nmass_valid],
                      Project=meta(ref.train)$project[Nmass_valid],
                      measured=meta(ref.train)$Nmass[Nmass_valid],
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

EWT_CVmodel<-plsr(meta(ref.train)$EWT~as.matrix(ref.train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_CVmodel <- selectNcomp(EWT_CVmodel, method = "onesigma", plot = FALSE)
EWT_valid <- which(!is.na(meta(ref.train)$EWT))
EWT_pred<-data.frame(ID=meta(ref.train)$sample_id[EWT_valid],
                     Species=meta(ref.train)$species[EWT_valid],
                     Project=meta(ref.train)$project[EWT_valid],
                     measured=meta(ref.train)$EWT[EWT_valid],
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

LDMC_CVmodel<-plsr(meta(ref.train)$LDMC~as.matrix(ref.train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LDMC_CVmodel <- selectNcomp(LDMC_CVmodel, method = "onesigma", plot = FALSE)
LDMC_valid <- which(!is.na(meta(ref.train)$LDMC))
LDMC_pred<-data.frame(ID=meta(ref.train)$sample_id[LDMC_valid],
                      Species=meta(ref.train)$species[LDMC_valid],
                      Project=meta(ref.train)$project[LDMC_valid],
                      measured=meta(ref.train)$LDMC[LDMC_valid],
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

LMA_CVmodel<-plsr(meta(ref.train)$LMA~as.matrix(ref.train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_CVmodel <- selectNcomp(LMA_CVmodel, method = "onesigma", plot = FALSE)
LMA_valid <- which(!is.na(meta(ref.train)$LMA))
LMA_pred<-data.frame(ID=meta(ref.train)$sample_id[LMA_valid],
                      Species=meta(ref.train)$species[LMA_valid],
                      Project=meta(ref.train)$project[LMA_valid],
                      measured=meta(ref.train)$LMA[LMA_valid],
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

chlA_mass_CVmodel<-plsr(meta(ref.train)$chlA_mass~as.matrix(ref.train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_mass_CVmodel <- selectNcomp(chlA_mass_CVmodel, method = "onesigma", plot = FALSE)
chlA_mass_valid <- which(!is.na(meta(ref.train)$chlA_mass))
chlA_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[chlA_mass_valid],
                     Species=meta(ref.train)$species[chlA_mass_valid],
                     Project=meta(ref.train)$project[chlA_mass_valid],
                     measured=meta(ref.train)$chlA_mass[chlA_mass_valid],
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

chlB_mass_CVmodel<-plsr(meta(ref.train)$chlB_mass~as.matrix(ref.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_mass_CVmodel <- selectNcomp(chlB_mass_CVmodel, method = "onesigma", plot = FALSE)
chlB_mass_valid <- which(!is.na(meta(ref.train)$chlB_mass))
chlB_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[chlB_mass_valid],
                            Species=meta(ref.train)$species[chlB_mass_valid],
                            Project=meta(ref.train)$project[chlB_mass_valid],
                            measured=meta(ref.train)$chlB_mass[chlB_mass_valid],
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

car_mass_CVmodel<-plsr(meta(ref.train)$car_mass~as.matrix(ref.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_mass_CVmodel <- selectNcomp(car_mass_CVmodel, method = "onesigma", plot = FALSE)
car_mass_valid <- which(!is.na(meta(ref.train)$car_mass))
car_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[car_mass_valid],
                            Species=meta(ref.train)$species[car_mass_valid],
                            Project=meta(ref.train)$project[car_mass_valid],
                            measured=meta(ref.train)$car_mass[car_mass_valid],
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

Al_mass_CVmodel<-plsr(meta(ref.train)$Al_mass~as.matrix(ref.train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Al_mass_CVmodel <- selectNcomp(Al_mass_CVmodel, method = "onesigma", plot = FALSE)
Al_mass_valid <- which(!is.na(meta(ref.train)$Al_mass))
Al_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[Al_mass_valid],
                          Species=meta(ref.train)$species[Al_mass_valid],
                          Project=meta(ref.train)$project[Al_mass_valid],
                          measured=meta(ref.train)$Al_mass[Al_mass_valid],
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

Ca_mass_CVmodel<-plsr(meta(ref.train)$Ca_mass~as.matrix(ref.train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Ca_mass_CVmodel <- selectNcomp(Ca_mass_CVmodel, method = "onesigma", plot = FALSE)
Ca_mass_valid <- which(!is.na(meta(ref.train)$Ca_mass))
Ca_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[Ca_mass_valid],
                          Species=meta(ref.train)$species[Ca_mass_valid],
                          Project=meta(ref.train)$project[Ca_mass_valid],
                          measured=meta(ref.train)$Ca_mass[Ca_mass_valid],
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

Cu_mass_CVmodel<-plsr(meta(ref.train)$Cu_mass~as.matrix(ref.train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cu_mass_CVmodel <- selectNcomp(Cu_mass_CVmodel, method = "onesigma", plot = FALSE)
Cu_mass_valid <- which(!is.na(meta(ref.train)$Cu_mass))
Cu_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[Cu_mass_valid],
                          Species=meta(ref.train)$species[Cu_mass_valid],
                          Project=meta(ref.train)$project[Cu_mass_valid],
                          measured=meta(ref.train)$Cu_mass[Cu_mass_valid],
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

Fe_mass_CVmodel<-plsr(meta(ref.train)$Fe_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Fe_mass_CVmodel <- selectNcomp(Fe_mass_CVmodel, method = "onesigma", plot = FALSE)
Fe_mass_valid <- which(!is.na(meta(ref.train)$Fe_mass))
Fe_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[Fe_mass_valid],
                         Species=meta(ref.train)$species[Fe_mass_valid],
                         Project=meta(ref.train)$project[Fe_mass_valid],
                         measured=meta(ref.train)$Fe_mass[Fe_mass_valid],
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

K_mass_CVmodel<-plsr(meta(ref.train)$K_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_mass_CVmodel <- selectNcomp(K_mass_CVmodel, method = "onesigma", plot = FALSE)
K_mass_valid <- which(!is.na(meta(ref.train)$K_mass))
K_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[K_mass_valid],
                         Species=meta(ref.train)$species[K_mass_valid],
                         Project=meta(ref.train)$project[K_mass_valid],
                         measured=meta(ref.train)$K_mass[K_mass_valid],
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

Mg_mass_CVmodel<-plsr(meta(ref.train)$Mg_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mg_mass_CVmodel <- selectNcomp(Mg_mass_CVmodel, method = "onesigma", plot = FALSE)
Mg_mass_valid <- which(!is.na(meta(ref.train)$Mg_mass))
Mg_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[Mg_mass_valid],
                         Species=meta(ref.train)$species[Mg_mass_valid],
                         Project=meta(ref.train)$project[Mg_mass_valid],
                         measured=meta(ref.train)$Mg_mass[Mg_mass_valid],
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

Mn_mass_CVmodel<-plsr(meta(ref.train)$Mn_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Mn_mass_CVmodel <- selectNcomp(Mn_mass_CVmodel, method = "onesigma", plot = FALSE)
Mn_mass_valid <- which(!is.na(meta(ref.train)$Mn_mass))
Mn_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[Mn_mass_valid],
                         Species=meta(ref.train)$species[Mn_mass_valid],
                         Project=meta(ref.train)$project[Mn_mass_valid],
                         measured=meta(ref.train)$Mn_mass[Mn_mass_valid],
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

Na_mass_CVmodel<-plsr(meta(ref.train)$Na_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Na_mass_CVmodel <- selectNcomp(Na_mass_CVmodel, method = "onesigma", plot = FALSE)
Na_mass_valid <- which(!is.na(meta(ref.train)$Na_mass))
Na_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[Na_mass_valid],
                         Species=meta(ref.train)$species[Na_mass_valid],
                         Project=meta(ref.train)$project[Na_mass_valid],
                         measured=meta(ref.train)$Na_mass[Na_mass_valid],
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

P_mass_CVmodel<-plsr(meta(ref.train)$P_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_P_mass_CVmodel <- selectNcomp(P_mass_CVmodel, method = "onesigma", plot = FALSE)
P_mass_valid <- which(!is.na(meta(ref.train)$P_mass))
P_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[P_mass_valid],
                         Species=meta(ref.train)$species[P_mass_valid],
                         Project=meta(ref.train)$project[P_mass_valid],
                         measured=meta(ref.train)$P_mass[P_mass_valid],
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

Zn_mass_CVmodel<-plsr(meta(ref.train)$Zn_mass~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Zn_mass_CVmodel <- selectNcomp(Zn_mass_CVmodel, method = "onesigma", plot = FALSE)
Zn_mass_valid <- which(!is.na(meta(ref.train)$Zn_mass))
Zn_mass_pred<-data.frame(ID=meta(ref.train)$sample_id[Zn_mass_valid],
                         Species=meta(ref.train)$species[Zn_mass_valid],
                         Project=meta(ref.train)$project[Zn_mass_valid],
                         measured=meta(ref.train)$Zn_mass[Zn_mass_valid],
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

VIP2.df<-data.frame(solubles_mass=VIP(solubles_mass_CVmodel)[ncomp_solubles_mass_CVmodel,],
                    hemicellulose_mass=VIP(hemicellulose_mass_CVmodel)[ncomp_hemicellulose_mass_CVmodel,],
                    cellulose_mass=VIP(cellulose_mass_CVmodel)[ncomp_cellulose_mass_CVmodel,],
                    lignin_mass=VIP(lignin_mass_CVmodel)[ncomp_lignin_mass_CVmodel,],
                    wavelength=400:2400)

VIP3.df<-data.frame(chlA_mass=VIP(chlA_mass_CVmodel)[ncomp_chlA_mass_CVmodel,],
                    chlB_mass=VIP(chlB_mass_CVmodel)[ncomp_chlB_mass_CVmodel,],
                    car_mass=VIP(car_mass_CVmodel)[ncomp_car_mass_CVmodel,],
                    wavelength=400:2400)

VIP4.df<-data.frame(Al_mass=VIP(Al_mass_CVmodel)[ncomp_Al_mass_CVmodel,],
                    Ca_mass=VIP(Ca_mass_CVmodel)[ncomp_Ca_mass_CVmodel,],
                    Cu_mass=VIP(Cu_mass_CVmodel)[ncomp_Cu_mass_CVmodel,],
                    Fe_mass=VIP(Fe_mass_CVmodel)[ncomp_Fe_mass_CVmodel,],
                    K_mass=VIP(K_mass_CVmodel)[ncomp_K_mass_CVmodel,],
                    wavelength=400:2400)

VIP5.df<-data.frame(Mg_mass=VIP(Mg_mass_CVmodel)[ncomp_Mg_mass_CVmodel,],
                    Mn_mass=VIP(Mn_mass_CVmodel)[ncomp_Mn_mass_CVmodel,],
                    Na_mass=VIP(Na_mass_CVmodel)[ncomp_Na_mass_CVmodel,],
                    P_mass=VIP(P_mass_CVmodel)[ncomp_P_mass_CVmodel,],
                    Zn_mass=VIP(Zn_mass_CVmodel)[ncomp_Zn_mass_CVmodel,],
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
  scale_color_manual(labels=c(expression("C"[mass]),
                              expression("N"[mass]),
                              "LMA","LDMC","EWT"),
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

# pdf("Images/VIP_ref.pdf",height=12,width=6,onefile = F)
# VIP1.plot / VIP2.plot / VIP3.plot / VIP4.plot / VIP5.plot
# dev.off()

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
  
  n.cal.spec<-nrow(ref.train)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-ref.train[train.jack]
  val.jack<-ref.train[test.jack]
  
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
                                   percentRMSE=percentRMSD(meta(val.jack)$solubles_mass,solubles_mass.jack.val.pred,0.025,0.975),
                                   bias=mean(solubles_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$solubles_mass,na.rm=T))
  
  hemicellulose_mass.jack.val.pred<-as.vector(predict(hemicellulose_mass.jack,newdata=as.matrix(val.jack),
                                                      ncomp=ncomp_hemicellulose_mass_CVmodel)[,,1])
  hemicellulose_mass.jack.val.fit<-lm(hemicellulose_mass.jack.val.pred~meta(val.jack)$hemicellulose_mass)
  hemicellulose_mass.jack.stats[[i]]<-c(R2=summary(hemicellulose_mass.jack.val.fit)$r.squared,
                                        RMSE=RMSD(meta(val.jack)$hemicellulose_mass,hemicellulose_mass.jack.val.pred),
                                        percentRMSE=percentRMSD(meta(val.jack)$hemicellulose_mass,hemicellulose_mass.jack.val.pred,0.025,0.975),
                                        bias=mean(hemicellulose_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$hemicellulose_mass,na.rm=T))
  
  cellulose_mass.jack.val.pred<-as.vector(predict(cellulose_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_cellulose_mass_CVmodel)[,,1])
  cellulose_mass.jack.val.fit<-lm(cellulose_mass.jack.val.pred~meta(val.jack)$cellulose_mass)
  cellulose_mass.jack.stats[[i]]<-c(R2=summary(cellulose_mass.jack.val.fit)$r.squared,
                                    RMSE=RMSD(meta(val.jack)$cellulose_mass,cellulose_mass.jack.val.pred),
                                    percentRMSE=percentRMSD(meta(val.jack)$cellulose_mass,cellulose_mass.jack.val.pred,0.025,0.975),
                                    bias=mean(cellulose_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$cellulose_mass,na.rm=T))
  
  lignin_mass.jack.val.pred<-as.vector(predict(lignin_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_lignin_mass_CVmodel)[,,1])
  lignin_mass.jack.val.fit<-lm(lignin_mass.jack.val.pred~meta(val.jack)$lignin_mass)
  lignin_mass.jack.stats[[i]]<-c(R2=summary(lignin_mass.jack.val.fit)$r.squared,
                                 RMSE=RMSD(meta(val.jack)$lignin_mass,lignin_mass.jack.val.pred),
                                 percentRMSE=percentRMSD(meta(val.jack)$lignin_mass,lignin_mass.jack.val.pred,0.025,0.975),
                                 bias=mean(lignin_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$lignin_mass,na.rm=T))
  
  chlA_mass.jack.val.pred<-as.vector(predict(chlA_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlA_mass_CVmodel)[,,1])
  chlA_mass.jack.val.fit<-lm(chlA_mass.jack.val.pred~meta(val.jack)$chlA_mass)
  chlA_mass.jack.stats[[i]]<-c(R2=summary(chlA_mass.jack.val.fit)$r.squared,
                               RMSE=RMSD(meta(val.jack)$chlA_mass,chlA_mass.jack.val.pred),
                               percentRMSE=percentRMSD(meta(val.jack)$chlA_mass,chlA_mass.jack.val.pred,0.025,0.975),
                               bias=mean(chlA_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlA_mass,na.rm=T))
  
  chlB_mass.jack.val.pred<-as.vector(predict(chlB_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlB_mass_CVmodel)[,,1])
  chlB_mass.jack.val.fit<-lm(chlB_mass.jack.val.pred~meta(val.jack)$chlB_mass)
  chlB_mass.jack.stats[[i]]<-c(R2=summary(chlB_mass.jack.val.fit)$r.squared,
                               RMSE=RMSD(meta(val.jack)$chlB_mass,chlB_mass.jack.val.pred),
                               percentRMSE=percentRMSD(meta(val.jack)$chlB_mass,chlB_mass.jack.val.pred,0.025,0.975),
                               bias=mean(chlB_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlB_mass,na.rm=T))
  
  car_mass.jack.val.pred<-as.vector(predict(car_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_car_mass_CVmodel)[,,1])
  car_mass.jack.val.fit<-lm(car_mass.jack.val.pred~meta(val.jack)$car_mass)
  car_mass.jack.stats[[i]]<-c(R2=summary(car_mass.jack.val.fit)$r.squared,
                              RMSE=RMSD(meta(val.jack)$car_mass,car_mass.jack.val.pred),
                              percentRMSE=percentRMSD(meta(val.jack)$car_mass,car_mass.jack.val.pred,0.025,0.975),
                              bias=mean(car_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$car_mass,na.rm=T))
  
  Cmass.jack.val.pred<-as.vector(predict(Cmass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Cmass_CVmodel)[,,1])
  Cmass.jack.val.fit<-lm(Cmass.jack.val.pred~meta(val.jack)$Cmass)
  Cmass.jack.stats[[i]]<-c(R2=summary(Cmass.jack.val.fit)$r.squared,
                           RMSE=RMSD(meta(val.jack)$Cmass,Cmass.jack.val.pred),
                           percentRMSE=percentRMSD(meta(val.jack)$Cmass,Cmass.jack.val.pred,0.025,0.975),
                           bias=mean(Cmass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Cmass,na.rm=T))
  
  Nmass.jack.val.pred<-as.vector(predict(Nmass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Nmass_CVmodel)[,,1])
  Nmass.jack.val.fit<-lm(Nmass.jack.val.pred~meta(val.jack)$Nmass)
  Nmass.jack.stats[[i]]<-c(R2=summary(Nmass.jack.val.fit)$r.squared,
                           RMSE=RMSD(meta(val.jack)$Nmass,Nmass.jack.val.pred),
                           percentRMSE=percentRMSD(meta(val.jack)$Nmass,Nmass.jack.val.pred,0.025,0.975),
                           bias=mean(Nmass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Nmass,na.rm=T))
  
  EWT.jack.val.pred<-as.vector(predict(EWT.jack,newdata=as.matrix(val.jack),ncomp=ncomp_EWT_CVmodel)[,,1])
  EWT.jack.val.fit<-lm(EWT.jack.val.pred~meta(val.jack)$EWT)
  EWT.jack.stats[[i]]<-c(R2=summary(EWT.jack.val.fit)$r.squared,
                         RMSE=RMSD(meta(val.jack)$EWT,EWT.jack.val.pred),
                         percentRMSE=percentRMSD(meta(val.jack)$EWT,EWT.jack.val.pred,0.025,0.975),
                         bias=mean(EWT.jack.val.pred,na.rm=T)-mean(meta(val.jack)$EWT,na.rm=T))
  
  LMA.jack.val.pred<-as.vector(predict(LMA.jack,newdata=as.matrix(val.jack),ncomp=ncomp_LMA_CVmodel)[,,1])
  LMA.jack.val.fit<-lm(LMA.jack.val.pred~meta(val.jack)$LMA)
  LMA.jack.stats[[i]]<-c(R2=summary(LMA.jack.val.fit)$r.squared,
                         RMSE=RMSD(meta(val.jack)$LMA,LMA.jack.val.pred),
                         percentRMSE=percentRMSD(meta(val.jack)$LMA,LMA.jack.val.pred,0.025,0.975),
                         bias=mean(LMA.jack.val.pred,na.rm=T)-mean(meta(val.jack)$LMA,na.rm=T))
  
  LDMC.jack.val.pred<-as.vector(predict(LDMC.jack,newdata=as.matrix(val.jack),ncomp=ncomp_LDMC_CVmodel)[,,1])
  LDMC.jack.val.fit<-lm(LDMC.jack.val.pred~meta(val.jack)$LDMC)
  LDMC.jack.stats[[i]]<-c(R2=summary(LDMC.jack.val.fit)$r.squared,
                          RMSE=RMSD(meta(val.jack)$LDMC,LDMC.jack.val.pred),
                          percentRMSE=percentRMSD(meta(val.jack)$LDMC,LDMC.jack.val.pred,0.025,0.975),
                          bias=mean(LDMC.jack.val.pred,na.rm=T)-mean(meta(val.jack)$LDMC,na.rm=T))
  
  Al_mass.jack.val.pred<-as.vector(predict(Al_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Al_mass_CVmodel)[,,1])
  Al_mass.jack.val.fit<-lm(Al_mass.jack.val.pred~meta(val.jack)$Al_mass)
  Al_mass.jack.stats[[i]]<-c(R2=summary(Al_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Al_mass,Al_mass.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Al_mass,Al_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Al_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Al_mass,na.rm=T))
  
  Ca_mass.jack.val.pred<-as.vector(predict(Ca_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Ca_mass_CVmodel)[,,1])
  Ca_mass.jack.val.fit<-lm(Ca_mass.jack.val.pred~meta(val.jack)$Ca_mass)
  Ca_mass.jack.stats[[i]]<-c(R2=summary(Ca_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Ca_mass,Ca_mass.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Ca_mass,Ca_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Ca_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Ca_mass,na.rm=T))
  
  Cu_mass.jack.val.pred<-as.vector(predict(Cu_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Cu_mass_CVmodel)[,,1])
  Cu_mass.jack.val.fit<-lm(Cu_mass.jack.val.pred~meta(val.jack)$Cu_mass)
  Cu_mass.jack.stats[[i]]<-c(R2=summary(Cu_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Cu_mass,Cu_mass.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Cu_mass,Cu_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Cu_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Cu_mass,na.rm=T))
  
  Fe_mass.jack.val.pred<-as.vector(predict(Fe_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Fe_mass_CVmodel)[,,1])
  Fe_mass.jack.val.fit<-lm(Fe_mass.jack.val.pred~meta(val.jack)$Fe_mass)
  Fe_mass.jack.stats[[i]]<-c(R2=summary(Fe_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Fe_mass,Fe_mass.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Fe_mass,Fe_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Fe_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Fe_mass,na.rm=T))
  
  K_mass.jack.val.pred<-as.vector(predict(K_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_K_mass_CVmodel)[,,1])
  K_mass.jack.val.fit<-lm(K_mass.jack.val.pred~meta(val.jack)$K_mass)
  K_mass.jack.stats[[i]]<-c(R2=summary(K_mass.jack.val.fit)$r.squared,
                            RMSE=RMSD(meta(val.jack)$K_mass,K_mass.jack.val.pred),
                            percentRMSE=percentRMSD(meta(val.jack)$K_mass,K_mass.jack.val.pred,0.025,0.975),
                            bias=mean(K_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$K_mass,na.rm=T))
  
  Mg_mass.jack.val.pred<-as.vector(predict(Mg_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Mg_mass_CVmodel)[,,1])
  Mg_mass.jack.val.fit<-lm(Mg_mass.jack.val.pred~meta(val.jack)$Mg_mass)
  Mg_mass.jack.stats[[i]]<-c(R2=summary(Mg_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Mg_mass,Mg_mass.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Mg_mass,Mg_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Mg_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Mg_mass,na.rm=T))
  
  Mn_mass.jack.val.pred<-as.vector(predict(Mn_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Mn_mass_CVmodel)[,,1])
  Mn_mass.jack.val.fit<-lm(Mn_mass.jack.val.pred~meta(val.jack)$Mn_mass)
  Mn_mass.jack.stats[[i]]<-c(R2=summary(Mn_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Mn_mass,Mn_mass.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Mn_mass,Mn_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Mn_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Mn_mass,na.rm=T))
  
  Na_mass.jack.val.pred<-as.vector(predict(Na_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Na_mass_CVmodel)[,,1])
  Na_mass.jack.val.fit<-lm(Na_mass.jack.val.pred~meta(val.jack)$Na_mass)
  Na_mass.jack.stats[[i]]<-c(R2=summary(Na_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Na_mass,Na_mass.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Na_mass,Na_mass.jack.val.pred,0.025,0.975),
                             bias=mean(Na_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Na_mass,na.rm=T))
  
  P_mass.jack.val.pred<-as.vector(predict(P_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_P_mass_CVmodel)[,,1])
  P_mass.jack.val.fit<-lm(P_mass.jack.val.pred~meta(val.jack)$P_mass)
  P_mass.jack.stats[[i]]<-c(R2=summary(P_mass.jack.val.fit)$r.squared,
                            RMSE=RMSD(meta(val.jack)$P_mass,P_mass.jack.val.pred),
                            percentRMSE=percentRMSD(meta(val.jack)$P_mass,P_mass.jack.val.pred,0.025,0.975),
                            bias=mean(P_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$P_mass,na.rm=T))
  
  Zn_mass.jack.val.pred<-as.vector(predict(Zn_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Zn_mass_CVmodel)[,,1])
  Zn_mass.jack.val.fit<-lm(Zn_mass.jack.val.pred~meta(val.jack)$Zn_mass)
  Zn_mass.jack.stats[[i]]<-c(R2=summary(Zn_mass.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$Zn_mass,Zn_mass.jack.val.pred),
                             percentRMSE=percentRMSD(meta(val.jack)$Zn_mass,Zn_mass.jack.val.pred,0.025,0.975),
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

solubles_mass.jack.pred<-apply.coefs(solubles_mass.jack.coefs,as.matrix(ref.test))
solubles_mass.jack.stat<-t(apply(solubles_mass.jack.pred,1,
                                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
solubles_mass.jack.df<-data.frame(pred.mean=solubles_mass.jack.stat[,1],
                                  pred.low=solubles_mass.jack.stat[,2],
                                  pred.high=solubles_mass.jack.stat[,3],
                                  Measured=meta(ref.test)$solubles_mass,
                                  ncomp=ncomp_solubles_mass_CVmodel,
                                  Species=meta(ref.test)$species,
                                  Project=meta(ref.test)$project,
                                  functional.group=meta(ref.test)$functional.group,
                                  ID=meta(ref.test)$sample_id)

solubles_mass.val.plot<-ggplot(solubles_mass.jack.df,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(25,90),ylim=c(25,90))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)

hemicellulose_mass.jack.pred<-apply.coefs(hemicellulose_mass.jack.coefs,as.matrix(ref.test))
hemicellulose_mass.jack.stat<-t(apply(hemicellulose_mass.jack.pred,1,
                                      function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
hemicellulose_mass.jack.df<-data.frame(pred.mean=hemicellulose_mass.jack.stat[,1],
                                       pred.low=hemicellulose_mass.jack.stat[,2],
                                       pred.high=hemicellulose_mass.jack.stat[,3],
                                       Measured=meta(ref.test)$hemicellulose_mass,
                                       ncomp=ncomp_hemicellulose_mass_CVmodel,
                                       Species=meta(ref.test)$species,
                                       Project=meta(ref.test)$project,
                                       functional.group=meta(ref.test)$functional.group,
                                       ID=meta(ref.test)$sample_id)

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

cellulose_mass.jack.pred<-apply.coefs(cellulose_mass.jack.coefs,as.matrix(ref.test))
cellulose_mass.jack.stat<-t(apply(cellulose_mass.jack.pred,1,
                                  function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
cellulose_mass.jack.df<-data.frame(pred.mean=cellulose_mass.jack.stat[,1],
                                   pred.low=cellulose_mass.jack.stat[,2],
                                   pred.high=cellulose_mass.jack.stat[,3],
                                   Measured=meta(ref.test)$cellulose_mass,
                                   ncomp=ncomp_cellulose_mass_CVmodel,
                                   Species=meta(ref.test)$species,
                                   Project=meta(ref.test)$project,
                                   functional.group=meta(ref.test)$functional.group,
                                   ID=meta(ref.test)$sample_id)

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

lignin_mass.jack.pred<-apply.coefs(lignin_mass.jack.coefs,as.matrix(ref.test))
lignin_mass.jack.stat<-t(apply(lignin_mass.jack.pred,1,
                               function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
lignin_mass.jack.df<-data.frame(pred.mean=lignin_mass.jack.stat[,1],
                                pred.low=lignin_mass.jack.stat[,2],
                                pred.high=lignin_mass.jack.stat[,3],
                                Measured=meta(ref.test)$lignin_mass,
                                ncomp=ncomp_lignin_mass_CVmodel,
                                Species=meta(ref.test)$species,
                                Project=meta(ref.test)$project,
                                functional.group=meta(ref.test)$functional.group,
                                ID=meta(ref.test)$sample_id)

lignin_mass.val.plot<-ggplot(lignin_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-2.5,22.5),ylim=c(-2.5,22.5))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)

chlA_mass.jack.pred<-apply.coefs(chlA_mass.jack.coefs,as.matrix(ref.test))
chlA_mass.jack.stat<-t(apply(chlA_mass.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_mass.jack.df<-data.frame(pred.mean=chlA_mass.jack.stat[,1],
                              pred.low=chlA_mass.jack.stat[,2],
                              pred.high=chlA_mass.jack.stat[,3],
                              Measured=meta(ref.test)$chlA_mass,
                              ncomp=ncomp_chlA_mass_CVmodel,
                              Species=meta(ref.test)$species,
                              Project=meta(ref.test)$project,
                              functional.group=meta(ref.test)$functional.group,
                              ID=meta(ref.test)$sample_id)

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

chlB_mass.jack.pred<-apply.coefs(chlB_mass.jack.coefs,as.matrix(ref.test))
chlB_mass.jack.stat<-t(apply(chlB_mass.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlB_mass.jack.df<-data.frame(pred.mean=chlB_mass.jack.stat[,1],
                              pred.low=chlB_mass.jack.stat[,2],
                              pred.high=chlB_mass.jack.stat[,3],
                              Measured=meta(ref.test)$chlB_mass,
                              ncomp=ncomp_chlB_mass_CVmodel,
                              Species=meta(ref.test)$species,
                              Project=meta(ref.test)$project,
                              functional.group=meta(ref.test)$functional.group,
                              ID=meta(ref.test)$sample_id)

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

car_mass.jack.pred<-apply.coefs(car_mass.jack.coefs,as.matrix(ref.test))
car_mass.jack.stat<-t(apply(car_mass.jack.pred,1,
                            function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
car_mass.jack.df<-data.frame(pred.mean=car_mass.jack.stat[,1],
                             pred.low=car_mass.jack.stat[,2],
                             pred.high=car_mass.jack.stat[,3],
                             Measured=meta(ref.test)$car_mass,
                             ncomp=ncomp_car_mass_CVmodel,
                             Species=meta(ref.test)$species,
                             Project=meta(ref.test)$project,
                             functional.group=meta(ref.test)$functional.group,
                             ID=meta(ref.test)$sample_id)

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

Cmass.jack.pred<-apply.coefs(Cmass.jack.coefs,as.matrix(ref.test))
Cmass.jack.stat<-t(apply(Cmass.jack.pred,1,
                         function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cmass.jack.df<-data.frame(pred.mean=Cmass.jack.stat[,1],
                          pred.low=Cmass.jack.stat[,2],
                          pred.high=Cmass.jack.stat[,3],
                          Measured=meta(ref.test)$Cmass,
                          ncomp=ncomp_Cmass_CVmodel,
                          Species=meta(ref.test)$species,
                          Project=meta(ref.test)$project,
                          functional.group=meta(ref.test)$functional.group,
                          ID=meta(ref.test)$sample_id)

Cmass.val.plot<-ggplot(Cmass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(36,55),ylim=c(36,55))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured C"[mass]*" (%)"),
       x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)

Nmass.jack.pred<-apply.coefs(Nmass.jack.coefs,as.matrix(ref.test))
Nmass.jack.stat<-t(apply(Nmass.jack.pred,1,
                         function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Nmass.jack.df<-data.frame(pred.mean=Nmass.jack.stat[,1],
                          pred.low=Nmass.jack.stat[,2],
                          pred.high=Nmass.jack.stat[,3],
                          Measured=meta(ref.test)$Nmass,
                          ncomp=ncomp_Nmass_CVmodel,
                          Species=meta(ref.test)$species,
                          Project=meta(ref.test)$project,
                          functional.group=meta(ref.test)$functional.group,
                          ID=meta(ref.test)$sample_id)

Nmass.val.plot<-ggplot(Nmass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,5.5),ylim=c(0,5.5))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured N"[mass]*" (%)"),
       x=expression("Predicted N"[mass]*" (%)"))+
  guides(color=F)

EWT.jack.pred<-apply.coefs(EWT.jack.coefs,as.matrix(ref.test))
EWT.jack.stat<-t(apply(EWT.jack.pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT.jack.df<-data.frame(pred.mean=EWT.jack.stat[,1],
                        pred.low=EWT.jack.stat[,2],
                        pred.high=EWT.jack.stat[,3],
                        Measured=meta(ref.test)$EWT,
                        ncomp=ncomp_EWT_CVmodel,
                        Species=meta(ref.test)$species,
                        Project=meta(ref.test)$project,
                        functional.group=meta(ref.test)$functional.group,
                        ID=meta(ref.test)$sample_id)

EWT.val.plot<-ggplot(EWT.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.068),ylim=c(0,0.068))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (cm)",x="Predicted EWT (cm)")+
  guides(color=F)

LMA.jack.pred<-apply.coefs(LMA.jack.coefs,as.matrix(ref.test))
LMA.jack.stat<-t(apply(LMA.jack.pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA.jack.df<-data.frame(pred.mean=LMA.jack.stat[,1],
                        pred.low=LMA.jack.stat[,2],
                        pred.high=LMA.jack.stat[,3],
                        Measured=meta(ref.test)$LMA,
                        ncomp=ncomp_LMA_CVmodel,
                        Species=meta(ref.test)$species,
                        Project=meta(ref.test)$project,
                        functional.group=meta(ref.test)$functional.group,
                        ID=meta(ref.test)$sample_id)

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

LDMC.jack.pred<-apply.coefs(LDMC.jack.coefs,as.matrix(ref.test))
LDMC.jack.stat<-t(apply(LDMC.jack.pred,1,
                        function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LDMC.jack.df<-data.frame(pred.mean=LDMC.jack.stat[,1],
                         pred.low=LDMC.jack.stat[,2],
                         pred.high=LDMC.jack.stat[,3],
                         Measured=meta(ref.test)$LDMC,
                         ncomp=ncomp_LDMC_CVmodel,
                         Species=meta(ref.test)$species,
                         Project=meta(ref.test)$project,
                         functional.group=meta(ref.test)$functional.group,
                         ID=meta(ref.test)$sample_id)

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

Al_mass.jack.pred<-apply.coefs(Al_mass.jack.coefs,as.matrix(ref.test))
Al_mass.jack.stat<-t(apply(Al_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Al_mass.jack.df<-data.frame(pred.mean=Al_mass.jack.stat[,1],
                            pred.low=Al_mass.jack.stat[,2],
                            pred.high=Al_mass.jack.stat[,3],
                            Measured=meta(ref.test)$Al_mass,
                            ncomp=ncomp_Al_mass_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

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

Ca_mass.jack.pred<-apply.coefs(Ca_mass.jack.coefs,as.matrix(ref.test))
Ca_mass.jack.stat<-t(apply(Ca_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Ca_mass.jack.df<-data.frame(pred.mean=Ca_mass.jack.stat[,1],
                            pred.low=Ca_mass.jack.stat[,2],
                            pred.high=Ca_mass.jack.stat[,3],
                            Measured=meta(ref.test)$Ca_mass,
                            ncomp=ncomp_Ca_mass_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

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

Cu_mass.jack.pred<-apply.coefs(Cu_mass.jack.coefs,as.matrix(ref.test))
Cu_mass.jack.stat<-t(apply(Cu_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cu_mass.jack.df<-data.frame(pred.mean=Cu_mass.jack.stat[,1],
                            pred.low=Cu_mass.jack.stat[,2],
                            pred.high=Cu_mass.jack.stat[,3],
                            Measured=meta(ref.test)$Cu_mass,
                            ncomp=ncomp_Cu_mass_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

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

Fe_mass.jack.pred<-apply.coefs(Fe_mass.jack.coefs,as.matrix(ref.test))
Fe_mass.jack.stat<-t(apply(Fe_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Fe_mass.jack.df<-data.frame(pred.mean=Fe_mass.jack.stat[,1],
                            pred.low=Fe_mass.jack.stat[,2],
                            pred.high=Fe_mass.jack.stat[,3],
                            Measured=meta(ref.test)$Fe_mass,
                            ncomp=ncomp_Fe_mass_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

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

K_mass.jack.pred<-apply.coefs(K_mass.jack.coefs,as.matrix(ref.test))
K_mass.jack.stat<-t(apply(K_mass.jack.pred,1,
                          function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_mass.jack.df<-data.frame(pred.mean=K_mass.jack.stat[,1],
                           pred.low=K_mass.jack.stat[,2],
                           pred.high=K_mass.jack.stat[,3],
                           Measured=meta(ref.test)$K_mass,
                           ncomp=ncomp_K_mass_CVmodel,
                           Species=meta(ref.test)$species,
                           Project=meta(ref.test)$project,
                           functional.group=meta(ref.test)$functional.group,
                           ID=meta(ref.test)$sample_id)

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

Mg_mass.jack.pred<-apply.coefs(Mg_mass.jack.coefs,as.matrix(ref.test))
Mg_mass.jack.stat<-t(apply(Mg_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mg_mass.jack.df<-data.frame(pred.mean=Mg_mass.jack.stat[,1],
                            pred.low=Mg_mass.jack.stat[,2],
                            pred.high=Mg_mass.jack.stat[,3],
                            Measured=meta(ref.test)$Mg_mass,
                            ncomp=ncomp_Mg_mass_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

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

Mn_mass.jack.pred<-apply.coefs(Mn_mass.jack.coefs,as.matrix(ref.test))
Mn_mass.jack.stat<-t(apply(Mn_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Mn_mass.jack.df<-data.frame(pred.mean=Mn_mass.jack.stat[,1],
                            pred.low=Mn_mass.jack.stat[,2],
                            pred.high=Mn_mass.jack.stat[,3],
                            Measured=meta(ref.test)$Mn_mass,
                            ncomp=ncomp_Mn_mass_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

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

Na_mass.jack.pred<-apply.coefs(Na_mass.jack.coefs,as.matrix(ref.test))
Na_mass.jack.stat<-t(apply(Na_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Na_mass.jack.df<-data.frame(pred.mean=Na_mass.jack.stat[,1],
                            pred.low=Na_mass.jack.stat[,2],
                            pred.high=Na_mass.jack.stat[,3],
                            Measured=meta(ref.test)$Na_mass,
                            ncomp=ncomp_Na_mass_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

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

P_mass.jack.pred<-apply.coefs(P_mass.jack.coefs,as.matrix(ref.test))
P_mass.jack.stat<-t(apply(P_mass.jack.pred,1,
                          function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
P_mass.jack.df<-data.frame(pred.mean=P_mass.jack.stat[,1],
                           pred.low=P_mass.jack.stat[,2],
                           pred.high=P_mass.jack.stat[,3],
                           Measured=meta(ref.test)$P_mass,
                           ncomp=ncomp_P_mass_CVmodel,
                           Species=meta(ref.test)$species,
                           Project=meta(ref.test)$project,
                           functional.group=meta(ref.test)$functional.group,
                           ID=meta(ref.test)$sample_id)

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

Zn_mass.jack.pred<-apply.coefs(Zn_mass.jack.coefs,as.matrix(ref.test))
Zn_mass.jack.stat<-t(apply(Zn_mass.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Zn_mass.jack.df<-data.frame(pred.mean=Zn_mass.jack.stat[,1],
                            pred.low=Zn_mass.jack.stat[,2],
                            pred.high=Zn_mass.jack.stat[,3],
                            Measured=meta(ref.test)$Zn_mass,
                            ncomp=ncomp_Zn_mass_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

Zn_mass.val.plot<-ggplot(Zn_mass.jack.df,aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.04,0.35),ylim=c(-0.04,0.35))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Zn (mg g"^-1*")"),
       x=expression("Predicted Zn (mg g"^-1*")"))+
  guides(color=F)

all.jack.coef.list<-list(solubles_mass=solubles_mass.jack.coefs,
                         hemicellulose_mass=hemicellulose_mass.jack.coefs,
                         cellulose_mass=cellulose_mass.jack.coefs,
                         lignin_mass=lignin_mass.jack.coefs,
                         chlA_mass=chlA_mass.jack.coefs,
                         chlB_mass=chlB_mass.jack.coefs,
                         car_mass=car_mass.jack.coefs,
                         Nmass=Nmass.jack.coefs,
                         Cmass=Cmass.jack.coefs,
                         EWT=EWT.jack.coefs,
                         LDMC=LDMC.jack.coefs,
                         LMA=LMA.jack.coefs,
                         Al_mass=Al_mass.jack.coefs,
                         Ca_mass=Ca_mass.jack.coefs,
                         Cu_mass=Cu_mass.jack.coefs,
                         Fe_mass=Fe_mass.jack.coefs,
                         K_mass=K_mass.jack.coefs,
                         Mg_mass=Mg_mass.jack.coefs,
                         Mn_mass=Mn_mass.jack.coefs,
                         Na_mass=Na_mass.jack.coefs,
                         P_mass=P_mass.jack.coefs,
                         Zn_mass=Zn_mass.jack.coefs)
saveRDS(all.jack.coef.list,"SavedResults/all_jack_coefs_list_ref.rds")

all.jack.df.list<-list(solubles_mass=solubles_mass.jack.df,
                       hemicellulose_mass=hemicellulose_mass.jack.df,
                       cellulose_mass=cellulose_mass.jack.df,
                       lignin_mass=lignin_mass.jack.df,
                       chlA_mass=chlA_mass.jack.df,
                       chlB_mass=chlB_mass.jack.df,
                       car_mass=car_mass.jack.df,
                       Nmass=Nmass.jack.df,
                       Cmass=Cmass.jack.df,
                       EWT=EWT.jack.df,
                       LDMC=LDMC.jack.df,
                       LMA=LMA.jack.df,
                       Al_mass=Al_mass.jack.df,
                       Ca_mass=Ca_mass.jack.df,
                       Cu_mass=Cu_mass.jack.df,
                       Fe_mass=Fe_mass.jack.df,
                       K_mass=K_mass.jack.df,
                       Mg_mass=Mg_mass.jack.df,
                       Mn_mass=Mn_mass.jack.df,
                       Na_mass=Na_mass.jack.df,
                       P_mass=P_mass.jack.df,
                       Zn_mass=Zn_mass.jack.df)
saveRDS(all.jack.df.list,"SavedResults/all_jack_df_list_ref.rds")

all.jack.stats.list<-list(solubles_mass=solubles_mass.jack.stats,
                       hemicellulose_mass=hemicellulose_mass.jack.stats,
                       cellulose_mass=cellulose_mass.jack.stats,
                       lignin_mass=lignin_mass.jack.stats,
                       chlA_mass=chlA_mass.jack.stats,
                       chlB_mass=chlB_mass.jack.stats,
                       car_mass=car_mass.jack.stats,
                       Nmass=Nmass.jack.stats,
                       Cmass=Cmass.jack.stats,
                       EWT=EWT.jack.stats,
                       LDMC=LDMC.jack.stats,
                       LMA=LMA.jack.stats,
                       Al_mass=Al_mass.jack.stats,
                       Ca_mass=Ca_mass.jack.stats,
                       Cu_mass=Cu_mass.jack.stats,
                       Fe_mass=Fe_mass.jack.stats,
                       K_mass=K_mass.jack.stats,
                       Mg_mass=Mg_mass.jack.stats,
                       Mn_mass=Mn_mass.jack.stats,
                       Na_mass=Na_mass.jack.stats,
                       P_mass=P_mass.jack.stats,
                       Zn_mass=Zn_mass.jack.stats)
saveRDS(all.jack.stats.list,"SavedResults/all_jack_stats_list_ref.rds")

pdf("Images/val_plots_ref1.pdf",width = 16,height = 20)
(EWT.val.plot + LDMC.val.plot + LMA.val.plot) / 
  (cellulose_mass.val.plot + solubles_mass.val.plot + Cmass.val.plot) / 
  (hemicellulose_mass.val.plot + Nmass.val.plot + chlB_mass.val.plot) / 
  (chlA_mass.val.plot + lignin_mass.val.plot + car_mass.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Images/val_plots_ref2.pdf",width = 16,height = 20)
(Ca_mass.val.plot + K_mass.val.plot + Zn_mass.val.plot) / 
  (P_mass.val.plot + Mg_mass.val.plot + Na_mass.val.plot) / 
  (Fe_mass.val.plot + Mn_mass.val.plot + Al_mass.val.plot) / 
  (Cu_mass.val.plot + guide_area() + guide_area()) &
  plot_layout(guides="collect")
dev.off()
