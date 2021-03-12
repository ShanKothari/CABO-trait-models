setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(caret)
library(pls)
library(lmodel2)
library(reshape2)
source("Scripts/VIP.R")

spec.traits<-readRDS("ProcessedSpectra/all_spectra_and_traits.rds")

##########################################
## to dos

## fix RMSD calculations
## try Type II regression?
## remove Dessain data and instead use as a validation data set

#########################################
## define functions

## RMSD between predicted and observed values
RMSD<-function(measured,predicted){
  not.na<-which(!is.na(measured) & !is.na(predicted))
  return(sqrt(sum((measured-predicted)^2,na.rm=T)/(length(not.na)-1)))
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

##########################################
## divide into training and testing

train.sample <- createDataPartition(
  y = meta(spec.traits)$project,
  p = .8,
  list = FALSE
)

test.sample<-setdiff(1:nrow(as.matrix(spec.traits)),train.sample)

spec.train<-spec.traits[train.sample,]
spec.test<-spec.traits[test.sample,]

###########################################
## start building first-pass models using
## K-fold cross validation to get the
## optimal number of components

EWT_CVmodel<-plsr(meta(spec.train)$EWT~as.matrix(spec.train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_CVmodel <- selectNcomp(EWT_CVmodel, method = "onesigma", plot = TRUE)
EWT_valid <- which(!is.na(meta(spec.train)$EWT))
EWT_pred<-data.frame(ID=meta(spec.train)$sample_id[EWT_valid],
                     Species=meta(spec.train)$species[EWT_valid],
                     Project=meta(spec.train)$project[EWT_valid],
                     measured=meta(spec.train)$EWT[EWT_valid],
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

solubles_mass_CVmodel<-plsr(meta(spec.train)$solubles_mass~as.matrix(spec.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_solubles_mass_CVmodel <- selectNcomp(solubles_mass_CVmodel, method = "onesigma", plot = TRUE)
solubles_mass_valid <- which(!is.na(meta(spec.train)$solubles_mass))
solubles_mass_pred<-data.frame(ID=meta(spec.train)$sample_id[solubles_mass_valid],
                         Species=meta(spec.train)$species[solubles_mass_valid],
                         Project=meta(spec.train)$project[solubles_mass_valid],
                         measured=meta(spec.train)$solubles_mass[solubles_mass_valid],
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

hemicellulose_mass_CVmodel<-plsr(meta(spec.train)$hemicellulose_mass~as.matrix(spec.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_hemicellulose_mass_CVmodel <- selectNcomp(hemicellulose_mass_CVmodel, method = "onesigma", plot = TRUE)
hemicellulose_mass_valid <- which(!is.na(meta(spec.train)$hemicellulose_mass))
hemicellulose_mass_pred<-data.frame(ID=meta(spec.train)$sample_id[hemicellulose_mass_valid],
                         Species=meta(spec.train)$species[hemicellulose_mass_valid],
                         Project=meta(spec.train)$project[hemicellulose_mass_valid],
                         measured=meta(spec.train)$hemicellulose_mass[hemicellulose_mass_valid],
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

cellulose_mass_CVmodel<-plsr(meta(spec.train)$cellulose_mass~as.matrix(spec.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_cellulose_mass_CVmodel <- selectNcomp(cellulose_mass_CVmodel, method = "onesigma", plot = TRUE)
cellulose_mass_valid <- which(!is.na(meta(spec.train)$cellulose_mass))
cellulose_mass_pred<-data.frame(ID=meta(spec.train)$sample_id[cellulose_mass_valid],
                         Species=meta(spec.train)$species[cellulose_mass_valid],
                         Project=meta(spec.train)$project[cellulose_mass_valid],
                         measured=meta(spec.train)$cellulose_mass[cellulose_mass_valid],
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

lignin_mass_CVmodel<-plsr(meta(spec.train)$lignin_mass~as.matrix(spec.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_lignin_mass_CVmodel <- selectNcomp(lignin_mass_CVmodel, method = "onesigma", plot = TRUE)
lignin_mass_valid <- which(!is.na(meta(spec.train)$lignin_mass))
lignin_mass_pred<-data.frame(ID=meta(spec.train)$sample_id[lignin_mass_valid],
                         Species=meta(spec.train)$species[lignin_mass_valid],
                         Project=meta(spec.train)$project[lignin_mass_valid],
                         measured=meta(spec.train)$lignin_mass[lignin_mass_valid],
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

Cmass_CVmodel<-plsr(meta(spec.train)$Cmass~as.matrix(spec.train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cmass_CVmodel <- selectNcomp(Cmass_CVmodel, method = "onesigma", plot = TRUE)
Cmass_valid <- which(!is.na(meta(spec.train)$Cmass))
Cmass_pred<-data.frame(ID=meta(spec.train)$sample_id[Cmass_valid],
                            Species=meta(spec.train)$species[Cmass_valid],
                            Project=meta(spec.train)$project[Cmass_valid],
                            measured=meta(spec.train)$Cmass[Cmass_valid],
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

Nmass_CVmodel<-plsr(meta(spec.train)$Nmass~as.matrix(spec.train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Nmass_CVmodel <- selectNcomp(Nmass_CVmodel, method = "onesigma", plot = TRUE)
Nmass_valid <- which(!is.na(meta(spec.train)$Nmass))
Nmass_pred<-data.frame(ID=meta(spec.train)$sample_id[Nmass_valid],
                      Species=meta(spec.train)$species[Nmass_valid],
                      Project=meta(spec.train)$project[Nmass_valid],
                      measured=meta(spec.train)$Nmass[Nmass_valid],
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

# Narea_CVmodel<-plsr(meta(spec.train)$Narea~as.matrix(spec.train),
#                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
# ncomp_Narea_CVmodel <- selectNcomp(Narea_CVmodel, method = "onesigma", plot = TRUE)
# Narea_valid <- which(!is.na(meta(spec.train)$Narea))
# Narea_pred<-data.frame(ID=meta(spec.train)$sample_id[Narea_valid],
#                        Species=meta(spec.train)$species[Narea_valid],
#                        Project=meta(spec.train)$project[Narea_valid],
#                        measured=meta(spec.train)$Narea[Narea_valid],
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
# Nnorm_CVmodel<-plsr(meta(spec.train)$Nnorm~as.matrix(spec.train),
#                     ncomp=30,method = "oscorespls",validation="CV",segments=10)
# ncomp_Nnorm_CVmodel <- selectNcomp(Nnorm_CVmodel, method = "onesigma", plot = TRUE)
# Nnorm_valid <- which(!is.na(meta(spec.train)$Nnorm))
# Nnorm_pred<-data.frame(ID=meta(spec.train)$sample_id[Nnorm_valid],
#                        Species=meta(spec.train)$species[Nnorm_valid],
#                        Project=meta(spec.train)$project[Nnorm_valid],
#                        measured=meta(spec.train)$Nnorm[Nnorm_valid],
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

LDMC_CVmodel<-plsr(meta(spec.train)$LDMC~as.matrix(spec.train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LDMC_CVmodel <- selectNcomp(LDMC_CVmodel, method = "onesigma", plot = TRUE)
LDMC_valid <- which(!is.na(meta(spec.train)$LDMC))
LDMC_pred<-data.frame(ID=meta(spec.train)$sample_id[LDMC_valid],
                      Species=meta(spec.train)$species[LDMC_valid],
                      Project=meta(spec.train)$project[LDMC_valid],
                      measured=meta(spec.train)$LDMC[LDMC_valid],
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

LMA_CVmodel<-plsr(meta(spec.train)$LMA~as.matrix(spec.train),
                   ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_CVmodel <- selectNcomp(LMA_CVmodel, method = "onesigma", plot = TRUE)
LMA_valid <- which(!is.na(meta(spec.train)$LMA))
LMA_pred<-data.frame(ID=meta(spec.train)$sample_id[LMA_valid],
                      Species=meta(spec.train)$species[LMA_valid],
                      Project=meta(spec.train)$project[LMA_valid],
                      measured=meta(spec.train)$LMA[LMA_valid],
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

chlA_fresh_CVmodel<-plsr(meta(spec.train)$chlA_fresh~as.matrix(spec.train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_fresh_CVmodel <- selectNcomp(chlA_fresh_CVmodel, method = "onesigma", plot = TRUE)
chlA_fresh_valid <- which(!is.na(meta(spec.train)$chlA_fresh))
chlA_fresh_pred<-data.frame(ID=meta(spec.train)$sample_id[chlA_fresh_valid],
                     Species=meta(spec.train)$species[chlA_fresh_valid],
                     Project=meta(spec.train)$project[chlA_fresh_valid],
                     measured=meta(spec.train)$chlA_fresh[chlA_fresh_valid],
                     val_pred=chlA_fresh_CVmodel$validation$pred[,,ncomp_chlA_fresh_CVmodel])
ggplot(chlA_fresh_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting chlA from fresh-leaf spectra")

chlB_fresh_CVmodel<-plsr(meta(spec.train)$chlB_fresh~as.matrix(spec.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlB_fresh_CVmodel <- selectNcomp(chlB_fresh_CVmodel, method = "onesigma", plot = TRUE)
chlB_fresh_valid <- which(!is.na(meta(spec.train)$chlB_fresh))
chlB_fresh_pred<-data.frame(ID=meta(spec.train)$sample_id[chlB_fresh_valid],
                            Species=meta(spec.train)$species[chlB_fresh_valid],
                            Project=meta(spec.train)$project[chlB_fresh_valid],
                            measured=meta(spec.train)$chlB_fresh[chlB_fresh_valid],
                            val_pred=chlB_fresh_CVmodel$validation$pred[,,ncomp_chlB_fresh_CVmodel])
ggplot(chlB_fresh_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting chlB from fresh-leaf spectra")

car_fresh_CVmodel<-plsr(meta(spec.train)$car_fresh~as.matrix(spec.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_car_fresh_CVmodel <- selectNcomp(car_fresh_CVmodel, method = "onesigma", plot = TRUE)
car_fresh_valid <- which(!is.na(meta(spec.train)$car_fresh))
car_fresh_pred<-data.frame(ID=meta(spec.train)$sample_id[car_fresh_valid],
                            Species=meta(spec.train)$species[car_fresh_valid],
                            Project=meta(spec.train)$project[car_fresh_valid],
                            measured=meta(spec.train)$car_fresh[car_fresh_valid],
                            val_pred=car_fresh_CVmodel$validation$pred[,,ncomp_car_fresh_CVmodel])
ggplot(car_fresh_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting carotenoids from fresh-leaf spectra")

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

EWT.jack.coefs<-list()
solubles_mass.jack.coefs<-list()
hemicellulose_mass.jack.coefs<-list()
cellulose_mass.jack.coefs<-list()
lignin_mass.jack.coefs<-list()
chlA_fresh.jack.coefs<-list()
chlB_fresh.jack.coefs<-list()
car_fresh.jack.coefs<-list()
Cmass.jack.coefs<-list()
Nmass.jack.coefs<-list()
LMA.jack.coefs<-list()
LDMC.jack.coefs<-list()

EWT.jack.stats<-list()
solubles_mass.jack.stats<-list()
hemicellulose_mass.jack.stats<-list()
cellulose_mass.jack.stats<-list()
lignin_mass.jack.stats<-list()
chlA_fresh.jack.stats<-list()
chlB_fresh.jack.stats<-list()
car_fresh.jack.stats<-list()
Cmass.jack.stats<-list()
Nmass.jack.stats<-list()
LMA.jack.stats<-list()
LDMC.jack.stats<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n.cal.spec<-nrow(spec.train)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-spec.train[train.jack]
  val.jack<-spec.train[test.jack]
  
  EWT.jack<-plsr(meta(calib.jack)$EWT~as.matrix(calib.jack),
                 ncomp=30,method = "oscorespls",validation="none")
  solubles_mass.jack<-plsr(meta(calib.jack)$solubles_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  hemicellulose_mass.jack<-plsr(meta(calib.jack)$hemicellulose_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  cellulose_mass.jack<-plsr(meta(calib.jack)$cellulose_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  lignin_mass.jack<-plsr(meta(calib.jack)$lignin_mass~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  chlA_fresh.jack<-plsr(meta(calib.jack)$chlA_fresh~as.matrix(calib.jack),
                           ncomp=30,method = "oscorespls",validation="none")
  chlB_fresh.jack<-plsr(meta(calib.jack)$chlB_fresh~as.matrix(calib.jack),
                           ncomp=30,method = "oscorespls",validation="none")
  car_fresh.jack<-plsr(meta(calib.jack)$car_fresh~as.matrix(calib.jack),
                           ncomp=30,method = "oscorespls",validation="none")
  Cmass.jack<-plsr(meta(calib.jack)$Cmass~as.matrix(calib.jack),
                          ncomp=30,method = "oscorespls",validation="none")
  Nmass.jack<-plsr(meta(calib.jack)$Nmass~as.matrix(calib.jack),
                          ncomp=30,method = "oscorespls",validation="none")
  LMA.jack<-plsr(meta(calib.jack)$LMA~as.matrix(calib.jack),
                         ncomp=30,method = "oscorespls",validation="none")
  LDMC.jack<-plsr(meta(calib.jack)$LDMC~as.matrix(calib.jack),
                          ncomp=30,method = "oscorespls",validation="none")

  EWT.jack.val.pred<-as.vector(predict(EWT.jack,newdata=as.matrix(val.jack),ncomp=ncomp_EWT_CVmodel)[,,1])
  EWT.jack.val.fit<-lm(EWT.jack.val.pred~meta(val.jack)$EWT)
  EWT.jack.stats[[i]]<-c(R2=summary(EWT.jack.val.fit)$r.squared,
                         RMSE=sqrt(c(crossprod(EWT.jack.val.fit$residuals))/length(EWT.jack.val.fit$residuals)),
                         max.val=max(meta(val.jack)$EWT,na.rm=T),
                         min.val=min(meta(val.jack)$EWT,na.rm=T),
                         bias=mean(EWT.jack.val.pred,na.rm=T)-mean(meta(val.jack)$EWT,na.rm=T))
    
  solubles_mass.jack.val.pred<-as.vector(predict(solubles_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_solubles_mass_CVmodel)[,,1])
  solubles_mass.jack.val.fit<-lm(solubles_mass.jack.val.pred~meta(val.jack)$solubles_mass)
  solubles_mass.jack.stats[[i]]<-c(R2=summary(solubles_mass.jack.val.fit)$r.squared,
                             RMSE=sqrt(c(crossprod(solubles_mass.jack.val.fit$residuals))/length(solubles_mass.jack.val.fit$residuals)),
                             max.val=max(meta(val.jack)$solubles_mass,na.rm=T),
                             min.val=min(meta(val.jack)$solubles_mass,na.rm=T),
                             bias=mean(solubles_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$solubles_mass,na.rm=T))
  
  hemicellulose_mass.jack.val.pred<-as.vector(predict(hemicellulose_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_hemicellulose_mass_CVmodel)[,,1])
  hemicellulose_mass.jack.val.fit<-lm(hemicellulose_mass.jack.val.pred~meta(val.jack)$hemicellulose_mass)
  hemicellulose_mass.jack.stats[[i]]<-c(R2=summary(hemicellulose_mass.jack.val.fit)$r.squared,
                             RMSE=sqrt(c(crossprod(hemicellulose_mass.jack.val.fit$residuals))/length(hemicellulose_mass.jack.val.fit$residuals)),
                             max.val=max(meta(val.jack)$hemicellulose_mass,na.rm=T),
                             min.val=min(meta(val.jack)$hemicellulose_mass,na.rm=T),
                             bias=mean(hemicellulose_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$hemicellulose_mass,na.rm=T))
  
  cellulose_mass.jack.val.pred<-as.vector(predict(cellulose_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_cellulose_mass_CVmodel)[,,1])
  cellulose_mass.jack.val.fit<-lm(cellulose_mass.jack.val.pred~meta(val.jack)$cellulose_mass)
  cellulose_mass.jack.stats[[i]]<-c(R2=summary(cellulose_mass.jack.val.fit)$r.squared,
                             RMSE=sqrt(c(crossprod(cellulose_mass.jack.val.fit$residuals))/length(cellulose_mass.jack.val.fit$residuals)),
                             max.val=max(meta(val.jack)$cellulose_mass,na.rm=T),
                             min.val=min(meta(val.jack)$cellulose_mass,na.rm=T),
                             bias=mean(cellulose_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$cellulose_mass,na.rm=T))
  
  lignin_mass.jack.val.pred<-as.vector(predict(lignin_mass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_lignin_mass_CVmodel)[,,1])
  lignin_mass.jack.val.fit<-lm(lignin_mass.jack.val.pred~meta(val.jack)$lignin_mass)
  lignin_mass.jack.stats[[i]]<-c(R2=summary(lignin_mass.jack.val.fit)$r.squared,
                             RMSE=sqrt(c(crossprod(lignin_mass.jack.val.fit$residuals))/length(lignin_mass.jack.val.fit$residuals)),
                             max.val=max(meta(val.jack)$lignin_mass,na.rm=T),
                             min.val=min(meta(val.jack)$lignin_mass,na.rm=T),
                             bias=mean(lignin_mass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$lignin_mass,na.rm=T))
  
  chlA_fresh.jack.val.pred<-as.vector(predict(chlA_fresh.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlA_fresh_CVmodel)[,,1])
  chlA_fresh.jack.val.fit<-lm(chlA_fresh.jack.val.pred~meta(val.jack)$chlA_fresh)
  chlA_fresh.jack.stats[[i]]<-c(R2=summary(chlA_fresh.jack.val.fit)$r.squared,
                                   RMSE=sqrt(c(crossprod(chlA_fresh.jack.val.fit$residuals))/length(chlA_fresh.jack.val.fit$residuals)),
                                   max.val=max(meta(val.jack)$chlA_fresh,na.rm=T),
                                   min.val=min(meta(val.jack)$chlA_fresh,na.rm=T),
                                   bias=mean(chlA_fresh.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlA_fresh,na.rm=T))
  
  chlB_fresh.jack.val.pred<-as.vector(predict(chlB_fresh.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlB_fresh_CVmodel)[,,1])
  chlB_fresh.jack.val.fit<-lm(chlB_fresh.jack.val.pred~meta(val.jack)$chlB_fresh)
  chlB_fresh.jack.stats[[i]]<-c(R2=summary(chlB_fresh.jack.val.fit)$r.squared,
                                   RMSE=sqrt(c(crossprod(chlB_fresh.jack.val.fit$residuals))/length(chlB_fresh.jack.val.fit$residuals)),
                                   max.val=max(meta(val.jack)$chlB_fresh,na.rm=T),
                                   min.val=min(meta(val.jack)$chlB_fresh,na.rm=T),
                                   bias=mean(chlB_fresh.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlB_fresh,na.rm=T))
  
  car_fresh.jack.val.pred<-as.vector(predict(car_fresh.jack,newdata=as.matrix(val.jack),ncomp=ncomp_car_fresh_CVmodel)[,,1])
  car_fresh.jack.val.fit<-lm(car_fresh.jack.val.pred~meta(val.jack)$car_fresh)
  car_fresh.jack.stats[[i]]<-c(R2=summary(car_fresh.jack.val.fit)$r.squared,
                                   RMSE=sqrt(c(crossprod(car_fresh.jack.val.fit$residuals))/length(car_fresh.jack.val.fit$residuals)),
                                   max.val=max(meta(val.jack)$car_fresh,na.rm=T),
                                   min.val=min(meta(val.jack)$car_fresh,na.rm=T),
                                   bias=mean(car_fresh.jack.val.pred,na.rm=T)-mean(meta(val.jack)$car_fresh,na.rm=T))
  
  Cmass.jack.val.pred<-as.vector(predict(Cmass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Cmass_CVmodel)[,,1])
  Cmass.jack.val.fit<-lm(Cmass.jack.val.pred~meta(val.jack)$Cmass)
  Cmass.jack.stats[[i]]<-c(R2=summary(Cmass.jack.val.fit)$r.squared,
                                  RMSE=sqrt(c(crossprod(Cmass.jack.val.fit$residuals))/length(Cmass.jack.val.fit$residuals)),
                                  max.val=max(meta(val.jack)$Cmass,na.rm=T),
                                  min.val=min(meta(val.jack)$Cmass,na.rm=T),
                                  bias=mean(Cmass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Cmass,na.rm=T))
  
  Nmass.jack.val.pred<-as.vector(predict(Nmass.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Nmass_CVmodel)[,,1])
  Nmass.jack.val.fit<-lm(Nmass.jack.val.pred~meta(val.jack)$Nmass)
  Nmass.jack.stats[[i]]<-c(R2=summary(Nmass.jack.val.fit)$r.squared,
                                  RMSE=sqrt(c(crossprod(Nmass.jack.val.fit$residuals))/length(Nmass.jack.val.fit$residuals)),
                                  max.val=max(meta(val.jack)$Nmass,na.rm=T),
                                  min.val=min(meta(val.jack)$Nmass,na.rm=T),
                                  bias=mean(Nmass.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Nmass,na.rm=T))
  
  LMA.jack.val.pred<-as.vector(predict(LMA.jack,newdata=as.matrix(val.jack),ncomp=ncomp_LMA_CVmodel)[,,1])
  LMA.jack.val.fit<-lm(LMA.jack.val.pred~meta(val.jack)$LMA)
  LMA.jack.stats[[i]]<-c(R2=summary(LMA.jack.val.fit)$r.squared,
                                 RMSE=sqrt(c(crossprod(LMA.jack.val.fit$residuals))/length(LMA.jack.val.fit$residuals)),
                                 max.val=max(meta(val.jack)$LMA,na.rm=T),
                                 min.val=min(meta(val.jack)$LMA,na.rm=T),
                                 bias=mean(LMA.jack.val.pred,na.rm=T)-mean(meta(val.jack)$LMA,na.rm=T))
  
  LDMC.jack.val.pred<-as.vector(predict(LDMC.jack,newdata=as.matrix(val.jack),ncomp=ncomp_LDMC_CVmodel)[,,1])
  LDMC.jack.val.fit<-lm(LDMC.jack.val.pred~meta(val.jack)$LDMC)
  LDMC.jack.stats[[i]]<-c(R2=summary(LDMC.jack.val.fit)$r.squared,
                                  RMSE=sqrt(c(crossprod(LDMC.jack.val.fit$residuals))/length(LDMC.jack.val.fit$residuals)),
                                  max.val=max(meta(val.jack)$LDMC,na.rm=T),
                                  min.val=min(meta(val.jack)$LDMC,na.rm=T),
                                  bias=mean(LDMC.jack.val.pred,na.rm=T)-mean(meta(val.jack)$LDMC,na.rm=T))

  EWT.jack.coefs[[i]]<-as.vector(coef(EWT.jack,ncomp=ncomp_EWT_CVmodel,intercept=TRUE))
  solubles_mass.jack.coefs[[i]]<-as.vector(coef(solubles_mass.jack,ncomp=ncomp_solubles_mass_CVmodel,intercept=TRUE))
  hemicellulose_mass.jack.coefs[[i]]<-as.vector(coef(hemicellulose_mass.jack,ncomp=ncomp_hemicellulose_mass_CVmodel,intercept=TRUE))
  cellulose_mass.jack.coefs[[i]]<-as.vector(coef(cellulose_mass.jack,ncomp=ncomp_cellulose_mass_CVmodel,intercept=TRUE))
  lignin_mass.jack.coefs[[i]]<-as.vector(coef(lignin_mass.jack,ncomp=ncomp_lignin_mass_CVmodel,intercept=TRUE))
  chlA_fresh.jack.coefs[[i]]<-as.vector(coef(chlA_fresh.jack,ncomp=ncomp_chlA_fresh_CVmodel,intercept=TRUE))
  chlB_fresh.jack.coefs[[i]]<-as.vector(coef(chlB_fresh.jack,ncomp=ncomp_chlB_fresh_CVmodel,intercept=TRUE))
  car_fresh.jack.coefs[[i]]<-as.vector(coef(car_fresh.jack,ncomp=ncomp_car_fresh_CVmodel,intercept=TRUE))
  Cmass.jack.coefs[[i]]<-as.vector(coef(Cmass.jack,ncomp=ncomp_Cmass_CVmodel,intercept=TRUE))
  Nmass.jack.coefs[[i]]<-as.vector(coef(Nmass.jack,ncomp=ncomp_Nmass_CVmodel,intercept=TRUE))
  LMA.jack.coefs[[i]]<-as.vector(coef(LMA.jack,ncomp=ncomp_LMA_CVmodel,intercept=TRUE))
  LDMC.jack.coefs[[i]]<-as.vector(coef(LDMC.jack,ncomp=ncomp_LDMC_CVmodel,intercept=TRUE))
  
}

EWT.jack.pred<-EWT.jack.pred<-apply.coefs(EWT.jack.coefs,as.matrix(spec.test))
EWT.jack.stat<-t(apply(EWT.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
EWT.jack.df<-data.frame(pred.mean=EWT.jack.stat[,1],
                        pred.low=EWT.jack.stat[,1]-1.96*EWT.jack.stat[,2],
                        pred.high=EWT.jack.stat[,1]+1.96*EWT.jack.stat[,2],
                        Measured=meta(spec.test)$EWT,
                        Species=meta(spec.test)$species,
                        Project=meta(spec.test)$project,
                        ID=meta(spec.test)$sample_id)

EWT.val.plot<-ggplot(EWT.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.068),ylim=c(0,0.068))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (cm)",x="Predicted EWT (cm)")+
  guides(color=F)

solubles_mass.jack.pred<-apply.coefs(solubles_mass.jack.coefs,as.matrix(spec.test))
solubles_mass.jack.stat<-t(apply(solubles_mass.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
solubles_mass.jack.df<-data.frame(pred.mean=solubles_mass.jack.stat[,1],
                            pred.low=solubles_mass.jack.stat[,1]-1.96*solubles_mass.jack.stat[,2],
                            pred.high=solubles_mass.jack.stat[,1]+1.96*solubles_mass.jack.stat[,2],
                            Measured=meta(spec.test)$solubles_mass,
                            Species=meta(spec.test)$species,
                            Project=meta(spec.test)$project,
                            ID=meta(spec.test)$sample_id)

solubles_mass.val.plot<-ggplot(solubles_mass.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(30,90),ylim=c(30,90))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)

hemicellulose_mass.jack.pred<-apply.coefs(hemicellulose_mass.jack.coefs,as.matrix(spec.test))
hemicellulose_mass.jack.stat<-t(apply(hemicellulose_mass.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
hemicellulose_mass.jack.df<-data.frame(pred.mean=hemicellulose_mass.jack.stat[,1],
                            pred.low=hemicellulose_mass.jack.stat[,1]-1.96*hemicellulose_mass.jack.stat[,2],
                            pred.high=hemicellulose_mass.jack.stat[,1]+1.96*hemicellulose_mass.jack.stat[,2],
                            Measured=meta(spec.test)$hemicellulose_mass,
                            Species=meta(spec.test)$species,
                            Project=meta(spec.test)$project,
                            ID=meta(spec.test)$sample_id)

hemicellulose_mass.val.plot<-ggplot(hemicellulose_mass.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,35),ylim=c(0,35))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)

cellulose_mass.jack.pred<-apply.coefs(cellulose_mass.jack.coefs,as.matrix(spec.test))
cellulose_mass.jack.stat<-t(apply(cellulose_mass.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
cellulose_mass.jack.df<-data.frame(pred.mean=cellulose_mass.jack.stat[,1],
                            pred.low=cellulose_mass.jack.stat[,1]-1.96*cellulose_mass.jack.stat[,2],
                            pred.high=cellulose_mass.jack.stat[,1]+1.96*cellulose_mass.jack.stat[,2],
                            Measured=meta(spec.test)$cellulose_mass,
                            Species=meta(spec.test)$species,
                            Project=meta(spec.test)$project,
                            ID=meta(spec.test)$sample_id)

cellulose_mass.val.plot<-ggplot(cellulose_mass.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,35),ylim=c(0,35))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)

lignin_mass.jack.pred<-apply.coefs(lignin_mass.jack.coefs,as.matrix(spec.test))
lignin_mass.jack.stat<-t(apply(lignin_mass.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
lignin_mass.jack.df<-data.frame(pred.mean=lignin_mass.jack.stat[,1],
                            pred.low=lignin_mass.jack.stat[,1]-1.96*lignin_mass.jack.stat[,2],
                            pred.high=lignin_mass.jack.stat[,1]+1.96*lignin_mass.jack.stat[,2],
                            Measured=meta(spec.test)$lignin_mass,
                            Species=meta(spec.test)$species,
                            Project=meta(spec.test)$project,
                            ID=meta(spec.test)$sample_id)

lignin_mass.val.plot<-ggplot(lignin_mass.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-5,20),ylim=c(-5,20))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)

chlA_fresh.jack.pred<-apply.coefs(chlA_fresh.jack.coefs,as.matrix(spec.test))
chlA_fresh.jack.stat<-t(apply(chlA_fresh.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
chlA_fresh.jack.df<-data.frame(pred.mean=chlA_fresh.jack.stat[,1],
                                pred.low=chlA_fresh.jack.stat[,1]-1.96*chlA_fresh.jack.stat[,2],
                                pred.high=chlA_fresh.jack.stat[,1]+1.96*chlA_fresh.jack.stat[,2],
                                Measured=meta(spec.test)$chlA_fresh,
                                Species=meta(spec.test)$species,
                                Project=meta(spec.test)$project,
                                ID=meta(spec.test)$sample_id)

chlA_fresh.val.plot<-ggplot(chlA_fresh.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-5,20),ylim=c(-5,20))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl a (mg g"^-1*")"),x=expression("Predicted Chl b (mg g"^-1*")"))+
  guides(color=F)

chlB_fresh.jack.pred<-apply.coefs(chlB_fresh.jack.coefs,as.matrix(spec.test))
chlB_fresh.jack.stat<-t(apply(chlB_fresh.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
chlB_fresh.jack.df<-data.frame(pred.mean=chlB_fresh.jack.stat[,1],
                                pred.low=chlB_fresh.jack.stat[,1]-1.96*chlB_fresh.jack.stat[,2],
                                pred.high=chlB_fresh.jack.stat[,1]+1.96*chlB_fresh.jack.stat[,2],
                                Measured=meta(spec.test)$chlB_fresh,
                                Species=meta(spec.test)$species,
                                Project=meta(spec.test)$project,
                                ID=meta(spec.test)$sample_id)

chlB_fresh.val.plot<-ggplot(chlB_fresh.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-5,20),ylim=c(-5,20))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl b (mg g"^-1*")"),x=expression("Predicted Chl b (mg g"^-1*")"))+
  guides(color=F)

car_fresh.jack.pred<-apply.coefs(car_fresh.jack.coefs,as.matrix(spec.test))
car_fresh.jack.stat<-t(apply(car_fresh.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
car_fresh.jack.df<-data.frame(pred.mean=car_fresh.jack.stat[,1],
                                pred.low=car_fresh.jack.stat[,1]-1.96*car_fresh.jack.stat[,2],
                                pred.high=car_fresh.jack.stat[,1]+1.96*car_fresh.jack.stat[,2],
                                Measured=meta(spec.test)$car_fresh,
                                Species=meta(spec.test)$species,
                                Project=meta(spec.test)$project,
                                ID=meta(spec.test)$sample_id)

car_fresh.val.plot<-ggplot(car_fresh.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-5,20),ylim=c(-5,20))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),x=expression("Predicted carotenoids (mg g"^-1*")"))+
  guides(color=F)

Cmass.jack.pred<-apply.coefs(Cmass.jack.coefs,as.matrix(spec.test))
Cmass.jack.stat<-t(apply(Cmass.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
Cmass.jack.df<-data.frame(pred.mean=Cmass.jack.stat[,1],
                         pred.low=Cmass.jack.stat[,1]-1.96*Cmass.jack.stat[,2],
                         pred.high=Cmass.jack.stat[,1]+1.96*Cmass.jack.stat[,2],
                         Measured=meta(spec.test)$Cmass,
                         Species=meta(spec.test)$species,
                         Project=meta(spec.test)$project,
                         ID=meta(spec.test)$sample_id)

Cmass.val.plot<-ggplot(Cmass.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(36,55),ylim=c(36,55))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured C"[mass]*" (%)"),x=expression("Predicted C"[mass]*" (%)"))+
  guides(color=F)

Nmass.jack.pred<-apply.coefs(Nmass.jack.coefs,as.matrix(spec.test))
Nmass.jack.stat<-t(apply(Nmass.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
Nmass.jack.df<-data.frame(pred.mean=Nmass.jack.stat[,1],
                         pred.low=Nmass.jack.stat[,1]-1.96*Nmass.jack.stat[,2],
                         pred.high=Nmass.jack.stat[,1]+1.96*Nmass.jack.stat[,2],
                         Measured=meta(spec.test)$Nmass,
                         Species=meta(spec.test)$species,
                         Project=meta(spec.test)$project,
                         ID=meta(spec.test)$sample_id)

Nmass.val.plot<-ggplot(Nmass.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,6),ylim=c(0,6))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured N"[mass]*" (%)"),x=expression("Predicted N"[mass]*" (%)"))+  guides(color=F)

LMA.jack.pred<-apply.coefs(LMA.jack.coefs,as.matrix(spec.test))
LMA.jack.stat<-t(apply(LMA.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
LMA.jack.df<-data.frame(pred.mean=LMA.jack.stat[,1],
                        pred.low=LMA.jack.stat[,1]-1.96*LMA.jack.stat[,2],
                        pred.high=LMA.jack.stat[,1]+1.96*LMA.jack.stat[,2],
                        Measured=meta(spec.test)$LMA,
                        Species=meta(spec.test)$species,
                        Project=meta(spec.test)$project,
                        ID=meta(spec.test)$sample_id)

LMA.val.plot<-ggplot(LMA.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.025,0.35),ylim=c(-0.025,0.35))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured LMA (kg m"^-2*")"),x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)

LDMC.jack.pred<-apply.coefs(LDMC.jack.coefs,as.matrix(spec.test))
LDMC.jack.stat<-t(apply(LDMC.jack.pred,1,function(obs) c(mean(obs),sd(obs))))
LDMC.jack.df<-data.frame(pred.mean=LDMC.jack.stat[,1],
                         pred.low=LDMC.jack.stat[,1]-1.96*LDMC.jack.stat[,2],
                         pred.high=LDMC.jack.stat[,1]+1.96*LDMC.jack.stat[,2],
                         Measured=meta(spec.test)$LDMC,
                         Species=meta(spec.test)$species,
                         Project=meta(spec.test)$project,
                         ID=meta(spec.test)$sample_id)

LDMC.val.plot<-ggplot(LDMC.jack.df,aes(y=Measured,x=pred.mean,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray",alpha=0.7)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,600),ylim=c(0,600))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)

all.jack.coef.list<-list(EWT=EWT.jack.coefs,
                         solubles_mass=solubles_mass.jack.coefs,
                         hemicellulose_mass=hemicellulose_mass.jack.coefs,
                         cellulose_mass=cellulose_mass.jack.coefs,
                         lignin_mass=lignin_mass.jack.coefs,
                         chlA_fresh=chlA_fresh.jack.coefs,
                         chlB_fresh=chlB_fresh.jack.coefs,
                         car_fresh=car_fresh.jack.coefs,
                         Nmass=Nmass.jack.coefs,
                         Cmass=Cmass.jack.coefs,
                         LDMC=LDMC.jack.coefs,
                         LMA=LMA.jack.coefs)
saveRDS(all.jack.coef.list,"SavedResults/all_jack_coefs_list.rds")

all.jack.df.list<-list(EWT=EWT.jack.df,
                       solubles_mass=solubles_mass.jack.df,
                       hemicellulose_mass=hemicellulose_mass.jack.df,
                       cellulose_mass=cellulose_mass.jack.df,
                       lignin_mass=lignin_mass.jack.df,
                       chlA_fresh=chlA_fresh.jack.df,
                       chlB_fresh=chlB_fresh.jack.df,
                       car_fresh=car_fresh.jack.df,
                       Nmass=Nmass.jack.df,
                       Cmass=Cmass.jack.df,
                       LDMC=LDMC.jack.df,
                       LMA=LMA.jack.df)
saveRDS(all.jack.df.list,"SavedResults/all_jack_df_list.rds")

pdf("Images/val_plots.pdf",width = 10,height = 9)
EWT.val.plot
solubles_mass.val.plot
hemicellulose_mass.val.plot
cellulose_mass.val.plot
lignin_mass.val.plot
Nmass.val.plot
Cmass.val.plot
LMA.val.plot
LDMC.val.plot
dev.off()
