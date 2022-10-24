setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(pls)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

## download VIP.R from Bjorn-Helge Mevik's website:
## https://mevik.net/work/software/pls.html
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

## sqrt-transform

LMA_sqrt_CVmodel<-plsr(sqrt(meta(ref.train)$LMA)~as.matrix(ref.train),
                       ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_sqrt_CVmodel <- selectNcomp(LMA_sqrt_CVmodel, method = "onesigma", plot = FALSE)
LMA_sqrt_valid <- which(!is.na(meta(ref.train)$LMA))
LMA_sqrt_pred<-data.frame(ID=meta(ref.train)$sample_id[LMA_sqrt_valid],
                          Species=meta(ref.train)$species[LMA_sqrt_valid],
                          Project=meta(ref.train)$project[LMA_sqrt_valid],
                          measured=meta(ref.train)$LMA[LMA_sqrt_valid],
                          val_pred=LMA_sqrt_CVmodel$validation$pred[,,ncomp_LMA_sqrt_CVmodel]^2)

Nmass_sqrt_CVmodel<-plsr(sqrt(meta(ref.train)$Nmass)~as.matrix(ref.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Nmass_sqrt_CVmodel <- selectNcomp(Nmass_sqrt_CVmodel, method = "onesigma", plot = FALSE)
Nmass_sqrt_valid <- which(!is.na(meta(ref.train)$Nmass))
Nmass_sqrt_pred<-data.frame(ID=meta(ref.train)$sample_id[Nmass_sqrt_valid],
                            Species=meta(ref.train)$species[Nmass_sqrt_valid],
                            Project=meta(ref.train)$project[Nmass_sqrt_valid],
                            measured=meta(ref.train)$Nmass[Nmass_sqrt_valid],
                            val_pred=Nmass_sqrt_CVmodel$validation$pred[,,ncomp_Nmass_sqrt_CVmodel]^2)

Cmass_sqrt_CVmodel<-plsr(sqrt(meta(ref.train)$Cmass)~as.matrix(ref.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cmass_sqrt_CVmodel <- selectNcomp(Cmass_sqrt_CVmodel, method = "onesigma", plot = FALSE)
Cmass_sqrt_valid <- which(!is.na(meta(ref.train)$Cmass))
Cmass_sqrt_pred<-data.frame(ID=meta(ref.train)$sample_id[Cmass_sqrt_valid],
                            Species=meta(ref.train)$species[Cmass_sqrt_valid],
                            Project=meta(ref.train)$project[Cmass_sqrt_valid],
                            measured=meta(ref.train)$Cmass[Cmass_sqrt_valid],
                            val_pred=Cmass_sqrt_CVmodel$validation$pred[,,ncomp_Cmass_sqrt_CVmodel]^2)

chlA_mass_sqrt_CVmodel<-plsr(sqrt(meta(ref.train)$chlA_mass)~as.matrix(ref.train),
                             ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_mass_sqrt_CVmodel <- selectNcomp(chlA_mass_sqrt_CVmodel, method = "onesigma", plot = FALSE)
chlA_mass_sqrt_valid <- which(!is.na(meta(ref.train)$chlA_mass))
chlA_mass_sqrt_pred<-data.frame(ID=meta(ref.train)$sample_id[chlA_mass_sqrt_valid],
                                Species=meta(ref.train)$species[chlA_mass_sqrt_valid],
                                Project=meta(ref.train)$project[chlA_mass_sqrt_valid],
                                measured=meta(ref.train)$chlA_mass[chlA_mass_sqrt_valid],
                                val_pred=chlA_mass_sqrt_CVmodel$validation$pred[,,ncomp_chlA_mass_sqrt_CVmodel]^2)

K_mass_sqrt_CVmodel<-plsr(sqrt(meta(ref.train)$K_mass)~as.matrix(ref.train),
                          ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_mass_sqrt_CVmodel <- selectNcomp(K_mass_sqrt_CVmodel, method = "onesigma", plot = FALSE)
K_mass_sqrt_valid <- which(!is.na(meta(ref.train)$K_mass))
K_mass_sqrt_pred<-data.frame(ID=meta(ref.train)$sample_id[K_mass_sqrt_valid],
                             Species=meta(ref.train)$species[K_mass_sqrt_valid],
                             Project=meta(ref.train)$project[K_mass_sqrt_valid],
                             measured=meta(ref.train)$K_mass[K_mass_sqrt_valid],
                             val_pred=K_mass_sqrt_CVmodel$validation$pred[,,ncomp_K_mass_sqrt_CVmodel]^2)

## log-transform

LMA_log_CVmodel<-plsr(log(meta(ref.train)$LMA)~as.matrix(ref.train),
                      ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_LMA_log_CVmodel <- selectNcomp(LMA_log_CVmodel, method = "onesigma", plot = FALSE)
LMA_log_valid <- which(!is.na(meta(ref.train)$LMA))
LMA_log_pred<-data.frame(ID=meta(ref.train)$sample_id[LMA_log_valid],
                         Species=meta(ref.train)$species[LMA_log_valid],
                         Project=meta(ref.train)$project[LMA_log_valid],
                         measured=meta(ref.train)$LMA[LMA_log_valid],
                         val_pred=exp(LMA_log_CVmodel$validation$pred[,,ncomp_LMA_log_CVmodel]))

Nmass_log_CVmodel<-plsr(log(meta(ref.train)$Nmass)~as.matrix(ref.train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Nmass_log_CVmodel <- selectNcomp(Nmass_log_CVmodel, method = "onesigma", plot = FALSE)
Nmass_log_valid <- which(!is.na(meta(ref.train)$Nmass))
Nmass_log_pred<-data.frame(ID=meta(ref.train)$sample_id[Nmass_log_valid],
                           Species=meta(ref.train)$species[Nmass_log_valid],
                           Project=meta(ref.train)$project[Nmass_log_valid],
                           measured=meta(ref.train)$Nmass[Nmass_log_valid],
                           val_pred=exp(Nmass_log_CVmodel$validation$pred[,,ncomp_Nmass_log_CVmodel]))

Cmass_log_CVmodel<-plsr(log(meta(ref.train)$Cmass)~as.matrix(ref.train),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cmass_log_CVmodel <- selectNcomp(Cmass_log_CVmodel, method = "onesigma", plot = FALSE)
Cmass_log_valid <- which(!is.na(meta(ref.train)$Cmass))
Cmass_log_pred<-data.frame(ID=meta(ref.train)$sample_id[Cmass_log_valid],
                           Species=meta(ref.train)$species[Cmass_log_valid],
                           Project=meta(ref.train)$project[Cmass_log_valid],
                           measured=meta(ref.train)$Cmass[Cmass_log_valid],
                           val_pred=exp(Cmass_log_CVmodel$validation$pred[,,ncomp_Cmass_log_CVmodel]))

chlA_mass_log_CVmodel<-plsr(log(meta(ref.train)$chlA_mass)~as.matrix(ref.train),
                            ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_chlA_mass_log_CVmodel <- selectNcomp(chlA_mass_log_CVmodel, method = "onesigma", plot = FALSE)
chlA_mass_log_valid <- which(!is.na(meta(ref.train)$chlA_mass))
chlA_mass_log_pred<-data.frame(ID=meta(ref.train)$sample_id[chlA_mass_log_valid],
                               Species=meta(ref.train)$species[chlA_mass_log_valid],
                               Project=meta(ref.train)$project[chlA_mass_log_valid],
                               measured=meta(ref.train)$chlA_mass[chlA_mass_log_valid],
                               val_pred=exp(chlA_mass_log_CVmodel$validation$pred[,,ncomp_chlA_mass_log_CVmodel]))

K_mass_log_CVmodel<-plsr(log(meta(ref.train)$K_mass)~as.matrix(ref.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_K_mass_log_CVmodel <- selectNcomp(K_mass_log_CVmodel, method = "onesigma", plot = FALSE)
K_mass_log_valid <- which(!is.na(meta(ref.train)$K_mass))
K_mass_log_pred<-data.frame(ID=meta(ref.train)$sample_id[K_mass_log_valid],
                            Species=meta(ref.train)$species[K_mass_log_valid],
                            Project=meta(ref.train)$project[K_mass_log_valid],
                            measured=meta(ref.train)$K_mass[K_mass_log_valid],
                            val_pred=exp(K_mass_log_CVmodel$validation$pred[,,ncomp_K_mass_log_CVmodel]))

#######################################################
## jackknife analyses

chlA_mass.sqrt.jack.coefs<-list()
Nmass.sqrt.jack.coefs<-list()
Cmass.sqrt.jack.coefs<-list()
LMA.sqrt.jack.coefs<-list()
K_mass.sqrt.jack.coefs<-list()

chlA_mass.log.jack.coefs<-list()
Nmass.log.jack.coefs<-list()
Cmass.log.jack.coefs<-list()
LMA.log.jack.coefs<-list()
K_mass.log.jack.coefs<-list()

chlA_mass.sqrt.jack.stats<-list()
Nmass.sqrt.jack.stats<-list()
Cmass.sqrt.jack.stats<-list()
LMA.sqrt.jack.stats<-list()
K_mass.sqrt.jack.stats<-list()

chlA_mass.log.jack.stats<-list()
Nmass.log.jack.stats<-list()
Cmass.log.jack.stats<-list()
LMA.log.jack.stats<-list()
K_mass.log.jack.stats<-list()

nreps<-100

for(i in 1:nreps){
  print(i)
  
  n.cal.spec<-nrow(ref.train)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-ref.train[train.jack]
  val.jack<-ref.train[test.jack]
  
  chlA_mass.sqrt.jack<-plsr(sqrt(meta(calib.jack)$chlA_mass)~as.matrix(calib.jack),
                            ncomp=30,method = "oscorespls",validation="none")
  Nmass.sqrt.jack<-plsr(sqrt(meta(calib.jack)$Nmass)~as.matrix(calib.jack),
                        ncomp=30,method = "oscorespls",validation="none")
  Cmass.sqrt.jack<-plsr(sqrt(meta(calib.jack)$Cmass)~as.matrix(calib.jack),
                        ncomp=30,method = "oscorespls",validation="none")
  LMA.sqrt.jack<-plsr(sqrt(meta(calib.jack)$LMA)~as.matrix(calib.jack),
                      ncomp=30,method = "oscorespls",validation="none")
  K_mass.sqrt.jack<-plsr(sqrt(meta(calib.jack)$K_mass)~as.matrix(calib.jack),
                         ncomp=30,method = "oscorespls",validation="none")
  
  chlA_mass.log.jack<-plsr(log(meta(calib.jack)$chlA_mass)~as.matrix(calib.jack),
                           ncomp=30,method = "oscorespls",validation="none")
  Nmass.log.jack<-plsr(log(meta(calib.jack)$Nmass)~as.matrix(calib.jack),
                       ncomp=30,method = "oscorespls",validation="none")
  Cmass.log.jack<-plsr(log(meta(calib.jack)$Cmass)~as.matrix(calib.jack),
                       ncomp=30,method = "oscorespls",validation="none")
  LMA.log.jack<-plsr(log(meta(calib.jack)$LMA)~as.matrix(calib.jack),
                     ncomp=30,method = "oscorespls",validation="none")
  K_mass.log.jack<-plsr(log(meta(calib.jack)$K_mass)~as.matrix(calib.jack),
                        ncomp=30,method = "oscorespls",validation="none")
  
  chlA_mass.sqrt.jack.val.pred<-as.vector(predict(chlA_mass.sqrt.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlA_mass_sqrt_CVmodel)[,,1])^2
  chlA_mass.sqrt.jack.val.fit<-lm(chlA_mass.sqrt.jack.val.pred~meta(val.jack)$chlA_mass)
  chlA_mass.sqrt.jack.stats[[i]]<-c(R2=summary(chlA_mass.sqrt.jack.val.fit)$r.squared,
                                    RMSE=RMSD(meta(val.jack)$chlA_mass,chlA_mass.sqrt.jack.val.pred),
                                    perRMSE=percentRMSD(meta(val.jack)$chlA_mass,chlA_mass.sqrt.jack.val.pred,0.025,0.975),
                                    bias=mean(chlA_mass.sqrt.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlA_mass,na.rm=T))
  
  Nmass.sqrt.jack.val.pred<-as.vector(predict(Nmass.sqrt.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Nmass_sqrt_CVmodel)[,,1])^2
  Nmass.sqrt.jack.val.fit<-lm(Nmass.sqrt.jack.val.pred~meta(val.jack)$Nmass)
  Nmass.sqrt.jack.stats[[i]]<-c(R2=summary(Nmass.sqrt.jack.val.fit)$r.squared,
                                RMSE=RMSD(meta(val.jack)$Nmass,Nmass.sqrt.jack.val.pred),
                                perRMSE=percentRMSD(meta(val.jack)$Nmass,Nmass.sqrt.jack.val.pred,0.025,0.975),
                                bias=mean(Nmass.sqrt.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Nmass,na.rm=T))
  
  Cmass.sqrt.jack.val.pred<-as.vector(predict(Cmass.sqrt.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Cmass_sqrt_CVmodel)[,,1])^2
  Cmass.sqrt.jack.val.fit<-lm(Cmass.sqrt.jack.val.pred~meta(val.jack)$Cmass)
  Cmass.sqrt.jack.stats[[i]]<-c(R2=summary(Cmass.sqrt.jack.val.fit)$r.squared,
                                RMSE=RMSD(meta(val.jack)$Cmass,Cmass.sqrt.jack.val.pred),
                                perRMSE=percentRMSD(meta(val.jack)$Cmass,Cmass.sqrt.jack.val.pred,0.025,0.975),
                                bias=mean(Cmass.sqrt.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Cmass,na.rm=T))
  
  LMA.sqrt.jack.val.pred<-as.vector(predict(LMA.sqrt.jack,newdata=as.matrix(val.jack),ncomp=ncomp_LMA_sqrt_CVmodel)[,,1])^2
  LMA.sqrt.jack.val.fit<-lm(LMA.sqrt.jack.val.pred~meta(val.jack)$LMA)
  LMA.sqrt.jack.stats[[i]]<-c(R2=summary(LMA.sqrt.jack.val.fit)$r.squared,
                              RMSE=RMSD(meta(val.jack)$LMA,LMA.sqrt.jack.val.pred),
                              perRMSE=percentRMSD(meta(val.jack)$LMA,LMA.sqrt.jack.val.pred,0.025,0.975),
                              bias=mean(LMA.sqrt.jack.val.pred,na.rm=T)-mean(meta(val.jack)$LMA,na.rm=T))
  
  K_mass.sqrt.jack.val.pred<-as.vector(predict(K_mass.sqrt.jack,newdata=as.matrix(val.jack),ncomp=ncomp_K_mass_sqrt_CVmodel)[,,1])^2
  K_mass.sqrt.jack.val.fit<-lm(K_mass.sqrt.jack.val.pred~meta(val.jack)$K_mass)
  K_mass.sqrt.jack.stats[[i]]<-c(R2=summary(K_mass.sqrt.jack.val.fit)$r.squared,
                                 RMSE=RMSD(meta(val.jack)$K_mass,K_mass.sqrt.jack.val.pred),
                                 perRMSE=percentRMSD(meta(val.jack)$K_mass,K_mass.sqrt.jack.val.pred,0.025,0.975),
                                 bias=mean(K_mass.sqrt.jack.val.pred,na.rm=T)-mean(meta(val.jack)$K_mass,na.rm=T))
  
  chlA_mass.log.jack.val.pred<-exp(as.vector(predict(chlA_mass.log.jack,newdata=as.matrix(val.jack),ncomp=ncomp_chlA_mass_log_CVmodel)[,,1]))
  chlA_mass.log.jack.val.fit<-lm(chlA_mass.log.jack.val.pred~meta(val.jack)$chlA_mass)
  chlA_mass.log.jack.stats[[i]]<-c(R2=summary(chlA_mass.log.jack.val.fit)$r.squared,
                                   RMSE=RMSD(meta(val.jack)$chlA_mass,chlA_mass.log.jack.val.pred),
                                   perRMSE=percentRMSD(meta(val.jack)$chlA_mass,chlA_mass.log.jack.val.pred,0.025,0.975),
                                   bias=mean(chlA_mass.log.jack.val.pred,na.rm=T)-mean(meta(val.jack)$chlA_mass,na.rm=T))
  
  Nmass.log.jack.val.pred<-exp(as.vector(predict(Nmass.log.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Nmass_log_CVmodel)[,,1]))
  Nmass.log.jack.val.fit<-lm(Nmass.log.jack.val.pred~meta(val.jack)$Nmass)
  Nmass.log.jack.stats[[i]]<-c(R2=summary(Nmass.log.jack.val.fit)$r.squared,
                               RMSE=RMSD(meta(val.jack)$Nmass,Nmass.log.jack.val.pred),
                               perRMSE=percentRMSD(meta(val.jack)$Nmass,Nmass.log.jack.val.pred,0.025,0.975),
                               bias=mean(Nmass.log.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Nmass,na.rm=T))
  
  Cmass.log.jack.val.pred<-exp(as.vector(predict(Cmass.log.jack,newdata=as.matrix(val.jack),ncomp=ncomp_Cmass_log_CVmodel)[,,1]))
  Cmass.log.jack.val.fit<-lm(Cmass.log.jack.val.pred~meta(val.jack)$Cmass)
  Cmass.log.jack.stats[[i]]<-c(R2=summary(Cmass.log.jack.val.fit)$r.squared,
                               RMSE=RMSD(meta(val.jack)$Cmass,Cmass.log.jack.val.pred),
                               perRMSE=percentRMSD(meta(val.jack)$Cmass,Cmass.log.jack.val.pred,0.025,0.975),
                               bias=mean(Cmass.log.jack.val.pred,na.rm=T)-mean(meta(val.jack)$Cmass,na.rm=T))
  
  LMA.log.jack.val.pred<-exp(as.vector(predict(LMA.log.jack,newdata=as.matrix(val.jack),ncomp=ncomp_LMA_log_CVmodel)[,,1]))
  LMA.log.jack.val.fit<-lm(LMA.log.jack.val.pred~meta(val.jack)$LMA)
  LMA.log.jack.stats[[i]]<-c(R2=summary(LMA.log.jack.val.fit)$r.squared,
                             RMSE=RMSD(meta(val.jack)$LMA,LMA.log.jack.val.pred),
                             perRMSE=percentRMSD(meta(val.jack)$LMA,LMA.log.jack.val.pred,0.025,0.975),
                             bias=mean(LMA.log.jack.val.pred,na.rm=T)-mean(meta(val.jack)$LMA,na.rm=T))
  
  K_mass.log.jack.val.pred<-exp(as.vector(predict(K_mass.log.jack,newdata=as.matrix(val.jack),ncomp=ncomp_K_mass_log_CVmodel)[,,1]))
  K_mass.log.jack.val.fit<-lm(K_mass.log.jack.val.pred~meta(val.jack)$K_mass)
  K_mass.log.jack.stats[[i]]<-c(R2=summary(K_mass.log.jack.val.fit)$r.squared,
                                RMSE=RMSD(meta(val.jack)$K_mass,K_mass.log.jack.val.pred),
                                perRMSE=percentRMSD(meta(val.jack)$K_mass,K_mass.log.jack.val.pred,0.025,0.975),
                                bias=mean(K_mass.log.jack.val.pred,na.rm=T)-mean(meta(val.jack)$K_mass,na.rm=T))
  
  chlA_mass.sqrt.jack.coefs[[i]]<-as.vector(coef(chlA_mass.sqrt.jack,ncomp=ncomp_chlA_mass_sqrt_CVmodel,intercept=TRUE))
  Nmass.sqrt.jack.coefs[[i]]<-as.vector(coef(Nmass.sqrt.jack,ncomp=ncomp_Nmass_sqrt_CVmodel,intercept=TRUE))
  Cmass.sqrt.jack.coefs[[i]]<-as.vector(coef(Cmass.sqrt.jack,ncomp=ncomp_Cmass_sqrt_CVmodel,intercept=TRUE))
  LMA.sqrt.jack.coefs[[i]]<-as.vector(coef(LMA.sqrt.jack,ncomp=ncomp_LMA_sqrt_CVmodel,intercept=TRUE))
  K_mass.sqrt.jack.coefs[[i]]<-as.vector(coef(K_mass.sqrt.jack,ncomp=ncomp_K_mass_sqrt_CVmodel,intercept=TRUE))
  
  chlA_mass.log.jack.coefs[[i]]<-as.vector(coef(chlA_mass.log.jack,ncomp=ncomp_chlA_mass_log_CVmodel,intercept=TRUE))
  Nmass.log.jack.coefs[[i]]<-as.vector(coef(Nmass.log.jack,ncomp=ncomp_Nmass_log_CVmodel,intercept=TRUE))
  Cmass.log.jack.coefs[[i]]<-as.vector(coef(Cmass.log.jack,ncomp=ncomp_Cmass_log_CVmodel,intercept=TRUE))
  LMA.log.jack.coefs[[i]]<-as.vector(coef(LMA.log.jack,ncomp=ncomp_LMA_log_CVmodel,intercept=TRUE))
  K_mass.log.jack.coefs[[i]]<-as.vector(coef(K_mass.log.jack,ncomp=ncomp_K_mass_log_CVmodel,intercept=TRUE))
  
}

chlA_mass.sqrt.jack.pred<-apply.coefs(chlA_mass.sqrt.jack.coefs,as.matrix(ref.test))^2
chlA_mass.sqrt.jack.stat<-t(apply(chlA_mass.sqrt.jack.pred,1,
                                  function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_mass.sqrt.jack.df<-data.frame(pred.mean=chlA_mass.sqrt.jack.stat[,1],
                                   pred.low=chlA_mass.sqrt.jack.stat[,2],
                                   pred.high=chlA_mass.sqrt.jack.stat[,3],
                                   Measured=meta(ref.test)$chlA_mass,
                                   ncomp=ncomp_chlA_mass_sqrt_CVmodel,
                                   Species=meta(ref.test)$species,
                                   Project=meta(ref.test)$project,
                                   functional.group=meta(ref.test)$functional.group,
                                   ID=meta(ref.test)$sample_id)

Nmass.sqrt.jack.pred<-apply.coefs(Nmass.sqrt.jack.coefs,as.matrix(ref.test))^2
Nmass.sqrt.jack.stat<-t(apply(Nmass.sqrt.jack.pred,1,
                              function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Nmass.sqrt.jack.df<-data.frame(pred.mean=Nmass.sqrt.jack.stat[,1],
                               pred.low=Nmass.sqrt.jack.stat[,2],
                               pred.high=Nmass.sqrt.jack.stat[,3],
                               Measured=meta(ref.test)$Nmass,
                               ncomp=ncomp_Nmass_sqrt_CVmodel,
                               Species=meta(ref.test)$species,
                               Project=meta(ref.test)$project,
                               functional.group=meta(ref.test)$functional.group,
                               ID=meta(ref.test)$sample_id)

Cmass.sqrt.jack.pred<-apply.coefs(Cmass.sqrt.jack.coefs,as.matrix(ref.test))^2
Cmass.sqrt.jack.stat<-t(apply(Cmass.sqrt.jack.pred,1,
                              function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cmass.sqrt.jack.df<-data.frame(pred.mean=Cmass.sqrt.jack.stat[,1],
                               pred.low=Cmass.sqrt.jack.stat[,2],
                               pred.high=Cmass.sqrt.jack.stat[,3],
                               Measured=meta(ref.test)$Cmass,
                               ncomp=ncomp_Cmass_sqrt_CVmodel,
                               Species=meta(ref.test)$species,
                               Project=meta(ref.test)$project,
                               functional.group=meta(ref.test)$functional.group,
                               ID=meta(ref.test)$sample_id)

LMA.sqrt.jack.pred<-apply.coefs(LMA.sqrt.jack.coefs,as.matrix(ref.test))^2
LMA.sqrt.jack.stat<-t(apply(LMA.sqrt.jack.pred,1,
                            function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA.sqrt.jack.df<-data.frame(pred.mean=LMA.sqrt.jack.stat[,1],
                             pred.low=LMA.sqrt.jack.stat[,2],
                             pred.high=LMA.sqrt.jack.stat[,3],
                             Measured=meta(ref.test)$LMA,
                             ncomp=ncomp_LMA_sqrt_CVmodel,
                             Species=meta(ref.test)$species,
                             Project=meta(ref.test)$project,
                             functional.group=meta(ref.test)$functional.group,
                             ID=meta(ref.test)$sample_id)

K_mass.sqrt.jack.pred<-apply.coefs(K_mass.sqrt.jack.coefs,as.matrix(ref.test))^2
K_mass.sqrt.jack.stat<-t(apply(K_mass.sqrt.jack.pred,1,
                               function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_mass.sqrt.jack.df<-data.frame(pred.mean=K_mass.sqrt.jack.stat[,1],
                                pred.low=K_mass.sqrt.jack.stat[,2],
                                pred.high=K_mass.sqrt.jack.stat[,3],
                                Measured=meta(ref.test)$K_mass,
                                ncomp=ncomp_K_mass_sqrt_CVmodel,
                                Species=meta(ref.test)$species,
                                Project=meta(ref.test)$project,
                                functional.group=meta(ref.test)$functional.group,
                                ID=meta(ref.test)$sample_id)

chlA_mass.log.jack.pred<-exp(apply.coefs(chlA_mass.log.jack.coefs,as.matrix(ref.test)))
chlA_mass.log.jack.stat<-t(apply(chlA_mass.log.jack.pred,1,
                                 function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
chlA_mass.log.jack.df<-data.frame(pred.mean=chlA_mass.log.jack.stat[,1],
                                  pred.low=chlA_mass.log.jack.stat[,2],
                                  pred.high=chlA_mass.log.jack.stat[,3],
                                  Measured=meta(ref.test)$chlA_mass,
                                  ncomp=ncomp_chlA_mass_log_CVmodel,
                                  Species=meta(ref.test)$species,
                                  Project=meta(ref.test)$project,
                                  functional.group=meta(ref.test)$functional.group,
                                  ID=meta(ref.test)$sample_id)

Nmass.log.jack.pred<-exp(apply.coefs(Nmass.log.jack.coefs,as.matrix(ref.test)))
Nmass.log.jack.stat<-t(apply(Nmass.log.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Nmass.log.jack.df<-data.frame(pred.mean=Nmass.log.jack.stat[,1],
                              pred.low=Nmass.log.jack.stat[,2],
                              pred.high=Nmass.log.jack.stat[,3],
                              Measured=meta(ref.test)$Nmass,
                              ncomp=ncomp_Nmass_log_CVmodel,
                              Species=meta(ref.test)$species,
                              Project=meta(ref.test)$project,
                              functional.group=meta(ref.test)$functional.group,
                              ID=meta(ref.test)$sample_id)

Cmass.log.jack.pred<-exp(apply.coefs(Cmass.log.jack.coefs,as.matrix(ref.test)))
Cmass.log.jack.stat<-t(apply(Cmass.log.jack.pred,1,
                             function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
Cmass.log.jack.df<-data.frame(pred.mean=Cmass.log.jack.stat[,1],
                              pred.low=Cmass.log.jack.stat[,2],
                              pred.high=Cmass.log.jack.stat[,3],
                              Measured=meta(ref.test)$Cmass,
                              ncomp=ncomp_Cmass_log_CVmodel,
                              Species=meta(ref.test)$species,
                              Project=meta(ref.test)$project,
                              functional.group=meta(ref.test)$functional.group,
                              ID=meta(ref.test)$sample_id)

LMA.log.jack.pred<-exp(apply.coefs(LMA.log.jack.coefs,as.matrix(ref.test)))
LMA.log.jack.stat<-t(apply(LMA.log.jack.pred,1,
                           function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
LMA.log.jack.df<-data.frame(pred.mean=LMA.log.jack.stat[,1],
                            pred.low=LMA.log.jack.stat[,2],
                            pred.high=LMA.log.jack.stat[,3],
                            Measured=meta(ref.test)$LMA,
                            ncomp=ncomp_LMA_log_CVmodel,
                            Species=meta(ref.test)$species,
                            Project=meta(ref.test)$project,
                            functional.group=meta(ref.test)$functional.group,
                            ID=meta(ref.test)$sample_id)

K_mass.log.jack.pred<-exp(apply.coefs(K_mass.log.jack.coefs,as.matrix(ref.test)))
K_mass.log.jack.stat<-t(apply(K_mass.log.jack.pred,1,
                              function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
K_mass.log.jack.df<-data.frame(pred.mean=K_mass.log.jack.stat[,1],
                               pred.low=K_mass.log.jack.stat[,2],
                               pred.high=K_mass.log.jack.stat[,3],
                               Measured=meta(ref.test)$K_mass,
                               ncomp=ncomp_K_mass_log_CVmodel,
                               Species=meta(ref.test)$species,
                               Project=meta(ref.test)$project,
                               functional.group=meta(ref.test)$functional.group,
                               ID=meta(ref.test)$sample_id)

all.sqrt.jack.coef.list<-list(chlA=chlA_mass.sqrt.jack.coefs,
                              N=Nmass.sqrt.jack.coefs,
                              C=Cmass.sqrt.jack.coefs,
                              LMA=LMA.sqrt.jack.coefs,
                              K=K_mass.sqrt.jack.coefs)

all.log.jack.coef.list<-list(chlA=chlA_mass.log.jack.coefs,
                             N=Nmass.log.jack.coefs,
                             C=Cmass.log.jack.coefs,
                             LMA=LMA.log.jack.coefs,
                             K=K_mass.log.jack.coefs)

saveRDS(all.sqrt.jack.coef.list,"SavedResults/all_jack_coefs_list_ref_sqrt.rds")
saveRDS(all.log.jack.coef.list,"SavedResults/all_jack_coefs_list_ref_log.rds")

all.sqrt.jack.df.list<-list(chlA=chlA_mass.sqrt.jack.df,
                            N=Nmass.sqrt.jack.df,
                            C=Cmass.sqrt.jack.df,
                            LMA=LMA.sqrt.jack.df,
                            K=K_mass.sqrt.jack.df)

all.log.jack.df.list<-list(chlA=chlA_mass.log.jack.df,
                           N=Nmass.log.jack.df,
                           C=Cmass.log.jack.df,
                           LMA=LMA.log.jack.df,
                           K=K_mass.log.jack.df)

saveRDS(all.sqrt.jack.df.list,"SavedResults/all_jack_df_list_ref_sqrt.rds")
saveRDS(all.log.jack.df.list,"SavedResults/all_jack_df_list_ref_log.rds")

all.sqrt.jack.stats.list<-list(chlA=chlA_mass.sqrt.jack.stats,
                               N=Nmass.sqrt.jack.stats,
                               C=Cmass.sqrt.jack.stats,
                               LMA=LMA.sqrt.jack.stats,
                               K=K_mass.sqrt.jack.stats)

all.log.jack.stats.list<-list(chlA=chlA_mass.log.jack.stats,
                              N=Nmass.log.jack.stats,
                              C=Cmass.log.jack.stats,
                              LMA=LMA.log.jack.stats,
                              K=K_mass.log.jack.stats)

saveRDS(all.sqrt.jack.stats.list,"SavedResults/all_jack_stats_list_ref_sqrt.rds")
saveRDS(all.log.jack.stats.list,"SavedResults/all_jack_stats_list_ref_log.rds")