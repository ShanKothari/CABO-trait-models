setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(pls)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

source("Scripts/CABO-trait-models/00 useful_functions.R")

####################################
## the purpose of this script is to add on the corrected values of EWT
## without creating a new training/testing data split
## script 03 has been modified so that the trait data
## incorporates the corrected values of EWT in the file
## all_ref_and_traits.rds
## we'll use those values to create a new column in the
## training and testing splits

all.ref<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")

ref.train<-readRDS("ProcessedSpectra/ref_train.rds")
ref.test<-readRDS("ProcessedSpectra/ref_test.rds")

## attach correct EWT values to training and testing data
meta(ref.train)$EWT_actual<-meta(all.ref)$EWT[match(meta(ref.train)$sample_id,
                                                    meta(all.ref)$sample_id)]
meta(ref.test)$EWT_actual<-meta(all.ref)$EWT[match(meta(ref.test)$sample_id,
                                                   meta(all.ref)$sample_id)]

## compare old and corrected EWT values
all.ref.old<-combine(ref.train,ref.test)
summary(lm(EWT_actual~EWT,data=meta(all.ref.old)))
with(meta(ref.train),quantile(EWT_actual/EWT,
                              probs=c(0.025,0.25,0.5,0.75,0.975),
                              na.rm=T))

###############################################
## rerun model development and internal validation
## reflectance first

## note that many variable names are shared between
## the sections for reflectance, transmittance, and
## absorptance-based models, so it is important to
## clear the environment (e.g. rm(list=ls()) )
## between sections

EWT_actual_CVmodel<-plsr(meta(ref.train)$EWT_actual~as.matrix(ref.train),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_actual_CVmodel <- selectNcomp(EWT_actual_CVmodel, method = "onesigma", plot = FALSE)
EWT_actual_valid <- which(!is.na(meta(ref.train)$EWT_actual))
EWT_actual_pred<-data.frame(ID=meta(ref.train)$sample_id[EWT_actual_valid],
                     Species=meta(ref.train)$species[EWT_actual_valid],
                     Project=meta(ref.train)$project[EWT_actual_valid],
                     measured=meta(ref.train)$EWT_actual[EWT_actual_valid],
                     val_pred=EWT_actual_CVmodel$validation$pred[,,ncomp_EWT_actual_CVmodel])
ggplot(EWT_actual_pred,aes(y=measured,x=val_pred,color=Project))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  ##  coord_cartesian(xlim=c(38,55),ylim=c(38,55))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting %EWT from fresh-leaf spectra")

EWT_actual.jack.coefs<-list()
EWT_actual.jack.stats<-list()
nreps<-100

for(i in 1:nreps){
  print(i)
  
  n.cal.spec<-nrow(ref.train)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-ref.train[train.jack]
  val.jack<-ref.train[test.jack]

  EWT_actual.jack<-plsr(meta(calib.jack)$EWT_actual~as.matrix(calib.jack),
                 ncomp=30,method = "oscorespls",validation="none")

  EWT_actual.jack.val.pred<-as.vector(predict(EWT_actual.jack,newdata=as.matrix(val.jack),ncomp=ncomp_EWT_actual_CVmodel)[,,1])
  EWT_actual.jack.val.fit<-lm(EWT_actual.jack.val.pred~meta(val.jack)$EWT_actual)
  EWT_actual.jack.stats[[i]]<-c(R2=summary(EWT_actual.jack.val.fit)$r.squared,
                         RMSE=RMSD(meta(val.jack)$EWT_actual,EWT_actual.jack.val.pred),
                         perRMSE=percentRMSD(meta(val.jack)$EWT_actual,EWT_actual.jack.val.pred,0.025,0.975),
                         bias=mean(EWT_actual.jack.val.pred,na.rm=T)-mean(meta(val.jack)$EWT_actual,na.rm=T))
  
  EWT_actual.jack.coefs[[i]]<-as.vector(coef(EWT_actual.jack,ncomp=ncomp_EWT_actual_CVmodel,intercept=TRUE))

}

EWT_actual.jack.pred<-apply.coefs(EWT_actual.jack.coefs,as.matrix(ref.test))
EWT_actual.jack.stat<-t(apply(EWT_actual.jack.pred,1,
                       function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_actual.jack.df<-data.frame(pred.mean=EWT_actual.jack.stat[,1],
                        pred.low=EWT_actual.jack.stat[,2],
                        pred.high=EWT_actual.jack.stat[,3],
                        Measured=meta(ref.test)$EWT_actual,
                        ncomp=ncomp_EWT_actual_CVmodel,
                        Species=meta(ref.test)$species,
                        Project=meta(ref.test)$project,
                        functional.group=meta(ref.test)$functional.group,
                        ID=meta(ref.test)$sample_id)

saveRDS(EWT_actual.jack.coefs,"SavedResults/EWT_corrected_jack_coefs_ref.rds")
saveRDS(EWT_actual.jack.df,"SavedResults/EWT_corrected_jack_df_ref.rds")
saveRDS(EWT_actual.jack.stats,"SavedResults/EWT_corrected_jack_stats_ref.rds")

rm(list=ls())

############################################
## rerun model development and internal validation
## transmittance

all.trans<-readRDS("ProcessedSpectra/all_trans_and_traits.rds")
trans.train<-readRDS("ProcessedSpectra/trans_train.rds")
trans.test<-readRDS("ProcessedSpectra/trans_test.rds")

meta(trans.train)$EWT_actual<-meta(all.trans)$EWT[match(meta(trans.train)$sample_id,
                                                        meta(all.trans)$sample_id)]
meta(trans.test)$EWT_actual<-meta(all.trans)$EWT[match(meta(trans.test)$sample_id,
                                                       meta(all.trans)$sample_id)]

###############################################
## rerun model development and internal validation
## absorptance

all.abs<-readRDS("ProcessedSpectra/all_abs_and_traits.rds")
abs.train<-readRDS("ProcessedSpectra/abs_train.rds")
abs.test<-readRDS("ProcessedSpectra/abs_test.rds")

meta(abs.train)$EWT_actual<-meta(all.abs)$EWT[match(meta(abs.train)$sample_id,
                                                        meta(all.abs)$sample_id)]
meta(abs.test)$EWT_actual<-meta(all.abs)$EWT[match(meta(abs.test)$sample_id,
                                                       meta(all.abs)$sample_id)]


###############################################
## external validation