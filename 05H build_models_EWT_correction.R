setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(pls)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

####################################
## the purpose of this script is to add on the corrected values of EWT
## without creating a new training/testing data split
## script 03 has been modified so that the trait data
## incorporates the corrected values of EWT in the file
## all_ref_and_traits.rds
## we'll use those values to create a new column in the
## training and testing splits

source("Scripts/CABO-trait-models/00 useful_functions.R")
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
with(meta(all.ref.old),quantile(EWT_actual/EWT,
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

source("Scripts/CABO-trait-models/00 useful_functions.R")
all.trans<-readRDS("ProcessedSpectra/all_trans_and_traits.rds")
trans.train<-readRDS("ProcessedSpectra/trans_train.rds")
trans.test<-readRDS("ProcessedSpectra/trans_test.rds")

meta(trans.train)$EWT_actual<-meta(all.trans)$EWT[match(meta(trans.train)$sample_id,
                                                        meta(all.trans)$sample_id)]
meta(trans.test)$EWT_actual<-meta(all.trans)$EWT[match(meta(trans.test)$sample_id,
                                                       meta(all.trans)$sample_id)]

EWT_actual_CVmodel<-plsr(meta(trans.train)$EWT_actual~as.matrix(trans.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_actual_CVmodel <- selectNcomp(EWT_actual_CVmodel, method = "onesigma", plot = FALSE)
EWT_actual_valid <- which(!is.na(meta(trans.train)$EWT_actual))
EWT_actual_pred<-data.frame(ID=meta(trans.train)$sample_id[EWT_actual_valid],
                            Species=meta(trans.train)$species[EWT_actual_valid],
                            Project=meta(trans.train)$project[EWT_actual_valid],
                            measured=meta(trans.train)$EWT_actual[EWT_actual_valid],
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
  
  n.cal.spec<-nrow(trans.train)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-trans.train[train.jack]
  val.jack<-trans.train[test.jack]
  
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

EWT_actual.jack.pred<-apply.coefs(EWT_actual.jack.coefs,as.matrix(trans.test))
EWT_actual.jack.stat<-t(apply(EWT_actual.jack.pred,1,
                              function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_actual.jack.df<-data.frame(pred.mean=EWT_actual.jack.stat[,1],
                               pred.low=EWT_actual.jack.stat[,2],
                               pred.high=EWT_actual.jack.stat[,3],
                               Measured=meta(trans.test)$EWT_actual,
                               ncomp=ncomp_EWT_actual_CVmodel,
                               Species=meta(trans.test)$species,
                               Project=meta(trans.test)$project,
                               functional.group=meta(trans.test)$functional.group,
                               ID=meta(trans.test)$sample_id)

saveRDS(EWT_actual.jack.coefs,"SavedResults/EWT_corrected_jack_coefs_trans.rds")
saveRDS(EWT_actual.jack.df,"SavedResults/EWT_corrected_jack_df_trans.rds")
saveRDS(EWT_actual.jack.stats,"SavedResults/EWT_corrected_jack_stats_trans.rds")

rm(list=ls())

###############################################
## rerun model development and internal validation
## absorptance

source("Scripts/CABO-trait-models/00 useful_functions.R")
all.abs<-readRDS("ProcessedSpectra/all_abs_and_traits.rds")
abs.train<-readRDS("ProcessedSpectra/abs_train.rds")
abs.test<-readRDS("ProcessedSpectra/abs_test.rds")

meta(abs.train)$EWT_actual<-meta(all.abs)$EWT[match(meta(abs.train)$sample_id,
                                                        meta(all.abs)$sample_id)]
meta(abs.test)$EWT_actual<-meta(all.abs)$EWT[match(meta(abs.test)$sample_id,
                                                       meta(all.abs)$sample_id)]

EWT_actual_CVmodel<-plsr(meta(abs.train)$EWT_actual~as.matrix(abs.train),
                         ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_EWT_actual_CVmodel <- selectNcomp(EWT_actual_CVmodel, method = "onesigma", plot = FALSE)
EWT_actual_valid <- which(!is.na(meta(abs.train)$EWT_actual))
EWT_actual_pred<-data.frame(ID=meta(abs.train)$sample_id[EWT_actual_valid],
                            Species=meta(abs.train)$species[EWT_actual_valid],
                            Project=meta(abs.train)$project[EWT_actual_valid],
                            measured=meta(abs.train)$EWT_actual[EWT_actual_valid],
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
  
  n.cal.spec<-nrow(abs.train)
  train.jack<-sample(1:n.cal.spec,floor(0.7*n.cal.spec))
  test.jack<-setdiff(1:n.cal.spec,train.jack)
  
  calib.jack<-abs.train[train.jack]
  val.jack<-abs.train[test.jack]
  
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

EWT_actual.jack.pred<-apply.coefs(EWT_actual.jack.coefs,as.matrix(abs.test))
EWT_actual.jack.stat<-t(apply(EWT_actual.jack.pred,1,
                              function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_actual.jack.df<-data.frame(pred.mean=EWT_actual.jack.stat[,1],
                               pred.low=EWT_actual.jack.stat[,2],
                               pred.high=EWT_actual.jack.stat[,3],
                               Measured=meta(abs.test)$EWT_actual,
                               ncomp=ncomp_EWT_actual_CVmodel,
                               Species=meta(abs.test)$species,
                               Project=meta(abs.test)$project,
                               functional.group=meta(abs.test)$functional.group,
                               ID=meta(abs.test)$sample_id)

saveRDS(EWT_actual.jack.coefs,"SavedResults/EWT_corrected_jack_coefs_abs.rds")
saveRDS(EWT_actual.jack.df,"SavedResults/EWT_corrected_jack_df_abs.rds")
saveRDS(EWT_actual.jack.stats,"SavedResults/EWT_corrected_jack_stats_abs.rds")

rm(list=ls())

###############################################
## plotting internal validation output

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

EWT_actual.jack.df.ref<-readRDS("SavedResults/EWT_corrected_jack_df_ref.rds")
EWT_actual.jack.df.trans<-readRDS("SavedResults/EWT_corrected_jack_df_trans.rds")
EWT_actual.jack.df.abs<-readRDS("SavedResults/EWT_corrected_jack_df_abs.rds")

all.EWT<-c(EWT_actual.jack.df.ref$Measured,
           EWT_actual.jack.df.ref$pred.mean,
           EWT_actual.jack.df.trans$pred.mean,
           EWT_actual.jack.df.abs$pred.mean)
EWT_upper<-max(all.EWT,na.rm=T)+0.02
EWT_lower<-min(all.EWT,na.rm=T)-0.02

EWT_actual.ref.val.plot<-ggplot(EWT_actual.jack.df.ref,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

EWT_actual.trans.val.plot<-ggplot(EWT_actual.jack.df.trans,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

EWT_actual.abs.val.plot<-ggplot(EWT_actual.jack.df.abs,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)",
       color="Functional group")+
  scale_color_manual(values=colorBlind)

pdf("Images/EWT_corrected_val_plot.pdf",width=16,height=6)
  (EWT_actual.ref.val.plot + EWT_actual.trans.val.plot + EWT_actual.abs.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

###############################################
## external validation

source("Scripts/CABO-trait-models/00 useful_functions.R")
EWT_actual.jack.coefs.ref<-readRDS("SavedResults/EWT_corrected_jack_coefs_ref.rds")

## read in processed data
LOPEX<-readRDS("IndependentValidationData/LOPEX/LOPEX_processed.rds")
ANGERS<-readRDS("IndependentValidationData/ANGERS/ANGERS_processed.rds")
Dessain.spec<-readRDS("IndependentValidationData/Dessain_processed.rds")

LOPEX_ANGERS<-spectrolab::combine(LOPEX,ANGERS)
all_val_ref<-spectrolab::combine(LOPEX_ANGERS,Dessain.spec)

EWT_pred<-apply.coefs(EWT_actual.jack.coefs.ref,val.spec = all_val_ref)
EWT_stat<-t(apply(EWT_pred,1,
                  function(obs) c(mean(obs),quantile(obs,probs=c(0.025,0.975)))))
EWT_pred_df<-data.frame(Measured=meta(all_val_ref)$EWT,
                        pred.mean=EWT_stat[,1],
                        pred.low=EWT_stat[,2],
                        pred.high=EWT_stat[,3],
                        dataset=meta(all_val_ref)$dataset)
EWT_all<-with(EWT_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
EWT_upper<-max(EWT_all,na.rm=T)+0.03
EWT_lower<-min(EWT_all,na.rm=T)-0.03

EWT_ind_val<-ggplot(data=EWT_pred_df,
                    aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  scale_color_manual(values=colorBlind)
