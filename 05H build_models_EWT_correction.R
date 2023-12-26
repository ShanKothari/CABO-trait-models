setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(pls)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

####################################
## the purpose of this section is to add on the corrected values of EWT
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
## save coefficients

EWT_actual.jack.coefs.ref<-readRDS("SavedResults/EWT_corrected_jack_coefs_ref.rds")
EWT_actual.jack.coefs.trans<-readRDS("SavedResults/EWT_corrected_jack_coefs_trans.rds")
EWT_actual.jack.coefs.abs<-readRDS("SavedResults/EWT_corrected_jack_coefs_abs.rds")

## save coefficients
write.coefs<-function(obj,path,filename){
  coef_mat<-matrix(unlist(obj),nrow=length(obj),byrow=T)
  colnames(coef_mat)<-c("intercept",400:2400)
  write.csv(coef_mat,
            paste(path,filename,".csv",sep=""),
            row.names=F)
}

write.coefs(obj=EWT_actual.jack.coefs.ref,
            path="ModelCoefficients/EWTCorrectedModels/",
            filename="EWTFresh_ref")

write.coefs(obj=EWT_actual.jack.coefs.trans,
            path="ModelCoefficients/EWTCorrectedModels/",
            filename="EWTFresh_trans")

write.coefs(obj=EWT_actual.jack.coefs.abs,
            path="ModelCoefficients/EWTCorrectedModels/",
            filename="EWTFresh_abs")

###############################################
## plotting internal validation output

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

EWT_actual.jack.df.ref<-readRDS("SavedResults/EWT_corrected_jack_df_ref.rds")
EWT_actual.jack.df.trans<-readRDS("SavedResults/EWT_corrected_jack_df_trans.rds")
EWT_actual.jack.df.abs<-readRDS("SavedResults/EWT_corrected_jack_df_abs.rds")

source("Scripts/CABO-trait-models/06A plotting_models_ref.R")
source("Scripts/CABO-trait-models/06B plotting_models_combined.R")

all.EWT<-c(EWT_actual.jack.df.ref$Measured,
           EWT_actual.jack.df.ref$pred.mean,
           EWT_actual.jack.df.trans$pred.mean,
           EWT_actual.jack.df.abs$pred.mean)
EWT_upper<-max(all.EWT,na.rm=T)+0.02
EWT_lower<-min(all.EWT,na.rm=T)-0.02

EWT_actual.ref.val.plot1<-ggplot(EWT_actual.jack.df.ref,
                                 aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

pdf("Images/val_plots_ref1_corrected.pdf",width = 16,height = 19)
(LMA.val.plot+LDMC.val.plot+EWT_actual.ref.val.plot1)/
  (Nmass.val.plot+Cmass.val.plot+solubles_mass.val.plot)/
  (hemicellulose_mass.val.plot+cellulose_mass.val.plot+lignin_mass.val.plot)/
  (chlA_mass.val.plot+chlB_mass.val.plot+car_mass.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

## larger text
EWT_actual.ref.val.plot2<-ggplot(EWT_actual.jack.df.ref,
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

pdf("Images/val_plots_comp1_corrected.pdf",width=16,height=17)
(LMA.ref.val.plot + LMA.trans.val.plot + LMA.abs.val.plot) /
  (LDMC.ref.val.plot + LDMC.trans.val.plot + LDMC.abs.val.plot) /
  (EWT_actual.ref.val.plot2 + EWT_actual.trans.val.plot + EWT_actual.abs.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

summary(lm(Measured~pred.mean,data=EWT_actual.jack.df.abs))
with(EWT_actual.jack.df.abs,
     RMSD(measured = Measured,predicted = pred.mean))
with(EWT_actual.jack.df.abs,
     percentRMSD(measured = Measured,predicted = pred.mean,min=0.025,max=0.975))

###############################################
## external validation

source("Scripts/CABO-trait-models/00 useful_functions.R")
EWT_actual.jack.coefs.ref<-readRDS("SavedResults/EWT_corrected_jack_coefs_ref.rds")
ind_val_list<-readRDS("SavedResults/ind_val_list.rds")

## read in processed data
LOPEX<-readRDS("IndependentValidationData/LOPEX/LOPEX_processed.rds")
ANGERS<-readRDS("IndependentValidationData/ANGERS/ANGERS_processed.rds")
Dessain.spec<-readRDS("IndependentValidationData/Dessain_processed.rds")

LOPEX_ANGERS<-spectrolab::combine(LOPEX,ANGERS)
all_val_ref<-spectrolab::combine(LOPEX_ANGERS,Dessain.spec)

solubles_pred_df<-ind_val_list$solubles_mass
solubles_all<-with(solubles_pred_df,c(pred.low[!is.na(Measured)],
                                      pred.high[!is.na(Measured)],
                                      Measured))
solubles_upper<-max(solubles_all,na.rm=T)+1
solubles_lower<-min(solubles_all,na.rm=T)-1

hemicellulose_pred_df<-ind_val_list$hemicellulose_mass
hemicellulose_all<-with(hemicellulose_pred_df,c(pred.low[!is.na(Measured)],
                                                pred.high[!is.na(Measured)],
                                                Measured))
hemicellulose_upper<-max(hemicellulose_all,na.rm=T)+1
hemicellulose_lower<-min(hemicellulose_all,na.rm=T)-1

cellulose_pred_df<-ind_val_list$cellulose_mass
cellulose_all<-with(cellulose_pred_df,c(pred.low[!is.na(Measured)],
                                        pred.high[!is.na(Measured)],
                                        Measured))
cellulose_upper<-max(cellulose_all,na.rm=T)+1
cellulose_lower<-min(cellulose_all,na.rm=T)-1

lignin_pred_df<-ind_val_list$lignin_mass
lignin_all<-with(lignin_pred_df,c(pred.low[!is.na(Measured)],
                                  pred.high[!is.na(Measured)],
                                  Measured))
lignin_upper<-max(lignin_all,na.rm=T)+1
lignin_lower<-min(lignin_all,na.rm=T)-1

Nmass_pred_df<-ind_val_list$Nmass
Nmass_all<-with(Nmass_pred_df,c(pred.low[!is.na(Measured)],
                                pred.high[!is.na(Measured)],
                                Measured))
Nmass_upper<-max(Nmass_all,na.rm=T)+0.1
Nmass_lower<-min(Nmass_all,na.rm=T)-0.1

Cmass_pred_df<-ind_val_list$Cmass
Cmass_all<-with(Cmass_pred_df,c(pred.low[!is.na(Measured)],
                                pred.high[!is.na(Measured)],
                                Measured))
Cmass_upper<-max(Cmass_all,na.rm=T)+1
Cmass_lower<-min(Cmass_all,na.rm=T)-1

LMA_pred_df<-ind_val_list$LMA
LMA_all<-with(LMA_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
LMA_upper<-max(LMA_all,na.rm=T)+0.02
LMA_lower<-min(LMA_all,na.rm=T)-0.02

LDMC_pred_df<-ind_val_list$LDMC
LDMC_all<-with(LDMC_pred_df,c(pred.low[!is.na(Measured)],
                              pred.high[!is.na(Measured)],
                              Measured))
LDMC_upper<-max(LDMC_all,na.rm=T)+10
LDMC_lower<-min(LDMC_all,na.rm=T)-10

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

chlA_pred_df<-ind_val_list$chlA
chlA_all<-with(chlA_pred_df,c(pred.low[!is.na(Measured)],
                              pred.high[!is.na(Measured)],
                              Measured))
chlA_upper<-max(chlA_all,na.rm=T)+0.5
chlA_lower<-min(chlA_all,na.rm=T)-0.5

chlB_pred_df<-ind_val_list$chlB
chlB_all<-with(chlB_pred_df,c(pred.low[!is.na(Measured)],
                              pred.high[!is.na(Measured)],
                              Measured))
chlB_upper<-max(chlB_all,na.rm=T)+0.3
chlB_lower<-min(chlB_all,na.rm=T)-0.3

car_pred_df<-ind_val_list$car
car_all<-with(car_pred_df,c(pred.low[!is.na(Measured)],
                            pred.high[!is.na(Measured)],
                            Measured))
car_upper<-max(car_all,na.rm=T)+0.2
car_lower<-min(car_all,na.rm=T)-0.2

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

solubles_ind_val<-ggplot(data=solubles_pred_df,
                         aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),
                  ylim=c(solubles_lower,solubles_upper))+
  scale_color_manual(values=colorBlind)

hemicellulose_ind_val<-ggplot(data=hemicellulose_pred_df,
                              aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),
                  ylim=c(hemicellulose_lower,hemicellulose_upper))+
  scale_color_manual(values=colorBlind)

cellulose_ind_val<-ggplot(data=cellulose_pred_df,
                          aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),
                  ylim=c(cellulose_lower,cellulose_upper))+
  scale_color_manual(values=colorBlind)

lignin_ind_val<-ggplot(data=lignin_pred_df,
                       aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),
                  ylim=c(lignin_lower,lignin_upper))+
  scale_color_manual(values=colorBlind)

Nmass_ind_val<-ggplot(data=Nmass_pred_df,
                      aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  coord_cartesian(xlim=c(Nmass_lower,Nmass_upper),
                  ylim=c(Nmass_lower,Nmass_upper))+
  scale_color_manual(values=colorBlind)

Cmass_ind_val<-ggplot(data=Cmass_pred_df,
                      aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  coord_cartesian(xlim=c(Cmass_lower,Cmass_upper),
                  ylim=c(Cmass_lower,Cmass_upper))+
  scale_color_manual(values=colorBlind)

LMA_ind_val<-ggplot(data=LMA_pred_df,
                    aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),
                  ylim=c(LMA_lower,LMA_upper))+
  scale_color_manual(values=colorBlind)

LDMC_ind_val<-ggplot(data=LDMC_pred_df,
                     aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),
                  ylim=c(LDMC_lower,LDMC_upper))+
  scale_color_manual(values=colorBlind)

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
  scale_color_manual(values=colorBlind)+
  guides(color=guide_legend(title="Dataset"))

chlA_ind_val<-ggplot(data=chlA_pred_df,
                     aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),
                  ylim=c(chlA_lower,chlA_upper))+
  scale_color_manual(values=colorBlind)

chlB_ind_val<-ggplot(data=chlB_pred_df,
                     aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),
                  ylim=c(chlB_lower,chlB_upper))+
  scale_color_manual(values=colorBlind)

car_ind_val<-ggplot(data=car_pred_df,
                    aes(x=pred.mean,y=Measured,color=dataset))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm",se=F)+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"))+
  coord_cartesian(xlim=c(car_lower,car_upper),
                  ylim=c(car_lower,car_upper))+
  scale_color_manual(values=colorBlind)

pdf("Images/ind_val_plots1_corrected.pdf",width = 16,height = 19, onefile=F)
ggarrange(plotlist=list(solubles_ind_val,hemicellulose_ind_val,
                        cellulose_ind_val,lignin_ind_val,
                        Nmass_ind_val,Cmass_ind_val,
                        LMA_ind_val,LDMC_ind_val,EWT_ind_val,
                        chlA_ind_val,chlB_ind_val,car_ind_val),
          common.legend = T,legend = "bottom",
          nrow=4,ncol=3)
dev.off()

summary(lm(Measured~pred.mean,data=EWT_pred_df[EWT_pred_df$dataset=="ANGERS",]))
with(EWT_pred_df[EWT_pred_df$dataset=="ANGERS",],
     RMSD(measured = Measured,predicted = pred.mean))
with(EWT_pred_df[EWT_pred_df$dataset=="ANGERS",],
     percentRMSD(measured = Measured,predicted = pred.mean,min=0.025,max=0.975))

##########################################
## plot comparisons between fresh and rehydrated EWT

library(reshape2)

all.ref<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")
meta(all.ref)$EWT_rehydrated<-with(meta(all.ref),(1/(LDMC/1000)-1)*LMA)
EWT_df<-meta(all.ref)[,c("sample_id","EWT","EWT_rehydrated")]
EWT_df_long<-melt(EWT_df, id.vars="sample_id")
levels(EWT_df_long$variable)<-c("Fresh EWT","Rehydrated EWT")

EWT_corr_density<-ggplot(data=EWT_df_long,
                         aes(x=value,color=variable))+
  geom_density(size=1.5,bounds=c(0,Inf))+
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=30),
        legend.position = c(0.7,0.7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  labs(y="Density",x="EWT (mm)")

EWT_corr_plot<-ggplot(data=meta(all.ref),
                      aes(x=EWT_rehydrated,y=EWT))+
  geom_point(size=2)+
  geom_abline(intercept=0,slope=1,size=2,linetype="dashed")+
  theme_bw()+
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=c(0,0.7),ylim=c(0,0.7))+
  labs(y="Rehydrated EWT (mm)",x="Fresh EWT (mm)")
