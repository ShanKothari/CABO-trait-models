setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(caret)
library(pls)
library(lmodel2)
library(reshape2)
library(hypervolume)
library(GSAR)
library(statip)
library(lme4)
library(FNN)

spec.traits<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")
## the Pardo dataset has no trait data (yet)
spec.traits<-spec.traits[-which(meta(spec.traits)$project=="2019-Pardo-MSc-UdeM")]

#########################################
## define RMSD function

RMSD<-function(measured,predicted){
  not.na<-which(!is.na(measured) & !is.na(predicted))
  return(sqrt(sum((measured-predicted)^2,na.rm=T)/(length(not.na)-1)))
}

percentRMSD<-function(measured,predicted,min,max){
  RMSD_data<-RMSD(measured,predicted)
  range<-quantile(measured,probs=max)-quantile(measured,probs=min)
  return(RMSD_data/range)
}

#########################################
## processing spectral data

## spectral dimensionality reduction
spec.traits.ref=as.data.frame(spec.traits,metadata=F)[,as.character(400:2400)]
colnames(spec.traits.ref)<-paste("X",colnames(spec.traits.ref),sep="")

spec.pca<-prcomp(~.,spec.traits.ref,scale.=F)
spec.pca.scale<-prcomp(~.,spec.traits.ref,scale.=T)

meta(spec.traits)$PC1.spec<-spec.pca$x[,1]
meta(spec.traits)$PC2.spec<-spec.pca$x[,2]
meta(spec.traits)$PC3.spec<-spec.pca$x[,3]
meta(spec.traits)$PC4.spec<-spec.pca$x[,4]
meta(spec.traits)$PC5.spec<-spec.pca$x[,5]

meta(spec.traits)$PC1.spec.scale<-spec.pca.scale$x[,1]
meta(spec.traits)$PC2.spec.scale<-spec.pca.scale$x[,2]
meta(spec.traits)$PC3.spec.scale<-spec.pca.scale$x[,3]
meta(spec.traits)$PC4.spec.scale<-spec.pca.scale$x[,4]
meta(spec.traits)$PC5.spec.scale<-spec.pca.scale$x[,5]

meta(spec.traits)$PC1.spec.scale.norm<-meta(spec.traits)$PC1.spec.scale/sd(meta(spec.traits)$PC1.spec.scale)
meta(spec.traits)$PC2.spec.scale.norm<-meta(spec.traits)$PC2.spec.scale/sd(meta(spec.traits)$PC2.spec.scale)
meta(spec.traits)$PC3.spec.scale.norm<-meta(spec.traits)$PC3.spec.scale/sd(meta(spec.traits)$PC3.spec.scale)
meta(spec.traits)$PC4.spec.scale.norm<-meta(spec.traits)$PC4.spec.scale/sd(meta(spec.traits)$PC4.spec.scale)
meta(spec.traits)$PC5.spec.scale.norm<-meta(spec.traits)$PC5.spec.scale/sd(meta(spec.traits)$PC5.spec.scale)

## trait dimensionality reduction

traits.ref=meta(spec.traits)[,c("LDMC","EWT","Cmass","Nmass","LMA",
                                "solubles_mass","hemicellulose_mass",
                                "cellulose_mass","lignin_mass",
                                "chlA_mass","chlB_mass","car_mass")]

traits.pca.scale<-prcomp(~.,traits.ref,scale.=T,na.action=na.exclude)
meta(spec.traits)$PC1.traits.scale<-traits.pca.scale$x[,1]
meta(spec.traits)$PC2.traits.scale<-traits.pca.scale$x[,2]
meta(spec.traits)$PC3.traits.scale<-traits.pca.scale$x[,3]
meta(spec.traits)$PC4.traits.scale<-traits.pca.scale$x[,4]
meta(spec.traits)$PC5.traits.scale<-traits.pca.scale$x[,5]

meta(spec.traits)$PC1.traits.scale.norm<-meta(spec.traits)$PC1.traits.scale/sd(meta(spec.traits)$PC1.traits.scale)
meta(spec.traits)$PC2.traits.scale.norm<-meta(spec.traits)$PC2.traits.scale/sd(meta(spec.traits)$PC2.traits.scale)
meta(spec.traits)$PC3.traits.scale.norm<-meta(spec.traits)$PC3.traits.scale/sd(meta(spec.traits)$PC3.traits.scale)
meta(spec.traits)$PC4.traits.scale.norm<-meta(spec.traits)$PC4.traits.scale/sd(meta(spec.traits)$PC4.traits.scale)
meta(spec.traits)$PC5.traits.scale.norm<-meta(spec.traits)$PC5.traits.scale/sd(meta(spec.traits)$PC5.traits.scale)

## find projects with more than 30 samples (all but CABO-General)
common.proj<-names(table(meta(spec.traits)$project))[table(meta(spec.traits)$project)>30]

########################################
## spectral overlap analyses

spec_overlap<-list()
spec_scale_overlap<-list()
spec_scale_norm_overlap<-list()
spec_pca_cols<-c("PC1.spec","PC2.spec","PC3.spec","PC4.spec","PC5.spec")
spec_pca_scale_cols<-c("PC1.spec.scale","PC2.spec.scale",
                       "PC3.spec.scale","PC4.spec.scale",
                       "PC5.spec.scale")
spec_pca_scale_norm_cols<-c("PC1.spec.scale.norm","PC2.spec.scale.norm",
                       "PC3.spec.scale.norm","PC4.spec.scale.norm",
                       "PC5.spec.scale.norm")

for(i in common.proj){
  print(i)
  spec.traits.proj<-spec.traits[meta(spec.traits)$project==i]
  spec.traits.other<-spec.traits[meta(spec.traits)$project!=i]
  
  # spec_proj_KDE<-hypervolume(meta(spec.traits.proj)[,spec_pca_cols[1:3]])
  # spec_other_KDE<-hypervolume(meta(spec.traits.other)[,spec_pca_cols[1:3]])
  # spec_set<-hypervolume_set(spec_proj_KDE,spec_other_KDE,check.memory = F)
  # 
  # spec_scale_proj_KDE<-hypervolume(meta(spec.traits.proj)[,spec_pca_scale_cols[1:3]])
  # spec_scale_other_KDE<-hypervolume(meta(spec.traits.other)[,spec_pca_scale_cols[1:3]])
  # spec_scale_set<-hypervolume_set(spec_scale_proj_KDE,spec_scale_other_KDE,check.memory = F)

  spec_scale_norm_proj_KDE<-hypervolume(meta(spec.traits.proj)[,spec_pca_scale_norm_cols[1:3]])
  spec_scale_norm_other_KDE<-hypervolume(meta(spec.traits.other)[,spec_pca_scale_norm_cols[1:3]])
  spec_scale_norm_set<-hypervolume_set(spec_scale_norm_proj_KDE,spec_scale_norm_other_KDE,check.memory = F)
  
  # spec_overlap[[i]]<-hypervolume_overlap_statistics(spec_set)
  # spec_scale_overlap[[i]]<-hypervolume_overlap_statistics(spec_scale_set)
  spec_scale_norm_overlap[[i]]<-hypervolume_overlap_statistics(spec_scale_norm_set)
}

spec_overlap_df<-data.frame(dataset=names(spec_scale_norm_overlap),
                            # jaccard=unlist(lapply(spec_overlap,function(x) x[["jaccard"]])),
                            # sorensen=unlist(lapply(spec_overlap,function(x) x[["sorensen"]])),
                            # unique_frac=unlist(lapply(spec_overlap,function(x) x[["frac_unique_1"]])),
                            # jaccard_scale=unlist(lapply(spec_scale_overlap,function(x) x[["jaccard"]])),
                            # sorensen_scale=unlist(lapply(spec_scale_overlap,function(x) x[["sorensen"]])),
                            # unique_frac_scale=unlist(lapply(spec_scale_overlap,function(x) x[["frac_unique_1"]])),
                            jaccard_scale_norm=unlist(lapply(spec_scale_norm_overlap,function(x) x[["jaccard"]])),
                            sorensen_scale_norm=unlist(lapply(spec_scale_norm_overlap,function(x) x[["sorensen"]])),
                            unique_frac_scale_norm=unlist(lapply(spec_scale_norm_overlap,function(x) x[["frac_unique_1"]])))

#######################################
## spectral distance analyses

spec_distance<-list()

for(i in common.proj){
  print(i)
  data_dummy<-vector(mode="numeric",length=nrow(meta(spec.traits)))
  data_dummy[meta(spec.traits)$project==i]<-1
  data_dummy[meta(spec.traits)$project!=i]<-2
  
  spec_pca_select<-t(as.matrix(meta(spec.traits)[,spec_pca_cols]))
  spec_pca_scale_select<-t(as.matrix(meta(spec.traits)[,spec_pca_scale_cols]))
  spec_pca_scale_norm_select<-t(as.matrix(meta(spec.traits)[,spec_pca_scale_norm_cols]))
  
  KS_spec<-KStest(spec_pca_select,group=data_dummy,pvalue.only=F)
  KS_spec_scale<-KStest(spec_pca_scale_select,group=data_dummy,pvalue.only=F)
  KS_spec_scale_norm<-KStest(spec_pca_scale_norm_select,group=data_dummy,pvalue.only=F)
  spec_distance[[i]]<-c(KS_spec$statistic,KS_spec_scale$statistic,KS_spec_scale_norm$statistic)
}

spec_distance_df<-data.frame(dataset=names(spec_distance),
                             KS_stat=unlist(lapply(spec_distance,function(x) x[[1]])),
                             KS_stat_scale=unlist(lapply(spec_distance,function(x) x[[2]])),
                             KS_stat_scale_norm=unlist(lapply(spec_distance,function(x) x[[3]])))

saveRDS(spec_distance,"SavedResults/spec_distance.rds")

#######################################
## trait overlap analyses

trait_overlap<-list()
traits_pca_scale_cols<-c("PC1.traits.scale","PC2.traits.scale",
                         "PC3.traits.scale","PC4.traits.scale",
                         "PC5.traits.scale")
traits_pca_scale_norm_cols<-c("PC1.traits.scale.norm","PC2.traits.scale.norm",
                         "PC3.traits.scale.norm","PC4.traits.scale.norm",
                         "PC5.traits.scale.norm")
spec.traits.na.omit<-spec.traits[!is.na(meta(spec.traits)$PC1.traits.scale)]

for(i in common.proj){
  print(i)
  traits.proj<-meta(spec.traits.na.omit)[meta(spec.traits.na.omit)$project==i,traits_pca_scale_cols[1:3]]
  traits.other<-meta(spec.traits.na.omit)[meta(spec.traits.na.omit)$project!=i,traits_pca_scale_cols[1:3]]
  
  ## hypervolume approach
  trait_proj_KDE<-hypervolume(traits.proj)
  trait_other_KDE<-hypervolume(traits.other)
  trait_set<-hypervolume_set(trait_proj_KDE,trait_other_KDE,check.memory = F)
  trait_overlap[[i]]<-hypervolume_overlap_statistics(trait_set)
}

trait_overlap_df<-data.frame(dataset=names(trait_overlap),
                              jaccard=unlist(lapply(trait_overlap,function(x) x[["jaccard"]])),
                              sorensen=unlist(lapply(trait_overlap,function(x) x[["sorensen"]])),
                              unique_frac=unlist(lapply(trait_overlap,function(x) x[["frac_unique_1"]])))

########################################
## trait distance analyses

traits_distance<-list()

for(i in common.proj){
  print(i)
  data_dummy<-vector(mode="numeric",length=nrow(meta(spec.traits.na.omit)))
  data_dummy[meta(spec.traits.na.omit)$project==i]<-1
  data_dummy[meta(spec.traits.na.omit)$project!=i]<-2
  
  traits_pca_scale_select<-t(as.matrix(meta(spec.traits.na.omit)[,traits_pca_scale_cols]))
  traits_pca_scale_norm_select<-t(as.matrix(meta(spec.traits.na.omit)[,traits_pca_scale_norm_cols]))
  
  KS_traits_scale<-KStest(traits_pca_scale_select,group=data_dummy,pvalue.only=F)
  KS_traits_scale_norm<-KStest(traits_pca_scale_norm_select,group=data_dummy,pvalue.only=F)
  
  traits_distance[[i]]<-c(KS_traits_scale$statistic,KS_traits_scale_norm$statistic)
}

traits_distance_df<-data.frame(dataset=names(traits_distance),
                             KS_stat_scale=unlist(lapply(traits_distance,function(x) x[[1]])),
                             KS_stat_scale_norm=unlist(lapply(traits_distance,function(x) x[[2]])))

saveRDS(traits_distance,"SavedResults/traits_distance.rds")

########################################
## leave one project out' trait prediction analyses

LMA_pred_list<-list()
LDMC_pred_list<-list()
Nmass_pred_list<-list()
lignin_mass_pred_list<-list()
chlA_mass_pred_list<-list()
EWT_pred_list<-list()

LMA_dist<-list()
LDMC_dist<-list()
Nmass_dist<-list()
lignin_mass_dist<-list()
chlA_mass_dist<-list()
EWT_dist<-list()

LMA_KL<-list()
LDMC_KL<-list()
Nmass_KL<-list()
lignin_mass_KL<-list()
chlA_mass_KL<-list()
EWT_KL<-list()

for(i in common.proj){
  print(i)
  spec.traits.proj<-spec.traits[meta(spec.traits)$project==i]
  spec.traits.other<-spec.traits[meta(spec.traits)$project!=i]
  
  LMA_cal<-plsr(meta(spec.traits.other)$LMA~as.matrix(spec.traits.other),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_LMA<-selectNcomp(LMA_cal, method = "onesigma", plot = F)
  LMA_val<-predict(LMA_cal,newdata=as.matrix(spec.traits.proj),ncomp=ncomp_LMA)
  LMA_pred_list[[i]]<-data.frame(measured=meta(spec.traits.proj)$LMA,
                                 val_pred=as.vector(LMA_val))
  LMA_dist[[i]]<-hellinger(na.omit(meta(spec.traits.proj)$LMA),
                           na.omit(meta(spec.traits.other)$LMA),
                           lower=-0.2,
                           upper=0.6)
  LMA_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$LMA),
                             na.omit(meta(spec.traits.proj)$LMA))
  
  LDMC_cal<-plsr(meta(spec.traits.other)$LDMC~as.matrix(spec.traits.other),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_LDMC<-selectNcomp(LDMC_cal, method = "onesigma", plot = F)
  LDMC_val<-predict(LDMC_cal,newdata=as.matrix(spec.traits.proj),ncomp=ncomp_LDMC)
  LDMC_pred_list[[i]]<-data.frame(measured=meta(spec.traits.proj)$LDMC,
                                 val_pred=as.vector(LDMC_val))
  LDMC_dist[[i]]<-hellinger(na.omit(meta(spec.traits.proj)$LDMC),
                           na.omit(meta(spec.traits.other)$LDMC),
                           lower=-200,
                           upper=800)
  LDMC_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$LDMC),
                             na.omit(meta(spec.traits.proj)$LDMC))
  
  Nmass_cal<-plsr(meta(spec.traits.other)$Nmass~as.matrix(spec.traits.other),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_Nmass<-selectNcomp(Nmass_cal, method = "onesigma", plot = F)
  Nmass_val<-predict(Nmass_cal,newdata=as.matrix(spec.traits.proj),ncomp=ncomp_Nmass)
  Nmass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.proj)$Nmass,
                                 val_pred=as.vector(Nmass_val))
  Nmass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.proj)$Nmass),
                           na.omit(meta(spec.traits.other)$Nmass),
                           lower=-1,
                           upper=8)
  Nmass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$Nmass),
                             na.omit(meta(spec.traits.proj)$Nmass))
  
  lignin_mass_cal<-plsr(meta(spec.traits.other)$lignin_mass~as.matrix(spec.traits.other),
                  ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_lignin_mass<-selectNcomp(lignin_mass_cal, method = "onesigma", plot = F)
  lignin_mass_val<-predict(lignin_mass_cal,newdata=as.matrix(spec.traits.proj),ncomp=ncomp_lignin_mass)
  lignin_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.proj)$lignin_mass,
                                   val_pred=as.vector(lignin_mass_val))
  lignin_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.proj)$lignin_mass),
                             na.omit(meta(spec.traits.other)$lignin_mass),
                             lower=-5,
                             upper=30)
  lignin_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$lignin_mass),
                                     na.omit(meta(spec.traits.proj)$lignin_mass))
  
  chlA_mass_cal<-plsr(meta(spec.traits.other)$chlA_mass~as.matrix(spec.traits.other),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_chlA_mass<-selectNcomp(chlA_mass_cal, method = "onesigma", plot = F)
  chlA_mass_val<-predict(chlA_mass_cal,newdata=as.matrix(spec.traits.proj),ncomp=ncomp_chlA_mass)
  chlA_mass_pred_list[[i]]<-data.frame(measured=meta(spec.traits.proj)$chlA_mass,
                                         val_pred=as.vector(chlA_mass_val))
  chlA_mass_dist[[i]]<-hellinger(na.omit(meta(spec.traits.proj)$chlA_mass),
                                   na.omit(meta(spec.traits.other)$chlA_mass),
                                   lower=-5,
                                   upper=30)
  chlA_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$chlA_mass),
                             na.omit(meta(spec.traits.proj)$chlA_mass))
  
  EWT_cal<-plsr(meta(spec.traits.other)$EWT~as.matrix(spec.traits.other),
                        ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_EWT<-selectNcomp(EWT_cal, method = "onesigma", plot = F)
  EWT_val<-predict(EWT_cal,newdata=as.matrix(spec.traits.proj),ncomp=ncomp_EWT)
  EWT_pred_list[[i]]<-data.frame(measured=meta(spec.traits.proj)$EWT,
                                         val_pred=as.vector(EWT_val))
  EWT_dist[[i]]<-hellinger(na.omit(meta(spec.traits.proj)$EWT),
                                   na.omit(meta(spec.traits.other)$EWT),
                                   lower=-0.02,
                                   upper=0.09)
  EWT_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$EWT),
                             na.omit(meta(spec.traits.proj)$EWT))
  
}

LMA_sum_df<-data.frame(dataset=names(LMA_pred_list),
                       RMSD=unlist(lapply(LMA_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                       R2=unlist(lapply(LMA_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                       hellinger=unlist(LMA_dist),
                       KL=unlist(lapply(LMA_KL,function(x) x[[5]])),
                       trait="LMA")

LDMC_sum_df<-data.frame(dataset=names(LDMC_pred_list),
                       RMSD=unlist(lapply(LDMC_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                       R2=unlist(lapply(LDMC_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                       hellinger=unlist(LDMC_dist),
                       KL=unlist(lapply(LDMC_KL,function(x) x[[5]])),
                       trait="LDMC")

Nmass_sum_df<-data.frame(dataset=names(Nmass_pred_list),
                        RMSD=unlist(lapply(Nmass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                        R2=unlist(lapply(Nmass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                        hellinger=unlist(Nmass_dist),
                        KL=unlist(lapply(Nmass_KL,function(x) x[[5]])),
                        trait="Nmass")

lignin_mass_sum_df<-data.frame(dataset=names(lignin_mass_pred_list),
                         RMSD=unlist(lapply(lignin_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                         R2=unlist(lapply(lignin_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                         hellinger=unlist(lignin_mass_dist),
                         KL=unlist(lapply(lignin_mass_KL,function(x) x[[5]])),
                         trait="lignin_mass")

chlA_mass_sum_df<-data.frame(dataset=names(chlA_mass_pred_list),
                               RMSD=unlist(lapply(chlA_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                               R2=unlist(lapply(chlA_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                               hellinger=unlist(chlA_mass_dist),
                             KL=unlist(lapply(chlA_mass_KL,function(x) x[[5]])),
                               trait="chlA_mass")

EWT_sum_df<-data.frame(dataset=names(EWT_pred_list),
                       RMSD=unlist(lapply(EWT_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                       R2=unlist(lapply(EWT_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                       hellinger=unlist(EWT_dist),
                       KL=unlist(lapply(EWT_KL,function(x) x[[5]])),
                       trait="EWT")

summary_df<-do.call(rbind,list(LMA_sum_df,LDMC_sum_df,Nmass_sum_df,
                               lignin_mass_sum_df,chlA_mass_sum_df,
                               EWT_sum_df))

summary_df$KS_stat<-spec_distance_df$KS_stat[match(summary_df$dataset,spec_distance_df$dataset)]
summary_df$KS_stat_scale<-spec_distance_df$KS_stat_scale[match(summary_df$dataset,spec_distance_df$dataset)]
summary_df$KS_stat_scale_norm<-spec_distance_df$KS_stat_scale_norm[match(summary_df$dataset,spec_distance_df$dataset)]

# summary_df$jaccard<-spec_overlap_df$jaccard[match(summary_df$dataset,spec_overlap_df$dataset)]
# summary_df$unique_frac<-spec_overlap_df$unique_frac[match(summary_df$dataset,spec_overlap_df$dataset)]
summary_df$jaccard_scale_norm<-spec_overlap_df$jaccard_scale_norm[match(summary_df$dataset,spec_overlap_df$dataset)]
summary_df$unique_frac_scale_norm<-spec_overlap_df$unique_frac_scale_norm[match(summary_df$dataset,spec_overlap_df$dataset)]


ggplot(data=summary_df,aes(x=KL,y=R2,color=trait))+
  geom_point()+geom_smooth(method="lm",se=F)

summary(lmer(R2~KL+(1|trait)+(1|dataset),data=summary_df))


## bringing in randomized validation data
val_models<-readRDS("SavedResults/all_jack_df_list.rds")
summary_df$R2_val<-NA
summary_df$R2_val[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,data=val_models[["EWT"]]))$r.squared
summary_df$R2_val[summary_df$trait=="LMA"]<-summary(lm(Measured~pred.mean,val_models[["LMA"]]))$r.squared
summary_df$R2_val[summary_df$trait=="LDMC"]<-summary(lm(Measured~pred.mean,val_models[["LDMC"]]))$r.squared
summary_df$R2_val[summary_df$trait=="Nmass"]<-summary(lm(Measured~pred.mean,val_models[["Nmass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="lignin_mass"]<-summary(lm(Measured~pred.mean,val_models[["lignin_mass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="chlA_mass"]<-summary(lm(Measured~pred.mean,val_models[["chlA_mass"]]))$r.squared

ggplot(data=summary_df,aes(x=trait,y=R2))+
  geom_violin()+
  geom_point(data=summary_df,aes(x=trait,y=R2_val),color="red")+
  theme_bw()

##########################################
## model transferral among growth forms

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
                           lower=-0.2,
                           upper=0.6)
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
                           lower=-0.02,
                           upper=0.09)
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
                            lower=-200,
                            upper=800)
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
                             lower=-1,
                             upper=8)
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
                             lower=20,
                             upper=80)
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
                                   lower=0,
                                   upper=100)
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
                                   lower=-5,
                                   upper=50)
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
                                   lower=-5,
                                   upper=50)
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
                                   lower=-5,
                                   upper=40)
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
                                  lower=-5,
                                  upper=30)
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
                                 lower=-5,
                                 upper=15)
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
                                 lower=-5,
                                 upper=10)
  car_mass_KL[[i]]<-KL.divergence(na.omit(meta(spec.traits.other)$car_mass),
                             na.omit(meta(spec.traits.fg)$car_mass))
  
}

LMA_sum_df<-data.frame(functional.group=names(LMA_pred_list),
                       RMSD=unlist(lapply(LMA_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                       R2=unlist(lapply(LMA_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                       hellinger=unlist(LMA_dist),
                       KL=unlist(lapply(LMA_KL,function(x) x[[5]])),
                       trait="LMA")

EWT_sum_df<-data.frame(functional.group=names(EWT_pred_list),
                       RMSD=unlist(lapply(EWT_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                       R2=unlist(lapply(EWT_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                       hellinger=unlist(EWT_dist),
                       KL=unlist(lapply(EWT_KL,function(x) x[[5]])),
                       trait="EWT")

LDMC_sum_df<-data.frame(functional.group=names(LDMC_pred_list),
                        RMSD=unlist(lapply(LDMC_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                        R2=unlist(lapply(LDMC_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                        hellinger=unlist(LDMC_dist),
                        KL=unlist(lapply(LDMC_KL,function(x) x[[5]])),
                        trait="LDMC")

Nmass_sum_df<-data.frame(functional.group=names(Nmass_pred_list),
                         RMSD=unlist(lapply(Nmass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                         R2=unlist(lapply(Nmass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                         hellinger=unlist(Nmass_dist),
                         KL=unlist(lapply(Nmass_KL,function(x) x[[5]])),
                         trait="Nmass")

## not sure why KL divergence for C mass comes up infinite
## for K<8 or so
Cmass_sum_df<-data.frame(functional.group=names(Cmass_pred_list),
                         RMSD=unlist(lapply(Cmass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                         R2=unlist(lapply(Cmass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                         hellinger=unlist(Cmass_dist),
                         KL=unlist(lapply(Cmass_KL,function(x) x[[10]])),
                         trait="Cmass")

solubles_mass_sum_df<-data.frame(functional.group=names(solubles_mass_pred_list),
                                 RMSD=unlist(lapply(solubles_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                                 R2=unlist(lapply(solubles_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                                 hellinger=unlist(solubles_mass_dist),
                                 KL=unlist(lapply(solubles_mass_KL,function(x) x[[5]])),
                                 trait="solubles_mass")

hemicellulose_mass_sum_df<-data.frame(functional.group=names(hemicellulose_mass_pred_list),
                                      RMSD=unlist(lapply(hemicellulose_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                                      R2=unlist(lapply(hemicellulose_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                                      hellinger=unlist(hemicellulose_mass_dist),
                                      KL=unlist(lapply(hemicellulose_mass_KL,function(x) x[[5]])),
                                      trait="hemicellulose_mass")

cellulose_mass_sum_df<-data.frame(functional.group=names(cellulose_mass_pred_list),
                                  RMSD=unlist(lapply(cellulose_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                                  R2=unlist(lapply(cellulose_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                                  hellinger=unlist(cellulose_mass_dist),
                                  KL=unlist(lapply(cellulose_mass_KL,function(x) x[[5]])),
                                  trait="cellulose_mass")

lignin_mass_sum_df<-data.frame(functional.group=names(lignin_mass_pred_list),
                               RMSD=unlist(lapply(lignin_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                               R2=unlist(lapply(lignin_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                               hellinger=unlist(lignin_mass_dist),
                               KL=unlist(lapply(lignin_mass_KL,function(x) x[[5]])),
                               trait="lignin_mass")

chlA_mass_sum_df<-data.frame(functional.group=names(chlA_mass_pred_list),
                             RMSD=unlist(lapply(chlA_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                             R2=unlist(lapply(chlA_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                             hellinger=unlist(chlA_mass_dist),
                             KL=unlist(lapply(chlA_mass_KL,function(x) x[[5]])),
                             trait="chlA_mass")

chlB_mass_sum_df<-data.frame(functional.group=names(chlB_mass_pred_list),
                             RMSD=unlist(lapply(chlB_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                             R2=unlist(lapply(chlB_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                             hellinger=unlist(chlB_mass_dist),
                             KL=unlist(lapply(chlB_mass_KL,function(x) x[[5]])),
                             trait="chlB_mass")

car_mass_sum_df<-data.frame(functional.group=names(car_mass_pred_list),
                            RMSD=unlist(lapply(car_mass_pred_list,function(x) RMSD(x$measured,x$val_pred))),
                            R2=unlist(lapply(car_mass_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2)),
                            hellinger=unlist(car_mass_dist),
                            KL=unlist(lapply(car_mass_KL,function(x) x[[5]])),
                            trait="car_mass")

summary_df<-do.call(rbind,list(LMA_sum_df,LDMC_sum_df,EWT_sum_df,
                               Nmass_sum_df,Cmass_sum_df,solubles_mass_sum_df,
                               hemicellulose_mass_sum_df,cellulose_mass_sum_df,
                               lignin_mass_sum_df,chlA_mass_sum_df,
                               chlB_mass_sum_df,car_mass_sum_df))

## bringing in randomized validation data
val_models<-readRDS("SavedResults/all_jack_df_list.rds")
summary_df$R2_val<-NA
summary_df$R2_val[summary_df$trait=="LMA"]<-summary(lm(Measured~pred.mean,val_models[["LMA"]]))$r.squared
summary_df$R2_val[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,data=val_models[["EWT"]]))$r.squared
summary_df$R2_val[summary_df$trait=="LDMC"]<-summary(lm(Measured~pred.mean,val_models[["LDMC"]]))$r.squared
summary_df$R2_val[summary_df$trait=="Nmass"]<-summary(lm(Measured~pred.mean,val_models[["Nmass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="Cmass"]<-summary(lm(Measured~pred.mean,val_models[["Cmass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="solubles_mass"]<-summary(lm(Measured~pred.mean,val_models[["solubles_mass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="hemicellulose_mass"]<-summary(lm(Measured~pred.mean,val_models[["hemicellulose_mass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="cellulose_mass"]<-summary(lm(Measured~pred.mean,val_models[["cellulose_mass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="lignin_mass"]<-summary(lm(Measured~pred.mean,val_models[["lignin_mass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="chlA_mass"]<-summary(lm(Measured~pred.mean,val_models[["chlA_mass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="chlB_mass"]<-summary(lm(Measured~pred.mean,val_models[["chlB_mass"]]))$r.squared
summary_df$R2_val[summary_df$trait=="car_mass"]<-summary(lm(Measured~pred.mean,val_models[["car_mass"]]))$r.squared

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
# summary_df$R2_ANGERS[summary_df$trait=="LDMC"]<-summary(lm(Measured~pred.mean,
#                                                           LDMC_pred_df[which(LDMC_pred_df$dataset=="ANGERS"),]))$r.squared

EWT_pred_df<-ind_val_list[["EWT"]]
summary_df$R2_Dessain[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,
                                                           EWT_pred_df[which(EWT_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,
                                                         EWT_pred_df[which(EWT_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="EWT"]<-summary(lm(Measured~pred.mean,
                                                          EWT_pred_df[which(EWT_pred_df$dataset=="ANGERS"),]))$r.squared

solubles_mass_pred_df<-ind_val_list[["solubles_mass"]]
summary_df$R2_Dessain[summary_df$trait=="solubles_mass"]<-summary(lm(Measured~pred.mean,
                                                           solubles_mass_pred_df[which(solubles_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="solubles_mass"]<-summary(lm(Measured~pred.mean,
                                                         solubles_mass_pred_df[which(solubles_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="solubles_mass"]<-summary(lm(Measured~pred.mean,
                                                          solubles_mass_pred_df[which(solubles_mass_pred_df$dataset=="ANGERS"),]))$r.squared

hemicellulose_mass_pred_df<-ind_val_list[["hemicellulose_mass"]]
summary_df$R2_Dessain[summary_df$trait=="hemicellulose_mass"]<-summary(lm(Measured~pred.mean,
                                                           hemicellulose_mass_pred_df[which(hemicellulose_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="hemicellulose_mass"]<-summary(lm(Measured~pred.mean,
                                                         hemicellulose_mass_pred_df[which(hemicellulose_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="hemicellulose_mass"]<-summary(lm(Measured~pred.mean,
                                                          hemicellulose_mass_pred_df[which(hemicellulose_mass_pred_df$dataset=="ANGERS"),]))$r.squared

cellulose_mass_pred_df<-ind_val_list[["cellulose_mass"]]
summary_df$R2_Dessain[summary_df$trait=="cellulose_mass"]<-summary(lm(Measured~pred.mean,
                                                           cellulose_mass_pred_df[which(cellulose_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="cellulose_mass"]<-summary(lm(Measured~pred.mean,
                                                         cellulose_mass_pred_df[which(cellulose_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="cellulose_mass"]<-summary(lm(Measured~pred.mean,
                                                          cellulose_mass_pred_df[which(cellulose_mass_pred_df$dataset=="ANGERS"),]))$r.squared

lignin_mass_pred_df<-ind_val_list[["lignin_mass"]]
summary_df$R2_Dessain[summary_df$trait=="lignin_mass"]<-summary(lm(Measured~pred.mean,
                                                           lignin_mass_pred_df[which(lignin_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="lignin_mass"]<-summary(lm(Measured~pred.mean,
                                                         lignin_mass_pred_df[which(lignin_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="lignin_mass"]<-summary(lm(Measured~pred.mean,
                                                          lignin_mass_pred_df[which(lignin_mass_pred_df$dataset=="ANGERS"),]))$r.squared

Nmass_pred_df<-ind_val_list[["Nmass"]]
summary_df$R2_Dessain[summary_df$trait=="Nmass"]<-summary(lm(Measured~pred.mean,
                                                           Nmass_pred_df[which(Nmass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="Nmass"]<-summary(lm(Measured~pred.mean,
                                                         Nmass_pred_df[which(Nmass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="Nmass"]<-summary(lm(Measured~pred.mean,
                                                          Nmass_pred_df[which(Nmass_pred_df$dataset=="ANGERS"),]))$r.squared

Cmass_pred_df<-ind_val_list[["Cmass"]]
summary_df$R2_Dessain[summary_df$trait=="Cmass"]<-summary(lm(Measured~pred.mean,
                                                           Cmass_pred_df[which(Cmass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="Cmass"]<-summary(lm(Measured~pred.mean,
                                                         Cmass_pred_df[which(Cmass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="Cmass"]<-summary(lm(Measured~pred.mean,
                                                          Cmass_pred_df[which(Cmass_pred_df$dataset=="ANGERS"),]))$r.squared

chlA_mass_pred_df<-ind_val_list[["chlA_mass"]]
summary_df$R2_Dessain[summary_df$trait=="chlA_mass"]<-summary(lm(Measured~pred.mean,
                                                           chlA_mass_pred_df[which(chlA_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="chlA_mass"]<-summary(lm(Measured~pred.mean,
                                                         chlA_mass_pred_df[which(chlA_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="chlA_mass"]<-summary(lm(Measured~pred.mean,
                                                          chlA_mass_pred_df[which(chlA_mass_pred_df$dataset=="ANGERS"),]))$r.squared

chlB_mass_pred_df<-ind_val_list[["chlB_mass"]]
summary_df$R2_Dessain[summary_df$trait=="chlB_mass"]<-summary(lm(Measured~pred.mean,
                                                           chlB_mass_pred_df[which(chlB_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="chlB_mass"]<-summary(lm(Measured~pred.mean,
                                                         chlB_mass_pred_df[which(chlB_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="chlB_mass"]<-summary(lm(Measured~pred.mean,
                                                          chlB_mass_pred_df[which(chlB_mass_pred_df$dataset=="ANGERS"),]))$r.squared


car_mass_pred_df<-ind_val_list[["car_mass"]]
summary_df$R2_Dessain[summary_df$trait=="car_mass"]<-summary(lm(Measured~pred.mean,
                                                           car_mass_pred_df[which(car_mass_pred_df$dataset=="Dessain"),]))$r.squared
summary_df$R2_LOPEX[summary_df$trait=="car_mass"]<-summary(lm(Measured~pred.mean,
                                                         car_mass_pred_df[which(car_mass_pred_df$dataset=="LOPEX"),]))$r.squared
summary_df$R2_ANGERS[summary_df$trait=="car_mass"]<-summary(lm(Measured~pred.mean,
                                                          car_mass_pred_df[which(car_mass_pred_df$dataset=="ANGERS"),]))$r.squared

pdf("Images/R2_summary.pdf",height = 6,width=8,onefile=F)
ggplot(data=summary_df,aes(x=trait,y=R2))+
  geom_point(color="blue",size=2)+
  geom_point(data=summary_df,aes(x=trait,y=R2_Dessain),color="black",size=2)+
  geom_point(data=summary_df,aes(x=trait,y=R2_LOPEX),color="black",size=2)+
  geom_point(data=summary_df,aes(x=trait,y=R2_ANGERS),color="black",size=2)+
  geom_point(data=summary_df,aes(x=trait,y=R2_val),color="red",size=2)+
  theme_bw()+
  scale_x_discrete(labels=expression("LMA","LDMC","EWT","N"[mass],
                            "C"[mass],"sol","hemi","cell",
                            "lign","Chl "~italic("a"),
                            "Chl "~italic("b"),"car"))+
  theme(text=element_text(size=20),
        axis.text.x = element_text(angle = 90,hjust=1,vjust=0))+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Trait",y=expression("R"^2))
dev.off()

RMSD_df<-data.frame(functional.group=LMA_sum_df$functional.group,
                    LMA=LMA_sum_df$RMSD,
                    LDMC=LDMC_sum_df$RMSD,
                    Nmass=Nmass_sum_df$RMSD,
                    lignin_mass=lignin_mass_sum_df$RMSD,
                    chlA_mass=chlA_mass_sum_df$RMSD,
                    EWT=EWT_sum_df$RMSD)

R2_df<-data.frame(functional.group=LMA_sum_df$functional.group,
                    LMA=LMA_sum_df$R2,
                    LDMC=LDMC_sum_df$R2,
                    Nmass=Nmass_sum_df$R2,
                    lignin_mass=lignin_mass_sum_df$R2,
                    chlA_mass=chlA_mass_sum_df$R2,
                    EWT=EWT_sum_df$R2)

ggplot(data=summary_df,aes(x=hellinger,y=R2,
                           color=trait))+
  geom_point(size=2,aes(shape=functional.group))+
  geom_smooth(method="lm",se=F)
