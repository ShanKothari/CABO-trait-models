setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(caret)
library(pls)
library(lmodel2)
library(reshape2)
library(statip)
## library(hypervolume)

spec.traits<-readRDS("ProcessedSpectra/all_spectra_and_traits.rds")
## the Pardo dataset has no trait data (yet)
spec.traits<-spec.traits[-which(meta(spec.traits)$project=="2019-Pardo-MSc-UdeM")]
## removing the Dessain dataset since some of these analyses
## might be sensitive to the exact errors involved
## and Dessain was measured with a different instrument
spec.traits<-spec.traits[-which(meta(spec.traits)$project=="2017-Dessain-MSc")]

#########################################
## define RMSD function

RMSD<-function(measured,predicted){
  not.na<-which(!is.na(measured) & !is.na(predicted))
  return(sqrt(sum((measured-predicted)^2,na.rm=T)/(length(not.na)-1)))
}

########################################
## spectral overlap analyses

## find projects with more than 50 samples (all but CABO-General)
common.proj<-names(table(meta(spec.traits)$project))[table(meta(spec.traits)$project)>50]

spec_overlap<-list()
bands_select<-c(531,680,760,1000,1440,1790)

for(i in common.proj){
  print(i)
  spec.traits.proj<-spec.traits[meta(spec.traits)$project==i]
  spec.traits.other<-spec.traits[meta(spec.traits)$project!=i]
  
  ## hypervolume approach
  spec_proj_KDE<-hypervolume(as.matrix(spec.traits.proj[,bands_select]))
  spec_other_KDE<-hypervolume(as.matrix(spec.traits.other[,bands_select]))
  spec_set<-hypervolume_set(spec_proj_KDE,spec_other_KDE,check.memory = F)
  spec_overlap[[i]]<-hypervolume_overlap_statistics(spec_set)
}

spec_jaccard<-unlist(lapply(spec_overlap,function(x) x[["jaccard"]]))
spec_sorensen<-unlist(lapply(spec_overlap,function(x) x[["sorensen"]]))
spec_unique_frac<-unlist(lapply(spec_overlap,function(x) x[["frac_unique_1"]]))

#######################################
## trait overlap analyses

trait_overlap<-list()
traits_select<-c("Nmass","EWT","lignin_mass","chlA_fresh","LMA")

for(i in common.proj){
  print(i)
  traits.proj<-na.omit(as.data.frame(meta(spec.traits)[,traits_select])[meta(spec.traits)$project==i,])
  traits.other<-na.omit(as.data.frame(meta(spec.traits)[,traits_select])[meta(spec.traits)$project!=i,])
  
  ## hypervolume approach
  trait_proj_KDE<-hypervolume(traits.proj)
  trait_other_KDE<-hypervolume(traits.other)
  trait_set<-hypervolume_set(trait_proj_KDE,trait_other_KDE,check.memory = F)
  trait_overlap[[i]]<-hypervolume_overlap_statistics(trait_set)
}

trait_jaccard<-unlist(lapply(trait_overlap,function(x) x[["jaccard"]]))
trait_sorensen<-unlist(lapply(trait_overlap,function(x) x[["sorensen"]]))
trait_unique_frac<-unlist(lapply(trait_overlap,function(x) x[["frac_unique_1"]]))

########################################
## leave one project out' trait prediction analyses

LDMC_pred_list<-list()

for(i in common.proj){
  print(i)
  spec.traits.proj<-spec.traits[meta(spec.traits)$project==i]
  spec.traits.other<-spec.traits[meta(spec.traits)$project!=i]
  
  LDMC_cal<-plsr(meta(spec.traits.other)$LDMC~as.matrix(spec.traits.other),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)
  ncomp_LDMC<-selectNcomp(LDMC_cal, method = "onesigma", plot = F)
  LDMC_val<-predict(LDMC_cal,newdata=as.matrix(spec.traits.proj),ncomp=ncomp_LDMC)
  LDMC_pred_list[[i]]<-data.frame(measured=meta(spec.traits.proj)$LDMC,
                                 val_pred=as.vector(LDMC_val))

}

LDMC_RMSD_pred<-unlist(lapply(LDMC_pred_list,function(x) RMSD(x$measured,x$val_pred)))
LDMC_R2_pred<-unlist(lapply(LDMC_pred_list,function(x) cor(x$measured,x$val_pred,use="c")^2))
