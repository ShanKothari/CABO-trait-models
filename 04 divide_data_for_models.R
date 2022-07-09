setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(caret)

#########################################
## this script is meant to ensure that the
## training and testing data for PLSR trait
## modeling analyses using reflectance,
## transmittance, and absorptance are the same

## read spectra and traits
ref.traits<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")
trans.traits<-readRDS("ProcessedSpectra/all_trans_and_traits.rds")
abs.traits<-readRDS("ProcessedSpectra/all_abs_and_traits.rds")

## remove Pardo project, which currently has no trait or summary data
ref.traits<-ref.traits[which(meta(ref.traits)$project!="2019-Pardo-MSc-UdeM")]
trans.traits<-trans.traits[which(meta(trans.traits)$project!="2019-Pardo-MSc-UdeM")]
abs.traits<-abs.traits[which(meta(abs.traits)$project!="2019-Pardo-MSc-UdeM")]

## save data for archiving
# write.csv(as.data.frame(ref.traits),"ProcessedSpectra/ref_spec.csv",row.names = F)
# write.csv(as.data.frame(trans.traits),"ProcessedSpectra/trans_spec.csv",row.names = F)
# write.csv(as.data.frame(abs.traits),"ProcessedSpectra/abs_spec.csv",row.names = F)

###################################
## divide up data

## check that IDs are the same and in the same order
if(meta(ref.traits)$sample_id != meta(trans.traits)$sample_id ||
   meta(ref.traits)$sample_id != meta(abs.traits)$sample_id){
  print("sample ids not the same")
}

## create division
train.sample <- createDataPartition(
  y = meta(ref.traits)$functional.group,
  p = .75,
  list = FALSE
)

test.sample<-setdiff(1:nrow(as.matrix(ref.traits)),train.sample)

ref.train<-ref.traits[train.sample,]
ref.test<-ref.traits[test.sample,]

trans.train<-trans.traits[train.sample,]
trans.test<-trans.traits[test.sample,]

abs.train<-abs.traits[train.sample,]
abs.test<-abs.traits[test.sample,]

## save divided up data
saveRDS(ref.train,"ProcessedSpectra/ref_train.rds")
saveRDS(ref.test,"ProcessedSpectra/ref_test.rds")
saveRDS(trans.train,"ProcessedSpectra/trans_train.rds")
saveRDS(trans.test,"ProcessedSpectra/trans_test.rds")
saveRDS(abs.train,"ProcessedSpectra/abs_train.rds")
saveRDS(abs.test,"ProcessedSpectra/abs_test.rds")
