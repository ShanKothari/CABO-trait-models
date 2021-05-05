setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(caret)

#########################################
## this script is meant to ensure that the
## training and testing data for PLSR trait
## modeling analyses using reflectance and
## transmittance are the same

## read spectra and traits
ref.traits<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")
trans.traits<-readRDS("ProcessedSpectra/all_trans_and_traits.rds")
abs.traits<-readRDS("ProcessedSpectra/all_abs_and_traits.rds")

## check that IDs are the same and in the same order
if(meta(all.ref)$sample_id != meta(all.trans)$sample_id ||
   meta(all.ref)$sample_id != meta(all.abs)$sample_id){
  stop("sample ids not the same")
}

## divide up data
train.sample <- createDataPartition(
  y = meta(ref.traits)$project,
  p = .8,
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
