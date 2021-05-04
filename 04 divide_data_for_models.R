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

## this should return integer(0) if the two data sets have
## spectra for the sample samples
which(!(meta(ref.traits)$sample_id %in% meta(trans.traits)$sample_id))

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

## save divided up data
saveRDS(ref.train,"ProcessedSpectra/ref_train.rds")
saveRDS(ref.test,"ProcessedSpectra/ref_test.rds")
saveRDS(trans.train,"ProcessedSpectra/trans_train.rds")
saveRDS(trans.test,"ProcessedSpectra/trans_test.rds")