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
with(meta(ref.train),median(EWT_actual/EWT,na.rm=T))
