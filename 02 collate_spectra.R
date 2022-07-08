setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels/")

library(spectrolab)

##################################
## read in and process reflectance data

BeauchampRioux.ref.df<-read.csv("ProcessedSpectra/BeauchampRioux_ref_processed.csv")
colnames(BeauchampRioux.ref.df)<-gsub("X","",colnames(BeauchampRioux.ref.df))
BeauchampRioux.ref.spec<-as_spectra(BeauchampRioux.ref.df,name_idx=1)
meta(BeauchampRioux.ref.spec)$sample_id<-gsub("spec_","",names(BeauchampRioux.ref.spec))
BeauchampRioux.ref.spec<-BeauchampRioux.ref.spec[,400:2400]

Blanchard.ref.df<-read.csv("ProcessedSpectra/Blanchard_ref_processed.csv")
colnames(Blanchard.ref.df)<-gsub("X","",colnames(Blanchard.ref.df))
Blanchard.ref.spec<-as_spectra(Blanchard.ref.df,name_idx=1)
meta(Blanchard.ref.spec)$sample_id<-gsub("spec_","",names(Blanchard.ref.spec))
Blanchard.ref.spec<-Blanchard.ref.spec[,400:2400]

Boucherville2018.ref.df<-read.csv("ProcessedSpectra/Boucherville2018_ref_processed.csv")
colnames(Boucherville2018.ref.df)<-gsub("X","",colnames(Boucherville2018.ref.df))
Boucherville2018.ref.spec<-as_spectra(Boucherville2018.ref.df,name_idx=1)
meta(Boucherville2018.ref.spec)$sample_id<-gsub("spec_","",names(Boucherville2018.ref.spec))
Boucherville2018.ref.spec<-Boucherville2018.ref.spec[,400:2400]

Boucherville2019.ref.df<-read.csv("ProcessedSpectra/Boucherville2019_ref_processed.csv")
colnames(Boucherville2019.ref.df)<-gsub("X","",colnames(Boucherville2019.ref.df))
Boucherville2019.ref.spec<-as_spectra(Boucherville2019.ref.df,name_idx=1)
meta(Boucherville2019.ref.spec)$sample_id<-gsub("spec_","",names(Boucherville2019.ref.spec))
Boucherville2019.ref.spec<-Boucherville2019.ref.spec[,400:2400]

CABOGeneral2019.ref.df<-read.csv("ProcessedSpectra/CABOGeneral2019_ref_processed.csv")
colnames(CABOGeneral2019.ref.df)<-gsub("X","",colnames(CABOGeneral2019.ref.df))
CABOGeneral2019.ref.spec<-as_spectra(CABOGeneral2019.ref.df,name_idx=1)
meta(CABOGeneral2019.ref.spec)$sample_id<-gsub("spec_","",names(CABOGeneral2019.ref.spec))
CABOGeneral2019.ref.spec<-CABOGeneral2019.ref.spec[,400:2400]
CABOGeneral2019.ref.spec<-CABOGeneral2019.ref.spec[-which(names(CABOGeneral2019.ref.spec) %in% names(Blanchard.ref.spec))]

CABOGeneralOther.ref.df<-read.csv("ProcessedSpectra/CABOGeneralOther_ref_processed.csv")
colnames(CABOGeneralOther.ref.df)<-gsub("X","",colnames(CABOGeneralOther.ref.df))
CABOGeneralOther.ref.spec<-as_spectra(CABOGeneralOther.ref.df,name_idx=1)
meta(CABOGeneralOther.ref.spec)$sample_id<-gsub("spec_","",names(CABOGeneralOther.ref.spec))
CABOGeneralOther.ref.spec<-CABOGeneralOther.ref.spec[,400:2400]

Crofts.ref.df<-read.csv("ProcessedSpectra/Crofts_ref_processed.csv")
colnames(Crofts.ref.df)<-gsub("X","",colnames(Crofts.ref.df))
Crofts.ref.spec<-as_spectra(Crofts.ref.df,name_idx=1)
meta(Crofts.ref.spec)$sample_id<-gsub("spec_","",names(Crofts.ref.spec))
Crofts.ref.spec<-Crofts.ref.spec[,400:2400]

Girard.ref.df<-read.csv("ProcessedSpectra/Girard_ref_processed.csv")
colnames(Girard.ref.df)<-gsub("X","",colnames(Girard.ref.df))
Girard.ref.spec<-as_spectra(Girard.ref.df,name_idx=1)
meta(Girard.ref.spec)$sample_id<-gsub("spec_","",names(Girard.ref.spec))
Girard.ref.spec<-Girard.ref.spec[,400:2400]

Hacker2019.ref.df<-read.csv("ProcessedSpectra/Hacker2019_ref_processed.csv")
colnames(Hacker2019.ref.df)<-gsub("X","",colnames(Hacker2019.ref.df))
Hacker2019.ref.spec<-as_spectra(Hacker2019.ref.df,name_idx=1)
meta(Hacker2019.ref.spec)$sample_id<-gsub("spec_","",names(Hacker2019.ref.spec))
Hacker2019.ref.spec<-Hacker2019.ref.spec[,400:2400]

Pardo.ref.df<-read.csv("ProcessedSpectra/Pardo_ref_processed.csv")
colnames(Pardo.ref.df)<-gsub("X","",colnames(Pardo.ref.df))
Pardo.ref.spec<-as_spectra(Pardo.ref.df,name_idx=1)
meta(Pardo.ref.spec)$sample_id<-gsub("spec_","",names(Pardo.ref.spec))
Pardo.ref.spec<-Pardo.ref.spec[,400:2400]

PhragmitesTemporal.ref.df<-read.csv("ProcessedSpectra/PhragmitesTemporal_ref_processed.csv")
colnames(PhragmitesTemporal.ref.df)<-gsub("X","",colnames(PhragmitesTemporal.ref.df))
PhragmitesTemporal.ref.spec<-as_spectra(PhragmitesTemporal.ref.df,name_idx=1)
meta(PhragmitesTemporal.ref.spec)$sample_id<-gsub("spec_","",names(PhragmitesTemporal.ref.spec))
PhragmitesTemporal.ref.spec<-PhragmitesTemporal.ref.spec[,400:2400]

Warren.ref.df<-read.csv("ProcessedSpectra/Warren_ref_processed.csv")
colnames(Warren.ref.df)<-gsub("X","",colnames(Warren.ref.df))
Warren.ref.spec<-as_spectra(Warren.ref.df,name_idx=1)
meta(Warren.ref.spec)$sample_id<-gsub("spec_","",names(Warren.ref.spec))
Warren.ref.spec<-Warren.ref.spec[,400:2400]

all.ref.spec<-Reduce(combine, list(BeauchampRioux.ref.spec, Blanchard.ref.spec, Boucherville2018.ref.spec,
                               Boucherville2019.ref.spec, CABOGeneral2019.ref.spec, CABOGeneralOther.ref.spec,
                               Crofts.ref.spec, Girard.ref.spec, Hacker2019.ref.spec,
                               Pardo.ref.spec, PhragmitesTemporal.ref.spec, Warren.ref.spec))

#####################################################
## read in and process transmittance data

BeauchampRioux.trans.df<-read.csv("ProcessedSpectra/BeauchampRioux_trans_processed.csv")
colnames(BeauchampRioux.trans.df)<-gsub("X","",colnames(BeauchampRioux.trans.df))
BeauchampRioux.trans.spec<-as_spectra(BeauchampRioux.trans.df,name_idx=1)
meta(BeauchampRioux.trans.spec)$sample_id<-gsub("spec_","",names(BeauchampRioux.trans.spec))
BeauchampRioux.trans.spec<-BeauchampRioux.trans.spec[,400:2400]

Blanchard.trans.df<-read.csv("ProcessedSpectra/Blanchard_trans_processed.csv")
colnames(Blanchard.trans.df)<-gsub("X","",colnames(Blanchard.trans.df))
Blanchard.trans.spec<-as_spectra(Blanchard.trans.df,name_idx=1)
meta(Blanchard.trans.spec)$sample_id<-gsub("spec_","",names(Blanchard.trans.spec))
Blanchard.trans.spec<-Blanchard.trans.spec[,400:2400]

Boucherville2018.trans.df<-read.csv("ProcessedSpectra/Boucherville2018_trans_processed.csv")
colnames(Boucherville2018.trans.df)<-gsub("X","",colnames(Boucherville2018.trans.df))
Boucherville2018.trans.spec<-as_spectra(Boucherville2018.trans.df,name_idx=1)
meta(Boucherville2018.trans.spec)$sample_id<-gsub("spec_","",names(Boucherville2018.trans.spec))
Boucherville2018.trans.spec<-Boucherville2018.trans.spec[,400:2400]

Boucherville2019.trans.df<-read.csv("ProcessedSpectra/Boucherville2019_trans_processed.csv")
colnames(Boucherville2019.trans.df)<-gsub("X","",colnames(Boucherville2019.trans.df))
Boucherville2019.trans.spec<-as_spectra(Boucherville2019.trans.df,name_idx=1)
meta(Boucherville2019.trans.spec)$sample_id<-gsub("spec_","",names(Boucherville2019.trans.spec))
Boucherville2019.trans.spec<-Boucherville2019.trans.spec[,400:2400]

CABOGeneral2019.trans.df<-read.csv("ProcessedSpectra/CABOGeneral2019_trans_processed.csv")
colnames(CABOGeneral2019.trans.df)<-gsub("X","",colnames(CABOGeneral2019.trans.df))
CABOGeneral2019.trans.spec<-as_spectra(CABOGeneral2019.trans.df,name_idx=1)
meta(CABOGeneral2019.trans.spec)$sample_id<-gsub("spec_","",names(CABOGeneral2019.trans.spec))
CABOGeneral2019.trans.spec<-CABOGeneral2019.trans.spec[,400:2400]
CABOGeneral2019.trans.spec<-CABOGeneral2019.trans.spec[-which(names(CABOGeneral2019.trans.spec) %in% names(Blanchard.trans.spec))]

CABOGeneralOther.trans.df<-read.csv("ProcessedSpectra/CABOGeneralOther_trans_processed.csv")
colnames(CABOGeneralOther.trans.df)<-gsub("X","",colnames(CABOGeneralOther.trans.df))
CABOGeneralOther.trans.spec<-as_spectra(CABOGeneralOther.trans.df,name_idx=1)
meta(CABOGeneralOther.trans.spec)$sample_id<-gsub("spec_","",names(CABOGeneralOther.trans.spec))
CABOGeneralOther.trans.spec<-CABOGeneralOther.trans.spec[,400:2400]

Crofts.trans.df<-read.csv("ProcessedSpectra/Crofts_trans_processed.csv")
colnames(Crofts.trans.df)<-gsub("X","",colnames(Crofts.trans.df))
Crofts.trans.spec<-as_spectra(Crofts.trans.df,name_idx=1)
meta(Crofts.trans.spec)$sample_id<-gsub("spec_","",names(Crofts.trans.spec))
Crofts.trans.spec<-Crofts.trans.spec[,400:2400]

Girard.trans.df<-read.csv("ProcessedSpectra/Girard_trans_processed.csv")
colnames(Girard.trans.df)<-gsub("X","",colnames(Girard.trans.df))
Girard.trans.spec<-as_spectra(Girard.trans.df,name_idx=1)
meta(Girard.trans.spec)$sample_id<-gsub("spec_","",names(Girard.trans.spec))
Girard.trans.spec<-Girard.trans.spec[,400:2400]

Hacker2019.trans.df<-read.csv("ProcessedSpectra/Hacker2019_trans_processed.csv")
colnames(Hacker2019.trans.df)<-gsub("X","",colnames(Hacker2019.trans.df))
Hacker2019.trans.spec<-as_spectra(Hacker2019.trans.df,name_idx=1)
meta(Hacker2019.trans.spec)$sample_id<-gsub("spec_","",names(Hacker2019.trans.spec))
Hacker2019.trans.spec<-Hacker2019.trans.spec[,400:2400]

Pardo.trans.df<-read.csv("ProcessedSpectra/Pardo_trans_processed.csv")
colnames(Pardo.trans.df)<-gsub("X","",colnames(Pardo.trans.df))
Pardo.trans.spec<-as_spectra(Pardo.trans.df,name_idx=1)
meta(Pardo.trans.spec)$sample_id<-gsub("spec_","",names(Pardo.trans.spec))
Pardo.trans.spec<-Pardo.trans.spec[,400:2400]

PhragmitesTemporal.trans.df<-read.csv("ProcessedSpectra/PhragmitesTemporal_trans_processed.csv")
colnames(PhragmitesTemporal.trans.df)<-gsub("X","",colnames(PhragmitesTemporal.trans.df))
PhragmitesTemporal.trans.spec<-as_spectra(PhragmitesTemporal.trans.df,name_idx=1)
meta(PhragmitesTemporal.trans.spec)$sample_id<-gsub("spec_","",names(PhragmitesTemporal.trans.spec))
PhragmitesTemporal.trans.spec<-PhragmitesTemporal.trans.spec[,400:2400]

Warren.trans.df<-read.csv("ProcessedSpectra/Warren_trans_processed.csv")
colnames(Warren.trans.df)<-gsub("X","",colnames(Warren.trans.df))
Warren.trans.spec<-as_spectra(Warren.trans.df,name_idx=1)
meta(Warren.trans.spec)$sample_id<-gsub("spec_","",names(Warren.trans.spec))
Warren.trans.spec<-Warren.trans.spec[,400:2400]

all.trans.spec<-Reduce(combine, list(BeauchampRioux.trans.spec, Blanchard.trans.spec, Boucherville2018.trans.spec,
                                   Boucherville2019.trans.spec, CABOGeneral2019.trans.spec, CABOGeneralOther.trans.spec,
                                   Crofts.trans.spec, Girard.trans.spec, Hacker2019.trans.spec,
                                   Pardo.trans.spec, PhragmitesTemporal.trans.spec, Warren.trans.spec))

#####################################################
## attach summary data

Fulcrum.summary<-read.csv("SummaryData/leaf_spectra.csv")
Fulcrum.summary$measurement.date<-as.POSIXlt(Fulcrum.summary$date_measured,format="%Y-%m-%d")
Fulcrum.sub<-data.frame(sample.id=Fulcrum.summary$sample_id,
                        species=Fulcrum.summary$scientific_name,
                        project=Fulcrum.summary$project,
                        site=Fulcrum.summary$site_id,
                        measurement.date=Fulcrum.summary$measurement.date)

Warren.summary<-read.csv("SummaryData/WarrenSummary.csv")
Fulcrum.sub$stage<-Warren.summary$Stage[match(Fulcrum.summary$sample_id,Warren.summary$Bulk.sample.ID)]

site.summary<-read.csv("SummaryData/sites.csv")
Fulcrum.sub$latitude<-site.summary$latitude[match(Fulcrum.sub$site,site.summary$site_id)]
Fulcrum.sub$longitude<-site.summary$longitude[match(Fulcrum.sub$site,site.summary$site_id)]

sp_split<-strsplit(as.character(Fulcrum.sub$species),split=" ")
Fulcrum.sub$latin.genus<-unlist(lapply(sp_split,function(entry) entry[[1]]))
Fulcrum.sub$latin.species<-unlist(lapply(sp_split,function(entry) entry[[2]]))

## missing full summary data for: Pardo (but we'll drop those samples later)

meta(all.ref.spec)$species<-Fulcrum.sub$species[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.ref.spec)$latin.genus<-Fulcrum.sub$latin.genus[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.ref.spec)$latin.species<-Fulcrum.sub$latin.species[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.ref.spec)$measurement.date<-Fulcrum.sub$measurement.date[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.ref.spec)$project<-Fulcrum.sub$project[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.ref.spec)$site<-Fulcrum.sub$site[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.ref.spec)$stage<-Fulcrum.sub$stage[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.ref.spec)$latitude<-Fulcrum.sub$latitude[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.ref.spec)$longitude<-Fulcrum.sub$longitude[match(meta(all.ref.spec)$sample_id,Fulcrum.sub$sample.id)]

meta(all.trans.spec)$species<-Fulcrum.sub$species[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.trans.spec)$latin.genus<-Fulcrum.sub$latin.genus[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.trans.spec)$latin.species<-Fulcrum.sub$latin.species[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.trans.spec)$measurement.date<-Fulcrum.sub$measurement.date[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.trans.spec)$project<-Fulcrum.sub$project[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.trans.spec)$site<-Fulcrum.sub$site[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.trans.spec)$stage<-Fulcrum.sub$stage[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.trans.spec)$latitude<-Fulcrum.sub$latitude[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]
meta(all.trans.spec)$longitude<-Fulcrum.sub$longitude[match(meta(all.trans.spec)$sample_id,Fulcrum.sub$sample.id)]

## remove certain data from Beauchamp Rioux dataset
## flagged for having red or yellow leaves (or being otherwise non-standard)
all.ref.spec<-all.ref.spec[-which(meta(all.ref.spec)$sample_id %in% c("11851575","13221003","22405244",
                                                          "21854888","21854774","21854585",
                                                          "10449089","10450505","10452131",
                                                          "10453164","10454172","10461643",
                                                          "10465061","10466192","10466735",
                                                          "10468015"))]

all.trans.spec<-all.trans.spec[-which(meta(all.trans.spec)$sample_id %in% c("11851575","13221003","22405244",
                                                                      "21854888","21854774","21854585",
                                                                      "10449089","10450505","10452131",
                                                                      "10453164","10454172","10461643",
                                                                      "10465061","10466192","10466735",
                                                                      "10468015"))]

## remove a few missing spectra from Blanchard project?
## not sure why they're missing though!
all.ref.spec<-all.ref.spec[-which(is.na(rowSums(as.matrix(all.ref.spec))))]
all.trans.spec<-all.trans.spec[-which(is.na(rowSums(as.matrix(all.trans.spec))))]

## check that spectra are in the same order
sum(meta(all.ref.spec)$sample_id==meta(all.trans.spec)$sample_id)==nrow(all.ref.spec)

all.abs.mat<-1-as.matrix(all.ref.spec)-as.matrix(all.trans.spec)
all.abs.spec<-spectra(all.abs.mat,bands=400:2400,names=names(all.ref.spec))
meta(all.abs.spec)<-meta(all.ref.spec)

saveRDS(all.ref.spec,file = "ProcessedSpectra/all_ref.rds")
saveRDS(all.trans.spec,file = "ProcessedSpectra/all_trans.rds")
saveRDS(all.abs.spec,file = "ProcessedSpectra/all_abs.rds")
