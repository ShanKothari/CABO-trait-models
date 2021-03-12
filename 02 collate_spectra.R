setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels/")

library(spectrolab)

BeauchampRioux.df<-read.csv("ProcessedSpectra/BeauchampRioux_spec_processed.csv")
colnames(BeauchampRioux.df)<-gsub("X","",colnames(BeauchampRioux.df))
BeauchampRioux.spec<-as_spectra(BeauchampRioux.df,name_idx=1)
meta(BeauchampRioux.spec)$sample_id<-gsub("spec_","",names(BeauchampRioux.spec))
BeauchampRioux.spec<-BeauchampRioux.spec[,400:2400]

Blanchard.df<-read.csv("ProcessedSpectra/Blanchard_spec_processed.csv")
colnames(Blanchard.df)<-gsub("X","",colnames(Blanchard.df))
Blanchard.spec<-as_spectra(Blanchard.df,name_idx=1)
meta(Blanchard.spec)$sample_id<-gsub("spec_","",names(Blanchard.spec))
Blanchard.spec<-Blanchard.spec[,400:2400]

Boucherville2018.df<-read.csv("ProcessedSpectra/Boucherville2018_spec_processed.csv")
colnames(Boucherville2018.df)<-gsub("X","",colnames(Boucherville2018.df))
Boucherville2018.spec<-as_spectra(Boucherville2018.df,name_idx=1)
meta(Boucherville2018.spec)$sample_id<-gsub("spec_","",names(Boucherville2018.spec))
Boucherville2018.spec<-Boucherville2018.spec[,400:2400]

Boucherville2019.df<-read.csv("ProcessedSpectra/Boucherville2019_spec_processed.csv")
colnames(Boucherville2019.df)<-gsub("X","",colnames(Boucherville2019.df))
Boucherville2019.spec<-as_spectra(Boucherville2019.df,name_idx=1)
meta(Boucherville2019.spec)$sample_id<-gsub("spec_","",names(Boucherville2019.spec))
Boucherville2019.spec<-Boucherville2019.spec[,400:2400]

CABOGeneral2019.df<-read.csv("ProcessedSpectra/CABOGeneral2019_spec_processed.csv")
colnames(CABOGeneral2019.df)<-gsub("X","",colnames(CABOGeneral2019.df))
CABOGeneral2019.spec<-as_spectra(CABOGeneral2019.df,name_idx=1)
meta(CABOGeneral2019.spec)$sample_id<-gsub("spec_","",names(CABOGeneral2019.spec))
CABOGeneral2019.spec<-CABOGeneral2019.spec[,400:2400]
CABOGeneral2019.spec<-CABOGeneral2019.spec[-which(names(CABOGeneral2019.spec) %in% names(Blanchard.spec))]

CABOGeneralOther.df<-read.csv("ProcessedSpectra/CABOGeneralOther_spec_processed.csv")
colnames(CABOGeneralOther.df)<-gsub("X","",colnames(CABOGeneralOther.df))
CABOGeneralOther.spec<-as_spectra(CABOGeneralOther.df,name_idx=1)
meta(CABOGeneralOther.spec)$sample_id<-gsub("spec_","",names(CABOGeneralOther.spec))
CABOGeneralOther.spec<-CABOGeneralOther.spec[,400:2400]

Crofts.df<-read.csv("ProcessedSpectra/Crofts_spec_processed.csv")
colnames(Crofts.df)<-gsub("X","",colnames(Crofts.df))
Crofts.spec<-as_spectra(Crofts.df,name_idx=1)
meta(Crofts.spec)$sample_id<-gsub("spec_","",names(Crofts.spec))
Crofts.spec<-Crofts.spec[,400:2400]

Dessain.df<-read.csv("ProcessedSpectra/Dessain_spec_processed.csv")
colnames(Dessain.df)<-gsub("X","",colnames(Dessain.df))
Dessain.spec<-as_spectra(Dessain.df,name_idx=1)
meta(Dessain.spec)$sample_id<-gsub("spec_","",names(Dessain.spec))
Dessain.spec<-Dessain.spec[,400:2400]

Girard.df<-read.csv("ProcessedSpectra/Girard_spec_processed.csv")
colnames(Girard.df)<-gsub("X","",colnames(Girard.df))
Girard.spec<-as_spectra(Girard.df,name_idx=1)
meta(Girard.spec)$sample_id<-gsub("spec_","",names(Girard.spec))
Girard.spec<-Girard.spec[,400:2400]

Hacker.df<-read.csv("ProcessedSpectra/Hacker_spec_processed.csv")
colnames(Hacker.df)<-gsub("X","",colnames(Hacker.df))
Hacker.spec<-as_spectra(Hacker.df,name_idx=1)
meta(Hacker.spec)$sample_id<-gsub("spec_","",names(Hacker.spec))
Hacker.spec<-Hacker.spec[,400:2400]

Pardo.df<-read.csv("ProcessedSpectra/Pardo_spec_processed.csv")
colnames(Pardo.df)<-gsub("X","",colnames(Pardo.df))
Pardo.spec<-as_spectra(Pardo.df,name_idx=1)
meta(Pardo.spec)$sample_id<-gsub("spec_","",names(Pardo.spec))
Pardo.spec<-Pardo.spec[,400:2400]

PhragmitesTemporal.df<-read.csv("ProcessedSpectra/PhragmitesTemporal_spec_processed.csv")
colnames(PhragmitesTemporal.df)<-gsub("X","",colnames(PhragmitesTemporal.df))
PhragmitesTemporal.spec<-as_spectra(PhragmitesTemporal.df,name_idx=1)
meta(PhragmitesTemporal.spec)$sample_id<-gsub("spec_","",names(PhragmitesTemporal.spec))
PhragmitesTemporal.spec<-PhragmitesTemporal.spec[,400:2400]

Warren.df<-read.csv("ProcessedSpectra/Warren_spec_processed.csv")
colnames(Warren.df)<-gsub("X","",colnames(Warren.df))
Warren.spec<-as_spectra(Warren.df,name_idx=1)
meta(Warren.spec)$sample_id<-gsub("spec_","",names(Warren.spec))
Warren.spec<-Warren.spec[,400:2400]

all.spec<-Reduce(combine, list(BeauchampRioux.spec, Blanchard.spec, Boucherville2018.spec,
                               Boucherville2019.spec, CABOGeneral2019.spec, CABOGeneralOther.spec,
                               Crofts.spec, Dessain.spec, Girard.spec, Hacker.spec,
                               Pardo.spec, PhragmitesTemporal.spec, Warren.spec))

Fulcrum.summary<-read.csv("SummaryData/leaf_spectra.csv")
Fulcrum.sub<-data.frame(sample.id=Fulcrum.summary$sample_id,
                        species=Fulcrum.summary$scientific_name,
                        project=Fulcrum.summary$project)
PhragmitesTemporal.summary<-read.csv("SummaryData/leaf_spectra_phragmites_temporal.csv")
PhragmitesTemporal.sub<-data.frame(sample.id=PhragmitesTemporal.summary$sample_id,
                        species=PhragmitesTemporal.summary$scientific_name,
                        project=PhragmitesTemporal.summary$project)
Dessain.summary<-read.csv("SummaryData/DessainHerbarium.csv")
Dessain.sub<-data.frame(sample.id=Dessain.summary$parentEventID,
                        species=Dessain.summary$scientificName,
                        project="2017-Dessain-MSc")
Warren.summary<-read.csv("SummaryData/WarrenSummary.csv")
Warren.sub<-data.frame(sample.id=Warren.summary$Bulk.sample.ID,
                       species=Warren.summary$Species,
                       project="SWA-Warren")
allmeta.sub<-do.call(rbind,args=list(Fulcrum.sub,Dessain.sub,
                                  PhragmitesTemporal.sub,Warren.sub))

## missing full summary data for: Pardo

meta(all.spec)$species<-allmeta.sub$species[match(meta(all.spec)$sample_id,allmeta.sub$sample.id)]
meta(all.spec)$project<-allmeta.sub$project[match(meta(all.spec)$sample_id,allmeta.sub$sample.id)]

## this is temporary
meta(all.spec)$project<-as.character(meta(all.spec)$project)
meta(all.spec)$project[meta(all.spec)$sample_id %in% meta(Pardo.spec)$sample_id]<-"2019-Pardo-MSc-UdeM"

## remove data from Beauchamp Rioux dataset
## of red or yellow leaves (or otherwise non-standard)
all.spec<-all.spec[-which(meta(all.spec)$sample_id %in% c("11851575","13221003","22405244",
                                                                   "21854888","21854774","21854585",
                                                                   "10449089","10450505","10452131",
                                                                   "10453164","10454172","10461643",
                                                                   "10465061","10466192","10466735",
                                                                   "10468015"))]

## remove a few missing spectra from Blanchard project?
## not sure why they're missing though!
all.spec<-all.spec[-which(is.na(rowSums(as.matrix(all.spec))))]

saveRDS(all.spec,file = "ProcessedSpectra/all_spectra.rds")
