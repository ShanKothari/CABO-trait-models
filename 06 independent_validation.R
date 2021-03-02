setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)

################################################
## LOPEX

LOPEX_ref<-t(read.csv("IndependentValidationData/LOPEX/lopex1993reflectance.csv",row.names = "Sample_."))
LOPEX_ref<-spectra(LOPEX_ref,
                   bands=colnames(LOPEX_ref),
                   names=rownames(LOPEX_ref))

LOPEX_traits<-read.csv("IndependentValidationData/LOPEX/lopex1993metadata.csv")
LOPEX_traits[which(LOPEX_traits== -999,arr.ind=T)]<-NA
meta(LOPEX_ref)$Nmass<-LOPEX_traits$Nitrogen....DW.
meta(LOPEX_ref)$Cmass<-LOPEX_traits$Carbon....DW.
meta(LOPEX_ref)$EWT<-LOPEX_traits$Equivalent.Water.Thickness..g.cm2.
meta(LOPEX_ref)$LMA<-(LOPEX_traits$Dry.Weight..g./1000)/(LOPEX_traits$Leaf.Area..cm2./10000)
meta(LOPEX_ref)$LDMC<-(LOPEX_traits$Dry.Weight..g./LOPEX_traits$Fresh.Weight..g.)*1000
meta(LOPEX_ref)$cellulose1<-LOPEX_traits$Cellulose_1....DW.
meta(LOPEX_ref)$cellulose2<-LOPEX_traits$Cellulose_2....DW.
meta(LOPEX_ref)$dataset<-"LOPEX"
## measures of lignin are highly inconsistent, so we ignore them here

bad_spectra_LOPEX<-c("X0176","X0177","X0178","X0179","X0180",
               "X0196","X0197","X0198","X0199","X0200",
               "X0321","X0322","X0323","X0324","X0325")
LOPEX_ref<-LOPEX_ref[-which(names(LOPEX_ref) %in% bad_spectra_LOPEX)]

LOPEX_ref<-match_sensors(LOPEX_ref,splice_at = 862,fixed_sensor = 2,interpolate_wvl = 5)
## need to do smoothing here

################################################
## ANGERS

ANGERS_ref<-t(read.csv("IndependentValidationData/ANGERS/angers2003reflectance.csv",row.names = "Sample_."))
ANGERS_ref<-spectra(ANGERS_ref,
                   bands=colnames(ANGERS_ref),
                   names=rownames(ANGERS_ref))

ANGERS_traits<-read.csv("IndependentValidationData/ANGERS/angers2003metadata.csv")
ANGERS_traits[which(ANGERS_traits== -999,arr.ind=T)]<-NA
meta(ANGERS_ref)$EWT<-ANGERS_traits$Equivalent.Water.Thickness..g.cm2.
meta(ANGERS_ref)$LMA<-ANGERS_traits$Leaf.mass.per.area..g.cm2.*10
meta(ANGERS_ref)$dataset<-"ANGERS"

bad_spectra_ANGERS<-c("X0178","X0179","X0184","X0185","X0196",
               "X0197","X0241","X0250","X0254","X0257",
               "X0258","X0269")
ANGERS_ref<-ANGERS_ref[-which(names(ANGERS_ref) %in% bad_spectra_ANGERS)]

## smooth here

#################################################
## combine

LOPEX_ref<-LOPEX_ref[,400:2400]
ANGERS_ref<-ANGERS_ref[,400:2400]
all_val_ref<-combine(LOPEX_ref,ANGERS_ref)

meta(all_val_ref)$Nmass_pred<-predict(Nmass_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_Nmass_CVmodel)[,,1]
meta(all_val_ref)$Cmass_pred<-predict(Cmass_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_Cmass_CVmodel)[,,1]
meta(all_val_ref)$EWT_pred<-predict(EWT_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_EWT_CVmodel)[,,1]
meta(all_val_ref)$LMA_pred<-predict(LMA_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_LMA_CVmodel)[,,1]
meta(all_val_ref)$LDMC_pred<-predict(LDMC_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_LDMC_CVmodel)[,,1]
meta(all_val_ref)$cellulose_pred<-predict(cellulose_mass_CVmodel,newdata = as.matrix(all_val_ref[,400:2400]),ncomp = ncomp_cellulose_mass_CVmodel)[,,1]

ggplot(data=meta(all_val_ref),aes(x=EWT_pred,y=EWT,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured EWT",x="Predicted EWT")+
  coord_cartesian(xlim=c(0,0.085),ylim=c(0,0.085))

ggplot(data=meta(all_val_ref),aes(x=LDMC_pred,y=LDMC*1000,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured LDMC",x="Predicted LDMC")

ggplot(data=meta(all_val_ref),aes(x=LMA_pred,y=LMA*1000,color=dataset))+
  geom_point()+
  theme_bw()+
  theme(text = element_text(size=20))+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_smooth(method="lm")+
  labs(y="Measured LMA",x="Predicted LMA")
