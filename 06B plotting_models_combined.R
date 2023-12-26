# setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(ggplot2)
library(RColorBrewer)
library(patchwork)
source("Scripts/CABO-trait-models/00 useful_functions.R")

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

#######################################
## read data

all.jack.df.list.ref<-readRDS("SavedResults/all_jack_df_list_ref.rds")
all.jack.df.list.trans<-readRDS("SavedResults/all_jack_df_list_trans.rds")
all.jack.df.list.abs<-readRDS("SavedResults/all_jack_df_list_abs.rds")

#######################################
## axis limits for validation plots

all.solubles<-c(all.jack.df.list.ref$sol$Measured,
                all.jack.df.list.ref$sol$pred.mean,
                all.jack.df.list.trans$sol$pred.mean,
                all.jack.df.list.abs$sol$pred.mean)
solubles_upper<-max(all.solubles,na.rm=T)+3
solubles_lower<-min(all.solubles,na.rm=T)-3

all.hemicellulose<-c(all.jack.df.list.ref$hemi$Measured,
                     all.jack.df.list.ref$hemi$pred.mean,
                     all.jack.df.list.trans$hemi$pred.mean,
                     all.jack.df.list.abs$hemi$pred.mean)
hemicellulose_upper<-max(all.hemicellulose,na.rm=T)+3
hemicellulose_lower<-min(all.hemicellulose,na.rm=T)-3

all.cellulose<-c(all.jack.df.list.ref$cell$Measured,
                 all.jack.df.list.ref$cell$pred.mean,
                 all.jack.df.list.trans$cell$pred.mean,
                 all.jack.df.list.abs$cell$pred.mean)
cellulose_upper<-max(all.cellulose,na.rm=T)+2
cellulose_lower<-min(all.cellulose,na.rm=T)-2

all.lignin<-c(all.jack.df.list.ref$lign$Measured,
              all.jack.df.list.ref$lign$pred.mean,
              all.jack.df.list.trans$lign$pred.mean,
              all.jack.df.list.abs$lign$pred.mean)
lignin_upper<-max(all.lignin,na.rm=T)+2
lignin_lower<-min(all.lignin,na.rm=T)-2

all.perN<-c(all.jack.df.list.ref$N$Measured,
            all.jack.df.list.ref$N$pred.mean,
            all.jack.df.list.trans$N$pred.mean,
            all.jack.df.list.abs$N$pred.mean)
perN_upper<-max(all.perN,na.rm=T)+0.2
perN_lower<-min(all.perN,na.rm=T)-0.2

all.perC<-c(all.jack.df.list.ref$C$Measured,
            all.jack.df.list.ref$C$pred.mean,
            all.jack.df.list.trans$C$pred.mean,
            all.jack.df.list.abs$C$pred.mean)
perC_upper<-max(all.perC,na.rm=T)+2
perC_lower<-min(all.perC,na.rm=T)-2

all.LMA<-c(all.jack.df.list.ref$LMA$Measured,
           all.jack.df.list.ref$LMA$pred.mean,
           all.jack.df.list.trans$LMA$pred.mean,
           all.jack.df.list.abs$LMA$pred.mean)
LMA_upper<-max(all.LMA,na.rm=T)+0.02
LMA_lower<-min(all.LMA,na.rm=T)-0.02

all.LDMC<-c(all.jack.df.list.ref$LDMC$Measured,
            all.jack.df.list.ref$LDMC$pred.mean,
            all.jack.df.list.trans$LDMC$pred.mean,
            all.jack.df.list.abs$LDMC$pred.mean)
LDMC_upper<-max(all.LDMC,na.rm=T)+30
LDMC_lower<-min(all.LDMC,na.rm=T)-30

all.EWT<-c(all.jack.df.list.ref$EWT$Measured,
           all.jack.df.list.ref$EWT$pred.mean,
           all.jack.df.list.trans$EWT$pred.mean,
           all.jack.df.list.abs$EWT$pred.mean)
EWT_upper<-max(all.EWT,na.rm=T)+0.02
EWT_lower<-min(all.EWT,na.rm=T)-0.02

all.chlA<-c(all.jack.df.list.ref$chlA$Measured,
            all.jack.df.list.ref$chlA$pred.mean,
            all.jack.df.list.trans$chlA$pred.mean,
            all.jack.df.list.abs$chlA$pred.mean)
chlA_upper<-max(all.chlA,na.rm=T)+1
chlA_lower<-min(all.chlA,na.rm=T)-1

all.chlB<-c(all.jack.df.list.ref$chlB$Measured,
            all.jack.df.list.ref$chlB$pred.mean,
            all.jack.df.list.trans$chlB$pred.mean,
            all.jack.df.list.abs$chlB$pred.mean)
chlB_upper<-max(all.chlB,na.rm=T)+0.4
chlB_lower<-min(all.chlB,na.rm=T)-0.4

all.car<-c(all.jack.df.list.ref$car$Measured,
           all.jack.df.list.ref$car$pred.mean,
           all.jack.df.list.trans$car$pred.mean,
           all.jack.df.list.abs$car$pred.mean)
car_upper<-max(all.car,na.rm=T)+0.2
car_lower<-min(all.car,na.rm=T)-0.2

all.Al<-c(all.jack.df.list.ref$Al$Measured,
          all.jack.df.list.ref$Al$pred.mean[!is.na(all.jack.df.list.ref$Al$Measured)],
          all.jack.df.list.trans$Al$pred.mean[!is.na(all.jack.df.list.ref$Al$Measured)],
          all.jack.df.list.abs$Al$pred.mean[!is.na(all.jack.df.list.ref$Al$Measured)])
Al_upper<-max(all.Al,na.rm=T)+0.02
Al_lower<-min(all.Al,na.rm=T)-0.02

all.Ca<-c(all.jack.df.list.ref$Ca$Measured,
          all.jack.df.list.ref$Ca$pred.mean[!is.na(all.jack.df.list.ref$Ca$Measured)],
          all.jack.df.list.trans$Ca$pred.mean[!is.na(all.jack.df.list.ref$Ca$Measured)],
          all.jack.df.list.abs$Ca$pred.mean[!is.na(all.jack.df.list.ref$Ca$Measured)])
Ca_upper<-max(all.Ca,na.rm=T)+2
Ca_lower<-min(all.Ca,na.rm=T)-2

all.Cu<-c(all.jack.df.list.ref$Cu$Measured,
          all.jack.df.list.ref$Cu$pred.mean[!is.na(all.jack.df.list.ref$Cu$Measured)],
          all.jack.df.list.trans$Cu$pred.mean[!is.na(all.jack.df.list.ref$Cu$Measured)],
          all.jack.df.list.abs$Cu$pred.mean[!is.na(all.jack.df.list.ref$Cu$Measured)])
Cu_upper<-max(all.Cu,na.rm=T)+0.003
Cu_lower<-min(all.Cu,na.rm=T)-0.003

all.Fe<-c(all.jack.df.list.ref$Fe$Measured,
          all.jack.df.list.ref$Fe$pred.mean[!is.na(all.jack.df.list.ref$Fe$Measured)],
          all.jack.df.list.trans$Fe$pred.mean[!is.na(all.jack.df.list.ref$Fe$Measured)],
          all.jack.df.list.abs$Fe$pred.mean[!is.na(all.jack.df.list.ref$Fe$Measured)])
Fe_upper<-max(all.Fe,na.rm=T)+0.03
Fe_lower<-min(all.Fe,na.rm=T)-0.03

all.K<-c(all.jack.df.list.ref$K$Measured,
         all.jack.df.list.ref$K$pred.mean[!is.na(all.jack.df.list.ref$K$Measured)],
         all.jack.df.list.trans$K$pred.mean[!is.na(all.jack.df.list.ref$K$Measured)],
         all.jack.df.list.abs$K$pred.mean[!is.na(all.jack.df.list.ref$K$Measured)])
K_upper<-max(all.K,na.rm=T)+3
K_lower<-min(all.K,na.rm=T)-3

all.Mg<-c(all.jack.df.list.ref$Mg$Measured,
          all.jack.df.list.ref$Mg$pred.mean[!is.na(all.jack.df.list.ref$Mg$Measured)],
          all.jack.df.list.trans$Mg$pred.mean[!is.na(all.jack.df.list.ref$Mg$Measured)],
          all.jack.df.list.abs$Mg$pred.mean[!is.na(all.jack.df.list.ref$Mg$Measured)])
Mg_upper<-max(all.Mg,na.rm=T)+0.5
Mg_lower<-min(all.Mg,na.rm=T)-0.5

all.Mn<-c(all.jack.df.list.ref$Mn$Measured,
          all.jack.df.list.ref$Mn$pred.mean[!is.na(all.jack.df.list.ref$Mn$Measured)],
          all.jack.df.list.trans$Mn$pred.mean[!is.na(all.jack.df.list.ref$Mn$Measured)],
          all.jack.df.list.abs$Mn$pred.mean[!is.na(all.jack.df.list.ref$Mn$Measured)])
Mn_upper<-max(all.Mn,na.rm=T)+0.1
Mn_lower<-min(all.Mn,na.rm=T)-0.1

all.Na<-c(all.jack.df.list.ref$Na$Measured,
          all.jack.df.list.ref$Na$pred.mean[!is.na(all.jack.df.list.ref$Na$Measured)],
          all.jack.df.list.trans$Na$pred.mean[!is.na(all.jack.df.list.ref$Na$Measured)],
          all.jack.df.list.abs$Na$pred.mean[!is.na(all.jack.df.list.ref$Na$Measured)])
Na_upper<-max(all.Na,na.rm=T)+0.5
Na_lower<-min(all.Na,na.rm=T)-0.5

all.P<-c(all.jack.df.list.ref$P$Measured,
         all.jack.df.list.ref$P$pred.mean[!is.na(all.jack.df.list.ref$P$Measured)],
         all.jack.df.list.trans$P$pred.mean[!is.na(all.jack.df.list.ref$P$Measured)],
         all.jack.df.list.abs$P$pred.mean[!is.na(all.jack.df.list.ref$P$Measured)])
P_upper<-max(all.P,na.rm=T)+0.5
P_lower<-min(all.P,na.rm=T)-0.5

all.Zn<-c(all.jack.df.list.ref$Zn$Measured,
          all.jack.df.list.ref$Zn$pred.mean[!is.na(all.jack.df.list.ref$Zn$Measured)],
          all.jack.df.list.trans$Zn$pred.mean[!is.na(all.jack.df.list.ref$Zn$Measured)],
          all.jack.df.list.abs$Zn$pred.mean[!is.na(all.jack.df.list.ref$Zn$Measured)])
Zn_upper<-max(all.Zn,na.rm=T)+0.05
Zn_lower<-min(all.Zn,na.rm=T)-0.05

######################################################
## reflectance plots

solubles_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$sol,
                                   aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),ylim=c(solubles_lower,solubles_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)+
  ggtitle("Reflectance")+
  scale_color_manual(values=colorBlind)

hemicellulose_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$hemi,
                                        aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),ylim=c(hemicellulose_lower,hemicellulose_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

cellulose_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$cell,
                                    aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),ylim=c(cellulose_lower,cellulose_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

lignin_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$lign,
                                 aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),ylim=c(lignin_lower,lignin_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$chlA,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Reflectance")+
  scale_color_manual(values=colorBlind)

chlB_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$chlB,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),ylim=c(chlB_lower,chlB_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

car_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$car,
                              aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_lower,car_upper),ylim=c(car_lower,car_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"),
       color="Functional group")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cmass.ref.val.plot<-ggplot(all.jack.df.list.ref$C,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perC_lower,perC_upper),ylim=c(perC_lower,perC_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Nmass.ref.val.plot<-ggplot(all.jack.df.list.ref$N,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)+
  ggtitle("Reflectance")+
  scale_color_manual(values=colorBlind)

EWT.ref.val.plot<-ggplot(all.jack.df.list.ref$EWT,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

LMA.ref.val.plot<-ggplot(all.jack.df.list.ref$LMA,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)+
  ggtitle("Reflectance")+
  scale_color_manual(values=colorBlind)

LDMC.ref.val.plot<-ggplot(all.jack.df.list.ref$LDMC,
                          aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Al_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$Al,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Al_lower,Al_upper),ylim=c(Al_lower,Al_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Al (mg g"^-1*")"),x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Reflectance")+
  scale_color_manual(values=colorBlind)

Ca_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$Ca,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),ylim=c(Ca_lower,Ca_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Ca (mg g"^-1*")"),
       x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cu_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$Cu,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),ylim=c(Cu_lower,Cu_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Cu (mg g"^-1*")"),
       x=expression("Predicted Cu (mg g"^-1*")"),
       color="Functional group")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Fe_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$Fe,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),ylim=c(Fe_lower,Fe_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Fe (mg g"^-1*")"),
       x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Reflectance")+
  scale_color_manual(values=colorBlind)

K_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$K,
                            aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mg_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$Mg,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),ylim=c(Mg_lower,Mg_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Mg (mg g"^-1*")"),
       x=expression("Predicted Mg (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mn_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$Mn,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Mn (mg g"^-1*")"),
       x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Reflectance")+
  scale_color_manual(values=colorBlind)

Na_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$Na,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Na_lower,Na_upper),ylim=c(Na_lower,Na_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Na (mg g"^-1*")"),
       x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

P_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$P,
                            aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(P_lower,P_upper),ylim=c(P_lower,P_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured P (mg g"^-1*")"),
       x=expression("Predicted P (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Zn_mass.ref.val.plot<-ggplot(all.jack.df.list.ref$Zn,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),ylim=c(Zn_lower,Zn_upper))+
  theme(text = element_text(size=25))+
  labs(y=expression("Measured Zn (mg g"^-1*")"),
       x=expression("Predicted Zn (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

#########################################
## transmittance plots

solubles_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$sol,
                                   aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),ylim=c(solubles_lower,solubles_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)+
  ggtitle("Transmittance")+
  scale_color_manual(values=colorBlind)

hemicellulose_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$hemi,
                                        aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),ylim=c(hemicellulose_lower,hemicellulose_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

cellulose_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$cell,
                                    aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),ylim=c(cellulose_lower,cellulose_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

lignin_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$lign,
                                 aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),ylim=c(lignin_lower,lignin_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$chlA,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Transmittance")+
  scale_color_manual(values=colorBlind)

chlB_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$chlB,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),ylim=c(chlB_lower,chlB_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

car_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$car,
                              aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_lower,car_upper),ylim=c(car_lower,car_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"),
       color="Functional group")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cmass.trans.val.plot<-ggplot(all.jack.df.list.trans$C,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perC_lower,perC_upper),ylim=c(perC_lower,perC_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Nmass.trans.val.plot<-ggplot(all.jack.df.list.trans$N,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)+
  ggtitle("Transmittance")+
  scale_color_manual(values=colorBlind)

EWT.trans.val.plot<-ggplot(all.jack.df.list.trans$EWT,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

LMA.trans.val.plot<-ggplot(all.jack.df.list.trans$LMA,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)+
  ggtitle("Transmittance")+
  scale_color_manual(values=colorBlind)

LDMC.trans.val.plot<-ggplot(all.jack.df.list.trans$LDMC,
                          aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Al_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$Al,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Al_lower,Al_upper),ylim=c(Al_lower,Al_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Al (mg g"^-1*")"),x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Transmittance")+
  scale_color_manual(values=colorBlind)

Ca_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$Ca,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),ylim=c(Ca_lower,Ca_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Ca (mg g"^-1*")"),
       x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cu_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$Cu,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),ylim=c(Cu_lower,Cu_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Cu (mg g"^-1*")"),
       x=expression("Predicted Cu (mg g"^-1*")"),
       color="Functional group")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Fe_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$Fe,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),ylim=c(Fe_lower,Fe_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Fe (mg g"^-1*")"),
       x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Transmittance")+
  scale_color_manual(values=colorBlind)

K_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$K,
                            aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mg_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$Mg,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),ylim=c(Mg_lower,Mg_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Mg (mg g"^-1*")"),
       x=expression("Predicted Mg (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mn_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$Mn,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Mn (mg g"^-1*")"),
       x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Transmittance")+
  scale_color_manual(values=colorBlind)

Na_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$Na,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Na_lower,Na_upper),ylim=c(Na_lower,Na_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Na (mg g"^-1*")"),
       x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

P_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$P,
                            aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(P_lower,P_upper),ylim=c(P_lower,P_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured P (mg g"^-1*")"),
       x=expression("Predicted P (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Zn_mass.trans.val.plot<-ggplot(all.jack.df.list.trans$Zn,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),ylim=c(Zn_lower,Zn_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Zn (mg g"^-1*")"),
       x=expression("Predicted Zn (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

# pdf("Images/val_plots_trans1.pdf",width = 16,height = 20)
# (EWT.val.plot + LDMC.val.plot + LMA.val.plot) / 
#   (cellulose_mass.val.plot + solubles_mass.val.plot + Cmass.val.plot) / 
#   (hemicellulose_mass.val.plot + Nmass.val.plot + chlB_mass.val.plot) / 
#   (chlA_mass.val.plot + lignin_mass.val.plot + car_mass.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()
# 
# pdf("Images/val_plots_trans2.pdf",width = 16,height = 20)
# (Ca_mass.val.plot + K_mass.val.plot + Zn_mass.val.plot) / 
#   (P_mass.val.plot + Mg_mass.val.plot + Na_mass.val.plot) / 
#   (Fe_mass.val.plot + Mn_mass.val.plot + Al_mass.val.plot) / 
#   (Cu_mass.val.plot + guide_area() + guide_area()) &
#   plot_layout(guides="collect")
# dev.off()

############################################
## absorptance plots

solubles_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$sol,
                                   aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),ylim=c(solubles_lower,solubles_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)+
  ggtitle("Absorptance")+
  scale_color_manual(values=colorBlind)

hemicellulose_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$hemi,
                                        aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),ylim=c(hemicellulose_lower,hemicellulose_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured hemicellulose (%)",x="Predicted hemicellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

cellulose_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$cell,
                                    aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),ylim=c(cellulose_lower,cellulose_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

lignin_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$lign,
                                 aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),ylim=c(lignin_lower,lignin_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)",
       color="Functional group")+
  scale_color_manual(values=colorBlind)

chlA_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$chlA,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Absorptance")+
  scale_color_manual(values=colorBlind)

chlB_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$chlB,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),ylim=c(chlB_lower,chlB_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

car_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$car,
                              aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_lower,car_upper),ylim=c(car_lower,car_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

Cmass.abs.val.plot<-ggplot(all.jack.df.list.abs$C,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perC_lower,perC_upper),ylim=c(perC_lower,perC_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Nmass.abs.val.plot<-ggplot(all.jack.df.list.abs$N,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(perN_lower,perN_upper),ylim=c(perN_lower,perN_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)+
  ggtitle("Absorptance")+
  scale_color_manual(values=colorBlind)

EWT.abs.val.plot<-ggplot(all.jack.df.list.abs$EWT,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=25),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)",
       color="Functional group")+
  scale_color_manual(values=colorBlind)

LMA.abs.val.plot<-ggplot(all.jack.df.list.abs$LMA,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)+
  ggtitle("Absorptance")+
  scale_color_manual(values=colorBlind)

LDMC.abs.val.plot<-ggplot(all.jack.df.list.abs$LDMC,
                          aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Al_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$Al,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Al_lower,Al_upper),ylim=c(Al_lower,Al_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Al (mg g"^-1*")"),
       x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Absorptance")+
  scale_color_manual(values=colorBlind)

Ca_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$Ca,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),ylim=c(Ca_lower,Ca_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Ca (mg g"^-1*")"),
       x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cu_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$Cu,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),ylim=c(Cu_lower,Cu_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Cu (mg g"^-1*")"),
       x=expression("Predicted Cu (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

Fe_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$Fe,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),ylim=c(Fe_lower,Fe_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Fe (mg g"^-1*")"),
       x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Absorptance")+
  scale_color_manual(values=colorBlind)

K_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$K,
                            aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mg_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$Mg,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),ylim=c(Mg_lower,Mg_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Mg (mg g"^-1*")"),
       x=expression("Predicted Mg (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

Mn_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$Mn,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Mn (mg g"^-1*")"),
       x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)+
  ggtitle("Absorptance")+
  scale_color_manual(values=colorBlind)

Na_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$Na,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Na_lower,Na_upper),ylim=c(Na_lower,Na_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Na (mg g"^-1*")"),
       x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

P_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$P,
                            aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(P_lower,P_upper),ylim=c(P_lower,P_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured P (mg g"^-1*")"),
       x=expression("Predicted P (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

Zn_mass.abs.val.plot<-ggplot(all.jack.df.list.abs$Zn,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),ylim=c(Zn_lower,Zn_upper))+
  theme(text = element_text(size=25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Zn (mg g"^-1*")"),
       x=expression("Predicted Zn (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

# pdf("Images/val_plots_abs1.pdf",width = 16,height = 20)
# (EWT.val.plot + LDMC.val.plot + LMA.val.plot) / 
#   (cellulose_mass.val.plot + solubles_mass.val.plot + Cmass.val.plot) / 
#   (hemicellulose_mass.val.plot + Nmass.val.plot + chlB_mass.val.plot) / 
#   (chlA_mass.val.plot + lignin_mass.val.plot + car_mass.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()
# 
# pdf("Images/val_plots_abs2.pdf",width = 16,height = 20)
# (Ca_mass.val.plot + K_mass.val.plot + Zn_mass.val.plot) / 
#   (P_mass.val.plot + Mg_mass.val.plot + Na_mass.val.plot) / 
#   (Fe_mass.val.plot + Mn_mass.val.plot + Al_mass.val.plot) / 
#   (Cu_mass.val.plot + guide_area() + guide_area()) &
#   plot_layout(guides="collect")
# dev.off()

###########################################
## combined plots

# pdf("Images/val_plots_comp1.pdf",width=16,height=17)
# (LMA.ref.val.plot + LMA.trans.val.plot + LMA.abs.val.plot) / 
#   (LDMC.ref.val.plot + LDMC.trans.val.plot + LDMC.abs.val.plot) / 
#   (EWT.ref.val.plot + EWT.trans.val.plot + EWT.abs.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()

# pdf("Images/val_plots_comp2.pdf",width=16,height=17)
# (Nmass.ref.val.plot + Nmass.trans.val.plot + Nmass.abs.val.plot) / 
#   (Cmass.ref.val.plot + Cmass.trans.val.plot + Cmass.abs.val.plot) / 
#   (P_mass.ref.val.plot + P_mass.trans.val.plot + P_mass.abs.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()

# pdf("Images/val_plots_comp3.pdf",width=16,height=21)
# (solubles_mass.ref.val.plot + solubles_mass.trans.val.plot + solubles_mass.abs.val.plot) / 
#   (hemicellulose_mass.ref.val.plot + hemicellulose_mass.trans.val.plot + hemicellulose_mass.abs.val.plot) / 
#   (cellulose_mass.ref.val.plot + cellulose_mass.trans.val.plot + cellulose_mass.abs.val.plot) / 
#   (lignin_mass.ref.val.plot + lignin_mass.trans.val.plot + lignin_mass.abs.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()

# pdf("Images/val_plots_comp4.pdf",width=16,height=17)
# (chlA_mass.ref.val.plot + chlA_mass.trans.val.plot + chlA_mass.abs.val.plot) / 
#   (chlB_mass.ref.val.plot + chlB_mass.trans.val.plot + chlB_mass.abs.val.plot) / 
#   (car_mass.ref.val.plot + car_mass.trans.val.plot + car_mass.abs.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()

# pdf("Images/val_plots_comp5.pdf",width=16,height=17)
# (Al_mass.ref.val.plot + Al_mass.trans.val.plot + Al_mass.abs.val.plot) / 
#   (Ca_mass.ref.val.plot + Ca_mass.trans.val.plot + Ca_mass.abs.val.plot) / 
#   (Cu_mass.ref.val.plot + Cu_mass.trans.val.plot + Cu_mass.abs.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()

# pdf("Images/val_plots_comp6.pdf",width=16,height=17)
# (Fe_mass.ref.val.plot + Fe_mass.trans.val.plot + Fe_mass.abs.val.plot) / 
#   (K_mass.ref.val.plot + K_mass.trans.val.plot + K_mass.abs.val.plot) / 
#   (Mg_mass.ref.val.plot + Mg_mass.trans.val.plot + Mg_mass.abs.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()

# pdf("Images/val_plots_comp7.pdf",width=16,height=17)
# (Mn_mass.ref.val.plot + Mn_mass.trans.val.plot + Mn_mass.abs.val.plot) / 
#   (Na_mass.ref.val.plot + Na_mass.trans.val.plot + Na_mass.abs.val.plot) / 
#   (Zn_mass.ref.val.plot + Zn_mass.trans.val.plot + Zn_mass.abs.val.plot) &
#   plot_layout(guides="collect") & theme(legend.position = "bottom")
# dev.off()

#########################################
## summary statistics

all.jack.df.list.trans<-readRDS("SavedResults/all_jack_df_list_trans.rds")

## validation R2, RMSE, and %RMSE
summ<-data.frame(ncomp=unlist(lapply(all.jack.df.list.trans,function(x) x$ncomp[1])),
                 r2=round(unlist(lapply(all.jack.df.list.trans,function(x) summary(lm(Measured~pred.mean,data=x))$r.squared)),3),
                 rmse=signif(unlist(lapply(all.jack.df.list.trans,function(x) RMSD(x$Measured,x$pred.mean))),3),
                 perrmse=signif(unlist(lapply(all.jack.df.list.trans,function(x) percentRMSD(x$Measured,x$pred.mean,0.025,0.975)))*100,3))
# write.csv(summ,"SavedResults/plsr_summ.csv")

# RMSD(all.jack.df.list.ref.area$N$Measured,all.jack.df.list.ref$N$pred.mean*all.jack.df.list.ref$LMA$pred.mean/1000)*1000
# percentRMSD(all.jack.df.list.ref.area$N$Measured*1000,all.jack.df.list.ref$N$pred.mean*all.jack.df.list.ref$LMA$pred.mean,0.025,0.975)*100
# summary(lm(all.jack.df.list.ref.area$N$Measured*1000~I(all.jack.df.list.ref$N$pred.mean*all.jack.df.list.ref$LMA$pred.mean)))
