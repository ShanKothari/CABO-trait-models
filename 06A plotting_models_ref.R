setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
source("Scripts/CABO-trait-models/00 useful_functions.R")

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

all.jack.df.list.ref<-readRDS("SavedResults/all_jack_df_list_ref.rds")
all.jack.df.list.ref.area<-readRDS("SavedResults/all_jack_df_list_ref_area.rds")

##################################
## axis limits for mass-based models

solubles_all<-with(all.jack.df.list.ref$sol,c(pred.low,pred.high,Measured))
solubles_upper<-max(solubles_all,na.rm=T)+1
solubles_lower<-min(solubles_all,na.rm=T)-1

hemicellulose_all<-with(all.jack.df.list.ref$hemi,c(pred.low,pred.high,Measured))
hemicellulose_upper<-max(hemicellulose_all,na.rm=T)+1
hemicellulose_lower<-min(hemicellulose_all,na.rm=T)-1

cellulose_all<-with(all.jack.df.list.ref$cell,c(pred.low,pred.high,Measured))
cellulose_upper<-max(cellulose_all,na.rm=T)+1
cellulose_lower<-min(cellulose_all,na.rm=T)-1

lignin_all<-with(all.jack.df.list.ref$lign,c(pred.low,pred.high,Measured))
lignin_upper<-max(lignin_all,na.rm=T)+0.5
lignin_lower<-min(lignin_all,na.rm=T)-0.5

chlA_all<-with(all.jack.df.list.ref$chlA,c(pred.low,pred.high,Measured))
chlA_upper<-max(chlA_all,na.rm=T)+1
chlA_lower<-min(chlA_all,na.rm=T)-1

chlB_all<-with(all.jack.df.list.ref$chlB,c(pred.low,pred.high,Measured))
chlB_upper<-max(chlB_all,na.rm=T)+0.5
chlB_lower<-min(chlB_all,na.rm=T)-0.5

car_all<-with(all.jack.df.list.ref$car,c(pred.low,pred.high,Measured))
car_upper<-max(car_all,na.rm=T)+0.3
car_lower<-min(car_all,na.rm=T)-0.3

C_all<-with(all.jack.df.list.ref$C,c(pred.low,pred.high,Measured))
C_upper<-max(C_all,na.rm=T)+1
C_lower<-min(C_all,na.rm=T)-1

N_all<-with(all.jack.df.list.ref$N,c(pred.low,pred.high,Measured))
N_upper<-max(N_all,na.rm=T)+0.1
N_lower<-min(N_all,na.rm=T)-0.1

EWT_all<-with(all.jack.df.list.ref$EWT,c(pred.low,pred.high,Measured))
EWT_upper<-max(EWT_all,na.rm=T)+0.03
EWT_lower<-min(EWT_all,na.rm=T)-0.03

LMA_all<-with(all.jack.df.list.ref$LMA,c(pred.low,pred.high,Measured))
LMA_upper<-max(LMA_all,na.rm=T)+0.02
LMA_lower<-min(LMA_all,na.rm=T)-0.02

LDMC_all<-with(all.jack.df.list.ref$LDMC,c(pred.low,pred.high,Measured))
LDMC_upper<-max(LDMC_all,na.rm=T)+10
LDMC_lower<-min(LDMC_all,na.rm=T)-10

Al_all<-with(all.jack.df.list.ref$Al,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Al_upper<-max(Al_all,na.rm=T)+0.02
Al_lower<-min(Al_all,na.rm=T)-0.02

Ca_all<-with(all.jack.df.list.ref$Ca,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Ca_upper<-max(Ca_all,na.rm=T)+1
Ca_lower<-min(Ca_all,na.rm=T)-1

Cu_all<-with(all.jack.df.list.ref$Cu,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Cu_upper<-max(Cu_all,na.rm=T)+0.003
Cu_lower<-min(Cu_all,na.rm=T)-0.003

Fe_all<-with(all.jack.df.list.ref$Fe,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Fe_upper<-max(Fe_all,na.rm=T)+0.01
Fe_lower<-min(Fe_all,na.rm=T)-0.01

K_all<-with(all.jack.df.list.ref$K,c(pred.low[!is.na(Measured)],
                             pred.high[!is.na(Measured)],
                             Measured))
K_upper<-max(K_all,na.rm=T)+1
K_lower<-min(K_all,na.rm=T)-1

Mg_all<-with(all.jack.df.list.ref$Mg,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Mg_upper<-max(Mg_all,na.rm=T)+0.5
Mg_lower<-min(Mg_all,na.rm=T)-0.5

Mn_all<-with(all.jack.df.list.ref$Mn,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Mn_upper<-max(Mn_all,na.rm=T)+0.05
Mn_lower<-min(Mn_all,na.rm=T)-0.05

Na_all<-with(all.jack.df.list.ref$Na,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Na_upper<-max(Na_all,na.rm=T)+0.2
Na_lower<-min(Na_all,na.rm=T)-0.2

P_all<-with(all.jack.df.list.ref$P,c(pred.low[!is.na(Measured)],
                             pred.high[!is.na(Measured)],
                             Measured))
P_upper<-max(P_all,na.rm=T)+0.5
P_lower<-min(P_all,na.rm=T)-0.5

Zn_all<-with(all.jack.df.list.ref$Zn,c(pred.low[!is.na(Measured)],
                               pred.high[!is.na(Measured)],
                               Measured))
Zn_upper<-max(Zn_all,na.rm=T)+0.03
Zn_lower<-min(Zn_all,na.rm=T)-0.03

###########################################
## plotting mass-based models

solubles_mass.val.plot<-ggplot(all.jack.df.list.ref$sol,
                               aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(solubles_lower,solubles_upper),
                  ylim=c(solubles_lower,solubles_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured solubles (%)",x="Predicted solubles (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

hemicellulose_mass.val.plot<-ggplot(all.jack.df.list.ref$hemi,
                                    aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(hemicellulose_lower,hemicellulose_upper),
                  ylim=c(hemicellulose_lower,hemicellulose_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured hemicellulose (%)",
       x="Predicted hemicellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

cellulose_mass.val.plot<-ggplot(all.jack.df.list.ref$cell,
                                aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(cellulose_lower,cellulose_upper),
                  ylim=c(cellulose_lower,cellulose_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured cellulose (%)",x="Predicted cellulose (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

lignin_mass.val.plot<-ggplot(all.jack.df.list.ref$lign,
                             aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(lignin_lower,lignin_upper),
                  ylim=c(lignin_lower,lignin_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured lignin (%)",x="Predicted lignin (%)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_mass.val.plot<-ggplot(all.jack.df.list.ref$chlA,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),
                  ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlB_mass.val.plot<-ggplot(all.jack.df.list.ref$chlB,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlB_lower,chlB_upper),
                  ylim=c(chlB_lower,chlB_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured Chl"~italic("b")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("b")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

car_mass.val.plot<-ggplot(all.jack.df.list.ref$car,
                          aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(car_lower,car_upper),
                  ylim=c(car_lower,car_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured carotenoids (mg g"^-1*")"),
       x=expression("Predicted carotenoids (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

Cmass.val.plot<-ggplot(all.jack.df.list.ref$C,
                       aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(C_lower,C_upper),
                  ylim=c(C_lower,C_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured C (%)"),
       x=expression("Predicted C (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Nmass.val.plot<-ggplot(all.jack.df.list.ref$N,
                       aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(N_lower,N_upper),
                  ylim=c(N_lower,N_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

EWT.val.plot<-ggplot(all.jack.df.list.ref$EWT,
                     aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(EWT_lower,EWT_upper),
                  ylim=c(EWT_lower,EWT_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y="Measured EWT (mm)",x="Predicted EWT (mm)")+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

LMA.val.plot<-ggplot(all.jack.df.list.ref$LMA,
                     aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),
                  ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

LDMC.val.plot<-ggplot(all.jack.df.list.ref$LDMC,
                      aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LDMC_lower,LDMC_upper),
                  ylim=c(LDMC_lower,LDMC_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured LDMC (mg g"^-1*")"),
       x=expression("Predicted LDMC (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Al_mass.val.plot<-ggplot(all.jack.df.list.ref$Al,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Al_lower,Al_upper),
                  ylim=c(Al_lower,Al_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Al (mg g"^-1*")"),
       x=expression("Predicted Al (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Ca_mass.val.plot<-ggplot(all.jack.df.list.ref$Ca,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Ca_lower,Ca_upper),
                  ylim=c(Ca_lower,Ca_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Ca (mg g"^-1*")"),
       x=expression("Predicted Ca (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cu_mass.val.plot<-ggplot(all.jack.df.list.ref$Cu,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Cu_lower,Cu_upper),
                  ylim=c(Cu_lower,Cu_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Cu (mg g"^-1*")"),
       x=expression("Predicted Cu (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Fe_mass.val.plot<-ggplot(all.jack.df.list.ref$Fe,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Fe_lower,Fe_upper),
                  ylim=c(Fe_lower,Fe_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Fe (mg g"^-1*")"),
       x=expression("Predicted Fe (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

K_mass.val.plot<-ggplot(all.jack.df.list.ref$K,
                        aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),
                  ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mg_mass.val.plot<-ggplot(all.jack.df.list.ref$Mg,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mg_lower,Mg_upper),
                  ylim=c(Mg_lower,Mg_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Mg (mg g"^-1*")"),
       x=expression("Predicted Mg (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mn_mass.val.plot<-ggplot(all.jack.df.list.ref$Mn,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Mn_lower,Mn_upper),
                  ylim=c(Mn_lower,Mn_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Mn (mg g"^-1*")"),
       x=expression("Predicted Mn (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Na_mass.val.plot<-ggplot(all.jack.df.list.ref$Na,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Na_lower,Na_upper),
                  ylim=c(Na_lower,Na_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Na (mg g"^-1*")"),
       x=expression("Predicted Na (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

P_mass.val.plot<-ggplot(all.jack.df.list.ref$P,
                        aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(P_lower,P_upper),
                  ylim=c(P_lower,P_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured P (mg g"^-1*")"),
       x=expression("Predicted P (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Zn_mass.val.plot<-ggplot(all.jack.df.list.ref$Zn,
                         aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(Zn_lower,Zn_upper),
                  ylim=c(Zn_lower,Zn_upper))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured Zn (mg g"^-1*")"),
       x=expression("Predicted Zn (mg g"^-1*")"),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

pdf("Images/val_plots_ref1.pdf",width = 16,height = 19)
(LMA.val.plot+LDMC.val.plot+EWT.val.plot)/
  (Nmass.val.plot+Cmass.val.plot+solubles_mass.val.plot)/
  (hemicellulose_mass.val.plot+cellulose_mass.val.plot+lignin_mass.val.plot)/
  (chlA_mass.val.plot+chlB_mass.val.plot+car_mass.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Images/val_plots_ref2.pdf",width = 16,height = 19,onefile=F)
ggarrange(Al_mass.val.plot, Ca_mass.val.plot, Cu_mass.val.plot,
          Fe_mass.val.plot,K_mass.val.plot,Mg_mass.val.plot,
          Mn_mass.val.plot, Na_mass.val.plot,P_mass.val.plot,
          Zn_mass.val.plot,ncol=3, nrow=4, common.legend = TRUE, legend="bottom")
dev.off()

#################################################
## plotting area-based models

solubles_area.val.plot<-ggplot(all.jack.df.list.ref.area$sol,
                               aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-3,25),ylim=c(-3,25))+
  theme(text = element_text(size=20),
        legend.position = c(0, 0.25))+
  labs(y=expression("Measured solubles (mg cm"^-2*")"),
       x=expression("Measured solubles (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

hemicellulose_area.val.plot<-ggplot(all.jack.df.list.ref.area$hemi,
                                    aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-1,5.5),ylim=c(-1,5.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured hemicellulose (mg cm"^-2*")"),
       x=expression("Measured hemicellulose (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

cellulose_area.val.plot<-ggplot(all.jack.df.list.ref.area$cell,
                                aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-1,7),ylim=c(-1,7))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured cellulose (mg cm"^-2*")"),
       x=expression("Measured cellulose (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

lignin_area.val.plot<-ggplot(all.jack.df.list.ref.area$lign,
                             aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.5,4.5),ylim=c(-0.5,4.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression("Measured lignin (mg cm"^-2*")"),
       x=expression("Measured lignin (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_area.val.plot<-ggplot(all.jack.df.list.ref.area$chlA,
                           aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,250),ylim=c(0,250))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression(paste("Measured Chl"~italic("a")~"(",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Chl"~italic("a")~"(",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlB_area.val.plot<-ggplot(all.jack.df.list.ref.area$chlB,
                           aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,75),ylim=c(0,75))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression(paste("Measured Chl"~italic("b")~"(",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Chl"~italic("b")~"(",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

car_area.val.plot<-ggplot(all.jack.df.list.ref.area$car,
                          aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,50),ylim=c(0,50))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25))+
  labs(y=expression(paste("Measured carotenoids (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted carotenoids (",mu,"g cm"^-2*")")),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

Carea.val.plot<-ggplot(all.jack.df.list.ref.area$C,
                       aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,19),ylim=c(0,19))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured C (mg cm"^-2*")"),
       x=expression("Predicted C (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Narea.val.plot<-ggplot(all.jack.df.list.ref.area$N,
                       aes(y=Measured*1000,x=pred.mean*1000,color=functional.group))+
  geom_errorbarh(aes(y=Measured*1000,xmin=pred.low*1000,xmax=pred.high*1000),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.42),ylim=c(0,0.42))+
  theme(text = element_text(size=20))+
  labs(y=expression("Measured N (mg cm"^-2*")"),
       x=expression("Predicted N (mg cm"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Al_area.val.plot<-ggplot(all.jack.df.list.ref.area$Al,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.1,1.8),ylim=c(-0.1,1.8))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Al (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Al (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Ca_area.val.plot<-ggplot(all.jack.df.list.ref.area$Ca,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-50,300),ylim=c(-50,300))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Ca (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Ca (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Cu_area.val.plot<-ggplot(all.jack.df.list.ref.area$Cu,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-0.025,0.25),ylim=c(-0.025,0.25))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Cu (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Cu (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Fe_area.val.plot<-ggplot(all.jack.df.list.ref.area$Fe,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,2),ylim=c(0,2))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Fe (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Fe (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

K_area.val.plot<-ggplot(all.jack.df.list.ref.area$K,
                        aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,150),ylim=c(0,150))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured K (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted K (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mg_area.val.plot<-ggplot(all.jack.df.list.ref.area$Mg,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,60),ylim=c(0,60))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Mg (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Mg (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Mn_area.val.plot<-ggplot(all.jack.df.list.ref.area$Mn,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-5,12),ylim=c(-5,12))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Mn (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Mn (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Na_area.val.plot<-ggplot(all.jack.df.list.ref.area$Na,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-10,60),ylim=c(-10,60))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Na (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Na (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

P_area.val.plot<-ggplot(all.jack.df.list.ref.area$P,
                        aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,45),ylim=c(0,45))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured P (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted P (",mu,"g cm"^-2*")")))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Zn_area.val.plot<-ggplot(all.jack.df.list.ref.area$Zn,
                         aes(y=Measured*10^6,x=pred.mean*10^6,color=functional.group))+
  geom_errorbarh(aes(y=Measured*10^6,xmin=pred.low*10^6,xmax=pred.high*10^6),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(-1,2.5),ylim=c(-1,2.5))+
  theme(text = element_text(size=20))+
  labs(y=expression(paste("Measured Zn (",mu,"g cm"^-2*")")),
       x=expression(paste("Predicted Zn (",mu,"g cm"^-2*")")),
       color="Functional group")+
  scale_color_manual(values=colorBlind)

pdf("Images/val_plots_ref_area1.pdf",width = 16,height = 15)
(Narea.val.plot+Carea.val.plot+solubles_area.val.plot)/
  (hemicellulose_area.val.plot+cellulose_area.val.plot+lignin_area.val.plot)/
  (chlA_area.val.plot+chlB_area.val.plot+car_area.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

pdf("Images/val_plots_ref_area2.pdf",width = 16,height = 19,onefile=F)
ggarrange(Al_area.val.plot, Ca_area.val.plot, Cu_area.val.plot,
          Fe_area.val.plot,K_area.val.plot,Mg_area.val.plot,
          Mn_area.val.plot, Na_area.val.plot,P_area.val.plot,
          Zn_area.val.plot,ncol=3, nrow=4, common.legend = TRUE, legend="bottom")
dev.off()

####################################
## comparing trait transformations

all.sqrt.jack.df.list.ref<-readRDS("SavedResults/all_jack_df_list_ref_sqrt.rds")
all.log.jack.df.list.ref<-readRDS("SavedResults/all_jack_df_list_ref_log.rds")

Nmass.sqrt.val.plot<-ggplot(all.sqrt.jack.df.list.ref$N,
                       aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(N_lower,N_upper),
                  ylim=c(N_lower,N_upper))+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

Nmass.log.val.plot<-ggplot(all.log.jack.df.list.ref$N,
                            aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(N_lower,N_upper),
                  ylim=c(N_lower,N_upper))+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured N (%)"),
       x=expression("Predicted N (%)"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_mass.sqrt.val.plot<-ggplot(all.sqrt.jack.df.list.ref$chlA,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),
                  ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

chlA_mass.log.val.plot<-ggplot(all.log.jack.df.list.ref$chlA,
                           aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(chlA_lower,chlA_upper),
                  ylim=c(chlA_lower,chlA_upper))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.25),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured Chl"~italic("a")~"(mg g"^-1*")"),
       x=expression("Predicted Chl"~italic("a")~"(mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

LMA.sqrt.val.plot<-ggplot(all.sqrt.jack.df.list.ref$LMA,
                     aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),
                  ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)+
  ggtitle("Square root-transformed")

LMA.log.val.plot<-ggplot(all.log.jack.df.list.ref$LMA,
                     aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(LMA_lower,LMA_upper),
                  ylim=c(LMA_lower,LMA_upper))+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)+
  ggtitle("Log-transformed")

K_mass.sqrt.val.plot<-ggplot(all.sqrt.jack.df.list.ref$K,
                        aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),
                  ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

K_mass.log.val.plot<-ggplot(all.log.jack.df.list.ref$K,
                        aes(y=Measured,x=pred.mean,color=functional.group))+
  geom_errorbarh(aes(y=Measured,xmin=pred.low,xmax=pred.high),
                 color="gray")+
  geom_point(size=2,alpha=0.7)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(K_lower,K_upper),
                  ylim=c(K_lower,K_upper))+
  theme(text = element_text(size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  labs(y=expression("Measured K (mg g"^-1*")"),
       x=expression("Predicted K (mg g"^-1*")"))+
  guides(color=F)+
  scale_color_manual(values=colorBlind)

pdf("Images/val_plots_comp_transform.pdf",width=16,height=21)
((LMA.val.plot + ggtitle("Untransformed")) + LMA.sqrt.val.plot + LMA.log.val.plot) / 
  (Nmass.val.plot + Nmass.sqrt.val.plot + Nmass.log.val.plot) / 
  (chlA_mass.val.plot + chlA_mass.sqrt.val.plot + chlA_mass.log.val.plot) /
  (K_mass.val.plot + K_mass.sqrt.val.plot + K_mass.log.val.plot) &
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()
