setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
source("Scripts/CABO-trait-models/00 useful_functions.R")

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

all.jack.coefs.list.ref<-readRDS("SavedResults/all_jack_coefs_list_ref.rds")
all.jack.coefs.list.ref.summ<-lapply(all.jack.coefs.list.ref,function(coef.list) {
  coef.mat<-matrix(unlist(coef.list),nrow=length(coef.list),byrow=T)
  return(coef.mat)
  })

coef.plot.df.list<-lapply(all.jack.coefs.list.ref.summ,function(coef.mat){
  coef.means<-colMeans(coef.mat)
  coef.025<-apply(coef.mat,2,quantile,probs=0.025)
  coef.975<-apply(coef.mat,2,quantile,probs=0.975)
  coef.summ.df<-data.frame(mean=coef.means,
                           per.025=coef.025,
                           per.975=coef.975)
  return(coef.summ.df)
})

sol_label<-paste0("Intercept = ",signif(coef.plot.df.list$sol$mean[1],digits=3),
                 " (",signif(coef.plot.df.list$sol$per.025[1],digits=3),
                 ", ",signif(coef.plot.df.list$sol$per.975[1],digits=3),")")
sol_coefs<-ggplot(coef.plot.df.list$sol[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$sol[-1,]),label=sol_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$sol[-1,]*1.05),label="Solubles")

hemi_label<-paste0("Intercept = ",signif(coef.plot.df.list$hemi$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$hemi$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$hemi$per.975[1],digits=3),")")
hemi_coefs<-ggplot(coef.plot.df.list$hemi[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$hemi[-1,]),label=hemi_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$hemi[-1,]*1.05),label="Hemicellulose")

cell_label<-paste0("Intercept = ",signif(coef.plot.df.list$cell$mean[1],digits=3),
                   " (",signif(coef.plot.df.list$cell$per.025[1],digits=3),
                   ", ",signif(coef.plot.df.list$cell$per.975[1],digits=3),")")
cell_coefs<-ggplot(coef.plot.df.list$cell[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$cell[-1,]),label=cell_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$cell[-1,]*1.05),label="Cellulose")

lign_label<-paste0("Intercept = ",signif(coef.plot.df.list$lign$mean[1],digits=3),
                   " (",signif(coef.plot.df.list$lign$per.025[1],digits=3),
                   ", ",signif(coef.plot.df.list$lign$per.975[1],digits=3),")")
lign_coefs<-ggplot(coef.plot.df.list$lign[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$lign[-1,]),label=lign_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$lign[-1,]*1.05),label="Lignin")

LMA_label<-paste0("Intercept = ",signif(coef.plot.df.list$LMA$mean[1],digits=3),
                   " (",signif(coef.plot.df.list$LMA$per.025[1],digits=3),
                   ", ",signif(coef.plot.df.list$LMA$per.975[1],digits=3),")")
LMA_coefs<-ggplot(coef.plot.df.list$LMA[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$LMA[-1,]),label=LMA_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$LMA[-1,]*1.05),label="LMA")

LDMC_label<-paste0("Intercept = ",signif(coef.plot.df.list$LDMC$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$LDMC$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$LDMC$per.975[1],digits=3),")")
LDMC_coefs<-ggplot(coef.plot.df.list$LDMC[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$LDMC[-1,]),label=LDMC_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$LDMC[-1,]*1.05),label="LDMC")

EWT_label<-paste0("Intercept = ",signif(coef.plot.df.list$EWT$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$EWT$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$EWT$per.975[1],digits=3),")")
EWT_coefs<-ggplot(coef.plot.df.list$EWT[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$EWT[-1,]),label=EWT_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$EWT[-1,]*1.05),label="EWT")

N_label<-paste0("Intercept = ",signif(coef.plot.df.list$N$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$N$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$N$per.975[1],digits=3),")")
N_coefs<-ggplot(coef.plot.df.list$N[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$N[-1,]),label=N_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$N[-1,]*1.05),label="%N")

C_label<-paste0("Intercept = ",signif(coef.plot.df.list$C$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$C$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$C$per.975[1],digits=3),")")
C_coefs<-ggplot(coef.plot.df.list$C[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$C[-1,]),label=C_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$C[-1,]*1.05),label="%C")

chlA_label<-paste0("Intercept = ",signif(coef.plot.df.list$chlA$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$chlA$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$chlA$per.975[1],digits=3),")")
chlA_coefs<-ggplot(coef.plot.df.list$chlA[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$chlA[-1,]),label=chlA_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$chlA[-1,]*1.05),label=expression("Chl"~italic("a")))

chlB_label<-paste0("Intercept = ",signif(coef.plot.df.list$chlB$mean[1],digits=3),
                   " (",signif(coef.plot.df.list$chlB$per.025[1],digits=3),
                   ", ",signif(coef.plot.df.list$chlB$per.975[1],digits=3),")")
chlB_coefs<-ggplot(coef.plot.df.list$chlB[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$chlB[-1,]),label=chlB_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$chlB[-1,]*1.05),label=expression("Chl"~italic("b")))

car_label<-paste0("Intercept = ",signif(coef.plot.df.list$car$mean[1],digits=3),
                   " (",signif(coef.plot.df.list$car$per.025[1],digits=3),
                   ", ",signif(coef.plot.df.list$car$per.975[1],digits=3),")")
car_coefs<-ggplot(coef.plot.df.list$car[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$car[-1,]),label=car_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$car[-1,]*1.05),label="Carotenoids")

Al_label<-paste0("Intercept = ",signif(coef.plot.df.list$Al$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$Al$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$Al$per.975[1],digits=3),")")
Al_coefs<-ggplot(coef.plot.df.list$Al[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$Al[-1,]),label=Al_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$Al[-1,]*1.05),label="Al")

Ca_label<-paste0("Intercept = ",signif(coef.plot.df.list$Ca$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$Ca$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$Ca$per.975[1],digits=3),")")
Ca_coefs<-ggplot(coef.plot.df.list$Ca[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$Ca[-1,]),label=Ca_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$Ca[-1,]*1.05),label="Ca")

Cu_label<-paste0("Intercept = ",signif(coef.plot.df.list$Cu$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$Cu$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$Cu$per.975[1],digits=3),")")
Cu_coefs<-ggplot(coef.plot.df.list$Cu[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$Cu[-1,]),label=Cu_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$Cu[-1,]*1.05),label="Cu")

Fe_label<-paste0("Intercept = ",signif(coef.plot.df.list$Fe$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$Fe$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$Fe$per.975[1],digits=3),")")
Fe_coefs<-ggplot(coef.plot.df.list$Fe[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$Fe[-1,]),label=Fe_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$Fe[-1,]*1.05),label="Fe")


K_label<-paste0("Intercept = ",signif(coef.plot.df.list$K$mean[1],digits=3),
                 " (",signif(coef.plot.df.list$K$per.025[1],digits=3),
                 ", ",signif(coef.plot.df.list$K$per.975[1],digits=3),")")
K_coefs<-ggplot(coef.plot.df.list$K[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$K[-1,]),label=K_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$K[-1,]*1.05),label="K")

Mg_label<-paste0("Intercept = ",signif(coef.plot.df.list$Mg$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$Mg$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$Mg$per.975[1],digits=3),")")
Mg_coefs<-ggplot(coef.plot.df.list$Mg[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$Mg[-1,]),label=Mg_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$Mg[-1,]*1.05),label="Mg")

Mn_label<-paste0("Intercept = ",signif(coef.plot.df.list$Mn$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$Mn$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$Mn$per.975[1],digits=3),")")
Mn_coefs<-ggplot(coef.plot.df.list$Mn[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$Mn[-1,]),label=Mn_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$Mn[-1,]*1.05),label="Mn")

Na_label<-paste0("Intercept = ",signif(coef.plot.df.list$Na$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$Na$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$Na$per.975[1],digits=3),")")
Na_coefs<-ggplot(coef.plot.df.list$Na[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$Na[-1,]),label=Na_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$Na[-1,]*1.05),label="Na")

P_label<-paste0("Intercept = ",signif(coef.plot.df.list$P$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$P$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$P$per.975[1],digits=3),")")
P_coefs<-ggplot(coef.plot.df.list$P[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$P[-1,]),label=P_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$P[-1,]*1.05),label="P")

Zn_label<-paste0("Intercept = ",signif(coef.plot.df.list$Zn$mean[1],digits=3),
                  " (",signif(coef.plot.df.list$Zn$per.025[1],digits=3),
                  ", ",signif(coef.plot.df.list$Zn$per.975[1],digits=3),")")
Zn_coefs<-ggplot(coef.plot.df.list$Zn[-1,])+
  geom_ribbon(aes(x=400:2400,
                  ymin = per.025,
                  ymax = per.975),
              alpha = 0.5,fill = "blue")+
  geom_line(aes(x=400:2400,
                y=mean),
            size=1,color="black")+
  theme_bw()+
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"))+
  labs(x="Wavelength (nm)",y="PLSR coefficient")+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2410))+
  annotate(geom="text",hjust=0, x=400,size=6,
           y=min(coef.plot.df.list$Zn[-1,]),label=Zn_label)+
  annotate(geom="text",hjust=0, x=400,size=8,
           y=max(coef.plot.df.list$Zn[-1,]*1.05),label="Zn")

pdf("Images/PLSR_coef_plot1.pdf",width=20,height=12)
plot_grid(plotlist = list(LMA_coefs,EWT_coefs,LDMC_coefs,
                          N_coefs,C_coefs,sol_coefs,
                          hemi_coefs,cell_coefs,lign_coefs,
                          chlA_coefs,chlB_coefs,car_coefs),
          ncol = 3, align = "v")
dev.off()

pdf("Images/PLSR_coef_plot2.pdf",width=20,height=12)
plot_grid(plotlist = list(Al_coefs,Ca_coefs,Cu_coefs,
                          Fe_coefs,K_coefs,Mg_coefs,
                          Mn_coefs,Na_coefs,P_coefs,
                          Zn_coefs,NULL,NULL),
          ncol = 3, align = "v")
dev.off()
