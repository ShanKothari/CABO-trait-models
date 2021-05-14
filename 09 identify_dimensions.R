setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(factoextra)
library(dplyr)
library(hsdar)
library(prospectr)
library(vegan)
library(intrinsicDimension)
library(lle)
library(ider)

ref.traits<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")

## why is there negative reflectance in some of these spectra? from the Dessain data
ref.traits<-ref.traits[-which(apply(as.matrix(ref.traits),1,function(x) sum(x<0))>0)]

#################################################
## factor analysis of spectra

## sample no more than three of the same species
## to avoid over-weighting common species in importance
ref.traits.samp<-subset_by(ref.traits,by=as.vector(meta(ref.traits)$species),
                            n_min = 1,n_max = 3)
ref.traits.ref<-as.data.frame(ref.traits.samp)[,as.character(400:2400)]
colnames(ref.traits.ref)<-paste("X",colnames(ref.traits.ref),sep="")
ref.pca<-prcomp(~.,data=ref.traits.ref,scale.=T)
meta(ref.traits.samp)$PC1<-ref.pca$x[,1]
meta(ref.traits.samp)$PC2<-ref.pca$x[,2]
meta(ref.traits.samp)$PC3<-ref.pca$x[,3]
meta(ref.traits.samp)$PC4<-ref.pca$x[,4]

## average spectrum
# plot(400:2400,ref.pca$center,type="l")
## scree plot
# fviz_eig(ref.pca)

## reconstructing the spectrum
# recon.unscaled<-ref.pca$x %*% t(ref.pca$rotation)
# recon<-sapply(1:ncol(recon.unscaled),function(i) recon.unscaled[,i]*ref.pca$scale[i]+ref.pca$center[i])

## here we rescale all PCA axes to a standard deviation of 1
ref.pca.scale<-ref.pca
ref.pca.scale$x<-sapply(1:ncol(ref.pca$x),function(i) ref.pca$x[,i]/ref.pca$sdev[i])
ref.pca.scale$rotation<-sapply(1:ncol(ref.pca$rotation),function(i) ref.pca$rotation[,i]*ref.pca$sdev[i])

## we can again reconstruct the spectrum with the same procedure
# recon.unscaled2<-ref.pca.scale$x %*% t(ref.pca.scale$rotation)
# recon2<-sapply(1:ncol(recon.unscaled2),function(i) recon.unscaled2[,i]*ref.pca.scale$scale[i]+ref.pca.scale$center[i])

## matrix with the specified PC at -1, 0, or 1 and all other PCs at 0
PC<-1
trial.matrix<-matrix(data=0,nrow=3,ncol=ncol(ref.pca.scale$x))
trial.matrix[,PC]<-c(-1,0,1)
recon.trial.unscaled<-trial.matrix %*% t(ref.pca.scale$rotation)
recon.trial<-sapply(1:ncol(recon.trial.unscaled),function(i) recon.trial.unscaled[,i]*ref.pca.scale$scale[i]+ref.pca.scale$center[i])

plot(400:2400,recon.trial[1,],type="l",col="red",cex.lab=2,
     xlab="Wavelength",ylab="Reflectance",ylim=c(0,0.55))
points(400:2400,recon.trial[2,],type="l",col="black")
points(400:2400,recon.trial[3,],type="l",col="blue")

recon.trial.lib<-speclib(spectra=recon.trial,wavelength=400:2400)
ref.inv<-list()
for(i in 1:dim(recon.trial.lib)[1]){ref.inv[[i]]<-PROSPECTinvert(recon.trial.lib[i])}

###################################################
## continuum removal

ref.traits.cont<-continuumRemoval(ref.traits.ref,wav=400:2400,
                                   interpol="linear",method="substraction")
ref.pca.cont<-prcomp(~.,data=data.frame(ref.traits.cont[,which(apply(ref.traits.cont,2,sd)!=0)]),scale.=T)

plot(spectrolab::spectra(ref.traits.cont,bands=400:2400,names=names(ref.traits.samp)))

################################
## nonlinear dimensionality reduction

ref.traits.norm<-apply(as.matrix(ref.traits.samp),2,function(x) (x-mean(x))/sd(x))
ref.traits.dist<-dist(as.matrix(ref.traits.norm),method = "manhattan")

## isomap
ref.isomap<-isomap(ref.traits.dist,ndim = 10,k=3)

colnames(ref.isomap$points)<-paste("D",1:ncol(ref.isomap$points),sep="")
ref.isomap.coords<-as.data.frame(ref.isomap$points)
ref.isomap.coords$functional.group<-as.factor(meta(ref.traits.samp)$functional.group)
ref.isomap.coords$LMA<-meta(ref.traits.samp)$LMA
ref.isomap.coords$EWT<-meta(ref.traits.samp)$EWT
ref.isomap.coords$Nmass<-meta(ref.traits.samp)$Nmass
ref.isomap.coords$chlA_mass<-meta(ref.traits.samp)$chlA_mass
ref.isomap.coords$R800<-ref.traits.samp[,800]

ggplot(data=ref.isomap.coords,
       aes(x=D1,y=D2,color=EWT))+
  geom_point()+
  scale_color_viridis_c()+
  theme_bw()

ggplot(data=ref.isomap.coords,
       aes(x=D3,y=D4,color=chlA_mass))+
  geom_point()+
  scale_color_viridis_c()+
  theme_bw()

plot(ref.isomap,col=as.numeric(functional.group),cex=2)
legend(150,0,levels(functional.group),col=1:length(functional.group),pch=1)

## NMDS
ref.nmds<-metaMDS(ref.traits.norm,distance="euclidean",k=7,
                    autotransform = F)
plot(ref.isomap,col=as.numeric(functional.group),cex=2)
legend(150,0,levels(functional.group),col=1:length(functional.group),pch=1)


#################################
## estimating intrinsic dimensionality

traits.samp<-meta(ref.traits.samp)[,c("LDMC","EWT","LMA","Cmass","Nmass",
                                      "solubles_mass","hemicellulose_mass",
                                      "cellulose_mass","lignin_mass","chlA_mass",
                                      "chlB_mass","car_mass")]
traits.norm<-apply(traits.samp,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm = T))


maxLikGlobalDimEst(ref.traits.norm,k=3)

maxLikGlobalDimEst(na.omit(traits.norm),k=3)

###########################################
## LLE (also with dimensionality estimate)

ref.lle<-lle(ref.traits.norm,m=7,k=3,id = T)
traits.lle<-lle(na.omit(traits.norm),m=7,k=3,id = T)

ref.lle.coords<-as.data.frame(ref.lle$Y)
colnames(ref.lle.coords)<-paste("D",1:ncol(ref.lle.coords),sep="")
ref.lle.coords$functional.group<-as.factor(meta(ref.traits.samp)$functional.group)
ref.lle.coords$LMA<-meta(ref.traits.samp)$LMA
ref.lle.coords$EWT<-meta(ref.traits.samp)$EWT
ref.lle.coords$R800<-ref.traits.samp[,800]

ggplot(as.data.frame(ref.lle.coords),
       aes(x=D1,y=D2,color=functional.group))+
  geom_point()+
  theme_bw()

ggplot(as.data.frame(ref.lle.coords),
       aes(x=D2,y=D3,color=R800))+
  scale_colour_viridis_c()+
  geom_point()+
  theme_bw()

#########################################
## does nonlinear dimensionality reduction work
## on simulated spectra?

trait.quant<-c(0.025,0.2,0.4,0.6,0.8,0.975)
N<-seq(1,2,by=0.2)
Cw<-sample(meta(ref.traits.samp)$EWT,probs=trait.quant,na.rm=T)
Cab<-quantile((meta(ref.traits.samp)$chlA_area+meta(ref.traits.samp)$chlB_area)*1000000,
              probs=trait.quant,na.rm=T)
Cm<-quantile(meta(ref.traits.samp)$LMA/10,probs=trait.quant,na.rm=T)
param.list<-list(N=N,Cw=Cw,Cab=Cab,Cm=Cm)
prospect.params<-expand.grid(param.list)

prospect.params$Car<-prospect.params$Cab/6.24 ## mean chl:car ratio
prospect.params$Anth<-0
prospect.params$Cbrown<-0

prospect.sim<-t(apply(prospect.params,1,
                      function(x) {
                        spec.sim<-PROSPECT(N=x["N"],Cab=x["Cab"],Car=x["Car"],
                                           Anth=x["Anth"],Cbrown=x["Cbrown"],
                                           Cm=x["Cm"],Cw=x["Cw"])
                        spec<-spec.sim@spectra@spectra_ma
                        return(spec)
                      }))

prospect.sim.norm<-apply(prospect.sim,2,function(x) (x-mean(x))/sd(x))

maxLikGlobalDimEst(prospect.sim,k=10)
maxLikGlobalDimEst(prospect.sim.norm,k=10)
convU(x = prospect.sim)


prospect.lle<-lle(prospect.sim,m=7,k=4,id = T)
prospect.isomap<-isomap(dist(prospect.sim,method="euclidean"),
                        ndim = 10,k=8)
plot(prospect.isomap$points[,2]~prospect.params$Cw)
