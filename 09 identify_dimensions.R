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
convU(x = ref.traits.norm)
convU(x = na.omit(traits.norm))
## Grassberger-Procaccia; unsure about correct choice of p or k1/k2
corint(x=ref.traits.norm,p=20)
corint(x=na.omit(traits.norm),p=20)
## ML; not sure about correct choice of p or k1/k2
lbmle(x = ref.traits.norm, k1=5, k2=10)
lbmle(x = na.omit(traits.norm), k1=5, k2=10)
mada(x = ref.traits.norm, k = NULL, comb = "average", maxDim = 10)
mada(x = na.omit(traits.norm), k = NULL, comb = "average", maxDim = 10)
nni(x = ref.traits.norm, k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = na.omit(traits.norm), k1 = 5, k2 = 30, eps = 0.01, p = NULL)

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
## dimensionality of fresh, pressed, and ground spectra

fresh.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/fresh_spec_EL_agg.rds")
pressed.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/pressed_spec_EL_agg.rds")
ground.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/ground_spec_EL_agg.rds")

fresh.spec.norm<-apply(as.matrix(fresh.spec),2,function(x) (x-mean(x))/sd(x))
pressed.spec.norm<-apply(as.matrix(pressed.spec),2,function(x) (x-mean(x))/sd(x))
ground.spec.norm<-apply(as.matrix(ground.spec),2,function(x) (x-mean(x))/sd(x))

traits.sub<-apply(meta(fresh.spec)[,c("LDMC","EWT","LMA","C","N",
                                      "solubles","hemicellulose",
                                      "cellulose","lignin","chlA",
                                      "chlB","car","Ca","K","P")],2,
                  function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))

maxLikGlobalDimEst(as.matrix(fresh.spec),k=3)
maxLikGlobalDimEst(as.matrix(pressed.spec),k=3)
maxLikGlobalDimEst(as.matrix(ground.spec),k=3)
maxLikGlobalDimEst(na.omit(traits.sub),k=3)

convU(as.matrix(fresh.spec),maxDim=10)
convU(as.matrix(pressed.spec),maxDim=10)
convU(as.matrix(ground.spec),maxDim=10)
convU(na.omit(traits.sub),maxDim=10)

corint(x=as.matrix(fresh.spec),p=20)
corint(x=as.matrix(pressed.spec),p=20)
corint(x=as.matrix(ground.spec),p=20)
corint(x=na.omit(traits.sub),p=20)

## with k1=k2=k, this yields the same result as
## maxLikGlobalDimEst
lbmle(x = as.matrix(fresh.spec), k1=5, k2=10)
lbmle(x = as.matrix(pressed.spec), k1=5, k2=10)
lbmle(x = as.matrix(ground.spec), k1=5, k2=10)

mada(x = as.matrix(fresh.spec), k = NULL, comb = "average", maxDim = 10)
mada(x = as.matrix(pressed.spec), k = NULL, comb = "average", maxDim = 10)
mada(x = as.matrix(ground.spec), k = NULL, comb = "average", maxDim = 10)
mada(x = na.omit(traits.sub), k = NULL, comb = "average", maxDim = 10)

nni(x = as.matrix(fresh.spec), k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = as.matrix(pressed.spec), k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = as.matrix(ground.spec), k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = na.omit(traits.sub), k1 = 5, k2 = 30, eps = 0.01, p = NULL)

pack(x = as.matrix(fresh.spec))
pack(x = as.matrix(pressed.spec))
pack(x = as.matrix(ground.spec))
pack(x = na.omit(traits.sub))

side(x = as.matrix(fresh.spec),maxDim=10)
side(x = as.matrix(pressed.spec),maxDim=10)
side(x = as.matrix(ground.spec),maxDim=10)
side(x = na.omit(traits.sub),maxDim=10)

lle.fresh<-lle(as.matrix(fresh.spec),m=10,k=10,id = T)
lle.pressed<-lle(as.matrix(pressed.spec),m=10,k=10,id = T)
lle.ground<-lle(as.matrix(ground.spec),m=10,k=10,id = T)
lle.traits<-lle(na.omit(traits.sub),m=10,k=10,id = T)

isomap.fresh<-isomap(dist(as.matrix(fresh.spec),method = "euclidean"),
                     ndim = 10,k=4)
isomap.pressed<-isomap(dist(as.matrix(pressed.spec),method = "euclidean"),
                     ndim = 10,k=4)
isomap.ground<-isomap(dist(as.matrix(ground.spec),method = "euclidean"),
                     ndim = 10,k=4)
isomap.traits<-isomap(dist(na.omit(traits.sub),method = "euclidean"),
                     ndim = 10,k=4)

## try stressplot() with MDS in vegan

#########################################
## does nonlinear dimensionality reduction work
## on simulated spectra?

N<-runif(n = 1000,min=1,max=2.2)
Cw<-sample(na.omit(meta(ref.traits.samp)$EWT),
                           size = 1000,replace=T)
Cab<-sample((na.omit(meta(ref.traits.samp)$chlA_area+meta(ref.traits.samp)$chlB_area)*1000000),
                            size=1000,replace=T)
Cm<-sample(na.omit(meta(ref.traits.samp)$LMA/10),
                           size=1000,replace=T)
prospect.params<-data.frame(N=N,Cw=Cw,Cab=Cab,Cm=Cm)
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
convU(x = prospect.sim.norm)
## Grassberger-Procaccia; unsure about correct choice of p or k1/k2
corint(x=prospect.sim,p=20)
corint(x=prospect.sim.norm,p=20)
## ML; not sure about correct choice of p or k1/k2
lbmle(x = prospect.sim, k1=5, k2=10)
lbmle(x = prospect.sim.norm, k1=5, k2=10)
mada(x = prospect.sim, k = NULL, comb = "average", maxDim = 10)
mada(x = prospect.sim.norm, k = NULL, comb = "average", maxDim = 10)
nni(x = prospect.sim, k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = prospect.sim.norm, k1 = 5, k2 = 30, eps = 0.01, p = NULL)

prospect.lle<-lle(prospect.sim,m=7,k=4,id = T)
prospect.isomap<-isomap(dist(prospect.sim.norm,method="euclidean"),
                        ndim = 10,k=3)
plot(prospect.isomap$points[,3]~prospect.params$Cab)
plot(prospect.isomap$points[,1]~prospect.isomap$points[,2])
