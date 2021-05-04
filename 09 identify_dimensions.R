setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(factoextra)
library(dplyr)
library(hsdar)
library(prospectr)

spec.traits<-readRDS("ProcessedSpectra/all_spectra_and_traits.rds")

## why is there negative reflectance in some of these spectra? from the Dessain data
spec.traits<-spec.traits[-which(apply(as.matrix(spec.traits),1,function(x) sum(x<0))>0)]

#################################################
## factor analysis of spectra

## sample no more than three of the same species
## to avoid over-weighting common species in importance
spec.traits.samp<-subset_by(spec.traits,by=as.vector(meta(spec.traits)$species),
                            n_min = 1,n_max = 3)
spec.traits.ref<-as.data.frame(spec.traits.samp)[,as.character(400:2400)]
colnames(spec.traits.ref)<-paste("X",colnames(spec.traits.ref),sep="")
spec.pca<-prcomp(~.,data=spec.traits.ref,scale.=T)
meta(spec.traits.samp)$PC1<-spec.pca$x[,1]
meta(spec.traits.samp)$PC2<-spec.pca$x[,2]
meta(spec.traits.samp)$PC3<-spec.pca$x[,3]
meta(spec.traits.samp)$PC4<-spec.pca$x[,4]

## average spectrum
# plot(400:2400,spec.pca$center,type="l")
## scree plot
# fviz_eig(spec.pca)

## reconstructing the spectrum
# recon.unscaled<-spec.pca$x %*% t(spec.pca$rotation)
# recon<-sapply(1:ncol(recon.unscaled),function(i) recon.unscaled[,i]*spec.pca$scale[i]+spec.pca$center[i])

## here we rescale all PCA axes to a standard deviation of 1
spec.pca.scale<-spec.pca
spec.pca.scale$x<-sapply(1:ncol(spec.pca$x),function(i) spec.pca$x[,i]/spec.pca$sdev[i])
spec.pca.scale$rotation<-sapply(1:ncol(spec.pca$rotation),function(i) spec.pca$rotation[,i]*spec.pca$sdev[i])

## we can again reconstruct the spectrum with the same procedure
# recon.unscaled2<-spec.pca.scale$x %*% t(spec.pca.scale$rotation)
# recon2<-sapply(1:ncol(recon.unscaled2),function(i) recon.unscaled2[,i]*spec.pca.scale$scale[i]+spec.pca.scale$center[i])

## matrix with the specified PC at -1, 0, or 1 and all other PCs at 0
PC<-1
trial.matrix<-matrix(data=0,nrow=3,ncol=ncol(spec.pca.scale$x))
trial.matrix[,PC]<-c(-1,0,1)
recon.trial.unscaled<-trial.matrix %*% t(spec.pca.scale$rotation)
recon.trial<-sapply(1:ncol(recon.trial.unscaled),function(i) recon.trial.unscaled[,i]*spec.pca.scale$scale[i]+spec.pca.scale$center[i])

plot(400:2400,recon.trial[1,],type="l",col="red",cex.lab=2,
     xlab="Wavelength",ylab="Reflectance",ylim=c(0,0.55))
points(400:2400,recon.trial[2,],type="l",col="black")
points(400:2400,recon.trial[3,],type="l",col="blue")

recon.trial.lib<-speclib(spectra=recon.trial,wavelength=400:2400)
spec.inv<-list()
for(i in 1:dim(recon.trial.lib)[1]){spec.inv[[i]]<-PROSPECTinvert(recon.trial.lib[i])}

###################################################
## continuum removal

spec.traits.cont<-continuumRemoval(spec.traits.ref,wav=400:2400,
                                   interpol="linear",method="substraction")
spec.pca.cont<-prcomp(~.,data=data.frame(spec.traits.cont[,which(apply(spec.traits.cont,2,sd)!=0)]),scale.=T)

plot(spectrolab::spectra(spec.traits.cont,bands=400:2400,names=names(spec.traits.samp)))
