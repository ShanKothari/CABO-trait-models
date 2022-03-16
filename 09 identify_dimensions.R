setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(dplyr)
library(hsdar)
library(truncnorm)
library(vegan)
library(intrinsicDimension)
library(lle)
library(ider)
library(paran)
library(pls)
library(reshape2)
library(patchwork)
library(ggfortify)
library(RDRToolbox)

ref.traits<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")
ref.traits<-ref.traits[which(meta(ref.traits)$project!="2019-Pardo-MSc-UdeM")]

meta(ref.traits)$chl_mass<-meta(ref.traits)$chlA_mass+meta(ref.traits)$chlB_mass

## area basis, in g/cm^2
meta(ref.traits)$chl_area<-meta(ref.traits)$chl_mass*meta(ref.traits)$LMA/10000
chl_area_norm<-lm(log(chl_area)~log(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$chl_norm<-resid(chl_area_norm)

meta(ref.traits)$logLMA<-log(meta(ref.traits)$LMA)
meta(ref.traits)$logEWT<-log(meta(ref.traits)$EWT)
meta(ref.traits)$logchl<-log(meta(ref.traits)$chl_mass)
meta(ref.traits)$logcar<-log(meta(ref.traits)$car_mass)
meta(ref.traits)$logN<-log(meta(ref.traits)$Nmass)

## remove samples with NAs for key traits
## so that sampling picks out samples with all trait values present
na.traits<-data.frame(which(is.na(meta(ref.traits)[c("EWT","LMA","Cmass","Nmass",
                               "hemicellulose_mass","cellulose_mass","lignin_mass",
                               "chl_mass","car_mass")]),arr.ind=T))
ref.traits<-ref.traits[-unique(na.traits$row),]

## sample no more than three of the same species
## to avoid over-weighting common species in importance
ref.traits.samp<-subset_by(ref.traits,by=as.vector(meta(ref.traits)$species),
                           n_min = 1,n_max = 10)

traits.samp<-meta(ref.traits.samp)[,c("EWT","LMA","Cmass","Nmass",
                                      "hemicellulose_mass","cellulose_mass","lignin_mass",
                                      "chl_mass","car_mass")]

traits.trans.samp<-meta(ref.traits.samp)[,c("logEWT","logLMA","Cmass","logN",
                                            "hemicellulose_mass","cellulose_mass",
                                            "lignin_mass","logchl","logcar")]

traits.osnas.samp<-meta(ref.traits.samp)[,c("logEWT","logLMA","Cnorm","Nnorm",
                                            "hemicellulose_norm","cellulose_norm",
                                            "lignin_norm","chl_norm","car_norm")]

## normalize all to mean 0, sd 1
ref.traits.norm<-apply(as.matrix(ref.traits.samp),2,function(x) (x-mean(x))/sd(x))
traits.norm<-apply(traits.samp,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm = T))
traits.trans.norm<-apply(traits.trans.samp,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm = T))
traits.osnas.norm<-apply(traits.osnas.samp,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm = T))

## choose bands to downsample to
select.bands<-c(seq(from=400,to=780,by=20),
                seq(from=800,to=1350,by=50),
                seq(from=1400,to=2400,by=25))

#################################################
## factor analysis of spectra

ref.traits.ref<-as.matrix(ref.traits.samp)
colnames(ref.traits.ref)<-paste("X",colnames(ref.traits.ref),sep="")
ref.pca<-prcomp(~.,data=as.data.frame(ref.traits.ref),scale.=T)
meta(ref.traits.samp)$PC1<-ref.pca$x[,1]
meta(ref.traits.samp)$PC2<-ref.pca$x[,2]
meta(ref.traits.samp)$PC3<-ref.pca$x[,3]
meta(ref.traits.samp)$PC4<-ref.pca$x[,4]

## average spectrum
# plot(400:2400,ref.pca$center,type="l")
## scree plot
screeplot(ref.pca)
screeplot(prcomp(~.,data=as.data.frame(ref.traits.ref[,select.bands-399]),scale.=T))
screeplot(prcomp(~.,data=as.data.frame(traits.trans.norm),scale.=T))

## reconstructing the spectrum
# recon.unscaled<-ref.pca$x %*% t(ref.pca$rotation)
# recon<-sapply(1:ncol(recon.unscaled),function(i) recon.unscaled[,i]*ref.pca$scale[i]+ref.pca$center[i])

# ## here we rescale all PCA axes to a standard deviation of 1
# ref.pca.scale<-ref.pca
# ref.pca.scale$x<-sapply(1:ncol(ref.pca$x),function(i) ref.pca$x[,i]/ref.pca$sdev[i])
# ref.pca.scale$rotation<-sapply(1:ncol(ref.pca$rotation),function(i) ref.pca$rotation[,i]*ref.pca$sdev[i])
# 
# ## matrix with the specified PC at -1, 0, or 1 and all other PCs at 0
# PC<-1
# trial.matrix<-matrix(data=0,nrow=3,ncol=ncol(ref.pca.scale$x))
# trial.matrix[,PC]<-c(-1,0,1)
# recon.trial.unscaled<-trial.matrix %*% t(ref.pca.scale$rotation)
# recon.trial<-sapply(1:ncol(recon.trial.unscaled),function(i) recon.trial.unscaled[,i]*ref.pca.scale$scale[i]+ref.pca.scale$center[i])
# 
# plot(400:2400,recon.trial[1,],type="l",col="red",cex.lab=2,
#      xlab="Wavelength",ylab="Reflectance",ylim=c(0,0.55))
# points(400:2400,recon.trial[2,],type="l",col="black")
# points(400:2400,recon.trial[3,],type="l",col="blue")
# 
# recon.trial.lib<-speclib(spectra=recon.trial,wavelength=400:2400)
# ref.inv<-list()
# for(i in 1:dim(recon.trial.lib)[1]){ref.inv[[i]]<-PROSPECTinvert(recon.trial.lib[i])}

################################
## nonlinear dimensionality reduction

ref.traits.dist<-dist(as.matrix(ref.traits.norm)[,select.bands-399],method = "manhattan")

## isomap
ref.isomap<-isomap(ref.traits.dist,ndim = 9,k=5)
# ref.isomap.proc<-ref.isomap

colnames(ref.isomap$points)<-paste("D",1:ncol(ref.isomap$points),sep="")
ref.isomap.coords<-as.data.frame(ref.isomap$points)
# ref.isomap.proc.coords<-data.frame(procrustes(X = traits.trans.norm, Y = ref.isomap.coords)$Yrot)
# ref.isomap.proc$points<-ref.isomap.proc.coords

ref.isomap.envfit.12<-envfit(ref.isomap,traits.trans.norm,choices = 1:2)
scores.12.df<-as.data.frame(vegan::scores(ref.isomap.envfit.12, "vectors")) * ordiArrowMul(ref.isomap.envfit.12)
ref.isomap.envfit.34<-envfit(ref.isomap,traits.trans.norm,choices = 3:4)
scores.34.df<-as.data.frame(vegan::scores(ref.isomap.envfit.34, "vectors")) * ordiArrowMul(ref.isomap.envfit.34)
rownames(scores.12.df)<-c("EWT","LMA","C","N","hemi","cell","lign","chl","car")
rownames(scores.34.df)<-c("EWT","LMA","C","N","hemi","cell","lign","chl","car")

ref.isomap.coords$functional.group<-as.factor(meta(ref.traits.samp)$functional.group)
ref.isomap.coords$LMA<-meta(ref.traits.samp)$LMA
ref.isomap.coords$EWT<-meta(ref.traits.samp)$EWT
ref.isomap.coords$Nmass<-meta(ref.traits.samp)$Nmass
ref.isomap.coords$chl<-meta(ref.traits.samp)$chl_mass
ref.isomap.coords$hemi<-meta(ref.traits.samp)$hemicellulose_mass
ref.isomap.coords$cell<-meta(ref.traits.samp)$cellulose_mass
ref.isomap.coords$lignin<-meta(ref.traits.samp)$lignin_mass
ref.isomap.coords$car<-meta(ref.traits.samp)$car_mass

ggplot(data=ref.isomap.coords,
       aes(x=D1,y=D2,color=functional.group))+
  geom_point(size=3)+
  geom_segment(aes(x = 0, y = 0, xend = D1*125, yend = D2*125),
               data = scores.12.df, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = scores.12.df, aes(x = D1*135, y = D2*135), 
            label = row.names(scores.12.df), colour = "navy", fontface = "bold") +
#  scale_color_viridis_c()+
  theme_bw()

ggplot(data=ref.isomap.coords,
       aes(x=D3,y=D4,color=functional.group))+
  geom_point(size=3)+
  geom_segment(aes(x = 0, y = 0, xend = D3*125, yend = D4*125),
               data = scores.34.df, size =1, alpha = 0.5, colour = "grey30") +
  geom_text(data = scores.34.df, aes(x = D3*135, y = D4*135), 
            label = row.names(scores.34.df), colour = "navy", fontface = "bold") +
  #  scale_color_viridis_c()+
  theme_bw()

# isomap.cca<-CCorA(Y = ref.isomap.coords[,c("D1","D2","D3","D4","D5")],
#                   X = data.frame(traits.trans.norm))
# biplot(isomap.cca)

# ref.isomap.proc.coords$functional.group<-as.factor(meta(ref.traits.samp)$functional.group)
# ref.isomap.proc.coords$logLMA<-data.frame(traits.trans.norm)$logLMA
# ref.isomap.proc.coords$logEWT<-data.frame(traits.trans.norm)$logEWT
# ref.isomap.proc.coords$logN<-data.frame(traits.trans.norm)$logN
# ref.isomap.proc.coords$Cmass<-data.frame(traits.trans.norm)$Cmass
# ref.isomap.proc.coords$logchl<-data.frame(traits.trans.norm)$logchl
# ref.isomap.proc.coords$hemi<-data.frame(traits.trans.norm)$hemicellulose_mass
# ref.isomap.proc.coords$cell<-data.frame(traits.trans.norm)$cellulose_mass
# ref.isomap.proc.coords$lignin<-data.frame(traits.trans.norm)$lignin_mass
# ref.isomap.proc.coords$R800<-ref.traits.samp[,800]
# 
# ref.isomap.envfit.12<-envfit(ref.isomap.proc,traits.trans.norm,choices = 1:2)
# scores.12.df<-as.data.frame(vegan::scores(ref.isomap.envfit.12, "vectors")) * ordiArrowMul(ref.isomap.envfit.12)
# ref.isomap.envfit.34<-envfit(ref.isomap.proc,traits.trans.norm,choices = 3:4)
# scores.34.df<-as.data.frame(vegan::scores(ref.isomap.envfit.34, "vectors")) * ordiArrowMul(ref.isomap.envfit.34)
# rownames(scores.12.df)<-c("EWT","LMA","C","N","hemi","cell","lign","chl","car")
# rownames(scores.34.df)<-c("EWT","LMA","C","N","hemi","cell","lign","chl","car")
# 
# ggplot(data=ref.isomap.proc.coords,
#        aes(x=X1,y=X2,color=logEWT))+
#   geom_point(size=3)+
# #  scale_color_viridis_c()+
#   theme_bw()+
#   geom_segment(aes(x = 0, y = 0, xend = X1*2, yend = X2*2),
#                data = scores.12.df, size =1, alpha = 0.5, colour = "grey30") +
#   geom_text(data = scores.12.df, aes(x = X1*2.1, y = X2*2.1), 
#             label = row.names(scores.12.df), colour = "navy", fontface = "bold")
# 
# ggplot(data=ref.isomap.proc.coords,
#        aes(x=X3,y=X4,color=logLMA))+
#   geom_point(size=3)+
# #  scale_color_discrete()+
#   theme_bw()+
#   geom_segment(aes(x = 0, y = 0, xend = X3*1.5, yend = X4*1.5),
#                data = scores.34.df, size =1, alpha = 0.5, colour = "grey30") +
#   geom_text(data = scores.34.df, aes(x = X3*1.6, y = X4*1.6), 
#             label = row.names(scores.34.df), colour = "navy", fontface = "bold")
# 
# isomap.proc.cca<-CCorA(Y = ref.isomap.proc[,c("X1","X2","X3","X4","X5")],
#                   X = data.frame(traits.trans.norm))
# biplot(isomap.cca)

## NMDS
ref.nmds<-metaMDS(ref.traits.norm[,select.bands-399],distance="manhattan",k=10,
                    autotransform = F)
stressplot(ref.nmds)

#################################
## estimating intrinsic dimensionality

## Grassberger-Procaccia; unsure about correct choice of p or k1/k2
corint(x=ref.traits.norm[,select.bands-399],p=20)
corint(x=na.omit(traits.norm),p=20)
## ML; not sure about correct choice of p or k1/k2
lbmle(x = ref.traits.norm[,select.bands-399], k1=5, k2=10)
lbmle(x = na.omit(traits.norm), k1=5, k2=10)
mada(x = ref.traits.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
mada(x = na.omit(traits.norm), k = NULL, comb = "average", maxDim = 10)
side(x = ref.traits.norm[,select.bands-399], comb = "average", maxDim = 10)
side(x = na.omit(traits.norm), comb = "average", maxDim = 10)
nni(x = ref.traits.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = na.omit(traits.norm), k1 = 5, k2 = 30, eps = 0.01, p = NULL)
dancoDimEst(data = ref.traits.norm[,select.bands-399], k=5, D=10, ver = "DANCo")
dancoDimEst(data = na.omit(traits.norm), k=5, D=10, ver = "DANCo")
dancoDimEst(data = ref.traits.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")
dancoDimEst(data = na.omit(traits.norm), k=5, D=10, ver = "MIND_MLi")

ref.isomap<-isomap(dist(ref.traits.norm[,select.bands-399],method="manhattan"),ndim = 2,k = 5)
trait.isomap<-isomap(dist(na.omit(traits.norm),method="manhattan"),ndim = 2,k = 5)

calc_k(ref.traits.norm[,select.bands-399], m=10, kmin=1, kmax=20, plotres=TRUE, 
       parallel=FALSE, cpus=2, iLLE=FALSE)
calc_k(na.omit(traits.norm), m=10, kmin=1, kmax=20, plotres=TRUE, 
       parallel=FALSE, cpus=2, iLLE=FALSE) 
ref.lle<-lle(ref.traits.norm[,select.bands-399],m=10,k=8,id = T,iLLE=T)
traits.lle<-lle(na.omit(traits.norm),m=10,k=8,id = T,iLLE=T)

#########################################
## dimensionality of fresh, pressed, and ground spectra

fresh.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/fresh_spec_EL_agg.rds")
pressed.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/pressed_spec_EL_agg.rds")
ground.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/ground_spec_EL_agg.rds")

meta(fresh.spec)$chl<-meta(fresh.spec)$chlA+meta(fresh.spec)$chlB

## sample no more than ten of the same species
## to avoid over-weighting common species in importance
fresh.spec.samp<-subset_by(fresh.spec,by=as.vector(meta(fresh.spec)$Species),
                           n_min = 1,n_max = 10)
pressed.spec.samp<-subset_by(pressed.spec,by=as.vector(meta(pressed.spec)$Species),
                             n_min = 1,n_max = 10)
ground.spec.samp<-subset_by(ground.spec,by=as.vector(meta(ground.spec)$Species),
                            n_min = 1,n_max = 10)

fresh.spec.norm<-apply(as.matrix(fresh.spec.samp),2,function(x) (x-mean(x))/sd(x))
pressed.spec.norm<-apply(as.matrix(pressed.spec.samp),2,function(x) (x-mean(x))/sd(x))
ground.spec.norm<-apply(as.matrix(ground.spec.samp),2,function(x) (x-mean(x))/sd(x))

traits.sub<-apply(meta(fresh.spec.samp)[,c("LDMC","EWT","LMA","C","N",
                                      "hemicellulose","cellulose","lignin",
                                      "chl","car","Ca","K","P")],2,
                  function(x) (x-mean(x,na.rm=T))/sd(x,na.rm=T))

colnames(fresh.spec.norm)<-paste("X",colnames(fresh.spec.norm),sep="")
colnames(pressed.spec.norm)<-paste("X",colnames(pressed.spec.norm),sep="")
colnames(ground.spec.norm)<-paste("X",colnames(ground.spec.norm),sep="")

screeplot(prcomp(~.,data=as.data.frame(fresh.spec.norm[,select.bands-399]),scale.=T))
screeplot(prcomp(~.,data=as.data.frame(pressed.spec.norm[,select.bands-399]),scale.=T))
screeplot(prcomp(~.,data=as.data.frame(ground.spec.norm[,select.bands-399]),scale.=T))

corint(x=fresh.spec.norm[,select.bands-399],p=20)
corint(x=pressed.spec.norm[,select.bands-399],p=20)
corint(x=ground.spec.norm[,select.bands-399],p=20)
corint(x=na.omit(traits.sub),p=20)

## with k1=k2=k, this yields the same result as
## maxLikGlobalDimEst
lbmle(x = fresh.spec.norm[,select.bands-399], k1=5, k2=10)
lbmle(x = pressed.spec.norm[,select.bands-399], k1=5, k2=10)
lbmle(x = ground.spec.norm[,select.bands-399], k1=5, k2=10)
lbmle(x = na.omit(traits.sub), k1=5, k2=10)

mada(x = fresh.spec.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
mada(x = pressed.spec.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
mada(x = ground.spec.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
mada(x = na.omit(traits.sub), k = NULL, comb = "average", maxDim = 10)

side(x = fresh.spec.norm[,select.bands-399],maxDim=10)
side(x = pressed.spec.norm[,select.bands-399],maxDim=10)
side(x = ground.spec.norm[,select.bands-399],maxDim=10)
side(x = na.omit(traits.sub),maxDim=10)

nni(x = fresh.spec.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = pressed.spec.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = ground.spec.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = na.omit(traits.sub), k1 = 5, k2 = 30, eps = 0.01, p = NULL)

dancoDimEst(data = fresh.spec.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")
dancoDimEst(data = pressed.spec.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")
dancoDimEst(data = ground.spec.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")
dancoDimEst(data = na.omit(traits.sub), k=5, D=10, ver = "MIND_MLi")

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

N<-runif(n = 1000,min=1,max=2)
Cw<-sample(na.omit(meta(ref.traits.samp)$EWT/10),
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
colnames(prospect.sim.norm)<-paste("X",400:2500,sep="")
screeplot(prcomp(~.,data=as.data.frame(prospect.sim.norm[,select.bands-399]),scale.=T))

## Grassberger-Procaccia; unsure about correct choice of p or k1/k2
corint(x=prospect.sim.norm[,select.bands-399],p=20)
## ML; not sure about correct choice of p or k1/k2
lbmle(x = prospect.sim.norm[,select.bands-399], k1=5, k2=10)
mada(x = prospect.sim.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
side(x = prospect.sim.norm[,select.bands-399],maxDim=10)
nni(x = prospect.sim.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
dancoDimEst(data = prospect.sim.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")

calc_k(prospect.sim.norm[,select.bands-399], m=10, kmin=1, kmax=20, plotres=TRUE, 
       parallel=FALSE, cpus=2, iLLE=FALSE) 
prospect.lle<-lle(prospect.sim.norm[,select.bands-399],m=4,k=6,id = T)
prospect.isomap<-isomap(dist(prospect.sim.norm,method="euclidean"),
                        ndim = 10,k=3)
plot(prospect.isomap$points[,3]~prospect.params$Cab)
plot(prospect.isomap$points[,1]~prospect.isomap$points[,2])

############################################
## trait covariance simulation

FG.sizes<-c(50,50,50,25,25)
FG<-c(rep("FG1",50),rep("FG2",50),rep("FG3",50),rep("FG4",25),rep("FG5",25))

Cm1<-rtruncnorm(50,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.125),sd=0.001)
Cm2<-rtruncnorm(50,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.375),sd=0.0015)
Cm3<-rtruncnorm(50,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.625),sd=0.002)
Cm4<-rtruncnorm(25,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.75),sd=0.0025)
Cm5<-rtruncnorm(25,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.875),sd=0.003)
Cm<-c(Cm1,Cm2,Cm3,Cm4,Cm5)

N<-1.5
Cw<-0.007
Cab<-40
prospect.params<-data.frame(N=N,Cw=Cw,Cab=Cab,Cm=Cm,FG=FG)
prospect.params$Car<-prospect.params$Cab/6.24 ## mean chl:car ratio
prospect.params$Anth<-0
prospect.params$Cbrown<-0

colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

Cm_density<-ggplot(data=prospect.params,
                   aes(x=Cm,color=FG))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=15))+
  labs(y="Density",x=expression("LMA (g cm"^-2*")"),tag="A")+
  coord_cartesian(xlim=c(0,0.021))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

prospect.sim<-t(apply(prospect.params,1,
                      function(x) {
                        spec.sim<-PROSPECT(N=x["N"],Cab=x["Cab"],Car=x["Car"],
                                           Anth=x["Anth"],Cbrown=x["Cbrown"],
                                           Cm=x["Cm"],Cw=x["Cw"])
                        spec<-spec.sim@spectra@spectra_ma
                        return(spec)
                      }))

prospect.spec<-spectrolab::spectra(prospect.sim,bands=400:2500,names=1:200)
prospect.spec.df<-data.frame(sample_id=names(prospect.spec),
                             FG=prospect.params$FG,
                             as.matrix(prospect.spec))

prospect.spec.long<-melt(prospect.spec.df,id.vars = c("sample_id","FG"))
prospect.spec.long$variable<-as.numeric(gsub(pattern = "X",replacement = "",x = prospect.spec.long$variable))

prospect_fg_spec_plot<-ggplot(prospect.spec.long,
                              aes(x=variable,y=value,color=FG))+
  stat_summary(fun=median,na.rm=T,geom="line",size=1)+
  theme_bw()+theme(text=element_text(size=15),
                   plot.margin=unit(c(0,0.2,0,0),"in"))+
  labs(tag = "B",x="Wavelength (nm)",y="Reflectance")+
  guides(color=guide_legend("Functional group",nrow = 2))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.5))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2510))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])

protoT1<-prospect.params$Cm+rnorm(200,mean=0,sd=0.003)
prospect.params$T1<-(max(prospect.params$Cm)-min(prospect.params$Cm))*(protoT1-min(protoT1))/(max(protoT1)-min(protoT1))+min(prospect.params$Cm)
prospect.params$FG.mean<-unlist(apply(prospect.params,1,function(x) median(prospect.params[["Cm"]][prospect.params[["FG"]]==x[["FG"]]])))
protoT2<-prospect.params$FG.mean+rnorm(200,mean=0,sd=0.001)
prospect.params$T2<-(max(prospect.params$Cm)-min(prospect.params$Cm))*(protoT2-min(protoT2))/(max(protoT2)-min(protoT2))+min(prospect.params$Cm)

T1_Cm_plot<-ggplot(prospect.params,aes(y=T1,x=Cm,color=FG))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
#  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.02),ylim=c(0,0.02))+
  theme(text = element_text(size=15),
        legend.position = c(0.8, 0.2))+
  labs(tag="C",y="FT1",x=expression("LMA (g cm"^-2*")"))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

T2_Cm_plot<-ggplot(prospect.params,aes(y=T2,x=Cm,color=FG))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  #  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.02),ylim=c(0,0.02))+
  theme(text = element_text(size=15),
        legend.position = c(0.8, 0.2))+
  labs(tag="D",y="FT2",x=expression("LMA (g cm"^-2*")"))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

Cm_model<-plsr(prospect.params$Cm~as.matrix(prospect.sim),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_Cm_model <- selectNcomp(Cm_model, method = "onesigma", plot = FALSE)
Cm_pred <- data.frame(measured=prospect.params$Cm,
                      val_pred=Cm_model$validation$pred[,,ncomp_Cm_model],
                      FG=prospect.params$FG)

ggplot(Cm_pred,aes(y=measured,x=val_pred,color=FG))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  coord_cartesian(xlim=c(0,0.02),ylim=c(0,0.02))+
  labs(y="Measured",x="Predicted")+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])

T1_model<-plsr(prospect.params$T1~as.matrix(prospect.sim),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_T1_model <- selectNcomp(T1_model, method = "onesigma", plot = FALSE)
T1_pred <- data.frame(measured=prospect.params$T1,
                      val_pred=T1_model$validation$pred[,,ncomp_T1_model],
                      FG=prospect.params$FG)

T1_plot<-ggplot(T1_pred,aes(y=measured,x=val_pred,color=FG))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.02),ylim=c(0,0.02))+
  theme(text = element_text(size=15),
        legend.position = c(0.8, 0.2))+
  labs(tag="E",y="Measured FT1",x="Predicted FT1")+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

T2_model<-plsr(prospect.params$T2~as.matrix(prospect.sim),
               ncomp=30,method = "oscorespls",validation="CV",segments=10)
ncomp_T2_model <- selectNcomp(T2_model, method = "onesigma", plot = FALSE)
T2_pred <- data.frame(measured=prospect.params$T2,
                      val_pred=T2_model$validation$pred[,,ncomp_T2_model],
                      FG=prospect.params$FG)

T2_plot<-ggplot(T2_pred,aes(y=measured,x=val_pred,color=FG))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.02),ylim=c(0,0.02))+
  theme(text = element_text(size=15),
        legend.position = c(0.8, 0.2))+
  labs(tag="F",y="Measured FT2",x="Predicted FT2")+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

pdf("GrimeReview/prospect_covar.pdf",height=12,width=11)
Cm_density+prospect_fg_spec_plot+
  T1_Cm_plot+T2_Cm_plot+
  T1_plot+T2_plot+
  plot_layout(ncol=2,guides="collect",heights=c(0.5,1,1)) &
  theme(legend.position = "bottom")
dev.off()

source("Scripts/VIP.R")
plot(400:2500,VIP(T2_model)[ncomp_T2_model,])
