setwd("C:/Users/kotha020/Dropbox/PostdocProjects/FreshLeafModels")

library(spectrolab)
library(hsdar)
library(truncnorm)
library(vegan)
library(intrinsicDimension)
library(lle)
library(ider)
library(pls)
library(reshape2)
library(patchwork)
library(dimRed)
library(ggrepel)

ref.traits<-readRDS("ProcessedSpectra/all_ref_and_traits.rds")
## remove the Pardo project, which has no trait data
ref.traits<-ref.traits[which(meta(ref.traits)$project!="2019-Pardo-MSc-UdeM")]

## color-blind friendly palette
colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

## sum chlA and chlB to get total chl
meta(ref.traits)$chl_mass<-meta(ref.traits)$chlA_mass+meta(ref.traits)$chlB_mass
## area basis, in g/cm^2
meta(ref.traits)$chl_area<-meta(ref.traits)$chl_mass*meta(ref.traits)$LMA/10000

###################################
## calculate normalization-independent traits
## sensu Osnas et al. 2013 Science and 2018 PNAS

Narea_norm<-lm(log10(Narea)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Nnorm<-resid(Narea_norm)

Carea_norm<-lm(log10(Carea)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Cnorm<-resid(Carea_norm)

solubles_area_norm<-lm(log10(solubles_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$solubles_norm<-resid(solubles_area_norm)

hemicellulose_area_norm<-lm(log10(hemicellulose_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$hemicellulose_norm<-resid(hemicellulose_area_norm)

cellulose_area_norm<-lm(log10(cellulose_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$cellulose_norm<-resid(cellulose_area_norm)

lignin_area_norm<-lm(log10(lignin_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$lignin_norm<-resid(lignin_area_norm)

chlA_area_norm<-lm(log10(chlA_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$chlA_norm<-resid(chlA_area_norm)

chlB_area_norm<-lm(log10(chlB_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$chlB_norm<-resid(chlB_area_norm)

chl_area_norm<-lm(log(chl_area)~log(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$chl_norm<-resid(chl_area_norm)

car_area_norm<-lm(log10(car_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$car_norm<-resid(car_area_norm)

meta(ref.traits)$Al_area_omit<-meta(ref.traits)$Al_area
meta(ref.traits)$Al_area_omit[which(meta(ref.traits)$Al_area_omit==0)]<-NA
Al_area_norm<-lm(log10(Al_area_omit)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Al_norm<-resid(Al_area_norm)
meta(ref.traits)$Al_area_omit<-NULL

Ca_area_norm<-lm(log10(Ca_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Ca_norm<-resid(Ca_area_norm)

meta(ref.traits)$Cu_area_omit<-meta(ref.traits)$Cu_area
meta(ref.traits)$Cu_area_omit[which(meta(ref.traits)$Cu_area_omit==0)]<-NA
Cu_area_norm<-lm(log10(Cu_area_omit)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Cu_norm<-resid(Cu_area_norm)
meta(ref.traits)$Cu_area_omit<-NULL

Fe_area_norm<-lm(log10(Fe_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Fe_norm<-resid(Fe_area_norm)

K_area_norm<-lm(log10(K_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$K_norm<-resid(K_area_norm)

Mg_area_norm<-lm(log10(Mg_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Mg_norm<-resid(Mg_area_norm)

meta(ref.traits)$Mn_area_omit<-meta(ref.traits)$Mn_area
meta(ref.traits)$Mn_area_omit[which(meta(ref.traits)$Mn_area_omit==0)]<-NA
Mn_area_norm<-lm(log10(Mn_area_omit)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Mn_norm<-resid(Mn_area_norm)
meta(ref.traits)$Mn_area_omit<-NULL

meta(ref.traits)$Na_area_omit<-meta(ref.traits)$Na_area
meta(ref.traits)$Na_area_omit[which(meta(ref.traits)$Na_area_omit==0)]<-NA
Na_area_norm<-lm(log10(Na_area_omit)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Na_norm<-resid(Na_area_norm)
meta(ref.traits)$Na_area_omit<-NULL

P_area_norm<-lm(log10(P_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$P_norm<-resid(P_area_norm)

Zn_area_norm<-lm(log10(Zn_area)~log10(0.1*LMA),data=meta(ref.traits),na.action=na.exclude)
meta(ref.traits)$Zn_norm<-resid(Zn_area_norm)

################################################
## PCA of traits for trait modeling paper (not Grime Review)

trait.pca<-prcomp(~.,meta(ref.traits)[,c("LMA","LDMC","EWT","Cmass","Nmass",
                                      "solubles_mass","hemicellulose_mass",
                                      "cellulose_mass","lignin_mass","chlA_mass",
                                      "chlB_mass","car_mass")],
                  scale.=T,na.action = na.exclude)

trait.pca.val <- data.frame(functional.group=meta(ref.traits)$functional.group,trait.pca$x)
trait.pca.loadings<-data.frame(variables = rownames(trait.pca$rotation), trait.pca$rotation)
trait.pca.loadings$variables<-gsub("_mass","",trait.pca.loadings$variables)
trait.pca.loadings$variables<-gsub("mass","",trait.pca.loadings$variables)
trait.pca.perc<-trait.pca$sdev^2/sum(trait.pca$sdev^2)*100

trait.pca.plot<-ggplot(trait.pca.val, aes(x = PC1, y = PC2, color = functional.group)) +
  geom_point(size = 2,alpha=0.6) +
  geom_segment(data = trait.pca.loadings, aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10),
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black",size=1) +
  annotate("text",
           x = trait.pca.loadings$PC1*12+c(0,0,0,0,0.3,0,0,0,0,-0.4,0.3,0),
           y = trait.pca.loadings$PC2*12+c(0,0,0,0,-0.4,0,0,0,0,0.3,-0.1,0),
           label = trait.pca.loadings$variables)+
  theme_bw()+theme(text=element_text(size=15))+
  guides(color=guide_legend("Functional group"))+
  coord_fixed(ratio=trait.pca.perc[2]/trait.pca.perc[1])+
  labs(x=paste("PC1 (",round(trait.pca.perc[1],1),"% variance)",sep=""),
       y=paste("PC2 (",round(trait.pca.perc[2],1),"% variance)",sep=""))+
  scale_color_manual(values=colorBlind)

## log-transform LMA and EWT to reduce skewness
meta(ref.traits)$logLMA<-log(meta(ref.traits)$LMA)
meta(ref.traits)$logEWT<-log(meta(ref.traits)$EWT)

## using normalization-independent traits
## and leaving out some traits strongly correlated with others
traitnorm.pca<-prcomp(~.,meta(ref.traits)[,c("logLMA","LDMC","logEWT","Cnorm","Nnorm",
                                             "hemicellulose_norm","cellulose_norm",
                                             "lignin_norm","chlA_norm")],
                      scale.=T,na.action = na.exclude)

traitnorm.pca.val <- data.frame(functional.group=meta(ref.traits)$functional.group,traitnorm.pca$x)
traitnorm.pca.loadings<-data.frame(variables = rownames(traitnorm.pca$rotation), traitnorm.pca$rotation)
traitnorm.pca.loadings$variables<-gsub("_norm","",traitnorm.pca.loadings$variables)
traitnorm.pca.loadings$variables<-gsub("norm","",traitnorm.pca.loadings$variables)
traitnorm.pca.perc<-traitnorm.pca$sdev^2/sum(traitnorm.pca$sdev^2)*100

traitnorm.pca.plot<-ggplot(traitnorm.pca.val, aes(x = PC1*2, y = -PC2*2, color = functional.group)) +
  geom_point(size = 2,alpha=0.6) +
  geom_segment(data = traitnorm.pca.loadings, aes(x = 0, y = 0, xend = PC1*10, yend = -PC2*10),
               arrow = arrow(length = unit(1/2, "picas")),
               color = "black",size=1) +
  annotate("text",
           x = traitnorm.pca.loadings$PC1*12+c(0,0,0,-0.1,0,-0.2,-0.2,0,-0.1),
           y = -traitnorm.pca.loadings$PC2*12+c(-0.1,0,0,0.3,-0.2,-0.2,0.3,-0.2,0),
           label = traitnorm.pca.loadings$variables)+
  theme_bw()+theme(text=element_text(size=15))+
  guides(color=guide_legend("Functional group"))+
  coord_fixed(ratio=traitnorm.pca$sdev[2]^2/traitnorm.pca$sdev[1]^2)+
  labs(x=paste("PC1 (",round(traitnorm.pca.perc[1],1),"% variance)",sep=""),
       y=paste("PC2 (",round(traitnorm.pca.perc[2],1),"% variance)",sep=""))+
  scale_color_manual(values=colorBlind)

pdf("Images/pca_plot.pdf")
trait.pca.plot+traitnorm.pca.plot+plot_layout(ncol=1,guides="collect") &
  theme(legend.position = "bottom")
dev.off()

#################################################
## preparing data for dimensionality analyses in Grime Review

## log-transform LMA and EWT to reduce skewness
meta(ref.traits)$logLMA<-log(meta(ref.traits)$LMA)
meta(ref.traits)$logEWT<-log(meta(ref.traits)$EWT)

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

## subsetting traits of interest
traits.samp<-meta(ref.traits.samp)[,c("EWT","LMA","Cmass","Nmass",
                                      "hemicellulose_mass","cellulose_mass","lignin_mass",
                                      "chl_mass","car_mass")]

## and the normalization-independent version
## also using the log-transformed EWT and LMA
traits.ni.samp<-meta(ref.traits.samp)[,c("logEWT","logLMA","Cnorm","Nnorm",
                                         "hemicellulose_norm","cellulose_norm",
                                         "lignin_norm","chl_norm","car_norm")]

## normalize all to mean 0, sd 1
ref.traits.norm<-apply(as.matrix(ref.traits.samp),2,function(x) (x-mean(x))/sd(x))
traits.norm<-apply(traits.samp,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm = T))
traits.ni.norm<-apply(traits.ni.samp,2,function(x) (x-mean(x,na.rm=T))/sd(x,na.rm = T))
## note that normalization-independent traits should already have mean ~0
## but not exactly because we did the renormalization prior to choosing
## 10 samples from species

## choose bands to downsample to
select.bands<-c(seq(from=400,to=780,by=20),
                seq(from=800,to=1350,by=50),
                seq(from=1400,to=2400,by=25))

#################################################
## factor analysis of spectra

## PCA scree plots
## because the scaling is done by the PCA function
## we can use the non-z-standardized data
ref.traits.mat<-as.matrix(ref.traits.samp)
colnames(ref.traits.mat)<-paste("X",colnames(ref.traits.mat),sep="")
screeplot(prcomp(~.,data=as.data.frame(ref.traits.mat[,select.bands-399]),scale.=T))
screeplot(prcomp(~.,data=as.data.frame(traits.ni.samp),scale.=T))

################################
## nonlinear dimensionality reduction

ref.traits.dist<-dist(as.matrix(ref.traits.norm)[,select.bands-399],method = "manhattan")

## isomap
ref.isomap<-isomap(ref.traits.dist,ndim = 9,k=5)
colnames(ref.isomap$points)<-paste("D",1:ncol(ref.isomap$points),sep="")
ref.isomap.coords<-as.data.frame(ref.isomap$points)

## procrustes
ref.isomap.proc<-ref.isomap
ref.isomap.proc$points<-procrustes(Y = ref.isomap.coords,X = traits.ni.norm)$Yrot
colnames(ref.isomap.proc$points)<-paste("D",1:ncol(ref.isomap.proc$points),sep="")
ref.isomap.proc.coords<-as.data.frame(ref.isomap.proc$points)

## plotting without Procrustes analysis with traits overlain
# ref.isomap.envfit.12<-envfit(ref.isomap,traits.ni.norm,choices = 1:2)
# scores.12.df<-as.data.frame(vegan::scores(ref.isomap.envfit.12, "vectors"))
# ref.isomap.envfit.34<-envfit(ref.isomap,traits.ni.norm,choices = 3:4)
# scores.34.df<-as.data.frame(vegan::scores(ref.isomap.envfit.34, "vectors"))
# rownames(scores.12.df)<-c("EWT","LMA","C","N","hemi","cell","lign","chl","car")
# rownames(scores.34.df)<-c("EWT","LMA","C","N","hemi","cell","lign","chl","car")
# 
# ref.isomap.coords$functional.group<-as.factor(meta(ref.traits.samp)$functional.group)
# 
# isomap.12<-ggplot(data=ref.isomap.coords,
#        aes(x=D1,y=D2,color=functional.group))+
#   geom_point(size=3)+
#   geom_segment(aes(x = 0, y = 0, xend = D1*150, yend = D2*150),
#                data = scores.12.df, size =1, alpha = 0.5, colour = "grey30") +
#   geom_text(data = scores.12.df, aes(x = D1*160, y = D2*160),
#             label = row.names(scores.12.df), colour = "black", fontface = "bold",
#             size=7) +
#   theme_bw()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         text = element_text(size=20))+
#   scale_color_manual(values=colorBlind)+
#   guides(color=guide_legend("Functional group"))+
#   labs(x="Spectral axis 1",y="Spectral axis 2")
# 
# isomap.34<-ggplot(data=ref.isomap.coords,
#        aes(x=D3,y=D4,color=functional.group))+
#   geom_point(size=3)+
#   geom_segment(aes(x = 0, y = 0, xend = D3*110, yend = D4*110),
#                data = scores.34.df, size =1, alpha = 0.5, colour = "grey30") +
#   geom_text(data = scores.34.df, aes(x = D3*120, y = D4*120),
#             label = row.names(scores.34.df), colour = "black", fontface = "bold",
#             size=7) +
#   theme_bw()+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         text = element_text(size=20))+
#   scale_color_manual(values=colorBlind)+
#   guides(color=guide_legend("Functional group"))+
#   labs(x="Spectral axis 3",y="Spectral axis 4")
# 
# # canonical correspondence analysis with traits
# isomap.cca<-CCorA(Y = ref.isomap.coords[,c("D1","D2","D3","D4","D5")],
#                   X = data.frame(traits.ni.norm))
# biplot(isomap.cca)
# 
# pdf("GrimeReview/isomap_no_procrustes_horizontal2.pdf",width=7,height=7)
# isomap.34+
#   theme(legend.position = "bottom")
# dev.off()

## plotting with Procrustes analysis
# add functional group for visualization
ref.isomap.proc.coords$functional.group<-as.factor(meta(ref.traits.samp)$functional.group)

# extract scores matrices
ref.isomap.envfit.12<-envfit(ref.isomap.proc,traits.ni.norm,choices = 1:2)
scores.12.df<-as.data.frame(vegan::scores(ref.isomap.envfit.12, "vectors")) * ordiArrowMul(ref.isomap.envfit.12)
ref.isomap.envfit.34<-envfit(ref.isomap.proc,traits.ni.norm,choices = 3:4)
scores.34.df<-as.data.frame(vegan::scores(ref.isomap.envfit.34, "vectors")) * ordiArrowMul(ref.isomap.envfit.34)
rownames(scores.12.df)<-c("EWT","LMA","C","N","hemi","cell","lign","chl","car")
rownames(scores.34.df)<-c("EWT","LMA","C","N","hemi","cell","lign","chl","car")

isomap.12<-ggplot(data=ref.isomap.proc.coords,
                  aes(x=D1,y=D2,color=functional.group))+
  geom_point(size=3)+
  geom_segment(aes(x = 0, y = 0, xend = D1*2, yend = D2*2),
               data = scores.12.df, size =1.5, alpha = 0.5, colour = "black") +
  geom_text_repel(data = scores.12.df, aes(x = D1*2.4, y = D2*2.4), size=7,
            label = row.names(scores.12.df), colour = "black", fontface = "bold",
            max.overlaps = 15,force = 2,force_pull = 3) +
  theme_bw()+coord_fixed()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=15))+
  scale_color_manual(values=colorBlind)+
  guides(color=guide_legend("Functional group"))+
  labs(x="Spectral Axis 1",y="Spectral Axis 2")

isomap.34<-ggplot(data=ref.isomap.proc.coords,
                  aes(x=D3,y=D4,color=functional.group))+
  geom_point(size=3)+
  geom_segment(aes(x = 0, y = 0, xend = D3, yend = D4),
               data = scores.34.df, size =1.5, alpha = 0.5, colour = "black") +
  geom_text(data = scores.34.df, aes(x = D3*1.2, y = D4*1.2), size=7,
            label = row.names(scores.34.df), colour = "black", fontface = "bold",
            nudge_x = c(0,-0.1,-0.05,0,-0.15,0,0.05,0.05,-0.05),
            nudge_y = c(0,0,0,-0.05,0,0,0,0,0)) +
  theme_bw()+coord_fixed()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=15))+
  scale_color_manual(values=colorBlind)+
  guides(color=guide_legend("Functional group"))+
  labs(x="Spectral Axis 3",y="Spectral Axis 4")

pdf("GrimeReview/isomap_horizontal.pdf",width=11.5,height=6.5)
isomap.12+isomap.34+
  plot_layout(ncol=2,guides="collect") &
  theme(legend.position = "bottom")
dev.off()

#################################
## estimating intrinsic dimensionality

## Grassberger-Procaccia correlation integral (1983)
corint(x=ref.traits.norm[,select.bands-399],p=20)
corint(x=na.omit(traits.ni.norm),p=20)
## Levina-Bickel maximum likelihood (2004)
lbmle(x = ref.traits.norm[,select.bands-399], k1=5, k2=10)
lbmle(x = na.omit(traits.ni.norm), k1=5, k2=10)
## Manifold-adaptive from Farahmand et al. 2007
mada(x = ref.traits.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
mada(x = na.omit(traits.ni.norm), k = NULL, comb = "average", maxDim = 10)
## Nearest-neighbor from Pettis et al. 1979
nni(x = ref.traits.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = na.omit(traits.ni.norm), k1 = 5, k2 = 30, eps = 0.01, p = NULL)
## minimum neighbor distance-maximum likelihood from Rozza et al.  (2012)
dancoDimEst(data = ref.traits.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")
dancoDimEst(data = na.omit(traits.ni.norm), k=5, D=10, ver = "MIND_MLi")

ref.isomap<-isomap(dist(ref.traits.norm[,select.bands-399],method="manhattan"),ndim = 10,k = 5)
trait.isomap<-isomap(dist(na.omit(traits.ni.norm),method="manhattan"),ndim = 10,k = 5)

calc_k(ref.traits.norm[,select.bands-399], m=10, kmin=1, kmax=20, plotres=TRUE, 
       parallel=FALSE, cpus=2, iLLE=FALSE)
calc_k(na.omit(traits.ni.norm), m=10, kmin=1, kmax=20, plotres=TRUE, 
       parallel=FALSE, cpus=2, iLLE=FALSE) 
ref.lle<-lle(ref.traits.norm[,select.bands-399],m=10,k=8,id = T,iLLE=T)
traits.lle<-lle(na.omit(traits.norm),m=10,k=8,id = T,iLLE=T)

#########################################
## dimensionality of fresh, pressed, and ground spectra

fresh.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/fresh_spec_all.rds")
pressed.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/pressed_spec_all.rds")
ground.spec<-readRDS("../../TraitModels2018/HerbariumPaper/ProcessedSpectralData/ground_spec_all.rds")

## sample no more than ten of the same species
## to avoid over-weighting common species in importance
fresh.spec.samp<-subset_by(fresh.spec,by=as.vector(meta(fresh.spec)$Species),
                           n_min = 1,n_max = 10)
pressed.spec.samp<-subset_by(pressed.spec,by=as.vector(meta(pressed.spec)$Species),
                             n_min = 1,n_max = 10)
ground.spec.samp<-subset_by(ground.spec,by=as.vector(meta(ground.spec)$Species),
                            n_min = 1,n_max = 10)

## z-standardize spectral bands
fresh.spec.norm<-apply(as.matrix(fresh.spec.samp),2,function(x) (x-mean(x))/sd(x))
pressed.spec.norm<-apply(as.matrix(pressed.spec.samp),2,function(x) (x-mean(x))/sd(x))
ground.spec.norm<-apply(as.matrix(ground.spec.samp),2,function(x) (x-mean(x))/sd(x))

colnames(fresh.spec.norm)<-paste("X",colnames(fresh.spec.norm),sep="")
colnames(pressed.spec.norm)<-paste("X",colnames(pressed.spec.norm),sep="")
colnames(ground.spec.norm)<-paste("X",colnames(ground.spec.norm),sep="")

## the dimensionality estimation methods are all mentioned above
## with names/references

screeplot(prcomp(~.,data=as.data.frame(fresh.spec.norm[,select.bands-399]),scale.=T))
screeplot(prcomp(~.,data=as.data.frame(pressed.spec.norm[,select.bands-399]),scale.=T))
screeplot(prcomp(~.,data=as.data.frame(ground.spec.norm[,select.bands-399]),scale.=T))

corint(x=fresh.spec.norm[,select.bands-399],p=20)
corint(x=pressed.spec.norm[,select.bands-399],p=20)
corint(x=ground.spec.norm[,select.bands-399],p=20)

lbmle(x = fresh.spec.norm[,select.bands-399], k1=5, k2=10)
lbmle(x = pressed.spec.norm[,select.bands-399], k1=5, k2=10)
lbmle(x = ground.spec.norm[,select.bands-399], k1=5, k2=10)

mada(x = fresh.spec.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
mada(x = pressed.spec.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
mada(x = ground.spec.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)

nni(x = fresh.spec.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = pressed.spec.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
nni(x = ground.spec.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)

dancoDimEst(data = fresh.spec.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")
dancoDimEst(data = pressed.spec.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")
dancoDimEst(data = ground.spec.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")

lle.fresh<-lle(as.matrix(fresh.spec),m=10,k=10,id = T)
lle.pressed<-lle(as.matrix(pressed.spec),m=10,k=10,id = T)
lle.ground<-lle(as.matrix(ground.spec),m=10,k=10,id = T)

isomap.fresh<-isomap(dist(as.matrix(fresh.spec),method = "manhattan"),
                     ndim = 10,k=5)
isomap.pressed<-isomap(dist(as.matrix(pressed.spec),method = "manhattan"),
                     ndim = 10,k=5)
isomap.ground<-isomap(dist(as.matrix(ground.spec),method = "manhattan"),
                     ndim = 10,k=5)

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

corint(x=prospect.sim.norm[,select.bands-399],p=20)
lbmle(x = prospect.sim.norm[,select.bands-399], k1=5, k2=10)
mada(x = prospect.sim.norm[,select.bands-399], k = NULL, comb = "average", maxDim = 10)
nni(x = prospect.sim.norm[,select.bands-399], k1 = 5, k2 = 30, eps = 0.01, p = NULL)
dancoDimEst(data = prospect.sim.norm[,select.bands-399], k=5, D=10, ver = "MIND_MLi")

calc_k(prospect.sim.norm[,select.bands-399], m=10, kmin=1, kmax=20, plotres=TRUE, 
       parallel=FALSE, cpus=2, iLLE=FALSE) 
prospect.lle<-lle(prospect.sim.norm[,select.bands-399],m=4,k=6,id = T)
prospect.isomap<-isomap(dist(prospect.sim.norm,method="manhattan"),
                        ndim = 10,k=3)

############################################
## trait covariance simulation
## NOTE: this section is functionally independent of the section
## immediately above and some of the variable names overlap

## pick the sizes of functional groups
FG.sizes<-c(50,50,50,25,25)
FG<-c(rep("FG1",50),rep("FG2",50),rep("FG3",50),rep("FG4",25),rep("FG5",25))

## generate LMA distribution for each functional group
Cm1<-rtruncnorm(50,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.125),sd=0.001)
Cm2<-rtruncnorm(50,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.375),sd=0.0015)
Cm3<-rtruncnorm(50,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.625),sd=0.002)
Cm4<-rtruncnorm(25,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.75),sd=0.0025)
Cm5<-rtruncnorm(25,a=0.0008,mean=quantile(meta(ref.traits.samp)$LMA/10,probs=0.875),sd=0.003)
Cm<-c(Cm1,Cm2,Cm3,Cm4,Cm5)

## set constant values for the other functional groups
N<-1.5
Cw<-0.007
Cab<-40
prospect.params<-data.frame(N=N,Cw=Cw,Cab=Cab,Cm=Cm,FG=FG)
prospect.params$Car<-prospect.params$Cab/6.24 ## mean chl:car ratio
prospect.params$Anth<-0
prospect.params$Cbrown<-0

## palette
colorBlind  <- c("#E69F00","#009E73","#56B4E9","#F0E442",
                 "#0072B2","#CC79A7","#D55E00","#999999")

## plot functional group LMA distributions
Cm_density<-ggplot(data=prospect.params,
                   aes(x=Cm,color=FG))+
  geom_density(size=1.5)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=15))+
  labs(y="Density",x=expression("LMA (g cm"^-2*")"),tag="(a)")+
  coord_cartesian(xlim=c(0,0.021))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

## simulate spectra with PROSPECT
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

## plot spectra per functional group
prospect_fg_spec_plot<-ggplot(prospect.spec.long,
                              aes(x=variable,y=value,color=FG))+
  stat_summary(fun=median,na.rm=T,geom="line",size=1)+
  theme_bw()+theme(text=element_text(size=15),
                   plot.margin=unit(c(0,0.2,0,0),"in"))+
  labs(tag = "(b)",x="Wavelength (nm)",y="Reflectance")+
  guides(color=guide_legend("Functional group",nrow = 2))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,0.5))+
  scale_x_continuous(expand = c(0, 0),limits=c(390,2510))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])

## generate other traits
## create T1 by adding noise and then constraining it to the same range as LMA
protoT1<-prospect.params$Cm+rnorm(200,mean=0,sd=0.003)
prospect.params$T1<-(max(prospect.params$Cm)-min(prospect.params$Cm))*(protoT1-min(protoT1))/(max(protoT1)-min(protoT1))+min(prospect.params$Cm)

## T2 is based on functional group means rather than the individual values
prospect.params$FG.mean<-unlist(apply(prospect.params,1,function(x) median(prospect.params[["Cm"]][prospect.params[["FG"]]==x[["FG"]]])))
protoT2<-prospect.params$FG.mean+rnorm(200,mean=0,sd=0.001)
prospect.params$T2<-(max(prospect.params$Cm)-min(prospect.params$Cm))*(protoT2-min(protoT2))/(max(protoT2)-min(protoT2))+min(prospect.params$Cm)

## plot relationship between synthetic traits and LMA
T1_Cm_plot<-ggplot(prospect.params,aes(y=T1,x=Cm,color=FG))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
#  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.02),ylim=c(0,0.02))+
  theme(text = element_text(size=15),
        legend.position = c(0.8, 0.2))+
  labs(tag="(c)",y="Functional Trait 1",x=expression("LMA (g cm"^-2*")"))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

T2_Cm_plot<-ggplot(prospect.params,aes(y=T2,x=Cm,color=FG))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  #  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.02),ylim=c(0,0.02))+
  theme(text = element_text(size=15),
        legend.position = c(0.8, 0.2))+
  labs(tag="(d)",y="Functional Trait 2",x=expression("LMA (g cm"^-2*")"))+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

## just to show we can predict LMA from the spectrum
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

## can we predict T1 from the spectrum?
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
  labs(tag="(e)",y="Functional Trait 1 (measured)",x="Functional Trait 1 (predicted)")+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

## can we predict T2 from the spectrum?
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
  labs(tag="(f)",y="Functional Trait 2 (measured)",x="Functional Trait 2 (predicted)")+
  scale_color_manual(values=colorBlind[c(1,2,4,5,6)])+
  guides(color=F)

pdf("GrimeReview/prospect_covar.pdf",height=12,width=11)
Cm_density+prospect_fg_spec_plot+
  T1_Cm_plot+T2_Cm_plot+
  T1_plot+T2_plot+
  plot_layout(ncol=2,guides="collect",heights=c(0.5,1,1)) &
  theme(legend.position = "bottom")
dev.off()

## just to see what the VIP looks like...
source("Scripts/VIP.R")
plot(400:2500,VIP(T2_model)[ncomp_T2_model,])
