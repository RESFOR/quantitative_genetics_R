# Title: Standard, Spatial and competition individual-tree mixed models / Spatial and Competition diagnosis using breedR R-package:
#        an example using total height trait at 30 age for the RES-FOR Red Earth (REDE) site.
# Author: Eduardo Pablo Cappa
# Date: "February 21th, 2019"
# Revised: "April 10th, 2019"
# Revised to publish: "May 25th, 2021"

# PACKAGES
##############################################################################################
wants <- c("plyr", "fields","reshape2","psych","ggplot2","scatterplot3d",
           "breedR","viridis","lattice","pedigreemm")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
suppressPackageStartupMessages(sapply(wants, require, character.only=TRUE))

# THE breedR PACKAGE

# The breedR package on the web
# http://famuvie.github.io/breedR/

# breedR tutorials
# https://github.com/famuvie/breedR/wiki/Overview

# breedR development site 
# https://github.com/famuvie/breedR

# Joint training workshop on breedR and phenotypic plasticity June 30 – July 3 2015 Jaca, Spain
# http://famuvie.github.io/breedR/workshop/
# Some slides about Genetic models, spatial autocorrelation and competition
# http://famuvie.github.io/breedR/workshop/day1/ECappa_theory.pdf

# Some papers related to spatial and/or competiton models using breedR package and diagnosis tools
# 1. https://link.springer.com/article/10.1007/s13595-019-0836-9
# 2. https://link.springer.com/article/10.1007/s11056-018-9682-0
# 3. https://academic.oup.com/forestscience/article/65/5/570/5556922
# 4. https://link.springer.com/article/10.1007/s11295-016-1061-4


##############################################################################################
# 1. READ THE DATA
data.s3<-read.table(file="data.s3.txt", head=TRUE)

# Number of the column where start the phenotypic information 
start <- 14


##############################################################################################
# 2. STRUCTURE AND VISUALIZATION OF THE DATA
# Dimension of the data
dim(data.s3)

# Structure of the data
str(data.s3)
data.s3[,c(start:ncol(data.s3))] <- sapply(data.s3[,c(start:ncol(data.s3))], as.numeric)
data.s3$site=as.factor(data.s3$site);data.s3$proc=as.factor(data.s3$proc);data.s3$rep=as.factor(data.s3$rep)
data.s3$tier=as.factor(data.s3$tier);data.s3$tree=as.factor(data.s3$tree)
str(data.s3)

# Create the plot effects
#data.s3 <- transform(data.s3, plot= factor(as.factor(mum):rep))

# Data visualization
head(data.s3, n=10);tail(data.s3, n=10)


##########################################################################################################
# 3. CHOOSE THE TARGET TRAIT
data.s3$phe<-data.s3$HT30/100

##########################################################################################################
# 4.STANDARD ANALISIS FAMILY MODEL: The family additive model with design effects (replicates and plots)
res.std_flia <- remlf90(fixed = phe ~ proc,
                       random = ~ rep + plot + factor(mum),
                       data = data.s3, 
                       method = 'em')
summary(res.std_flia)

#Narrow-sense individual-tree heritability (4*Var_flia / (Var_Flia + Var_ Residual))
round(with(res.std_flia, 4*var["factor(mum)",1] / sum(var["factor(mum)",1]+var["Residual",1])),4)


##############################################################################################
# 5. STANDARD ANALYSIS INDIVIDUAL-TREE MIXED MDOEL: The additive individual-tree model with design effects (replicates and plots)
# Pedigre A-matrix
ped <- build_pedigree(c('self', 'dad', 'mum'), data = data.s3)
A <- as(getA(ped),"matrix")
Ai <- solve(A)
dim(A)

# Incidence Matrix for generic effect (Breeding values)
Z1 <- as(data.s3$self, 'indMatrix')
Z <- as(Z1,"indMatrix")
dim(Z)

# The model
res.blk <- remlf90(fixed = phe ~ proc,
                   random= ~ plot,
                   generic = list(A = list(incidence = Z, precision = Ai)),
                   spatial = list(model = 'blocks', 
                                  coord = data.s3[, c('row','col')], 
                                  id = "rep"),
                   data = data.s3, 
                   method = 'em')
summary(res.blk)

##############################################################################################
# 5.1. STANDARD ANALYSIS INDIVIDUAL-TREE MIXED MDOEL without design effects (replicates and plots)
res.stdSB <- remlf90(fixed = phe ~ proc, 
                     genetic = list(model = 'add_animal', 
                                    pedigree = data.s3[, c('self','dad','mum')], 
                                    id = 'self'), 
                     data = data.s3, 
                     method = 'em')
summary(res.stdSB)
# Residuals
Residuals_stdSB<-residuals(res.stdSB)


##############################################################################################
# 6. DIAGNOSIS
# 6.1. Plots the phenotypes and residuals by rows and columns (2-D plots)
# Using breedR package
coordinates(res.stdSB) <- data.s3[, c('row','col')]
# If you would like change the default colors
breedR.setOption(col.seq = c('yellow', 'red'))
#breedR.setOption(col.div = c('yellow', 'red'))
# Phenotype and Residuals plots
plot(res.stdSB, type = 'phenotype')
plot(res.stdSB, type = 'residuals')

# Using virdis package
vp = plot(res.stdSB, 'phenotype') 
vp + scale_fill_viridis(name = 'Phenotype',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns')

vr = plot(res.stdSB, 'residuals') 
vr + scale_fill_viridis(name = 'Residuals',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns')

# 6.2.Variograms of residuals
(variogram_stdSB<-variogram(res.stdSB))
variogram(res.stdSB, plot = 'isotropic')
variogram(res.stdSB, plot = 'perspective')
variogram(res.stdSB, plot = 'perspective', R=15)

# 6.3.Residuals against row and column position
xyplot(Residuals_stdSB ~ data.s3$col | data.s3$row,as.table=2,strip = strip.custom(strip.names = TRUE, strip.levels = TRUE),
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.loess(x,y, span= 0.6,degree = 2)
       })

xyplot(Residuals_stdSB ~ data.s3$row | data.s3$col,as.table=2,strip = strip.custom(strip.names = TRUE, strip.levels = TRUE),
       panel = function(x, y) {
         panel.xyplot(x, y)
         panel.loess(x,y, span= 0.6,degree = 2)
       })

####################################################################################
# 7. SPATIAL ANALISIS using individual-tree mixed model with two-dimensional B-spline: The B-Spline model
# BE CAREFULLY!!!!!!!!!!!
# The spatial parameters n.konts for row and column should be optimized for each particular data set
# See tuning spatial parameters here: https://github.com/famuvie/breedR/wiki/Overview#exercise--tuning-spatial-parameters

res.spl <- remlf90(fixed = phe ~ proc, 
                   random = ~ rep + plot,
                   generic = list(A = list(incidence = Z, precision = Ai)),
                   spatial = list(model = 'splines', 
                                  coord = data.s3[, c('row','col')], 
                                  n.knots = c(7, 7)),
                   data = data.s3, 
                   method = 'em')
summary(res.spl)


##############################################################################################
# 8. SPATIAL ANALISIS using individual-tree mixed model separable kronecker product of First order 
# Autoregressive processes on the rows and the colums: The autoregressive model
# BE CAREFULLY!!!!!!!!!!!
# The spatial parameters rho for row and column should be optimized or each particular data set
# See tuning spatial parameters here: https://github.com/famuvie/breedR/wiki/Overview#exercise--tuning-spatial-parameters

res.ar <- remlf90(fixed = phe ~ proc,
                   random = ~ rep + plot,
                   generic = list(A = list(incidence = Z, precision = Ai)),
                   spatial = list(model = 'AR', 
                                  coord = data.s3[, c('row','col')],
                                  rho = c(0.924982,0.956806)),
                   data = data.s3, 
                   method = 'em')
summary(res.ar)


##############################################################################################
# 9. Plots comparing the diferent spatial models - we preserve the scale by using compare.plots() -
# 9.1. Comparison of variograms: isotropic variograms to be compared
variogram_blk<-variogram(res.blk)
variogram_spl<-variogram(res.spl)
variogram_ar<-variogram(res.ar)

isotropic_variograms <- rbind(cbind(model = "res.stdSB", variogram_stdSB[["isotropic"]]),
                              cbind(model = "res.blk", variogram_blk[["isotropic"]]),
                              cbind(model = "res.spl", variogram_spl[["isotropic"]]),
                              cbind(model = "res.ar", variogram_ar[["isotropic"]]))

# Then, you can plot the variograms together under the same scale
ggplot(isotropic_variograms, aes(distance, variogram)) + geom_point() + geom_line() +
  stat_smooth(se = FALSE, method = 'auto') + facet_wrap(~ model)

# 9.2. Comparison of residuals
compare.plots(
  list(`Individual-tree model only` = plot(res.stdSB, 'residuals') + scale_fill_viridis(name = 'Residuals',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'),
       `Individual-tree/blocks model` = plot(res.blk, 'residuals') + scale_fill_viridis(name = 'Residuals',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'),
       `Individual-tree/splines model` = plot(res.spl, 'residuals')+scale_fill_viridis(name = 'Residuals',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'),
       `Individual-tree/AR model` = plot(res.ar, 'residuals')      +scale_fill_viridis(name = 'Residuals',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'))
)

# 9.3. Comparison of spatial components
compare.plots(
  list(Blocks = plot(res.blk, type = 'spatial')  + scale_fill_viridis(name = 'Spatial',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'),
       Splines = plot(res.spl, type = 'spatial') +scale_fill_viridis(name = 'Spatial',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'),
       AR1xAR1 = plot(res.ar, type = 'spatial')  +scale_fill_viridis(name = 'Spatial',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'))
)

# 9.4. Prediction of the spatial effect in unobserved locations. 
# The type fullspatial fills the holes (when possible)
compare.plots(
  list(Blocks = plot(res.blk, type = 'fullspatial')  + scale_fill_viridis(name = 'Full_Spatial',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'),
       Splines = plot(res.spl, type = 'fullspatial') +scale_fill_viridis(name = 'Full_Spatial',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'),
       AR1xAR1 = plot(res.ar, type = 'fullspatial')  +scale_fill_viridis(name = 'Full_Spatial',  option = "B") + theme_bw() + xlab('Rows') + ylab('Columns'))
)


##############################################################################################
# 10. Spearman correlations and plots of breeding values (BVs) for parents and offspring between standard and spatial analyses.
# 10.1. Number of parents and offspring
(nparents<-c(max(data.s3$mum,data.s3$dad)))
(nprogenies<-c(max(data.s3$self)) - nparents)

# 10.2. BVs
BVs_blk_parents  <-res.blk$ranef$A[[1]]$value[1:nparents]
BVs_blk_progenies<-res.blk$ranef$A[[1]]$value[nparents+1:nprogenies]
BVs_spl_parents  <-res.spl$ranef$A[[1]]$value[1:nparents]
BVs_spl_progenies<-res.spl$ranef$A[[1]]$value[nparents+1:nprogenies]
BVs_ar_parents   <-res.ar$ranef$A[[1]]$value[1:nparents]
BVs_ar_progenies <-res.ar$ranef$A[[1]]$value[nparents+1:nprogenies]

# 10.3. Spearman correlations of BVs between the different models
a=cor(BVs_blk_parents,BVs_ar_parents,method = c("spearman"))
b=cor(BVs_blk_parents,BVs_spl_parents,method = c("spearman"))
c=cor(BVs_spl_parents,BVs_ar_parents,method = c("spearman"))
d=cor(BVs_blk_progenies,BVs_ar_progenies,method = c("spearman"))
e=cor(BVs_blk_progenies,BVs_spl_progenies,method = c("spearman"))
f=cor(BVs_spl_progenies,BVs_ar_progenies,method = c("spearman"))

All_cor=round(c(a,b,c,d,e,f),3)

Cor_Final_Results <- matrix(All_cor,nrow=2,byrow=TRUE)
rownames (Cor_Final_Results) <- c("Parents","Progenies")
colnames (Cor_Final_Results) <- c("BLK-AR","BLK-SPL","SPL-AR")
(Cor_Final_Results <- as.table(Cor_Final_Results))

# 10.4. Plots of BVs between the different standard (BLK) and spatial (SPL and AR) models
par(mfrow = c(2, 3))
plot(BVs_blk_parents,BVs_ar_parents)
abline(lm(BVs_ar_parents~BVs_blk_parents), col="red")
plot(BVs_blk_parents,BVs_spl_parents)
abline(lm(BVs_spl_parents~BVs_blk_parents), col="red")
plot(BVs_spl_parents,BVs_ar_parents)
abline(lm(BVs_ar_parents~BVs_spl_parents), col="red")
plot(BVs_blk_progenies,BVs_ar_progenies)
abline(lm(BVs_ar_progenies~BVs_blk_progenies), col="red")
plot(BVs_blk_progenies,BVs_spl_progenies)
abline(lm(BVs_spl_progenies~BVs_blk_progenies), col="red")
plot(BVs_spl_progenies,BVs_ar_progenies)
abline(lm(BVs_ar_progenies~BVs_spl_progenies), col="red")


##############################################################################################
# 11. THEORETICAL ACCURACIES OF BVs
# 11.1. SEs
SEs_blk_parents  <-res.blk$ranef$A[[1]]$s.e.[1:nparents]
SEs_blk_progenies<-res.blk$ranef$A[[1]]$s.e.[nparents+1:nprogenies]
SEs_spl_parents  <-res.spl$ranef$A[[1]]$s.e.[1:nparents]
SEs_spl_progenies<-res.spl$ranef$A[[1]]$s.e.[nparents+1:nprogenies]
SEs_ar_parents   <-res.ar$ranef$A[[1]]$s.e.[1:nparents]
SEs_ar_progenies <-res.ar$ranef$A[[1]]$s.e.[nparents+1:nprogenies]

# 11.2  PEVs: prediction error variance
SE2s_blk_parents  =SEs_blk_parents *SEs_blk_parents 
SE2s_blk_progenies=SEs_blk_progenies*SEs_blk_progenies
SE2s_spl_parents  =SEs_spl_parents *SEs_spl_parents 
SE2s_spl_progenies=SEs_spl_progenies*SEs_spl_progenies
SE2s_ar_parents   =SEs_ar_parents  *SEs_ar_parents  
SE2s_ar_progenies =SEs_ar_progenies*SEs_ar_progenies

# 11.3 Fi
Fi   = diag(A) - 1
Fipa =Fi[1:nparents]
Fipr =Fi[nparents+1:nprogenies]

# 11.4  Accuracies (ACCs)
ACC_blk_parents  =sqrt(1-(SE2s_blk_parents  /((1+Fipa)*with(res.blk, var["A",1]))))
ACC_blk_progenies=sqrt(1-(SE2s_blk_progenies/((1+Fipr)*with(res.blk, var["A",1]))))
ACC_spl_parents  =sqrt(1-(SE2s_spl_parents  /((1+Fipa)*with(res.spl, var["A",1]))))
ACC_spl_progenies=sqrt(1-(SE2s_spl_progenies/((1+Fipr)*with(res.spl, var["A",1]))))
ACC_ar_parents   =sqrt(1-(SE2s_ar_parents   /((1+Fipa)*with(res.ar, var["A",1]) )))
ACC_ar_progenies =sqrt(1-(SE2s_ar_progenies /((1+Fipr)*with(res.ar, var["A",1]) )))

All_acc=round(c(mean(ACC_blk_parents,na.rm = TRUE),
                mean(ACC_spl_parents,na.rm = TRUE),
                mean(ACC_ar_parents,na.rm = TRUE), 
                mean(ACC_blk_progenies,na.rm = TRUE),
                mean(ACC_spl_progenies,na.rm = TRUE),
                mean(ACC_ar_progenies,na.rm = TRUE)),3)

Acc_Final_Results <- matrix(All_acc,nrow=2,byrow=TRUE)
rownames (Acc_Final_Results) <- c("Parents","Progenies")
colnames (Acc_Final_Results) <- c("BLK","SPL","AR")
(Acc_Final_Results <- as.table(Acc_Final_Results))

# 11.5. Plots of ACCs between the different models
par(mfrow = c(2, 3))
plot(ACC_blk_parents,  ACC_ar_parents)
abline(lm(ACC_ar_parents~ACC_blk_parents), col="red")
plot(ACC_blk_parents,  ACC_spl_parents)
abline(lm(ACC_spl_parents~ACC_blk_parents), col="red")
plot(ACC_spl_parents,  ACC_ar_parents)
abline(lm(ACC_ar_parents~ACC_spl_parents), col="red")
plot(ACC_blk_progenies,ACC_ar_progenies)
abline(lm(ACC_ar_progenies~ACC_blk_progenies), col="red")
plot(ACC_blk_progenies,ACC_spl_progenies)
abline(lm(ACC_spl_progenies~ACC_blk_progenies), col="red")
plot(ACC_spl_progenies,ACC_ar_progenies)
abline(lm(ACC_ar_progenies~ACC_spl_progenies), col="red")

# 11.6 OUTPUT files BVs and ACCs
BLUP_ACC_Pa = cbind(BVs_blk_parents,ACC_blk_parents,BVs_spl_parents,ACC_spl_parents,BVs_ar_parents,ACC_ar_parents)
BLUP_ACC_Of = cbind(BVs_blk_progenies,ACC_blk_progenies,BVs_spl_progenies,ACC_spl_progenies,BVs_ar_progenies,ACC_ar_progenies)  

write.table(BLUP_ACC_Pa, file = "BVs_ACCs_Parent.txt")
write.table(BLUP_ACC_Of, file = "BVs_ACCs_Offspring.txt")


######################################################################################################
# 12. COMPETITION ANALISIS
# BE CAREFULLY!!!!!!!!!!!
# The spatial parameters rho should be optimized
# See tuning spatial parameters here: https://github.com/famuvie/breedR/wiki/Overview#exercise--tuning-spatial-parameters

res.comp <- remlf90(fixed = phe ~ proc,
                    random = ~ rep + plot,
                    genetic = list(model = c('comp'),
                                   pedigree = data.s3[, c('self','dad','mum')],
                                   id = 'self',
                                   coord = data.s3[, c('row', 'col')],
                                   competition_decay = 1, # IC decay 1/distance
                                   pec = list(present = TRUE)), #envirmonetal compettion effect
                    spatial = list(model = 'AR', 
                                   coord = data.s3[, c('row','col')],
                                   rho = c(0.924982,0.956806)),
                    data = data.s3,
                    method = 'em'
                    #debug = F
                    )
summary(res.comp)
# Direct and competition additive correlation
(Correlation.com = cov2cor(res.comp$var$genetic))


###################################################################################################
# 13. VARIOGRAMS FOR ALL MODELS
variogram(res.stdSB)
variogram(res.blk)
variogram(res.spl)
variogram(res.ar)
variogram(res.comp)

###################################################################################################
# 14. ADDITONAL DIAGNOSIS of competition effects:
# Plot of phenotypic and residual values after fitting genetic effects plotted against means 
# of the 8 nearest neighbour trees. Means are weight by IC factor.
# Firs at all I remove the trees with NA from the phenoipic data and the matrix Zc
data.s3na=subset(data.s3,!is.na(phe))
# Desing additive genetic competition matrix
Zc <- as.matrix(model.matrix(res.comp)$"genetic_competition")[1:length(data.s3$phe),(nparents+1):(nparents+nprogenies)]
rownames(Zc)=colnames(Zc)=data.s3$self
# Desing additive genetic competition matrix without NA phenotypic values
Zcna =Zc[which(rownames(Zc)%in%data.s3na$self),which(rownames(Zc)%in%data.s3na$self)]

# Phenotype values
phe_X<-as.matrix(data.s3na$phe)
# Phenotypic means of the 8 nearest neighbour trees
Nei_mean_pheno<- Zcna%*%phe_X
# Residual means of the 8 nearest neighbour trees
resid=as.matrix(Residuals_stdSB)
rownames(resid)=data.s3$self
resid.na=subset(resid,!is.na(resid[,1]))
Nei_mean_res<- Zcna%*%resid.na

# 14.1 Plot of the phenotypic values of each tree against the phenotypic means of the neighbouring trees.
# Means are weight by IC factor.
par(mfrow = c(1, 2))
plot(Nei_mean_pheno,phe_X)
abline(lm(phe_X~Nei_mean_pheno), col="red")
cor_pheno = cor(phe_X,Nei_mean_pheno)
cor_pheno = round(cor_pheno,3)
text(10,13,paste("r =", cor_pheno))

# 14. 2. Plot of the residual values of each tree after fitting genetic effects against the phenotypic means of the neighbouring trees.
#Means are weight by IC factor.
plot(Nei_mean_pheno,resid.na)
abline(lm(resid.na~Nei_mean_pheno), col="blue")
cor_resid = cor(resid.na,Nei_mean_pheno)
cor_resid = round(cor_resid,3)
text(10,3,paste("r =", cor_resid))


############################################################################################################
# 15. TABLE WITH THE RESULTS OF THE FITTED SPATIAL MODELS: LogL, AIC, VAriance components and heritabili
# 15.1 Fit of each model: AIC
AIC.blk<-round(res.blk$fit$AIC,0)
AIC.spl<-round(res.spl$fit$AIC,0)
AIC.ar<-round(res.ar$fit$AIC,0)
AIC.comp<-round(res.comp$fit$AIC,0)

# 14.2 Variance components
Genetic.blk     <-round(with(res.blk, var["A",1]),2)
Replication.blk <-round(with(res.blk, var["spatial",1]),2)
Plot.blk <-round(with(res.blk, var["plot",1]),2)
Residual.blk    <-round(with(res.blk, var["Residual",1]),2)

Replication.spl <-round(with(res.spl, var["rep",1]),2)
Plot.spl <-round(with(res.spl, var["plot",1]),2)
Genetic.spl     <-round(with(res.spl, var["A",1]),2)
Spatial.spl     <-round(with(res.spl, var["spatial",1]),2)
Residual.spl    <-round(with(res.spl, var["Residual",1]),2)

Replication.ar <-round(with(res.ar, var["rep",1]),2)
Plot.ar <-round(with(res.ar, var["plot",1]),2)
Genetic.ar     <-round(with(res.ar, var["A",1]),2)
Spatial.ar     <-round(with(res.ar, var["spatial",1]),2)
Residual.ar    <-round(with(res.ar, var["Residual",1]),2)

Replication.comp<-round(with(res.comp, var$"rep"[1,1]),2)
Plot.comp<-round(with(res.comp, var$"plot"[1,1]),2)
Genetic.comp    <-round(res.comp$var$genetic[1,1],2)
Competition.comp<-round(res.comp$var$genetic[2,2],2)
Correlation.comp<-round(Correlation.com[1,2],2)
Spatial.comp    <-round(with(res.comp, var$"spatial"[1,1]),2)
Env.comp        <-round(with(res.comp, var$"pec"[1,1]),2)
Residual.comp   <-round(with(res.comp, var$"Residual"[1,1]),2)

# 15.3 Narrow-sense individual-tree heritabilities for each model
(h2N_blk<- with(res.blk, var["A",1] / sum(var["A",1]+var["Residual",1])))
(h2N_spl<- with(res.spl, var["A",1] / sum(var["A",1]+var["Residual",1])))
(h2N_ar<- with(res.ar, var["A",1]   / sum(var["A",1]+var["Residual",1])))

All_par<-c(
  AIC.blk,AIC.spl,AIC.ar,AIC.comp,
  Genetic.blk,Genetic.spl,Genetic.ar,Genetic.comp,
  "-","-","-",Competition.comp,
  "-","-","-",Correlation.comp,
  Replication.blk,Replication.spl,Replication.ar,Replication.comp,
  Plot.blk,Plot.spl,Plot.ar,Plot.comp,
  "-",Spatial.spl,Spatial.ar,Spatial.comp,
  "-","-","-",Env.comp,
  Residual.blk,Residual.spl,Residual.ar,Residual.comp,
  round(h2N_blk,2),round(h2N_spl,2),round(h2N_ar,2),"-"
)

Par_Final_Results <- matrix(All_par,ncol=4,byrow=TRUE)
rownames (Par_Final_Results) <- c("AIC","Genetic","Competition","Correlation","Replication","Plot","Spatial","Env.Comp","Residual", "Heritability")
colnames (Par_Final_Results) <- c("STR","SPL","AR1","COMP.AR")
(Par_Final_Results <- as.table(Par_Final_Results))


############################################################################################################
# 16. Estimated provenance effects (BLUEs) for each model
BLUEs_proc = cbind(res.blk$fixed$proc[[1]]$value, 
                   res.spl$fixed$proc[[1]]$value, 
                   res.ar$fixed$proc[[1]]$value,
                   res.comp$fixed$proc[[1]]$value)
colnames(BLUEs_proc)=c("STR","SPL","AR1","COMP.AR")
rownames(BLUEs_proc)=rownames(res.blk$fixed$proc[[1]])
BLUEs_proc=as.table((BLUEs_proc))
round(BLUEs_proc,1)

par(mfrow = c(1, 1))
barplot(t(BLUEs_proc), 
        beside=TRUE, 
        legend = colnames(BLUEs_proc),
        main="Estimated Provenance effects (BLUEs) by model",
        args.legend = list(x = "topright", bty = "n", inset=c(0.06, 0.06))
        )


############################################################################################################
# 17. ADJUSTED DESIGN AND/OR SPATIAL PHENOTYPES
# 17.1. Design and Spatial phenotype adjustment (by repplicates, plots, and residual spatial -AR model- effects)
# Replicate effects from AR spatial model
rep.phe <- model.matrix(res.ar)$rep %*% ranef(res.ar)$rep
# Plot effects from AR spatial model
plot.phe <- model.matrix(res.ar)$plot %*% ranef(res.ar)$plot
# Residual spatial effects from AR spatial model
sp.phe <- model.matrix(res.ar)$spatial %*%  ranef(res.ar)$spatial
# Design and Spatial phenotype adjusted phenotype
data.s3$phe_Spatial_adj = data.s3$phe - as.vector(rep.phe) - as.vector(plot.phe) - as.vector(sp.phe)

# 17.2. Design adjustment (by replicates and plot effects) from the Block Model phenotype for the target traits
# Replicate effects from standard model
repb.phe <- model.matrix(res.blk)$spatial %*% ranef(res.blk)$spatial
# Plot effects from standard model
plotb.phe <- model.matrix(res.blk)$plot %*% ranef(res.blk)$plot
# Design phenotype adjustment 
data.s3$phe_Design_adj = data.s3$phe - as.vector(repb.phe) - as.vector(plotb.phe)

# PLOTs Raw vs. Spatial and Design adjusted phenotype (colored by replication)
par(mfrow = c(1, 2))
P.phe = plot(data.s3$phe_Spatial_adj, data.s3$phe, xlab  = "Spatially adjusted phenotype" , ylab = "Raw Phenotype", main = "Spatial Adjusted vs. Raw phenotype", col=data.s3$rep)
legend("topleft", bty ="n", legend=levels(factor(data.s3$rep)), text.col=seq_along(levels(factor(data.s3$rep))))

P.phe = plot(data.s3$phe_Design_adj, data.s3$phe, xlab  = "Design adjusted phenotype" , ylab = "Raw Phenotype", main = "Design adjusted vs. Raw phenotype", col=data.s3$rep)
legend("topleft", bty ="n", legend=levels(factor(data.s3$rep)), text.col=seq_along(levels(factor(data.s3$rep))))

# OUTPUT FILE with Spatiala and Design, and Design adjusted phenotype
write.table(cbind(data.s3$self, data.s3$phe, data.s3$phe_Spatial_adj,data.s3$phe_Design_adj), file = "Adjusted Phenotypes.txt", row.names = FALSE, col.names = c("Self", "Raw_Pheno", "Spatial_Adj", "Design_Adj"))
