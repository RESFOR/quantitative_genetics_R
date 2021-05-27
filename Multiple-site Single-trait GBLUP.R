# Title: Multi-site single-trait GBLUP individual-tree mixed models using breedR R-package: an example using the HT30 trait
#        and the 4 sites (JUDY, VIRG, SWAN, and TIME) for the genotyped 1,490 RES-FOR lodgepole pine trees.
#        The phenotypes were adjusted by design effects and scaled to mean 0 and variance 1. 
# Author: Eduardo Pablo Cappa (cappa.eduardo@inta.gob.ar)
# Copyright @ Dr. Eduardo Pable Cappa
# Revised: September 10th 2019"
# Revised: "November 27th, 2019" I inlcuded the WWD trait
# Revised to publish: "May 26th, 2021"

# PACKAGES
##############################################################################################
wants <- c("plyr", "fields","reshape2","psych","ggplot2","scatterplot3d",
           "breedR","viridis","lattice","pedigreemm","corrplot", "matrixcalc",
           "network","stats","data.table")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
suppressPackageStartupMessages(sapply(wants, require, character.only=TRUE))


##############################################################################################
# 1. READ THE DATA
data<-read.table(file="data.txt", head=TRUE)

# Number of the column where start the fenotyic information 
start <- 15
# Defined number of traits below analysis
# For this analysis I will keep just 5 traits (DBH30, HT30, NSWGR36, WD and MFA) from a total of 20 in the data file
ntrait = 5


###################################################################################
# 2. STRUCTURE AND VISUALIZATION OF THE DATA
# Dimension of the data
dim(data)

# Structure of the data
data[,c(start:ncol(data))] <- sapply(data[,c(start:ncol(data))], as.numeric)
data$site=as.factor(data$site);data$proc=as.factor(data$proc);data$rep=as.factor(data$rep)
data$tree=as.factor(data$tree);data$stake=as.factor(data$stake);data$set=as.factor(data$set)
data$plot=as.factor(data$plot);data$row=as.integer(data$row);data$col=as.integer(data$col)
str(data)

# Replication nested in sites
# data <- transform(data, sitrep= factor(site:rep))

# Sets nested in Replication 
data <- transform(data, repset= factor(rep:set))

# Create the plot effects
#data <- transform(data, plot= factor(as.factor(mum):rep:set))

# Data visualization
head(data, n=10);tail(data, n=10)


###############################################################################################################################
# 3. Data for Multi-site analyses (reshaping)
# In this case, I reshaped each trait in a separately file in a list.
DATA = list()
for (i in start:(ntrait+start-1)){
DATA[[i]] <-reshape(data[,c(1:6,8,i,ncol(data))], idvar = c("ORDEN","self","dad","mum","proc","rep","repset"), timevar = "site", direction = "wide")
}
# Print a few lines of the 5 studied traits
head(DATA[[15]]);tail(DATA[[15]]) # first trait HT30
head(DATA[[16]]);tail(DATA[[16]]) # DBH30
head(DATA[[17]]);tail(DATA[[17]]) # NSWGR36
head(DATA[[18]]);tail(DATA[[18]]) # Wood Density (WWD)
head(DATA[[19]]);tail(DATA[[19]]) # MFA


##########################################################################################################
# 4. PEDIGREE (A) AND GENOMIC (G) RELATIONSHIP MATRICES
# Pedigre A-matrix
ped <- build_pedigree(c('self', 'dad', 'mum'), data = data)
A <- as(getA(ped),"matrix")
Ai <- solve(A)
Ai[1:5,1:5];dim(Ai)

# Genomic G-matrix
G <- as.matrix(readRDS(file="G_matrix.RDS"))
G[1:5,1:5];dim(G)
# Order the row and column of the G-matrix equal to the dataset  using the codg column. The codg is the RES-FOR ID of each tree
G<-G[order(match(rownames(G),data$codg)),order(match(colnames(G),data$codg))]
G[1:5,1:5];dim(G)
# Check the order of the trees from the data file and the G-matrix
all.equal(as.character(data$codg),rownames(G))
all.equal(as.character(data$codg),colnames(G))
# Inverser of G-matrix
Gi <- solve(G)
Gi[1:5,1:5];dim(Gi)

# Incidence Matrix for generic effect (Breeding values of the A-matrix)
Z1 <- as(data$self, 'indMatrix')
Z <- as(Z1,"indMatrix")
dim(Z)

# Incidence Matrix for generic effect (Breeding values of the G-matrix)
Zg <- as(c(1:nrow(data)), 'indMatrix')
dim(Zg)


############################################################################################################
### 5. MULTI-SITE SINGLE-TRAIT INDIVIDUAL-TREE MIXED MODEL using genomic relationships (GBLUP model)
# Strategy: Using the Expectation-Maximization (EM) algorithm followed by one iteration with the Average Information (AI) algorithm 
#           to compute the approximated standard errors of the heritabilities and (co)variance components
# Initial genetic (Gg_ini) and residual (Rg_ini) (co)variance matrices
start_time <- Sys.time()

initial_covs = list(Gg_ini = matrix(c(0.5571,	0.3303,	0.4171,	0.2827,
                                      0.3303,	0.4084,	0.3801,	0.2693,
                                      0.4171,	0.3801,	0.5258,	0.4016,
                                      0.2827,	0.2693,	0.4016,	0.3267),4,4),
                    
                    Rg_ini = matrix(c(0.3287,	0.0000,	0.0000,	0.0000,
                                      0.0000,	0.5158,	0.0000,	0.0000,
                                      0.0000,	0.0000,	0.3999,	0.0000,
                                      0.0000,	0.0000,	0.0000,	0.5626),4,4))

# Using as example the trait HT30 in DATA[[15]]
# The GBLUP model using REML-EM algorithm
breedR.G = remlf90(fixed = cbind(DATA[[15]][,8],DATA[[15]][,9],DATA[[15]][,10],DATA[[15]][,11]) ~  - 1 + proc, 
                   generic = list(G = list(incidence = Zg, precision = Gi,var.ini = initial_covs$Gg_ini)),
                   data = DATA[[15]], 
                   method = 'em',
                   var.ini = list(residual=initial_covs$Rg_ini),
                   progsf90.options=c('use_yams','fact_once memory'))
summary(breedR.G)
# Use cov2cor() to compute genetic and residual correlations
round(cov2cor(breedR.G$var[[1]]),4)
#round(cov2cor(breedR.G$var[[2]]),4)

# Heretabilities for each site (4 in total, 4 sites)
h2_1  <- 'G_2_2_1_1/(G_2_2_1_1+R_1_1)'; h2_2  <- 'G_2_2_2_2/(G_2_2_2_2+R_2_2)';       
h2_3  <- 'G_2_2_3_3/(G_2_2_3_3+R_3_3)'; h2_4  <- 'G_2_2_4_4/(G_2_2_4_4+R_4_4)'
# Genetic correlations
r12<-'G_2_2_1_2/(G_2_2_1_1*G_2_2_2_2)**0.5'; r13<-'G_2_2_1_3/(G_2_2_1_1*G_2_2_3_3)**0.5' ;r14<-'G_2_2_1_4/(G_2_2_1_1*G_2_2_4_4)**0.5'
r23 <- 'G_2_2_2_3/(G_2_2_2_2*G_2_2_3_3)**0.5'; r24 <- 'G_2_2_2_4/(G_2_2_2_2*G_2_2_4_4)**0.5'; r34 <- 'G_2_2_3_4/(G_2_2_3_3*G_2_2_4_4)**0.5'

# The GBLUP model using REML-EM algorithm
breedR.G.ai = remlf90(fixed = cbind(DATA[[15]][,8],DATA[[15]][,9],DATA[[15]][,10],DATA[[15]][,11]) ~  - 1 + proc, 
                      generic = list(G = list(incidence = Zg, precision = Gi,var.ini = breedR.G$var$G)),
                      data = DATA[[15]], 
                      method = 'ai',
                      var.ini = list(residual=breedR.G$var$Residual),
                      progsf90.options=c('use_yams','fact_once memory','maxrounds 1',
                                         paste('se_covar_function h2_1', h2_1), paste('se_covar_function h2_2', h2_2),
                                         paste('se_covar_function h2_3', h2_3), paste('se_covar_function h2_4', h2_4),
                                         paste('se_covar_function r12', r12),paste('se_covar_function r13', r13),
                                         paste('se_covar_function r14', r14),paste('se_covar_function r23', r23),
                                         paste('se_covar_function r24', r24),paste('se_covar_function r34', r34)))
summary(breedR.G.ai)
  
# Save the ABLUP and GBLUP results
saveRDS(list(GBLUP.em=breedR.G, GBLUP.ai=breedR.G.ai),paste("Multiple-site Single-trait ABLUP and GBLUP",".RDS",sep=""))

end_time <- Sys.time()
print(end_time - start_time)

##############################################################################################
# 6.  Breeding values (BVs) for each site 
# BLUP of BVs and its standard errors (s.e.) for parents and progenies and for each trait
BV.G<-breedR.G$ranef$G

BV1.G<-cbind(data$self, data$site, BV.G[[1]], BV.G[[2]], BV.G[[3]], BV.G[[4]])
colnames(BV1.G)[1] = "self"
colnames(BV1.G)[2] = "site"
# BVs of the tree from the site where the tree come
BV2.G.s1 = subset(BV1.G, site == 1)[,c(1:2,3:4)]; BV2.G.s2 = subset(BV1.G, site == 2)[,c(1:2,5:6)]
BV2.G.s3 = subset(BV1.G, site == 3)[,c(1:2,7:8)]; BV2.G.s4 = subset(BV1.G, site == 4)[,c(1:2,9:10)]
names(BV2.G.s1) <- c("self","site","BV_G.HT30","Acc_G.HT30"); names(BV2.G.s2) <- c("self","site","BV_G.HT30","Acc_G.HT30")
names(BV2.G.s3) <- c("self","site","BV_G.HT30","Acc_G.HT30"); names(BV2.G.s4) <- c("self","site","BV_G.HT30","Acc_G.HT30")
GBLUP_Of = rbind(BV2.G.s1, BV2.G.s2, BV2.G.s3, BV2.G.s4)

head(GBLUP_Of);tail(GBLUP_Of)

write.table(GBLUP_Of, file = "GBLUP_BVs_SE_HT30.txt", row.names = FALSE)
