# Title: Multi-trait multi-site ABLUP and HBLUP individual-tree mixed models using breedR R-package: an example using three traits (DBH30, HT30, NSWGR36) 
#        and 4 sites (JUDY, VIRG, SWAN, and TIME) from the RES-FOR lodgepole pine dataset.
# Author: Eduardo Pablo Cappa
# Date: "Abril 29th, 2019"
# Revised: "October 09th, 2019"  Cheking code
# Revised to publish: "May 25th, 2021"

# PACKAGES
##############################################################################################
wants <- c("plyr", "fields","reshape2","psych","ggplot2","scatterplot3d",
           "breedR","viridis","lattice","pedigreemm","data.table","matrixcalc")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
suppressPackageStartupMessages(sapply(wants, require, character.only=TRUE))

##############################################################################################
# 1. READ THE DATA
data<-read.table(file="data.txt", head=TRUE)

# Number of the column where start the fenotyic information 
start <- 16

###################################################################################


# 2. STRUCTURE AND VISUALIZATION OF THE DATA
# Dimension of the data
dim(data)

# Structure of the data
#data[,c(5:(start-1))] <- sapply(data[,c(5:(start-1))], as.factor)
data[,c(start:ncol(data))] <- sapply(data[,c(start:ncol(data))], as.numeric)
data$site=as.factor(data$site);data$proc=as.factor(data$proc);data$rep=as.factor(data$rep)
data$set=as.factor(data$set);data$tree=as.factor(data$tree);data$stake=as.factor(data$stake)
data$test=as.factor(data$test);data$plot=as.factor(data$plot)
data$row=as.integer(data$row);data$col=as.integer(data$col)
str(data)

# Replication nested in sites
data <- transform(data, sitrep= factor(site:rep))

# Sets nested in Replication 
data <- transform(data, repset= factor(rep:set))

# Create the plot effects
#data <- transform(data, plot= factor(as.factor(mum):rep:set))

# Data visualization
head(data, n=10);tail(data, n=10)


##########################################################################################################
# 3. Data for Multi-trait Multiple-site analyses (reshaping)
# Spatial adjusted traits below analysis: DBH30c, HT30c, NSWGR3604c
dat = as.data.table(data)
head(dat)

data = dcast(setDT(dat), ORDEN+self+dad+mum+site+proc+codg+rep+set+test+row+col+stake+tree+plot  ~ site, value.var = c("HT30c","DBH30c","NSWGR3604c"))
data[1:5,];data[(nrow(data)-10):nrow(data),];dim(data)

data = as.data.frame(data)
class(data);head(data);dim(data)


##########################################################################################################
# 4. PEDIGREE (A) AND COMBINED PEDIGREE-GENOMIC (H) RELATIONSHIP MATRICES
# Pedigre A-matrix
ped <- build_pedigree(c('self', 'dad', 'mum'), data = data)
A <- as(getA(ped),"matrix");dim(A)
Ai <- solve(A)
Ai[1:5,1:5];dim(Ai)

# (Inverse) Combined pedigree-genomic H-matrix from blupf90 family programs:
# http://nce.ads.uga.edu/wiki/doku.php?id=readme.pregsf90
Hi <- as.matrix(readRDS(file="Hinv8K.RDS"))
Hi[1:5,1:5];dim(Hi)

# Checking trees in A- and H-matrices
all.equal(rownames(Ai),rownames(Hi))

# Incidence Matrix for generic effect (Breeding values)
Z1 <- as(data$self, 'indMatrix')
Z <- as(Z1,"indMatrix")
dim(Z)

# Check the Hi matrix for positive definite (PD) 
eigen <- eigen(Hi)
(eigen <- head(sort(eigen$values)))  # values must to be > 0

############################################################################################################
### 5. STANDARD MULTI-TRAIT MULTI-SITE INDIVIDUAL-TREE MIXED MODEL usnig pedigree relationships (ABLUP model)
# Strategy: using the Expectation-Maximization (EM) algorithm followed by one iteration with the Average Information (AI) algorithm 
#           to compute the approximated standard errors of the heritabilities and (co)variance components
# Initial genetic (GA_ini) and residual (RA_ini) (co)variance matrices
Ga_ini = matrix(c(1.916,	1.222,	1.16,	1.257,	0.5195,	0.246,	0.1853,	0.3315,	-0.3071,	-0.1946,	-0.3549,	-0.2493,
                  1.222,	2.121,	1.448,	1.32,	0.2754,	0.457,	0.3921,	0.4061,	-0.2173,	-0.269,	-0.2236,	-0.2203,
                  1.16,	1.448,	1.464,	1.235,	0.3012,	0.3437,	0.4115,	0.3812,	-0.2909,	-0.3196,	-0.2991,	-0.2739,
                  1.257,	1.32,	1.235,	1.504,	0.1975,	0.1612,	0.1782,	0.3655,	-0.3391,	-0.3222,	-0.3415,	-0.2677,
                  0.5195,	0.2754,	0.3012,	0.1975,	0.3325,	0.2472,	0.204,	0.2099,	-0.1244,	-0.02782,	-0.1125,	-0.08652,
                  0.246,	0.457,	0.3437,	0.1612,	0.2472,	0.436,	0.2643,	0.2269,	-0.04835,	-0.0109,	-0.02192,	-0.02774,
                  0.1853,	0.3921,	0.4115,	0.1782,	0.204,	0.2643,	0.3102,	0.2361,	-0.1401,	-0.1147,	-0.124,	-0.1275,
                  0.3315,	0.4061,	0.3812,	0.3655,	0.2099,	0.2269,	0.2361,	0.316,	-0.1904,	-0.1254,	-0.1734,	-0.1585,
                  -0.3071,	-0.2173,	-0.2909,	-0.3391,	-0.1244,	-0.04835,	-0.1401,	-0.1904,	0.3396,	0.3116,	0.3491,	0.3076,
                  -0.1946,	-0.269,	-0.3196,	-0.3222,	-0.02782,	-0.0109,	-0.1147,	-0.1254,	0.3116,	0.3907,	0.3663,	0.3424,
                  -0.3549,	-0.2236,	-0.2991,	-0.3415,	-0.1125,	-0.02192,	-0.124,	-0.1734,	0.3491,	0.3663,	0.416,	0.3635,
                  -0.2493,	-0.2203,	-0.2739,	-0.2677,	-0.08652,	-0.02774,	-0.1275,	-0.1585,	0.3076,	0.3424,	0.3635,	0.3521),12,12)
# Non-zero covariance betweeen traits whithin sites and zero covariance between sites. 
Ra_ini = matrix(c(6.487,	0,	0,	0,	1.804,	0,	0,	0,	-0.01132,	0,	0,	0,
                  0,	5.586,	0,	0,	0,	1.662,	0,	0,	0,	-0.1669,	0,	0,
                  0,	0,	4.338,	0,	0,	0,	0.8702,	0,	0,	0,	-0.1768,	0,
                  0,	0,	0,	5.285,	0,	0,	0,	1.268,	0,	0,	0,	-0.179,
                  1.804,	0,	0,	0,	1.061,	0,	0,	0,	0.02311,	0,	0,	0,
                  0,	1.662,	0,	0,	0,	1.103,	0,	0,	0,	-0.08942,	0,	0,
                  0,	0,	0.8702,	0,	0,	0,	0.469,	0,	0,	0,	-0.04814,	0,
                  0,	0,	0,	1.268,	0,	0,	0,	0.6588,	0,	0,	0,	-0.03263,
                  -0.01132,	0,	0,	0,	0.02311,	0,	0,	0,	0.4235,	0,	0,	0,
                  0,	-0.1669,	0,	0,	0,	-0.08942,	0,	0,	0,	0.2934,	0,	0,
                  0,	0,	-0.1768,	0,	0,	0,	-0.04814,	0,	0,	0,	0.3285,	0,
                  0,	0,	0,	-0.179,	0,	0,	0,	-0.03263,	0,	0,	0,	0.3842),12,12)

(initial_covs <- list(G_a = Ga_ini,
                      Residual_a = Ra_ini))
is.positive.definite(initial_covs$G_a, tol=1e-8)
is.positive.definite(initial_covs$Residual_a, tol=1e-8)

ptm <- Sys.time()
# The ABLUP model using REML-EM algorithm
res.spa.A <- remlf90(fixed   = cbind(DBH30c_1,DBH30c_2,DBH30c_3,DBH30c_4,
                                     HT30c_1, HT30c_2, HT30c_3, HT30c_4,
                                     NSWGR3604c_1, NSWGR3604c_2, NSWGR3604c_3, NSWGR3604c_4) ~ -1 + proc,
                     generic = list(A = list(incidence = Z, precision = Ai,var.ini = initial_covs$G_a)),
                     data = data,
                     method = 'em',
                     var.ini = list(residual = initial_covs$Residual_a),
                     progsf90.options=c('use_yams','fact_once memory'))
summary(res.spa.A)

# Use cov2cor() to compute genetic and residual correlations
round(cov2cor(res.spa.A$var[[1]]),4)
round(cov2cor(res.spa.A$var[[2]]),4)

# The ABLUP model using REML-AI algorithm (just 1 iteration)
# Heretabilities for each trait-site combiantion (12 in total, 3 traits and 4 sites)
h2_1  <- 'G_2_2_1_1/(G_2_2_1_1+R_1_1)';       h2_2  <- 'G_2_2_2_2/(G_2_2_2_2+R_2_2)';       h2_3  <- 'G_2_2_3_3/(G_2_2_3_3+R_3_3)'
h2_4  <- 'G_2_2_4_4/(G_2_2_4_4+R_4_4)';       h2_5  <- 'G_2_2_5_5/(G_2_2_5_5+R_5_5)';       h2_6  <- 'G_2_2_6_6/(G_2_2_6_6+R_6_6)'
h2_7  <- 'G_2_2_7_7/(G_2_2_7_7+R_7_7)';       h2_8  <- 'G_2_2_8_8/(G_2_2_8_8+R_8_8)';       h2_9  <- 'G_2_2_9_9/(G_2_2_9_9+R_9_9)'
h2_10 <- 'G_2_2_10_10/(G_2_2_10_10+R_10_10)'; h2_11 <- 'G_2_2_11_11/(G_2_2_11_11+R_11_11)'; h2_12 <- 'G_2_2_12_12/(G_2_2_12_12+R_12_12)'

# Genetic correlations
r12<-'G_2_2_1_2/(G_2_2_1_1*G_2_2_2_2)**0.5'; r13<-'G_2_2_1_3/(G_2_2_1_1*G_2_2_3_3)**0.5'; r14<-'G_2_2_1_4/(G_2_2_1_1*G_2_2_4_4)**0.5'
r15<-'G_2_2_1_5/(G_2_2_1_1*G_2_2_5_5)**0.5'; r16<-'G_2_2_1_6/(G_2_2_1_1*G_2_2_6_6)**0.5'; r17<-'G_2_2_1_7/(G_2_2_1_1*G_2_2_7_7)**0.5'
r18<-'G_2_2_1_8/(G_2_2_1_1*G_2_2_8_8)**0.5'; r19<-'G_2_2_1_9/(G_2_2_1_1*G_2_2_9_9)**0.5'; r110<- 'G_2_2_1_10/(G_2_2_1_1*G_2_2_10_10)**0.5'
r111<- 'G_2_2_1_11/(G_2_2_1_1*G_2_2_11_11)**0.5'; r112<- 'G_2_2_1_12/(G_2_2_1_1*G_2_2_12_12)**0.5'

r23 <- 'G_2_2_2_3/(G_2_2_2_2*G_2_2_3_3)**0.5'; r24 <- 'G_2_2_2_4/(G_2_2_2_2*G_2_2_4_4)**0.5';    r25 <- 'G_2_2_2_5/(G_2_2_2_2*G_2_2_5_5)**0.5'
r26 <- 'G_2_2_2_6/(G_2_2_2_2*G_2_2_6_6)**0.5'; r27 <- 'G_2_2_2_7/(G_2_2_2_2*G_2_2_7_7)**0.5';    r28 <- 'G_2_2_2_8/(G_2_2_2_2*G_2_2_8_8)**0.5'
r29 <- 'G_2_2_2_9/(G_2_2_2_2*G_2_2_9_9)**0.5'; r210<- 'G_2_2_2_10/(G_2_2_2_2*G_2_2_10_10)**0.5'; r211<- 'G_2_2_2_11/(G_2_2_2_2*G_2_2_11_11)**0.5'
r212<- 'G_2_2_2_12/(G_2_2_2_2*G_2_2_12_12)**0.5'

r34 <- 'G_2_2_3_4/(G_2_2_3_3*G_2_2_4_4)**0.5';    r35 <- 'G_2_2_3_5/(G_2_2_3_3*G_2_2_5_5)**0.5';    r36 <- 'G_2_2_3_6/(G_2_2_3_3*G_2_2_6_6)**0.5'
r37 <- 'G_2_2_3_7/(G_2_2_3_3*G_2_2_7_7)**0.5';    r38 <- 'G_2_2_3_8/(G_2_2_3_3*G_2_2_8_8)**0.5';    r39 <- 'G_2_2_3_9/(G_2_2_3_3*G_2_2_9_9)**0.5'
r310<- 'G_2_2_3_10/(G_2_2_3_3*G_2_2_10_10)**0.5'; r311<- 'G_2_2_3_11/(G_2_2_3_3*G_2_2_11_11)**0.5'; r312<- 'G_2_2_3_12/(G_2_2_3_3*G_2_2_12_12)**0.5'

r45 <- 'G_2_2_4_5/(G_2_2_4_4*G_2_2_5_5)**0.5'; r46 <- 'G_2_2_4_6/(G_2_2_4_4*G_2_2_6_6)**0.5'; r47 <- 'G_2_2_4_7/(G_2_2_4_4*G_2_2_7_7)**0.5'
r48 <- 'G_2_2_4_8/(G_2_2_4_4*G_2_2_8_8)**0.5'; r49 <- 'G_2_2_4_9/(G_2_2_4_4*G_2_2_9_9)**0.5'; r410<- 'G_2_2_4_10/(G_2_2_4_4*G_2_2_10_10)**0.5'
r411<- 'G_2_2_4_11/(G_2_2_4_4*G_2_2_11_11)**0.5'; r412<- 'G_2_2_4_12/(G_2_2_4_4*G_2_2_12_12)**0.5'

r56 <- 'G_2_2_5_6/(G_2_2_5_5*G_2_2_6_6)**0.5'; r57 <- 'G_2_2_5_7/(G_2_2_5_5*G_2_2_7_7)**0.5'; r58 <- 'G_2_2_5_8/(G_2_2_5_5*G_2_2_8_8)**0.5'
r59 <- 'G_2_2_5_9/(G_2_2_5_5*G_2_2_9_9)**0.5'; r510<- 'G_2_2_5_10/(G_2_2_5_5*G_2_2_10_10)**0.5'; r511<- 'G_2_2_5_11/(G_2_2_5_5*G_2_2_11_11)**0.5'
r512<- 'G_2_2_5_12/(G_2_2_5_5*G_2_2_12_12)**0.5'

r67 <- 'G_2_2_6_7/(G_2_2_6_6*G_2_2_7_7)**0.5'; r68 <- 'G_2_2_6_8/(G_2_2_6_6*G_2_2_8_8)**0.5'; r69 <- 'G_2_2_6_9/(G_2_2_6_6*G_2_2_9_9)**0.5'
r610<- 'G_2_2_6_10/(G_2_2_6_6*G_2_2_10_10)**0.5'; r611<- 'G_2_2_6_11/(G_2_2_6_6*G_2_2_11_11)**0.5'; r612<- 'G_2_2_6_12/(G_2_2_6_6*G_2_2_12_12)**0.5'

r78 <- 'G_2_2_7_8/(G_2_2_7_7*G_2_2_8_8)**0.5'; r79 <- 'G_2_2_7_9/(G_2_2_7_7*G_2_2_9_9)**0.5'; r710<- 'G_2_2_7_10/(G_2_2_7_7*G_2_2_10_10)**0.5'
r711<- 'G_2_2_7_11/(G_2_2_7_7*G_2_2_11_11)**0.5'; r712<- 'G_2_2_7_12/(G_2_2_7_7*G_2_2_12_12)**0.5'

r89 <- 'G_2_2_8_9/(G_2_2_8_8*G_2_2_9_9)**0.5'; r810<- 'G_2_2_8_10/(G_2_2_8_8*G_2_2_10_10)**0.5'; r811<- 'G_2_2_8_11/(G_2_2_8_8*G_2_2_11_11)**0.5'
r812<- 'G_2_2_8_12/(G_2_2_8_8*G_2_2_12_12)**0.5'

r910<- 'G_2_2_9_10/(G_2_2_9_9*G_2_2_10_10)**0.5'; r911<- 'G_2_2_9_11/(G_2_2_9_9*G_2_2_11_11)**0.5'; r912<- 'G_2_2_9_12/(G_2_2_9_9*G_2_2_12_12)**0.5'

r1011<- 'G_2_2_10_11/(G_2_2_10_10*G_2_2_11_11)**0.5'; r1012<- 'G_2_2_10_12/(G_2_2_10_10*G_2_2_12_12)**0.5'

r1112<- 'G_2_2_11_12/(G_2_2_11_11*G_2_2_12_12)**0.5'

res.spa.A.ai <- remlf90(fixed   = cbind(DBH30c_1,DBH30c_2,DBH30c_3,DBH30c_4,
                                        HT30c_1, HT30c_2, HT30c_3, HT30c_4,
                                        NSWGR3604c_1, NSWGR3604c_2, NSWGR3604c_3, NSWGR3604c_4) ~ -1 + proc,
                        generic = list(A = list(incidence = Z, precision = Ai,var.ini = res.spa.A$var[[1]])),
                        data = data,
                        method = 'ai',
                        var.ini = list(residual = res.spa.A$var[[2]]),
                        progsf90.options=c('use_yams','fact_once memory','maxrounds 1',
                                           paste('se_covar_function h2_1', h2_1),  paste('se_covar_function h2_2', h2_2),
                                           paste('se_covar_function h2_3', h2_3),  paste('se_covar_function h2_4', h2_4),
                                           paste('se_covar_function h2_5', h2_5),  paste('se_covar_function h2_6', h2_6),
                                           paste('se_covar_function h2_7', h2_7),  paste('se_covar_function h2_8', h2_8),
                                           paste('se_covar_function h2_9', h2_9),  paste('se_covar_function h2_10', h2_10),                                           
                                           paste('se_covar_function h2_11', h2_11),paste('se_covar_function h2_12', h2_12),
                                           paste('se_covar_function r12  ', r12  ),paste('se_covar_function r13  ', r13  ),paste('se_covar_function r14  ', r14  ),
                                           paste('se_covar_function r15  ', r15  ),paste('se_covar_function r16  ', r16  ),paste('se_covar_function r17  ', r17  ),
                                           paste('se_covar_function r18  ', r18  ),paste('se_covar_function r19  ', r19  ),paste('se_covar_function r110 ', r110 ),
                                           paste('se_covar_function r111 ', r111 ),paste('se_covar_function r112 ', r112 ),paste('se_covar_function r23  ', r23  ),
                                           paste('se_covar_function r24  ', r24  ),paste('se_covar_function r25  ', r25  ),paste('se_covar_function r26  ', r26  ),
                                           paste('se_covar_function r27  ', r27  ),paste('se_covar_function r28  ', r28  ),paste('se_covar_function r29  ', r29  ),
                                           paste('se_covar_function r210 ', r210 ),paste('se_covar_function r211 ', r211 ),paste('se_covar_function r212 ', r212 ),
                                           paste('se_covar_function r34  ', r34  ),paste('se_covar_function r35  ', r35  ),paste('se_covar_function r36  ', r36  ),
                                           paste('se_covar_function r37  ', r37  ),paste('se_covar_function r38  ', r38  ),paste('se_covar_function r39  ', r39  ),
                                           paste('se_covar_function r310 ', r310 ),paste('se_covar_function r311 ', r311 ),paste('se_covar_function r312 ', r312 ),
                                           paste('se_covar_function r45  ', r45  ),paste('se_covar_function r46  ', r46  ),paste('se_covar_function r47  ', r47  ),
                                           paste('se_covar_function r48  ', r48  ),paste('se_covar_function r49  ', r49  ),paste('se_covar_function r410 ', r410 ),
                                           paste('se_covar_function r411 ', r411 ),paste('se_covar_function r412 ', r412 ),paste('se_covar_function r56  ', r56  ),
                                           paste('se_covar_function r57  ', r57  ),paste('se_covar_function r58  ', r58  ),paste('se_covar_function r59  ', r59  ),
                                           paste('se_covar_function r510 ', r510 ),paste('se_covar_function r511 ', r511 ),paste('se_covar_function r512 ', r512 ),
                                           paste('se_covar_function r67  ', r67  ),paste('se_covar_function r68  ', r68  ),paste('se_covar_function r69  ', r69  ),
                                           paste('se_covar_function r610 ', r610 ),paste('se_covar_function r611 ', r611 ),paste('se_covar_function r612 ', r612 ),
                                           paste('se_covar_function r78  ', r78  ),paste('se_covar_function r79  ', r79  ),paste('se_covar_function r710 ', r710 ),
                                           paste('se_covar_function r711 ', r711 ),paste('se_covar_function r712 ', r712 ),paste('se_covar_function r89  ', r89  ),
                                           paste('se_covar_function r810 ', r810 ),paste('se_covar_function r811 ', r811 ),paste('se_covar_function r812 ', r812 ),
                                           paste('se_covar_function r910 ', r910 ),paste('se_covar_function r911 ', r911 ),paste('se_covar_function r912 ', r912 ),
                                           paste('se_covar_function r1011', r1011),paste('se_covar_function r1012', r1012),paste('se_covar_function r1112', r1112)))

summary(res.spa.A.ai)

(time <- (Sys.time() - ptm))

##############################################################################################
# 5.1. Spearman correlations and plots of breeding values (BVs) for parents and offspring between regular and spatial analyses.
# 5.1.1. Number of parents and offspring
(nparents<-c(max(data$mum,data$dad)))
(nprogenies<-c(max(data$self)) - nparents)

# 5.1.2. Breeding values (BVs) for each trait-site combination
# DBH30S1
BVs_spa_DBH301_parents  <-res.spa.A$ranef$A$DBH30c_1$value[1:nparents]
BVs_spa_DBH301_progenies<-res.spa.A$ranef$A$DBH30c_1$value[nparents+1:nprogenies]
# HT30S1
BVs_spa_HT301_parents  <-res.spa.A$ranef$A$HT30c_1$value[1:nparents]
BVs_spa_HT301_progenies<-res.spa.A$ranef$A$HT30c_1$value[nparents+1:nprogenies]
# WGR36S1
BVs_spa_WGR361_parents  <-res.spa.A$ranef$A$NSWGR3604c_1$value[1:nparents]
BVs_spa_WGR361_progenies<-res.spa.A$ranef$A$NSWGR3604c_1$value[nparents+1:nprogenies]

# DBH30S2
BVs_spa_DBH302_parents  <-res.spa.A$ranef$A$DBH30c_2$value[1:nparents]
BVs_spa_DBH302_progenies<-res.spa.A$ranef$A$DBH30c_2$value[nparents+1:nprogenies]
# HT30S2
BVs_spa_HT302_parents  <-res.spa.A$ranef$A$HT30c_2$value[1:nparents]
BVs_spa_HT302_progenies<-res.spa.A$ranef$A$HT30c_2$value[nparents+1:nprogenies]
# WGR36S2
BVs_spa_WGR362_parents  <-res.spa.A$ranef$A$NSWGR3604c_2$value[1:nparents]
BVs_spa_WGR362_progenies<-res.spa.A$ranef$A$NSWGR3604c_2$value[nparents+1:nprogenies]

# DBH30S3
BVs_spa_DBH303_parents  <-res.spa.A$ranef$A$DBH30c_3$value[1:nparents]
BVs_spa_DBH303_progenies<-res.spa.A$ranef$A$DBH30c_3$value[nparents+1:nprogenies]
# HT30S3
BVs_spa_HT303_parents  <-res.spa.A$ranef$A$HT30c_3$value[1:nparents]
BVs_spa_HT303_progenies<-res.spa.A$ranef$A$HT30c_3$value[nparents+1:nprogenies]
# WGR36S3
BVs_spa_WGR363_parents  <-res.spa.A$ranef$A$NSWGR3604c_3$value[1:nparents]
BVs_spa_WGR363_progenies<-res.spa.A$ranef$A$NSWGR3604c_3$value[nparents+1:nprogenies]

# DBH30S4
BVs_spa_DBH304_parents  <-res.spa.A$ranef$A$DBH30c_4$value[1:nparents]
BVs_spa_DBH304_progenies<-res.spa.A$ranef$A$DBH30c_4$value[nparents+1:nprogenies]
# HT30S4
BVs_spa_HT304_parents  <-res.spa.A$ranef$A$HT30c_4$value[1:nparents]
BVs_spa_HT304_progenies<-res.spa.A$ranef$A$HT30c_4$value[nparents+1:nprogenies]
# WGR36S4
BVs_spa_WGR364_parents  <-res.spa.A$ranef$A$NSWGR3604c_4$value[1:nparents]
BVs_spa_WGR364_progenies<-res.spa.A$ranef$A$NSWGR3604c_4$value[nparents+1:nprogenies]

##############################################################################################
# 5.1.3. THEORETICAL ACCURACIES OF BVs
# 5.1.3.1. SEs
# DBH30S1
SEs_spa_DBH301_parents   <-res.spa.A$ranef$A$DBH30c_1$s.e.[1:nparents]
SEs_spa_DBH301_progenies <-res.spa.A$ranef$A$DBH30c_1$s.e.[nparents+1:nprogenies]
# HT30cS1
SEs_spa_HT301_parents   <-res.spa.A$ranef$A$HT30c_1$s.e.[1:nparents]
SEs_spa_HT301_progenies <-res.spa.A$ranef$A$HT30c_1$s.e.[nparents+1:nprogenies]
# WGR36S1
SEs_spa_WGR361_parents   <-res.spa.A$ranef$A$NSWGR3604c_1$s.e.[1:nparents]
SEs_spa_WGR361_progenies <-res.spa.A$ranef$A$NSWGR3604c_1$s.e.[nparents+1:nprogenies]

# DBH30cS2
SEs_spa_DBH302_parents   <-res.spa.A$ranef$A$DBH30c_2$s.e.[1:nparents]
SEs_spa_DBH302_progenies <-res.spa.A$ranef$A$DBH30c_2$s.e.[nparents+1:nprogenies]
# HT30cS2
SEs_spa_HT302_parents   <-res.spa.A$ranef$A$HT30c_2$s.e.[1:nparents]
SEs_spa_HT302_progenies <-res.spa.A$ranef$A$HT30c_2$s.e.[nparents+1:nprogenies]
# WGR36S2
SEs_spa_WGR362_parents   <-res.spa.A$ranef$A$NSWGR3604c_2$s.e.[1:nparents]
SEs_spa_WGR362_progenies <-res.spa.A$ranef$A$NSWGR3604c_2$s.e.[nparents+1:nprogenies]

# DBH30cS3
SEs_spa_DBH303_parents   <-res.spa.A$ranef$A$DBH30c_3$s.e.[1:nparents]
SEs_spa_DBH303_progenies <-res.spa.A$ranef$A$DBH30c_3$s.e.[nparents+1:nprogenies]
# HT30cS3
SEs_spa_HT303_parents   <-res.spa.A$ranef$A$HT30c_3$s.e.[1:nparents]
SEs_spa_HT303_progenies <-res.spa.A$ranef$A$HT30c_3$s.e.[nparents+1:nprogenies]
# WGR36S3
SEs_spa_WGR363_parents   <-res.spa.A$ranef$A$NSWGR3604c_3$s.e.[1:nparents]
SEs_spa_WGR363_progenies <-res.spa.A$ranef$A$NSWGR3604c_3$s.e.[nparents+1:nprogenies]

# DBH30cS4
SEs_spa_DBH304_parents   <-res.spa.A$ranef$A$DBH30c_4$s.e.[1:nparents]
SEs_spa_DBH304_progenies <-res.spa.A$ranef$A$DBH30c_4$s.e.[nparents+1:nprogenies]
# HT30cS4
SEs_spa_HT304_parents   <-res.spa.A$ranef$A$HT30c_4$s.e.[1:nparents]
SEs_spa_HT304_progenies <-res.spa.A$ranef$A$HT30c_4$s.e.[nparents+1:nprogenies]
# WGR36S4
SEs_spa_WGR364_parents   <-res.spa.A$ranef$A$NSWGR3604c_4$s.e.[1:nparents]
SEs_spa_WGR364_progenies <-res.spa.A$ranef$A$NSWGR3604c_4$s.e.[nparents+1:nprogenies]


# 5.1.3.2.  Prediction error variances (PEVs)
# DBH30S1
SE2s_spa_DBH301_parents   = SEs_spa_DBH301_parents  *SEs_spa_DBH301_parents  
SE2s_spa_DBH301_progenies = SEs_spa_DBH301_progenies*SEs_spa_DBH301_progenies
# HT30S1
SE2s_spa_HT301_parents    = SEs_spa_HT301_parents   *SEs_spa_HT301_parents   
SE2s_spa_HT301_progenies  = SEs_spa_HT301_progenies *SEs_spa_HT301_progenies
# WGR30S1
SE2s_spa_WGR361_parents    = SEs_spa_WGR361_parents   *SEs_spa_WGR361_parents   
SE2s_spa_WGR361_progenies  = SEs_spa_WGR361_progenies *SEs_spa_WGR361_progenies

# DBH30S2
SE2s_spa_DBH302_parents   = SEs_spa_DBH302_parents  *SEs_spa_DBH302_parents  
SE2s_spa_DBH302_progenies = SEs_spa_DBH302_progenies*SEs_spa_DBH302_progenies
# HT30S2
SE2s_spa_HT302_parents    = SEs_spa_HT302_parents   *SEs_spa_HT302_parents   
SE2s_spa_HT302_progenies  = SEs_spa_HT302_progenies *SEs_spa_HT302_progenies
# WGR30S2
SE2s_spa_WGR362_parents    = SEs_spa_WGR362_parents   *SEs_spa_WGR362_parents   
SE2s_spa_WGR362_progenies  = SEs_spa_WGR362_progenies *SEs_spa_WGR362_progenies

# DBH30S3
SE2s_spa_DBH303_parents   = SEs_spa_DBH303_parents  *SEs_spa_DBH303_parents  
SE2s_spa_DBH303_progenies = SEs_spa_DBH303_progenies*SEs_spa_DBH303_progenies
# HT30S3
SE2s_spa_HT303_parents    = SEs_spa_HT303_parents   *SEs_spa_HT303_parents   
SE2s_spa_HT303_progenies  = SEs_spa_HT303_progenies *SEs_spa_HT303_progenies
# WGR30S3
SE2s_spa_WGR363_parents    = SEs_spa_WGR363_parents   *SEs_spa_WGR363_parents   
SE2s_spa_WGR363_progenies  = SEs_spa_WGR363_progenies *SEs_spa_WGR363_progenies

# DBH30S4
SE2s_spa_DBH304_parents   = SEs_spa_DBH304_parents  *SEs_spa_DBH304_parents  
SE2s_spa_DBH304_progenies = SEs_spa_DBH304_progenies*SEs_spa_DBH304_progenies
# HT30S4
SE2s_spa_HT304_parents    = SEs_spa_HT304_parents   *SEs_spa_HT304_parents   
SE2s_spa_HT304_progenies  = SEs_spa_HT304_progenies *SEs_spa_HT304_progenies
# WGR30S4
SE2s_spa_WGR364_parents    = SEs_spa_WGR364_parents   *SEs_spa_WGR364_parents   
SE2s_spa_WGR364_progenies  = SEs_spa_WGR364_progenies *SEs_spa_WGR364_progenies

# 5.1.3.3. Fi
Fi   = diag(A) - 1
Fipa =Fi[1:nparents]
Fipr =Fi[nparents+1:nprogenies]
summary(Fipr)

# 5.1.3.4. Accuracies (ACCs)
# DBH30S1
ACC_spa_DBH301_parents   =sqrt(1-(SE2s_spa_DBH301_parents  /((1+Fipa)* res.spa.A$var$A[1,1])))
ACC_spa_DBH301_progenies =sqrt(1-(SE2s_spa_DBH301_progenies/((1+Fipr)* res.spa.A$var$A[1,1])))
# HT30S1
ACC_spa_HT301_parents    =sqrt(1-(SE2s_spa_HT301_parents  /((1+Fipa)* res.spa.A$var$A[5,5])))
ACC_spa_HT301_progenies  =sqrt(1-(SE2s_spa_HT301_progenies/((1+Fipr)* res.spa.A$var$A[5,5])))
# WGR30S1
ACC_spa_WGR361_parents    =sqrt(1-(SE2s_spa_WGR361_parents  /((1+Fipa)* res.spa.A$var$A[9,9])))
ACC_spa_WGR361_progenies  =sqrt(1-(SE2s_spa_WGR361_progenies/((1+Fipr)* res.spa.A$var$A[9,9])))

# DBH30S2
ACC_spa_DBH302_parents   =sqrt(1-(SE2s_spa_DBH302_parents  /((1+Fipa)* res.spa.A$var$A[2,2])))
ACC_spa_DBH302_progenies =sqrt(1-(SE2s_spa_DBH302_progenies/((1+Fipr)* res.spa.A$var$A[2,2])))
# HT30S2
ACC_spa_HT302_parents    =sqrt(1-(SE2s_spa_HT302_parents  /((1+Fipa)* res.spa.A$var$A[6,6])))
ACC_spa_HT302_progenies  =sqrt(1-(SE2s_spa_HT302_progenies/((1+Fipr)* res.spa.A$var$A[6,6])))
# WGR30S2
ACC_spa_WGR362_parents    =sqrt(1-(SE2s_spa_WGR362_parents  /((1+Fipa)* res.spa.A$var$A[10,10])))
ACC_spa_WGR362_progenies  =sqrt(1-(SE2s_spa_WGR362_progenies/((1+Fipr)* res.spa.A$var$A[10,10])))

# DBH30S3
ACC_spa_DBH303_parents   =sqrt(1-(SE2s_spa_DBH303_parents  /((1+Fipa)* res.spa.A$var$A[3,3])))
ACC_spa_DBH303_progenies =sqrt(1-(SE2s_spa_DBH303_progenies/((1+Fipr)* res.spa.A$var$A[3,3])))
# HT30S3
ACC_spa_HT303_parents    =sqrt(1-(SE2s_spa_HT303_parents  /((1+Fipa)* res.spa.A$var$A[7,7])))
ACC_spa_HT303_progenies  =sqrt(1-(SE2s_spa_HT303_progenies/((1+Fipr)* res.spa.A$var$A[7,7])))
# WGR30S3
ACC_spa_WGR363_parents    =sqrt(1-(SE2s_spa_WGR363_parents  /((1+Fipa)* res.spa.A$var$A[11,11])))
ACC_spa_WGR363_progenies  =sqrt(1-(SE2s_spa_WGR363_progenies/((1+Fipr)* res.spa.A$var$A[11,11])))

# DBH30S4
ACC_spa_DBH304_parents   =sqrt(1-(SE2s_spa_DBH304_parents  /((1+Fipa)* res.spa.A$var$A[4,4])))
ACC_spa_DBH304_progenies =sqrt(1-(SE2s_spa_DBH304_progenies/((1+Fipr)* res.spa.A$var$A[4,4])))
# HT30S4
ACC_spa_HT304_parents    =sqrt(1-(SE2s_spa_HT304_parents  /((1+Fipa)* res.spa.A$var$A[8,8])))
ACC_spa_HT304_progenies  =sqrt(1-(SE2s_spa_HT304_progenies/((1+Fipr)* res.spa.A$var$A[8,8])))
# WGR30S4
ACC_spa_WGR364_parents    =sqrt(1-(SE2s_spa_WGR364_parents  /((1+Fipa)* res.spa.A$var$A[12,12])))
ACC_spa_WGR364_progenies  =sqrt(1-(SE2s_spa_WGR364_progenies/((1+Fipr)* res.spa.A$var$A[12,12])))

# 5.1.3.5. OUTPUT BVs and ACCs
# Parents
BLUP_ACC_Pa_site1= cbind(BVs_spa_DBH301_parents,ACC_spa_DBH301_parents,BVs_spa_HT301_parents,ACC_spa_HT301_parents,
                         BVs_spa_WGR361_parents, ACC_spa_WGR361_parents)

BLUP_ACC_Pa_site2 = cbind(BVs_spa_DBH302_parents, ACC_spa_DBH302_parents, BVs_spa_HT302_parents,ACC_spa_HT302_parents,
                          BVs_spa_WGR362_parents, ACC_spa_WGR362_parents)

BLUP_ACC_Pa_site3 = cbind(BVs_spa_DBH303_parents, ACC_spa_DBH303_parents, BVs_spa_HT303_parents,ACC_spa_HT303_parents,
                          BVs_spa_WGR363_parents, ACC_spa_WGR363_parents)

BLUP_ACC_Pa_site4 = cbind(BVs_spa_DBH304_parents, ACC_spa_DBH304_parents, BVs_spa_HT304_parents,ACC_spa_HT304_parents,
                          BVs_spa_WGR364_parents, ACC_spa_WGR364_parents)

BLUP_Pa_Final1 = cbind(c(1:242), BLUP_ACC_Pa_site1,BLUP_ACC_Pa_site2,BLUP_ACC_Pa_site3,BLUP_ACC_Pa_site4)
colnames(BLUP_Pa_Final1)[1] = "self"; colnames(BLUP_Pa_Final1)[2] = "site"
head(BLUP_Pa_Final1)

BLUP_Pa_Final1 = as.data.frame(BLUP_Pa_Final1)

BLUP_Pa_FinalS1 = subset(BLUP_Pa_Final1)[,c(1,2:7)]
BLUP_Pa_FinalS2 = subset(BLUP_Pa_Final1)[,c(1,8:13)]
BLUP_Pa_FinalS3 = subset(BLUP_Pa_Final1)[,c(1,14:19)]
BLUP_Pa_FinalS4 = subset(BLUP_Pa_Final1)[,c(1,20:25)]
head(BLUP_Pa_FinalS1); head(BLUP_Pa_FinalS2); head(BLUP_Pa_FinalS3); head(BLUP_Pa_FinalS4)

BLUP_Pa_FinalS1 = cbind(BLUP_Pa_FinalS1[,1], rep(1,nparents),BLUP_Pa_FinalS1[,-1])
BLUP_Pa_FinalS2 = cbind(BLUP_Pa_FinalS2[,1], rep(2,nparents),BLUP_Pa_FinalS2[,-1])
BLUP_Pa_FinalS3 = cbind(BLUP_Pa_FinalS3[,1], rep(3,nparents),BLUP_Pa_FinalS3[,-1])
BLUP_Pa_FinalS4 = cbind(BLUP_Pa_FinalS4[,1], rep(4,nparents),BLUP_Pa_FinalS4[,-1])

names(BLUP_Pa_FinalS1) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Pa_FinalS2) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Pa_FinalS3) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Pa_FinalS4) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")

BLUP_Pa_Final = rbind(BLUP_Pa_FinalS1,BLUP_Pa_FinalS2,BLUP_Pa_FinalS3,BLUP_Pa_FinalS4)
write.table(BLUP_Pa_Final, file = "ABLUP_BVs_ACCs_Parent.txt", row.names = FALSE)

# Offspring
BLUP_ACC_Of_site1= cbind(BVs_spa_DBH301_progenies,ACC_spa_DBH301_progenies, BVs_spa_HT301_progenies, ACC_spa_HT301_progenies,
                         BVs_spa_WGR361_progenies,ACC_spa_WGR361_progenies)

BLUP_ACC_Of_site2 = cbind(BVs_spa_DBH302_progenies,ACC_spa_DBH302_progenies, BVs_spa_HT302_progenies, ACC_spa_HT302_progenies,
                          BVs_spa_WGR362_progenies,ACC_spa_WGR362_progenies)

BLUP_ACC_Of_site3 = cbind(BVs_spa_DBH303_progenies,ACC_spa_DBH303_progenies, BVs_spa_HT303_progenies, ACC_spa_HT303_progenies,
                          BVs_spa_WGR363_progenies,ACC_spa_WGR363_progenies)

BLUP_ACC_Of_site4 = cbind(BVs_spa_DBH304_progenies,ACC_spa_DBH304_progenies, BVs_spa_HT304_progenies, ACC_spa_HT304_progenies,
                          BVs_spa_WGR364_progenies,ACC_spa_WGR364_progenies)

BLUP_Of_Final1 = cbind(data$self, data$site, BLUP_ACC_Of_site1, BLUP_ACC_Of_site2, BLUP_ACC_Of_site3, BLUP_ACC_Of_site4)
colnames(BLUP_Of_Final1)[1] = "self"; colnames(BLUP_Of_Final1)[2] = "site"
head(BLUP_Of_Final1)

BLUP_Of_Final1 = as.data.frame(BLUP_Of_Final1)

BLUP_Of_FinalS1 = subset(BLUP_Of_Final1, site == 1)[,c(1:2,3:8)]
BLUP_Of_FinalS2 = subset(BLUP_Of_Final1, site == 2)[,c(1:2,9:14)]
BLUP_Of_FinalS3 = subset(BLUP_Of_Final1, site == 3)[,c(1:2,15:20)]
BLUP_Of_FinalS4 = subset(BLUP_Of_Final1, site == 4)[,c(1:2,21:26)]
head(BLUP_Of_FinalS1);head(BLUP_Of_FinalS2);head(BLUP_Of_FinalS3);head(BLUP_Of_FinalS4)

names(BLUP_Of_FinalS1) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Of_FinalS2) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Of_FinalS3) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Of_FinalS4) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")

BLUP_Of_Final = rbind(BLUP_Of_FinalS1,BLUP_Of_FinalS2,BLUP_Of_FinalS3,BLUP_Of_FinalS4)
write.table(BLUP_Of_Final, file = "ABLUP_BVs_ACCs_Offspring.txt", row.names = FALSE)

############################################################################################################
# 5.1.4. Estimated provenance effects (BLUEs)
(BLUEs_A_proc = res.spa.A$fixed$proc)

write.table(BLUEs_A_proc, file = "ABLUEs_proc.txt")


############################################################################################################
### 6. STANDARD MULTI-TRAIT MULTI-SITE INDIVIDUAL-TREE MIXED MODEL using combined pedigree-genomic relationships (HBLUP model)
# Strategy: using the Expectation-Maximization (EM) algorithm followed by one iteration with the Average Information (AI) algorithm 
#           to compute the approximated standard errors of the heritabilities and (co)variance components
# Initial (co)variances genetic (GH_ini) and residual (RH_ini) 
GH_ini = matrix(c(1.75,	1.224,	1.148,	1.166,	0.5447,	0.2376,	0.2063,	0.3452,	-0.2996,	-0.1836,	-0.3643,	-0.2727,
                 1.224,	1.861,	1.335,	1.266,	0.2814,	0.386,	0.342,	0.3978,	-0.2055,	-0.2568,	-0.2464,	-0.2296,
                 1.148,	1.335,	1.323,	1.173,	0.3182,	0.3292,	0.3662,	0.3668,	-0.2705,	-0.2934,	-0.2987,	-0.2678,
                 1.166,	1.266,	1.173,	1.354,	0.1917,	0.1436,	0.1641,	0.3335,	-0.3054,	-0.3108,	-0.3346,	-0.2652,
                 0.5447,	0.2814,	0.3182,	0.1917,	0.3302,	0.2377,	0.2013,	0.1994,	-0.1202,	-0.009187,	-0.1025,	-0.07811,
                 0.2376,	0.386,	0.3292,	0.1436,	0.2377,	0.3782,	0.2525,	0.2069,	-0.0523,	-0.003267,	-0.0125,	-0.01802,
                 0.2063,	0.342,	0.3662,	0.1641,	0.2013,	0.2525,	0.2923,	0.2222,	-0.1284,	-0.09402,	-0.1067,	-0.1147,
                 0.3452,	0.3978,	0.3668,	0.3335,	0.1994,	0.2069,	0.2222,	0.3029,	-0.1738,	-0.1074,	-0.1633,	-0.1538,
                 -0.2996,	-0.2055,	-0.2705,	-0.3054,	-0.1202,	-0.0523,	-0.1284,	-0.1738,	0.2996,	0.2691,	0.3165,	0.2839,
                 -0.1836,	-0.2568,	-0.2934,	-0.3108,	-0.009187,	-0.003267,	-0.09402,	-0.1074,	0.2691,	0.3237,	0.3271,	0.3073,
                 -0.3643,	-0.2464,	-0.2987,	-0.3346,	-0.1025,	-0.0125,	-0.1067,	-0.1633,	0.3165,	0.3271,	0.3915,	0.3494,
                 -0.2727,	-0.2296,	-0.2678,	-0.2652,	-0.07811,	-0.01802,	-0.1147,	-0.1538,	0.2839,	0.3073,	0.3494,	0.3423),12,12)

RH_ini = matrix(c(6.634,	0,	0,	0,	1.785,	0,	0,	0,	-0.01944,	0,	0,	0,
                0,	5.817,	0,	0,	0,	1.726,	0,	0,	0,	-0.1759,	0,	0,
                0,	0,	4.462,	0,	0,	0,	0.9083,	0,	0,	0,	-0.175,	0,
                0,	0,	0,	5.423,	0,	0,	0,	1.298,	0,	0,	0,	-0.1796,
                1.785,	0,	0,	0,	1.063,	0,	0,	0,	0.01897,	0,	0,	0,
                0,	1.726,	0,	0,	0,	1.159,	0,	0,	0,	-0.09382,	0,	0,
                0,	0,	0.9083,	0,	0,	0,	0.4852,	0,	0,	0,	-0.06207,	0,
                0,	0,	0,	1.298,	0,	0,	0,	0.6712,	0,	0,	0,	-0.03609,
                -0.01944,	0,	0,	0,	0.01897,	0,	0,	0,	0.4601,	0,	0,	0,
                0,	-0.1759,	0,	0,	0,	-0.09382,	0,	0,	0,	0.3546,	0,	0,
                0,	0,	-0.175,	0,	0,	0,	-0.06207,	0,	0,	0,	0.3502,	0,
                0,	0,	0,	-0.1796,	0,	0,	0,	-0.03609,	0,	0,	0,	0.3913),12,12)

(initial_covs2 <- list(G_h = GH_ini,
                      Residual_h = RH_ini))
is.positive.definite(initial_covs2$G_h, tol=1e-8)
is.positive.definite(initial_covs2$Residual_h, tol=1e-8)

start_time <- Sys.time()
# The HBLUP model using REML-EM algorithm
res.spa.H <- remlf90(fixed   = cbind(DBH30c_1,DBH30c_2,DBH30c_3,DBH30c_4,
                                     HT30c_1, HT30c_2, HT30c_3, HT30c_4,
                                     NSWGR3604c_1, NSWGR3604c_2, NSWGR3604c_3, NSWGR3604c_4) ~ -1 + proc,
                     generic = list(H = list(incidence = Z, precision = Hi,var.ini = initial_covs2$G_h)),
                     data = data,
                     method = 'em',
                     var.ini = list(residual = initial_covs2$Residual_h),
                     progsf90.options=c('use_yams','fact_once memory')) #test = initial_covs$test,

summary(res.spa.H)
# Use cov2cor() to compute correlations
round(cov2cor(res.spa.H$var$H),4)
round(cov2cor(res.spa.H$var$Residual),4)

# The HBLUP model using REML-AI algorithm (just 1 iteration)
# Heretabilities for each trait-site combiantion (12 in total, 3 traits and 4 sites)
h2_H_1  <- 'G_2_2_1_1/(G_2_2_1_1+R_1_1)';       h2_H_2  <- 'G_2_2_2_2/(G_2_2_2_2+R_2_2)';       h2_H_3  <- 'G_2_2_3_3/(G_2_2_3_3+R_3_3)'
h2_H_4  <- 'G_2_2_4_4/(G_2_2_4_4+R_4_4)';       h2_H_5  <- 'G_2_2_5_5/(G_2_2_5_5+R_5_5)';       h2_H_6  <- 'G_2_2_6_6/(G_2_2_6_6+R_6_6)'
h2_H_7  <- 'G_2_2_7_7/(G_2_2_7_7+R_7_7)';       h2_H_8  <- 'G_2_2_8_8/(G_2_2_8_8+R_8_8)';       h2_H_9  <- 'G_2_2_9_9/(G_2_2_9_9+R_9_9)'
h2_H_10 <- 'G_2_2_10_10/(G_2_2_10_10+R_10_10)'; h2_H_11 <- 'G_2_2_11_11/(G_2_2_11_11+R_11_11)'; h2_H_12 <- 'G_2_2_12_12/(G_2_2_12_12+R_12_12)'

# Genetic correlatios
r_H_12<-'G_2_2_1_2/(G_2_2_1_1*G_2_2_2_2)**0.5'; r_H_13<-'G_2_2_1_3/(G_2_2_1_1*G_2_2_3_3)**0.5'; r_H_14<-'G_2_2_1_4/(G_2_2_1_1*G_2_2_4_4)**0.5'
r_H_15<-'G_2_2_1_5/(G_2_2_1_1*G_2_2_5_5)**0.5'; r_H_16<-'G_2_2_1_6/(G_2_2_1_1*G_2_2_6_6)**0.5'; r_H_17<-'G_2_2_1_7/(G_2_2_1_1*G_2_2_7_7)**0.5'
r_H_18<-'G_2_2_1_8/(G_2_2_1_1*G_2_2_8_8)**0.5'; r_H_19<-'G_2_2_1_9/(G_2_2_1_1*G_2_2_9_9)**0.5'; r_H_110<- 'G_2_2_1_10/(G_2_2_1_1*G_2_2_10_10)**0.5'
r_H_111<- 'G_2_2_1_11/(G_2_2_1_1*G_2_2_11_11)**0.5'; r_H_112<- 'G_2_2_1_12/(G_2_2_1_1*G_2_2_12_12)**0.5'

r_H_23 <- 'G_2_2_2_3/(G_2_2_2_2*G_2_2_3_3)**0.5'; r_H_24 <- 'G_2_2_2_4/(G_2_2_2_2*G_2_2_4_4)**0.5'; r_H_25 <- 'G_2_2_2_5/(G_2_2_2_2*G_2_2_5_5)**0.5'
r_H_26 <- 'G_2_2_2_6/(G_2_2_2_2*G_2_2_6_6)**0.5'; r_H_27 <- 'G_2_2_2_7/(G_2_2_2_2*G_2_2_7_7)**0.5'; r_H_28 <- 'G_2_2_2_8/(G_2_2_2_2*G_2_2_8_8)**0.5'
r_H_29 <- 'G_2_2_2_9/(G_2_2_2_2*G_2_2_9_9)**0.5'; r_H_210<- 'G_2_2_2_10/(G_2_2_2_2*G_2_2_10_10)**0.5'; r_H_211<- 'G_2_2_2_11/(G_2_2_2_2*G_2_2_11_11)**0.5'
r_H_212<- 'G_2_2_2_12/(G_2_2_2_2*G_2_2_12_12)**0.5'

r_H_34 <- 'G_2_2_3_4/(G_2_2_3_3*G_2_2_4_4)**0.5'; r_H_35 <- 'G_2_2_3_5/(G_2_2_3_3*G_2_2_5_5)**0.5'; r_H_36 <- 'G_2_2_3_6/(G_2_2_3_3*G_2_2_6_6)**0.5'
r_H_37 <- 'G_2_2_3_7/(G_2_2_3_3*G_2_2_7_7)**0.5'; r_H_38 <- 'G_2_2_3_8/(G_2_2_3_3*G_2_2_8_8)**0.5'; r_H_39 <- 'G_2_2_3_9/(G_2_2_3_3*G_2_2_9_9)**0.5'
r_H_310<- 'G_2_2_3_10/(G_2_2_3_3*G_2_2_10_10)**0.5'; r_H_311<- 'G_2_2_3_11/(G_2_2_3_3*G_2_2_11_11)**0.5'; r_H_312<- 'G_2_2_3_12/(G_2_2_3_3*G_2_2_12_12)**0.5'

r_H_45 <- 'G_2_2_4_5/(G_2_2_4_4*G_2_2_5_5)**0.5'; r_H_46 <- 'G_2_2_4_6/(G_2_2_4_4*G_2_2_6_6)**0.5'; r_H_47 <- 'G_2_2_4_7/(G_2_2_4_4*G_2_2_7_7)**0.5'
r_H_48 <- 'G_2_2_4_8/(G_2_2_4_4*G_2_2_8_8)**0.5'; r_H_49 <- 'G_2_2_4_9/(G_2_2_4_4*G_2_2_9_9)**0.5'; r_H_410<- 'G_2_2_4_10/(G_2_2_4_4*G_2_2_10_10)**0.5'
r_H_411<- 'G_2_2_4_11/(G_2_2_4_4*G_2_2_11_11)**0.5'; r_H_412<- 'G_2_2_4_12/(G_2_2_4_4*G_2_2_12_12)**0.5'

r_H_56 <- 'G_2_2_5_6/(G_2_2_5_5*G_2_2_6_6)**0.5'; r_H_57 <- 'G_2_2_5_7/(G_2_2_5_5*G_2_2_7_7)**0.5'; r_H_58 <- 'G_2_2_5_8/(G_2_2_5_5*G_2_2_8_8)**0.5'
r_H_59 <- 'G_2_2_5_9/(G_2_2_5_5*G_2_2_9_9)**0.5'; r_H_510<- 'G_2_2_5_10/(G_2_2_5_5*G_2_2_10_10)**0.5'; r_H_511<- 'G_2_2_5_11/(G_2_2_5_5*G_2_2_11_11)**0.5'
r_H_512<- 'G_2_2_5_12/(G_2_2_5_5*G_2_2_12_12)**0.5'

r_H_67 <- 'G_2_2_6_7/(G_2_2_6_6*G_2_2_7_7)**0.5'; r_H_68 <- 'G_2_2_6_8/(G_2_2_6_6*G_2_2_8_8)**0.5'; r_H_69 <- 'G_2_2_6_9/(G_2_2_6_6*G_2_2_9_9)**0.5'
r_H_610<- 'G_2_2_6_10/(G_2_2_6_6*G_2_2_10_10)**0.5'; r_H_611<- 'G_2_2_6_11/(G_2_2_6_6*G_2_2_11_11)**0.5'; r_H_612<- 'G_2_2_6_12/(G_2_2_6_6*G_2_2_12_12)**0.5'

r_H_78 <- 'G_2_2_7_8/(G_2_2_7_7*G_2_2_8_8)**0.5'; r_H_79 <- 'G_2_2_7_9/(G_2_2_7_7*G_2_2_9_9)**0.5'; r_H_710<- 'G_2_2_7_10/(G_2_2_7_7*G_2_2_10_10)**0.5'
r_H_711<- 'G_2_2_7_11/(G_2_2_7_7*G_2_2_11_11)**0.5'; r_H_712<- 'G_2_2_7_12/(G_2_2_7_7*G_2_2_12_12)**0.5'

r_H_89 <- 'G_2_2_8_9/(G_2_2_8_8*G_2_2_9_9)**0.5'; r_H_810<- 'G_2_2_8_10/(G_2_2_8_8*G_2_2_10_10)**0.5'; r_H_811<- 'G_2_2_8_11/(G_2_2_8_8*G_2_2_11_11)**0.5'
r_H_812<- 'G_2_2_8_12/(G_2_2_8_8*G_2_2_12_12)**0.5'

r_H_910<- 'G_2_2_9_10/(G_2_2_9_9*G_2_2_10_10)**0.5'; r_H_911<- 'G_2_2_9_11/(G_2_2_9_9*G_2_2_11_11)**0.5'; r_H_912<- 'G_2_2_9_12/(G_2_2_9_9*G_2_2_12_12)**0.5'

r_H_1011<- 'G_2_2_10_11/(G_2_2_10_10*G_2_2_11_11)**0.5'; r_H_1012<- 'G_2_2_10_12/(G_2_2_10_10*G_2_2_12_12)**0.5'

r_H_1112<- 'G_2_2_11_12/(G_2_2_11_11*G_2_2_12_12)**0.5'

res.spa.H.ai <- remlf90(fixed   = cbind(DBH30c_1,DBH30c_2,DBH30c_3,DBH30c_4,
                                        HT30c_1, HT30c_2, HT30c_3, HT30c_4,
                                        NSWGR3604c_1, NSWGR3604c_2, NSWGR3604c_3, NSWGR3604c_4) ~ -1 + proc,
                        generic = list(H = list(incidence = Z, precision = Hi,var.ini = res.spa.H$var$H)),
                        data = data,
                        method = 'ai',
                        var.ini = list(residual = res.spa.H$var$Residual),
                        progsf90.options=c('use_yams','fact_once memory','maxrounds 1',
                                           paste('se_covar_H__function h2_H_1', h2_H_1), paste('se_covar_H__function h2_H_2', h2_H_2),
                                           paste('se_covar_H__function h2_H_3', h2_H_3), paste('se_covar_H__function h2_H_4', h2_H_4),
                                           paste('se_covar_H__function h2_H_5', h2_H_5), paste('se_covar_H__function h2_H_6', h2_H_6),
                                           paste('se_covar_H__function h2_H_7', h2_H_7), paste('se_covar_H__function h2_H_8', h2_H_8),
                                           paste('se_covar_H__function h2_H_9', h2_H_9), paste('se_covar_H__function h2_H_10', h2_H_10),                                           
                                           paste('se_covar_H__function h2_H_11', h2_H_11), paste('se_covar_H__function h2_H_12', h2_H_12),
                                           paste('se_covar_H__function r_H_12  ', r_H_12  ),paste('se_covar_H__function r_H_13  ', r_H_13  ),paste('se_covar_H__function r_H_14  ', r_H_14  ),
                                           paste('se_covar_H__function r_H_15  ', r_H_15  ),paste('se_covar_H__function r_H_16  ', r_H_16  ),paste('se_covar_H__function r_H_17  ', r_H_17  ),
                                           paste('se_covar_H__function r_H_18  ', r_H_18  ),paste('se_covar_H__function r_H_19  ', r_H_19  ),paste('se_covar_H__function r_H_110 ', r_H_110 ),
                                           paste('se_covar_H__function r_H_111 ', r_H_111 ),paste('se_covar_H__function r_H_112 ', r_H_112 ),paste('se_covar_H__function r_H_23  ', r_H_23  ),
                                           paste('se_covar_H__function r_H_24  ', r_H_24  ),paste('se_covar_H__function r_H_25  ', r_H_25  ),paste('se_covar_H__function r_H_26  ', r_H_26  ),
                                           paste('se_covar_H__function r_H_27  ', r_H_27  ),paste('se_covar_H__function r_H_28  ', r_H_28  ),paste('se_covar_H__function r_H_29  ', r_H_29  ),
                                           paste('se_covar_H__function r_H_210 ', r_H_210 ),paste('se_covar_H__function r_H_211 ', r_H_211 ),paste('se_covar_H__function r_H_212 ', r_H_212 ),
                                           paste('se_covar_H__function r_H_34  ', r_H_34  ),paste('se_covar_H__function r_H_35  ', r_H_35  ),paste('se_covar_H__function r_H_36  ', r_H_36  ),
                                           paste('se_covar_H__function r_H_37  ', r_H_37  ),paste('se_covar_H__function r_H_38  ', r_H_38  ),paste('se_covar_H__function r_H_39  ', r_H_39  ),
                                           paste('se_covar_H__function r_H_310 ', r_H_310 ),paste('se_covar_H__function r_H_311 ', r_H_311 ),paste('se_covar_H__function r_H_312 ', r_H_312 ),
                                           paste('se_covar_H__function r_H_45  ', r_H_45  ),paste('se_covar_H__function r_H_46  ', r_H_46  ),paste('se_covar_H__function r_H_47  ', r_H_47  ),
                                           paste('se_covar_H__function r_H_48  ', r_H_48  ),paste('se_covar_H__function r_H_49  ', r_H_49  ),paste('se_covar_H__function r_H_410 ', r_H_410 ),
                                           paste('se_covar_H__function r_H_411 ', r_H_411 ),paste('se_covar_H__function r_H_412 ', r_H_412 ),paste('se_covar_H__function r_H_56  ', r_H_56  ),
                                           paste('se_covar_H__function r_H_57  ', r_H_57  ),paste('se_covar_H__function r_H_58  ', r_H_58  ),paste('se_covar_H__function r_H_59  ', r_H_59  ),
                                           paste('se_covar_H__function r_H_510 ', r_H_510 ),paste('se_covar_H__function r_H_511 ', r_H_511 ),paste('se_covar_H__function r_H_512 ', r_H_512 ),
                                           paste('se_covar_H__function r_H_67  ', r_H_67  ),paste('se_covar_H__function r_H_68  ', r_H_68  ),paste('se_covar_H__function r_H_69  ', r_H_69  ),
                                           paste('se_covar_H__function r_H_610 ', r_H_610 ),paste('se_covar_H__function r_H_611 ', r_H_611 ),paste('se_covar_H__function r_H_612 ', r_H_612 ),
                                           paste('se_covar_H__function r_H_78  ', r_H_78  ),paste('se_covar_H__function r_H_79  ', r_H_79  ),paste('se_covar_H__function r_H_710 ', r_H_710 ),
                                           paste('se_covar_H__function r_H_711 ', r_H_711 ),paste('se_covar_H__function r_H_712 ', r_H_712 ),paste('se_covar_H__function r_H_89  ', r_H_89  ),
                                           paste('se_covar_H__function r_H_810 ', r_H_810 ),paste('se_covar_H__function r_H_811 ', r_H_811 ),paste('se_covar_H__function r_H_812 ', r_H_812 ),
                                           paste('se_covar_H__function r_H_910 ', r_H_910 ),paste('se_covar_H__function r_H_911 ', r_H_911 ),paste('se_covar_H__function r_H_912 ', r_H_912 ),
                                           paste('se_covar_H__function r_H_1011', r_H_1011),paste('se_covar_H__function r_H_1012', r_H_1012),paste('se_covar_H__function r_H_1112', r_H_1112)))

summary(res.spa.H.ai)

end_time <- Sys.time()
print(end_time - start_time)

##############################################################################################
# 6.1. Spearman correlations and plots of breeding values (BVs) for parents and offspring between regular and spatial analyses.
# 6.1.1. Number of parents and offspring
#parents<-as.matrix(summary(as.data.frame(get_pedigree(res.blk))$dam))
(nparents<-c(max(data$mum,data$dad)))
(nprogenies<-c(max(data$self)) - nparents)

# 6.1.2. Breeding values (BVs) for each trait-site combination
# DBH30S1
BVs_spa_DBH301_parents  <-res.spa.H$ranef$H$DBH30c_1$value[1:nparents]
BVs_spa_DBH301_progenies<-res.spa.H$ranef$H$DBH30c_1$value[nparents+1:nprogenies]
# HT30S1
BVs_spa_HT301_parents  <-res.spa.H$ranef$H$HT30c_1$value[1:nparents]
BVs_spa_HT301_progenies<-res.spa.H$ranef$H$HT30c_1$value[nparents+1:nprogenies]
# WGR36S1
BVs_spa_WGR361_parents  <-res.spa.H$ranef$H$NSWGR3604c_1$value[1:nparents]
BVs_spa_WGR361_progenies<-res.spa.H$ranef$H$NSWGR3604c_1$value[nparents+1:nprogenies]

# DBH30S2
BVs_spa_DBH302_parents  <-res.spa.H$ranef$H$DBH30c_2$value[1:nparents]
BVs_spa_DBH302_progenies<-res.spa.H$ranef$H$DBH30c_2$value[nparents+1:nprogenies]
# HT30S2
BVs_spa_HT302_parents  <-res.spa.H$ranef$H$HT30c_2$value[1:nparents]
BVs_spa_HT302_progenies<-res.spa.H$ranef$H$HT30c_2$value[nparents+1:nprogenies]
# WGR36S2
BVs_spa_WGR362_parents  <-res.spa.H$ranef$H$NSWGR3604c_2$value[1:nparents]
BVs_spa_WGR362_progenies<-res.spa.H$ranef$H$NSWGR3604c_2$value[nparents+1:nprogenies]

# DBH30S3
BVs_spa_DBH303_parents  <-res.spa.H$ranef$H$DBH30c_3$value[1:nparents]
BVs_spa_DBH303_progenies<-res.spa.H$ranef$H$DBH30c_3$value[nparents+1:nprogenies]
# HT30S3
BVs_spa_HT303_parents  <-res.spa.H$ranef$H$HT30c_3$value[1:nparents]
BVs_spa_HT303_progenies<-res.spa.H$ranef$H$HT30c_3$value[nparents+1:nprogenies]
# WGR36S3
BVs_spa_WGR363_parents  <-res.spa.H$ranef$H$NSWGR3604c_3$value[1:nparents]
BVs_spa_WGR363_progenies<-res.spa.H$ranef$H$NSWGR3604c_3$value[nparents+1:nprogenies]

# DBH30S4
BVs_spa_DBH304_parents  <-res.spa.H$ranef$H$DBH30c_4$value[1:nparents]
BVs_spa_DBH304_progenies<-res.spa.H$ranef$H$DBH30c_4$value[nparents+1:nprogenies]
# HT30S4
BVs_spa_HT304_parents  <-res.spa.H$ranef$H$HT30c_4$value[1:nparents]
BVs_spa_HT304_progenies<-res.spa.H$ranef$H$HT30c_4$value[nparents+1:nprogenies]
# WGR36S4
BVs_spa_WGR364_parents  <-res.spa.H$ranef$H$NSWGR3604c_4$value[1:nparents]
BVs_spa_WGR364_progenies<-res.spa.H$ranef$H$NSWGR3604c_4$value[nparents+1:nprogenies]

##############################################################################################
# 6.1.3. THEORETICAL ACCURACIES OF BVs
# 6.1.3.1. SEs
# DBH30S1
SEs_spa_DBH301_parents   <-res.spa.H$ranef$H$DBH30c_1$s.e.[1:nparents]
SEs_spa_DBH301_progenies <-res.spa.H$ranef$H$DBH30c_1$s.e.[nparents+1:nprogenies]
# HT30cS1
SEs_spa_HT301_parents   <-res.spa.H$ranef$H$HT30c_1$s.e.[1:nparents]
SEs_spa_HT301_progenies <-res.spa.H$ranef$H$HT30c_1$s.e.[nparents+1:nprogenies]
# WGR36S1
SEs_spa_WGR361_parents   <-res.spa.H$ranef$H$NSWGR3604c_1$s.e.[1:nparents]
SEs_spa_WGR361_progenies <-res.spa.H$ranef$H$NSWGR3604c_1$s.e.[nparents+1:nprogenies]

# DBH30cS2
SEs_spa_DBH302_parents   <-res.spa.H$ranef$H$DBH30c_2$s.e.[1:nparents]
SEs_spa_DBH302_progenies <-res.spa.H$ranef$H$DBH30c_2$s.e.[nparents+1:nprogenies]
# HT30cS2
SEs_spa_HT302_parents   <-res.spa.H$ranef$H$HT30c_2$s.e.[1:nparents]
SEs_spa_HT302_progenies <-res.spa.H$ranef$H$HT30c_2$s.e.[nparents+1:nprogenies]
# WGR36S2
SEs_spa_WGR362_parents   <-res.spa.H$ranef$H$NSWGR3604c_2$s.e.[1:nparents]
SEs_spa_WGR362_progenies <-res.spa.H$ranef$H$NSWGR3604c_2$s.e.[nparents+1:nprogenies]

# DBH30cS3
SEs_spa_DBH303_parents   <-res.spa.H$ranef$H$DBH30c_3$s.e.[1:nparents]
SEs_spa_DBH303_progenies <-res.spa.H$ranef$H$DBH30c_3$s.e.[nparents+1:nprogenies]
# HT30cS3
SEs_spa_HT303_parents   <-res.spa.H$ranef$H$HT30c_3$s.e.[1:nparents]
SEs_spa_HT303_progenies <-res.spa.H$ranef$H$HT30c_3$s.e.[nparents+1:nprogenies]
# WGR36S3
SEs_spa_WGR363_parents   <-res.spa.H$ranef$H$NSWGR3604c_3$s.e.[1:nparents]
SEs_spa_WGR363_progenies <-res.spa.H$ranef$H$NSWGR3604c_3$s.e.[nparents+1:nprogenies]

# DBH30cS4
SEs_spa_DBH304_parents   <-res.spa.H$ranef$H$DBH30c_4$s.e.[1:nparents]
SEs_spa_DBH304_progenies <-res.spa.H$ranef$H$DBH30c_4$s.e.[nparents+1:nprogenies]
# HT30cS4
SEs_spa_HT304_parents   <-res.spa.H$ranef$H$HT30c_4$s.e.[1:nparents]
SEs_spa_HT304_progenies <-res.spa.H$ranef$H$HT30c_4$s.e.[nparents+1:nprogenies]
# WGR36S4
SEs_spa_WGR364_parents   <-res.spa.H$ranef$H$NSWGR3604c_4$s.e.[1:nparents]
SEs_spa_WGR364_progenies <-res.spa.H$ranef$H$NSWGR3604c_4$s.e.[nparents+1:nprogenies]


# 6.1.3.2.  Prediction error variances (PEVs)
# DBH30S1
SE2s_spa_DBH301_parents   = SEs_spa_DBH301_parents  *SEs_spa_DBH301_parents  
SE2s_spa_DBH301_progenies = SEs_spa_DBH301_progenies*SEs_spa_DBH301_progenies
# HT30S1
SE2s_spa_HT301_parents    = SEs_spa_HT301_parents   *SEs_spa_HT301_parents   
SE2s_spa_HT301_progenies  = SEs_spa_HT301_progenies *SEs_spa_HT301_progenies
# WGR30S1
SE2s_spa_WGR361_parents    = SEs_spa_WGR361_parents   *SEs_spa_WGR361_parents   
SE2s_spa_WGR361_progenies  = SEs_spa_WGR361_progenies *SEs_spa_WGR361_progenies

# DBH30S2
SE2s_spa_DBH302_parents   = SEs_spa_DBH302_parents  *SEs_spa_DBH302_parents  
SE2s_spa_DBH302_progenies = SEs_spa_DBH302_progenies*SEs_spa_DBH302_progenies
# HT30S2
SE2s_spa_HT302_parents    = SEs_spa_HT302_parents   *SEs_spa_HT302_parents   
SE2s_spa_HT302_progenies  = SEs_spa_HT302_progenies *SEs_spa_HT302_progenies
# WGR30S2
SE2s_spa_WGR362_parents    = SEs_spa_WGR362_parents   *SEs_spa_WGR362_parents   
SE2s_spa_WGR362_progenies  = SEs_spa_WGR362_progenies *SEs_spa_WGR362_progenies

# DBH30S3
SE2s_spa_DBH303_parents   = SEs_spa_DBH303_parents  *SEs_spa_DBH303_parents  
SE2s_spa_DBH303_progenies = SEs_spa_DBH303_progenies*SEs_spa_DBH303_progenies
# HT30S3
SE2s_spa_HT303_parents    = SEs_spa_HT303_parents   *SEs_spa_HT303_parents   
SE2s_spa_HT303_progenies  = SEs_spa_HT303_progenies *SEs_spa_HT303_progenies
# WGR30S3
SE2s_spa_WGR363_parents    = SEs_spa_WGR363_parents   *SEs_spa_WGR363_parents   
SE2s_spa_WGR363_progenies  = SEs_spa_WGR363_progenies *SEs_spa_WGR363_progenies

# DBH30S4
SE2s_spa_DBH304_parents   = SEs_spa_DBH304_parents  *SEs_spa_DBH304_parents  
SE2s_spa_DBH304_progenies = SEs_spa_DBH304_progenies*SEs_spa_DBH304_progenies
# HT30S4
SE2s_spa_HT304_parents    = SEs_spa_HT304_parents   *SEs_spa_HT304_parents   
SE2s_spa_HT304_progenies  = SEs_spa_HT304_progenies *SEs_spa_HT304_progenies
# WGR30S4
SE2s_spa_WGR364_parents    = SEs_spa_WGR364_parents   *SEs_spa_WGR364_parents   
SE2s_spa_WGR364_progenies  = SEs_spa_WGR364_progenies *SEs_spa_WGR364_progenies

# 6.1.3.3. Fi
Fi_H   = diag(solve(Hi)) - 1
Fipa =Fi_H[1:nparents]
Fipr =Fi_H[nparents+1:nprogenies]
summary(Fipr)

# 6.1.3.4. Accuracies (ACCs)
# DBH30S1
ACC_spa_DBH301_parents   =sqrt(1-(SE2s_spa_DBH301_parents  /((1+Fipa)* res.spa.H$var$H[1,1])))
ACC_spa_DBH301_progenies =sqrt(1-(SE2s_spa_DBH301_progenies/((1+Fipr)* res.spa.H$var$H[1,1])))
# HT30S1
ACC_spa_HT301_parents    =sqrt(1-(SE2s_spa_HT301_parents  /((1+Fipa)* res.spa.H$var$H[5,5])))
ACC_spa_HT301_progenies  =sqrt(1-(SE2s_spa_HT301_progenies/((1+Fipr)* res.spa.H$var$H[5,5])))
# WGR30S1
ACC_spa_WGR361_parents    =sqrt(1-(SE2s_spa_WGR361_parents  /((1+Fipa)* res.spa.H$var$H[9,9])))
ACC_spa_WGR361_progenies  =sqrt(1-(SE2s_spa_WGR361_progenies/((1+Fipr)* res.spa.H$var$H[9,9])))

# DBH30S2
ACC_spa_DBH302_parents   =sqrt(1-(SE2s_spa_DBH302_parents  /((1+Fipa)* res.spa.H$var$H[2,2])))
ACC_spa_DBH302_progenies =sqrt(1-(SE2s_spa_DBH302_progenies/((1+Fipr)* res.spa.H$var$H[2,2])))
# HT30S2
ACC_spa_HT302_parents    =sqrt(1-(SE2s_spa_HT302_parents  /((1+Fipa)* res.spa.H$var$H[6,6])))
ACC_spa_HT302_progenies  =sqrt(1-(SE2s_spa_HT302_progenies/((1+Fipr)* res.spa.H$var$H[6,6])))
# WGR30S2
ACC_spa_WGR362_parents    =sqrt(1-(SE2s_spa_WGR362_parents  /((1+Fipa)* res.spa.H$var$H[10,10])))
ACC_spa_WGR362_progenies  =sqrt(1-(SE2s_spa_WGR362_progenies/((1+Fipr)* res.spa.H$var$H[10,10])))

# DBH30S3
ACC_spa_DBH303_parents   =sqrt(1-(SE2s_spa_DBH303_parents  /((1+Fipa)* res.spa.H$var$H[3,3])))
ACC_spa_DBH303_progenies =sqrt(1-(SE2s_spa_DBH303_progenies/((1+Fipr)* res.spa.H$var$H[3,3])))
# HT30S3
ACC_spa_HT303_parents    =sqrt(1-(SE2s_spa_HT303_parents  /((1+Fipa)* res.spa.H$var$H[7,7])))
ACC_spa_HT303_progenies  =sqrt(1-(SE2s_spa_HT303_progenies/((1+Fipr)* res.spa.H$var$H[7,7])))
# WGR30S3
ACC_spa_WGR363_parents    =sqrt(1-(SE2s_spa_WGR363_parents  /((1+Fipa)* res.spa.H$var$H[11,11])))
ACC_spa_WGR363_progenies  =sqrt(1-(SE2s_spa_WGR363_progenies/((1+Fipr)* res.spa.H$var$H[11,11])))

# DBH30S4
ACC_spa_DBH304_parents   =sqrt(1-(SE2s_spa_DBH304_parents  /((1+Fipa)* res.spa.H$var$H[4,4])))
ACC_spa_DBH304_progenies =sqrt(1-(SE2s_spa_DBH304_progenies/((1+Fipr)* res.spa.H$var$H[4,4])))
# HT30S4
ACC_spa_HT304_parents    =sqrt(1-(SE2s_spa_HT304_parents  /((1+Fipa)* res.spa.H$var$H[8,8])))
ACC_spa_HT304_progenies  =sqrt(1-(SE2s_spa_HT304_progenies/((1+Fipr)* res.spa.H$var$H[8,8])))
# WGR30S4
ACC_spa_WGR364_parents    =sqrt(1-(SE2s_spa_WGR364_parents  /((1+Fipa)* res.spa.H$var$H[12,12])))
ACC_spa_WGR364_progenies  =sqrt(1-(SE2s_spa_WGR364_progenies/((1+Fipr)* res.spa.H$var$H[12,12])))

# 6.1.3.5. OUTPUT BVs and ACCs
# Parents
BLUP_ACC_Pa_site1= cbind(BVs_spa_DBH301_parents, ACC_spa_DBH301_parents, BVs_spa_HT301_parents, ACC_spa_HT301_parents,
                         BVs_spa_WGR361_parents, ACC_spa_WGR361_parents)

BLUP_ACC_Pa_site2 = cbind(BVs_spa_DBH302_parents, ACC_spa_DBH302_parents, BVs_spa_HT302_parents, ACC_spa_HT302_parents,
                          BVs_spa_WGR362_parents, ACC_spa_WGR362_parents)

BLUP_ACC_Pa_site3 = cbind(BVs_spa_DBH303_parents, ACC_spa_DBH303_parents, BVs_spa_HT303_parents, ACC_spa_HT303_parents,
                          BVs_spa_WGR363_parents, ACC_spa_WGR363_parents)

BLUP_ACC_Pa_site4 = cbind(BVs_spa_DBH304_parents, ACC_spa_DBH304_parents, BVs_spa_HT304_parents, ACC_spa_HT304_parents,
                          BVs_spa_WGR364_parents, ACC_spa_WGR364_parents)

BLUP_Pa_Final1 = cbind(c(1:242), BLUP_ACC_Pa_site1,BLUP_ACC_Pa_site2,BLUP_ACC_Pa_site3,BLUP_ACC_Pa_site4)
colnames(BLUP_Pa_Final1)[1] = "self"; colnames(BLUP_Pa_Final1)[2] = "site"
head(BLUP_Pa_Final1);dim(BLUP_Pa_Final1)

BLUP_Pa_Final1 = as.data.frame(BLUP_Pa_Final1)

BLUP_Pa_FinalS1 = subset(BLUP_Pa_Final1)[,c(1,2:7)]
BLUP_Pa_FinalS2 = subset(BLUP_Pa_Final1)[,c(1,8:13)]
BLUP_Pa_FinalS3 = subset(BLUP_Pa_Final1)[,c(1,14:19)]
BLUP_Pa_FinalS4 = subset(BLUP_Pa_Final1)[,c(1,20:25)]
head(BLUP_Pa_FinalS1); head(BLUP_Pa_FinalS2); head(BLUP_Pa_FinalS3); head(BLUP_Pa_FinalS4)

BLUP_Pa_FinalS1 = cbind(BLUP_Pa_FinalS1[,1], rep(1,242),BLUP_Pa_FinalS1[,-1])
BLUP_Pa_FinalS2 = cbind(BLUP_Pa_FinalS2[,1], rep(2,242),BLUP_Pa_FinalS2[,-1])
BLUP_Pa_FinalS3 = cbind(BLUP_Pa_FinalS3[,1], rep(3,242),BLUP_Pa_FinalS3[,-1])
BLUP_Pa_FinalS4 = cbind(BLUP_Pa_FinalS4[,1], rep(4,242),BLUP_Pa_FinalS4[,-1])

names(BLUP_Pa_FinalS1) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Pa_FinalS2) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Pa_FinalS3) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Pa_FinalS4) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")

BLUP_Pa_Final = rbind(BLUP_Pa_FinalS1,BLUP_Pa_FinalS2,BLUP_Pa_FinalS3,BLUP_Pa_FinalS4)
write.table(BLUP_Pa_Final, file = "HBLUP_BVs_ACCs_Parents.txt", row.names = FALSE)

# Offspring
BLUP_ACC_Of_site1= cbind(BVs_spa_DBH301_progenies, ACC_spa_DBH301_progenies, BVs_spa_HT301_progenies,ACC_spa_HT301_progenies,
                         BVs_spa_WGR361_progenies, ACC_spa_WGR361_progenies)

BLUP_ACC_Of_site2 = cbind(BVs_spa_DBH302_progenies, ACC_spa_DBH302_progenies, BVs_spa_HT302_progenies,ACC_spa_HT302_progenies,
                          BVs_spa_WGR362_progenies, ACC_spa_WGR362_progenies)

BLUP_ACC_Of_site3 = cbind(BVs_spa_DBH303_progenies, ACC_spa_DBH303_progenies, BVs_spa_HT303_progenies,ACC_spa_HT303_progenies,
                          BVs_spa_WGR363_progenies, ACC_spa_WGR363_progenies)

BLUP_ACC_Of_site4 = cbind(BVs_spa_DBH304_progenies, ACC_spa_DBH304_progenies, BVs_spa_HT304_progenies,ACC_spa_HT304_progenies,
                          BVs_spa_WGR364_progenies, ACC_spa_WGR364_progenies)

BLUP_Of_Final1 = cbind(data$self,data$site, BLUP_ACC_Of_site1,BLUP_ACC_Of_site2,BLUP_ACC_Of_site3,BLUP_ACC_Of_site4)
colnames(BLUP_Of_Final1)[1] = "self"; colnames(BLUP_Of_Final1)[2] = "site"
head(BLUP_Of_Final1)

BLUP_Of_Final1 = as.data.frame(BLUP_Of_Final1)

BLUP_Of_FinalS1 = subset(BLUP_Of_Final1, site == 1)[,c(1:2,3:8)]
BLUP_Of_FinalS2 = subset(BLUP_Of_Final1, site == 2)[,c(1:2,9:14)]
BLUP_Of_FinalS3 = subset(BLUP_Of_Final1, site == 3)[,c(1:2,15:20)]
BLUP_Of_FinalS4 = subset(BLUP_Of_Final1, site == 4)[,c(1:2,21:26)]
head(BLUP_Of_FinalS1); head(BLUP_Of_FinalS2); head(BLUP_Of_FinalS3); head(BLUP_Of_FinalS4)

names(BLUP_Of_FinalS1) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Of_FinalS2) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Of_FinalS3) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")
names(BLUP_Of_FinalS4) <- c("self","site","BV_DBH30","Acc_DBH30","BV_HT30","Acc_HT30","BV_WGR36","Acc_WGR36")

BLUP_Of_Final = rbind(BLUP_Of_FinalS1,BLUP_Of_FinalS2,BLUP_Of_FinalS3,BLUP_Of_FinalS4)
write.table(BLUP_Of_Final, file = "HBLUP_BVs_ACCs_offspring.txt", row.names = FALSE)


############################################################################################################
# 5.1.4. Estimated provenance effects (BLUEs)
(BLUEs_H_proc = res.spa.H$fixed$proc)

write.table(BLUEs_H_proc, file = "HBLUEs_proc.txt")
