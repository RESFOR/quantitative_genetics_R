# Title: EXPLORATORY DATA ANALYSIS (EDA): an example using three traits (DBH30, HT30, WGR36) 
#        and 4 sites (JUDY, VIRG, SWAN, and TIME) from the RES-FOR lodgepole pine dataset.
# Author: Eduardo Pablo Cappa (cappa.eduardo@inta.gob.ar)
# Copyright @ Dr. Eduardo Pable Cappa
# Revised: "April 12th, 2019"
# Revised to publish: "May 28th, 2021"


# PACKAGES
################################
wants <- c("plyr", "fields","reshape2","psych","ggplot2","scatterplot3d")
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])
suppressPackageStartupMessages(sapply(wants, require, character.only=TRUE))

# 1. READ THE RAW DATA
data<-read.table(file="raw_data.txt", head=TRUE)

# Number of the column where start the fenotyic information 
start <- 16


##################################################################################
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

# Replicates nested in sites
data <- transform(data, sitrep= factor(site:rep))

# Sets nested in Replicates 
data <- transform(data, repset= factor(rep:set))

# Create the plot effects
#data <- transform(data, plot= factor(as.factor(mum):rep:set))

# Data visualization
head(data, n=10);tail(data, n=10)

# Data by site
data.s1=subset(data,data$site==1)
data.s2=subset(data,data$site==2)
data.s3=subset(data,data$site==3)
data.s4=subset(data,data$site==4)


#########################################################################################
# 3. CHOOSE THE TARGET TRAITs FOR THE EDA and SPLIT THE DATA FOR EACH SITE
# NSWGR3604 is the transformate WGR36 data into "normal score".
# See details of Normal Score in 
# Across sites
data$phe1<-data$DBH30
data$phe2<-data$HT30
data$phe3<-data$WGR36_04
data$phe4<-data$NSWGR3604

# For each site
data.s1$phe1<-data.s1$DBH30
data.s1$phe2<-data.s1$HT30
data.s1$phe3<-data.s1$WGR36_04
data.s1$phe4<-data.s1$NSWGR3604

data.s2$phe1<-data.s2$DBH30
data.s2$phe2<-data.s2$HT30
data.s2$phe3<-data.s2$WGR36_04
data.s2$phe4<-data.s2$NSWGR3604

data.s3$phe1<-data.s3$DBH30
data.s3$phe2<-data.s3$HT30
data.s3$phe3<-data.s3$WGR36_04
data.s3$phe4<-data.s3$NSWGR3604

data.s4$phe1<-data.s4$DBH30
data.s4$phe2<-data.s4$HT30
data.s4$phe3<-data.s4$WGR36_04
data.s4$phe4<-data.s4$NSWGR3604


########################################################################################
# 4. CHECKING INCONSISTENCIES IN THE DESIGN
# Number of total trees by replicates and site
(Total_tree = table(data[ ,c("rep","site")]))  

# Number of trees by set and site
table(data[ ,c("set","site")])  

# Number of trees by plot and site
head(table(data[ ,c("plot","site")]), n=30)  

# Number of flias (mum) by site
table(data[ ,c("mum","site")])  


########################################################################################
# 5. STUDY of SURVIVING TREES
# Number (and percentages) of surviving trees using recorded trees for DBH30
SV=as.data.frame(cbind(data$mum, as.factor(data$proc), as.factor(data$site), as.factor(data$rep), data$DBH30))
SV=na.omit(SV);colnames(SV) = c("mum","proc","site","rep","DBH30");SV$site=as.factor(SV$site);SV$proc=as.factor(SV$proc);SV$rep=as.factor(SV$rep)
SV$mum=as.factor(SV$mum)
head(SV)

# Number of survival trees by site
(SURV.s= table(SV[ ,c("site")])) 
# Survival rate (percentages) by site
(round((SURVAL.s = SURV.s/colSums(Total_tree)),3))

# Survival rate (percentages) by replicates and site
round((SURVAL = SURV/Total_tree),2) 

# Survival rate (percentages) by provenance and site
round(table(SV[ ,c("proc","site")])/table(data[ ,c("proc","site")]),3) 

# Number of survival trees by families (mum) and site
round(table(SV[ ,c("mum","site")])/table(data[ ,c("mum","site")]),3)  


#########################################################################################
# 6. SOME STATISTICS OF THE DATA
# 6.1. FOR CONTINOUS TRAITS 
# For all data and sites
describe(data)

# For the phenotypic data by site, replicates and provenance
describeBy(data$phe1, data$site,mat=T)[,-c(1,8)]
describeBy(data$phe2, data$site,mat=T)[,-c(1,8)]
describeBy(data$phe3, data$site,mat=T)[,-c(1,8)]
describeBy(data$phe4, data$site,mat=T)[,-c(1,8)]

describeBy(data$phe1, data$sitrep,mat=T)[,-c(1,8)]
describeBy(data$phe2, data$sitrep,mat=T)[,-c(1,8)]
describeBy(data$phe3, data$sitrep,mat=T)[,-c(1,8)]
describeBy(data$phe4, data$sitrep,mat=T)[,-c(1,8)]

describeBy(data$phe1, data$proc,mat=T)[,-c(1,8)]
describeBy(data$phe2, data$proc,mat=T)[,-c(1,8)]
describeBy(data$phe3, data$proc,mat=T)[,-c(1,8)]
describeBy(data$phe4, data$proc,mat=T)[,-c(1,8)]

# Phenotypic correlations and distribution for all traits and sites
pairs.panels(data[(start):(start+9)])
# by site
pairs.panels(data.s1[,(start+1):(start+8)], cex.cor=1)
pairs.panels(data.s2[,(start+1):(start+8)], cex.cor=1)
pairs.panels(data.s3[,(start+1):(start+8)], cex.cor=1)
pairs.panels(data.s4[,(start+1):(start+8)], cex.cor=1)

# Correlations and distribution for the target traits by site
pairs.panels(data.s1[,c("phe1","phe2","phe4")], main="Site JUDY", cex.cor=0.6)
pairs.panels(data.s2[,c("phe1","phe2","phe4")], main="Site VIRG", cex.cor=0.6)
pairs.panels(data.s3[,c("phe1","phe2","phe4")], main="Site SWAN", cex.cor=0.6)
pairs.panels(data.s4[,c("phe1","phe2","phe4")], main="Site TIME", cex.cor=0.6)

# Correaltions for the target traits
print("Correlations between the target traits across all sites")
round(cor(data[,c("phe1","phe2","phe3","phe4")], use = "na.or.complete"),3)
print("Correlations between the target traits for JUDY site")
round(cor(data.s1[,c("phe1","phe2","phe3","phe4")], use = "na.or.complete"),3)
print("Correlations between the target traits for VIRG site")
round(cor(data.s2[,c("phe1","phe2","phe3","phe4")], use = "na.or.complete"),3)
print("Correlations between the target traits for SWAN site")
round(cor(data.s3[,c("phe1","phe2","phe3","phe4")], use = "na.or.complete"),3)
print("Correlations between the target traits for TIME site")
round(cor(data.s3[,c("phe1","phe2","phe3","phe4")], use = "na.or.complete"),3)

# 6.2. FOR binary or categorical TRAITs
# Histograms for the WGR36 trait (4-point score)
# All sites
par(mfrow=c(1,1))
aa=barplot(table(data$phe3), ylim = c(0,6000), col = "blue", main="WGR36 for alll sites")
text(x = aa, y = table(data$phe3), label = table(data$phe3), pos = 3, cex = 0.8, col = "red")

# by sites
par(mfrow=c(2,2))
a=barplot(table(data.s1$phe3), col = "blue", main="WGR36 for JUDY")
text(x = a, y = table(data.s1$phe3), label = table(data.s1$phe3), pos = 3, cex = 0.8, col = "red")
b=barplot(table(data.s2$phe3), col = "blue", main="WGR36 for VIRG")
text(x = b, y = table(data.s2$phe3), label = table(data.s2$phe3), pos = 3, cex = 0.8, col = "red")
c=barplot(table(data.s3$phe3), col = "blue", main="WGR36 for SWAN")
text(x = c, y = table(data.s3$phe3), label = table(data.s3$phe3), pos = 3, cex = 0.8, col = "red")
d=barplot(table(data.s4$phe3), col = "blue", main="WGR36 for TIME")
text(x = d, y = table(data.s4$phe3), label = table(data.s4$phe3), pos = 3, cex = 0.8, col = "red")

# Number and percentages of trees in each category by sites
# Across all sites by site
par(mfrow=c(1,1))
(c.s=table(data[ , c("phe3","site")]))
(c.sp=round(prop.table(c.s, 2),3))
a=barplot(c.s , ylim= c(0,8000), ylab="number of observation",xlab= "sites" , main = "WGR36 by site", col = c("blue","yellow","red","black"))
legend('top', ncol = 5L, cex = 0.7, title = "Percentages by Sites",box.lty=0,
       legend = c('', '1', '2', '3', '4','JUDY', c.sp[1:4,1], 'VIRG', c.sp[1:4,2], 'SWAN', c.sp[1:4,3], 'TIME', c.sp[1:4,4]))

# Across all sites by replicates nested in site
par(mfrow=c(1,1))
(c.s2=table(data[ , c("phe3","sitrep")]))
(c.s2p=round(prop.table(c.s2, 2),3))
c=barplot(c.s2, ylim= c(0,1000), ylab="number of observation",xlab= "replicates nested in sites" , main = "WGR36 by replicates nested in sites", col = c("blue","yellow","red","black","brown","orange","green"))

# Across all sites by provenance
par(mfrow=c(1,1))
data.na<-data[!is.na(data$phe3),]
(c.p=table(data.na[ , c("phe3","proc")]))
(c.pp=round(prop.table(c.p, 2),3))
b=barplot(c.p , ylim= c(0,5000), ylab="number of observation",xlab= "provenances" , main = "WGR36 by provenances", col = c("blue","yellow","red","black","brown","orange","green"))

# Across all sites by provenances in each site
(c.s3=table(data[ , c("phe3","proc","site")]))
c.s3=table(data[ , c("phe3","proc","site")])
(c.s23=round(prop.table(c.s3, 2),3))


################################################################################################################
# 7. CHECKING PHENOTYPES
# Plots for all sites by trait
par(mfrow=c(1,3))
plot(data$phe1,col=data$site,ylab = "DBH30", xlab = "trees")
plot(data$phe2,col=data$site, ylab = "HT30", xlab = "trees")
plot(data$phe3,col=data$site, ylab = "WGR36", xlab = "trees")
axis(2, at = c(0.5,0.5), las=2)
legend(100,0.9,c("Site JUDY","Site VIRG","Site SWAN","Site TIME"), fill=c("black","red","green","blue"))

# Some plots of DBH30 and HT30 by site
par(mfrow = c(2, 2),oma = c(0, 0, 2, 0))
plot(data.s1$phe1, main = "Distribution");hist(data.s1$phe1,main = "Histogram");qqnorm(data.s1$phe1);boxplot(data.s1$phe1, main= "Boxplot")
mtext("JUDY Site - DBH30", outer = TRUE, cex = 1.5)
plot(data.s1$phe2, main = "Distribution");hist(data.s1$phe2,main = "Histogram");qqnorm(data.s1$phe2);boxplot(data.s1$phe2, main= "Boxplot")
mtext("JUDY Site - HT30", outer = TRUE, cex = 1.5)

plot(data.s2$phe1, main = "Distribution");hist(data.s2$phe1,main = "Histogram");qqnorm(data.s2$phe1);boxplot(data.s2$phe1, main= "Boxplot")
mtext("VIRG Site - DBH30", outer = TRUE, cex = 1.5)
plot(data.s2$phe2, main = "Distribution");hist(data.s2$phe2,main = "Histogram");qqnorm(data.s2$phe2);boxplot(data.s2$phe2, main= "Boxplot")
mtext("VIRG Site - HT30", outer = TRUE, cex = 1.5)

plot(data.s3$phe1, main = "Distribution");hist(data.s3$phe1,main = "Histogram");qqnorm(data.s3$phe1);boxplot(data.s3$phe1, main= "Boxplot")
mtext("SWAN Site - DBH30", outer = TRUE, cex = 1.5)
plot(data.s3$phe2, main = "Distribution");hist(data.s3$phe2,main = "Histogram");qqnorm(data.s3$phe2);boxplot(data.s3$phe2, main= "Boxplot")
mtext("SWAN Site - HT30", outer = TRUE, cex = 1.5)

plot(data.s4$phe1, main = "Distribution");hist(data.s4$phe1,main = "Histogram");qqnorm(data.s4$phe1);boxplot(data.s4$phe1, main= "Boxplot")
mtext("TIME Site - DBH30", outer = TRUE, cex = 1.5)
plot(data.s4$phe2, main = "Distribution");hist(data.s4$phe2,main = "Histogram");qqnorm(data.s4$phe2);boxplot(data.s4$phe2, main= "Boxplot")
mtext("TIME Site - HT30", outer = TRUE, cex = 1.5)

# Boxplots by provenances
# DBH30 and HT30
par(mfrow = c(2, 4))
boxplot(phe1 ~ proc, data = data.s1, col = rainbow(8), xlab = "Provenance", ylab = "DBH30", main="JUDY site")
boxplot(phe1 ~ proc, data = data.s2, col = rainbow(8), xlab = "Provenance", ylab = "DBH30", main="VIRG site")
boxplot(phe1 ~ proc, data = data.s3, col = rainbow(8), xlab = "Provenance", ylab = "DBH30", main="SWAN site")
boxplot(phe1 ~ proc, data = data.s4, col = rainbow(8), xlab = "Provenance", ylab = "DBH30", main="TIME site")
boxplot(phe2 ~ proc, data = data.s1, col = rainbow(8), xlab = "Provenance", ylab = "HT30")
boxplot(phe2 ~ proc, data = data.s2, col = rainbow(8), xlab = "Provenance", ylab = "HT30")
boxplot(phe2 ~ proc, data = data.s3, col = rainbow(8), xlab = "Provenance", ylab = "HT30")
boxplot(phe2 ~ proc, data = data.s4, col = rainbow(8), xlab = "Provenance", ylab = "HT30")


########################################################################################################################
# 8. PRELIMINAR SPATIAL DIAGNOSIS USING THE RAW PHENOTYPIC DATA
# Boxplots by replicates
par(mfrow = c(2, 4))
boxplot(phe1 ~ rep, data = data.s1, col = rainbow(8), xlab = "Replicate", ylab = "DBH30", main="JUDY site")
boxplot(phe1 ~ rep, data = data.s2, col = rainbow(8), xlab = "Replicate", ylab = "DBH30", main="VIRG site")
boxplot(phe1 ~ rep, data = data.s3, col = rainbow(8), xlab = "Replicate", ylab = "DBH30", main="SWAN site")
boxplot(phe1 ~ rep, data = data.s4, col = rainbow(8), xlab = "Replicate", ylab = "DBH30", main="TIME site")
boxplot(phe2 ~ rep, data = data.s1, col = rainbow(8), xlab = "Replicate", ylab = "HT30")
boxplot(phe2 ~ rep, data = data.s2, col = rainbow(8), xlab = "Replicate", ylab = "HT30")
boxplot(phe2 ~ rep, data = data.s3, col = rainbow(8), xlab = "Replicate", ylab = "HT30")
boxplot(phe2 ~ rep, data = data.s4, col = rainbow(8), xlab = "Replicate", ylab = "HT30")

# Plots of median of phenoypes along rows and column
# DBH30
ByRow_P<-describeBy(x=data.s1$phe1,group=data.s1$col, mat=TRUE)
ByCol_P<-describeBy(x=data.s1$phe1,group=data.s1$row, mat=TRUE)
ByRow_R<-describeBy(x=data.s1$phe1,group=data.s1$col, mat=TRUE)
ByCol_R<-describeBy(x=data.s1$phe1,group=data.s1$row, mat=TRUE)
par(mfrow = c(2, 2))
plot(as.numeric(ByRow_P$group1), ByRow_P$median, pch=19)
plot(as.numeric(ByCol_P$group1), ByCol_P$median, pch=19)
plot(as.numeric(ByRow_R$group1), ByRow_R$median, pch=19)
plot(as.numeric(ByCol_R$group1), ByCol_R$median, pch=19)
# HT30
ByRow_P<-describeBy(x=data.s1$phe2,group=data.s1$col, mat=TRUE)
ByCol_P<-describeBy(x=data.s1$phe2,group=data.s1$row, mat=TRUE)
ByRow_R<-describeBy(x=data.s1$phe2,group=data.s1$col, mat=TRUE)
ByCol_R<-describeBy(x=data.s1$phe2,group=data.s1$row, mat=TRUE)
par(mfrow = c(2, 2))
plot(as.numeric(ByRow_P$group1), ByRow_P$median, pch=19)
plot(as.numeric(ByCol_P$group1), ByCol_P$median, pch=19)
plot(as.numeric(ByRow_R$group1), ByRow_R$median, pch=19)
plot(as.numeric(ByCol_R$group1), ByCol_R$median, pch=19)

# Heatmaps (spatial patterns) of the phenotpic data: DBH30, HT30, and WGR36 for each of the 4 lodgepole sites 
par(mfrow = c(3, 4))
data.s1$rep1=as.integer(data.s1$rep)
ms1.bl=as.matrix(acast(data.s1, data.s1$row~data.s1$col, value.var="rep1"))
ms1.phe1=acast(data.s1, data.s1$row~data.s1$col, value.var="phe1")
image.plot(x=1:max(data.s1$row),y=1:max(data.s1$col),z=ms1.phe1,col=rainbow(100), main = "Site JUDY DBH30")
contour(1:max(data.s1$row), 1:max(data.s1$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s2$rep1=as.integer(data.s2$rep)
ms1.bl=as.matrix(acast(data.s2, data.s2$row~data.s2$col, value.var="rep1"))
ms1.phe1=acast(data.s2, data.s2$row~data.s2$col, value.var="phe1")
image.plot(x=1:max(data.s2$row),y=1:max(data.s2$col),z=ms1.phe1,col=rainbow(100), main = "Site VIRG DBH30")
contour(1:max(data.s2$row), 1:max(data.s2$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s3$rep1=as.integer(data.s3$rep)
ms1.bl=as.matrix(acast(data.s3, data.s3$row~data.s3$col, value.var="rep1"))
ms1.phe1=acast(data.s3, data.s3$row~data.s3$col, value.var="phe1")
image.plot(x=1:max(data.s3$row),y=1:max(data.s3$col),z=ms1.phe1,col=rainbow(100), main = "Site SWAN DBH30")
contour(1:max(data.s3$row), 1:max(data.s3$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s4$rep1=as.integer(data.s4$rep)
ms1.bl=as.matrix(acast(data.s4, data.s4$row~data.s4$col, value.var="rep1"))
ms1.phe1=acast(data.s4, data.s4$row~data.s4$col, value.var="phe1")
image.plot(x=1:max(data.s4$row),y=1:max(data.s4$col),z=ms1.phe1,col=rainbow(100), main = "Site TIME DBH30")
contour(1:max(data.s4$row), 1:max(data.s4$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)

data.s1$rep1=as.integer(data.s1$rep)
ms1.bl=as.matrix(acast(data.s1, data.s1$row~data.s1$col, value.var="rep1"))
ms1.phe2=acast(data.s1, data.s1$row~data.s1$col, value.var="phe2")
image.plot(x=1:max(data.s1$row),y=1:max(data.s1$col),z=ms1.phe2,col=rainbow(100), main = "Site JUDY HT30")
contour(1:max(data.s1$row), 1:max(data.s1$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s2$rep1=as.integer(data.s2$rep)
ms1.bl=as.matrix(acast(data.s2, data.s2$row~data.s2$col, value.var="rep1"))
ms1.phe2=acast(data.s2, data.s2$row~data.s2$col, value.var="phe2")
image.plot(x=1:max(data.s2$row),y=1:max(data.s2$col),z=ms1.phe2,col=rainbow(100), main = "Site VIRG HT30")
contour(1:max(data.s2$row), 1:max(data.s2$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s3$rep1=as.integer(data.s3$rep)
ms1.bl=as.matrix(acast(data.s3, data.s3$row~data.s3$col, value.var="rep1"))
ms1.phe2=acast(data.s3, data.s3$row~data.s3$col, value.var="phe2")
image.plot(x=1:max(data.s3$row),y=1:max(data.s3$col),z=ms1.phe2,col=rainbow(100), main = "Site SWAN HT30")
contour(1:max(data.s3$row), 1:max(data.s3$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s4$rep1=as.integer(data.s4$rep)
ms1.bl=as.matrix(acast(data.s4, data.s4$row~data.s4$col, value.var="rep1"))
ms1.phe2=acast(data.s4, data.s4$row~data.s4$col, value.var="phe2")
image.plot(x=1:max(data.s4$row),y=1:max(data.s4$col),z=ms1.phe2,col=rainbow(100), main = "Site TIME HT30")
contour(1:max(data.s4$row), 1:max(data.s4$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)

data.s1$rep1=as.integer(data.s1$rep)
ms1.bl=as.matrix(acast(data.s1, data.s1$row~data.s1$col, value.var="rep1"))
ms1.phe3=acast(data.s1, data.s1$row~data.s1$col, value.var="phe3")
image.plot(x=1:max(data.s1$row),y=1:max(data.s1$col),z=ms1.phe3, zlim= c(1,7), col=c("blue","yellow","red","black","brown","orange","green"), main = "Site JUDY WGR36")
contour(1:max(data.s1$row), 1:max(data.s1$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s2$rep1=as.integer(data.s2$rep)
ms1.bl=as.matrix(acast(data.s2, data.s2$row~data.s2$col, value.var="rep1"))
ms1.phe3=acast(data.s2, data.s2$row~data.s2$col, value.var="phe3")
image.plot(x=1:max(data.s2$row),y=1:max(data.s2$col),z=ms1.phe3, zlim= c(1,7), col=c("blue","yellow","red","black","brown","orange","green"), main = "Site VIRG WGR36")
contour(1:max(data.s2$row), 1:max(data.s2$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s3$rep1=as.integer(data.s3$rep)
ms1.bl=as.matrix(acast(data.s3, data.s3$row~data.s3$col, value.var="rep1"))
ms1.phe3=acast(data.s3, data.s3$row~data.s3$col, value.var="phe3")
image.plot(x=1:max(data.s3$row),y=1:max(data.s3$col),z=ms1.phe3, zlim= c(1,7), col=c("blue","yellow","red","black","brown","orange","green"), main = "Site SWAN WGR36")
contour(1:max(data.s3$row), 1:max(data.s3$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)
data.s4$rep1=as.integer(data.s4$rep)
ms1.bl=as.matrix(acast(data.s4, data.s4$row~data.s4$col, value.var="rep1"))
ms1.phe3=acast(data.s4, data.s4$row~data.s4$col, value.var="phe3")
image.plot(x=1:max(data.s4$row),y=1:max(data.s4$col),z=ms1.phe3, zlim= c(1,7), col=c("blue","yellow","red","black","brown","orange","green"), main = "Site TIME WGR36")
contour(1:max(data.s4$row), 1:max(data.s4$col), ms1.bl, col = "black", add = TRUE, lwd= 3, nlevel=6)

# 3D Distribution of the target DBH30 and HT30 traits
par(mfrow=c(1,1))
scatterplot3d(data.s1$row,data.s1$col,data.s1$phe1, pch=16, highlight.3d=TRUE,type="h",angle = 30,main = "Site JUDY DBH30")
scatterplot3d(data.s1$row,data.s1$col,data.s1$phe2, pch=16, highlight.3d=TRUE,type="h",angle = 30,main = "Site JUDY HT30")
scatterplot3d(data.s2$row,data.s2$col,data.s2$phe1, pch=16, highlight.3d=TRUE,type="h",angle = 30,main = "Site VIRG DBH30")
scatterplot3d(data.s2$row,data.s2$col,data.s2$phe2, pch=16, highlight.3d=TRUE,type="h",angle = 30,main = "Site VIRG HT30")
scatterplot3d(data.s3$row,data.s3$col,data.s3$phe1, pch=16, highlight.3d=TRUE,type="h",angle = 30,main = "Site SWAN DBH30")
scatterplot3d(data.s3$row,data.s3$col,data.s3$phe2, pch=16, highlight.3d=TRUE,type="h",angle = 30,main = "Site SWAN HT30")
scatterplot3d(data.s4$row,data.s4$col,data.s4$phe1, pch=16, highlight.3d=TRUE,type="h",angle = 30,main = "Site TIME DBH30")
scatterplot3d(data.s4$row,data.s4$col,data.s4$phe2, pch=16, highlight.3d=TRUE,type="h",angle = 30,main = "Site TIME HT30")
