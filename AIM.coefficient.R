#=====alternative - selection of ancestry informative markers (AIM) (Rosenberg et al. 2003)==============#
# Copyright @ Blaise Ratcliffe (b.ratcliffe@gmail.com)
# With advice of Dr. Jaroslav Klapste#

# GEN11 is an n by p matrix of samples with SNP genotypes coded as 0, 0.5, 1
# fam$ID is a vector of genetic groups of length n  
# AIM1 is an p x 2 matrix of SNP indices and AIM coefficients, respectively.

GEN12<-cbind(fam$ID,GEN11)  #add identifier of genetic group to genotypes

Pij<-NULL

for(j in 1:max(fam$ID)){
  
  print(j)
  
  GEN13<-GEN12[GEN12[,1]==j,]
  
  if(is.null(dim(GEN13))){
    CM<-GEN13
    Pij<-rbind(Pij,CM)
    
  }
  
  if(!is.null(dim(GEN13))){
    CM<-colMeans(GEN13)
    Pij<-rbind(Pij,CM)
  }
}

Pij[Pij==0]<-0.0000000001  #modify allelic frequency to be non-zero for fixed alleles at particular genetic group 
Pij[Pij==1]<-0.9999999999

AIM<-matrix(0,ncol(snp.aim1),2)
AIM[,1]<-c(1:ncol(snp.aim1))

for(k in 1:ncol(snp.aim1)){
  A1<-(-mean(GEN12[,(k+1)])*log(mean(GEN12[,(k+1)])))+sum((Pij[,(k+1)]/max(fam$ID))*log(Pij[,(k+1)]))
  
  A2<-(-(1-mean(GEN12[,(k+1)]))*log(1-mean(GEN12[,(k+1)])))+sum(((1-Pij[,(k+1)])/max(fam$ID))*log(1-Pij[,(k+1)]))
  AIM[k,2]<-A1+A2
}

AIM1<-as.data.frame(AIM[order(-AIM[,2]),],stringsAsFactors = F)
