###################################### 
setwd("~/Documents/HTP Workshop/Exercises")
library(asreml)
library(rrBLUP)
library(plyr)
asreml.options(maxit=500)

#function to refit a model until it converges
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-update(mod)
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

##read in data
datA<- read.csv('NDVIgf.EYTBWB5IR_01.blues.csv', row.names=1)
datB<- read.csv('NDVIveg.EYTBWB5IR_01.blues.csv', row.names=1)
datC<- read.csv('gryld.EYTBWB5IR_01.blues.csv', row.names=1)
datD<- read.csv('NDVIgf.EYTBWBHT_01.blues.csv', row.names=1)
datE<- read.csv('NDVIveg.EYTBWBHT_01.blues.csv', row.names=1)
datF<- read.csv('gryld.EYTBWBHT_01.blues.csv', row.names=1)

#combine into one data frame
dat<- rbind(datA, datB, datC, datD, datE, datF)

#change gid to a factor
dat$gid<- as.factor(as.character(dat$gid))

#order, trait within trial
dat<- dat[order(dat$trait),]
dat<- dat[order(dat$trial),]

#get the marker relationship matrix
load('Gmat.RData')

#keep the gids that are present in our dataset
gidset<- intersect(dat$gid, row.names(Gmat))
ixG<- match(gidset, row.names(Gmat))
Gmatsub<- Gmat[ixG, ixG]

#make matrix positive semi-definite
#Gmatsub<- as.matrix(nearPD(Gmatsub)$mat)

#remove any gids from phenotype data
#that are not in the marker data
dat<- droplevels.data.frame(dat[which(dat$gid %in% row.names(Gmatsub)),])

#sort data by trait and then trial
dat<- dat[order(dat$trait),]
dat<- dat[order(dat$gid),]
dat<- dat[order(dat$trial),]

#fit mixed model with genotypes as random, separate genetic variance per trait
wt<- 1/dat$std.error^2
mod<- asreml(fixed= predicted.value~1+trait+ trait:trial, 
             random= ~us(trait):vm(gid, Gmatsub), data=dat,
             family = asr_gaussian(dispersion = 1),
             #residual= ~dsum(~id(gid):us(trait)| trial),
             weights= wt,
             na.action = na.method(y='include', x='include'))

