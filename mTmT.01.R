###################################### 
setwd("~/Documents/HTP Workshop/Exercises")
library(asreml)
library(rrBLUP)
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
datA<- read.csv('multiTrait.EYTBWBHT_01.blues.csv', row.names=1)
datB<- read.csv('multiTrait.EYTBWB5IR_01.blues.csv', row.names=1)

#combine into one data frame
dat<- rbind(datA, datB)

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

#fit mixed model with genotypes as random, separate genetic variance per trait
wt<- 1/dat$std.error^2
mod<- asreml(fixed= predicted.value~1+trait_id+ trait_id:trial, 
              random= ~us(trait_id):vm(gid, Gmatsub), data=dat,
              family = asr_gaussian(dispersion = 1),
              weights= wt,
              na.action = na.method(y='include', x='include'))

summary(mod)$varcomp

###make the covariance matrix to get the
###genetic correlations between traits

#first get vector of the covariances
trtnm<- as.character(unique(dat$trait_id))
ntrait<- length(trtnm)
vc<- summary(mod)$varcomp
ncomb<- ntrait*ntrait-1
vc_vec<- vc[1:c(ncomb),'component']
names(vc_vec)<- gsub('trait_id:vm(gid, Gmatsub)!trait_id_', "", 
                 row.names(vc)[1:c(ncomb)], fixed=TRUE)

#assemble into a matrix
Gcov<- matrix(NA, ntrait, ntrait)
colnames(Gcov)<- trtnm
rownames(Gcov)<- trtnm
for(i in 1:nrow(Gcov)){
  for(j in 1:ncol(Gcov)){
    Gcov[i,j]<- vc_vec[paste(trtnm[i], trtnm[j], sep=":")]
  }
}
Gcov[upper.tri(Gcov)]<- Gcov[lower.tri(Gcov)]

#get correlation matrix
cov2cor(Gcov)

#get the blups for yield
pred<- predict(mod, classify='gid:trait_id')
blups<- pred$pvals
blupsGy<- blups[which(blups$trait=='GRYLD'),]

