setwd("~/Documents/HTP Workshop/Exercises")
library(asreml)
library(reshape)

#set asreml options
asreml.options(maxit=500)

#function to refit a model until it converges
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-update(mod)
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

dat<- read.csv('pheno.csv', row.names=1)

#change variables to factors
dat$gid<- as.factor(dat$gid)
dat$rep<- as.factor(dat$rep)
dat$block<- as.factor(dat$block)
dat$row<- as.factor(dat$row)
dat$col<- as.factor(dat$col)
dat$plot_no<- as.factor(dat$plot_no)

#Select a single trial
dat1<- dat[which(dat$trial=='EYTBWB5IR_01'),]
dat1<- droplevels.data.frame(dat1)

##########################################################
#  Traits with multiple dates and a trait with one date  #
##########################################################
#select which traits
sub<- droplevels.data.frame(dat1[which(dat1$trait_id %in% c('GRYLD', 'G-NDVI')),])

#for NDVI make the trait name different according to growth stage
sub$trait_id<- as.character(sub$trait_id) #change trait_id to character
sub[which(sub$trait_id=='G-NDVI'),'trait_id']<- 
  paste(sub[which(sub$trait_id=='G-NDVI'),'trait_id'], 
        sub[which(sub$trait_id=='G-NDVI'),'stage'],
        sep="")
sub$trait_id<- as.factor(sub$trait_id)#change trait_id back to factor

#create a trait-date variable
sub$trait_id2<- as.factor(paste(sub$trait_id, sub$phenotype_date, sep="_"))

#create a data frame with missing values in each row
datNA<- melt(cast(sub, plot_no~trait_id2, value='phenotype_value'))

#merge data with missing values with the original data.frame
#this is needed because rows of missing data are needed
#missing data should not be omitted, it needs to be there as NA
sub<- merge(sub, datNA, by=c('plot_no', 'trait_id2'), all=TRUE)

#check what traits were measured on what dates
smry<- cast(sub, trait_id~phenotype_date)
smry#notice that G-NDVIVeg was only measured on one date

#fit multivariate model with genotypes as random and iid
#this includes conditional factors for traits with multiple dates
#it assumes there is gid covariance between traits
#it assumes there is a residual covariance between all trait-date combinations
mod<- asreml(fixed= phenotype_value~1+trait_id+at(trait_id, 'GRYLD'):rep, 
             random= ~us(trait_id):gid+
               at(trait_id, c('G-NDVIVeg', 'G-NDVIGF')):phenotype_date+
               at(trait_id, c('G-NDVIVeg', 'G-NDVIGF')):gid:phenotype_date+
               at(trait_id, c('G-NDVIVeg', 'G-NDVIGF')):rep:phenotype_date+ 
               at(trait_id, c('G-NDVIVeg', 'G-NDVIGF')):block:rep:phenotype_date+
               at(trait_id, 'GRYLD'):block:rep,
             residual= ~id(plot_no):us(trait_id2),
             data=sub, na.action = na.method(y='include', x='include'))

#look at the variance component table
summary(mod)$varcomp

#fit multivariate model with genotypes as fixed
#this includes conditional factors for traits with multiple dates
#it fits a separate gid effect per trait
#it assumes there is a residual covariance between all trait-date combinations
mod<- asreml(fixed= phenotype_value~1+trait_id+at(trait_id, 'GRYLD'):rep+ 
               at(trait_id):gid, 
             random= ~at(trait_id, c('G-NDVIVeg', 'G-NDVIGF')):phenotype_date+
               at(trait_id, c('G-NDVIVeg', 'G-NDVIGF')):gid:phenotype_date+
               at(trait_id, c('G-NDVIVeg', 'G-NDVIGF')):rep:phenotype_date+ 
               at(trait_id, c('G-NDVIVeg', 'G-NDVIGF')):block:rep:phenotype_date+
               at(trait_id, 'GRYLD'):block:rep,
             residual= ~id(plot_no):us(trait_id2),
             data=sub, na.action = na.method(y='include', x='include'))

#look at the variance component table
summary(mod)$varcomp

#generate the predicted values for gid, these are the BLUPs
blues<- predict(mod, classify='gid:trait_id')$pvals

#add the trial name to the blues data.frame
blues<- data.frame(trial='EYTBWB5IR_01', blues)

#write the results to a file
write.csv(blues, file= "multiTrait.EYTBWB5IR_01.blues.csv")


