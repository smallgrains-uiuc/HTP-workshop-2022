setwd("~/Documents/HTP Workshop/Exercises")
library(asreml)

#read in the phenotypic data
dat<- read.csv('pheno.csv', row.names=1)

#look at the first six rows of the data
head(dat)

#display the structure of the dat object
str(dat)

#change variables to factors
dat$gid<- as.factor(dat$gid)
dat$rep<- as.factor(dat$rep)
dat$block<- as.factor(dat$block)
dat$row<- as.factor(dat$row)
dat$col<- as.factor(dat$col)

#Select a single trial
dat1<- dat[which(dat$trial=='EYTBWBHT_01'),]

#drop excess factor levels
dat1<- droplevels.data.frame(dat1)

############################################
#  Trait with one phenotyping date        #
############################################

#select a single trait
gy<- dat1[which(dat1$trait_id=='GRYLD'),]

#fit univariate model with genotypes as random and iid
mod<- asreml(fixed= phenotype_value~1+rep, 
             random= ~gid+block:rep, data=gy)

#look at the variance component estimates
summary(mod)$varcomp

#get the genetic variance from the model summary
Vg<- summary(mod)$varcomp['gid','component']

#generate the predicted values for gid, these are the BLUPs
blups<- predict(mod, classify='gid', ignore=c('(Intercept)', 'rep'))$pvals

#get the prediction error variance 
pev<- blups$std.error^2

#calcuate the reliability
#reliability is the square of the selection accuracy 
#this is equivalent to the 'line mean heritability'
rel<- 1-pev/Vg
mean(rel)

#now lets fit the model with genotypes as fixed
#this will let us get BLUEs that we could use 
#later in another model
mod<- asreml(fixed= phenotype_value~1+rep+gid, 
             random= ~block:rep, data=gy)
mod<- update(mod)

#generate the predicted values for gid, these are the BLUEs
blues<- predict(mod, classify='gid')$pvals

#look at the BLUEs
head(blues)

#add metadata to the BLUEs data.frame
blues<- data.frame(trait='GRYLD', trial='EYTBWBHT_01', blues)

#write the results to a csv file
#we will use this later on
write.csv(blues, file='gryld.EYTBWBHT_01.blues.csv')

############################################
#  Trait with multiple phenotyping dates   #
#      With analysis by growth stage       #
############################################
#select a single trait
ndvi1<- dat1[which(dat1$trait_id=='G-NDVI'),]

#select a stage
ndvi1_veg<- ndvi1[which(ndvi1$stage=='Veg'),]
#In this case, there is only one phenotyping date

#fit univariate model with genotypes as random and iid
mod<- asreml(fixed= phenotype_value~1+rep, 
             random= ~gid+block:rep,
             data=ndvi1_veg)

#get the genetic variance from the model summary
Vg<- summary(mod)$varcomp['gid','component']

#generate the predicted values for gid, these are the BLUPs
blups<- predict(mod, classify='gid', ignore=c('(Intercept)', 'rep'))$pvals

#get the prediction error variance 
pev<- blups$std.error^2

#calcuate the reliability
rel<- 1-pev/Vg
mean(rel)

#now lets fit the model with genotypes as fixed
#this will let us get BLUEs that we could use 
#later in another model
mod<- asreml(fixed= phenotype_value~1+rep+gid, 
             random= ~block:rep,
             data=ndvi1_veg)
mod<- update(mod)

#generate the predicted values for gid, these are the BLUEs
blues<- predict(mod, classify='gid')$pvals

#add metadata to the BLUEs data.frame
blues<- data.frame(trait='NDVIveg', trial='EYTBWBHT_01', blues)

#write the results to a csv file
#we will use this later on
write.csv(blues, file='NDVIveg.EYTBWBHT_01.blues.csv')

####Repeat for the grain filling stage
#select a stage
ndvi1_gf<- ndvi1[which(ndvi1$stage=='GF'),]

#fit univariate model with genotypes as random and iid
mod<- asreml(fixed= phenotype_value~1, 
             random= ~gid+phenotype_date+gid:phenotype_date+
               rep:phenotype_date+ block:rep:phenotype_date,
             data=ndvi1_gf)
mod<- update(mod)

#get the genetic variance from the model summary
Vg<- summary(mod)$varcomp['gid','component']

#generate the predicted values for gid, these are the BLUPs
blups<- predict(mod, classify='gid', ignore='(Intercept)')$pvals

#get the prediction error variance 
pev<- blups$std.error^2

#calcuate the reliability
rel<- 1-pev/Vg
mean(rel)

#now lets fit the model with genotypes as fixed
#this will let us get BLUEs that we could use 
#later in another model
mod<- asreml(fixed= phenotype_value~1+gid, 
             random= ~phenotype_date+gid:phenotype_date+
               rep:phenotype_date+ block:rep:phenotype_date,
              data=ndvi1_gf)

#generate the predicted values for gid, these are the BLUEs
blues<- predict(mod, classify='gid')$pvals

#add metadata to the BLUEs data.frame
blues<- data.frame(trait='NDVIgf', trial='EYTBWBHT_01', blues)

#write the results to a csv file
#we will use this later on
write.csv(blues, file='NDVIgf.EYTBWBHT_01.blues.csv')

