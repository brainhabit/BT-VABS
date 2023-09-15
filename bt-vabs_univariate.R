###########################################################################################################################################
#
# Script to test univariate assumptions for twin analysis of the BT-VABS paper investigating
# etiological structure underlying different domains of the Vineland instrument.
#
#
# G. Bussu, 28 November 2022 [adapted from a script tutorial provided by M. Taylor (Feb 2022)]
#
####################################
# The script needs as input the selected dataset in the wide format (twin_data.csv from the bt-vabs_prepdata script).
#
# The script gives as output results from the univariate model on Motor or Social-Communication scales compared to the saturated model,
# including tests of twin modelling assumptions.
#
###########################################################################################################################################

rm(list=ls())

### set a working directory (otherwise create a new R project):

setwd('C:/Users/giobu365/Documents/BT-VABS')
list.files()

### load OpenMx and functions:

require(OpenMx)
source('C:/Users/giobu365/Documents/twin_modelling/tutorial/miFunctions.R')

### prepare some data:

data <- read.csv(file='C:/Users/giobu365/Documents/BT-VABS/data/twin_data_VABS_NAincome.csv',header=T,sep=',')
names(data);dim(data)

### select variables for analysis:

var <- 'resSocComm' # name of variable
nv <- 1 # number of phenotypes
ntv <- nv*2 # total number of variables
(selVars <- paste(var,c(rep(1,nv),rep(2,nv)),sep=''))

# makes a list of
# variables as they appear 
# in the file
# OBS: all variables for twin 1 are specified first, followed by twin 2

### split the data into subsets by zygosity:

mz <- subset(data,zygosity=='MZ',selVars)
dz <- subset(data,zygosity=='DZ',selVars)

#=====
# Assumptions testing
#=====

### specify and fit a fully saturated model

# estimate the means (separately for each twin and zygosity group):

expMeanMZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.01,
                      labels=labFull('m_mz',1,ntv),name='ExpMeanMZ')
expMeanDZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.01,
                      labels=labFull('m_dz',1,ntv),name='ExpMeanDZ')

# estimate the variances and covariances:

expCovMZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=c(1,.5,1),
                     labels=labSymm('mz_cov',ntv),name='ExpCovMZ')
expCovDZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=c(1,.5,1),
                     labels=labSymm('dz_cov',ntv),name='ExpCovDZ')

# convert the covariances to correlations:

matI <- mxMatrix(type='Iden',nrow=ntv,ncol=ntv,name='I')

expCorMZ <- mxAlgebra(solve(sqrt(I*ExpCovMZ))%&%ExpCovMZ,name='ExpCorMZ')
expCorDZ <- mxAlgebra(solve(sqrt(I*ExpCovDZ))%&%ExpCovDZ,name='ExpCorDZ')

#expCorMZ <- mxAlgebra(solve(sqrt(I*ExpCovMZ))%*%ExpCovMZ%*%solve(sqrt(I*ExpCovMZ)),name='ExpCorMZ')  
#expCorDZ <- mxAlgebra(solve(sqrt(I*ExpCovDZ))%*%ExpCovDZ%*%solve(sqrt(I*ExpCovDZ)),name='ExpCorDZ')



# tell OpenMx where to find the observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# specify the objectives of the model (ie to estimate the means and
# covariances):

objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMeanMZ',
                             dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMeanDZ',
                             dimnames=selVars)

# specify that we will use maximum-likelihood to estimate parameters:

funcML <- mxFitFunctionML()

# specify which 'groups' are included in the model, and which parameters
# are to be estimated in each group:

modelMZ <- mxModel('MZ',expMeanMZ,expCovMZ,matI,expCorMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',expMeanDZ,expCovDZ,matI,expCorDZ,dataDZ,objDZ,funcML)

# combine all of the groups into one single objective:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

# specify which parameters Mx should estimate confidence intervals around:

ci <- mxCI(c('MZ.ExpMeanMZ','DZ.ExpMeanDZ',
             'MZ.ExpCovMZ','DZ.ExpCovDZ','MZ.ExpCorMZ','DZ.ExpCorDZ'))

# combine all of the model objects together into one model:

SatModel <- mxModel('Sat',modelMZ,modelDZ,multi,ci)

# and fit the model:

SatFit <- mxTryHard(SatModel,intervals=T)
summary(SatFit)

### check the results against observed statistics:

# put all the estimates from the model in an object:

est <- round(SatFit$output$confidenceIntervals,2)

# means:

ObsMeanMZ <- round(colMeans(mz[,1:2],na.rm=T),3)
ObsMeanDZ <- round(colMeans(dz[,1:2],na.rm=T),3)

EstMeanMZ <- c(est[1,2],est[2,2])
EstMeanDZ <- c(est[3,2],est[4,2])

# variance and covariance:

ObsCovMZ <- round(cov(mz[,1:2],use='complete'),3)
ObsCovDZ <- round(cov(dz[,1:2],use='complete'),3)

EstCovMZ <- rbind(cbind(est[5,2],est[6,2]),
                  cbind(est[6,2],est[7,2]))
EstCovDZ <- rbind(cbind(est[8,2],est[9,2]),
                  cbind(est[9,2],est[10,2]))

# check the results up against the observed data:

ObsMeanMZ;EstMeanMZ
ObsMeanDZ;EstMeanDZ

ObsCovMZ;EstCovMZ
ObsCovDZ;EstCovDZ

### test assumptions by fitting nested models

# test the assumption of equal means across twin order:

SatModel2 <- mxModel(SatFit,name='Sat2')
SatModel2 <- omxSetParameters(SatModel2,labels=labFull('m_mz',1,ntv),
                              free=T,values=1,newlabel='m_mz')
SatModel2 <- omxSetParameters(SatModel2,labels=labFull('m_dz',1,ntv),
                              free=T,values=1,newlabel='m_dz')
SatFit2 <- mxTryHard(SatModel2,intervals=T)
summary(SatFit2)

mxCompare(SatFit,SatFit2)

# test assumption of equal means across zygosity:

SatModel3 <- mxModel(SatFit2,name='Sat3')
SatModel3 <- omxSetParameters(SatModel3,labels=c('m_mz','m_dz'),
                              free=T,values=1,newlabel='m')
SatFit3 <- mxTryHard(SatModel3,intervals=T)
summary(SatFit3)

# test assumption of equal variances across twin order:

SatModel4 <- mxModel(SatFit3,name='Sat4')
SatModel4 <- omxSetParameters(SatModel4,labels=labDiag('mz_cov',ntv),
                              free=T,values=1.2,newlabel='v_mz')
SatModel4 <- omxSetParameters(SatModel4,labels=labDiag('dz_cov',ntv),
                              free=T,values=1.2,newlabel='v_dz')
SatFit4 <- mxTryHard(SatModel4,intervals=T)
summary(SatFit4)

# test assumption of equal variances across zygosity:

SatModel5 <- mxModel(SatFit4,name='Sat5')
SatModel5 <- omxSetParameters(SatModel5,labels=c('v_mz','v_dz'),free=T,
                              values=1.2,newlabel='v')
SatFit5 <- mxTryHard(SatModel5,intervals=T)
summary(SatFit5)

### compare the fit of the models to one another

mxCompare(SatFit,c(SatFit2,SatFit3,SatFit4,SatFit5))

#########################################################################################
#=====
# ACE model
#=====

### means:

expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=1,
                    label='m',name='ExpMean')

### path coefficients:

pathA <- mxMatrix(type='Full',nrow=nv,ncol=nv,free=T,values=1,label='a_1_1',
                  name='a')
pathC <- mxMatrix(type='Full',nrow=nv,ncol=nv,free=T,values=1,label='c_1_1',
                  name='c')
pathE <- mxMatrix(type='Full',nrow=nv,ncol=nv,free=T,values=1,label='e_1_1',
                  name='e')

### calculate variance components:

varA <- mxAlgebra(a^2,name='A')
varC <- mxAlgebra(c^2,name='C')
varE <- mxAlgebra(e^2,name='E')

### calculate the total variance:

varP <- mxAlgebra(A+C+E,name='V')

### convert the variance components to proportions:

estVC <- mxAlgebra(cbind(A/V,C/V,E/V),name='EstVC')

### expected variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

### observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

### objectives of the model:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',
                             dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',
                             dimnames=selVars)
funcML <- mxFitFunctionML()

### data groups:

# make a list of objects that are common to all groups (to make 
# specifying each group easier):

pars <- list(expMean,pathA,pathC,pathE,varA,varC,varE,varP)

# specify the two groups:

modelMZ <- mxModel('MZ',pars,expCovMZ,dataMZ,objMZ,funcML,estVC)
modelDZ <- mxModel('DZ',pars,expCovDZ,dataDZ,objDZ,funcML,estVC)

### combine the data groups:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

### specify the confidence intervals:

ci <- mxCI(c('MZ.EstVC'))

### combine all of the model objects:

ACEModel <- mxModel('ACE',modelMZ,modelDZ,multi,ci)

### fit the model:

ACEFit <- mxTryHard(ACEModel,intervals=T)
summary(ACEFit)

### compare the fit of the model to the saturated model:

mxCompare(SatFit,ACEFit)


#=====
# Fit nested models
#=====

### here we will drop individual parameters (ie fix them to 0) to test
### their statistical significance

# AE model:

AEModel <- mxModel(ACEFit,name='AE')
AEModel <- omxSetParameters(AEModel,label='c_1_1',free=F,values=0)
AEFit <- mxTryHard(AEModel,intervals=T)
summary(AEFit)

# CE model:

CEModel <- mxModel(ACEFit,name='CE')
CEModel <- omxSetParameters(CEModel,label='a_1_1',free=F,values=0)
CEFit <- mxTryHard(CEModel,intervals=T)
summary(CEFit)

# E model:

EModel <- mxModel(ACEFit,name='E')
EModel <- omxSetParameters(EModel,label='a_1_1',free=F,values=0)
EModel <- omxSetParameters(EModel,label='c_1_1',free=F,values=0)
EFit <- mxTryHard(EModel,intervals=T)
summary(EFit)

### compare the nested models to the ACE model:

mxCompare(ACEFit,c(AEFit,CEFit,EFit))

save.image("UNIV_Motor.RData")

