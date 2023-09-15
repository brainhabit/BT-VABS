###########################################################################################################################################
# Script to run multivariate twin analysis of the BT-VABS paper investigating
# etiological structure underlying different domains of the Vineland instrument.
#
#
# G. Bussu, 28 November 2022 [adapted from a script tutorial provided by M. Taylor (Feb 2022)]
#
####################################
# The script needs as input the selected dataset in the wide format (twin_data.csv from the bt-vabs_prepdata script).
#
# The script gives as output results from the bivariate model with correlated factor solution compared to the saturated model.
#
###########################################################################################################################################

rm(list=ls())

setwd('C:/Users/giobu365/Documents/BT-VABS')
list.files()

### load OpenMx and functions:
require(OpenMx)
source('C:/Users/giobu365/Documents/twin_modelling/tutorial/miFunctions.R')

mxOption(key="Number of Threads",
         value=parallel::detectCores())

### prepare some data:
data <- read.csv(file='C:/Users/giobu365/Documents/BT-VABS/data/twin_data_VABS.csv',header=T,sep=',')
names(data);dim(data)

Vars <- c('resMotor','resSocComm')
nv <- 2 # number of phenotypes
ntv <- nv*2 # number of measures/variables
(selVars <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep=''))


mz <- subset(data,zygosity=='MZ',c(selVars))
dz <- subset(data,zygosity=='DZ',c(selVars))

############## Fully Saturated model ##################################################################################

svCov <- c(1,rep(.5,3),1,rep(.5,2),1,.5,1)


# Create Matrices for Covariates and linear Regression Coefficients
expMeanMZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                      labels=labFull('m_mz',1,ntv),name='ExpMeanMZ')
expMeanDZ <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                      labels=labFull('m_dz',1,ntv),name='ExpMeanDZ')


expCovMZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=svCov,
                     labels=labSymm('mz_cov',ntv),name='ExpCovMZ')

expCovDZ <- mxMatrix(type='Symm',nrow=ntv,ncol=ntv,free=T,values=svCov,
                     labels=labSymm('dz_cov',ntv),name='ExpCovDZ')

matI <- mxMatrix(type='Iden',nrow=ntv,ncol=ntv,name='I')

expCorMZ <- mxAlgebra(solve(sqrt(I*ExpCovMZ))%&%ExpCovMZ,name='ExpCorMZ')
expCorDZ <- mxAlgebra(solve(sqrt(I*ExpCovDZ))%&%ExpCovDZ,name='ExpCorDZ')

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:
objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMeanMZ',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMeanDZ',dimnames=selVars)

funcML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups

# specify the submodels for each zygosity group:
modelMZ <- mxModel('MZ',expMeanMZ,expCovMZ,matI,expCorMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',expMeanDZ,expCovDZ,matI,expCorDZ,dataDZ,objDZ,funcML)

# combine submodels into one object to make sure they are evaluated together:
multi <- mxFitFunctionMultigroup(c('MZ','DZ'))
ci <- mxCI(c('MZ.ExpCorMZ','DZ.ExpCorDZ'))

# combine all model objects:
SatModel <- mxModel('Sat',modelMZ,modelDZ,multi,ci)

mxOption(model= SatModel, key="Number of Threads", value= (omxDetectCores() - 1))

SatFit <- mxTryHard(SatModel,intervals=T)
sumSatFit<-summary(SatFit)

####################################### Constrained Model (estimate twin correlation) ############################################################

CorModel <- mxModel(SatFit,name='Cor')

# equal means across twins and zygosity:
CorModel <- omxSetParameters(CorModel,labels=c(labFull('m_mz',1,ntv),
                                               labFull('m_dz',1,ntv)),
                             free=T,values=.01,newlabels=labFull('m',1,nv))

# equal variances:
CorModel <- omxSetParameters(CorModel,labels=c(labDiag('mz_cov',ntv),
                                               labDiag('dz_cov',ntv)),
                             free=T,values=1.2,newlabels=labDiag('v',nv))

# equal phenotypic correlation in twin 1 and twin 2:
CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_2_1','mz_cov_4_3',
                                               'dz_cov_2_1','dz_cov_4_3'),
                             free=T,values=.3,newlabel='rph')

# fix the CTCT to be the same between twins:
CorModel <- omxSetParameters(CorModel,labels=c('mz_cov_4_1','mz_cov_3_2'),
                             free=T,values=.3,newlabel='ctct_mz')
CorModel <- omxSetParameters(CorModel,labels=c('dz_cov_4_1','dz_cov_3_2'),
                             free=T,values=.3,newlabel='ctct_dz')

# fit the model
CorFit <- mxTryHard(CorModel,intervals=T)
sumCorFit<-summary(CorFit)

# compare with fully saturated to check assumptions
mxCompare(SatFit,CorFit)

########################################## correlated factors solution ####################################################

### full model:

# path coefficients:
coefpath<-c(1,.5,1)

pathA <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('a',nv),name='a')
pathC <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('c',nv),name='c')  
pathE <- mxMatrix(type='Lower',nrow=nv,ncol=nv,free=T,values=coefpath,
                  labels=labLower('e',nv),name='e')  

# calculate the variance components:

varA <- mxAlgebra(a%*%t(a),name='A')
varC <- mxAlgebra(c%*%t(c),name='C')
varE <- mxAlgebra(e%*%t(e),name='E')

# calculate the total variance and the standard deviation:

varP <- mxAlgebra(A+C+E,name='V')

matI <- mxMatrix(type='Iden',nrow=nv,ncol=nv,name='I')
isd <- mxAlgebra(solve(sqrt(I*V)),name='iSD')

# calculate phenotypic and etiological correlations:

corP <- mxAlgebra(solve(sqrt(I*V))%&%V,name='rPH')
corA <- mxAlgebra(solve(sqrt(I*A))%&%A,name='rA')
corC <- mxAlgebra(solve(sqrt(I*C))%&%C,name='rC')
corE <- mxAlgebra(solve(sqrt(I*E))%&%E,name='rE')

# means:

expMean <- mxMatrix(type='Full',nrow=1,ncol=ntv,free=T,values=0.001,
                    labels=labFull('m',1,nv),name='ExpMean')


# variance/covariance:

expCovMZ <- mxAlgebra(rbind(cbind(V,A+C),
                            cbind(A+C,V)),name='ExpCovMZ')
expCovDZ <- mxAlgebra(rbind(cbind(V,0.5%x%A+C),
                            cbind(0.5%x%A+C,V)),name='ExpCovDZ')

# convert the variance components to proportions:

estVC <- mxAlgebra(cbind(A/V,C/V,E/V),name='EstVC')

# observed data:

dataMZ <- mxData(mz,type='raw')
dataDZ <- mxData(dz,type='raw')

# objectives:

objMZ <- mxExpectationNormal(covariance='ExpCovMZ',means='ExpMean',dimnames=selVars)
objDZ <- mxExpectationNormal(covariance='ExpCovDZ',means='ExpMean',dimnames=selVars)

# method of estimation:

funcML <- mxFitFunctionML()

# specify submodels for each zygosity group:

pars <- list(pathA,pathC,pathE,varA,varC,varE,varP,matI,isd,corP,corA,corC,corE,expMean)

modelMZ <- mxModel('MZ',pars,estVC,expCovMZ,dataMZ,objMZ,funcML)
modelDZ <- mxModel('DZ',pars,estVC,expCovDZ,dataDZ,objDZ,funcML)

# specify that we need to evaluate all submodels together:

multi <- mxFitFunctionMultigroup(c('MZ','DZ'))

# confidence intervals:

ci <- mxCI(c('MZ.EstVC','MZ.rPH','MZ.rA','MZ.rC','MZ.rE'))

# combine all objects:

CorFac<- mxModel('Cor',modelMZ,modelDZ,multi,ci)

# fit the model:

FacFit <- mxTryHard(CorFac,intervals=T)
sumCorFac<-summary(FacFit)

# compare fit to saturated model:

mxCompare(SatFit,FacFit)


  save.image("bivariate_results.RData")

##################################################### PLOTS ################################################################################################

library(ggplot2)

### Twin correlations
rMZ<-sumCorFit$CI[c(3,8),2]
rMZlow<-sumCorFit$CI[c(3,8),1]
rMZhigh<-sumCorFit$CI[c(3,8),3]

rDZ<-sumCorFit$CI[c(19,24),2]
rDZlow<-sumCorFit$CI[c(19,24),1]
rDZhigh<-sumCorFit$CI[c(19,24),3]

type<-c(rep('MZ',2),rep('DZ',2))
type<-ordered(type,levels=c('MZ','DZ'))
pheno<-rep(c('Motor','SocComm'),2)
pheno<-ordered(pheno,levels=c('Motor','SocComm'))

twinplot<-data.frame(cbind(c(rMZ,rDZ),c(rMZlow,rDZlow),c(rMZhigh,rDZhigh),pheno,type))
names(twinplot)<-c('Correlation','Low','High','Phenotype','Type')
twinplot$Correlation<-as.numeric(twinplot$Correlation)
twinplot$Low<-as.numeric(twinplot$Low)
twinplot$High<-as.numeric(twinplot$High)
twinplot$Type<-as.factor(twinplot$Type)
twinplot$Phenotype<-as.factor(twinplot$Phenotype)
levels(twinplot$Phenotype)<-c('Motor','SocComm')
levels(twinplot$Type)<-c('MZ','DZ')

png('twincorr.png',width=2000,height=2000,res=300)
ggplot(twinplot, aes(fill=Type, y=Correlation, x=Phenotype)) +
  geom_bar(position="dodge", stat="identity") + geom_errorbar(position=position_dodge(width=0.8),aes(ymin = Low, ymax = High), width =.2, col = "black")+
  scale_fill_brewer(palette='Set2')+
  theme_bw() +
  theme(legend.position="right")+
  xlab("VABS scale score") + ylab('Correlation')
dev.off()

### Cross-twin cross-trait correlations + phenotypic correlation
rph<-sumCorFit$CI[2,2]
rph_low<-sumCorFit$CI[2,1]
rph_high<-sumCorFit$CI[2,3]
ctct_mz<-sumCorFit$CI[4,2]
ctct_mz_low<-sumCorFit$CI[4,1]
ctct_mz_high<-sumCorFit$CI[4,3]
ctct_dz<-sumCorFit$CI[20,2]
ctct_dz_low<-sumCorFit$CI[20,1]
ctct_dz_high<-sumCorFit$CI[20,3]
pheno<-rep('Motor-SocComm',3)

corrtype<-c('Rph','CTCT-MZ','CTCT-DZ')
corrplot<-data.frame(cbind(c(rph,ctct_mz,ctct_dz),c(rph_low,ctct_mz_low,ctct_dz_low),c(rph_high,ctct_mz_high,ctct_dz_high),pheno,corrtype))
names(corrplot)<-c('Correlation','Lower bound','Upper bound','Phenotype','Type')
corrplot$Type<-ordered(corrplot$Type,levels=c('Rph','CTCT-MZ','CTCT-DZ'))
corrplot$Correlation<-as.numeric(corrplot$Correlation)
corrplot$`Lower bound`<-as.numeric(corrplot$`Lower bound`)
corrplot$`Upper bound`<-as.numeric(corrplot$`Upper bound`)

png('CTCT.png',width=2000,height=2000,res=300)
ggplot(corrplot, aes(fill=Type, y=Correlation, x=Phenotype)) +
  geom_bar(position="dodge", stat="identity") + geom_errorbar(position=position_dodge(width=0.8),aes(ymin = `Lower bound`, ymax = `Upper bound`), width =.2, col = "black")+
  scale_fill_brewer(palette='Set2')+
  theme_bw() +
  theme(legend.position="right")+
  xlab("VABS phenotype") + ylab('Correlation')
dev.off()


### Variance components
CorFacVar<-data.frame(matrix(round(sumCorFac$CI[,2],digits=2),nrow=4))
names(CorFacVar)<-c('A','C','E','rPh','rA','rC','rE')

CorFacVar_low<-data.frame(matrix(round(sumCorFac$CI[,1],digits=2),nrow=4))
names(CorFacVar_low)<-c('A','C','E','rPh','rA','rC','rE')

CorFacVar_high<-data.frame(matrix(round(sumCorFac$CI[,3],digits=2),nrow=4))
names(CorFacVar_high)<-c('A','C','E','rPh','rA','rC','rE')

vartype<-c(rep('A',2),rep('C',2),rep('E',2))
pheno<-rep(c('Motor','SocComm'),3)

plotdata<-data.frame(cbind(c(CorFacVar_low$A[c(1,4)],CorFacVar_low$C[c(1,4)],CorFacVar_low$E[c(1,4)]),c(CorFacVar$A[c(1,4)],CorFacVar$C[c(1,4)],CorFacVar$E[c(1,4)]),c(CorFacVar_high$A[c(1,4)],CorFacVar_high$C[c(1,4)],CorFacVar_high$E[c(1,4)]),pheno,vartype))
plotdata$pheno<-ordered(pheno,level=c('Motor','SocComm'))
plotdata$vartype<-as.factor(vartype)
names(plotdata)<-c('CI_low','Variance','CI_high','Phenotype','Components')
plotdata$CI_low<-as.numeric(plotdata$CI_low)
plotdata$CI_high<-as.numeric(plotdata$CI_high)
plotdata$Variance<-as.numeric(plotdata$Variance)

png('variance.png',width=2000,height=2000,res=300)
ggplot(plotdata, aes(fill=Components, y=Variance, x=Phenotype)) +
  geom_bar(position="stack", stat="identity") +
  #geom_errorbar(position=position_dodge(width=0.8),aes(ymin = CI_low, ymax = CI_high), width =.2, col = "black")+
  scale_fill_brewer(palette='Set2')+
  theme_bw() +
  theme(legend.position="right") +
  xlab("VABS score") + ylab('Variance (%)')
dev.off()


