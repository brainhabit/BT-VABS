###########################################################################################################################################
# Script to run post-hoc power analysis based on Monte-Carlo simulations 
#
# Project: BT-VABS - investigating etiological structure underlying different domains of the Vineland instrument.
#
# Author: Giorgia Bussu
#
# Date: 12 February 2025
#
####################################
## Preparations
# Clean up
rm(list=ls())

# Load estimated twin model for post-hoc analysis
load("~/Downloads/btvabs/bivariate_results.RData")

# Load necessary libraries
library(OpenMx)
library(ggplot2)
library(MASS)

# Number of simulations
nSim <- 1046
sigLevel <- 0.05

# Storage for significant results
sigA1 <- numeric(nSim)
sigC1 <- numeric(nSim)
sigE1 <- numeric(nSim)
sigA2 <- numeric(nSim)
sigC2 <- numeric(nSim)
sigE2 <- numeric(nSim)
sig_corA <- numeric(nSim)
sig_corC <- numeric(nSim)
sig_corE <- numeric(nSim)

# n pairs and zygosity
nmz<-nrow(data[which(data$zygosity=='MZ'),])
ndz<-nrow(data[which(data$zygosity=='DZ'),])
npairs<-nmz+ndz

# Add variability to model estimates
varFactor <- 1

## Start simulation loop
for (i in 1:nSim) {
  
  set.seed(NULL)
  
  # Estimated parameters with variations
  estA1 <- max(min(0.21*varFactor+rnorm(1,0,0.05),1),0)
  estA2 <- max(min(0.12*varFactor+rnorm(1,0,0.05),1),0)
  estRA <- max(min(-0.39*varFactor+rnorm(1,0,0.05),1),-1)
  estC1 <- max(min(0.67*varFactor-rnorm(1,0,0.05),1),0)
  estC2 <- max(min(0.78*varFactor-rnorm(1,0,0.05),1),0)
  estRC <- max(min(0.45*varFactor+rnorm(1,0,0.05),1),-1)
  estE1 <- max(min(1-estA1-estC1,1),0)
  estE2 <- max(min(1-estA2-estC2,1),0)
  estRE <- max(min(0.32*varFactor+rnorm(1,0,0.05),1),-1)
  
  ## Simulate effects used to generate simulated datasets
  simA1 <- rnorm(npairs,mean=0,sd=sqrt(estA1))
  simC1 <- rnorm(npairs,mean=0,sd=sqrt(estC1))
  simE1 <- rnorm(2*npairs,mean=0,sd=sqrt(estE1))
  simA2 <- estRA*simA1 + rnorm(npairs, mean=0,sd=sqrt(estA2*(1-estRA*estRA)))
  simC2 <- estRC*simC1 + rnorm(npairs, mean=0,sd=sqrt(estC2*(1-estRC*estRC)))
  simE2 <- estRE*simE1 + rnorm(2*npairs,mean=0,sd=sqrt(estE2*(1-estRE*estRE)))
  
  # MZ-DZ twins: trait 1
  MZ_t1 <- simA1[1:nmz]+simC1[1:nmz]+simE1[1:nmz]
  MZ_t2 <- simA1[1:nmz]+simC1[1:nmz]+simE1[(nmz+1):(2*nmz)]
  DZ_t1 <- simA1[(nmz+1):(nmz+ndz)]+simC1[(nmz+1):(nmz+ndz)]+simE1[(2*nmz+1):(2*nmz+ndz)]
  DZ_t2 <- 0.5*simA1[(nmz+1):(nmz+ndz)]+simC1[(nmz+1):(nmz+ndz)]+simE1[(2*nmz+ndz+1):(2*nmz+2*ndz)]
  
  # trait 2
  MZ_t1_2 <- simA2[1:nmz]+simC2[1:nmz]+simE2[1:nmz]
  MZ_t2_2 <- simA2[1:nmz]+simC2[1:nmz]+simE2[(nmz+1):(2*nmz)]
  DZ_t1_2 <- simA2[(nmz+1):(nmz+ndz)]+simC2[(nmz+1):(nmz+ndz)]+simE2[(2*nmz+1):(2*nmz+ndz)]
  DZ_t2_2 <- 0.5*simA2[(nmz+1):(nmz+ndz)]+simC2[(nmz+1):(nmz+ndz)]+simE2[(2*nmz+ndz+1):(2*nmz+2*ndz)]
  
  ## simulated data
  simData <- data.frame(zygosity=c(rep('MZ',nmz),rep('DZ',ndz)),
  resMotor1=c(MZ_t1,DZ_t1),
  resMotor2=c(MZ_t2,DZ_t2),
  resSocComm1=c(MZ_t1_2,DZ_t1_2),
  resSocComm2=c(MZ_t2_2,DZ_t2_2))
  
  ## ru model
  mxOption(key="Number of Threads",
           value=parallel::detectCores())
  
  Vars <- c('resMotor','resSocComm')
  nv <- 2 # number of phenotypes
  ntv <- nv*2 # number of measures/variables
  (selVars <- paste(Vars,c(rep(1,nv),rep(2,nv)),sep=''))
  
  
  mz <- subset(simData,zygosity=='MZ',c(selVars))
  dz <- subset(simData,zygosity=='DZ',c(selVars))
  
  
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
  
  simFit <- mxTryHard(CorFac,intervals=T)
  
  # Test significance of the bivariate estimates
  if (simFit$output$status$code==0) {
    
    ## Test Power for Univariate Parameters
    # Create reduced models for A1, C1, E1 and A2, C2, E2
    
    # CE model for Variable 1 (fix A1 = 0)
    CE_model1 <- mxModel(simFit)
    CE_model1 <- omxSetParameters(CE_model1, labels = 'a_1_1', values = 0, free = FALSE)
    fit_CE1 <- try(mxTryHard(CE_model1), silent = TRUE)
    
    # AE model for Variable 1 (fix C1 = 0)
    AE_model1 <- mxModel(simFit)
    AE_model1 <- omxSetParameters(AE_model1, labels = 'c_1_1', values = 0, free = FALSE)
    fit_AE1 <- try(mxTryHard(AE_model1), silent = TRUE)
    
    # AC model for Variable 1 (fix E1 = 0)
    AC_model1 <- mxModel(simFit)
    AC_model1 <- omxSetParameters(AC_model1, labels = 'e_1_1', values = 0, free = FALSE)
    fit_AC1 <- try(mxTryHard(AC_model1), silent = TRUE)
    
    ## Repeat for Variable 2
    # CE model for Variable 2 (fix A2 = 0)
    CE_model2 <- mxModel(simFit)
    CE_model2 <- omxSetParameters(CE_model2, labels = 'a_2_2', values = 0, free = FALSE)
    fit_CE2 <- try(mxTryHard(CE_model2), silent = TRUE)
    
    # AE model for Variable 2 (fix C2 = 0)
    AE_model2 <- mxModel(simFit)
    AE_model2 <- omxSetParameters(AE_model2, labels = 'c_2_2', values = 0, free = FALSE)
    fit_AE2 <- try(mxTryHard(AE_model2), silent = TRUE)
    
    # AC model for Variable 2 (fix E2 = 0)
    AC_model2 <- mxModel(simFit)
    AC_model2 <- omxSetParameters(AC_model2, labels = 'e_2_2', values = 0, free = FALSE)
    fit_AC2 <- try(mxTryHard(AC_model2), silent = TRUE)
    
    ## Likelihood Ratio Tests for A1, C1, E1, A2, C2, E2
    if (fit_CE1$output$status$code==0) {
      sigA1[i] <- ifelse(mxCompare(simFit, fit_CE1)$p[2] < sigLevel, 1, 0)
    }else(sigA1[i]=1)
    if (fit_AE1$output$status$code==0) {
      sigC1[i] <- ifelse(mxCompare(simFit, fit_AE1)$p[2] < sigLevel, 1, 0)
    }else(sigC1[i]=1)
    if (fit_AC1$output$status$code==0) {
      sigE1[i] <- ifelse(mxCompare(simFit, fit_AC1)$p[2] < sigLevel, 1, 0)
    }else(sigE1[i]=1)
    if (fit_CE2$output$status$code==0) {
      sigA2[i] <- ifelse(mxCompare(simFit, fit_CE2)$p[2] < sigLevel, 1, 0)
    }else(sigA2[i]=1)
    if (fit_AE2$output$status$code==0) {
      sigC2[i] <- ifelse(mxCompare(simFit, fit_AE2)$p[2] < sigLevel, 1, 0)
    }else(sigC2[i]=1)
    if (fit_AC2$output$status$code==0) {
      sigE2[i] <- ifelse(mxCompare(simFit, fit_AC2)$p[2] < sigLevel, 1, 0)
    }else(sigE2[i]=1)
    
    ## Extract and Test Power for Correlations (rA, rC, rE)
    
    # fix rA = 0
    rA_model2 <- mxModel(simFit)
    rA_model2 <- omxSetParameters(rA_model2, labels = 'a_2_1', values = 0, free = FALSE)
    fit_rA2 <- try(mxRun(rA_model2), silent = TRUE)
    
    if (fit_rA2$output$status$code==0) {
      sig_corA[i] <- ifelse(mxCompare(simFit, fit_rA2)$p[2] < sigLevel, 1, 0)
    }else(sig_corA[i]=1)
    
    # fix rC = 0
    rC_model2 <- mxModel(simFit)
    rC_model2 <- omxSetParameters(rC_model2, labels = 'c_2_1', values = 0, free = FALSE)
    fit_rC2 <- try(mxRun(rC_model2), silent = TRUE)
    
    if (fit_rC2$output$status$code==0) {
      sig_corC[i] <- ifelse(mxCompare(simFit, fit_rC2)$p[2] < sigLevel, 1, 0)
    }else(sig_corC[i]=1)
    
    # fix rE = 0
    rE_model2 <- mxModel(simFit)
    rE_model2 <- omxSetParameters(rE_model2, labels = 'e_2_1', values = 0, free = FALSE)
    fit_rE2 <- try(mxRun(rE_model2), silent = TRUE)
    
    if (fit_rE2$output$status$code==0) {
      sig_corE[i] <- ifelse(mxCompare(simFit, fit_rE2)$p[2] < sigLevel, 1, 0)
    }else(sig_corE[i]=1)
  }
}

# Compute post hoc power for the univariate parameters
power_A1 <- mean(sigA1, na.rm = TRUE)
power_C1 <- mean(sigC1, na.rm = TRUE)
power_E1 <- mean(sigE1, na.rm = TRUE)
power_A2 <- mean(sigA2, na.rm = TRUE)
power_C2 <- mean(sigC2, na.rm = TRUE)
power_E2 <- mean(sigE2, na.rm = TRUE)

# Compute post hoc power for the correlations
power_corA <- mean(sig_corA, na.rm = TRUE)
power_corC <- mean(sig_corC, na.rm = TRUE)
power_corE <- mean(sig_corE, na.rm = TRUE)

# Print results
cat("Monte Carlo Post Hoc Power for A1 (Variable 1):", power_A1, "\n")
cat("Monte Carlo Post Hoc Power for C1 (Variable 1):", power_C1, "\n")
cat("Monte Carlo Post Hoc Power for E1 (Variable 1):", power_E1, "\n")
cat("Monte Carlo Post Hoc Power for A2 (Variable 2):", power_A2, "\n")
cat("Monte Carlo Post Hoc Power for C2 (Variable 2):", power_C2, "\n")
cat("Monte Carlo Post Hoc Power for E2 (Variable 2):", power_E2, "\n")
cat("Monte Carlo Post Hoc Power for corA (Genetic Correlation):", power_corA, "\n")
cat("Monte Carlo Post Hoc Power for corC (Shared Environmental Correlation):", power_corC, "\n")
cat("Monte Carlo Post Hoc Power for corE (Unique Environmental Correlation):", power_corE, "\n")

# Combine results for visualization
power_results <- data.frame(
  Parameter = ordered(as.factor(c('A1', 'C1', 'E1', 'A2', 'C2', 'E2', 'corA', 'corC', 'corE')),levels=c('A1', 'C1', 'E1', 'A2', 'C2', 'E2', 'corA', 'corC', 'corE')),
  Power = c(power_A1, power_C1, power_E1, power_A2, power_C2, power_E2, power_corA, power_corC, power_corE)
)

ggplot(power_results, aes(x = Parameter, y = Power, fill = Parameter)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  labs(title = "Post Hoc Power Analysis for Bivariate Model",
       x = "Parameter", y = "Power") +
  scale_fill_manual(values = c("steelblue", "forestgreen", "salmon", "gold", "purple", "skyblue", "pink", "orange", "yellow")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
