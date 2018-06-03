#R code for simulation of multiple cases sibships datasets and analysis of the data


library(mvtnorm)
library(lme4)
library(numDeriv)


#Function to define the relationship matrix
#parameter is family size
kinship<-function(sizeFamily){
  kinship<-matrix(0.5,nrow=sizeFamily,ncol=sizeFamily)+0.5*diag(sizeFamily)
  kinship
}

#Function to define the matrix between primary and secondary phenotype with as parameters:
#size : family size
#varianceA and Variance B genetic variances for both phenotypes
#coeff: sigma u
covPhenotypes<-function(size,varianceA,varianceB,coeff){
  matrix<-varianceA*varianceB*kinship(size)
  matrix<-matrix+coeff^2*diag(1,nrow=size)
  matrix
}

#Function to define the variance covariance for each phenotype:
#size : family size
#genetic,error and u: are genetic variance, residual variance and sigma u
Variance<-function(size,genetic,error,u){
  matrix<-kinship(size)*genetic^2+(error^2+u^2)*diag(size)
  matrix
}

#Function to define the variance covariance for secondary phenotypes for adhoc methods:
#size : family size
#genetic,error: are genetic variance, residual variance
Variance2<-function(size,genetic,error){
  matrix<-kinship(size)*genetic^2+(error^2)*diag(size)
  matrix
}

#Creation of the table of conditional genotype distributions of offspring given parental genotypes
possibleParentalGenotype<-as.matrix(expand.grid(rep(list(c(0:2)),2)))
AA<-c(1,1/2,0,1/2,1/4,0,0,0,0)
aA<-c(0,1/2,1,1/2,1/2,1/2,1,1/2,0)
aa<-c(0,0,0,0,1/4,1/2,0,1/2,1)
IBD<-cbind(possibleParentalGenotype,AA,aA,aa)

#Function to compute fixed effects
meanGenotype<-function(x,Intercept,Slope){
  mean<-Intercept + Slope*x
  mean
}

#function to obtain lower boundary of the integrals
lower<-function(x){
  lower<-0
  if(x==0){
    lower<--Inf
  }
  lower
}

#function to obtain upper boundary of the integrals
upper<-function(x){
  upper<-+Inf
  if(x==0){
    upper<-0
  }
  upper
}

##################################################################################################
#Function to obtain a vector corresponding to the allele of each parents
#parameter MAF: Minor allele frequency assumed in the global population
ParentalGenotypeSimulationBis<-function(MAF){
  parent1<-c(rbinom(1,1,prob=MAF),rbinom(1,1,prob=MAF)) #obtention of the alleles for parent 1 
  parent2<-c(rbinom(1,1,prob=MAF),rbinom(1,1,prob=MAF)) #obtention of the alleles for parent 2
  parents<-c(parent1,parent2)
  parents
}

#Function to obtain a vector corresponding to the genotype of siblings
#parameter MAF: Minor allele frequency assumed in the global population
#size family
OffspringsGenotypesSimulationBis<-function(MAF,sizeFamily){
  familyGenotype<-NULL
  parents<-ParentalGenotypeSimulationBis(MAF)#simulation of parental genotypes
  for(i in 1:sizeFamily){
    father<-rbinom(1,1,prob=0.5) #random obtention of 1 allele of parent 1 (0 first allele 1 second allele)
    mother<-rbinom(1,1,prob=0.5) #random obtention of 1 allele of parent 2 (0 first allele 1 second allele)
    genotypeOff<-sum(c(parents[father+1],parents[mother+3])) #sum to obtain the value of the genotype for 1 sibling 	
    familyGenotype<-c(familyGenotype,genotypeOff)
  }
  familyGenotype
}

#Function to simulate trait values for members of the same family
#parameters MAF, size family, index family
#parametersX : vector containing values of betas and variance parameters (genetic and residual) for the secondary phenotype
#parametersY : vector containing values of betas and the genetic variance for the primary phenotype
#parameterXY : sigma u
simulationDataFamilyBis<-function(MAF,sizeFamily,parametersX,parametersY,parameterXY,indexFamily){
  genotype<-OffspringsGenotypesSimulationBis(MAF,sizeFamily) #obtention of the siblings genotype
  meanX<-meanGenotype(genotype,Intercept=parametersX[1],Slope=parametersX[2]) #computation of the mean of the secondary phenotype
  meanY<-meanGenotype(genotype,Intercept=parametersY[1],Slope=parametersY[2]) #computation of the mean of the primary phenotype
  varX<-Variance(sizeFamily,genetic=parametersX[3],error=parametersX[4],u=parameterXY) #obtention of the variance-covariance matrix for the secondary phenotype
  varY<-Variance(sizeFamily,genetic=parametersY[3],error=1,u=parameterXY) #obtention of the variance-covariance matrix for the primary phenotype
  covariance<-covPhenotypes(sizeFamily,varianceA=parametersX[3],varianceB=parametersY[3],coeff=parameterXY) #obtention of the covariance matrix between phenotypes
  varPhenotypes<-rbind(cbind(varY,covariance),cbind(covariance,varX)) #construction of the variance-covariance matrix of the joint distribution of the secondary phenotype and the latent variable
  Phenotypes<-rmvnorm(1,mean=c(meanY,meanX),varPhenotypes) #Phenotypes simulation
  X<-Phenotypes[(sizeFamily+1):(2*sizeFamily)] #Obtention of the value of the secondary trait
  pij<-Phenotypes[1:(sizeFamily)] #Obtention of the value of the latent variable
  probY<-pnorm(pij) #transformation of the latent variable to obtain the probability to get the disease or not
  Y<-rbinom(sizeFamily,size=rep(1,sizeFamily),prob=probY) #simulation of the primary phenotype values
  IF<-rep(indexFamily,sizeFamily)
  IID<-c(1:sizeFamily)
  dataFamily<-as.data.frame(cbind(IF,IID,genotype,X,Y)) #Creation of the dataset for the family
  dataFamily
}


#Function to simulate the full dataset
#parameters MAF, size family, index family, number family
#parametersX : vector containing values of betas and variance parameters (genetic and residual) for the secondary phenotype
#parametersY : vector containing values of betas and the genetic variance for the primary phenotype
#parameterXY : sigma u
#threshold : value used for ascertainment
#generator: seed
simulationDatasetBis<-function(MAF,sizeFamily,parametersX,parametersY,numberFamily,parameterXY,threshold,generator){
  set.seed(generator)
  families<-NULL
  dataset<-NULL
  j<-1 #intialisation of the family index
  while(length(families) < numberFamily){
    familyDataset<-simulationDataFamilyBis(MAF,sizeFamily,parametersX,parametersY,parameterXY,j) #obtention of the dataset for family j
    if (sum(familyDataset$Y)>= threshold){
      dataset<-as.data.frame(rbind(dataset,familyDataset)) #if family have have more or the same number of case than the threshold she is integrated in the dataset
      families<-c(families,j) #vector containing the index of each family integrated 
    }
    j<-j+1
  }
  dataset
}


#######################Genotype probability function#########################

#Function to get possible parental genotypes given siblings genotype
#parameters: the family dataset, and the table of conditional genotype distributions of offspring given parental genotypes
possibleParentalGenotypes<-function(dataset,IBD){
  possibleGenotypes<-c(1:9)#vector indexing the 9 possible parental genotype combination
  genotypes<-unique(dataset$genotype)#vector with the unique genotype value for the marker in the offsprings ((0) or (0,1) or (0,2) or (1,2) or (0,1,2))
  for(j in 1:length(genotypes)){
    possibleParentalgenotype<-which(IBD[,3+genotypes[j]]>0) #Getting the rows from the table IBD defined before in order to get the possible parental combination for the genotype value of a sibling
    possibleGenotypes<-intersect(possibleGenotypes,possibleParentalgenotype) #intersecting possibleParentalgenotype and possibleGenotypes in order to keep only the possible parental combination
  }
  return(possibleGenotypes)
}

#Function to compute the genotype probability of a family
#parameters : dataset of one family
#IBD and MAF
genotypeProbabilityFamily<-function(dataset,IBD,MAF){
  HWE<-c((1-MAF)^2,2*MAF*(1-MAF),MAF^2)#defining probabilities under HW Equilibrium
  result<-NULL
  genotypeColumn<-3+dataset$genotype
  possibleGenotypes<-possibleParentalGenotypes(dataset,IBD)#obtaining the possible parental genotype given offsprings genotypes
  probabilityGivenGenotypes<-NULL
  for(j in 1:length(possibleGenotypes)){#computing probability for each possible parental combination
    value<-prod(IBD[possibleGenotypes[j],genotypeColumn])*HWE[IBD[possibleGenotypes[j],1]+1]*HWE[IBD[possibleGenotypes[j],2]+1]#obtention of the probability of siblings genotype given  parental genotypes
    probabilityGivenGenotypes<-c(probabilityGivenGenotypes,value)
  }
  prob<-sum(probabilityGivenGenotypes) #sum of the different probabilities obtained
  prob
}


#######################################################################################################################
####################################################################################################################################################



#Function to compute the Probability of the joint distribution of phenotypes given the genotype
#Parameters: dataset simulated for 1 family + added variable obtained in the first step of the likelihood function and
#variance parameters
CatalanoLikelihood<-function(dataset,sizeFamily,VarX,VarY,covYX,meanX,meanY){
  sizeFamily<-length(dataset[,1]) #size of the family
  lower<-sapply(dataset$Y,lower) #obtention of all lower boundaries for the family members
  upper<-sapply(dataset$Y,upper) #obtention of all upper boundaries for the family members
  meanYgivenX<-as.numeric(meanY[(dataset$genotype+1)]+covYX%*%solve(VarX)%*%(dataset$X-meanX[(dataset$genotype+1)])) #compute the mean of Y|X
  #print(meanYgivenX)
  VarYgivenX<-VarY-covYX%*%solve(VarX)%*%t(covYX) #compute the varicance-covariance matrix of Y|X
  prob<-dmvnorm(dataset$X,mean=meanX[(dataset$genotype+1)],sigma=VarX)*pmvnorm(lower,upper,mean=meanYgivenX,sigma=VarYgivenX,algorithm=Miwa()) #Computation of the joint density for one family
  prob
}


#Function to compute the denominator of the retrospective likelihood
#Parameters, parametersY: betas of the primary phenotype, IBD, MAF
#VarY: variance covariance matrix of the primary phenotype
denominatorBis<-function(dataset,parametersY,IBD,MAF,VarY,meanY){
  sizeFamily<-length(dataset[,1])
  lower<-sapply(dataset$Y,lower) #obtention of all lower boundaries for the family members
  upper<-sapply(dataset$Y,upper) #obtention of all upper boundaries for the family members	
  possibleOffspringsGenotype<-as.matrix(expand.grid(rep(list(c(0:2)),sizeFamily))) #obtention of all the possible genotype combination for offsprings
  value<-NULL
  for(j in 1:length(possibleOffspringsGenotype[,1])){ #for each possible genotype computation of the genotype probability
    dataset$genotype<-possibleOffspringsGenotype[j,] #replacing in the dataset the genotype value
    mean.Y<-meanY[(dataset$genotype+1)] #obtention of the mean of the primary phenotype with the genotype	
    ProbY<-pmvnorm(lower,upper,mean=mean.Y,sigma=VarY,algorithm=Miwa()) #computing the marginal probability of Y given the genotypes
    GenotypeProbability<-genotypeProbabilityFamily(dataset,IBD,MAF) #computing the genotype probability 
    value<-c(value,ProbY*GenotypeProbability)
  }
  denominatorValue<-sum(value) #sum of all the genotype probability to obtain the denominator value
  denominatorValue
}

likelihoodFamily<-function(x,dataset,VarY,VarX,CovYX,IBD,MAF,meanX,meanY,pattern,sizeFamily){
  family<-dataset[dataset[,1]==x,] #obtention of a subset of the dataset (one family)
  jointDensity<-CatalanoLikelihood(family,sizeFamily,VarX,VarY,CovYX,meanX,meanY) #Computation of the joint probability
  genotype<-genotypeProbabilityFamily(family,IBD,MAF) #Computation of the genotype probability
  numerator<-jointDensity*genotype #computation of the numerator
  denom<-pattern[(sum(family$Y)-threshold+1)] #obtention of the denominator value corresponding of the number of case in the family	
  resultsFamily<-numerator/denom #computation of the family value
  return(resultsFamily)
}

Likelihood.newer.algorithm<-function(parameters,dataset,IBD,MAF,sizeFamily,threshold){
    #print(parameters)
    family.id<-unique(dataset[,1])
    number.family<-length(family.id)
    parametersX<-parameters[1:4]
    parametersX[3:4]<-exp(parametersX[3:4])
    parametersY<-parameters[5:7]
    parametersY[3]<-exp(parametersY[3])
    parameterXY<-exp(parameters[8])
    possibleGenotype<-seq(0,2,by=1)
    meanX<-sapply(possibleGenotype,meanGenotype,Intercept=parametersX[1],Slope=parametersX[2])
    meanY<-sapply(possibleGenotype,meanGenotype,Intercept=parametersY[1],Slope=parametersY[2])
    VarX<-Variance(sizeFamily,genetic=parametersX[3],error=parametersX[4],u=parameterXY) #obtention of the variance-covariance matrix for the secondary phenotype
    VarY<-Variance(sizeFamily,genetic=parametersY[3],error=1,u=parameterXY) #obtention of the variance-covariance matrix for the primary phenotype
    covYX<-covPhenotypes(sizeFamily,varianceA=parametersX[3],varianceB=parametersY[3],coeff=parameterXY) #obtention of the covariance matrix between phenotypes
    pattern<-NULL
    for (k in threshold:sizeFamily){ #to increase speed computation, computation of the denimonator for each possible case
      Y<-rep(0,sizeFamily)
      Y[1:k]<-rep(1,k)
      t<-as.data.frame(cbind(fid=rep(1,sizeFamily),Y=Y,X=rep(1,sizeFamily),genotype=rep(1,sizeFamily))) #creation of a dataset for using function denominator only Y values are important
      denomina<-denominatorBis(t,parametersY,IBD,MAF,VarY,meanY) #computation of the denominator value
      pattern<-c(pattern,denomina)
    }
    res<-sapply(family.id,likelihoodFamily,dataset=dataset,VarY=VarY,VarX=VarX,CovYX=covYX,IBD=IBD,MAF=MAF,meanX=meanX,meanY=meanY,pattern=pattern,sizeFamily=sizeFamily)
    result<-sum(log(res))
}

Likelihood.newer.algorithm.NULL<-function(parameters,dataset,IBD,MAF,sizeFamily,threshold){
  #print(parameters)
  family.id<-unique(dataset[,1])
  number.family<-length(family.id)
  parametersX<-parameters[1:3]
  parametersX[2:3]<-exp(parametersX[2:3])
  parametersY<-parameters[4:6]
  parametersY[6]<-exp(parameters[6])
  parameterXY<-exp(parameters[7])
  possibleGenotype<-seq(0,2,by=1)
  meanX<-sapply(possibleGenotype,meanGenotype,Intercept=parametersX[1],Slope=0)
  #print(meanX)
  meanY<-sapply(possibleGenotype,meanGenotype,Intercept=parametersY[1],Slope=parametersY[2])
  VarX<-Variance(sizeFamily,genetic=parametersX[2],error=parametersX[3],u=parameterXY) #obtention of the variance-covariance matrix for the secondary phenotype
  VarY<-Variance(sizeFamily,genetic=parametersY[3],error=1,u=parameterXY) #obtention of the variance-covariance matrix for the primary phenotype
  covYX<-covPhenotypes(sizeFamily,varianceA=parametersX[2],varianceB=parametersY[3],coeff=parameterXY) #obtention of the covariance matrix between phenotypes
  pattern<-NULL
  for (k in threshold:sizeFamily){ #to increase speed computation, computation of the denimonator for each possible case
    Y<-rep(0,sizeFamily)
    Y[1:k]<-rep(1,k)
    t<-as.data.frame(cbind(fid=rep(1,sizeFamily),Y=Y,X=rep(1,sizeFamily),genotype=rep(1,sizeFamily))) #creation of a dataset for using function denominator only Y values are important
    denomina<-denominatorBis(t,parametersY,IBD,MAF,VarY,meanY) #computation of the denominator value
    pattern<-c(pattern,denomina)
  }
  res<-sapply(family.id,likelihoodFamily,dataset=dataset,VarY=VarY,VarX=VarX,CovYX=covYX,IBD=IBD,MAF=MAF,meanX=meanX,meanY=meanY,pattern=pattern,sizeFamily=sizeFamily)
  result<-sum(log(res))
}

#parameters for dataset simulation
Genetic.sd.x<-2
Genetic.sd.y<-sqrt(3)
beta.Y<-0.1
beta.X<-0
ascertainment<-2
sd.u<-sqrt(2)
ascertainment<-2
preval<-(-2.326348)
MAF<-0.3
sizeFamily<-5
numberFamily<-400
parametersX<-c(3.5,beta.X,Genetic.sd.x,sqrt(2))
parametersY<-c(preval[n],beta.Y,Genetic.sd.y)
numberFamily<-100
threshold<-ascertainment[n]



#dataset simulation
dataset<-simulationDatasetBis(MAF,sizeFamily,parametersX,parametersY,numberFamily,sd.u,threshold,m+3742) #creation of the dataset
write.table(dataset,"SimulatedDataset.txt",row.names =FALSE, sep="\t",quote=F)

#Computation of intial values
e<-glmer(Y~(1|IF)+genotype,data=dataset,family=binomial(link="probit")) #running glmer model in order to get primary phenotype initial values for optim
primaryFixed<-as.numeric(fixef(e))  #obtention of the intial values of the fixed effects of the primary phenotype
vcPrimary<-VarCorr(e)[[1]] #obtention of the initial value of the genetice parameter for the primary phenotype
f<-lmer(X~(1|IF)+genotype,data=dataset) #running lmer model in order to get secondary phenotype initial values for optim
secondaryFixed<-as.numeric(fixef(f)) #obtention of the intial values of the fixed effects of the secondary phenotype
vcSecondary <-as.numeric(VarCorr(f)[[1]]) #obtention of the initial value of the genetic parameter for the secondary phenotype
residuals<-attr(VarCorr(f),"sc") #obtention of the initial value of the residual parameter for the secondary phenotype

#optimization of the likelihood
res1<-try(optim(c(secondaryFixed[1],secondaryFixed[2],log(vcSecondary),log(residuals),primaryFixed[1],primaryFixed[2],log(sqrt(3)),0),Likelihood.newer.algorithm,dataset=dataset,IBD=IBD,MAF=MAF,sizeFamily=sizeFamily,threshold=threshold,method="BFGS",control=list(trace=1,fnscale=-1,parscale=c(1,0.1,0.1,0.1,1,0.1,0.1,1),maxit=200),hessian=TRUE)) #maximisation of the likelihood for the Full method
try(StdError<-sqrt(diag(solve(-res1$hessian)))) #obtention of the standard errors
#results are wriiten in a txt file. All variance paramters are log transforme. Use exponential function for getting the true variance parameters value.
try(write.table(t(c(m,res1$par[1],StdError[1],res1$par[2],StdError[2],res1$par[3],StdError[3],res1$par[4],StdError[4],res1$par[5],StdError[5],res1$par[6],StdError[6],res1$par[7],StdError[7],res1$par[8],StdError[8],res1$valueq)),file=paste("Results.txt",sep=""),append=TRUE,quote=FALSE,sep="\t",row.names = FALSE,col.names = FALSE)) #saving results obtained
