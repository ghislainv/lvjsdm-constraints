### Code pour travailler avec: fichier: Hmsc_structure_Simule_10thin_noreordering_3L_V2.RData ici: E:\Dossier Frederic\Articles\TPE\DATA_ANALYSIS\RData

path_alllvmodel<-"E:\\Dossier Frederic\\Articles\\TPE\\DATA_ANALYSIS\\RData\\Hmsc_structure_Simule_10thin_noreordering_3L_V2.RData"
path_nolvmodel<-"E:\\Dossier Frederic\\Articles\\TPE\\DATA_ANALYSIS\\RData\\Simule_10thin_noreordering_0L_V2.RData"

### seed for the pseudo-random number generator
seed<-5

library(Hmsc)
library(coda)
library(nimble)
mod_Hmsc_simulate_coda<-convertToCodaObject(mod_Hmsc_simulate)
names(mod_Hmsc_simulate_coda)


### code used in Ovaskainen & Abrego to to convergence check, p.48
######gelman.diag(mod_Hmsc_simulate_coda)
### does not work ???!!!
### OK corrigé/précisé p.75:
gelman.diag(mod_Hmsc_simulate_coda$Beta,multivariate=FALSE)$psrf
### non convergence des Omega:
 gelman.diag(mod_Hmsc_simulate_coda$Omega[[1]],multivariate=FALSE)$psrf



### code to get the target nlv (number of latent variables):
load(path_alllvmodel)
out.diag<-convertToCodaObject(mod_Hmsc_simulate)$Eta[[1]]
#starting from 1 nlv
nlv<-1
nlv.OK<-TRUE
while (nlv.OK)
{nlv<-nlv+1
name.lv<-paste0("factor",nlv)
nlv.OK<-sum(regexpr(name.lv,dimnames(out.diag[[1]])[[2]])>0)>0
}
nlv.ref<-nlv-1

### code for Remedy29 for Hmsc: to be applied to 0lv Hmsc model

### loading the output of the 0-lv JSDM 
### 
load(path_nolvmodel)

  ## putting the MCMC model output in an object called out.diag: here only Beta component because (i) 0lv; (ii) beta0 is in beta; (iii) no alphapar parameter
 out.diag<-convertToCodaObject(mod_Hmsc_simulate)$Beta
 ### to be adapted to 0lv model when available
 ###out.diag<-convertToCodaObject(get("mod_Hmsc_simulate",pos=2))$Beta

### converting it to a trully mcmc.list format

  # requiring coda package
  
  print("")
  print("################################################################")
  print("Diagnostics of species order:")

  ## putting the MCMC data in an object called ModelDataTMB
  #ModelDataTMB<-get("ModelConsts",pos=2)
  ## attaching this object in position 3
  #try(attach(ModelDataTMB,pos=3))
  #attached_3<-TRUE
  ### for Hmsc, we sioll recerate data for the code from the infos in mod_Hmsc_simulate

  ## if object idx_s that cpontains identifiers of the lower trianguler matrix of the nspecies*nespecies matrix does not exist, create it
  nspecies<-mod_Hmsc_simulate$ns
  ### formula valid in our unhierrachical study designs
  nsites<-mod_Hmsc_simulate$ny
  nvar<-mod_Hmsc_simulate$nc
  ### lower_triangular species*species matrix
  if (!exists("idx_s"))
  {
    temp<-matrix(0,ncol=nspecies,nrow=nspecies)
    idx_s <- lower.tri(temp)
  }


### method of species order correction n°29: Reordering species based on truemfcdiuncPCA of repeated quantile residuals of random parameter (without replacement) of model with 0lv
				### without using latent variable dnorm(CL,1)
      try({
      ## transform the mcmc.list out.diag into a matrix
      out.diag.matrix<-as.matrix(out.diag)
      
      
      #specification of the number of MCMCoutcomes taken into account to reorder species:
      Nrep<-1000
      Nrep<-min(Nrep,dim(out.diag.matrix)[1])
      
      #random selection of the 1000 values that will be used
	  ### ??? check whether with or without replacvement
      set.seed(seed)
      sample.ref<-sample(1:(dim(out.diag.matrix)[1]),Nrep,replace=FALSE)
      
      # transferring the 0-lv jsdm parameters corresponding to these outcomes in object samples.to.treat.EC
      samples.to.treat.EC<-out.diag.matrix[sample.ref,]
      
      
      
      
       nlv<-dim(mod_Hmsc_simulate_coda$Eta[[1]][[1]])[2]/nsites
	  #nlv<-0
                i<-1
                ### Nimble Code for the 0 lv Hmsc model (except for priors):
			ModelData<-list(Y=as.vector(mod_Hmsc_simulate$Y))
			ModelConsts<-list(nvar=nvar,nspecies=nspecies,nobs=nspecies*nsites,
					sd10=10,
					beta0=rep(1:nspecies,each=nsites),
					sites=rep(1:nsites,times=nspecies),
					env=cbind(rep(1,nsites),mod_Hmsc_simulate$XData),
					### specify the number of latent variables on which we want to reorder species##not active in the nimbleCode below
					nlv=nlv.ref)
			
			names.samples.EC<-dimnames(samples.to.treat.EC)[[2]]
				  ## changing names so that treatable by Nimble Model:
				  
				for (i in 1:nspecies)
					{ for (j in 1:nvar) {
						names.samples.EC[j+(i-1)*nvar]<-paste0("Beta[",j,",",i,"]")
						}
						
						}

logLiks<-matrix(0,ncol=Nrep,nrow=ModelConsts$nobs)
#indices.ref<-sapply(1:ModelConsts$nobs,function(x){sapply(1:ModelConsts$nvar,function(y){which(names.samples.EC==paste0("Beta[",y,",",ModelConsts$beta0[x],"]"))})})
#temp<-sapply(1:Nrep,function(x){sapply(1:ModelConsts$nobs,function(y){sum(samples.to.treat.EC[x,indices.ref[,y]]*ModelConsts$env[ModelConsts$sites[y],])})})

for (j in 1:Nrep)
{	print(j)
	for (i in 1:ModelConsts$nobs)
		{ for (k in 1:ModelConsts$nvar)
			{
			logLiks[i,j]<-logLiks[i,j]+samples.to.treat.EC[j,names.samples.EC==paste0("Beta[",k,",",ModelConsts$beta0[i],"]")]*ModelConsts$env[ModelConsts$sites[i],k]
			}
		}
}
	  		
      				  
      		
      			
      ## calculation of the best reordering for each MCMC output: one calculation for each set of MCMC outcomes
      
	  set.seed(seed+1)
      newspecies<-sapply(1:Nrep,function(x) {
            #selected.line<-sample.ref[x]
            
            ## for quantile residuals we need to have a vector of random uniform samples the same size as ModelData$Y
            randomizer<-runif(length(ModelData$Y))
            
            ### calculation of normalized quantile residuals as if we had a GLM with Bernoulli distribution and probit link function
            qresiduals<-qnorm(pbinom(ModelData$Y-1,1,pnorm(logLiks[,x]))+randomizer*dbinom(ModelData$Y,1,pnorm(logLiks[,x])))
            
            ### putting these residuals in a site*species matrix format
            ##qresidSpeciesMatrix<-sapply(tapply(1:length(ModelData$Y),Species,function(x){qresiduals[x]}),I,simplify=TRUE)
			qresidSpeciesMatrix<-sapply(tapply(1:length(ModelData$Y),ModelConsts$beta0,function(x){qresiduals[x]}),I,simplify=TRUE)
            
            #performing PCA of this residual matrix, without centering & scaling
            testPCA<-prcomp(qresidSpeciesMatrix,center=FALSE,scale=FALSE)
            
            #calculating newspecies, the index of the species - on the lambda scale???- to put on the diagonal
            newspecies<-NULL ### will be the number of the species -on the lambdascale, that is corresponding to species_order-  to put on the diagonal 
            
            for (i in 1:ModelConsts$nlv)
            {
            	to_substract<-rep(0,length(testPCA$sdev))
            	### ??? FG: explian why; and why this is from above: j>i???
            	 if (i<=(ModelConsts$nlv-1)) {
                	for (j in (i+1):ModelConsts$nlv)
                	{to_substract<-to_substract+(testPCA$rotation[,j]^2*testPCA$sdev[j]^2)}
                }
            	
            	### 22 FG: explain why the metric is testPCA$rotation[,i]^2*testPCA$sdev[i]^2
            	## chosing species to constrain on axis i by chosing the species with maximum metric: testPCA$rotation[,i]^2*testPCA$sdev[i]^2)-to_substract
            	testsp<-which.max((testPCA$rotation[,i]^2*testPCA$sdev[i]^2)-to_substract)[1]
            	
            	newspecies<-c(newspecies,testsp)
            	
            	### putting rotation component for the selected species to 0 so that thius species is not selected for subsequent axes
            	testPCA$rotation[testsp,]<-0
            	
            }
            
            
            newspecies})
      
      #printing the output, which should be a matrix of nlv rows and Nrep columns: each column contains the optimal species order for the first nlv axes for the MCMC outcome corresponding to the column.
      print(newspecies)
      
      ### forming the vector of Nrep characters that contain in charcater format the optimal reordering of species for the corresponding MCMC outcome; we then transform it into a factor vector and then summarize it, putting in first place the combination of species which is the most frequent; then taking the names of the first component of the summary
      bestcomb<- names(summary(as.factor(apply(newspecies,2,paste, collapse=",")))[1])
      
      ### then decomposing the best combination (as character) to the series of nlv species (in numbers); we put these in the vector: newspecies
      cuts<-gregexpr(",",bestcomb)[[1]]
      newspecies<-NULL
      for (i in 1:(ModelConsts$nlv-1))
      {
        temp<-as.double(substring(bestcomb,first=max(cuts[i-1]+1,1),last=cuts[i]-1))
        newspecies<-c(newspecies,temp)
      }
      temp<- as.double(substring(bestcomb,first=max(cuts[ModelConsts$nlv-1]+1,1),last=nchar(bestcomb)))
      newspecies<-c(newspecies,temp)
      
     
      print(newspecies)
      
	  
	  newspecies_29<-newspecies
      
     
      ###??? to be adapted to hmsc foprmat??? to put species in the new order for the new hmsc
     ###???
      
      
      })





### method of species order correction n°30: Reordering species based on mfcdiwoqruncPCA of repeated quantile residuals of random parameter (without replacement) of model with 0lv
		### with using latent variable dnorm(CL,1)

########################## ACTUALLY IMPOSSIBBLE WITH hmsc because latent variable just before Y loglik calculation not saved in hmsc


### method of species order correction n°3: Reordering species based on Rhat on lambdas, factor by factor

try({

load("E:\\Dossier Frederic\\Articles\\TPE\\DATA_ANALYSIS\\RData\\Hmsc_structure_Simule_10thin_noreordering_3L_V2.RData")
  ## putting the MCMC model output in an object called out.diag: here only Beta component because (i) 0lv; (ii) beta0 is in beta; (iii) no alphapar parameter
 out.diag<-convertToCodaObject(mod_Hmsc_simulate)$Lambda[[1]]
 
#preparing matrix with indices of the lower trinagular part: tempb
tempb<-NA*idx_s
tempb[idx_s]<-1:sum(idx_s)

## transform the mcmc.list out.diag into a matrix
      out.diag.matrix<-as.matrix(out.diag)
      
      
      #specification of the number of MCMCoutcomes taken into account to reorder species:
      Nrep<-1000
      Nrep<-min(Nrep,dim(out.diag.matrix)[1])
      
      #random selection of the 1000 values that will be used
	  ### ??? check whether with or without replacvement
      set.seed(seed)
      sample.ref<-sample(1:(dim(out.diag.matrix)[1]),Nrep,replace=FALSE)
      
      # transferring the 0-lv jsdm parameters corresponding to these outcomes in object samples.to.treat.EC
      samples.to.treat.EC<-out.diag.matrix[sample.ref,]
      
      
      
      
       
      nlv<-dim(mod_Hmsc_simulate_coda$Eta[[1]][[1]])[2]/nsites
	  #nlv<-0
	  
                i<-1
                
				
				#### Nimble Code for the 0 lv Hmsc model (except for priors):
			ModelData<-list(Y=as.vector(mod_Hmsc_simulate$Y))
			ModelConsts<-list(nvar=nvar,nspecies=nspecies,nobs=nspecies*nsites,
					sd10=10,
					beta0=rep(1:nspecies,each=nsites),
					sites=rep(1:nsites,times=nspecies),
					env=cbind(rep(1,nsites),mod_Hmsc_simulate$XData),
					### specify the number of latent variables on which we want to reorder species##not active in the nimbleCode below
					nlv=nlv.ref)
			

#calculating Rhat for each lambdapar and ordering it by level
diag_obs_Rhat<-lapply(1:ModelConsts$nlv,function(i){
temp_obs<-sapply(1:nspecies,function(x){gelman.diag(coda::as.mcmc.list(lapply(out.diag[,dimnames(out.diag[[1]])[[2]]==paste("Lambda1[sp_",x," (S",x,"), factor",i,"]",sep=""),drop=FALSE],function(x){coda::as.mcmc(x)})),autoburnin=FALSE,multivariate=FALSE)[[1]][1,2]})
names(temp_obs)<-1:nspecies
print(paste("Obs, level:",i))
print(sort(temp_obs,decreasing=TRUE))
sort(temp_obs,decreasing=TRUE)})



#calculating newspecies, the index of the species to put on the diagonal
newspecies<-NULL ### will be the number of the species -on the lambdascale, that is corresponding to species_order-  to put on the diagonal 
for (i in 1:ModelConsts$nlv)
{
	testsp<-as.double(names(which.max(diag_obs_Rhat[[i]])))
	while (is.element(testsp,newspecies))
	{
		diag_obs_Rhat[[i]][as.character(testsp)]<-0
		testsp<-as.double(names(which.max(diag_obs_Rhat[[i]])))
	}
	#if (testsp!=i)
	{
		newspecies<-c(newspecies,testsp)
	}

}
newspecies_3<-newspecies
})






