rm(list=ls())
source('~/Dropbox/EAWAG PostDoc/05_Spatial_mutualism/Spatial_mutualism/01_mutualism_code.R')
require(deSolve) ## for integrating ordinary differential equations
require(tidyverse) ## for efficient data manipulation & plotting
require(cowplot) ## for arranging plots in a grid
library(dplyr)
library(readr)


theme_set(theme_bw()) 


adj.mat<-function(data){
  #dat <- paste('network.csv',sep='')
  d <- read.csv(file=data,header=FALSE )
  dat<-as.matrix(d)
  dat[dat > 0] = 1
  dat<-apply(dat,2,as.numeric)
  return(dat)}



mat.comp<-function(matrix,degree.animals,degree.plants){
  Aspecies<- dim(matrix)[2]
  Plantspecies<- dim(matrix)[1]
  
  Amatrix<-matrix(runif(Aspecies^2, 0.001, 0.05), nrow=Aspecies, ncol = Aspecies)
  diag(Amatrix)<-1
  #diag(Amatrix)<-  diag(Amatrix) #/degree.animals
  
  Pmatrix<-matrix(runif(Plantspecies^2, 0.001, 0.05), nrow=Plantspecies, ncol = Plantspecies)
  diag(Pmatrix)<-1
  #diag(Pmatrix)<-2-  diag(Pmatrix)/degree.plants
  out<-return(list(Amatrix=Amatrix,Pmatrix=Pmatrix))
  
}

# measures connectance of a web network
Connectance<-function(web)
{
  return(sum(web)/(ncol(web)*nrow(web)))}


#computes the raw NODF
#input: web = mutualistic network
#output: raw NODF of the given network
nestedness_NODF <- function(web){
  web[web > 0] = 1
  SA <- nrow(web)
  SP <- ncol(web)
  N <- t(web) %*% web
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SP,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SP)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele == 0] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)] = 1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n1 <- sum(nes)
  
  N <- web %*% t(web)
  num <- N
  num[lower.tri(num,diag=TRUE)]=1
  den <- (matrix(1,nrow=SA,ncol=1)*diag(N))%*%matrix(1,nrow=1,ncol=SA)
  dele <- den - t(den)
  dele[lower.tri(dele,diag=TRUE)] <- 1
  num[dele ==0 ] <- 0
  den <- pmin(den,t(den))
  den[lower.tri(den,diag=TRUE)]=1
  nes <- num/den
  nes[lower.tri(nes,diag=TRUE)] = 0
  nes[is.na(nes)] <- 0
  n2 <- sum(nes)
  out <- 2*(n1 + n2) / (SA*(SA-1)+SP*(SP-1))
  return(out)
}


#ex_webs<-c("datasets/M_PL_011.csv")


#new_fact<-fact %>% filter(web == ex_webs)

# reading all the datasets
# calculating nestedness and connectance
mydir = 'datasets'
myfiles = list.files(path=mydir, pattern="*.csv", full.names=TRUE)
#myfiles<-myfile
newfiles<-myfiles[1:60]

datasets/M_PL_061_33.csv
#newfiles<-ex_webs
fact<- expand.grid(`Strength_mutualism`=c( 0.75), 
                   `web` ="datasets/M_PL_061_33.csv",
                   #`initial.trait`=ru
                   `individual.variation` = c("high","low"),
                   `interaction`=c("HOI","pairwise"),
                   `random_seed`=4127+(1:12)*100) %>%
  as_tibble %>%
  mutate(`Biomass.animal`=0,
         `Biomass.plant` =0,
         `plant.richness`=0,
         `animal.richness`=0,
         `richness`=0,
        `trait.matching` = 0)
model.t<-list()

set.seed(1234)


for(r in 1:nrow(fact)){

  g<-adj.mat(myfiles[which(myfiles == fact$web[r])]) #network web names
  g<-g[-1,-1] 
  Aspecies<- dim(g)[2] # no of animal species
  Plantspecies<- dim(g)[1] # no of plant species
  degree.animals<-degree.plants<-numeric()
  
  #degree of plants and anichmals
  for(i in 1:Plantspecies){
    degree.plants[i]<-sum(g[i,])} # degree of plants
  for(j in 1:Aspecies){
    degree.animals[j]<-sum(g[,j]) # degree of animals
  }
  
  ##control loop for selecting whether variation is high or low
  if(fact$individual.variation[r] == "low"){
    sig <-runif((Aspecies+Plantspecies),0.0001,0.005)
    } else if(fact$individual.variation[r] == "high"){
    sig <-runif((Aspecies+Plantspecies),0.05,0.1)}
  
  #higher-order interaction matrix for animals (intra+interHOIs)
  halpha_a<-array(NA,dim=c(Aspecies,Plantspecies,Aspecies))
  for (i in 1:Aspecies){
    halpha_a[i,,]<- rnorm(Aspecies*Plantspecies, 0, 0.01)
  }
  
  #higher-order interaction matrix for plants (intra+interHOIs)  
  halpha_p<-array(NA, dim=c(Plantspecies,Aspecies,Aspecies))
  for (i in 1:Plantspecies){
    halpha_p[i,,]<- rnorm(Aspecies*Aspecies, 0, 0.01)
    }           
                    
  h2<-runif((Aspecies+Plantspecies),0.4,0.4)
  
  ## vector of species trait standard deviations
  N <- runif( (Aspecies+Plantspecies) , 1,1)  ## initial species densities
  nainit<- N[1:Aspecies]
  npinit<-N[(Aspecies+1): (Aspecies+Plantspecies)]
  muinit <-runif((Aspecies+Plantspecies), -1,1)
  mainit<-muinit[1:Aspecies]
  mpinit<-muinit[(Aspecies+1): (Aspecies+Plantspecies)]
  
  Amatrix<-mat.comp(g,degree.animals,degree.plants)$Amatrix
  Pmatrix<-mat.comp(g,degree.animals,degree.plants)$Pmatrix
  gamma=fact$Strength_mutualism[r]#fact_lessvar$Strength_mutualism[r]
  mut.strength<-5
  nestedness<-nestedness_NODF(g)
  C<-Connectance(g)
  web.name<-myfiles[r]
  ba<-runif(Aspecies, -0.5,-0.05)
  bp<-runif(Plantspecies,-0.5,-0.05)
  dganimals<-degree.animals
  dgplants<-degree.plants
  
  
  ic <-c(nainit, npinit, mainit,mpinit)
  
  
  params <- list(time=time,matrix=g,sig=sig,Amatrix=Amatrix,
                 Pmatrix=Pmatrix,w=gamma,
                 mut.strength=mut.strength,m=muinit,C=C,nestedness=nestedness,
                 web.name=web.name,h2=h2, ba=ba,bp=bp,dganimals=dganimals,
                 dgplants=dgplants,halpha_p=halpha_p,halpha_a=halpha_a, model= fact$interaction[r])
  
  
  
  start.time =1500
  model.t<-lapply(1, Mcommunity,time=start.time,state=ic,
                  pars=params)

  
  plot_snapshot(Na = model.t[[1]]$Animals[1500,],
                Np = model.t[[1]]$Plants[1500,],
                m = c(model.t[[1]]$Animal.trait[1500,], model.t[[1]]$Plant.trait[1500,]),
                sigma =sig, moment=0, limits=c(-0.8, 0.8), res=1001)
  
  
  pbiomass<-sum(colMeans(model.t[[1]]$Plants[600:1500,]))
  abiomass<-sum(colMeans(model.t[[1]]$Animals[600:1500,]))
  
  fact$plant.richness[r] <-length(which(model.t[[1]]$Plants[1500,] > 0))
  fact$animal.richness[r] <-length(which(model.t[[1]]$Animals[1500,] > 0))
  fact$richness[r] <- fact$animal.richness[r] +   fact$plant.richness[r]
  fact$Biomass.animal[r] = abiomass
  fact$Biomass.plant[r] = pbiomass
  fact$CV.trait[r] = traitCV(mus = c(model.t[[1]]$Animal.trait[1500,], model.t[[1]]$Plant.trait[1500,]),
                             N = c(model.t[[1]]$Animals[1500,],model.t[[1]]$Plants[1500,]))
  fact$trait.matching[r] = trait.matching(mA = model.t[[1]]$Animal.trait[1500,],
                                         mP = model.t[[1]]$Plant.trait[1500,],
                                         adj.mat_1=g,
                                         gamma=gamma)

  print(r)
}
save(fact, file="Mutualism_data.RData")
