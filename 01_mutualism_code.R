rm(list=ls())


require(statmod)

#1. Web of interaction as matrix
#2. species trait variances
#3. m: species mean trait values that do not evolve
cutoff <- function(x) ifelse(x<1, (1*(x>0))*(x*x*x*(10+x*(-15+6*x))), 1)

traitCV <- function(mus,N){
  mus[N==0] <- NA
  v <- sort(mus[is.na(mus)==FALSE])
  d <- diff(v)
  cvs<-sd(d)/mean(d)

  return(cvs)
}


trait.matching<-function(mA,mP,adj.mat_1,gamma){
  tm<-numeric()
  for(i in 1:nrow(adj.mat_1)){
    tm[i] <- mean(adj.mat_1[i,]*exp(-(mA-mP[i])^2)/gamma)
    
  }
  return(tm=mean(tm))
}


plot_snapshot <- function(Na, Np, m, sigma, moment=0, limits=c(-0.8, 0.8), res=1001) {
  Sa <- length(Na) ## number of species
  Sp <- length(Np)
  ma<- m[1:(Sa)]
  mp<- m[Sa+1:Sp]
  sigma_a <-sigma[1:(Sa)]
  sigma_p <- sigma[Sa+1:Sp]
  traitaxis <- seq(limits[1], limits[2], l=res) ## sampling the trait axis
  #snap <- dat %>% filter(time==moment) %>% select(-time) ## time = moment
  traits_a <- expand.grid(species=1:Sa, trait=traitaxis) %>% as_tibble ## trait table
  traits_p <- expand.grid(species=Sa+1:Sp, trait=traitaxis) %>% as_tibble ## trait table
  
  traits_a["density"] <- 0 ## add column for population densities
  traits_p["density"] <- 0
  
  for (i in 1:Sa) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_a$density[(traits_a$species==i)] <- Na[i]*
      dnorm(traits_a$trait[(traits_a$species==i)], ma[i], sigma_a[i]) ## times density
  }
  traits_a$density[traits_a$density<max(traits_a$density)/1e3] <- NA
  
  for (i in (Sa+1):(Sa+Sp)) {
    #v <- snap %>% filter(species==i) %>% select(n, m, sigma)
    traits_p$density[(traits_p$species==i)] <- Np[-(Sa+Sp)+i+2]*dnorm(traits_p$trait[(traits_p$species==i)], 
                                                           mp[-(Sa+Sp)+i+2], sigma_p[-(Sa+Sp)+i+2]) ## times density
  }
  traits_p$density[traits_p$density<max(traits_p$density)/1e3] <- NA
  
  #landscape <- tibble(trait=traitaxis, r=0) %>% ## for plotting intrinsic rates
  #  mutate(r=1-trait^2/dat$theta[1]^2) %>% ## phenotype-specific growth rates
  #  mutate(r=ifelse(r<=0, NA, r)) %>%
   # mutate(r=r*max(traits$density, na.rm=TRUE)) ## scale growth rates for plotting
  traits<-data.frame(rbind(traits_a,traits_p), 
                     species_group=c(rep("Animals", nrow(traits_a)),
                                                 rep("Plants", nrow(traits_p))))
             
  ggplot(traits) + ## generate plot
    geom_line(aes(x=trait, y=density, colour=factor(species)), na.rm=TRUE) +
    geom_ribbon(aes(x=trait, ymin=0, ymax=density, fill=factor(species)),
                alpha=0.15, colour=NA)+facet_wrap(.~species_group, nrow = 2)+
    #geom_line(data=landscape, aes(x=trait, y=r), linetype="dashed",
     #         colour="darkred", alpha=0.5, na.rm=TRUE) +
    scale_x_continuous(name="trait value", limits=limits) +
    scale_y_continuous(name="density", limits=c(0, NA)) +
    theme(legend.position="none") %>%
    return }




gausquad.animals<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.animal){
  
  
  temp2<-dat2<-x2<-x3<-array(dim=c(points))
  if(mat == 0){
   return(list(G= 0, B = 0))
}
  else if(mat == 1){
  #nodes oir points in the abscissa where the integral will be evaluated numerically
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z'
  z2<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z''
  
  #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$ma,sigma =sigma$sa)$weights #pi(z')
  w2<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$mp,sigma =sigma$sp)$weights #pj(z'')
  
  
  #for the pairwise model however there are only two species interacting and hence i and j
  #or in other words the integral goes over z and z'
  for (i in 1: points){
    
    
    temp2[i]<- sum(np*(mut.strength/degree.animal)*exp(-(z1[i]- z2)^2/w^2)/(1+h*np*(mut.strength/degree.animal)*exp(-(z1[i]-z2)^2/w^2))*w2*w1[i])
    
      #sum(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])
    x2[i]<- sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
    dat2[i]<- sum((z1[i]-m$ma)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]))
    
     # x3[i] <-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
      
  }
  return(list(G= sum(temp2), 
              B = sum(dat2*np*(mut.strength/degree.animal)/(1+x2*(mut.strength/degree.animal)*np)) ))
  
  }
}


gausquad.plants<-function(m,sigma,w,h,np,na,mut.strength,points,mat,degree.plant){
  
  temp2<-dat2<-x3<-x4<-array(dim=c(points))
  
  if(mat == 0){
    
    return(list(G= 0, 
                B = 0))
    
  }
  else if (mat==1){
  #nodes oir points in the abscissa where the integral will be evaluated numerically
  z1<-gauss.quad.prob(points, dist = "normal", mu=m$mp, sigma =sigma$sp)$nodes #z'
  z2<-gauss.quad.prob(points, dist = "normal", mu=m$ma, sigma =sigma$sa)$nodes #z''
  
  #weights of the gaussian distribution given by mean trait value mu_i and its variance \sigma_i
  w1<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$mp,sigma =sigma$sp)$weights #pi(z')
  w2<-gauss.quad.prob(points, dist = "normal", 
                      mu=m$ma,sigma =sigma$sa)$weights #pj(z'')
  
  
  #for the pairwise model however there are only two species interacting and hence i and j
  #or in other words the integral goes over z and z'
  for (i in 1: points){
    temp2[i]<- sum(na*(mut.strength/degree.plant)*exp(-(z1[i]- z2)^2/w^2)/(1+h*na*(mut.strength/degree.plant)*exp(-(z1[i]-z2)^2/w^2))*w2*w1[i])
      
      #sum(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i])
    x3[i]<-sum(h*exp(-(z1[i]-z2)^2/w^2)*w2)
    dat2[i]<- sum((z1[i]-m$mp)*(exp(-(z1[i]- z2)^2/w^2)*w2*w1[i]))
    
   # x4[i]<-sum(h*w2*exp(-(z1[i]-z2)^2/w^2))
  }
#h=3
  return(list(G= sum(temp2), 
              B = sum(dat2*na*(mut.strength/degree.plant)/(1+x3*(mut.strength/degree.plant)*na))))
  }
  
}


eqs <- function(time, state, pars) {
  A <- dim(pars$matrix)[2]  ## number of animal species
  P <-dim(pars$matrix)[1]
  a <- state[1:A] ## species densities of animals
  p <- state[(A+1):(A+P)] ## species densities of plants
  s <- pars$sig ## species' trait standard deviations
  ma<-state[(A+P+1):(A+P+A)]
  mp<-state[(A+P+A+1):(A+P+A+P)]
  ## define g, where g[i] is the selection pressure on species i from growth
  alpha.a<-pars$Amatrix ## alpha matrix
  alpha.p<-pars$Pmatrix
  horder_a<-horder_p<-numeric()
  halpha_a<-pars$halpha_a
  halpha_p<-pars$halpha_p
  dt<-0.1
  #w<-pars$w
  aij<-bij<-matrix(0, nrow=A,ncol=P) 
  aji<-bji<-matrix(0, nrow=P,ncol=A) 
  Na<-muA<-matrix(0, nrow = time, ncol = A )
  Np<-muP<-matrix(0, nrow =time, ncol = P )
  Np[1,]<-state[(A+1):(A+P)]
  Na[1,]<-state[1:A]
  muA[1,]<-ma
  muP[1,]<-mp
  aj<-bj<-ai<-bi<-numeric()
  for (t in 1:(time-1)){
  for(r in 1:A){
    for(l in 1:P){
      #
      m.temp<-list(ma=muA[t,r],mp=muP[t,l])
      sigma1<-list(sa=s[r],sp=s[(A)+l])
      temp1<-gausquad.animals(m=m.temp,sigma=sigma1,w=pars$w,h=.5,np=Np[t,l],na=Na[t,r],
                             mut.strength=pars$mut.strength, points=7
                             ,mat=pars$matrix[l,r],degree.animal = pars$dganimals[r])
      aij[r,l] <-temp1$G
      bij[r,l]<-temp1$B
      
    }
    ai[r]<-sum(aij[r,])
    bi[r]<-sum(bij[r,])
    horder_a[r] <- t(Np[t,])%*%halpha_a[r,,]%*%Na[t,]
    
  }
  for(k in 1:P){
    for(m in 1:A){
      m2.temp<-list(ma=muA[t,m],mp=muP[t,k])
      sigma2<-list(sa=s[m],sp=s[(A+k)])
      temp2<-gausquad.plants(m=m2.temp,sigma=sigma2,w=pars$w,h=.5,np=Np[t,k],na=Na[t,m],
                             mut.strength=pars$mut.strength,
                             points=7,mat=pars$matrix[k,m], degree.plant =pars$dgplants[k])
      aji[k,m] <-temp2$G
      bji[k,m]<-temp2$B
    }
    aj[k]<-sum(aji[k,])
    bj[k]<-sum(bji[k,])
    
    horder_p[k] <- t(Na[t,])%*%halpha_p[k,,]%*%Na[t,]
  }
    
    
  #print(t)
    #a*sum(N[t,,i]*Disp.mat[k,])*dt - a*N[t,k,i]*dt +rnorm(1,0,0.02)*dt*N[t,k,i]
    
    if (pars$model == "HOI"){
    Na[t+1, ]<-Na[t,] + Na[t,]*(pars$ba-alpha.a%*%Na[t,]+ai+horder_a)*dt #+ rnorm(A, 0,sd=0.0)*Na[t,]*dt## density eqs
    Np[t+1, ]<-Np[t,] + Np[t,]*(pars$bp-alpha.p%*%Np[t,]+aj+horder_p)*dt #+ rnorm(P, 0,sd=0.0)*Np[t,]*dt ## trait mean eqs
    muA[t+1, ]<-muA[t,] +pars$h2[1:A]*(bij%*%Np[t,])*dt + rnorm(A, 0,sd=0.0)*dt ## trait mean eqs
    muP[t+1, ]<- muP[t,]+pars$h2[(A+1):(A+P)]*(bji%*%Na[t,])*dt + rnorm(P, 0,sd=0.0)*dt ## trait mean eqs
  
    Na[t+1,which(Na[t+1,] < 1e-4)]<-0
    Np[t+1,which(Np[t+1,] < 1e-4)]<-0
    } else if (pars$model == "pairwise"){
      Na[t+1, ]<-Na[t,] + Na[t,]*(pars$ba-alpha.a%*%Na[t,]+ai)*dt #+ rnorm(A, 0,sd=0.0)*Na[t,]*dt## density eqs
      Np[t+1, ]<-Np[t,] + Np[t,]*(pars$bp-alpha.p%*%Np[t,]+aj)*dt #+ rnorm(P, 0,sd=0.0)*Np[t,]*dt ## trait mean eqs
      muA[t+1, ]<-muA[t,] +pars$h2[1:A]*(bij%*%Np[t,])*dt + rnorm(A, 0,sd=0.0)*dt ## trait mean eqs
      muP[t+1, ]<- muP[t,]+pars$h2[(A+1):(A+P)]*(bji%*%Na[t,])*dt + rnorm(P, 0,sd=0.0)*dt ## trait mean eqs
      
      Na[t+1,which(Na[t+1,] < 1e-4)]<-0
      Np[t+1,which(Np[t+1,] < 1e-4)]<-0
      
    }

  } 
#ts.plot(Na,ylim=c(0,5))
#ts.plot(muA[1:100,])
  ## return equations by first flattening them back into a single vector
  output= list(Plants = Np[1:time,],Animals=Na[1:time,], Plant.trait = muP[1:time,], 
               Animal.trait=muA[1:time,],aij=aij,ba=pars$ba)
  return(output)
}

Mcommunity = function(iter, time, ...){
  set.seed(rnorm(1,as.numeric(Sys.time())-Sys.getpid(),10000)) 
  init = time
  replicate =try(eqs(time=init, ...))
  replicate$start = init
  replicate$iter = iter
  return(replicate)
}


## Organize simulation results into tidy table
## Input:
## - sol: output produced by the function ode()
## - pars: list of parameters, with the following elements:
##         $w: width of competition kernel
##         $theta: width of intrinsic growth function
##         $sigma: vector of species' intraspecific standard deviations
##         $h2: vector of heritabilities
## Output:
## - a tibble with columns: time, species, n (density), m (trait mean),
##   sigma, w, theta, and h2
organize_results <- function(sol, pars) {
  S <- length(pars$sigma) ## number of species
  A<-dim(pars$matrix)[2] # no. of animals
  P<-dim(pars$matrix)[1] # no. of plants
  temp<- sol %>% as.data.frame %>% as_tibble ## convert to tibble
  ## name the first column "time"
  temp<- temp %>% filter(time >= pars$cutoff.time)
  names(temp)[2:(A+1)] <- paste0("N_", 1:(A)) ## name abundance columns (n_k)
  names(temp)[1] <- "time"
  names(temp)[(A+2):(A+1+P)] <- paste0("P_", 1:P) ## name trait mean columns
  temp <- temp %>%
    gather("variable", "v", 2:ncol(temp)) %>% ## normalize the data
    separate(variable, c("type", "species"), sep="_") %>%
    #spread(type, v) %>% ## separate columns for animal densities n and plant densities m
    dplyr::select(time, type, species,v) %>% ## rearrange columns
    mutate(species=as.integer(species), Mut_strength=pars$mut.strength,
           Nestedness=pars$nestedness, Connectance=pars$C,
           theta=pars$theta,Web.name=pars$web.name) ## add params
  return(as_tibble(temp))
}


## Plot time series of densities, time series of trait values, and
## snapshot of the trait distributions at time = moment
## Input:
## - dat: data generated by organize_results()
## - moment: time at which trait distribution should be plotted
## - limits: a vector of two entries (x_low, x_high) for the x-axis limits
## - res: number of evenly spaced sampling points along the trait axis
##               for the trait distribution plot
## Output:
## - a ggplot2 plot with three panels in one column: abundance time series,
##   trait value time seties, and snapshot of trait distribution
plot_all <- function(dat, moment=0, limits=c(-0.6, 0.6), res=1001) {
  plot_grid(plot_density(dat), ncol=1, align="hv") %>%
    return
}
## Plot species densities through time
## Input:
## - dat: data generated by organize_results()
## Output:
## - a ggplot2 plot

plot_density<- function(dat) {
  dat %>%
    ggplot +
    geom_line(aes(x=time, y=v, colour = factor(species))) +
    scale_y_continuous(name="population density", limits=c(0, NA)) +
    theme(legend.position="none") + facet_wrap(.~type) %>%
    return
}



