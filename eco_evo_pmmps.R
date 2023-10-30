#### Preamble ####
rm(list=ls())
library(MASS)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(dplyr)
library(vioplot) 
pal_d <- brewer.pal(8,"Dark2")
pal_v <- viridis(3)

## Some useful functions
nnzero <- function(x){
  length(which(x!=0))
}

mysample <- function(x,size,replace=T,prob=NULL){
  if(length(x)==1) { rep(x,times=round(size)) }
  else { sample(x,size=round(size),replace=replace,prob=prob)}
}

## is.mature: Quick function to check whether an individual has reached age of maturity.
## Arguments:
#  x  : A vector of ages, can include NAs.
#  a  : Age at maturity.
#  Outputs: A boolean vector of length x.
is.mature <- function(x,a){
  x[is.na(x)] <- -1
  out <- ifelse(x>=a,TRUE,FALSE)
  return(out)
}

## ddfun: Quickly plot survival as a function of density for arbitrarily many m2 parameters
## stars: a matrix of dimensions N x reps of population densities, to plot along the lines.
ddfun <- function(m2,type="exponential",title=NA,stars=NA,xmax=2,pal=pal_v){
  xs <- seq(0,xmax,by=0.002)
  survs <- numeric(length(xs))
  plot(xs,rep(1,length(xs)),type="n",ylim=c(0,1),col=pal[1],xlab=bquote(italic(n[i]/K)),cex.lab=1.2,ylab="Survival probability",main=title)
  if(length(m2)>1){
    for(i in 1:length(m2)){
      if(type=="exponential"){
        survs <- dexp(xs,m2[i])/m2[i]
        lines(xs,survs,col=pal[i],lwd=2)
        if(sum(stars[i,],na.rm=T)>0){
          for(j in 1:length(stars[1,])){ # For each replicate, plot an open circle
            points(stars[i,j],dexp(stars[i,j],m2[i])/m2[i],pch=1,col=alpha(pal[i],0.5),cex=1.5)
          }
          points(mean(stars[i,]),dexp(mean(stars[i,]),m2[i])/m2[i],pch=20,col=pal[i],cex=1.8)
        }
      } else if(type=="linear"){
        survs <- 1-xs*m2[i]
        lines(xs,survs,col=pal[i],lwd=2)
        if(sum(stars,na.rm=T)>0){
          points(stars[i],1-xs*m2[i],pch=20,col=pal[i])
        }
      }
      #print(survs[101])
    }
  }
  abline(v=1,lty=2)
  legend("bottomleft",col=pal,pch=rep(20,length(m2)),legend=as.fractions(m2),title=expression(italic(gamma)[n]),bty="n",cex=1,pt.cex=1.3)
}
#####

#### Functions for wrangling and plotting results ####

### Popstats ###
##  Function to summarize output from simulations into a plottable format (takes some time for very long simulations)
##  Arguments:
#   local    : A simulation result (created locally)
#   saved    : A simulation result (from saved file)
#   max.age  : Number of age classes to bin individuals in for age distribution (individuals of age > max.age are binned into age=max.age)
##  Outputs
#   A list of two elements:
#   [[1]] EVOL: Storage array of format:
#   Row 1:      Mean breeding values
#   Row 2:      Sd breeding values
#   Row 3:(N+2):Frequencies of destination genes
#   Row (N+3):(N+2+max.age):Age distribution
#   per population for each Patch (2nd dim), recorded year (3rd dim) and Replicate (4th dim)
#   [[2]] LIAB: Storage array of format:
#   Row 1:      Mean liability
#   Row 2:      Sd liability
#   Row 3:      Fraction L>0
#   per population for each patch (2nd dim), recorded year (3rd dim) and Replicate (4th dim)
#   [[3]] AGE: Storage array of subpopulation mean age per patch (1st dim), recorded year (2nd dim) and replicate (3rd dim).
#   [[4]] EXT: Storage array of format
#   Row 1: Year of first extinction (pre-breeding population census of 0)
#   Row 2: Most recently recorded mean BV prior to extinction.
#   Row 3: Most recently recorded sd BV prior to extinction.
#   per site (2nd dim) and replicate (3rd dim)
popstats <- function(local=NULL,saved=NULL,max.age=10){
  if(is.null(local)){
    output <- readRDS(saved)
  } else {
    output <- local
  }
  
  popstorage <- output[[1]]
  dyn <- output[[2]]
  N <- dim(dyn)[2] # Number of subpopulations
  K <- dim(popstorage)[1]/N # Max number of individuals per subpopulation (carrying capacity)
  t <- dim(popstorage)[3] # Number of recorded timesteps
  
  reps <- dim(dyn)[3] # Number of replicate simulations
  rec.int <- dim(dyn)[1]/(2*t) # Interval of recording
  years <- (1:t)*rec.int # The years at which recordings took place
  
  ## Create storage matrices ##
  # Mean and sd gene values per pop per recording
  evol <- array(NA,dim=c(N+2+max.age,N,length(years),reps)) # Dims: 1) Stats (see below). 2) Patches. 3) Years. 4) Reps
  liab <- array(NA,dim=c(3,N,length(years),reps)) # Rows: 1) Mean L 2) Sd L 3) Fraction L>0 - per patch, year, rep.
  age <- array(NA,dim=c(N,length(years),reps)) # Subpopulation mean age. Dims: 1) Patch, 2) Years, 3) Rep
  ext <- array(NA,dim=c(3,N,reps)) # Subpop extinction times (1st row), most recent mean (2nd row) and sd (3rd row) BVs pre-extinction, per patch (cols) and rep (3rd dim)
  # Stats format for evol:
  # 1: Mean breeding value
  # 2: Sd breeding value
  # 3,4,5 (or "3:(N+2)"): Frequencies of destination genes
  # 6-16 (or "(N+3):(N+2+max.age)"): Age distribution
  t <- Sys.time()
  for(r in 1:reps){
    for(i in 1:length(years)){
      popstorage[which(popstorage[,3,i,r]>max.age),3,i,r] <- max.age # Final age class is "10 or older".
      for(j in 1:N){
        inds <- (1:K)+(j-1)*K
        evol[1,j,i,r] <- mean(popstorage[inds,1,i,r],na.rm=T) # Mean breeding value
        evol[2,j,i,r] <- sd(popstorage[inds,1,i,r],na.rm=T) # Sd breeding value
        liab[1,j,i,r] <- mean(popstorage[inds,6,i,r],na.rm=T) # Mean liability
        liab[2,j,i,r] <- sd(popstorage[inds,6,i,r],na.rm=T) # Sd liability
        liab[3,j,i,r] <- sum(popstorage[inds,6,i,r]>0,na.rm=T)/sum(!is.na(popstorage[inds,6,i,r])) # Fraction of migrants
        evol[3:(N+2),j,i,r] <- c(tabulate(popstorage[inds,2,i,r]),rep(0,times=N-length(tabulate(popstorage[inds,2,i,r]))))/sum(!is.na(popstorage[inds,1,i,r])) #Or divide by dyn[(i*2)-1,j,r]  # Destination genes distribution
        #print(length(c(tabulate(popstorage[(1:K)+(j-1)*K,3,i,r]),rep(0,times=max.age-length(tabulate(popstorage[(1:K)+(j-1)*K,3,i,r]))))))
        #print(length((N+3):(N+2+max.age))) # These are just for checking that everything runs smoothly.
        evol[(N+3):(N+2+max.age),j,i,r] <- c(tabulate(popstorage[inds,3,i,r])/sum(!is.na(popstorage[inds,1,i,r])),rep(0,times=max.age-length(tabulate(popstorage[inds,3,i,r])))) # Age distribution
        evol[is.infinite(evol)] <- NA # Find a better solution for this!
        liab[is.infinite(liab)] <- NA 
        age[j,i,r] <- mean(popstorage[((j-1)*K+1):(j*K),3,i,r],na.rm=T)
      }
    }
    for(j in 1:N){
    # Use dyn instead  ext[r,j] <- which(is.na(statstmp[[2]][3,3,,1]))[min(which(diff(which(is.na(statstmp[[2]][3,3,,1])))==1))]
      summers <- seq(1,length(dyn[,1,1]),by=2)
      if(any(dyn[summers,j,r]==0)){
        ext[1,j,r] <- min(which(dyn[summers,j,r]==0)) # Year of first extinction
        ext[2,j,r] <- evol[1,j,max(which(years<ext[1,j,r])),r] # Most recently recorded mean breeding value
        ext[3,j,r] <- evol[2,j,max(which(years<ext[1,j,r])),r]
      }
    }
    ##
    
    print(r)
  }
  
  print(Sys.time()-t)
  return(list(evol,liab,age,ext))
}

### Popplots ###
## Function to plot simulation results.
## Arguments:
#  local    : A stats file (generated by popstats) in the local environment
#  saved    : A stats file (generated by popstats) that can be read from disk
#  mainfile : A simulation result (generated by sim), created locally or read from file.
#  nview    : Number of replicates to choose when plotting fewer than all replicates. Default 5.
#  plot     : Shortnames for the different plots wanted: "all" (default),
#             "indBVs","indLs","meanBVs","meanLs","popsize","endpopsize","winteruse","agedist","Dest","meancorr","beforecorr","aftercorr",
#             "EndBVs","EndBVviol","EndLs","EndLiabviol","EndMeans","EndSDs","EndDests","EndDensity","avgBVs","avgBVviol","avgLs","avgLiabviol","avgMeans","avgSDs","avgDests","avgDensity",
#             "preshockBVs","preshockBVviol","preshockLs","preshockLviol","preshockMeans","preshockSDs","preshockDests","preshockDensity".
#  pal      : Color palette to use: "vir" (viridis, default), or "dark" (RColorbrewer's Dark2). Both should be colorblind friendly.
#  legend   : Setting FALSE can suppress some more annoying legends: EndDests, (add more here as they are implemented)
popplots <- function(local=NULL,saved="",mainfile=NULL,nview=5,plot="all",pal="vir",legend=T,...){
  tt <- Sys.time()
  if(is.null(local)){
    stats <- readRDS(saved)
  } else {
    stats <- local
  }
  evol <- stats[[1]]
  liab <- stats[[2]]
  
  if(is.null(mainfile)){
    tmp <- strsplit(saved," stats")[[1]][1]
    output <- readRDS(paste0(tmp,".R"))
  } else {
    output <- mainfile
  }
  
  popstorage <- output[[1]]
  dyn <- output[[2]]
  N <- dim(dyn)[2] # Number of subpopulations
  K <- dim(popstorage)[1]/N # Max number of individuals per subpopulation (carrying capacity)
  t <- dim(popstorage)[3] # Number of recorded timesteps
  
  reps <- dim(dyn)[3] # Number of replicate simulations
  rec.int <- dim(dyn)[1]/(2*t) # Interval of recording
  years <- (1:t)*rec.int # The years at which recordings took place
  max.age <- dim(evol)[1]-N-2 # Number of age classes individuals are binned into (for age distribution)
  
  pal <- switch(pal,"vir"=viridis(N),"dark"=brewer.pal(N,"Dark2"))
  sites <- switch(as.character(N),
                  "2"=c("North","South"),
                  "3"=c("North","Middle","South"),
                  "4"=c("North","Middle-N","Middle-S","South"),
                  "5"=c("North","Middle-N","Middle","Middle-S","South"),
                  "6"=c("North","","Middle-N","Middle-S","","South"),
                  "7"=c("North","","Middle-N","Middle","Middle-S","","South"),
                  "8"=c("North","","Middle-N","","","Middle-S","","South"),
                  "9"=c("North","","Middle-N","","Middle","","Middle-S","","South"))
  if(N>8){
    pal <- c(pal,brewerpal(8,"Set2"))
  }
  
  # Other info that will be useful to generate
  ncorrs <- choose(N,2) # Calculate how many combinations of patches there are
  viewr <- sample(reps,min(nview,reps)) # If there are too many reps to visualize easily, only plot 5 of them.
  shocks <- output[[3]]$shock.size
  print(Sys.time()-tt)

  # Mean gene values from random 3 replicates over time. Line type is replicate. Shaded bands are population SD.
  if(any(plot=="all") | any(plot=="indBVs")){
    par(mfrow=c(2,2),mar=c(4,4,2,1))
    
    plot(1:t,rep(0,t),type="n",xlab="Year",xaxt="n",ylab="Breeding values \U00B1 within-pop. sd",ylim=c(-5,5),cex.lab=1.2)
    axis(1,at=seq(0,t,by=1000/rec.int),labels=c(0,years[which(years%%1000==0)]))
    viewr <- sample(reps,size=min(reps,3),replace=F)
    for(r in 1:length(viewr)){
      for(j in 1:N){
        lowerbounds <- evol[1,j,,viewr[r]] - evol[2,j,,viewr[r]]
        upperbounds <- evol[1,j,,viewr[r]] + evol[2,j,,viewr[r]]
        polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],0.1),border=NA)
        lines(1:t,evol[1,j,,viewr[r]],col=pal[j],lty=r)
      }
    }
    if(legend) legend("bottomleft",legend=paste("Site",sites[1:N]),lty=1,col=pal[1:N],bg=alpha("White",alpha=0.5))
    abline(h=0,lty=2)
  }
  
  # Mean liabilities from 3 random replicates over time. Line type is replicate. Shaded bands are population SD.
  if(any(plot=="all") | any(plot=="indLs")){
    plot(1:t,rep(0,t),type="n",xlab="Year",xaxt="n",ylab="Liabilities \U00B1 within-pop. sd",ylim=c(-5,5),cex.lab=1.2)
    axis(1,at=seq(0,t,by=1000/rec.int),labels=c(0,years[which(years%%1000==0)]))
    for(r in 1:length(viewr)){
      for(j in 1:N){
        lowerbounds <- liab[1,j,,viewr[r]] - liab[2,j,,viewr[r]]
        upperbounds <- liab[1,j,,viewr[r]] + liab[2,j,,viewr[r]]
        polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],0.1),border=NA)
        lines(1:t,liab[1,j,,viewr[r]],col=pal[j],lty=r)
      }
    }
    if(legend) legend("bottomleft",legend=paste("Site",sites[1:N]),lty=1,lwd=1,col=pal[1:N],bg=alpha("white",alpha=0.5))
    abline(h=0,lty=2)
  }
  
  #Mean gene values across replicates over time. Shaded bands are standard error of the mean across all reps
  if(any(plot=="all") | any(plot=="meanBVs")){
#    if(any(plot!="all") & any(plot=="meanBVs")){
 #     par(mfrow=c(2,2),mar=c(4,4,2,1))
  #  }
    sigma_e <- output[[3]]$sigma_e
    plot(1:t,rep(0,t),type="n",xlab="Year",xaxt="n",ylab="Mean breeding values \U00B1 among-pop. SE",ylim=c(-5,5),cex.lab=1.2)
    axis(1,at=seq(0,t,by=1000/rec.int),labels=c(0,years[which(years%%1000==0)]))
    for(j in 1:N){
      means <- apply(evol[1,j,,],FUN=mean,MARGIN=1,na.rm=T)
      lowerbounds <- means - apply(evol[1,j,,],MARGIN=1,FUN=sd,na.rm=T)/sqrt(reps) # Standard error of the mean
      upperbounds <- means + apply(evol[1,j,,],MARGIN=1,FUN=sd,na.rm=T)/sqrt(reps) # Standard of the mean
      lines(1:t,apply(evol[1,j,,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2)
      polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],alpha=0.1),border=NA)
    }
    ys <- seq(-3*sigma_e,3*sigma_e,by=0.01)
    xs <- dnorm(ys,0,sigma_e)
    lines(xs*t/4,ys)
    abline(h=0,lty=2,lwd=2)
    if(legend) legend("bottomleft",title="Site",legend=sites[1:N],col=pal[1:N],lty=1,lwd=2,bty="n",seg.len=1)
  }
  
  # Mean liabilities across replicates over time. Shaded bands are standard error of the mean across all reps
  if(any(plot=="all") | any(plot=="meanLs")){
    plot(1:t,rep(0,t),type="n",xlab="Year",xaxt="n",ylab="Mean liabilities \U00B1 among pop. SE",ylim=c(-5,5),cex.lab=1.2)
    axis(1,at=seq(0,t,by=1000/rec.int),labels=c(0,years[which(years%%1000==0)]))
    for(j in 1:N){
      means <- apply(liab[1,j,,],FUN=mean,MARGIN=1,na.rm=T)
      #lowerbounds <- apply(liab[1,j,,],MARGIN=1,FUN=quantile,na.rm=T,probs=0.25) # 50 % quantile
      #upperbounds <- apply(liab[1,j,,],MARGIN=1,FUN=quantile,na.rm=T,probs=0.75) # 50 % quantile
      lowerbounds <- means - apply(liab[1,j,,],MARGIN=1,FUN=sd,na.rm=T)/sqrt(reps) # Standard error of the mean
      upperbounds <- means + apply(liab[1,j,,],MARGIN=1,FUN=sd,na.rm=T)/sqrt(reps) # Standard of the mean
      lines(1:t,apply(liab[1,j,,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2)
      polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],alpha=0.1),border=NA)
    }
    sigma_e <- output[[3]]$sigma_e
    ys <- seq(-3*sigma_e,3*sigma_e,by=0.01)
    xs <- dnorm(ys,0,sigma_e)
    lines(xs*t/4,ys)
    abline(h=0,lty=2,lwd=2)
    if(legend) legend("bottomleft",title="Site",legend=sites[1:N],col=pal[1:N],lty=1,lwd=2,bty="n",seg.len=1)
  }
  
  # Population sizes over time, means and SDs across all reps.
  if(any(plot=="all") | any(plot=="popsize")){
    par(mfrow=c(2,1),mar=c(4,4,1,0.5))
    plot(1:(rec.int*t),seq(1,K,length.out=rec.int*t),type="n",ylab="Local population size",xlab="Year",ylim=c(0,max(dyn,na.rm=T)),xaxt="n",cex.lab=1.2)
    axis(1,at=seq(0,rec.int*t,by=1000)-1,labels=c(0,years[which(years%%1000==0)]))
    for(j in 1:N){
      lines(1:(rec.int*t),apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=1,lty=2) # Summer
      lines(1:(rec.int*t),apply(dyn[seq(2,(rec.int*t*2),by=2),j,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=1) # Winter
      lowerbounds.b <- apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=mean,MARGIN=1,na.rm=T) - apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=sd,MARGIN=1,na.rm=T) # Summer
      upperbounds.b <- apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=mean,MARGIN=1,na.rm=T) + apply(dyn[seq(1,(rec.int*t*2)-1,by=2),j,],FUN=sd,MARGIN=1,na.rm=T)
      polygon(c(1:(rec.int*t),rev(1:(rec.int*t))),c(lowerbounds.b,rev(upperbounds.b)),col=alpha(pal[j],alpha=0.1),border=NA)
      lowerbounds.w <- apply(dyn[seq(2,rec.int*t*2,by=2),j,],FUN=mean,MARGIN=1,na.rm=T) - apply(dyn[seq(2,rec.int*t*2,by=2),j,],FUN=sd,MARGIN=1,na.rm=T) # Winter
      upperbounds.w <- apply(dyn[seq(2,rec.int*t*2,by=2),j,],FUN=mean,MARGIN=1,na.rm=T) + apply(dyn[seq(2,rec.int*t*2,by=2),j,],FUN=sd,MARGIN=1,na.rm=T)
      polygon(c(1:(rec.int*t),rev(1:(rec.int*t))),c(lowerbounds.w,rev(upperbounds.w)),col=alpha(pal[j],alpha=0.1),border=NA)
    }
    abline(h=K,lty=5) # Breeding season carrying capacity
    if(legend) legend("topleft",legend=c("Breeding","Non-breeding",paste(sites[1:N],"site")),lty=c(2,1,rep(0,N)),pch=c(NA,NA,rep(15,N)),lwd=2,col=c("Black","Black",pal[1:N]),bg=alpha("white",alpha=0.5))
  }
  
  # Population sizes +- sd across reps, zoom in on last 75 time steps. Change hashtags for 100 time steps. 
  if(any(plot=="all") | any(plot=="endpopsize") && rec.int*t>100){
    ts <- seq(rec.int*t*2-149,rec.int*t*2,by=2)
    plot(1:75,seq(1,K,length.out=75),type="n",ylab="Local population size",xlab="Year",ylim=c(0,max(dyn,na.rm=T)),xaxt="n",cex.lab=1.2)
    axis(1,at=c(1,25,50,75),labels=((rec.int*t)-25*(3:0))) 
    for(j in 1:N){
      lines(1:75,apply(dyn[ts,j,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2,lty=2) # Summer
      lines(1:75,apply(dyn[ts+1,j,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2) # Winter
      lowersum <- apply(dyn[ts,j,],FUN=mean,MARGIN=1,na.rm=T) - apply(dyn[ts,j,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps) # Summer
      uppersum <- apply(dyn[ts,j,],FUN=mean,MARGIN=1,na.rm=T) + apply(dyn[ts,j,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps) # Summer
      lowerwin <- apply(dyn[ts+1,j,],FUN=mean,MARGIN=1,na.rm=T) - apply(dyn[ts+1,j,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps) # Winter
      upperwin <- apply(dyn[ts+1,j,],FUN=mean,MARGIN=1,na.rm=T) + apply(dyn[ts+1,j,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps) # Winter
      polygon(c(1:75,75:1),c(lowersum,rev(uppersum)),col=alpha(pal[j],alpha=0.1),border=NA)
      polygon(c(1:75,75:1),c(lowerwin,rev(upperwin)),col=alpha(pal[j],alpha=0.1),border=NA)
      #for(r in viewr){
      #  lines(1:75,dyn[ts,j,r],col=alpha(pal[j],0.1),lty=2) # Summer
      #  lines(1:75,dyn[ts+1,j,r],col=alpha(pal[j],0.1)) # Winter
      #}
    }
    abline(h=K,lty=5,lwd=2) # Breeding season carrying capacity
    abline(v=25,lty=3,lwd=2) # Shock time
    if(legend){
      legend("topright",legend=c("Breeding", "Non-breeding",paste(sites[1:N],"site")),lty=c(2,1,rep(0,N)),pch=c(NA,NA,rep(15,N)),lwd=2,col=c("Black","Black",pal[1:N]),bg=alpha("white",alpha=0.5))
    }
  }
  
  # Frequencies of individuals overwintering where. Each point is a replicate - max 5 random replicates shown.
  if(any(plot=="all") | any(plot=="winteruse")){
    par(mfrow=c(N,1),mar=c(3.5,4,0.5,0.5))
    means <- array(NA,dim=c(N,N,length(years),reps)) # Dims: From, to, timestep, rep
    for(f in 1:N){
      inds <- ((f-1)*K+1):(f*K)
      for(y in 1:t){
        for(r in 1:length(viewr)){
          means[f,,y,viewr[r]] <- c(tabulate(popstorage[inds,5,y,viewr[r]]),rep(0,times=N-length(tabulate(popstorage[inds,5,y,viewr[r]]))))/sum(!is.na(popstorage[inds,5,y,viewr[r]]))
          if(sum(means[f,,y,viewr[r]],na.rm=T)>1.001 | sum(means[f,,y,viewr[r]],na.rm=T)<0.999){
            print(paste("Overwintering frequencies do not sum to 1, population",f,", year",y*rec.int,", rep",viewr[r]))
            print(sum(means[f,,y,viewr[r]],na.rm=T),digits=10)
          }
        }
      }
      plot(1:t,rep(1,t),type="n",ylim=c(0,1),ylab=paste("Freq of winter locations from",sites[f],"pop"),xlab="",xaxt="n")
      axis(1,at=seq(1,t+1,by=1000/rec.int)-1,labels=c(0,years[which(years%%1000==0)]))
      for(j in 1:N){
        for(r in 1:reps){
          points(jitter(1:t),means[f,j,,r],col=alpha(pal[j],0.5),pty=r)
        }
      }
    }
    if(legend) legend("bottomleft",title="Non-breeding season site",legend=sites[1:N],col=pal,pch=rep(15,N),bg=alpha("white",0.5))
    mtext("Year",side=1,line=2.5)
  }
  
  # Age distribution from up to 5 random replicates.
  if(any(plot=="all") | any(plot=="agedist")){ 
    par(mfrow=c(length(viewr),1),mar=c(3.5,4.2,0.5,0.3),oma=c(0,0,1.5,0))
    for(r in 1:length(viewr)){
      plot(1:max.age,rep(0.1,max.age),type="n",xlab="",ylab="Frequency",cex.lab=1.2,cex.axis=1.2,ylim=c(0,max(evol[(N+3):(N+2+max.age),,,],na.rm=T)))
      for(j in 1:N){
        lines((1:max.age)+0.1*(j-1),evol[(N+3):(N+2+max.age),j,t,viewr[r]],col=pal[j],type="h",lwd=2)
        abline(v=mean(popstorage[((j-1)*K+1):(j*K),3,t,viewr[r]],na.rm=T),lty=2,col=pal[j])
      }
      text(max.age/2,0.25,paste("Replicate",viewr[r]),cex=1.5)
    }
    mtext("Age class",side=1,line=2)
    mtext("   Age distribution, post shock",side=3,line=-0.3,cex=1.2,outer=T)
    if(legend) legend("topright",title="Home site",legend=sites[1:N],col=pal,pch=rep(15,N),bg=alpha("white",0.5))
  }
  
  # Destination gene evolution and migrant fraction over time.
  if(any(plot=="all") | any(plot=="Dest")){
    par(mfrow=c(N,1),mar=c(3.5,4,0.5,0.5))
    for(i in 1:N){
      plot(1:t,rep(1,t),type="n",ylim=c(0,1),ylab=paste("Freq of Dest. gene in",sites[i],"pop"),xlab="",xaxt="n",cex.lab=1+(4-N)/10)
      axis(1,at=seq(1,t+1,by=1000/rec.int)-1,labels=c(0,years[which(years%%1000==0)]))
      for(j in 1:N){
        lines(1:t,apply(evol[3+(j-1),i,1:t,],FUN=mean,MARGIN=1,na.rm=T),col=pal[j],lwd=2)
        lowerbounds <- apply(evol[3+j-1,i,1:t,],FUN=mean,MARGIN=1,na.rm=T) - apply(evol[3+j-1,i,1:t,],FUN=sd,MARGIN=1,na.rm=T)/reps
        upperbounds <- apply(evol[3+j-1,i,1:t,],FUN=mean,MARGIN=1,na.rm=T) + apply(evol[3+j-1,i,1:t,],FUN=sd,MARGIN=1,na.rm=T)/reps
        polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha(pal[j],alpha=0.1),border=NA)
      }
      lines(1:t,apply(liab[3,i,1:t,],FUN=mean,MARGIN=1,na.rm=T),col="Black",lty=2,lwd=2) # Proportion migrant
      lowerbounds <- apply(liab[3,i,1:t,],FUN=mean,MARGIN=1,na.rm=T) - apply(liab[3,i,1:t,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps)
      upperbounds <- apply(liab[3,i,1:t,],FUN=mean,MARGIN=1,na.rm=T) + apply(liab[3,i,1:t,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(reps)
      polygon(c(1:t,rev(1:t)),c(lowerbounds,rev(upperbounds)),col=alpha("Black",alpha=0.1),border=NA)
    }
    mtext("Year",side=1,line=2.5)
    legend("topleft",title="Home site",legend=sites[1:N],col=pal,pch=rep(15,N),bg=alpha("white",0.5),cex=1+(1-N)/20)
    if(legend) legend("topright",legend=c("Allele frequencies","Proportion migrant"),lty=1:2,lwd=2,bg=alpha("white",0.5))
  }
  
  # Correlations between population sizes
  if(any(plot=="all") | any(plot=="meancorr")){
    par(mfrow=c(3,1),mar=c(3,4,1.5,0.5))
    # Same kind of loop as for networks, plotting all combinations
    scorrs <- wcorrs <- matrix(NA,ncorrs,reps)
    names <- character(ncorrs)
    i <- counter <- 1
    summers <- seq(1,rec.int*t-1,by=2)
    winters <- summers+1
    plot(1:ncorrs,rep(0,ncorrs),type="n",xaxt="n",xlim=c(ifelse(N>3,1,0.5),ifelse(N>3,ncorrs,ncorrs+0.5)),ylim=c(-1,1),ylab="Correlation",xlab="",main="All years",bty="L")
    while(i<N){
      for(j in (i+1):N){
        polygon(c(counter-0.3,counter+0.3,counter+0.3,counter-0.3),c(-1.2,-1.2,1.2,1.2),col=alpha("Light grey",alpha=min(1,shocks[i]+shocks[j])),border=NA)
        scorrs[counter,] <- diag(cor(log(dyn[summers,i,]),log(dyn[summers,j,])))
        wcorrs[counter,] <- diag(cor(log(dyn[winters,i,]),log(dyn[winters,j,])))
        names[counter] <- paste0(LETTERS[i],LETTERS[j])
        vioplot(scorrs[counter,],add=T,at=counter-0.1,col=alpha("Dark green",0.2),wex=0.25,border="Dark green",lineCol=NA)
        vioplot(wcorrs[counter,],add=T,at=counter+0.1,col=alpha("Blue",0.2),wex=0.25,border="Blue",lineCol=NA)
        #points(counter+0.1*(1:length(viewr))-0.1*mean(c(1,length(viewr))),scorrs[counter,viewr],pch=2,col="Dark green")
        #points(counter+0.1*(1:length(viewr))-0.1*mean(c(1,length(viewr))),wcorrs[counter,viewr],pch=4,col="Blue")
        counter <- counter+1
      }
      i <- i+1
    }
    axis(1,at=1:ncorrs,labels=names,cex=1.2)
    abline(h=0,lty=2)
    #legend("bottomleft",legend=c("Summer","Winter"),col=c("Dark green", "Blue"),pch=c(2,4))
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),col=c("Dark green","Blue"),pt.bg=c(alpha("Dark green",0.2),alpha("Blue",0.2)),pch=22)
  }
  
  # Correlations between population sizes before shock
  if(any(plot=="all") | any(plot=="beforecorr")){
    # Same kind of loop as for networks, plotting all combinations
    scorrs <- wcorrs <- matrix(NA,ncorrs,reps)
    i <- counter <- 1
    summers <- seq(rec.int*(t-2)+1,rec.int*(t-1)-1,by=2)
    winters <- summers+1
    plot(1:ncorrs,rep(0,ncorrs),type="n",xaxt="n",xlim=c(ifelse(N>3,1,0.5),ifelse(N>3,ncorrs,ncorrs+0.5)),ylim=c(-1,1),ylab="Correlation",xlab="",main="50 years before shock",bty="L")
    while(i<N){
      for(j in (i+1):N){
        polygon(c(counter-0.3,counter+0.3,counter+0.3,counter-0.3),c(-1.2,-1.2,1.2,1.2),col=alpha("Light grey",alpha=min(1,shocks[i]+shocks[j])),border=NA)
        scorrs[counter,] <- diag(cor(log(dyn[summers,i,]),log(dyn[summers,j,])))
        wcorrs[counter,] <- diag(cor(log(dyn[winters,i,]),log(dyn[winters,j,])))
        names[counter] <- paste0(LETTERS[i],LETTERS[j])
        vioplot(scorrs[counter,],add=T,at=counter-0.1,col=alpha("Dark green",0.2),wex=0.25,border="Dark green",lineCol=NA)
        vioplot(wcorrs[counter,],add=T,at=counter+0.1,col=alpha("Blue",0.2),wex=0.25,border="Blue",lineCol=NA)
        #points(counter+0.1*(1:length(viewr))-0.1*mean(c(1,length(viewr))),scorrs[counter,viewr],pch=2,col="Dark green")
        #points(counter+0.1*(1:length(viewr))-0.1*mean(c(1,length(viewr))),wcorrs[counter,viewr],pch=4,col="Blue")
        counter <- counter+1
      }
      i <- i+1
    }
    axis(1,at=1:ncorrs,labels=names,cex=1.2)
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),col=c("Dark green","Blue"),pt.bg=c(alpha("Dark green",0.2),alpha("Blue",0.2)),pch=22)
  }
  
  # Correlations between population sizes after shock
  if(any(plot=="all") | any(plot=="aftercorr")){
    # Same kind of loop as for networks, plotting all combinations
    scorrs <- wcorrs <- matrix(NA,ncorrs,reps)
    names <- character(choose(N,2))
    i <- counter <- 1
    summers <- seq(rec.int*(t-1)+1,(rec.int*t)-1,by=2)
    winters <- summers+1
    plot(1:ncorrs,rep(0,ncorrs),type="n",xaxt="n",xlim=c(ifelse(N>3,1,0.5),ifelse(N>3,ncorrs,ncorrs+0.5)),ylim=c(-1,1),ylab="Correlation",xlab="",main="50 years after shock",bty="L")
    while(i<N){
      for(j in (i+1):N){
        polygon(c(counter-0.3,counter+0.3,counter+0.3,counter-0.3),c(-1.2,-1.2,1.2,1.2),col=alpha("Light grey",alpha=min(1,shocks[i]+shocks[j])),border=NA)
        scorrs[counter,] <- diag(cor(log(dyn[summers,i,]),log(dyn[summers,j,])))
        wcorrs[counter,] <- diag(cor(log(dyn[winters,i,]),log(dyn[winters,j,])))
        names[counter] <- paste0(LETTERS[i],LETTERS[j])
        vioplot(scorrs[counter,],add=T,at=counter-0.1,col=alpha("Dark green",0.2),wex=0.25,border="Dark green",lineCol=NA)
        vioplot(wcorrs[counter,],add=T,at=counter+0.1,col=alpha("Blue",0.2),wex=0.25,border="Blue",lineCol=NA)
        #points(counter-0.05*(1:length(viewr)),scorrs[counter,viewr],pch=2,col="Dark green")
        #points(counter+0.05*(1:length(viewr)),wcorrs[counter,viewr],pch=4,col="Blue")
        counter <- counter+1
      }
      i <- i+1
    }
    axis(1,at=1:ncorrs,labels=names,cex=1.2)
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),col=c("Dark green","Blue"),pt.bg=c(alpha("Dark green",0.2),alpha("Blue",0.2)),pch=22)
  }
  
  ### Average over all recordings of summary stats for last 1000 years.
  if(N>3) t <- t-2 # Might need to hashtag away...
  ts <-  (t-1000%/%rec.int):t # 20 recordings if rec.int=50
  
  if(any(plot=="all") | any(plot=="avgBVs")){
    par(mfrow=c(2,4),mar=c(4,4,2,1),oma=c(0,0,1.5,0))
    
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Breeding values",ylim=range(evol[1:2,,t,],na.rm=T),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),apply(evol[1,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),apply(evol[2,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8)
    }
    abline(h=0,lty=2)
    mtext("Average over recordings from last 1000 years",outer=T,cex=1.4,line=-0.5)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop BV",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="avgBVviol")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values",bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,na.rm=T,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  }
  
  if(any(plot=="all") | any(plot=="avgLs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Liabilities",ylim=c(range(liab[1:2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),apply(liab[1,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),pch=16,cex=1.8) # Means: filled points
      points(jitter(1:N,amount=0.25),apply(liab[2,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8) # SDs: open points
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop L",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="avgLiabviol")){ # Plots violins of liabilities in final year for viewr random replicates
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,6,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Liabilities")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,6,t,r]))){
          vioplot(popstorage[inds,6,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="avgMeans")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population means",ylim=c(range(liab[1:2,,t,])),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.2),apply(evol[1,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.2),apply(liab[1,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="avgSDs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population SDs",ylim=c(0,max(liab[2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.2),apply(evol[2,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.2),apply(liab[2,,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    if(legend) legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="avgDests")){
    plot(0:(N+1),0:(N+1),type="n",ylim=c(0,1),ylab="Frequency in subpopulation",xaxt="n",xlab="Subpopulation",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),apply(liab[3,,ts,r],FUN=mean,MARGIN=1,na.rm=T),pch=18,col=alpha("Black",0.5),cex=1.6) # proportion migrants
      for(j in 1:N){ # First plot all the points of 'going to A', then all 'going to B', etc.
        paltemp <- rep(pal[j],N)
        paltemp[j] <- NA
        points(jitter(1:N,amount=0.25),apply(evol[3+(j-1),,ts,r],FUN=mean,MARGIN=1,na.rm=T),col=alpha(paltemp[1:N],0.7),cex=1.2)
      }
    }
    if(legend) {
      legend("topleft",legend=c("Migrants","Destination alleles"),pch=c(18,1),pt.cex=c(1.6,1.2),bty="n",bg=alpha("white",0.5),col=alpha("Black",0.5),cex=0.9)
      if(N>3) legend("topright",legend=sites[1:N],col=pal[1:N],bg=alpha("White",0.5),cex=1+(1-N)/20,pch=1,title="Destination")
    }
  }
  
  ts <- seq(dim(dyn)[1]-2000,dim(dyn)[1]-1,by=2) # Average pop densities over every year
  if(any(plot=="all") | any(plot=="avgDensity")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Site",ylab="Local population size",xaxt="n",cex.lab=1.2,ylim=c(0,ifelse(grepl("extreme",saved,fixed=T),3,2)*K))
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),apply(dyn[ts,,r],FUN=mean,MARGIN=2,na.rm=T),pch=16,col=alpha(pal[1:N],0.5),cex=1.8) # Breeding 
      points(jitter(1:N,amount=0.25),apply(dyn[ts+1,,r],FUN=mean,MARGIN=2,na.rm=T),col=alpha(pal[1:N],0.5),cex=1.8) # Non-breeding
    }
    abline(h=K,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),pch=c(16,1),pt.cex=1.8,col=alpha("Black",0.5),bty="n")
  }
  
  ### Final time step summary stats
  if(any(plot=="all") | any(plot=="EndBVs")){
    par(mfrow=c(2,4),mar=c(4,4,2,1),oma=c(0,0,1.5,0))

    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Breeding values",ylim=range(evol[1:2,,t,],na.rm=T),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),evol[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    abline(h=0,lty=2)
    mtext("Post shock",outer=T,cex=1.4,line=-0.5)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop BV",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="EndBVviol")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values",bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,na.rm=T,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  }
  
  if(any(plot=="all") | any(plot=="EndBVbar")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab=expression(paste("Breeding values, ",italic(a))),bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      inds <- ((j-1)*K+1):(j*K)
      for(r in 1:reps){
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha("white",0),colMed=pal[j],colMed2=alpha(pal[j],0.3),border=NA,na.rm=T,lineCol=NA,rectCol=alpha(pal[j],0.3),cex=1.5,pchMed=1) #pchMed=21 if filled points
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="EndLs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Liabilities",ylim=c(range(liab[1:2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),liab[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8) # Means: filled points
      points(jitter(1:N,amount=0.25),liab[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8) # SDs: open points
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop L",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="EndLiabviol")){ # Plots violins of liabilities in final year for viewr random replicates
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,6,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Liabilities")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,6,t,r]))){
          vioplot(popstorage[inds,6,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="EndMeans")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population means",ylim=c(range(liab[1:2,,t,])),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),liab[1,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="EndSDs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population SDs",ylim=c(0,max(liab[2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[2,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),liab[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    if(legend) legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="EndDensity")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Site",ylab="Local population size",xaxt="n",cex.lab=1.2,ylim=c(0,ifelse(grepl("extreme",saved,fixed=T),3,2)*K))
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),dyn[(rec.int*t*2)-1,,r],pch=16,col=alpha(pal[1:N],0.5),cex=1.8)
      points(jitter(1:N,amount=0.25),dyn[rec.int*t*2,,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }
    abline(h=K,lty=2)
    
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),pch=c(16,1),pt.cex=1.8,col=alpha("Black",0.5),bty="n")
  }
  
  if(any(plot=="all") | any(plot=="EndDests")){
    plot(0:(N+1),0:(N+1),type="n",ylim=c(0,1),ylab="Frequency in subpopulation",xaxt="n",xlab="Subpopulation",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),liab[3,,t,r],pch=18,col=alpha("Black",0.5),cex=1.6) # proportion migrants
      for(j in 1:N){ # First plot all the points of 'going to A', then all 'going to B', etc.
        paltemp <- rep(pal[j],N)
        paltemp[j] <- NA
        points(jitter(1:N,amount=0.25),evol[3+(j-1),,t,r],col=alpha(paltemp[1:N],0.7),cex=1.2)
      }
    }
    if(legend) legend("topleft",legend=c("Migrants","Destination allele"),pch=c(18,1),pt.cex=c(1.6,1.2),bty="n",bg=alpha("white",0.5),col=alpha("Black",0.5),cex=0.9)
    #legend("topright",legend=sites[1:N],col=pal[1:N],bg=alpha("White",0.5),cex=1+(1-N)/20,pch=1,title="Destination")
  }
  
  ### End-but-pre-shock summary stats
  t <- dim(popstorage)[3]-2
  
  if(any(plot=="all") | any(plot=="preshockBVs")){
    par(mfrow=c(2,4),mar=c(4,4,2,1),oma=c(0,0,1.5,0))
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Breeding values",ylim=range(evol[1:2,,t,],na.rm=T),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),evol[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    abline(h=0,lty=2)
    mtext("Pre shock",outer=T,cex=1.4,line=-0.5)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop BV",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="preshockBVviol")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,na.rm=T,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  }
  
  if(any(plot=="all") | any(plot=="preshockBVviol")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values",bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,na.rm=T,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  }
  
  if(any(plot=="all") | any(plot=="preshockBVbar")){
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,1,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Breeding values",bty="L")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      inds <- ((j-1)*K+1):(j*K)
      for(r in 1:reps){
        if(any(!is.na(popstorage[inds,1,t,r]))){
          vioplot(popstorage[inds,1,t,r],add=T,at=jitter(j,amount=0.25),col=alpha("white",0),border=NA,na.rm=T,lineCol=NA,rectCol=alpha(pal[j],0.5))
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="preshockLs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Liabilities",ylim=c(range(liab[1:2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),liab[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8) # Means: filled points
      points(jitter(1:N,amount=0.25),liab[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8) # SDs: open points
    }    
    abline(h=0,lty=2)
    if(legend) legend("bottomleft",legend=c("Mean","Sd"),bty="n",pch=c(16,1),col=alpha("Black",0.5),title="Within-pop L",pt.cex=1.8)
  }
  
  if(any(plot=="all") | any(plot=="preshockLviol")){ # Plots violins of liabilities in final year for three random replicates
    plot(0:(N+1),1:(N+2),type="n",ylim=range(popstorage[,6,t,],na.rm=T),xaxt="n",cex.lab=1.2,xlab="Subpopulation",ylab="Liabilities")
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(j in 1:N){
      for(r in viewr){
        inds <- ((j-1)*K+1):(j*K)
        if(any(!is.na(popstorage[inds,6,t,r]))){
          vioplot(popstorage[inds,6,t,r],add=T,at=jitter(j,amount=0.25),col=alpha(pal[j],0.1),border=NA,lineCol=NA,rectCol=alpha("black",0.3))
        }
      }
    }
    abline(h=0,lty=2)
  } 
  
  if(any(plot=="all") | any(plot=="preshockMeans")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population means",ylim=c(range(liab[1:2,,t,])),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[1,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),liab[1,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    abline(h=0,lty=2)
    legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }
  
  if(any(plot=="all") | any(plot=="preshockSDs")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Subpopulation",ylab="Within-population SDs",ylim=c(0,max(liab[2,,t,],na.rm=T)),pch=16,xaxt="n",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),evol[2,,t,r],col=alpha(pal[1:N],0.5),pch=16,cex=1.8)
      points(jitter(1:N,amount=0.25),liab[2,,t,r],col=alpha(pal[1:N],0.5),cex=1.8)
    }    
    legend("bottomleft",legend=c("Breeding value","Liability"),pch=c(16,1),pt.cex=1.8,bty="n")
  }

  if(any(plot=="all") | any(plot=="preshockDensity")){
    plot(0:(N+1),0:(N+1),type="n",xlab="Site",ylab="Local population density",xaxt="n",cex.lab=1.2,ylim=c(0,2*K))
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),dyn[(rec.int*t*2)-1,,r],pch=16,col=alpha(pal[1:N],0.5),cex=1.8) # Breeding
      points(jitter(1:N,amount=0.25),dyn[rec.int*t*2,,r],col=alpha(pal[1:N],0.5),cex=1.8) # Non-breeding
    }
    abline(h=K,lty=2)
    if(legend) legend("bottomleft",legend=c("Breeding","Non-breeding"),pch=c(16,1),pt.cex=1.8,col=alpha("Black",0.5),bty="n")
  }
  
  if(any(plot=="all") | any(plot=="preshockDests")){
    plot(0:(N+1),0:(N+1),type="n",ylim=c(0,1),ylab="Frequency in subpopulation",xaxt="n",xlab="Subpopulation",cex.lab=1.2)
    axis(side=1,at=1:N,labels=sites[1:N],cex.lab=1)
    for(r in 1:reps){
      points(jitter(1:N,amount=0.25),liab[3,,t,r],pch=18,col=alpha("Black",0.5),cex=1.6) # proportion migrants
      for(j in 1:N){ # First plot all the points of 'going to A', then all 'going to B', etc.
        paltemp <- rep(pal[j],N)
        paltemp[j] <- NA
        points(jitter(1:N,amount=0.25),evol[3+(j-1),,t,r],col=alpha(paltemp[1:N],0.7),cex=1.2)
      }
    }
    if(legend){
      legend("topleft",legend=c("Migrants","Destination allele"),pch=c(18,1),pt.cex=c(1.6,1.2),bty="n",bg=alpha("white",0.5),col=alpha("Black",0.5),cex=0.9)
      if(N>3) legend("topright",legend=sites[1:N],col=pal[1:N],bg=alpha("White",0.5),cex=1+(1-N)/20,pch=1,title="Destination")
    }
  }

  # Age distribution from up to 3 random replicates before shocks.
  if(any(plot=="all") | any(plot=="agedist")){ 
    par(mfrow=c(length(viewr),1),mar=c(3.5,4.2,0.5,0.3),oma=c(0,0,1.5,0))
    for(r in 1:length(viewr)){
      plot(1:max.age,rep(0.1,max.age),type="n",xlab="",ylab="Frequency",cex.lab=1.2,cex.axis=1.2,ylim=c(0,max(evol[(N+3):(N+2+max.age),,,],na.rm=T)))
      for(j in 1:N){
        lines((1:max.age)+0.1*(j-1),evol[(N+3):(N+2+max.age),j,t,viewr[r]],col=pal[j],type="h",lwd=2)
        abline(v=mean(popstorage[((j-1)*K+1):(j*K),3,t,viewr[r]],na.rm=T),lty=2,col=pal[j])
      }
      text(max.age/2,0.25,paste("Replicate",viewr[r]),cex=1.5)
    }
    mtext("Age class",side=1,line=2)
    mtext("Age distribution, pre shock",side=3,line=-0.3,cex=1.2,outer=T)
    if(legend) legend("topright",title="Home site",legend=sites[1:N],col=pal,pch=rep(15,N))
  }
  
}

allbutpopsize <- c("indBVs", "indLs", "meanBVs", "meanLs","endpopsize","winteruse","agedist","Dest","meancorr","beforecorr","aftercorr",
                   "EndBVs", "EndBVviol", "EndLs", "EndLiabviol", "EndMeans", "EndSDs", "EndDests", "EndDensity","EndBVbar",
                   "avgBVs", "avgBVviol", "avgLs", "avgLiabviol", "avgMeans", "avgSDs", "avgDests", "avgDensity",
                   "preshockBVs", "preshockBVviol", "preshockLs", "preshockLviol", "preshockMeans", "preshockSDs", "preshockDests","preshockDensity","preshockBVbar")
shorter <- c("meanBVs","meanLs","endpopsize","preshockBVviol","preshockLviol","preshockDests","preshockDensity")
#####

#### Functions for individual-level database & analyses ####

### create.alldb(): Function for creating individual-level database, with one row per individual per year. Slow for large files!
## Arguments:
# data     : A file created by sim() - locally or read from file
# small    : Whether to subset to only a few individuals, in order for it to run faster (useful). Default NA.
#           If changing, use a number (e.g. s=5000), then it uses the s individuals with highest BVs.
## Outputs: A data frame with columns: ID, Year, BV, Dest, Age, Home, Winter, L, Mort, ind, max.age, problem
create.alldb <- function(data,small=NA){
  colnames(data[[4]]) <- c("ID","Year","BV","Dest","Age","Home","Winter","L","Mort")
  alldf <- data.frame("ID"=numeric(0),"Year"=numeric(0),"BV"=numeric(0),"Dest"=numeric(0),"Age"=numeric(0),"Home"=numeric(0),"Winter"=numeric(0),"L"=numeric(0),"Mort"=numeric(0),"Rep"=numeric(0))
  for(i in 1:dim(data[[4]])[3]){
    tmp <- as.data.frame(data[[4]][,,i]) %>% mutate("Rep"=i)
    alldf <- bind_rows(alldf,tmp)
  }
  # Make sure each individual gets its unique identity "ind".
  # First rank the individuals by breeding value.
  alldf <- mutate(alldf,ind=rank(BV,ties.method="min",na.last="keep"))
  
  alldf <- mutate(alldf,new=as.numeric(paste0(ind,ID))) %>% filter(new!="NANA") %>% filter(!is.na(new)) %>% group_by(new)
  alldf <- group_by(alldf,new) %>% mutate(max.age=max(Age,na.rm=T), problem=any(tabulate(Age)>1) | sum(Age==0)>1) %>% ungroup() %>% arrange(new) 
  alldf <- alldf %>% select(-ind) %>% rename(ind=new)
  #  if(!is.na(small)){
  #   alldf <- filter(alldf,ind %in% sample(unique(alldf$ind),size=small,replace=F)) # This should take less time! 
  #}
  #uniqueinds <- as.numeric(sort(unique(alldf$ind)))
  
  #  t <- Sys.time()
  # Then split up parents and offspring who have the same breeding values - use unique ID. Warning: this takes time! ca 15 minutes / 1 million observations?
  # for(i in uniqueinds){ # For each new identity
  #  counter <- i # For assigning new IDs
  # uniqueIDs <- unique(alldf[alldf$ind==i,1]) # Which old IDs have this identity?
  # Assumes that no individuals had the same slot and the same breeding value.
  #for(j in uniqueIDs){ ## CAN THIS LOOP BE SPED UP WITH GROUP_BY MULTIPLE CONDITIONS?? How to assign the separate number to each group?
  # alldf[alldf$ind %in% i & alldf$ID %in% j,11] <- counter
  #counter <- counter + 1
  #    }
  #   print(i)
  #}
  
  #print(Sys.time()-t)
  #alldf <- group_by(alldf,ind) %>% mutate(max.age=max(Age,na.rm=T), problem=any(tabulate(Age)>1) | sum(Age==0)>1) %>% ungroup() %>% arrange(ind) 
  return(alldf)
}

### polish.alldb(): Function for fixing problems in a database created by create.alldb. Returns a database with same dimensions as create.alldb().
polish.alldb <- function(newdf){
  t <- Sys.time()
  allinds <- 1:10*length(newdf$ind) # or max(newdf$ind,na.rm=T)? Huge when new IDs merged ind and ID. 
  free <- allinds[-na.omit(unique(newdf$ind))] # Find free slots to fill in when we split our individuals
  free <- c(free,tail(allinds,1)+1)
  cycle <- unique(newdf$ind[newdf$problem==T])
  cycle <- cycle[!is.na(cycle)]
  while(length(cycle)*10>length(free)){
    free <- c(free,tail(free,1)+1:length(cycle))
  }
  
  counter <- 1 
  for(i in cycle){ # For each problematic individual
    start <- counter # This is where the counter was when we started splitting this individual.
    n <- which(newdf$ind==i) # First identify the places this individual occupies
    morts1 <- which(diff(as.matrix(newdf[n,5]))!=1) # Where is the age in the next cell not 1 more than the age in the previous cell?
    morts2 <- which(newdf[newdf$ind %in% i,9]>0) # Registered death events
    morts3 <- sort(union(morts1,morts2))
    morts1 <- which(newdf[n,2]==max(newdf$Year,na.rm=T)) # Individual who reach the end of the simulation
    morts <- sort(union(morts1,morts3))
    newdf[n[1]:n[morts[1]],11] <- free[counter]
    counter <- counter+1
    if(length(morts)>1){
      for(j in 2:length(morts)){
        newdf[(n[morts[j-1]]+1):(n[morts[j]]),11] <- free[counter]
        counter <- counter+1
      }
    }
    #    for(j in free[start:(counter-1)]){
    #     newdf[newdf$ind %in% j,12] <- max(newdf[newdf$ind %in% j,5],na.rm=T) #Now assign new max.ages to all individuals split in this process.
    #  }
    print(i)
  }
  newdf <- group_by(newdf,ind) %>% mutate(max.age=max(Age,na.rm=T), problem=any(tabulate(Age)>1) | sum(Age==0)>1 ) %>% ungroup() %>% arrange(ind)
  newdf <- mutate(newdf,migrated=as.numeric(L>0))
  print(Sys.time()-t)
  return(newdf)
}

## Function linkstrength(): For calculating the strength of the link between two patches.
# x is the square matrix of transition probabilities (columns: from summer location; rows: to winter location)
# s1 and s2 are the subpopulation for which link strength is to be calculated. Must be as.numeric(subpopulation name). 0 < s1,s2 =< N.
linkstrength <- function(x,patch1,patch2){
  if(nrow(x)!=ncol(x)){
    print("x not a square matrix")
    break
  }
  if(any(colSums(x)>1.001) | any(colSums(x)<0.999)){
    print("Transition probabilities do not sum to 1")
    break
  }
  N <- dim(x)[1]
  tmp <- x[patch1,patch2]*x[patch2,patch2] + x[patch2,patch1]*x[patch1,patch1]
  for(i in (1:N)[-c(patch1,patch2)]){
    tmp <- tmp + x[patch1,i]*x[patch2,i]
  }
  
  return(tmp)
}

## Function trans.matrix(): For calculating the transition probabilities between patches, using polished data frame alldf.
# Output: A NxN transition matrix from (columns) each patch to (rows) each patch.
trans.matrix <- function(alldf){
  N <- max(alldf$Home)
  tmp <- matrix(NA,N,N) # Columns: From summer location; Rows: To winter location
  for(j in 1:N){
    tmp[,j] <- tabulate(alldf$Winter[alldf$Home==j])
    tmp[,j] <- tmp[,j]/sum(tmp[,j]) # Standardize based on number of observations from each patch?
  }
  return(tmp)
}

## Function for producing network plot and calculating link strengths (calling linkstrength() function). Returns the vector of link strengths and (if plot=T) a plot where edge color and thickness reflects link strength.
network <- function(x,plot=T,widthscale=10,main=NULL){
  N <- dim(x)[1]
  
  i <- counter <-  1
  links <- numeric(choose(N,2))
  names <- character(choose(N,2))
  if(plot){
    pal <- brewer.pal(8,"Dark2")
    origin <- (N+1)/2
    xs <- 1:N
    ys <- sqrt(((N+1)/2)^2-((1:N)-origin)^2)
    plot(xs,ys,xlim=c(0.5,N+0.5),ylim=c(min(ys)-0.3,max(ys)+0.3),xaxt="n",yaxt="n",main=main,xlab="Subpopulation",ylab="")
    axis(1,at=xs,labels=LETTERS[xs])
  }  
  while(i<N){
    for(j in (i+1):N){
      links[counter] <- linkstrength(x,i,j)
      names[counter] <- paste0(LETTERS[i],LETTERS[j]) #Finish this
      if(plot){
        lines(c(xs[i],xs[j]),c(ys[i],ys[j]),lwd=log(widthscale*links[counter]),col=grey((1-links[counter])^2))
      }
      counter <- counter+1
    }
    i <- i+1
  }
  if(plot){
    points(xs,ys,cex=4,pch=21,col=pal[1:N],bg=pal[1:N])
  }
  names(links) <- names
  return(links)
}


#####

#### SIMULATION ####

## Arguments
#  N         Number of patches (zones)
#  K         Carrying capacity per patch
#  T         End time
#  mu        Per-locus mutation rate. Default 0.01.
#  m1        Breeding season survival. A vector of length N. Default rep(0.9,N), but can use c(0.95,0.9,0.85). But be careful with overall lifespan!
#  m2        Strength of winter density-dependence (higher m2 -> lower survival prob). A vector of length N. Default type is exponential
#            Exponential: Survival=dexp(n_i/K,m2)/m2) Use e.g. c(1.5,1,0.5). But lines cross! see ddfun(). 
#            Linear: Survival=1-(n_i/K)*m2[i]. Use e.g. c(1/2,1/3,1/10)
#  m3        Mortality cost of migration. Default 0.01
#  maxClutch Maximum (female) clutch size. A vector of length N. Default rep(3,N).
#  age.mat   Age at first reproduction. Default 1.
#  sigma_e   Environmental component to autumn migration decision. Default 1. (Standard deviation of normal distr. from which values are chosen)
#  lambda_d  Demographic (population density) component to autumn migration decision. Default 0.1.
#  disp      Dispersal rate (between-patch movement). Default 0.01.
#  grain     Whether stochastic environmental cues for migration decisions are coarse-grained ("coarse", one for each subpop, default) or fine-grained ("fine", one fore each individual)
#  sigma_m1  Stochastic components to breeding season survival. A vector of length N. Default rep(0,N).
#  sigma_m2  Stochastic components to overwinter survival. A vector of length N. Default rep(0,N).
#  rec.int   Time interval for recording full population matrix. Default 50.
#  rec.all.t Number of final time steps over which to record full population matrix every year. Default 100.
#  rec.all.r Number of replicates over which to record full population matrices. Default 5.
#  shock.time When to introduce perturbations (shocks). Can be a vector. Default T-rec.all.t/2
#  shock.size Size of shocks (fraction of population) for each patch. A vector of length N. Default c(0,0.8,0).
#  distance  Whether all patches are equally easy to reach (FALSE, default) or whether farther away ones are costlier (TRUE).
#  sex       Whether reproduction is sexual (T, default) or asexual (F).
#  v0        Initial additive genetic variance
#  a0        Initial mean breeding value
#  popsaved  Whether to initiate simulation with an already saved evolved population. Default FALSE. If TRUE, load the whole simulation output file, using readRDS(), or local file.
#  ...       Additional arguments.
## Outputs:
#  A list of 4 elements
#  [[1]] popstorage: Storage matrix of format:
#  Rows:    Individuals
#  Columns: Individual-level information recorded at regular intervals:
#  1: BV    Breeding value for genetic component of liability. Migration occurs if BV>0
#  2: Dest  Gene for determining destination if migrating. Randomly sampled among all patches except home patch
#  3: Age   Current age. Starts at 0
#  4: Home  Each individual's breeding population
#  5: Location: The zone the individual is spending the winter in
#  6: L     Liability value. Individual migrated if L>0.
#  for each recorded time step (3rd dim) and replicate (4th dim) 
#  [[2]] dyn: Storage matrix of dimensions:
#  1: Time steps (NB! two recordings per year! Odd numbers are summer, even numbers winter)
#  2: Patches
#  3: Replicate simulations.
#  [[3]] pars: list of all function inputs.
#  [[4]] alldata: Storage matrix of format:
#  Rows:    Individuals
#  Rolumns: Individual-level information recorded every year for the last rec.all.t years:
#  1: Individual ID
#  2: Year
#  3: BV gene
#  4: Dest gene
#  5: Age
#  6: Home patch
#  7: Wintering location in year t
#  8: Liability in year t
#  9: Cause of mortality. 0: Still alive. 1: Breeding season mortality. 2: Winter mortality. 3: Migration mortality.
#  for each replicate simulation (3rd dim)
sim <- function(reps=10,N=3,K=1000,T=5000,m1=rep(0.9,N),m2=c(1/2,1/3,1/5),m3=0.01,maxClutch=rep(3,N),age.mat=1,mu=0.01,m_size=0.1,sigma_e=1,lambda_d=0.1,
                disp=0.01,grain="coarse",sigma_m1=rep(0,N),sigma_m2=rep(0,N),rec.int=50,rec.all.t=100,rec.all.r=5,shock.time=T-1-rec.all.t/2,
                shock.size=rep(NA,N),info="None",distance=F,sex=TRUE,v0=1,a0=0,popsaved=FALSE,...){
  t <- Sys.time()
  if(length(shock.size)!=N){
    print(paste("Argument shock.size must be a vector of length N"))
    break
  }
  if(length(m2)!=N){
    print(paste("Argument m2 must be a vector of length N"))
    break
  }
  if(length(sigma_m2)!=N){
    print(paste("Argument sigma_m2 must be a vector of length N"))
    break
  }
  if(length(m1)!=N){
    print(paste("Argument m1 must be a vector of length N"))
    break
  }
  if(length(sigma_m1)!=N){
    print(paste("Argument sigma_m1 must be a vector of length N"))
    break
  }
  if(length(maxClutch)!=N){
    print(paste("Argument maxClutch must be a vector of length N"))
    break
  }

  # Create storage vectors for recording data #
  dyn <- array(0,dim=c(T*2,N,reps)) # Dimensions: 1) Time steps (NB! two recordings per year! Odd numbers are summer, even numbers winter), 2) Patches, 3) Replicate simulations.
  colnames(dyn) <- LETTERS[1:N]
  popstorage <- array(NA,c(N*K,6,T%/%rec.int,reps)) # Dims: 1) Individuals, 2) Pop matrix contents (see below), 3) recording intervals, 4) replicates
  pars <- list(reps,N,K,T,m1,m2,m3,maxClutch,age.mat,mu,m_size,sigma_e,lambda_d,disp,grain,sigma_m1,sigma_m2,rec.int,rec.all.t,rec.all.r,shock.time,shock.size,info,distance,sex,v0,a0)
  names(pars) <- c("reps","N","K","maxTime","Breeding-season survival","Winter mort","Migration cost","maxClutch","age@mat","Mut rate","Mut size","sigma_e","lambda_d","Dispersal rate","grain",
                   "sigma_m1","sigma_m2","rec.int","rec.all.t","rec.all.r","shock.time","shock.size","info","distance","sex","v0","Initial mean BV")
  allstorage <- array(NA,c(N*K*rec.all.t*1.5,9,min(reps,rec.all.r)))

  for(r in 1:reps){
    #### Initialize population ####
    cue <- dens <- numeric(N*K) # Empty vectors for storing information on migration decision
    
    alldata <- matrix(NA,N*K*rec.all.t*1.5,9) 
    ### alldata matrix: 6 columns, 1 observation per row.
    # 1: Individual ID
    # 2: Year
    # 3: BV gene
    # 4: Dest gene
    # 5: Age
    # 6: Home patch
    # 7: Wintering location in year t
    # 8: Liability in year t
    # 9: Cause of mortality. 0: Still alive. 1: Breeding season mortality. 2: Winter mortality. 3: Migration mortality.
    colnames(alldata) <- c("ID","Year","BV","Dest","Age","Home","Winter","L","Mort")
    obs <- 1 # Counter for observations in Alldata
    
    n_off <- integer(1)
    inds <- fill <- rep(NA,K) # Initialize necessary storage vectors
    if(sex){
      mothers <- fathers <- rep(NA,K)
    } else {lucky <- rep(NA,K)}
    
    extinct <- N+1
    
    if(!isFALSE(popsaved)){ # If using a pre-saved population, 
      pop <- popsaved[[1]][,,dim(popsaved[[1]])[3],r]
    } else {
      pop <- matrix(NA,N*K,6)
      colnames(pop) <- c("BV","Dest","Age","Home","Winter","L")
      ### pop matrix: 6 columns, 1 individual per row ###
      # BV:     Breeding value for genetic component of liability. Migration occurs if BV>0
      # Dest:   Gene for determining destination if migrating. Randomly sampled among all patches except home patch
      # Age:    Current age. Starts at 0
      # Home:   Each individual's breeding population
      # Location:The zone the individual is spending the winter in
      # L:      Liability value. Individual migrated if L>0.
      colnames(pop) <- c("BV","Dest","Age","Pop","Zone","L") #
      
      pop[,4] <- rep(1:N,each=K) # Assign zones to each slot in pop matrix
      pop[,5] <- pop[,4] # Overwinter at home unless migrating.
      pop[,3] <- rpois(N*K,3) # Assign starting age distribution. This can be made fancier!
      pop[,1] <- rnorm(N*K,a0,sqrt(v0)) # Assign random breeding values, normally distributed around a0 with var=v0.
      
      # Assign random destination gene values - can't be the patch you breed in.
      for(i in 1:N){
        inds <- ((i-1)*K+1):(i*K) # IDs of individuals in this patch
        pop[inds,2] <- mysample((1:N)[-i],size=K,replace=T) # Assign random destination gene values 
      }
    }
    #####
    
    #### Start simulation ####
    
    for(i in 1:T){ # For every year
      pop[,3] <- pop[,3]+1 # Age increases for surviving individuals. Conveniently, NA+1=NA.
      
      if(i %% rec.int == 0){
        popstorage[,,i%/%rec.int,r] <- pop # Record individual traits
      }
      
      ### BREEDING SEASON ###
      
      for(j in (1:N)[-extinct]){ # For each patch
        fill[] <- NA # Reset parent vectors
        if(sex){
          mothers[] <- fathers[] <- NA
        } else {
          lucky[] <- NA
          removed <- 0}
        n_off <- 0  # Reset offspring counter
        
        inds[] <- ((j-1)*K+1):(j*K) # Slots in this patch
        alive <- which(!is.na(pop[inds,1])) + (j-1)*K # IDs of individuals alive 
        dyn[(i*2)-1,j,r] <- length(alive) # Pre-mortality and pre-breeding summer population census
        
        nmort <- rnorm(1,1-m1[j],sigma_m1[j])
        nmort <- min(1,nmort)
        nmort <- max(0,nmort)
        mort <- mysample(alive,size=round(nmort*length(alive)),replace=F) # Deterministic pre-breeding mortality due to seasonal suitability
        
        if(i>(T-rec.all.t) & rec.all.r >=r  & length(mort)>0){
          alldata[obs:(obs+length(mort)-1),1] <- as.numeric(paste0(r,mort))
          alldata[obs:(obs+length(mort)-1),2] <- i
          alldata[obs:(obs+length(mort)-1),3] <- pop[mort,1]
          alldata[obs:(obs+length(mort)-1),4] <- pop[mort,2]
          alldata[obs:(obs+length(mort)-1),5] <- pop[mort,3]
          alldata[obs:(obs+length(mort)-1),6] <- j
          alldata[obs:(obs+length(mort)-1),9] <- 1
          obs <- obs+length(mort)
        }
        
        pop[mort,] <- NA 
        
        ## Reproduction ##
        n_off <- sum(is.na(pop[inds,1])) # How many empty slots are there for offspring
        if(n_off==K){ # If everyone is dead
          print(paste("Subpopulation",j,"extinct in summer time",i))
          pop[inds,] <- NA # Shouldn't be necessary, but try 15.11 whether this is what causes ages to grow continue growing?
          extinct <- c(j,extinct)
          break
        }
        
        if(max(pop[inds,3],na.rm=T)>=age.mat & n_off>0){ #Check that there are breeders in the population
          n_off <- min((K-n_off)*maxClutch[j],n_off)
          
          if(sex){ # Sexual reproduction
            mothers[1:n_off] <- mysample(which(is.mature(pop[inds,3],a=age.mat))+(j-1)*K,size=n_off,replace=T) # Sample mothers
            fathers[1:n_off] <- mysample(which(is.mature(pop[inds,3],a=age.mat))+(j-1)*K,size=n_off,replace=T) # Sample fathers
            
            # Which slots get new offspring?
            fill[1:n_off] <- head(which(is.na(pop[inds,1])),n_off)+(j-1)*K # Empty slots in this patch
            pop[fill[1:n_off],3] <- 0 # Assign ages to all newborns
            pop[fill[1:n_off],4:5] <- j # Assign home and current patch to all newborns
            pop[fill[1:n_off],1] <- apply(cbind(pop[mothers[1:n_off],1],pop[fathers[1:n_off],1]),FUN=mean,MARGIN=1) # Midpoint of parents' BVs
            # Next, newborn's BV is drawn from a Gaussian around midparent BV with sd given by distance among parents' BVs
            pop[fill[1:n_off],1] <- pop[fill[1:n_off],1]+rnorm(n_off,0,v0)
            pop[fill[1:n_off],1] <- pop[fill[1:n_off],1]+rnorm(n_off,0,sqrt(apply(cbind(pop[mothers[1:n_off],1],pop[fathers[1:n_off],1]),FUN=sd,MARGIN=1)))
            pop[fill[1:n_off],2] <- ifelse(runif(n_off)>0.5,pop[mothers[1:n_off],2],pop[fathers[1:n_off],2]) # Mendelian segregation at the Destination locus
            
            if(any(is.na(pop[fill[1:n_off],1:5]))){
              print("oops")
              print(paste0("N_off=",n_off,", last parent ID=",tail(parents,1)))
            }
          } else { # Asexual reproduction
            if(maxClutch[j]==1){ # Simplify if max clutch size= 1.
              parents <- which(is.mature(pop[inds,3],a=age.mat))+(j-1)*K
              n_off <- min(n_off,parents)
              lucky[1:n_off] <- mysample(parents,size=n_off,replace=F)
            } else {
              lucky[1:n_off] <- mysample(which(is.mature(pop[inds,3],a=age.mat)) + (j-1)*K,size=n_off,replace=T) # IDs of parents who get offspring
              if(any(tabulate(lucky)>maxClutch[j])){ # Make sure no parent gets more than 3 offspring
                problem <- which(tabulate(lucky)>maxClutch[j]) # IDs of parents who got more than 3 offspring
                for(p in problem){ # For each problematic parent
                  while(tabulate(lucky)[p]>maxClutch[j]){ # As long as they have more than 3 offspring
                    lucky[1:(n_off-removed-1)] <- lucky[(1:(n_off-removed))[-head(which(lucky==p),1)]] # Remove one of their offspring
                    lucky[(n_off-removed):K] <- NA
                    removed <- removed+1 
                  }
                }
              }
              n_off <- n_off-removed
            }
            
            # Which slots get new offspring?
            fill[1:n_off] <- head(which(is.na(pop[inds,1])),n_off)+(j-1)*K
            pop[fill[1:n_off],3] <- 0 # Assign ages to all newborns
            pop[fill[1:n_off],4:5] <- j # Assign home and current patch to all newborns
            pop[fill[1:n_off],1] <- pop[lucky[1:n_off],1] # Newborns inherit their parent's genes
            pop[fill[1:n_off],2] <- pop[lucky[1:n_off],2] # Newborns inherit their parent's genes
            
          }
          
          # Mutation - only for Destination if sexual, BV is taken care of by the standard deviation of the random draw.
          if(N>2){
            mutants <- sample(fill[1:n_off],round(n_off*mu)) # A subset of newborns mutate their Destination
            for(m in mutants){
              pop[m,2] <- mysample((1:N)[-c(j,pop[m,2])],1)
            }
          }
          if(!sex){
            mutants <- sample(fill[1:n_off],round(n_off*mu)) # A subset of newborns mutate their BV.
            pop[mutants,1] <- rnorm(length(mutants),pop[mutants,1],m_size)
          }
        } else {
          print(paste("No breeding in patch",j,"in time",i))
        }
      }
      
      if(any(!is.na(pop))==FALSE){
        print(paste("Metapopulation extinct in summer time",i))
        break
      }
      if(length(extinct)==N+1){
        print(paste("Metapopulation extinct in summer time",i))
      }
      ### AUTUMN ###
      
      # Choose stochastic environmental cues
      if(grain=="coarse"){
        for(j in 1:N){
          cue[((j-1)*K+1):(j*K)] <- rnorm(1,0,sigma_e) # Choose a single cue for each patch
        }
      } else if(grain=="fine"){
        cue <- rnorm(N*K,0,sigma_e) # Choose a unique cue for each individual
      }
      
      # Work out density cue for each patch
      for(j in 1:N){
        inds <- ((j-1)*K+1):(j*K) # Slots in this patch
        dens[((j-1)*K+1):(j*K)] <- sum(is.na(pop[inds,1]))/K
      }  
      dens <- dens*lambda_d
      
      # Migrate if env val + pop dens cause liability to exceed threshold
      pop[,6] <- pop[,1] + cue + dens
      pop[which(pop[,6]>0),5] <- pop[which(pop[,6]>0),2] # Individuals who migrate go to patch determined by their Dest gene.
      cue[which(pop[,6]>0)] <- 0 # Individuals who migrate are 'released' from the bad weather in their home patch and get a neutral cue.
      
      if(i %% rec.int == 0){
        popstorage[,5:6,i%/%rec.int,r] <- pop[,5:6] # If recording pop matrix, track liabilities and where individuals overwintered.
      }
      
      # Deterministic mortality depending on seasonal suitability and density in post-migration patch
      for(j in 1:N){
        whos_there <- which(pop[,5]==j)
        localN <- length(whos_there)
        dyn[i*2,j,r] <- localN # Pre-mortality population census
        
        if(m2[j]==0){ # Calculate density-dependence term
          ddterm <- 0
        } else {
          ddterm <- 1-dexp(localN/K,m2[j])/m2[j]
        }
        nmort <- rnorm(1,ddterm,sigma_m2[j])
        nmort <- min(1,nmort)
        nmort <- max(0,nmort)
        if(i==shock.time){
          nmort <- max(c(nmort,shock.size[j]),na.rm=T)
        }

        # Mortality correlated with cue received? Hashtag away the one you're not using
        if(info=="m2"){
          mort <- mysample(whos_there,size=localN*nmort,replace=F,prob=exp(cue[whos_there])) # Cue gives information: Those who didn't migrate get worse weather.
        } else {
          mort <- mysample(whos_there,size=localN*nmort,replace=F) # Cue gives no information
        }
        #mort <- mysample(whos_there,size=min(K,round(localN*localN/K*m2[j])),replace=F)
        
        if(i>(T-rec.all.t) & rec.all.r >=r & length(mort)>0){
          alldata[obs:(obs+length(mort)-1),1] <- as.numeric(paste0(r,mort))
          alldata[obs:(obs+length(mort)-1),2] <- i 
          alldata[obs:(obs+length(mort)-1),3] <- pop[mort,1]
          alldata[obs:(obs+length(mort)-1),4] <- pop[mort,2]
          alldata[obs:(obs+length(mort)-1),5] <- pop[mort,3]
          alldata[obs:(obs+length(mort)-1),6] <- pop[mort,4]
          alldata[obs:(obs+length(mort)-1),7] <- j
          alldata[obs:(obs+length(mort)-1),8] <- pop[mort,6]
          alldata[obs:(obs+length(mort)-1),9] <- 2
          obs <- obs+length(mort)
        }
        
        pop[mort,] <- NA
      }
      ### WINTER (not yet?) ###
      
      if(any(!is.na(pop))==FALSE){
        print(paste("Metapopulation extinct in winter time",i))
        break
      }
      # Choose a stochastic env value for each site
       
      # Individuals still in breeding site migrate if env val + pop dens cause liability to exceed threshold
      
      # Deterministic mortality depending on seasonal suitability in post-migration patch
      
      ### End of winter
      
      # Dispersal: Some individuals from each patch "disperse" to a different patch.
      for(j in 1:N){
        inds <- ((j-1)*K+1):(j*K) # Slots in this patch
        free <- which(is.na(pop[inds,1])) + (j-1)*K # Find free slots
        maxdisp <- length(which(is.na(pop[,1]))[-inds])*disp/(N-1) # Max number of dispersers
        n_imm <- min(length(free),round(maxdisp)) # How many immigrants? Not more than free slots, nor as determined by dispersal parameter.
        if(n_imm>0){
          immigrants <- sample((1:(N*K))[-inds],size=n_imm,replace=F) # Choose random immigrants from other patches
          pop[free[1:n_imm],] <- pop[immigrants,] # Move their slots to this patch
          pop[free[1:n_imm],4] <- j # Immigrants get new home patch
          pop[free[1:n_imm],2] <- mysample((1:N)[-j],size=n_imm,replace=T) # Immigrants get new Dest gene
          pop[immigrants, ] <- NA # Remove their slots in other patches
        }
        if(j %in% extinct){
          extinct <- extinct[-which(extinct==j)]
        }
      }
      # Cost of migration: A fraction of everyone who migrated still alive are killed.
      if(distance){ # Are adjacent patches easier to reach than farther away ones?
        if(any(pop[,6]>0,na.rm=T)){
          mort <- mysample(which(pop[,6]>0),size=round(m3*sum(pop[,6]>0,na.rm=T)),replace=F,prob=abs(pop[which(pop[,6]>0),4]-pop[which(pop[,6]>0),2]))
        }
      } else{ # All patches are equally easy to reach
        mort <- mysample(which(pop[,6]>0),size=round(m3*sum(pop[,6]>0,na.rm=T)),replace=F)
      }
      
      # If recording individual-level data, track those who died from migration
      if(i>(T-rec.all.t) & rec.all.r >=r  & length(mort)>0){
        alldata[obs:(obs+length(mort)-1),1] <- as.numeric(paste0(r,mort))
        alldata[obs:(obs+length(mort)-1),2] <- i 
        alldata[obs:(obs+length(mort)-1),3] <- pop[mort,1]
        alldata[obs:(obs+length(mort)-1),4] <- pop[mort,2]
        alldata[obs:(obs+length(mort)-1),5] <- pop[mort,3]
        alldata[obs:(obs+length(mort)-1),6] <- pop[mort,4]
        alldata[obs:(obs+length(mort)-1),7] <- pop[mort,5]
        alldata[obs:(obs+length(mort)-1),8] <- pop[mort,6]
        alldata[obs:(obs+length(mort)-1),9] <- 3
        obs <- obs+length(mort)
      }
      
      # Kill individuals who died from migration
      pop[mort,] <- NA
      
      # If recording individual-level data, update all information on survived individuals.
      if(i>(T-rec.all.t) & rec.all.r >= r & length(alive)>0){
        alive <- which(!is.na(pop[,1])) 
        alldata[obs:(obs+length(alive)-1),1] <- as.numeric(paste0(r,alive))
        alldata[obs:(obs+length(alive)-1),2] <- i 
        alldata[obs:(obs+length(alive)-1),3] <- pop[alive,1]
        alldata[obs:(obs+length(alive)-1),4] <- pop[alive,2]
        alldata[obs:(obs+length(alive)-1),5] <- pop[alive,3]
        alldata[obs:(obs+length(alive)-1),6] <- pop[alive,4]
        alldata[obs:(obs+length(alive)-1),7] <- pop[alive,5]
        alldata[obs:(obs+length(alive)-1),8] <- pop[alive,6]
        alldata[obs:(obs+length(alive)-1),9] <- 0
        obs <- obs+length(alive)
      }
      # All survived individuals return to their home patch
      pop[,5] <- pop[,4]
      #print(i)
    }

    # End simulation
    #####
    
    if(r <= rec.all.r) {
      allstorage[,,r] <- alldata
    }
    print(r)
  }
  print(Sys.time()-t)
  return(list(popstorage,dyn,pars,allstorage))
}
#####

#### Run and save simulations####
tmp <- sim(T=10000,N=3,reps=20,lambda_d=0, sigma_e=0)
saveRDS(tmp,file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only winter gradient.R")
statstmp <- popstats(saved="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only winter gradient.R")
saveRDS(statstmp,file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only winter gradient stats.R") 

tmp <- sim(T=10000,N=3,reps=20,m1=c(0.95,0.9,0.85),lambda_d=0, sigma_e=0)
saveRDS(tmp,file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Both seasonality gradients.R")
statstmp <- popstats(saved="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Both seasonality gradients.R")
saveRDS(statstmp,file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Both seasonality gradients stats.R")

tmp <- sim(T=10000,N=3,reps=20,m1=c(0.95,0.9,0.85),m2=c(1/3,1/3,1/3),lambda_d=0, sigma_e=0)
saveRDS(tmp,file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only summer gradient.R")
statstmp <- popstats(local=tmp)
saveRDS(statstmp,file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only summer gradient stats.R")

tmp <- sim(T=10000,N=3,reps=20,m2=c(1/3,1/3,1/3),lambda_d=0, sigma_e=0)
saveRDS(tmp,file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\No seasonality gradient.R")
statstmp <- popstats(local=tmp)
saveRDS(statstmp,file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\No seasonality gradient stats.R")
#####

#### Create stats file from saved file ####
statstmp <- popstats(saved="M:\\Documents\\Migration\\sexresults\\sigma_m2=0.4\\No seasonality gradient.R")
saveRDS(statstmp,file="M:\\Documents\\Migration\\sexresults\\sigma_m2=0.4\\No seasonality gradient newstats.R")
#####

#### Create plots from local environment files (example) ####
test <- sim()
statstmp <- popstats(test)
popplots(local=statstmp, mainfile=test)
#####

#### Create plots from saved file ####
#set.seed(50) # If needing to retrieve the same 5 random violins each time.
pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\Figure S1 newlayout.pdf",width=7,height=9.5)
par(mfrow=c(3,2),mar=c(4.2,4.2,1.6,1),oma=c(0,0.3,1.7,0)) # Can be exported in 8x7.5
popplots(saved="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only winter gradient stats.R",plot="EndBVbar")
mtext("No non-genetic components to liability",line=0.3,cex=1)
mtext("(a)",2,line=2.5,las=1,padj=-8.5,cex=1)
mtext("   Alternative scenario",line=0,cex=1.25,outer=T)

popplots(saved="M:\\Documents\\Migration\\sexresults\\Fine-grain\\Only winter gradient stats.R",plot="EndBVbar")
mtext("(b)",2,line=2.5,las=1,padj=-8.5,cex=1)
mtext("Fine-grained environmental variation",line=0.5,cex=1)

popplots(saved="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only winter gradient stats.R",plot="avgDests")
mtext("(c)",2,line=2.5,las=1,padj=-8.5,cex=1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Fine-grain\\Only winter gradient stats.R",plot="avgDests",legend=F)
mtext("(d)",2,line=2.5,las=1,padj=-8.5,cex=1)

popplots(saved="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only winter gradient stats.R",plot="avgDensity")
mtext("(e)",2,line=2.5,las=1,padj=-8.5,cex=1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Fine-grain\\Only winter gradient stats.R",plot="avgDensity",legend=F)
mtext("(f)",2,line=2.5,las=1,padj=-8.5,cex=1)
dev.off()
#####

#### Figure showing differences in among-year variation in proportion migrant ####
par(mfrow=c(2,3),mar=c(4.2,4.2,2.8,1)) # Can be exported in 8x7.5
popplots(saved="M:\\Documents\\Migration\\sexresults\\Only winter gradient stats.R",plot="EndDests")
mtext("(a)",2,line=2.5,las=1,padj=-10.5,cex=1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only winter gradient stats.R",plot="EndDests")
mtext("No non-genetic components to liability",line=0.5,cex=1.2)
mtext("(b)",2,line=2.5,las=1,padj=-10.5,cex=1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Fine-grain\\Only winter gradient stats.R",plot="EndDests")
mtext("(c)",2,line=2.5,las=1,padj=-10.5,cex=1)
#####

#### Plotting proportion migrants at end time for different sigma_mw ####
pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\New fig 4 avgs.pdf",width=9,height=5.5)
par(mfrow=c(1,2),mar=c(2,4,2,0.5),oma=c(2,0,1.5,0))
plot(0.5:4.5,rep(0,5),type="n",xlab="",ylab="Frequency of migrants",xaxt="n",ylim=c(0,1),main="Density-dependent",cex.lab=1.2,cex.main=1.1)
axis(1,at=1:4,labels=c(0.1,0.2,0.3,0.4))
mtext(expression(paste("     Stochasticity in non-breeding season survival, ",italic(sigma)[mn])),side=1,outer=T,line=0.7,cex=1.2)
pal <- viridis(3) # or brewer.pal(8,"Dark2")
mtext("     Non-breeding season survival",side=3,outer=T,line=0.2,cex=1.3)
mtext("(a)",side=2,line=2.5,las=1,padj=-16)
statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.1\\No seasonality gradient stats.R")
dat <- matrix(NA,dim(statstmp[[2]])[2],dim(statstmp[[2]])[4]) # Rows: Patches. Columns: Reps.
for(i in 1:nrow(dat)){  # Extract site-specific means across time for each rep.
  dat[i,] <- apply(statstmp[[2]][3,i,(dim(statstmp[[2]])[3]-20):dim(statstmp[[2]])[3],],FUN=mean,MARGIN=2,na.rm=T)
}
vioplot(dat,at=c(0.8,1,1.2),col=alpha(pal[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))
points(c(0.8,1,1.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4) # Improve?

statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.2\\No seasonality gradient stats.R")
dat <- matrix(NA,dim(statstmp[[2]])[2],dim(statstmp[[2]])[4]) # Rows: Patches. Columns: Reps.
for(i in 1:nrow(dat)){  # Extract site-specific means across time for each rep.
  dat[i,] <- apply(statstmp[[2]][3,i,(dim(statstmp[[2]])[3]-20):dim(statstmp[[2]])[3],],FUN=mean,MARGIN=2,na.rm=T)
}
points(c(1.8,2,2.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)
vioplot(dat,at=c(1.8,2,2.2),col=alpha(pal[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))

statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.3\\No seasonality gradient stats.R")
dat <- matrix(NA,dim(statstmp[[2]])[2],dim(statstmp[[2]])[4]) # Rows: Patches. Columns: Reps.
for(i in 1:nrow(dat)){  # Extract site-specific means across time for each rep.
  dat[i,] <- apply(statstmp[[2]][3,i,(dim(statstmp[[2]])[3]-20):dim(statstmp[[2]])[3],],FUN=mean,MARGIN=2,na.rm=T)
}
points(c(2.8,3,3.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)
vioplot(dat,at=c(2.8,3,3.2),col=alpha(pal[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))

statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.4\\No seasonality gradient newstats.R") 
dat <- matrix(NA,dim(statstmp[[2]])[2],dim(statstmp[[2]])[4]) # Rows: Sites Columns: Reps.
for(i in 1:nrow(dat)){  # Extract site-specific means across time. Loop over sites.
  for(r in 1:ncol(dat)){ # Loop over reps
    yrs <- (statstmp[[4]][1,i,r]%/%50-20):(statstmp[[4]][1,i,r]%/%50)
    yrs <- yrs[yrs>0]
    dat[i,r] <- mean(statstmp[[2]][3,i,yrs,r])
  }
}
points(c(3.8,4,4.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)
vioplot(dat,at=c(3.8,4,4.2),col=alpha(pal[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))
legend("topleft",legend=c("North","Middle","South"),col=pal[1:3],title="Site",bty="n",pch=1)

plot(0.5:4.5,rep(0,5),type="n",xlab="",ylab="Frequency of migrants",xaxt="n",ylim=c(0,1),main="Density-independent",cex.lab=1.2,cex.main=1.1)
axis(1,at=1:4,labels=c(0.1,0.2,0.3,0.4))
mtext("(b)",side=2,line=2.5,las=1,padj=-16,cex=1)
statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.1\\No seasonality gradient, no winter DD stats.R")
dat <- matrix(NA,dim(statstmp[[2]])[2],dim(statstmp[[2]])[4]) # Rows: Patches. Columns: Reps.
for(i in 1:nrow(dat)){  # Extract site-specific means across time for each rep.
  dat[i,] <- apply(statstmp[[2]][3,i,(dim(statstmp[[2]])[3]-20):dim(statstmp[[2]])[3],],FUN=mean,MARGIN=2,na.rm=T)
}
vioplot(dat,at=c(0.8,1,1.2),col=alpha(pal[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))
points(c(0.8,1,1.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)

statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.2\\No seasonality gradient, no winter DD stats.R")
dat <- matrix(NA,dim(statstmp[[2]])[2],dim(statstmp[[2]])[4]) # Rows: Patches. Columns: Reps.
for(i in 1:nrow(dat)){  # Extract site-specific means across time for each rep.
  dat[i,] <- apply(statstmp[[2]][3,i,(dim(statstmp[[2]])[3]-20):dim(statstmp[[2]])[3],],FUN=mean,MARGIN=2,na.rm=T)
}
vioplot(dat,at=c(1.8,2,2.2),col=alpha(pal[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))
points(c(1.8,2,2.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)

statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.3\\No seasonality gradient, no winter DD stats.R")
dat <- matrix(NA,dim(statstmp[[2]])[2],dim(statstmp[[2]])[4]) # Rows: Patches. Columns: Reps.
for(i in 1:nrow(dat)){  # Extract site-specific means across time for each rep.
  dat[i,] <- apply(statstmp[[2]][3,i,(dim(statstmp[[2]])[3]-20):dim(statstmp[[2]])[3],],FUN=mean,MARGIN=2,na.rm=T)
}
points(c(2.8,3,3.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)
vioplot(dat,at=c(2.8,3,3.2),col=alpha(pal[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))

statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.4\\No seasonality gradient, no winter DD stats.R")
dat <- matrix(NA,dim(statstmp[[2]])[2],dim(statstmp[[2]])[4]) # Rows: Patches. Columns: Reps.
for(i in 1:nrow(dat)){  # Extract site-specific means across time for each rep.
  dat[i,] <- apply(statstmp[[2]][3,i,(dim(statstmp[[2]])[3]-20):dim(statstmp[[2]])[3],],FUN=mean,MARGIN=2,na.rm=T)
}
points(c(3.8,4,4.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)
vioplot(dat,at=c(3.8,4,4.2),col=alpha(pal[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))
dev.off()
#####

#### c=2 version of Fig. 4 ####
plot(0.5:3.5,rep(0,4),type="n",xlab=expression(italic(sigma)[mn]),ylab="Proportion migrants",xaxt="n",ylim=c(0,1),main="A) With winter density dependence",cex.lab=1.2)
axis(1,at=1:3,labels=c(0.1,0.2,0.3))
pal <- brewer.pal(8,"Dark2")
statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.1\\c=2, No seasonality gradient stats.R")
vioplot(statstmp[[2]][3,,1000,],at=c(0.8,1,1.2),col=alpha(pal_v[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))
points(c(0.8,1,1.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)
statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.2\\c=2, No seasonality gradient stats.R")
points(c(1.8,2,2.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)
vioplot(statstmp[[2]][3,,1000,],at=c(1.8,2,2.2),col=alpha(pal_v[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))
statstmp <- readRDS("M:\\Documents\\Migration\\sexresults\\sigma_m2=0.3\\c=2, No seasonality gradient stats.R")
points(c(2.8,3,3.2),apply(is.na(statstmp[[2]][3,,1000,]),FUN=mean,MARGIN=1),pch=4)
vioplot(statstmp[[2]][3,,1000,],at=c(2.8,3,3.2),col=alpha(pal_v[1:3],0.1),add=T,lineCol=NA,use.cols=F,wex=0.7,border=alpha(pal[1:3],0.3),rectCol=alpha("black",0.3))
legend("topleft",legend=c("North","Middle","South"),col=pal_v[1:3],title="Site",bty="n",pch=1)
#####

#### Fig showing distribution of winter mortalities ####
pdf(file="M:\\Documents\\Migration\\sexresults\\Fig. 2 density dependence.pdf",width=6,height=5)
par(mfrow=c(1,3),mar=c(4,4,1.5,0.5))

tmp <- readRDS("M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\No gradient, shock_B=0.8.R")
dyn <- tmp[[2]]
stars <- dyn[18000,,]/1000
ddfun(c(1/3,1/3,1/3),title="a)",stars=stars) # Maybe use a veery slight gradient instead?

tmp <- readRDS("M:\\Documents\\Migration\\sexresults\\Only winter gradient.R")
dyn <- tmp[[2]]
stars <- dyn[18000,,]/1000
ddfun(c(1/2,1/3,1/5),title="b)",stars=stars) # Moderate gradient

tmp <- readRDS("M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Extreme gradient, shock_B=0.8.R")
dyn <- tmp[[2]]
stars <- dyn[18000,,]/1000
ddfun(c(5,1,1/10),title="c)",stars=stars,xmax=3) # Extreme gradient
dev.off()
#####

#### Fig S3 - Effect of distance ####
pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\Fig S3 - effect of distance.pdf",width=8,height=8)
par(mfrow=c(2,2),mar=c(4,4,2,0.5),oma=c(0.3,3.2,2.2,0))
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient\\Moderate gradient, shock_C=0.8 stats.R",plot=c("avgDests"),legend=F)
mtext("      Distance travelled affects costs of migration",side=3,line=0.4,outer=T,cex=1.3)
mtext("   Non-breeding season suitability gradient",side=2,line=1.8,outer=T,cex=1.3)
mtext("(a)",2,line=2.5,las=1,padj=-9.5,cex=1.2)
mtext("No",3,line=0.6,cex=1.2)
mtext("Weak",2,line=4.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient w distance\\Moderate gradient, shock_C=0.8 stats.R",plot=c("avgDests"))
mtext("(b)",2,line=2.5,las=1,padj=-9.5,cex=1.2)
mtext("Yes",3,line=0.6,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient\\Extreme gradient, shock_C=0.8 stats.R",plot=c("avgDests"),legend=F)
mtext("(c)",2,line=2.5,las=1,padj=-9.5,cex=1.2)
mtext("Strong",2,line=4.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient w distance\\Extreme gradient, shock_C=0.8 stats.R",plot=c("avgDests"),legend=F)
mtext("(d)",2,line=2.5,las=1,padj=-9.5,cex=1.2)
dev.off()
#####

#### Popdyn plots shocks ####
pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\New Fig 5 morereps.pdf",width=8,height=9)
par(mfrow=c(3,2),mar=c(4,4,0.5,0.5),oma=c(0.3,2.3,3.7,0))
popplots(saved="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_A=0.8 stats.R",plot="endpopsize",legend=F)
mtext("(a)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
mtext("Weak",line=0.3)
mtext("     Geographical gradient in non-breeding season suitability",outer=T,line=2.4)
mtext("North site perturbed",side=2,line=4.5)
popplots(saved="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Extreme gradient, shock_A=0.8 stats.R",plot="endpopsize",legend=T)
mtext("Strong",line=0.3)
mtext("(b)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
popplots(saved="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_B=0.8 stats.R",plot="endpopsize",legend=F)
mtext("Middle site perturbed",side=2,line=4.5)
mtext("(c)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
popplots(saved="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Extreme gradient, shock_B=0.8 stats.R",plot="endpopsize",legend=F)
mtext("(d)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
popplots(saved="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_C=0.8 stats.R",plot="endpopsize",legend=F)
mtext("South site perturbed",side=2,line=4.5)
mtext("(e)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
popplots(saved="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Extreme gradient, shock_C=0.8 stats.R",plot="endpopsize",legend=F)
mtext("(f)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
dev.off()

pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\Fig S5 CIs viridis.pdf",width=8,height=9)
par(mfrow=c(3,2),mar=c(4,4,0.5,0.5),oma=c(0.3,2.3,3.7,0))
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient w distance\\Moderate gradient, shock_A=0.8 stats.R",plot="endpopsize",legend=F)
mtext("(a)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
mtext("Weak",line=0.3)
mtext("     Geographical gradient in non-breeding season suitability",outer=T,line=2.4)
mtext("North site perturbed",side=2,line=4.5)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient w distance\\Extreme gradient, shock_A=0.8 stats.R",plot="endpopsize",legend=T)
mtext("Strong",line=0.3)
mtext("(b)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient w distance\\Moderate gradient, shock_C=0.8 stats.R",plot="endpopsize",legend=F)
mtext("Middle site perturbed",side=2,line=4.5)
mtext("(c)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient w distance\\Extreme gradient, shock_C=0.8 stats.R",plot="endpopsize",legend=F)
mtext("(d)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient w distance\\Moderate gradient, shock_E=0.8 stats.R",plot="endpopsize",legend=F)
mtext("South site perturbed",side=2,line=4.5)
mtext("(e)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient w distance\\Extreme gradient, shock_E=0.8 stats.R",plot="endpopsize",legend=F)
mtext("(f)",2,line=2.3,las=1,padj=-7.9,cex=1.1)
dev.off()
#####

#### Plots of all patches individually perturbed ####
pdf(file=paste0("M:\\Documents\\Migration\\Shock exploration n=3\\Winter gradient\\",as.character(Sys.Date()),
                " Extreme gradient ABC shocks.pdf"))
par(mfrow=c(3,1),mar=c(4,4,2,1))
popplots(saved="M:\\Documents\\Migration\\Shock exploration n=3\\Winter gradient\\2023-01-15 Extreme gradient, shock_A=0.8 stats.R",plot="endpopsize")
mtext("North site perturbed",cex=0.9)
popplots(saved="M:\\Documents\\Migration\\Shock exploration n=3\\Winter gradient\\2023-01-15 Extreme gradient, shock_B=0.8 stats.R",plot="endpopsize")
mtext("Middle site perturbed",cex=0.9)
popplots(saved="M:\\Documents\\Migration\\Shock exploration n=3\\Winter gradient\\2023-01-15 Extreme gradient, shock_C=0.8 stats.R",plot="endpopsize")
mtext("South site perturbed",cex=0.9)
dev.off()
#####

#### c=2 popdyn plots ####
pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\Fig S6 CIs.pdf",width=8,height=9)
par(mfrow=c(3,2),mar=c(4,4,0.5,0.5),oma=c(0,2,4,0))
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Moderate gradient, c=2, shock_A=0.8 stats.R",plot="endpopsize",legend=F)
mtext("Weak",line=0.3)
mtext("       Geographical gradient in non-breeding season suitability",outer=T,line=2)
mtext("North site perturbed",side=2,line=4.5)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Extreme gradient, c=2, shock_A=0.8 stats.R",plot="endpopsize",legend=T)
mtext("Strong",line=0.3)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Moderate gradient, c=2, shock_B=0.8 stats.R",plot="endpopsize",legend=F)
mtext("Middle site perturbed",side=2,line=4.5)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Extreme gradient, c=2, shock_B=0.8 stats.R",plot="endpopsize",legend=F)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Moderate gradient, c=2, shock_C=0.8 stats.R",plot="endpopsize",legend=F)
mtext("South site perturbed",side=2,line=4.5)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Extreme gradient, c=2, shock_C=0.8 stats.R",plot="endpopsize",legend=F)
dev.off()
#####

#### fine-grain popdyn plots ####
pdf(file="M:\\Documents\\Migration\\sexresults\\Alt Fig 5 fine-grain N=5.pdf",width=8,height=9)
par(mfrow=c(3,2),mar=c(4,4,0.5,0.5),oma=c(0,2,4,0))
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient\\Extreme gradient, fine-grain, shock_A=0.8 stats.R",plot="endpopsize",legend=T)
mtext("Strong",line=0.3)
mtext("       Geographical gradient in winter mortality",outer=T,line=2)
mtext("North site perturbed",side=2,line=4.5)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient\\Moderate gradient, fine-grain, shock_A=0.8 stats.R",plot="endpopsize",legend=F)
mtext("Weak",line=0.3)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient\\Extreme gradient, fine-grain, shock_C=0.8 stats.R",plot="endpopsize",legend=F)
mtext("Middle site perturbed",side=2,line=4.5)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient\\Moderate gradient, fine-grain, shock_C=0.8 stats.R",plot="endpopsize",legend=F)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient\\Extreme gradient, fine-grain, shock_E=0.8 stats.R",plot="endpopsize",legend=F)
mtext("South site perturbed",side=2,line=4.5)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Shock exploration n=5\\Winter gradient\\Moderate gradient, fine-grain, shock_E=0.8 stats.R",plot="endpopsize",legend=F)
dev.off()
#####

#### New 6-panel Fig. 2 ####
pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\New fig 2 avgdests.pdf",width=9,height=7)
par(mfrow=c(2,3),mar=c(4,4.2,0.5,0.7),oma=c(0.2,0.1,3.8,0))
set.seed(1)
popplots(saved="M:\\Documents\\Migration\\sexresults\\No seasonality gradient stats.R",plot="avgDests")
mtext("(a)",side=2,line=2.3,las=1,padj=-8.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Only winter gradient moderate stats.R",plot="avgDests",legend=F)
mtext("(b)",side=2,line=2.3,las=1,padj=-8.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Only winter gradient extreme stats.R",plot="avgDests",legend=F)
mtext("(c)",side=2,line=2.3,las=1,padj=-8.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\No seasonality gradient stats.R",plot="EndBVbar")
mtext("(d)",side=2,line=2.3,las=1,padj=-8.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Only winter gradient moderate stats.R",plot="EndBVbar")
mtext("(e)",side=2,line=2.3,las=1,padj=-8.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Only winter gradient extreme stats.R",plot="EndBVbar")
mtext("(f)",side=2,line=2.3,las=1,padj=-8.5,cex=1.2)
mtext("   Gradient in non-breeding season suitability",side=3,line=2,cex=1.5,outer=T)
mtext("       None                                  Weak                                  Strong",side=3,line=-0.2,cex=1.4,outer=T)
dev.off()
#####

#### New 6-panel Fig. 3 ####
pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\New fig 3 bigavgs.pdf",width=9,height=7)
par(mfrow=c(2,3),mar=c(4,4,0.5,0.5),oma=c(0.1,2,3.9,0)) 
set.seed(1)

tmp <- readRDS("M:\\Documents\\Migration\\sexresults\\No seasonality gradient.R")
dyn <- tmp[[2]]
ts <- seq(dim(dyn)[[1]]-1000,dim(dyn)[[1]],by=2) # Plot average over last 50 years
stars <- matrix(NA,dim(dyn)[2],dim(dyn)[3])
for(i in 1:dim(dyn)[2]){ # For each patch
  stars[i,] <- apply(dyn[ts,i,],FUN=mean,MARGIN=2)/1000 
}
ddfun(c(1/3,1/3,1/3),stars=stars) # Maybe use a veery slight gradient instead?
mtext("(a)",side=2,line=2.5,las=1,padj=-8.5,cex=1.2)

tmp <- readRDS("M:\\Documents\\Migration\\sexresults\\Only winter gradient moderate.R")
dyn <- tmp[[2]]
stars <- matrix(NA,dim(dyn)[2],dim(dyn)[3])
for(i in 1:dim(dyn)[2]){ # For each patch
  stars[i,] <- apply(dyn[ts,i,],FUN=mean,MARGIN=2)/1000 
}
ddfun(c(1/2,1/3,1/5),stars=stars) # Moderate gradient
mtext("(b)",side=2,line=2.5,las=1,padj=-8.5,cex=1.2)

tmp <- readRDS("M:\\Documents\\Migration\\sexresults\\Only winter gradient extreme.R")
dyn <- tmp[[2]]
stars <- matrix(NA,dim(dyn)[2],dim(dyn)[3])
for(i in 1:dim(dyn)[2]){ # For each patch
  stars[i,] <- apply(dyn[ts,i,],FUN=mean,MARGIN=2)/1000 
}
ddfun(c(5,1,1/10),stars=stars,xmax=3) # Extreme gradient
mtext("(c)",side=2,line=2.5,las=1,padj=-8.5,cex=1.2)

popplots(saved="M:\\Documents\\Migration\\sexresults\\No seasonality gradient stats.R",plot="avgDensity") # Average over last 50 years
mtext("(d)",side=2,line=2.5,las=1,padj=-8.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Only winter gradient moderate stats.R",plot="avgDensity",legend=F)
mtext("(e)",side=2,line=2.5,las=1,padj=-8.5,cex=1.2)
popplots(saved="M:\\Documents\\Migration\\sexresults\\Only winter gradient extreme stats.R",plot="avgDensity",legend=F)
mtext("(f)",side=2,line=2.5,las=1,padj=-8.5,cex=1.2)
mtext("     Gradient in non-breeding season suitability",side=3,line=2,cex=1.5,outer=T)
mtext("       None                                  Weak                                  Strong",side=3,line=-0.2,cex=1.4,outer=T)
dev.off()
#####

#### Workflow for large simulations:
# 1: Run the simulation to element tmp and save tmp as file with name "x.R".
# 2: Run function popstats(saved="x") and create file statstmp
# 3: Save statstmp as file with name "x stats.R"
# 4: Run function popplots(saved="x stats.R") and use pdf() save to PDF with name "x plots.pdf"
# 5: Run function create.alldb(readRDS(file="x")) to temporary element alldb
# 6: Run polish.alldb(alldb) to temporary element alldb
# 7: Save alldb as file with name "x alldb.R"
 
#### Creating alldbs ####
alldb <- create.alldb(readRDS(file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_A=0.R"))
alldb <- polish.alldb(alldb)
saveRDS(alldb,file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_A=0.8 alldb.R")

alldb <- create.alldb(readRDS(file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_B=0.8.R"))
alldb <- polish.alldb(alldb)
saveRDS(alldb,file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_B=0.8 alldb.R")

alldb <- create.alldb(readRDS(file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_C=0.8.R"))
alldb <- polish.alldb(alldb)
saveRDS(alldb,file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, shock_C=0.8 alldb.R")
saveRDS(alldb,file="M:\\Documents\\Migration\\sexresults\\Fine-grain\\Only winter gradient alldb.R")
#####

#### Individual-level data wrangling and analysis

####Evol dynamics####
create.evoldyn <- function(alldb,years=9926:10000){
  alldb <- readRDS(file=alldb)
  N <- max(alldb$Home,na.rm=T)
  R <- max(alldb$Rep,na.rm=T)
  bvyears <- migryears <- array(NA,dim=c(length(years),N,R)) #Years,patches,reps
  dests <- array(0,dim=c(length(years),N,R,N)) #Years,patches,reps,dests
  alldb <- filter(alldb,Year %in% years)
  evoldyn <- alldb %>% group_by(Home,Year,Rep) %>% summarise(meanbv=mean(BV),sdbv=sd(BV),PropMigr=mean(migrated,na.rm=T))
  t <- Sys.time()
  for(y in 1:length(years)){
    for(n in 1:N){
      for(r in 1:R){
        bvyears[y,n,r] <- filter(evoldyn,Year==years[y] & Home==n & Rep==r)$meanbv
        migryears[y,n,r] <- filter(evoldyn,Year==years[y] & Home==n & Rep==r)$PropMigr
        tmp <- tabulate(filter(alldb,Year==years[y] & Home==n & Rep==r)$Dest)
        if(n==N){
          tmp <- c(tmp,0)
        }
        dests[y,n,r,] <- tmp
        dests[y,n,r,] <- dests[y,n,r,]/sum(dests[y,n,r,])
      }
    }
    print(years[y])
  }
  return(list(evoldyn,bvyears,migryears,dests,years))
}

tmp <- create.evoldyn("C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, c=2, shock_A=0.8 alldb.R")
saveRDS(tmp,file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, c=2, shock_A=0.8 evoldyn.R")
tmp <- create.evoldyn("C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, c=2, shock_B=0.8 alldb.R")
saveRDS(tmp,file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Moderate gradient, c=2, shock_B=0.8 evoldyn.R")
tmp <- create.evoldyn("C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Extreme gradient, c=2, shock_C=0.8 alldb.R")
saveRDS(tmp,file="C:\\Users\\thomhaa\\Local\\n=3 Winter gradient\\Extreme gradient, c=2, shock_C=0.8 evoldyn.R")
#####

#### Fig S4: Detailed dynamics following ECE in C. ####
tmp <- readRDS(file="M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Moderate gradient, shock_C=0.8 evoldyn.R")
evoldyn <- tmp[[1]]
bvyears <- tmp[[2]]
migryears <- tmp[[3]]
dests <- tmp[[4]]
years <- tmp[[5]]
N <- dim(bvyears)[2]
pdf("M:\\Documents\\Migration\\Publication figs Ecology\\Fig 6 CIs.pdf",width=7.5,height=8)
par(mfrow=c(2,2),mar=c(4,4,2,1),oma=c(0,0.4,0,0))
sites <- c("North","Middle","South")
pal <- viridis(N)

plot(1:length(years),rep(0,length(years)),ylim=range(evoldyn$meanbv,na.rm=T),xaxt="n",ylab=expression(paste("Breeding values, ",italic(a))),xlab="Year",type="l",lty=2,lwd=2,main="Metapopulation")
for(n in 1:dim(bvyears)[2]){
  #for(r in 1:dim(bvyears)[3]){
  #  lines(1:dim(bvyears)[1],bvyears[,n,r],col=alpha(pal[n],0.2))
  #}
  lowersum <- apply(bvyears[,n,],FUN=mean,MARGIN=1,na.rm=T) - apply(bvyears[,n,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(bvyears)[3]) # Summer
  uppersum <- apply(bvyears[,n,],FUN=mean,MARGIN=1,na.rm=T) + apply(bvyears[,n,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(bvyears)[3]) # Summer
  polygon(c(1:75,75:1),c(lowersum,rev(uppersum)),col=alpha(pal[n],alpha=0.1),border=NA)
  lines(1:dim(bvyears)[1],apply(bvyears[,n,],FUN=mean,MARGIN=1),lwd=2,col=pal[n]) # Average across reps
}
abline(v=which(years==9950),lty=3,lwd=2)
axis(1,at=c(1,(1:4)*25),labels=years[c(1,(1:4)*25)])
#legend("topleft",legend=sites[1:dim(bvyears)[2]],col=pal[1:dim(bvyears)[2]],pch=15,title="Site",bg=alpha("white",alpha=0.7),box.lty=0)
mtext("(a)",2,line=2.3,las=1,padj=-10.9,cex=1.1)

for(i in 1:dim(bvyears)[2]){ # For each home patch
  plot(1:length(years),rep(0,length(years)),ylim=c(0,1),xaxt="n",ylab="Frequency",xlab="Year",type="n",main=paste(sites[i],"subpopulation")) # 1 plot per subpopulation for this one?
  for(n in 1:dim(bvyears)[2]){ # For each going-to patch
    #for(r in 1:dim(bvyears)[3]){ # For each rep
    #  lines(1:dim(bvyears)[1],dests[,i,r,n],col=alpha(pal[n],0.3))
    #}
    lower <- apply(dests[,i,,n],FUN=mean,MARGIN=1,na.rm=T) - apply(dests[,i,,n],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(dests)[3]) # Summer
    upper <- apply(dests[,i,,n],FUN=mean,MARGIN=1,na.rm=T) + apply(dests[,i,,n],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(dests)[3]) # Summer
    lower[which(lower<0)] <- 0
    upper[which(upper>1)] <- 1
    polygon(c(1:dim(dests)[1],dim(dests)[1]:1),c(lower,rev(upper)),col=alpha(pal[n],alpha=0.1),border=NA)
    lines(1:dim(dests)[1],apply(dests[,i,,n],FUN=mean,MARGIN=1),lwd=2,col=pal[n]) # Average across reps
  }
  lines(1:dim(bvyears)[1],rep(0,dim(bvyears)[1]),lwd=2,col="White")
  #for(r in 1:dim(bvyears)[3]){
  #  lines(1:dim(bvyears)[1],migryears[,i,r],col=alpha("Black",0.1)) # Number of migrants
  #}
  lower <- apply(migryears[,i,],FUN=mean,MARGIN=1,na.rm=T) - apply(migryears[,i,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(dests)[3]) # Summer
  upper <- apply(migryears[,i,],FUN=mean,MARGIN=1,na.rm=T) + apply(migryears[,i,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(dests)[3]) # Summer
  lower[which(lower<0)] <- 0
  upper[which(upper>1)] <- 1
  polygon(c(1:dim(dests)[1],dim(dests)[1]:1),c(lower,rev(upper)),col=alpha("Black",alpha=0.1),border=NA)
  
  lines(1:dim(bvyears)[1],apply(migryears[,i,],FUN=mean,MARGIN=1),lwd=2,col="Black") # Mean number of migrants
  abline(v=which(years==9950),lty=3,lwd=2)
  mtext(paste0("(",letters[i+1],")"),2,line=2.3,las=1,padj=-10.9,cex=1.1)
  axis(1,at=c(1,(1:4)*25),labels=years[c(1,(1:4)*25)])
}
#legend("topright",legend=c(sites[1:dim(bvyears)[2]],"Migrants"),pch=15,col=c(pal[1:dim(bvyears)[2]],"Black"),box.lty=0,title="Destination alleles",bg=alpha("white",alpha=0.5))
dev.off()
#####

#### Fig S8: Detailed dynamics following ECE in C, f=2. ####
tmp <- readRDS(file="M:\\Documents\\Migration\\sexresults\\Shock exploration n=3\\Winter gradient\\Moderate gradient, c=2, shock_C=0.8 evoldyn.R")
evoldyn <- tmp[[1]]
bvyears <- tmp[[2]]
migryears <- tmp[[3]]
dests <- tmp[[4]]
years <- tmp[[5]]
N <- dim(bvyears)[2]
pdf("M:\\Documents\\Migration\\Publication figs Ecology\\Fig S7 CI.pdf",width=7.5,height=8)
par(mfrow=c(2,2),mar=c(4,4,2,1),oma=c(0,0.4,0,0))
sites <- c("North","Middle","South")
pal <- viridis(N)

plot(1:length(years),rep(0,length(years)),ylim=range(evoldyn$meanbv,na.rm=T),xaxt="n",ylab=expression(paste("Breeding values, ",italic(a))),xlab="Year",type="l",lty=2,lwd=2,main="Metapopulation")
for(n in 1:dim(bvyears)[2]){
  #for(r in 1:dim(bvyears)[3]){
  #  lines(1:dim(bvyears)[1],bvyears[,n,r],col=alpha(pal[n],0.2))
  #}
  lowersum <- apply(bvyears[,n,],FUN=mean,MARGIN=1,na.rm=T) - apply(bvyears[,n,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(bvyears)[3]) # Summer
  uppersum <- apply(bvyears[,n,],FUN=mean,MARGIN=1,na.rm=T) + apply(bvyears[,n,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(bvyears)[3]) # Summer
  polygon(c(1:75,75:1),c(lowersum,rev(uppersum)),col=alpha(pal[n],alpha=0.1),border=NA)
  
  lines(1:dim(bvyears)[1],apply(bvyears[,n,],FUN=mean,MARGIN=1),lwd=2,col=pal[n]) # Average across reps
}
abline(v=which(years==9950),lty=3,lwd=2)
axis(1,at=c(1,(1:4)*25),labels=years[c(1,(1:4)*25)])
#legend("topleft",legend=sites[1:dim(bvyears)[2]],col=pal[1:dim(bvyears)[2]],pch=15,title="Site",bg=alpha("white",alpha=0.7),box.lty=0)
mtext("(a)",2,line=2.3,las=1,padj=-10.9,cex=1.1)

for(i in 1:dim(bvyears)[2]){ # For each home patch
  plot(1:length(years),rep(0,length(years)),ylim=c(0,1),xaxt="n",ylab="Frequency",xlab="Year",type="n",main=paste(sites[i],"subpopulation")) # 1 plot per subpopulation for this one?
  for(n in 1:dim(bvyears)[2]){ # For each going-to patch
    #for(r in 1:dim(bvyears)[3]){ # For each rep
    #  lines(1:dim(bvyears)[1],dests[,i,r,n],col=alpha(pal[n],0.3))
    #}
    lower <- apply(dests[,i,,n],FUN=mean,MARGIN=1,na.rm=T) - apply(dests[,i,,n],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(dests)[3]) # Summer
    upper <- apply(dests[,i,,n],FUN=mean,MARGIN=1,na.rm=T) + apply(dests[,i,,n],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(dests)[3]) # Summer
    lower[which(lower<0)] <- 0
    upper[which(upper>1)] <- 1
    polygon(c(1:dim(dests)[1],dim(dests)[1]:1),c(lower,rev(upper)),col=alpha(pal[n],alpha=0.1),border=NA)
    lines(1:dim(dests)[1],apply(dests[,i,,n],FUN=mean,MARGIN=1),lwd=2,col=pal[n]) # Average across reps
  }
  lines(1:dim(bvyears)[1],rep(0,dim(bvyears)[1]),lwd=2,col="White")
  #for(r in 1:dim(bvyears)[3]){
  #  lines(1:dim(bvyears)[1],migryears[,i,r],col=alpha("Black",0.1)) # Number of migrants
  #}
  lower <- apply(migryears[,i,],FUN=mean,MARGIN=1,na.rm=T) - apply(migryears[,i,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(dests)[3]) # Summer
  upper <- apply(migryears[,i,],FUN=mean,MARGIN=1,na.rm=T) + apply(migryears[,i,],FUN=sd,MARGIN=1,na.rm=T)*qnorm(0.975)/sqrt(dim(dests)[3]) # Summer
  lower[which(lower<0)] <- 0
  upper[which(upper>1)] <- 1
  polygon(c(1:dim(dests)[1],dim(dests)[1]:1),c(lower,rev(upper)),col=alpha("Black",alpha=0.1),border=NA)

  lines(1:dim(bvyears)[1],apply(migryears[,i,],FUN=mean,MARGIN=1),lwd=2,col="Black") # Mean number of migrants
  abline(v=which(years==9950),lty=3,lwd=2)
  mtext(paste0("(",letters[i+1],")"),2,line=2.3,las=1,padj=-10.9,cex=1.1)
  axis(1,at=c(1,(1:4)*25),labels=years[c(1,(1:4)*25)])
}
#legend("topright",legend=c(sites[1:dim(bvyears)[2]],"Migrants"),pch=15,col=c(pal[1:dim(bvyears)[2]],"Black"),box.lty=0,title="Destination alleles",bg=alpha("white",alpha=0.5))
dev.off()
#####

#### Comparison of scenarios: grain & non-genetic liab components ####
pdf(file="M:\\Documents\\Migration\\Publication figs Ecology\\Fig S2 survdiffs.pdf",width=8.5,height=3.5)
par(mfrow=c(1,3),mar=c(4,4.1,3.9,0.7))
sites <- c("North","Middle","South")
plot(0,0,type="n",xlim=c(-0.3,0.3),ylim=c(0,15),xlab="",ylab="Density",main="(a) Baseline",bty="n",cex.lab=1.2,cex.main=1)
alldb <- readRDS(file="M:\\Documents\\Migration\\sexresults\\Only winter gradient alldb.R")
seldiffs <- numeric(100)
for(i in 3:1){
  tic <- 1
  for(y in seq(9901,10000,by=5)){
    for(r in 1:5){
      rtmp <- subset(alldb,Year==y & Rep==r & Home==i & migrated==0)
      rsurv <- mean(as.logical(rtmp$Mort),na.action=na.exclude)
      mtmp <- subset(alldb,Year==y & Rep==r & Home==i & migrated==1)
      msurv <- mean(as.logical(mtmp$Mort),na.action=na.exclude)
      seldiffs[tic] <- rsurv-msurv
      tic <- tic+1
      print(tic)
    }
  }
  hist(seldiffs,xlab="",freq=F,lty=1,angle=i*60,density=20,lwd=2,breaks=seq(-0.3,0.3,by=0.1*1/3),col=alpha(pal_v[i],0.3),border=alpha(pal_v[i],0.9),add=T)
  points(seldiffs,rep(-0.3,100),pch="|",col=alpha(pal_v[i],0.3),cex=1)
}
abline(v=0,lty=2,lwd=2)

alldb <- readRDS(file="M:\\Documents\\Migration\\sexresults\\Fine-grain\\Only winter gradient alldb.R")
plot(0,0,type="n",xlim=c(-0.3,0.3),ylim=c(0,15),xlab="Survival differentials",ylab="",main="(b) Fine-grained environment",bty="n",cex.lab=1.2,cex.main=1)
mtext("Model scenario",side=3,line=2.4)
for(i in 3:1){
  tic <- 1
  for(y in (1:20)*5){
    for(r in 1:5){
      rtmp <- subset(alldb,Year==y & Rep==r & Home==i & migrated==0)
      rsurv <- mean(as.logical(rtmp$Mort),na.action=na.exclude)
      mtmp <- subset(alldb,Year==y & Rep==r & Home==i & migrated==1)
      msurv <- mean(as.logical(mtmp$Mort),na.action=na.exclude)
      seldiffs[tic] <- rsurv-msurv
      tic <- tic+1
      print(tic)
    }
  }
  hist(seldiffs,main=bquote(paste(.(sites[i])," subpopulation")),xlab="Survival differentials",freq=FALSE,
       lty=1,angle=i*60,density=20,lwd=5,breaks=seq(-0.3,0.3,by=0.1*1/3),col=alpha(pal_v[i],0.3),border=alpha(pal_v[i],0.9),add=T)
  points(seldiffs,rep(-0.3,100),pch="|",col=alpha(pal_v[i],0.3),cex=1)
}
abline(v=0,lty=2,lwd=2)
legend("right",legend=sites,col=pal_v[1:3],pch=12,bty="n",title="Subpopulation",bg=alpha(pal_v[1:3],0.3))

alldb <- readRDS(file="M:\\Documents\\Migration\\sexresults\\lambda_d=0, sigma_e=0\\Only winter gradient alldb.R")
plot(0,0,type="n",xlim=c(-0.3,0.3),ylim=c(0,15),xlab="",ylab="",main="(c) No non-genetic components to liability",bty="n",cex.lab=1.2,cex.main=1)

for(i in 3:1){
  tic <- 1
  for(y in (1:20)*5){
    for(r in 1:5){
      rtmp <- subset(alldb,Year==y & Rep==r & Home==i & migrated==0)
      rsurv <- mean(as.logical(rtmp$Mort),na.action=na.exclude)
      mtmp <- subset(alldb,Year==y & Rep==r & Home==i & migrated==1)
      msurv <- mean(as.logical(mtmp$Mort),na.action=na.exclude)
      seldiffs[tic] <- rsurv-msurv
      tic <- tic+1
      print(tic)
    }
  }
  hist(seldiffs,main=bquote(paste(.(sites[i])," subpopulation")),xlab="Survival differentials",freq=F,add=T,
       lty=1,angle=i*60,density=20,lwd=5,breaks=seq(-0.3,0.3,by=0.1*1/3),col=alpha(pal_v[i],0.3),border=alpha(pal_v[i],0.9))
  points(seldiffs,rep(-0.3,100),pch="|",col=alpha(pal_v[i],0.3),cex=1)
}
abline(v=0,lty=2,lwd=2)
dev.off()
#####