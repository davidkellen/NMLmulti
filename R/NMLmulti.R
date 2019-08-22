
init <- function(ks,Ns){
  for (t in 1:length(ks)) {
    m <- t(rmultinom(1,Ns[t],MCMCpack::rdirichlet(1,rep(1/ks[t],ks[t])))) # sample y from dirichlet distribution with alpha = 1/k
    if (t==1) n <- m else n <- c(n,m)
  }
  return(n)
}


# function that returns running average and ESS (plus count)  
welford <- function(x,m1m2count){
    m1    <- m1m2count[1]
    m2    <- m1m2count[2]
    count <- m1m2count[3]
    
    for(iii in 1:length(x)){
      delta <- x[iii]-m1
      m1 <- m1 + delta/count
      m2 <- m2 + delta*(x[iii]-m1)
      count <- count + 1  
    }
    return(c(m1,m2,count))
  }
  

##########################################################
# function that fits a model and returns -log-likelihood #
##########################################################

fit_model <- function(chd,NN,fun,parl) return(nlminb(rnorm(parl),fun,chd=chd,NN=NN)$objective)

######################################
# function that puts it all together #
######################################

run_nml <- function(fun,parl,ks,Ns,batchsize=2000,burn=1000, thin=1,precision=0.2, cores=NULL){
  
  if(is.null(cores)) cores <- round(parallel::detectCores()*0.8)
  cl <- parallel::makePSOCKcluster(cores)
  
  NN <- rep(Ns,ks)
  parallel::clusterExport(cl, list("fit_model"))
  parallel::clusterExport(cl, list("fun","NN"), envir=environment())
  
    
  # burnin samples
  startvec <- list()
  for(ccc in 1:cores){ startvec[[ccc]] <- burnin(burn,init(ks,Ns),ks,Ns);print(ccc)}

  lbatch <- list()
  for(ccc in 1:cores) lbatch[[ccc]] <- gen_chain(ks,Ns,batchsize,thin,startvec[[ccc]])

  lfit <- parallel::parLapply(cl=cl, X=lbatch, function(x) exp(-apply(x,1,fit_model,NN=NN,fun=fun,parl=10)))
  bbatchh <- lapply(lfit,welford,m1m2count=c(0,0,1))
  estbatch <- matrix(unlist(bbatchh),ncol=3,byrow=TRUE)
  if(nrow(estbatch)> 1){poolvarmean <- mean((mean(estbatch[,1])-estbatch[,1])**2)}else{poolvarmean <- estbatch[1,2]/(estbatch[1,3]-1)}
  

  # compute variance and errors
  mm <- mean(estbatch[,1])
  vv <- poolvarmean
  # compute error
  err <- sqrt(poolvarmean)*3
  er1<-log(mm+err)-log(mm)
  if((mm-err)>0){ er2<-log(mm)-log(mm-err)} else {er2<-10.0}
  err<-er1+er2
  
  print(round(c("penalty"=log(mm),"error"=err,"count"=sum(estbatch[,3]-1)),4))

  
 while(err > precision){
  
  lastit <- list()
  for(ccc in 1:cores) lastit[[ccc]] <- lbatch[[ccc]][nrow(lbatch[[ccc]]),]
    
  lbatch <- list()
  for(ccc in 1:cores) lbatch[[ccc]] <- gen_chain(ks,Ns,batchsize,thin,lastit[[ccc]])
  
  lfit <- parLapply(cl=cl, X=lbatch, function(x) exp(-apply(x,1,fit_model,NN=NN,fun=fun,parl=10)))
  
  for(ccc in 1:cores) bbatchh[[ccc]] <- welford(lfit[[ccc]],m1m2count=estbatch[ccc,])
  estbatch <- matrix(unlist(bbatchh),ncol=3,byrow=TRUE)

  if(nrow(estbatch)> 1){poolvarmean <- mean((mean(estbatch[,1])-estbatch[,1])**2)}else{poolvarmean <- estbatch[1,2]/(estbatch[1,3]-1)}
  mm <- mean(estbatch[,1])
  vv <- poolvarmean
  
  # compute error
  err <- sqrt(poolvarmean)*3
  er1<-log(mm+err)-log(mm)
  if((mm-err)>0){ er2<-log(mm)-log(mm-err)} else {er2<-10.0}
  err<-er1+er2
  
  print(round(c("penalty"=log(mm),"error"=err,"count"=sum(estbatch[,3]-1)),4))
 } 
 
  parallel::stopCluster(cl)
  return(round(c("penalty"=log(mm),"error"=err,"count"=sum(estbatch[,3]-1)),4)) 

}  


