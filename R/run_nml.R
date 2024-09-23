

# function that generates initial sample for multinomial probabilities
# (replaces old version)
init <- function(ks,Ns){
  for (t in 1:length(ks)) {
    m <- t(rmultinom(1,Ns[t],MCMCpack::rdirichlet(1,rep(1,ks[t])))) # sample y from dirichlet distribution with alpha = 1
    if (t==1) n <- m else n <- c(n,m)
  }
  return(n)
}


# function that computes log-likelihood for ML of saturated model
# (new function)
satMaxLik <- function(x, ks, Ns, NN){
  subs <- 0
  y <- 0
  for(i in 1:length(ks)){
    subslow <- subs + 1
    subs <- subs + ks[i]
    y <- y + dmultinom(x[subslow:subs], Ns[i], x[subslow:subs]/NN[subslow:subs], log = T)
  }
  return(y)
}


# numerically stable function that computes log(exp(xa)+exp(xb))
# (new function)
logsum <- function(xa, xb){
  if(xa <= -Inf){return(xb)}else if(xb <= -Inf){return(xa)}
  if(xa > xb){
    temp <- xa + log1p(exp(xb-xa))
  } else {
    temp <- xb + log1p(exp(xa-xb))
  }
  return(temp)
}


# numerically stable function that computes log(exp(xa)-exp(xb))
# (new function)
logdif <- function(xa, xb){
  if(xa <= -Inf){return(xb)}else if(xb <= -Inf){return(xa)}
  if(xa > xb){
    temp <- xa + log1p(-(exp(xb-xa)))
  } else {
    temp <- xb + log1p(-(exp(xa-xb)))
  }
  return(temp)
}


# iterative version of logsum that takes a vector as argument
# (new function)
logsumVec <- function(logx) {
  y <- logsum(logx[1], logx[2])
  sumlen <- length(logx)
  if(sumlen <= 2){return(y)}else{
    for (i in 3:sumlen) {
      y <- logsum(y, logx[i])
    }
    return(y)
  }
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


fit_model <- function(chd,NN,fun,parl,fits){
  start_par_mat <- matrix(rnorm(parl*fits), nrow=fits)
  t_fits <- apply(start_par_mat,1,function(x) nlminb(x,dynGet(fun),chd=chd,NN=NN)$objective)
  return(return(min(t_fits)))
}

######################################
# function that puts it all together #
######################################
# (replaces old version)

run_nml <- function(
    fun,parl,ks,Ns,
    fits=2,batchsize=2000,burn=1000,
    thin=1,precision=0.2,
    cores=NULL, packages_multicore = NULL,
    G2=TRUE,
    stopcrit=c("betw_logmean","betw_meanlog","betw_max","within_max","all_max"),
    aggmethod=c("logmean","meanlog"))
{ # G2 is a new indicator: if T fun must return 1/2 * G^2, if F fun must return negative log-likelihood

  fun <- deparse(substitute(fun))
  
  aggmethod <- match.arg(aggmethod)
  stopcrit <- match.arg(stopcrit)
  if(is.null(cores)){ cores <- round(parallel::detectCores()*0.8) } else { cores <- round(cores)}
  if (cores == 1) {
    stop("'cores = 1' is curretnly not supported.", call. = TRUE)
  }
  
  cl <- parallel::makePSOCKcluster(cores)
  
  parallel::clusterEvalQ(cl=cl,library("Rcpp"))
  parallel::clusterEvalQ(cl=cl,library("NMLmulti"))
  if (!is.null(packages_multicore)) {
    for (p in packages_multicore) {
      parallel::clusterExport(cl, "p", envir=environment())
      parallel::clusterEvalQ(cl=cl,library(p, character.only = TRUE))
    }
  }
  
  NN <- rep(Ns,ks)
  parallel::clusterExport(cl, list("fit_model","fun","NN","fits","parl","Ns","ks", "burn", "batchsize", "thin"), envir=environment())   # adaptation for parallel sampling

  inits <- list()
  for(ccc in 1:cores){ inits[[ccc]] <- init(ks,Ns)}
  startvec <- parallel::parLapply(cl=cl, X=inits, function(x) burnin(burn,x,ks,Ns))                     # this allows parallel burins
  lbatch <- parallel::parLapply(cl=cl, X=startvec, function(x) gen_chain(ks,Ns,batchsize,thin,x))       # this allows parallel sampling


  lfit <- parallel::parLapply(cl=cl, X=lbatch, function(x) (-apply(x,1,fit_model,NN=NN,fits=fits, fun=fun,parl=parl)))        # this is transforms results not in probability space but remains in log-probability space



  if (G2) {                                                                                                                   # this checks if G^2 or negative log-likelihood was returned by fun and proceeds accordingly
    gfit <- Map("*", lfit, 0.5)
  } else {
    lfitsat <- parallel::parLapply(cl=cl, X=lbatch, function(x) (apply(x,1,satMaxLik,ks=ks, Ns=Ns, NN=NN)))                   # if neg-log-likelihood is returned by fun, the log-likelihood of the saturated model is computed
    gfit <- Map("-",lfit,lfitsat)                                                                                             # likeihood-ratio (LR) is computed
  }

  bbatchh <- lapply(gfit,logsumVec)                                             # this computes the the sum of of the LR of the samples
  bbatchh <- Map(c, bbatchh, lapply(gfit, length))
  estbatch <- matrix(unlist(bbatchh),ncol=2,byrow=TRUE)

  estbatchmean <- estbatch[,1]-log(estbatch[,2])                                # this computes the mean of of the LR of the samples

  SqDif_2_gfit <- lapply(gfit, function(x) numeric(length(x)))                  # this computes the standard errors of the estimate within cores
  for (j in 1:length(SqDif_2_gfit)) {
    for(i in 1:length(SqDif_2_gfit[[j]])){
      SqDif_2_gfit[[j]][i] <- 2*logdif(gfit[[j]][i], estbatchmean[j])
    }
  }
  SS_2_gfit <- lapply(SqDif_2_gfit, logsumVec)
  poolSEmean_within <- unlist(lapply(SqDif_2_gfit, function(x) (0.5*(logsumVec(x)-log(length(x)-1)))-(0.5*log(length(x))) ))

  errup_within <- apply(cbind(estbatchmean, poolSEmean_within), 1, function(x) logsum(x[1], x[2]+log(3)))
  errlo_within <- apply(cbind(estbatchmean, poolSEmean_within), 1, function(x) logdif(x[1], x[2]+log(3)))
  errint_within <- errup_within - errlo_within


  # compute variance and errors
  mm_1 <- mean(estbatchmean)                                                             # this computes the mean of log(T_z)

  if(cores>1){

    mm_2 <- logsumVec(estbatchmean) - log(length(estbatchmean))                          # this computes the log of mean(T_z)

    vv_1 <- sd(estbatchmean)                                                             # this computes the SD of the estimates of the different cores

    SS_2_estbatchmean <- numeric(length(estbatchmean))                                   # this computes the SD of T_z
    for(i in 1:length(estbatchmean)){
      SS_2_estbatchmean[i] <- 2*logdif(estbatchmean[i], mm_2)
    }

    vv_2 <- 0.5*(logsumVec(SS_2_estbatchmean)-log(length(estbatchmean)-1))

    # compute error
    erup_1 <- mm_1+vv_1*3
    erup_2 <- logsum(mm_2, vv_2+log(3))

    erlo_1 <- mm_1-vv_1*3
    erlo_2 <- logdif(mm_2, vv_2+log(3))

    errint_1 <- erup_1-erlo_1
    errint_2 <- erup_2-erlo_2

    if(stopcrit=="betw_max"){
      errint <- max(errint_1, errint_2)
    }else if(stopcrit=="betw_logmean"){
      errint <- errint_2
    }else if(stopcrit=="betw_meanlog"){
      errint <- errint_1
    }else if(stopcrit=="within_max"){
      errint <- max(errint_within)
    }else if(stopcrit=="all_max"){
      errint <- max(errint_1, errint_2, max(errint_within))
    }

    # compute penalty
    if(aggmethod=="logmean"){
      mm <- mm_2
    }else if(aggmethod=="meanlog"){
      mm <- mm_1
    }

    print(round(c("penalty (mean log)"=mm_1, "penalty (log mean)"=mm_2,"error (mean log)"=errint_1, "error (log mean)"=errint_2, "max error within"= max(errint_within), "count"=sum(estbatch[,2])),4))


  } else {
    # compute penalty
    mm <- mm_1

    # compute error
    errint_1 <- NA
    errint_2 <- NA
    errint <- errint_within

    print(round(c("penalty"=mm, "error within"=errint_within, "count"=sum(estbatch[,2])),4))
  }

  while(errint > precision){

    lastit <- list()
    for(ccc in 1:cores) lastit[[ccc]] <- lbatch[[ccc]][nrow(lbatch[[ccc]]),]

    lbatch <- parallel::parLapply(cl=cl, X=lastit, function(x) gen_chain(ks,Ns,batchsize,thin,x))

    lfit <- parallel::parLapply(cl=cl, X=lbatch, function(x) (-apply(x,1,fit_model,NN=NN,parl=parl,fits=fits,fun=fun)))

    if (G2) {
      gfit <- Map("*", lfit, 0.5)
      #gfit = lfit      # this would be accurate if fun returns 1/2 * G^2
    } else {
      lfitsat <- parallel::parLapply(cl=cl, X=lbatch, function(x) (apply(x,1,satMaxLik,ks=ks, Ns=Ns, NN=NN)))
      gfit <- Map("-",lfit,lfitsat)
    }

    estbatch_alt = estbatch

    bbatchh <- lapply(Map(c, estbatch_alt[,1], gfit), logsumVec)
    bbatchh <- Map(c, bbatchh, Map("+",lapply(gfit, length), estbatch_alt[,2]))
    estbatch <- matrix(unlist(bbatchh),ncol=2,byrow=TRUE)

    estbatchmean = estbatch[,1]-log(estbatch[,2])


    SS_2_gfit_old <- SS_2_gfit
    SqDif_2_gfit <- lapply(gfit, function(x) numeric(length(x)))
    for (j in 1:length(SqDif_2_gfit)) {
      for(i in 1:length(SqDif_2_gfit[[j]])){
        SqDif_2_gfit[[j]][i] <- 2*logdif(gfit[[j]][i], estbatchmean[j])
      }
    }
    SqDif_2_gfit <- Map(c, SS_2_gfit_old, SqDif_2_gfit)

    for (ccc in 1:cores) {
      SS_2_gfit[[ccc]] <- logsumVec(SqDif_2_gfit[[ccc]])
      poolSEmean_within[ccc] <- (0.5*(SS_2_gfit[[ccc]] - log(estbatch[ccc,2]-1))) - (0.5*log(estbatch[ccc,2]))
    }

    errup_within <- apply(cbind(estbatchmean, poolSEmean_within), 1, function(x) logsum(x[1], x[2]+log(3)))
    errlo_within <- apply(cbind(estbatchmean, poolSEmean_within), 1, function(x) logdif(x[1], x[2]+log(3)))
    errint_within <- errup_within - errlo_within



    # compute variance and errors
    mm_1 <- mean(estbatchmean)

    if(cores>1){

      mm_2 <- logsumVec(estbatchmean) - log(length(estbatchmean))

      vv_1 <- sd(estbatchmean)

      SS_2_estbatchmean <- numeric(length(estbatchmean))
      for(i in 1:length(estbatchmean)){
        SS_2_estbatchmean[i] <- 2*logdif(estbatchmean[i], mm_2)
      }

      vv_2 <- 0.5*(logsumVec(SS_2_estbatchmean)-log(length(estbatchmean)-1))

      # compute error
      erup_1 <- mm_1+vv_1*3
      erup_2 <- logsum(mm_2, vv_2+log(3))

      erlo_1 <- mm_1-vv_1*3
      erlo_2 <- logdif(mm_2, vv_2+log(3))

      errint_1 <- erup_1-erlo_1
      errint_2 <- erup_2-erlo_2

      if(stopcrit=="betw_max"){
        errint <- max(errint_1, errint_2)
      }else if(stopcrit=="betw_logmean"){
        errint <- errint_2
      }else if(stopcrit=="betw_meanlog"){
        errint <- errint_1
      }else if(stopcrit=="within_max"){
        errint <- max(errint_within)
      }else if(stopcrit=="all_max"){
        errint <- max(errint_1, errint_2, max(errint_within))
      }

      # compute penalty
      if(aggmethod=="logmean"){
        mm <- mm_2
      }else if(aggmethod=="meanlog"){
        mm <- mm_1
      }

      print(round(c("penalty (mean log)"=mm_1, "penalty (log mean)"=mm_2,"error (mean log)"=errint_1, "error (log mean)"=errint_2, "max error within"= max(errint_within), "count"=sum(estbatch[,2])),4))

    } else {
      # compute penalty
      mm <- mm_1

      # compute error
      errint_1 <- NA
      errint_2 <- NA
      errint <- errint_within

      print(round(c("penalty"=mm_2, "error within"= errint, "count"=sum(estbatch[,2])),4))

    }


  }

  parallel::stopCluster(cl)

  if(cores==1){
    return((list("penalty"=mm,"error_within"=errint,"count"=sum(estbatch[,2]-1),"stopping_criterion" = "within","aggregation_method" = "none")))
  } else {
    return((list("penalty"=mm,"error"=errint,"max_error_within"= max(errint_within),"count"=sum(estbatch[,2]-1), "stopping_criterion" = stopcrit,"aggregation_method" = aggmethod)))
  }

}
