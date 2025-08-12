# Nimble functions --------------

# Modified dDHMM function from nimble Ecology to get a weighted likelihood
#Weight is included in NUMi (for each unique life history)
dDHMMmy <- nimbleFunction(
  run = function(x = double(1),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 agei = double(1),
                 NUMi = double(),
                 len = double(),## length of x (needed as a separate param for rDHMM)
                 checkRowSums = double(0, default = 1),
                 log = integer(0, default = 0)) {
    if (length(init) != dim(probObs)[1]) stop("In dDHMM: Length of init does not match nrow of probObs in dDHMM.")
    if (length(x) != len) stop("In dDHMM: Length of x does not match len in dDHMM.")
    if (abs(sum(init) - 1) > 1e-6) stop("In dDHMM: Initial probabilities must sum to 1.")
    
    
    pi <- init # State probabilities at time t=1
    logL <- 0
    lengthX <- length(x)
    for (t in 1:lengthX) {
      Zpi <- pi # Vector of P(state) * P(observation class x[t] | state)
      if(t>1){Zpi <- probObs[, x[t],agei[t]] * pi}
      sumZpi <- sum(Zpi)    # Total P(observed as class x[t])
      logL <- logL + log(sumZpi)*NUMi  # Accumulate log probabilities through time
      if (t != lengthX) pi <- ((Zpi %*% probTrans[,,(t)])/sumZpi)[1, ] # State probabilities at t+1
    }
    
    returnType(double())
    if (log) return(logL)
    return(exp(logL))
  }
)

#This function has to be defined but should never be runned as it is not correct!!!
rDHMMmy <- nimbleFunction(
  run = function(n = integer(),    ## Observed capture (state) history
                 init = double(1),
                 probObs = double(3),
                 probTrans = double(3),
                 agei = double(1),
                 NUMi = double(),
                 len = double(),
                 checkRowSums = double(0, default = 1)) {
    nStates <- length(init)
    if (nStates != dim(probObs)[1]) stop("In rDHMM: Length of init does not match nrow of probObs in dDHMM.")
    if (abs(sum(init) - 1) > 1e-6) stop("In rDHMM: Initial probabilities must sum to 1.")
    if (checkRowSums) {
      transCheckPasses <- TRUE
      obsCheckPasses <- TRUE
      for (i in 1:dim(probObs)[1]) {
        thisCheckSum <- sum(probObs[i,,1])
        if (abs(thisCheckSum - 1) > 1e-6) {
          print("In rDHMM: Problem with sum(probObs[i,]) with i = ", i, ". The sum should be 1 but is ", thisCheckSum)
          obsCheckPasses <- FALSE
        }
      }
      if(!(transCheckPasses | obsCheckPasses))
        stop("In rDHMM: probTrans and probObs were not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!transCheckPasses)
        stop("In rDHMM: probTrans was not specified correctly.  Probabilities in each row (second dimension) must sum to 1.")
      if(!obsCheckPasses)
        stop("In rDHMM: probObs was not specified correctly. Probabilities in each row must sum to 1.")
    }
    
    returnType(double(1))
    ans <- numeric(len)
    
    trueState <- rcat(1, init)
    for (i in 1:len) {
      # Detect based on the true state
      ans[i] <- rcat(1, probObs[trueState,,1])
      # Transition to a new true state
      if (i != len) {
        trueState <- rcat(1, probTrans[trueState, , i])
      }
    }
    return(ans)
  })

registerDistributions(list(
  dDHMMmy = list(
    BUGSdist = 'dDHMMmy( init, probObs, probTrans, agei, NUMi, len, checkRowSums)',
    types    = c('value = double(1)', 'init = double(1)', 'probObs = double(3)', 
                 'probTrans = double(3)',
                 'agei = double(1)','NUMi = double()', 
                 'len = double()', 'checkRowSums = double(0)'))
))

# Other functions --------------

#Estimate Rhat
Rhatfun<-function(rb,nch,m,long){
  vari=apply(rb,c(2,3),var)
  mea=apply(rb,c(2,3),mean)
  meag=apply(rb,3,mean)
  
  W=apply(vari,2,mean)
  W1=apply(vari,2,var)
  B=apply(mea,2,var)
  cov1=cov2=rep(0,long)
  for (i in 1:long){
    cov1[i]=cov(vari[,i],y=(mea^2)[,i])
    cov2[i]=cov(vari[,i],y=mea[,i])
  }
  sig2=((m-1)/m)*W+B
  V=sqrt(sig2+B/3)^2
  varV=((m-1)/m)^2/3*W1+(4/3)^2*B^2+2*(m-1)*4/(9*m)*(cov1-2*meag*cov2)
  df=2*V^2/varV
  Rhat=abs((V/W*df)/(df-2))
  return(Rhat)
}

# Estimate mean, sd, 95%IC and Rhat of nimble outputs
sum_nim<-function(rb2){
  m=dim(rb2)[1]
  nch=dim(rb2)[2]
  long=dim(rb2)[3]
  sumres=matrix(NA,nrow=dim(rb2)[3],ncol=5)
  rownames(sumres)=dimnames(rb2)[[3]]
  colnames(sumres)=c('mean','sd','QI 2.5','QI 97.5','Rhat')
  sumres[,1]=   apply(rb2,3, mean, na.rm = T)
  sumres[,2]=   apply(rb2,3, sd)
  sumres[,3]=   apply(rb2,3, quantile,0.025)
  sumres[,4]=   apply(rb2,3,quantile,0.975)
  sumres[,5]=   Rhatfun(rb2,nch,m,long)
  return(sumres)
}