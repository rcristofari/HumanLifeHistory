# Dataset: one vector of (residual) age-at-death and one boolean vector with censorship.
# The algorithm is an extension to randomly right-censored data of the one proposed by Jiang & Kececioglu (1992)
# The module contains 4 main functions:
# SimWeibull(...)
# FitWeibull(...)
# PlotWeibull(...)
# CIWeibull(...)

# TODO list : multithread bootstrapping
# Graphical guess for multimodal data

# ChooseNWeibull(...)

##############################################################################################################
# CENSORED MIXTURE WEIBULL SIMULATOR
##############################################################################################################
# Simulate an N-population mixed-weibull dataset:

SimWeibull <- function(nSamp=1000, rightcensor=.2, leftcensor = 0, theta){
  
# Fix the missing degree of freedom on theta through the last weight in case of error:
weights <- theta[seq(1,length(theta),3)]
last_weight <- 1-sum(weights[-length(weights)])
theta[length(theta)-2] <- last_weight

nPop <- length(theta)/3
aad <- vector()
belonging <- matrix(ncol = nPop, nrow = nSamp, data = 0)
begin <- 1
for(p in 1:nPop){
  w <- theta[((p-1)*3+1)] ; a <- theta[((p-1)*3+2)] ; b <- theta[((p-1)*3+3)]
  nspop <- round(nSamp*w)
  belonging[begin:(begin+nspop-1),p] <- 1
  begin <- begin+nspop
  aad <- c(aad, rweibull(nspop, scale=a, shape=b))
}

censor <- sample(c(0,1,2), nSamp, replace=T, prob=c(rightcensor, (1-(rightcensor+leftcensor)), leftcensor))
for(a in 1:nSamp){
  if(censor[a]==0){
    aad[a]<-runif(1,(aad[a]/2),aad[a])
  } else if (censor[a]==2){
    aad[a]<-runif(1,aad[a],max(aad))}
}

return(list(aad=aad, censor=censor, belonging=belonging))}

####################################################################################
####################################################################################

load('~/Desktop/aad.RData')
load('~/Desktop/censor.RData')
simdata <- SimWeibull(nSamp=5000, leftcensor=0, rightcensor=0.4, theta = c(.5, 90, 4, .5, 27, 9))
aad <- as.numeric(simdata[['aad']])
censor <- as.numeric(simdata[['censor']])
belonging <- simdata[['belonging']]
                      
####################################################################################
# MAIN EXPECTATION-MAXIMISATION FITTING FUNCTION
####################################################################################

FitWeibull <- function(aad, censor, n=2, weights=NA, thres=1e-4, NR_thres=1e-6, verbose=TRUE, plot=TRUE){
  # Structure of the theta vector: {w1, a1, b1, w2, a2, b2... wi, ai, bi}
  # censor: 0 for right-censored, 1 for uncensored
  start_time <- proc.time()
  require(Rmpfr)
 
  if(is.na(weights)==T){ 
    censor <- censor[order(aad)]
    aad <- sort(aad)}
  
  # the simple, 2-parameter Weibull PDF function:
  pdf <- function(x,a,b){return((b/a)*((x/a)^(b-1)) * exp(-(x/a)^b))}
  # the simple, 2-parameter Weibull tail function:
  tf <- function(x,a,b){return(exp(-(x/a)^b))}

####################################################################################
# N-population Weibull log-likelihood function
LnL <- function(aad, censor, theta){
	# theta is a vector of model parameters {w1, a1, b1, w2, a2, b2... wi, ai, bi}
	# a are scales, b are shapes (opposite from other R functions but reflects the literature)
  
  LogPrTj <- function(t, c, theta){
  # Summed probability of failure time Tj over all components 'i'
  # the simple, 2-parameter Weibull PDF function for a single component:
  total <- 0
  if(c == 1) {
    for(k in 1:n){
      wi <- theta[((k-1)*3+1)]
      ai <- theta[((k-1)*3+2)]
      bi <- theta[((k-1)*3+3)]
      total <- total + log(wi*pdf(x=t, a=ai, b=bi))}
    return(total)
  } else if (c == 0){
    for(k in 1:n){
      wi <- theta[((k-1)*3+1)]
      ai <- theta[((k-1)*3+2)]
      bi <- theta[((k-1)*3+3)]
      total <- total + log(wi*tf(x=t, a=ai, b=bi))}
  } else if (c == 2){
    for(k in 1:n){
      wi <- theta[((k-1)*3+1)]
      ai <- theta[((k-1)*3+2)]
      bi <- theta[((k-1)*3+3)]
      total <- total + log(wi*cdf(x=t, a=ai, b=bi))}
  } else {
      stop("Invalid value in the censoring vector.")}
  return(total)}

  SumLogPrTj <- 0
  for(j in 1:length(aad)){
    SumLogPrTj <- SumLogPrTj + LogPrTj(aad[j], censor[j], theta)}
  
  # N <- length(censor)
  # R <- sum(censor)
  # Full likelihood function:
  # L <- (lfactorial(N) - lfactorial(N-R)) + SumLogPrTj
  # return(L)
  
  return(SumLogPrTj)}

##############################################################################################################
# For each time "tj", determine the probability of belonging to subp. "i" given a parameter estimate theta:
PRiTj <- function(t, c, i, theta){

  # the simple, 2-parameter Weibull PDF function:
  pdf <- function(x,a,b){return((b/a)*((x/a)^(b-1)) * exp(-(x/a)^b))}
  # the simple, 2-parameter Weibull tail function:
  tf <- function(x,a,b){return(exp(-(x/a)^b))}

  w <- theta[((i-1)*3+1)] ; a <- theta[((i-1)*3+2)] ; b <- theta[((i-1)*3+3)]

  total <- 0
  for(k in 1:n){
    wi <- theta[((k-1)*3+1)]
    ai <- theta[((k-1)*3+2)]
    bi <- theta[((k-1)*3+3)]
    total <- total + wi * (pdf(x=t, a=ai, b=bi) ^ c) * (tf(x=t, a=ai, b=bi) ^ (1 - c)) }
  return((w * (pdf(x=t, a=a, b=b) ^ c) * (tf(x=t, a=a, b=b) ^ (1 - c))) / total )}
  
##############################################################################################################
## FUNCTIONS TO USE DURING THE ITERATION PHASES

# Describe Bi given the belongings probs:
G <- function(b, i, aad, censor, probs) {
  # Using multiple-precision floating-point reliable numbers to get the G1/G3 ratio
  # in case the log goes over 709 to avoid failing on an Inf exponential
  
  # G1 <- function(b, i, aad, probs) sum(probs[,i] * log(aad) * aad^b) # goes do Inf
  
  G1 <- function(b, i, aad, probs){
    total <- 0
    for(j in 1:length(aad)){
      logpoint <- log(probs[j,i]) + log(log(aad[j])) + b*log(aad[j])
      if(logpoint > 709){
      total <- mpfr(total, precBits = 24)
      logpoint <- mpfr(logpoint, precBits = 24)}
      total <- total + exp(logpoint)}
    return(total)}
  
  G2 <- function(i, censor, probs) sum(probs[censor==1,i])
  
  # G3 <- function(b, i, aad, probs) sum(probs[,i] * aad^b) # goes do Inf
  
  G3 <- function(b, i, aad, probs){
    total <- 0
    for(j in 1:length(aad)){
      logpoint <- log(probs[j,i]) + b*log(aad[j])
      if(logpoint > 709){
        total <- mpfr(total, precBits = 24)
        logpoint <- mpfr(logpoint, precBits = 24)}
      total <- total + exp(logpoint)}
    return(total)}
  
  G13 <- function(b, i, aad, probs) as.numeric(G1(b, i, aad, probs) / G3(b, i, aad, probs))
  
  G4 <- function(i, aad, censor, probs) sum(probs[censor==1,i] * log(aad[censor==1]))
  G5 <- function(b, i, censor, probs) sum(probs[censor==1,i])/b
  # G <- (-1 * G1(b, i, aad, probs) * G2(i, censor, probs)) / G3(b, i, aad, probs) + G4(i, aad, censor, probs) + G5(b, i, censor, probs)
  G <- (-1 * G13(b, i, aad, probs) * G2(i, censor, probs)) + G4(i, aad, censor, probs) + G5(b, i, censor, probs)
return(G)}

# Derivative of G to be used in the Newton iteration:
g <- function(b, i, aad, censor, probs) {
  g1 <- function(b, i, aad, censor, probs) (-1/(b^2)) * sum(probs[censor==1,i])
  g2 <- function(i, aad, censor, probs) sum(probs[censor==1,i])
  
  # U <- v <-  function(b, i, aad, probs) sum(probs[,i] * log(aad) * aad^b)
  # u <- function(b, i, aad, probs) sum(probs[,i] * log(aad)^2 * aad^b)
  # V <- function(b, i, aad, probs) sum(probs[,i] * aad^b)

  U <- v <-  function(b, i, aad, probs){
    total <- 0
    for(j in 1:length(aad)){
      logpoint <- log(probs[j,i]) + log(log(aad[j])) + b*log(aad[j])
      if(logpoint > 709){
        total <- mpfr(total, precBits = 24)
        logpoint <- mpfr(logpoint, precBits = 24)}
      total <- total + exp(logpoint)}
    return(total)}
  
  u <-  function(b, i, aad, probs){
    total <- 0
    for(j in 1:length(aad)){
      logpoint <- log(probs[j,i]) + 2*log(log(aad[j])) + b*log(aad[j])
      if(logpoint > 709){
        total <- mpfr(total, precBits = 24)
        logpoint <- mpfr(logpoint, precBits = 24)}
      total <- total + exp(logpoint)}
    return(total)}  
  
  V <-  function(b, i, aad, probs){
    total <- 0
    for(j in 1:length(aad)){
      logpoint <- log(probs[j,i]) + b*log(aad[j])
      if(logpoint > 709){
        total <- mpfr(total, precBits = 24)
        logpoint <- mpfr(logpoint, precBits = 24)}
      total <- total + exp(logpoint)}
    return(total)} 
  
  uVUvV2 <- function(b, i, aad, probs) as.numeric((u(b, i, aad, probs)*V(b, i, aad, probs) - U(b, i, aad, probs)*v(b, i, aad, probs)) / V(b, i, aad, probs)^2)
  
  g <- g1(b, i, aad, censor, probs) - g2(i, aad, censor, probs) * uVUvV2(b, i, aad, probs)
  
  # g <- g1(b, i, aad, censor, probs) - g2(i, aad, censor, probs) * ((u(b, i, aad, probs) * V(b, i, aad, probs) - U(b, i, aad, probs) * v(b, i, aad, probs)) / V(b, i, aad, probs)^2)
return(g)}

# Newton-Raphson iteration to find the root of G, starting from a guess "B0":
# This converges very nicely provided the starting estimates are not too bad.
NR <- function(b, i, aad, censor, probs, nr_thres = NR_thres, verbose=F){
  B0 <- b
  B1 <- B0 - (G(B0, i, aad, censor, probs)/g(B0, i, aad, censor, probs))
  while(abs(B0-B1) >= nr_thres) {
    B0 <- B1
    B1 <- B0 - (G(B0, i, aad, censor, probs)/g(B0, i, aad, censor, probs))
    if(verbose==T) print(B1)
  }
  return(B1)}

##############################################################################################################
# Get an estimate of Ai given an updated estimate of Bi

Ai <- function(b, i, aad, censor, probs){
  A1 <- function(b, i, aad, probs) sum(probs[,i] * aad[]^b)
  A2 <- function(i, aad, censor, probs) sum(probs[censor==1,i])
  return((A1(b, i, aad, probs) / A2(i, aad, censor, probs))^(1/b))}

##############################################################################################################
# Get an estimate of Wi given a updated estimates of Ai and Bi sotred in theta
Wi <- function(i, aad, probs) sum(probs[,i])/length(aad)

####################################################################################
# EMPIRICAL QUANTILE FUNCTION GIVEN A BELONGING OBJECT
####################################################################################

#This ignores the censoring...
emp_quant <- function(p, i, belonging, c, a){
  belonging <- belonging[censor==1,]
  normalise_belonging <- function(x) x/sum(x)
  belonging <- apply(belonging, 2, normalise_belonging)
  with_aad <- data.frame(aad=aad[censor==1], belong = belonging[,i])
  with_aad <- with_aad[order(with_aad$aad),]
  cumbelong <- apply(with_aad, 2, cumsum)
  cumbelong[,1]<-sort(aad[censor==1])
  return(cumbelong[cumbelong[,2]>=p,][1,1])}

# This estimator is taken from Marks 2005, Estimation of Weibull parameters from common percentiles
b_quantile_estimate <- function(i, blg, aad, censor, level=.1){
  qtL <- level
  qtH = 1 - qtL
  L <- emp_quant(qtL, i, belonging=blg, a=aad, c=censor)
  H <- emp_quant(qtH, i, belonging=blg, a=aad, c=censor)
  b <-  log(log(qtL)/log(qtH)) * log(H/L)^(-1)
  return(b)
}


##############################################################################################################
# Initialise an Lx plot if desired
if(plot == TRUE){
  require(survival)
  sfit <- survfit(Surv(aad,censor,type="right")~1)
  lx <- data.frame('age'=summary(sfit)$time, 'lx'=(summary(sfit)[6]))
  plot(lx$age, lx$surv, cex=.05)
}

##############################################################################################################
# ARTICULATE EVERYTHING INTO THE ALGORITHM ITSELF
##############################################################################################################
# Initialise the belonging probabilities randomly after a Dirichlet distribution
if(is.na(weights)==T){
require(LearnBayes)
initial_belonging <- rdirichlet(length(aad), rep(.1, n))
initial_belonging <- initial_belonging[order(initial_belonging[,1]),]
} else { initial_belonging <- weights }
LogLikelihood <- vector()

# Make an empty vector Theta that will hold the initial estimates:
Init_Theta <- vector(length=3*n)

for(i in 1:n){
  b <- b_quantile_estimate(i, blg=initial_belonging, aad=aad, censor=censor, level=.1)
  # Get a first estimate of Bi:
  Init_Theta[(i-1)*3+3] <- NR(b, i, aad=aad, censor=censor, probs=initial_belonging, nr_thres = NR_thres)
  # Get the corresponding Ai:
  Init_Theta[(i-1)*3+2] <- Ai(b, i, aad=aad, censor=censor, probs=initial_belonging)
  # Get the corresponding weight:
  Init_Theta[(i-1)*3+1] <- Wi(i, aad=aad, probs=initial_belonging)
}

# Add in the first LnL value
LogLikelihood <- append(LogLikelihood, LnL(aad=aad, censor=censor, theta=Init_Theta))

# Proceed to the EM iterations:

# Make a matrix to receive the parameter estimates, and fill the first row:
thetas <- matrix(ncol=(3*n), nrow=1)
thetas[1,] <- Init_Theta

# Initialise the convergence criterion
converged <- FALSE
# Count the number of consecutive steps fulfilling convergence criterion:
stable_steps <- 0
k <- 1

while(converged == FALSE){
  this_theta <- vector(length=3*n)
  previous_theta <- thetas[k,]
  
  current_belonging <- matrix(ncol=n, nrow=length(aad))
  for(j in 1:length(aad)){
    for(i in 1:n){
      current_belonging[j,i] <- PRiTj(aad[j], censor[j], i, previous_theta)}}
  #plot(raster(current_belonging))
  
  for(i in 1:n){
    b <- b_quantile_estimate(i, blg=current_belonging, aad=aad, censor=censor, level=.1)
    # Get a first estimate of Bi:
    #this_theta[(i-1)*3+3] <- NR(b, i, aad, censor, theta=previous_theta, thres = thres)
    this_theta[(i-1)*3+3] <- NR(b, i, aad=aad, censor=censor, probs=current_belonging, nr_thres = NR_thres)
    # Get the corresponding Ai:
    #this_theta[(i-1)*3+2] <- Ai(b, i, aad, censor, theta=previous_theta)
    this_theta[(i-1)*3+2] <- Ai(b, i, aad, censor=censor, probs=current_belonging)
    # Get the corresponding weight:
    #this_theta[(i-1)*3+1] <- Wi(i, aad, censor, theta=previous_theta)
    this_theta[(i-1)*3+1] <- Wi(i, aad, probs=current_belonging)
  }

  LogLikelihood <- append(LogLikelihood, LnL(aad=aad, censor=censor, theta=this_theta))
  thetas <- rbind(thetas, this_theta)
  if(verbose==T) print(thetas)
  
  # Evaluate the convergence criterion:
  if(sum(abs(previous_theta - this_theta) < thres)){
    stable_steps <- stable_steps +1 } else { stable_steps <- 0 }
  if(stable_steps >= 3) converged <- TRUE
  
  if(plot == TRUE){
    formula <- paste(rep('%s * exp(-(x/%s)^%s) + ', length(this_theta)/3), sep='', collapse='')
    formula <- do.call(sprintf, c(list(formula), round(this_theta, 3)))
    formula <- substr(formula,1,nchar(formula)-3)
    funform <- function(x) eval(parse(text=formula))
    plot(lx$age, lx$surv, cex=.05)
    curve(funform, from=min(aad), to=max(aad), add=T, col='red')}

  k <- k+1}


# Format and return the output
estvector <- thetas[nrow(thetas),]
BestLogLikelihood <- LogLikelihood[length(LogLikelihood)]
estmatrix <- matrix(ncol = 3, nrow = length(estvector)/3)
for(i in 1:(length(estvector)/3)){
  estmatrix[i,1] <- estvector[(i-1)*3+1]
  estmatrix[i,2] <- estvector[(i-1)*3+2]
  estmatrix[i,3] <- estvector[(i-1)*3+3]}
colnames(estmatrix) <- c('Weight','Scale','Shape')
rownames(estmatrix) <- sprintf('SubPop_%s', 1:(length(estvector)/3))
iterations <- cbind(thetas, LogLikelihood)
wnames <- sprintf('Weight_%s', 1:(length(estvector)/3))
scnames <- sprintf('Scale_%s', 1:(length(estvector)/3))
shnames <- sprintf('Shape_%s', 1:(length(estvector)/3))
columnnames <- vector()
for(x in 1:(length(estvector)/3)){columnnames <- c(columnnames, wnames[x], scnames[x], shnames[x])}
columnnames <- c(columnnames, 'LnL')
colnames(iterations) <- columnnames
rownames(iterations) <- sprintf("Iter_%s", 1:nrow(iterations))
output <- list(BestLogLikelihood, estmatrix, iterations)
names(output) <- c('LogLikelihood', 'Estimates', 'Iterations')

print(paste('Convergence reached in', round((proc.time()[3] - start_time[3])/60), 'minutes, after', k, 'iterations.'))

return(output)
}

####################################################################################
# PLOT A MIXTURE WEIBULL SURVIVAL CURVE
####################################################################################

PlotWeibull <- function(aad, censor, theta, add=FALSE) {
  if(add==FALSE){
    require(survival)
    sfit <- survfit(Surv(aad,censor,type="right")~1)
    lx <- data.frame('age'=summary(sfit)$time, 'lx'=(summary(sfit)[6]))
    plot(lx$age, lx$surv, cex=.05)
    add <- TRUE}
  
  formula <- paste(rep('%s * exp(-(x/%s)^%s) + ', length(theta)/3), sep='', collapse='')
  formula <- do.call(sprintf, c(list(formula), round(theta, 3)))
  formula <- substr(formula,1,nchar(formula)-3)
  funform <- function(x) eval(parse(text=formula))
  #plot(lx$age, lx$surv, cex=.05)
  curve(funform, from=min(aad), to=max(aad), add=add, col='red')}

####################################################################################
# BOOTSTRAP ROUTINE
####################################################################################

CIWeibull <- function(method='nonparametric', theta=NA, n=2, aad=NA, censor=NA, nBs = 100, threads=1, verbose=F, plot=F, thres=1e-4, NR_thres=1e-6){
  # Type is (1) "parametric" ('p') or (2) "nonparametric" ('np')
  # For non-parametric bootstrapping, theta is not needed.
  # For parametric bootstrapping, n is not needed
  # For parametric bootstrapping, aad and censor are used to extract dataset shape
  ptm_total <- proc.time()
  
  # Some controls on passed arguments
  # XXXXX
  
  if(is.na(theta)==FALSE){
  n <- length(theta)/3}
  
  # Single-threaded case does not require parallel libraries for better portability
  if (threads == 1) {

  fits <- list()
  if(method=='nonparametric' | method=='np'){
    for(i in 1:nBs){
      print(paste("Non-parametric bootstrap replicate number", i))
      fit <- NA
      while(is.na(fit)==T){
      index <- sort(sample(1:length(aad), length(aad), replace=T))
      this_aad <- aad[index]
      this_censor <- censor[index]
      try(fit <- FitWeibull(this_aad, this_censor, n=n, weights=NA, thres=thres, NR_thres=NR_thres, verbose=verbose, plot=plot), silent=T)}
      fits[[i]] <- fit}
    
  } else if(method=='parametric' | method=='p'){
    for(i in 1:nBs){
      print(paste("Parametric bootstrap replicate number", i))
      fit <- NA
      while(is.na(fit)==T){
      simdata <- SimWeibull(nSamp=length(aad), rightcensor=((length(censor)-sum(censor))/length(censor)), theta = theta)
      try(fit <- FitWeibull(aad=as.numeric(simdata$aad), censor=as.numeric(simdata$censor), n=n, weights=NA, thres=thres, NR_thres=NR_thres, verbose=verbose, plot=plot), silent=T)}
      fits[[i]] <- fit}
    
  } else { print(paste(type, ': invalid bootstrapping type (only parametric, p or nonparametric, np)'))}

    
  # Multithreaded case uses de doParallel library:
  } else if (threads > 1) {
  
  #Setup the parallel cluster
  library("foreach")
  library("doParallel")
  cl <- makeCluster(3)
  registerDoParallel(cl)

  if(method=='nonparametric' | method=='np'){
  fits = foreach(b = 1:nBs, 
                    .export=c('FitWeibull')
  ) %dopar% {
        fit <- NA
        while(is.na(fit)==T){
          index <- sort(sample(1:length(aad), length(aad), replace=T))
          this_aad <- aad[index]
          this_censor <- censor[index]
          try(fit <- FitWeibull(this_aad, this_censor, n=n, weights=NA, thres=thres, NR_thres=NR_thres, verbose=verbose, plot=plot), silent=T)}
          return(fit)}
      
    } else if(method=='parametric' | method=='p'){
    fits = foreach(b = 1:nBs, 
                     .export=c('FitWeibull', 'SimWeibull')
    ) %dopar% {
        fit <- NA
        while(is.na(fit)==T){
          simdata <- SimWeibull(nSamp=length(aad), rightcensor=((length(censor)-sum(censor))/length(censor)), theta = theta)
          try(fit <- FitWeibull(aad=as.numeric(simdata$aad), censor=as.numeric(simdata$censor), n=n, weights=NA, thres=thres, NR_thres=NR_thres, verbose=verbose, plot=plot), silent=T)}
        return(fit)}

    } else { print(paste(type, ': invalid bootstrapping type (only parametric, p or nonparametric, np)'))}
    
  stopCluster(cl)}
  
  # Format the parameter estimates into a matrix:
  estimates <- matrix(nrow=nBs, ncol=n*3+1)
  for(b in 1:nBs){
    estimates[b,] <- c(fits[[b]]$Iterations[nrow(fits[[b]]$Iterations),], fits[[b]]$BestLogLikelihood)}
  wnames <- sprintf('Weight_%s', 1:n)
  scnames <- sprintf('Scale_%s', 1:n)
  shnames <- sprintf('Shape_%s', 1:n)
  columnnames <- vector()
  for(x in 1:n){columnnames <- c(columnnames, wnames[x], scnames[x], shnames[x])}
  columnnames <- c(columnnames, 'LnL')
  colnames(estimates) <- columnnames
  
  # Get median and quantile series to plot:
  ages <- 1:125
  
  # Cumulative Probability Function given theta
  CDF <- function(ages, theta){
    formula <- paste(rep('%s * exp(-(x/%s)^%s) + ', length(theta)/3), sep='', collapse='')
    formula <- do.call(sprintf, c(list(formula), round(theta, 3)))
    formula <- substr(formula,1,nchar(formula)-3)
    funform <- function(x) eval(parse(text=formula))
    return(funform(ages))}
  
  # Probability Density Function given theta 
  PDF <- function(ages, theta){
    formula <- paste(rep('%s * (%s/%s)*((x/%s)^(%s-1)) * exp(-(x/%s)^%s) + ', length(theta)/3), sep='', collapse='')
    params <- vector()
    for(i in 1:(length(theta)/3)){
      params[(i-1)*7+1] <- theta[(i-1)*3+1]
      params[(i-1)*7+2] <- theta[(i-1)*3+3]
      params[(i-1)*7+3] <- theta[(i-1)*3+2]
      params[(i-1)*7+4] <- theta[(i-1)*3+2]
      params[(i-1)*7+5] <- theta[(i-1)*3+3]
      params[(i-1)*7+6] <- theta[(i-1)*3+2]
      params[(i-1)*7+7] <- theta[(i-1)*3+3]}
    formula <- do.call(sprintf, c(list(formula), round(params, 3)))
    formula <- substr(formula,1,nchar(formula)-3)
    funform <- function(x) eval(parse(text=formula))
    return(funform(ages))}
  
  PDFs <- matrix(nrow=125, ncol=nBs)
  CDFs <- matrix(nrow=125, ncol=nBs)
  
  for(b in 1:nBs){
    PDFs[,b] <- PDF(ages, estimates[b,-length(estimates[b,])])  
    CDFs[,b] <- CDF(ages, estimates[b,-length(estimates[b,])])
  }
  
  # apply median, quant. 2.5%, 5%, 10%, 25%, 75%, 90%, 95%, 97.5%
  PDFquantiles <- apply(PDFs, 1, quantile, c(.025, .05, .1, .25, .5, .75, .90, .95, .975))
  CDFquantiles <- apply(CDFs, 1, quantile, c(.025, .05, .1, .25, .5, .75, .90, .95, .975))

  output <- list(estimates, PDFquantiles, CDFquantiles, fits)
  names(output) <- c('Estimates', 'PDF', 'CDF', 'Fits')
  
 print(paste("Total execution time: ", (proc.time() - ptm_total)[3]/60, " minutes"))
 return(output)}


####################################################################################
# PLOT CONFIDENCE INTERVALS FROM BOOTSTRAP
####################################################################################

PlotCIWeibull <- function(ci, aad=NA, censor=NA, type='surv') {
  require(ggplot2)
  # type is "surv" (the tail function) or "dens" (the density function)
  lx <- NA
  if(is.na(aad)==F && is.na(censor)==F){
    require(survival)
    sfit <- survfit(Surv(aad,censor,type="right")~1)
    lx <- data.frame('age'=summary(sfit)$time, 'lx'=(summary(sfit)[6]))
    plot(lx$age, lx$surv, cex=.05)}
  
  if(type=='dens'){
  plotdata <- as.data.frame(cbind(age=1:125, t(ci$PDF)))
  names(plotdata)<-c('age','q025','q050','q100','q250','q500','q750','q900','q950','q975')
  
  g <- ggplot(plotdata, aes(x=age)) +
    theme_bw() +
    geom_line(aes(y=q500))+
    geom_ribbon(aes(ymin=q025, ymax=q975), alpha=.5)
  
  } else if (type=='surv'){
    plotdata <- as.data.frame(cbind(age=1:125, t(ci$CDF)))
  names(plotdata)<-c('age','q025','q050','q100','q250','q500','q750','q900','q950','q975')
  
  g <- ggplot(plotdata, aes(x=age)) +
    theme_bw() +
    geom_line(aes(y=q500))+
    geom_ribbon(aes(ymin=q025, ymax=q975), alpha=.5)
  
  if(is.na(lx)==F){
  g <- g +
    geom_point(data=lx, aes(x=age, y=surv), size=.2, col='red')}
  
  }
  
  
  return(g)
  }


####################################################################################
# CHOOSE THE RIGHT NUMBER OF SUBPOPULATIONS BASED ON AIC
####################################################################################


