# Dataset: one vector of (residual) age-at-death and one boolean vector with censorship.
# The algorithm is an extension to randomly right-censored data of the one proposed by Jiang & Kececioglu (1992)
# The module contains 4 main functions:
# SimWeibull(...)
# FitWeibull(...)
# PlotWeibull(...)
# CIWeibull(...)
# ChooseNWeibull(...)

##############################################################################################################
# CENSORED MIXTURE WEIBULL SIMULATOR
##############################################################################################################
# Simulate an N-population mixed-weibull dataset:

SimWeibull <- function(nSamp=1000, percentcensor=.2, theta){
  
# Fix the missing degree of freedom on theta through the last weight in case of error:
weights <- theta[seq(1,length(theta),3)]
last_weight <- 1-sum(weights[-length(weights)])
theta[length(theta)-2] <- last_weight

nPop <- length(theta)/3
aad <- vector()
for(p in 1:nPop){
  w <- theta[((p-1)*3+1)] ; a <- theta[((p-1)*3+2)] ; b <- theta[((p-1)*3+3)]
  nspop <- round(nSamp*w)
  aad <- c(aad, rweibull(nspop, scale=a, shape=b))
}
censor <- sample(c(0,1), nSamp, replace=T, prob=c(percentcensor, (1-percentcensor)))
for(a in 1:nSamp){
  if(censor[a]==0){
    aad[a]<-runif(1,0,aad[a])}}

return(data.frame(aad, censor))}

####################################################################################
####################################################################################

sim <- SimWeibull(nSamp=5000, percentcensor = .1, theta = c(.6, 170, 4, .4, 57, 9))
aad <- as.numeric(sim[,1])
censor <- as.numeric(sim[,2])

####################################################################################
# MAIN EXPECTATION-MAXIMISATION FITTING FUNCTION
####################################################################################

FitWeibull <- function(aad, censor, n=2, weights=NA, thres=1e-4, NR_thres=1e-6, verbose=F, graphical=FALSE){
  # Structure of the theta vector: {w1, a1, b1, w2, a2, b2... wi, ai, bi} 
  require(LearnBayes)
  start_time <- proc.time()
  
  censor <- censor[order(aad)]
  aad <- sort(aad)
  
  
  # the simple, 2-parameter Weibull PDF function:
  f1 <- function(x,a,b){return((b/a)*((x/a)^(b-1)) * exp(-(x/a)^b))}
  # the 2-sub-population, 5-parameter mixed PDF:
  f <- function(x, theta){theta[1]*f1(x, a=theta[2], b=theta[3]) + (1-theta[1])*f1(x, a=theta[4], b=theta[5])}
  # the simple, 2-parameter Weibull CDF function:
  cdf <- function(x,a,b){return(1 - exp(-(x/a)^b))}
  # the simple, 2-parameter Weibull tail function:
  tf <- function(x,a,b){return(exp(-(x/a)^b))}
  # the 2-sub-population, 5-parameter mixed tail function:
  R <- function(x, theta){return(1 - (theta[1]*cdf(x, a=theta[2], b=theta[3]) + (1-theta[1])*cdf(x, a=theta[4], b=theta[5])))}
  
####################################################################################
# N-population Weibull log-likelihood function
LnL <- function(aad, censor, theta){
	# theta is a vector of model parameters {w1, a1, b1, w2, a2, b2... wi, ai, bi}
	# a are scales, b are shapes (opposite from other R functions but reflects the literature)
  nPop <- length(theta)/3
  
  LogPrTj <- function(t, c, theta){
  # Summed probability of failure time Tj over all components 'i'
  # the simple, 2-parameter Weibull PDF function for a single component:
  total <- 0
  if(c == 1) {
    for(k in 1:nPop){
      wi <- theta[((k-1)*3+1)]
      ai <- theta[((k-1)*3+2)]
      bi <- theta[((k-1)*3+3)]
      total <- total + wi*f1(x=t, a=ai, b=bi)}
    return(log(total))
  } else {
    for(k in 1:nPop){
      wi <- theta[((k-1)*3+1)]
      ai <- theta[((k-1)*3+2)]
      bi <- theta[((k-1)*3+3)]
      total <- total + wi*tf(x=t, a=ai, b=bi)}
    return(log(total))}
  }

  SumLogPrTj <- 0
  for(j in 1:length(aad)){
    SumLogPrTj <- SumLogPrTj + LogPrTj(aad[j], censor[j], theta)}
  
  n <- length(censor)
  r <- sum(censor)
  # Full likelihood function:
  L <- (lfactorial(n) - lfactorial(n-r)) + SumLogPrTj
  
  return(L)}

##############################################################################################################
# For each time "tj", determine the probability of belonging to subp. "i" given a parameter estimate theta:
PRiTj <- function(t, c, i, theta){
  w <- theta[((i-1)*3+1)] ; a <- theta[((i-1)*3+2)] ; b <- theta[((i-1)*3+3)]
  if(c == 1) {
    total <- 0
    for(k in 1:(length(theta)/3)){
      wi <- theta[((k-1)*3+1)]
      ai <- theta[((k-1)*3+2)]
      bi <- theta[((k-1)*3+3)]
      total <- total + wi*f1(x=t, a=ai, b=bi)}
    return(w * (f1(x=t, a=a, b=b)/total))
  } else {
    total <- 0
    for(k in 1:(length(theta)/3)){
      wi <- theta[((k-1)*3+1)]
      ai <- theta[((k-1)*3+2)]
      bi <- theta[((k-1)*3+3)]
      total <- total + wi*tf(x=t, a=ai, b=bi)}
    return(w * (tf(x=t, a=a, b=b)/total))   
  }}

##############################################################################################################
## FUNCTIONS TO USE DURING THE ITERATION PHASES

# Describe Bi given PRiTj:
G <- function(b, i, aad, censor, theta) {
  
  G1 <- function(b, i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + (PRiTj(aad[j], censor[j], i, theta) * log(aad[j]) * aad[j]^b)}
    return(total)}
  
  G2 <- function(i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta)}
    return(total)}

  G3 <- function(b, i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + (PRiTj(aad[j], censor[j], i, theta) * aad[j]^b)}
    return(total)}  

  G4 <- function(i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta) * log(aad[j])}
    return(total)}  

  G5 <- function(b, i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta)}
    return((1/b) * total)}

  G <- (-1 * G1(b, i, aad, censor, theta) * G2(i, aad, censor, theta)) / G3(b, i, aad, censor, theta) + G4(i, aad, censor, theta) + G5(b, i, aad, censor, theta)

  return(G)}

# Derivative of G to be used in the Newton iteration:
g <- function(b, i, aad, censor, theta) {
  
  g1 <- function(b, i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta)}
    total <- (-1/(b^2)) * total
    return(total)}
  
  g2 <- function(i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta)}
    return(total)}
  
  U <- v <-  function(b, i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta) * log(aad[j]) * aad[j]^b }
    return(total)}
  
  u <- function(b, i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta) * (log(aad[j]))^2 * aad[j]^b }
    return(total)}

  V <- function(b, i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta) * aad[j]^b }
    return(total)}
  
  g <- g1(b, i, aad, censor, theta) - g2(i, aad, censor, theta) * ((u(b, i, aad, censor, theta) * V(b, i, aad, censor, theta) - U(b, i, aad, censor, theta) * v(b, i, aad, censor, theta)) / V(b, i, aad, censor, theta)^2)
  
  return(g)}

# Newton-Raphson iteration to find the root of G, starting from a guess "B0":
# This converges very nicely provided the starting estimates are not too bad.
NR <- function(b, i, add, censor, theta, thres = NR_thres, verbose=F){
  B0 <- b
  B1 <- B0 - (G(B0, i, aad, censor, theta)/g(B0, i, aad, censor, theta))
  while(abs(B0-B1) >= thres) {
    B0 <- B1
    B1 <- B0 - (G(B0, i, aad, censor, theta)/g(B0, i, aad, censor, theta))
    if(verbose==T){
    print(B1)}
  }
  return(B1)
}

##############################################################################################################
# Get an estimate of Ai given an updated estimate of Bi

Ai <- function(b, i, add, censor, theta){
  
  A1 <- function(b, i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + (PRiTj(aad[j], censor[j], i, theta) * aad[j]^b)}
    return(total)}
  
  A2 <- function(i, aad, censor, theta) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + PRiTj(aad[j], censor[j], i, theta)}
    return(total)}
  
  return((A1(b, i, aad, censor, theta) / A2(i, aad, censor, theta))^(1/b))
  
}

##############################################################################################################
# Get an estimate of Wi given a updated estimates of Ai and Bi sotred in theta

Wi <- function(i, aad, censor, theta){
  total <- 0
  for(j in 1:length(aad)){
    total <- total + PRiTj(aad[j], censor[j], i, theta)}
  return(total/length(aad))}

##############################################################################################################
## FUNCTIONS TO USE DURING THE INITIAL PHASE (Random probabilities and no parameter estimates)

# Describe Bi given the random belongings RB:
Init_G <- function(b, i, aad, censor, RB) {
  
  G1 <- function(b, i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + (RB[j,i] * log(aad[j]) * aad[j]^b)}
    return(total)}
  
  G2 <- function(i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i]}
    return(total)}
  
  G3 <- function(b, i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + (RB[j,i] * aad[j]^b)}
    return(total)}  
  
  G4 <- function(i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i] * log(aad[j])}
    return(total)}  
  
  G5 <- function(b, i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i]}
    return((1/b) * total)}
  
  G <- (-1 * G1(b, i, aad, censor, RB) * G2(i, aad, censor, RB)) / G3(b, i, aad, censor, RB) + G4(i, aad, censor, RB) + G5(b, i, aad, censor, RB)
  
  return(G)}

# Derivative of G to be used in the Newton iteration:
Init_g <- function(b, i, aad, censor, RB) {
  
  g1 <- function(b, i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i]}
    total <- (-1/(b^2)) * total
    return(total)}
  
  g2 <- function(i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i]}
    return(total)}
  
  U <- v <-  function(b, i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i] * log(aad[j]) * aad[j]^b }
    return(total)}
  
  u <- function(b, i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i] * (log(aad[j]))^2 * aad[j]^b }
    return(total)}
  
  V <- function(b, i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i] * aad[j]^b }
    return(total)}
  
  g <- g1(b, i, aad, censor, RB) - g2(i, aad, censor, RB) * ((u(b, i, aad, censor, RB) * V(b, i, aad, censor, RB) - U(b, i, aad, censor, RB) * v(b, i, aad, censor, RB)) / V(b, i, aad, censor, RB)^2)
  
  return(g)}

# Newton-Raphson iteration to find the root of G, starting from a guess "B0":
# This converges very nicely provided the starting estimates are not too bad.
Init_NR <- function(b, i, add, censor, RB, thres = NR_thres, verbose=F){
  B0 <- b
  B1 <- B0 - (Init_G(B0, i, aad, censor, RB)/Init_g(B0, i, aad, censor, RB))
  while(abs(B0-B1) >= thres) {
    B0 <- B1
    B1 <- B0 - (Init_G(B0, i, aad, censor, RB)/Init_g(B0, i, aad, censor, RB))
    if(verbose==T){
    print(B1)}
  }
  return(B1)
}

##############################################################################################################
# Get an estimate of Ai given an updated estimate of Bi

Init_Ai <- function(b, i, add, censor, RB){
  
  A1 <- function(b, i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + (RB[j,i] * aad[j]^b)}
    return(total)}
  
  A2 <- function(i, aad, censor, RB) {
    total <- 0
    for(j in 1:length(aad)){
      total <- total + RB[j,i]}
    return(total)}
  
  return((A1(b, i, aad, censor, RB) / A2(i, aad, censor, RB))^(1/b))
  
}

##############################################################################################################
# Get an estimate of Wi given a updated estimates of Ai and Bi sotred in theta

Init_Wi <- function(i, aad, censor, RB){
  total <- 0
  for(j in 1:length(aad)){
    total <- total + RB[j,i]}
  return(total/length(aad))}

####################################################################################
# EMPIRICAL QUANTILE FUNCTION GIVEN A BELONGING OBJECT
####################################################################################

#This ignores the censoring...
emp_quant <- function(p, i, belonging, censor=censor, aad=aad){
  belonging <- belonging[censor==1,]
  normalise_belonging <- function(x) x/sum(x)
  belonging <- apply(belonging, 2, normalise_belonging)
  with_aad <- data.frame(aad=aad[censor==1], belong = belonging[,i])
  with_aad <- with_aad[order(with_aad$aad),]
  cumbelong <- apply(with_aad, 2, cumsum)
  cumbelong[,1]<-sort(aad[censor==1])
  return(cumbelong[cumbelong[,2]>=p,][1,1])}

# This estimator is taken from Marks 2005, Estimation of Weibull parameters from common percentiles
b_quantile_estimate <- function(i, belonging, aad=aad, censor=censor, level=.1){
  qtL <- level
  qtH = 1 - qtL
  L <- emp_quant(qtL, i, belonging=belonging, aad=aad, censor=censor)
  H <- emp_quant(qtH, i, belonging=belonging, aad=aad, censor=censor)
  b <-  log(log(qtL)/log(qtH)) * log(H/L)^(-1)
  return(b)
}


##############################################################################################################
# Initialise an Lx plot if desired
if(graphical == TRUE){
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
initial_belonging <- rdirichlet(length(aad), rep(.1, n))
initial_belonging <- initial_belonging[order(initial_belonging[,1]),]
} else { initial_belonging <- weights }
LogLikelihood <- vector()

# Make an empty vector Theta that will hold the initial estimates:
Init_Theta <- vector(length=3*n)

for(i in 1:n){
  b <- b_quantile_estimate(i, belonging=initial_belonging, aad=aad, censor=censor, level=.1)
  # Get a first estimate of Bi:
  Init_Theta[(i-1)*3+3] <- Init_NR(b, i, add=aad, censor=censor, RB=initial_belonging, thres = thres)
  # Get the corresponding Ai:
  Init_Theta[(i-1)*3+2] <- Init_Ai(b, i, add=aad, censor=censor, RB=initial_belonging)
  # Get the corresponding weight:
  Init_Theta[(i-1)*3+1] <- Init_Wi(i, aad=aad, censor=censor, RB=initial_belonging)
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

  for(i in 1:n){
    b <- b_quantile_estimate(i, belonging=current_belonging, aad=aad, censor=censor, level=.1)
    # Get a first estimate of Bi:
    this_theta[(i-1)*3+3] <- NR(b, i, add, censor, theta=previous_theta, thres = thres)
    # Get the corresponding Ai:
    this_theta[(i-1)*3+2] <- Ai(b, i, add, censor, theta=previous_theta)
    # Get the corresponding weight:
    this_theta[(i-1)*3+1] <- Wi(i, aad, censor, theta=previous_theta) 
  }

  LogLikelihood <- append(LogLikelihood, LnL(aad=aad, censor=censor, theta=this_theta))
  thetas <- rbind(thetas, this_theta)
  if(verbose==T) print(thetas)
  
  # Evaluate the convergence criterion:
  if(sum(abs(previous_theta - this_theta) < thres)){
    stable_steps <- stable_steps +1 } else { stable_steps <- 0 }
  if(stable_steps >= 3) converged <- TRUE
  
  if(graphical == TRUE){
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
  plot(lx$age, lx$surv, cex=.05)}
  formula <- paste(rep('%s * exp(-(x/%s)^%s) + ', length(this_theta)/3), sep='', collapse='')
  formula <- do.call(sprintf, c(list(formula), round(this_theta, 3)))
  formula <- substr(formula,1,nchar(formula)-3)
  funform <- function(x) eval(parse(text=formula))
  plot(lx$age, lx$surv, cex=.05)
  curve(funform, from=min(aad), to=max(aad), add=add, col='red')}

####################################################################################
# PARAMETRIC BOOTSTRAP ROUTINE
####################################################################################

CIWeibull <- function(method='nonparametric', theta=NA, n=2, aad=NA, censor=NA, nBs = 100, verbose=T, thres=1e-4, NR_thres=1e-6){
  # Type is (1) "parametric" or (2) "nonparametric"
  # For non-parametric bootstrapping, theta is not needed.
  # For parametric bootstrapping, n is not needed
  # For parametric bootstrapping, aad and censor are used to extract dataset shape
  
  # Some controls on passed arguments
  # XXXXX
  
  if(is.na(theta)==FALSE){
  n <- length(theta)/3}
  
  fits <- list()
  
  if(method=='nonparametric'){
    for(i in 1:nBs){
      fit <- NA
      while(is.na(fit)==T){
      index <- sort(sample(1:length(aad), length(aad), replace=T))
      this_aad <- aad[index]
      this_censor <- censor[index]
      fit <- FitWeibull(this_aad, this_censor, n=n, weights=NA, thres=thres, NR_thres=NR_thres)}
      fits[[i]] <- fit}
    
  } else if(method=='parametric'){
    for(i in 1:nBs){
      fit <- NA
      while(is.na(fit)==T){
      simdata <- SimWeibull(nSamp=length(aad), percentcensor=((length(censor)-sum(censor))/length(censor)), theta = theta)
      fit <- FitWeibull(aad=as.numeric(simdata$aad), censor=as.numeric(simdata$censor), n=n, weights=NA, thres=thres, NR_thres=NR_thres)}
      fits[[i]] <- fit}
    
  } else { print(paste(type, ': invalid bootstrapping type (parametric or nonparametric)'))}

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
  
 return(output)}


require(ggplot2)
pdfplot <- as.data.frame(cbind(age=1:125, t(ci$PDF)))
names(pdfplot)<-c('age','q025','q050','q100','q250','q500','q750','q900','q950','q975')
ggplot(pdfplot, aes(x=age)) +
  theme_bw() +
  geom_line(aes(y=q500))+
  geom_ribbon(aes(ymin=q025, ymax=q975), alpha=.5)

cdfplot <- as.data.frame(cbind(age=1:125, t(ci$CDF)))
names(cdfplot)<-c('age','q025','q050','q100','q250','q500','q750','q900','q950','q975')
ggplot(cdfplot, aes(x=age)) +
  theme_bw() +
  geom_line(aes(y=q500))+
  geom_ribbon(aes(ymin=q025, ymax=q975), alpha=.5)

####################################################################################
# CHOOSE THE RIGHT NUMBER OF SUBPOPULATIONS BASED ON AIC
####################################################################################


