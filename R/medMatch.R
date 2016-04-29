################
#' Calculate propensity score for matching.
#' 
#' Computes the propensity score for treatment based on covariates. Should not generally be called directly.
#' 
#' The propensity scores are computed by fitting a logistic regression model to the binary outcome \code{Y},
#'  with linear terms for each of the baseline covariates \code{W1, W2, ...}
#' 
#' @param dat A data frame containing columns for the treatment (\code{Z}), outcome (\code{Y}), 
#'  potential mediator (\code{S} and \code{S.cat}), and baseline covariates (\code{W1, W2, ...})
#' @return A vector of propensity scores of the same length as \code{Z}.
#' 
propScores <- function(dat) {
  w.ind <- !(colnames(dat) %in% c("Z","Y","S","S.cat"))
  fmla <- as.formula(paste0("Z~",paste0(colnames(dat[,w.ind]),collapse="+")))
  fitted(glm(fmla,family=binomial,data=dat))
}

################
#' Calculates matching-based bounds on the natural direct effect.
#' 
#' Internal function for calculating matching-based bounds for the natural direct effect.
#' 
#' This function is called by the bootstrap procedure in the main \code{\link{DeltaBounds}} function.
#' 
#' @param DAT A data frame containing columns for the treatment (\code{Z}), outcome (\code{Y}), 
#'  potential mediator (\code{S} and \code{S.cat}), and baseline covariates (\code{W1, W2, ...})
#' @param inds A vector of indicators for the bootstrap.
#' @param Gamma The sensitivity parameter controlling the odds ratio between the true and estimated odds of a discordant pair favoring the treatment.
#' @param Psi The sensitivity parameter controlling the odds ratio between the true and estimated odds of a pair being discordant.
#' @param effect The causal effect to be estimated. Can be \code{ATT} or \code{ATC}.
#' @return A vector of three elements containing: \enumerate{
#' \item The estimated lower bound
#' \item The estimated value
#' \item The estimated upper bound 
#' }
#' @export
getBounds <- function(DAT,inds,Gamma=1,Psi=Gamma,effect="ATT") {
  dat <- DAT[inds,]
  z <- dat$Z
  ##  w <- dat$W
  y <- dat$Y
  s.cat <- dat$S.cat
  
  ps <- rep(NA,nrow(dat))
  
  ## Compute the propensity scores
  for(s in 1:max(s.cat)) {
    ind <- which(s.cat==s)
    ps[ind] <- propScores(dat[ind,])
  }
  
  matched <- Matching::Matchby(y,Tr=z,by=s.cat,X=ps,est=effect,replace=TRUE,print.level=0,ties=FALSE,M=1) ## With or without replacement?
  MD <- matched$mdata
  
  IT <- matched$index.treated
  IC <- matched$index.control
  PAIRS <- cbind(IT,IC,y[IT],y[IC])
  
  #### Calculate the differences and perform the statistical test ####
  pairdiffs <- PAIRS[,3] - PAIRS[,4] ## Treated minus control differences
  if(all(pairdiffs==0)) {
    DeltaN.L <- DeltaN <- DeltaN.U <- NA
  } else {
    n <- length(pairdiffs)
    n.diffs <- sum(pairdiffs!=0)
    nT.diffs <- sum(pairdiffs==-1) ## Number favoring the treatment
    nC.diffs <- n.diffs - nT.diffs ## Number favoring the control
    
    pi.D <- n.diffs/n
    pi.D1 <- nT.diffs/n.diffs
    
    pL.D <- pi.D / (pi.D + Gamma*(1 - pi.D))
    pU.D <- pi.D / (pi.D + (1-pi.D)/Gamma)
    
    pL.D1 <- pi.D1 / (pi.D1 + Psi*(1 - pi.D1))
    pU.D1 <- pi.D1 / (pi.D1 + (1-pi.D1)/Psi)
    
    DeltaN.L <- -ifelse((2*pU.D1-1)>0,pU.D*(2*pU.D1 - 1),pL.D*(2*pU.D1-1))
    DeltaN.U <- -ifelse((2*pL.D1-1)>0,pL.D*(2*pL.D1 - 1),pU.D*(2*pL.D1-1))
    
    DeltaN <- mean(pairdiffs)
    #  DeltaN <- pi.D*(2*pi.D1-1)
  }
  
  c(DeltaN.L,DeltaN,DeltaN.U)
  #c(LB,ODDS,UB)
}

################
#' Calculates matching-based bounds on the natural direct effects.
#' 
#' This function estimates matching-based bounds on the two natural direct effects 
#' 
#' This function computes uncetainty intervals for natural direct effects which incorporate uncertainty in the point estimate
#' due to both possible unmeasured confounding and statistical uncertainty in estimation due to the matching procedure.
#' Currently, percentile-based bootstrap uncertainty intervals are returned as they appear to yield somewhat better coverage than
#' alternatives.
#' 
#' @param Z A vector of treatment indicators.
#' @param Y A vector of (binary) outcomes.
#' @param S A vector of biomarker values, i.e., potential mediators
#' @param W A vector or matrix of covariates to be used in the matching
#' @param S.cat (Optional) vector containing a coarsened version of \code{S}. There should be at least one treated and one untreated subject within each level of \code{S.cat}.
#'  By default (\code{S.cat=NULL}), \code{S} is split into quartiles.
#' @param Gamma The sensitivity parameter controlling the odds ratio between the true and estimated odds of a discordant pair favoring the treatment. Default is 1 (corresponding to no unmeasured confonding), values >= 1 are valid.
#' @param Psi The sensitivity parameter controlling the odds ratio between the true and estimated odds of a pair being discordant. Default is 1 (corresponding to no unmeasured confonding), values >= 1 are valid.
#' @param NDE Specifies which natural direct effect (NDE) to calculate. Valid options are \code{"Delta1"} (the default)
#'  and \code{"Delta0"}, corresponding to \eqn{E(Y^{1S(1)} - Y^{0S(1)})} and \eqn{E(Y^{1S(0)} - Y^{0S(0)})} respectively.
#' @param R Number of bootstrap samples used to estimate uncertainty intervals. Defaults to \code{500}.
#' @param conf Confidence level for the bootstrapped uncertainty intervals. Defaults to \code{0.95}.
#' @return A list containing:
#'  \describe{
#'  \item{\code{delta}}{The point estimate for the NDE.}
#'  \item{\code{CI}}{The confidence interval for the NDE, reflecting only stochastic uncertainty, ignoring uncertainty due to unmeasured confounding.}
#'  \item{\code{II}}{The ignorance interval for the NDE, reflecting only uncertainty due to unmeasured confounding as reflected by \code{Gamma} and \code{Psi}.}
#'  \item{\code{UI}}{The uncertainty interval for the NDE, reflecting both sources of uncertainty.}
#' }
#' @export
DeltaBounds <- function(Z,Y,S,W,S.cat=NULL,Gamma=1,Psi=1,NDE="Delta1",R=500,conf=0.95) {
  
  if(is.null(S.cat)) { 
    cut.pts <- c(-Inf,quantile(S,c(0.25,0.5,0.75)),Inf)
    S.cat <- cut(S,cut.pts,labels=FALSE)
  }
  sapply(levels(S.cat),function(l) {
    if(length(unique(S.cat[S.cat==l]))<=1) { stop('Coarsening of S resulted in at least one category with no variability in treatment assignment. 
                                                  If S.cat was not specified explicitly, please supply one.')}
  })
  dat <- data.frame(Z,Y,S,W,S.cat)
  if(NDE=="Delta1") { effect <- "ATT" }
  else {
    if(NDE=="Delta0") { effect <- "ATC" }
    else{ stop("Not a recognized type of Natural Direct effect. Valid types are 'Delta1' and 'Delta0'")}
  }
  if(Gamma < 1 | Psi < 1) { stop("Only sensitivity parameter values >= 1 are permitted.")}
  
  BOOT.W <- boot::boot(dat,
                       getBounds,
                       R=R,
                       Gamma=Gamma,
                       Psi=Psi,
                       effect=effect)
  
  BCI.LB <- boot::boot.ci(BOOT.W,
                          conf=conf,
                          type=c("norm","basic","perc"),
                          index=1)
  BCI.DELTA <- boot::boot.ci(BOOT.W,
                             conf=conf,
                             type=c("norm","basic","perc"),
                             index=2)
  BCI.UB <- boot::boot.ci(BOOT.W,
                          conf=conf,
                          type=c("norm","basic","perc"),
                          index=3)
  
  DELTA.L <- BOOT.W$t0[1]
  DELTA.MN <- BOOT.W$t0[2]
  DELTA.U <- BOOT.W$t0[3]

  return( list( delta=DELTA.MN,
                CI=BCI.DELTA$perc[4:5],
                II=c(DELTA.L,DELTA.U),
                UI=c(BCI.LB$perc[4],BCI.UB$perc[5]) )
  )
  
}

plotBounds <- function(Delta, CI, UI, Gamma.vals ) {
  ## In progress; need to make more general for Delta.1 or Delta.0
  
  df <-data.frame( Delta, CI = CI, UI = UI, Gamma.vals = factor(Gamma.vals) )
  
    g <- ggplot( df, aes( x = Gamma.vals, y = Delta ))
    g <- g + geom_point( size=5, shape = 18 ) +
          geom_errorbar( aes( ymin = UI.1, ymax = UI.2 ), width = 0.3, colour = gray(0.5)) +
          geom_linerange( aes( ymin = CI.1, ymax = CI.2 ), size=1.2) +
          geom_hline( y = 0 ) + 
          xlab( expression(Gamma == Psi) ) #+
    #ylab( expression(Delta^{N}*(1)) ) +
    #ggtitle( "Matching every treated subject" )
  
  return(g)
  

}
  
