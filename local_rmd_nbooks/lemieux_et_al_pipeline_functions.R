################################################################################
###                                                                          ###
###   Example code to implement temporal estimation for P. falciparum        ###
###   FUNCTIONS                                                              ###
###                                                                          ###
###   Based on Lemieux et al. 2009. ``Statistical estimation of cell-cycle   ###
###     progression and lineage commitment in Plasmodium falciparum reveals  ###
###     a homogeneous pattern of transcription in ex vivo culture,''         ###
###     PNAS, 106(18): 7559-7564.                                            ###
###                                                                          ###
###   For questions and comments, contact Avi Feller at avi.feller@gmail.com ###
###   Version: 26 May, 2009                                                  ###
###                                                                          ###
################################################################################


################################################################################
###
### Pre-processing functions
###

# Sync reference and test sets
# x is test set
# z is reference set
# Gene description must be "Name"
sync_data <- function(x, z){
  
  # xgenes <- x$?..Name !!NOTE - JR35 CHANGE
  xgenes <- x$Name
  zgenes <- z$Name
  
  these <- intersect(xgenes, zgenes)
  x.o <- x[xgenes %in% these,]
  z.o <- z[zgenes %in% these,]
  
  # xgenes <- as.character(x.o$?..Name)
  # xgenes <- as.character(x.o$..Name) !!NOTE - JR35 CHANGE
  xgenes <- as.character(x.o$Name)
  zgenes <- as.character(z.o$Name)
  
  # make sure genes are in the same order
  route <- rep(NA, nrow(x.o))
  for(i in 1:nrow(x.o)){
    route[i] <- which(xgenes == zgenes[i])
  }
  x.o <- x.o[route,]
  
  list(x.o, z.o)
  
}

# Make data ordinal by column
# 'use.name = T': removes first column
ordinal <- function(x, use.name = T){
  if(use.name == T){
    x <- x[,-c(1)]
  }
  
  these.names <- names(x)
  
  thismat <- NULL
  for(j in 1:ncol(x)){
    thisrank <- rank(x[,j], ties.method="average", na.last="keep")
    thismat <- cbind(thismat, thisrank)
  }
  
  colnames(thismat) <- these.names
  thismat
}

################################################################################
###
### Log-likelihood parameters
###

# Create log-likelihood paramaters
ll.par <- vector("list", 6)
names(ll.par) <- c("minT", "maxT", "window", "grid.unit", "ref", "hourly")

# Create "reference" grid
# Adjust grid parameters
ll.par$window <- 3
ll.par$minT <- 1
ll.par$maxT <- 48
ll.par$grid.unit <- 0.1

ll.par$ref <- seq(ll.par$minT - ll.par$window, ll.par$maxT + ll.par$window,
                  by = ll.par$grid.unit)
ll.par$hourly <- which(ll.par$ref %in% round(ll.par$minT:ll.par$maxT, digits = 1))

ll.par.default <- ll.par


################################################################################
###
### Smoothing Functions
###

# Higher level function for smoothing reference set
smooth.ref <- function(y, ll.par = ll.par.default, method = "spline",
                       span = 0.4, spar = 0.5){
  if(method == "spline"){
    mat <- t(apply(y, 1, function(y) this.smooth.spline(y, spar = spar, ll.par = ll.par)))
    return(mat)
  }
  if(method == "loess"){
    mat <- t(apply(y, 1, function(y) this.loess.smooth(y, span = span, ll.par = ll.par)))
    return(mat)
  }
}


# Wrapper function for loess smoothing

this.loess.smooth <- function(y, ll.par = ll.par.default, span = 0.4){
  ref <- ll.par$ref
  n <- length(y)
  x0 <- c( (-1*n + 1):(2*n) )
  new.y <- c(y, y, y)
  predict(loess(new.y ~ x0, span = span/3), ref)
}

# Wrapper function for spline smoothing
this.smooth.spline <- function(y, ll.par = ll.par.default, spar = 0.5){
  ref <- ll.par$ref
  n <- length(y)
  x0 <- c( (-1*n + 1):(2*n) )
  new.y <- c(y, y, y)
  predict(smooth.spline(x0, new.y, spar = spar), ref)$y
}


# Wrapper function for spline smoothing missing data
smooth.missing <- function(vec){
  vec.o <- na.omit(vec)
  is.missing <- which(is.na(vec))
  n <- length(vec.o)
  x0 <- 1:n
  this.smooth <- smooth.spline(x0, vec.o)
  
  # interpolate missing values
  interp <- rep(NA, sum(is.na(vec)))
  for( i in 1:sum(is.na(vec)) ){
    thisval <-  is.missing[i] - (i-1) - 0.5
    interp[i] <- predict(this.smooth, thisval)$y
  }
  
  vec[is.na(vec)] <- interp
  vec
}




################################################################################
###
### Likelihood Functions
###

# Function to compute likelihood with confidence intervals
full.ll <- function(test, z, sigma = 789.9056, B = 100,
                    sample.rate = 0.50, alpha = 0.05, bootstrap = T, correct.max = T,
                    ll.par = ll.par.default){
  
  # set-up
  ref <- ll.par$ref
  hourly <- ll.par$hourly
  n <- nrow(z)
  
  # get density values
  p <- apply(z, 2, function(y) dnorm(y - test, mean=0, sd=sigma))
  
  # take logs and correct for Inf/-Inf
  p <- log(p)
  p[p == Inf] <- NA
  p[p == -Inf] <- NA
  
  # compute log-likelihoods
  loglik <- colSums(p, na.rm = T)
  loglik.hourly <- loglik[hourly]
  
  # compute MLE
  mle <- ref[which(loglik == max(loglik))]
  
  # subsample confidence intervals
  if(bootstrap == T){
    
    # set-up
    these.ll <- matrix(NA, B, length(ref))
    sample.vec <- 1:n
    draw.rate <- n*sample.rate
    these.max <- rep(NA, B)
    
    # subsample
    for(i in 1:B){
      # compute log-likelihood for subsample of genes
      temp.ll <- colSums(p[ sample(sample.vec, draw.rate), ] , na.rm = T)
      
      # compute maxima in each subsamples
      these.max[i] <- ref[ which(temp.ll == max(temp.ll)) ]
    }
    
    # compute maxima in each subsamples
    #these.max <- apply(these.ll, 1,  function(y)  ref[ which(y == max(y)) ])
    
    # correct for CIs that wrap around; can be disabled
    if(correct.max){
      these.max[ these.max > ll.par$maxT ] <- these.max[ these.max > ll.par$maxT] -
        ll.par$maxT
    }
    
    # store CI
    this.ci <- quantile(these.max, probs = c(alpha/2, 1 - alpha/2) )
    
  }else{  # else, store NA
    this.ci <- NA
  }
  
  # output results
  output <- list(loglik, loglik.hourly, mle, this.ci)
  names(output) <- c("loglik", "loglik.hourly", "mle", "ci")
  
  output
}



## Wrapper function for full.ll
compute.ll <- function(x, z, sigma = 789.9056, B = 100,
                       sample.rate = 0.50, alpha = 0.05, bootstrap = T, correct.max = T,
                       ll.par = ll.par.default, ...){
  
  apply(x, 2, function(y){    full.ll(y,
                                      z,
                                      sigma = sigma,
                                      B = B,
                                      sample.rate = sample.rate,
                                      alpha = alpha,
                                      bootstrap = bootstrap,
                                      correct.max = correct.max,
                                      ll.par = ll.par, ...) })
}


## Convenience functions to extract parameters from ll
# Log-likelihoods
loglik <- function(ll){
  t(sapply(ll, function(y) y$loglik.hourly))
}

# Maximum Likelihood Estimates
mle <- function(ll){
  sapply(ll, function(y) y$mle)
}

# Confidence Intervals
ci <- function(ll){
  t(sapply(ll, function(y) y$ci))
}



################################################################################
###
### Plotting Functions
###

# Takes ll file and plots MLEs and CIs (if already computed)
# Plots rainbow colors by default
plot.mle <- function(ll, main = "MLEs of Parasite Age",
                     xlab = "HPI", ylab = "Samples", col = NULL, ...){
  
  # set-up
  this.mle <- mle(ll)
  this.ci <- ci(ll)
  maxT <- ncol( loglik(ll) )
  n <- length(this.mle)
  
  # change color of points
  if( is.null(col) ){
    this.col <- rainbow(n)
  }else{
    this.col <- col
  }
  
  # base plot
  plot(this.mle, 1:n, pch = 16, xlim = c(0, maxT), type = "n",
       yaxt = "n", xlab = xlab, ylab = ylab,
       main = main, ...)
  
  # plot CIs if not NA
  if(!all(is.na(this.ci))){
    arrows(this.ci[,1], 1:n, this.ci[,2], 1:n,
           length = 0.05, # width of the arrowhead
           angle = 90, # angle of the arrowhead
           code = 3 # arrowhead in both ends
    )
  }
  
  # plot MLEs
  points( this.mle, 1:length(this.mle), pch = 16, col = this.col )
}


plot.ll <- function(ll, main="log-Likelihood", col=NULL,
                    xlab="Time (hr)", ylab="log-Likelihood", ylim=NULL, ...){
  
  this.ll <- loglik(ll)
  
  x <- 1:ncol(this.ll)
  
  if(is.null(ylim)){
    this.ylim <- range(this.ll)
  }else{
    this.ylim <- ylim
  }
  plot(x, this.ll[1,], ylim=this.ylim, type="n", xlab=xlab, ylab=ylab,
       main=main, ...)
  
  if(is.null(col)){
    for(i in 1:nrow(this.ll)){
      lines(x, this.ll[i,], col=rainbow(nrow(this.ll))[i])
    }
  }else{
    for(i in 1:nrow(this.ll)){
      lines(x, this.ll[i,], col=col[i])
    }
  }
}


