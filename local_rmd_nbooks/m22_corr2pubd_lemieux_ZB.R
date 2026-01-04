##mle(ll)

################################################################################
###                                                                          ###
###   Example code to implement temporal estimation for P. falciparum        ###
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
### Source Functions
###

# setwd("/MLE/") !!NOTE - JR35 CHANGE
setwd("/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/")
source("mali22_code_lnk/local_rmd_nbooks/lemieux_et_al_pipeline_functions.R")


################################################################################
###
### Load Data
###


## Load reference set
# We use a cleaned version of the Bozdech et al. (2003) overview data set with
# NULL genes removed and duplicates averaged, yielding n = 3532 genes.
# Can also use 3D7 and Dd2 reference sets
z <- read.csv("data/raw/Pf/andrade_lemiux_MLE/bozdech_Hb3_clean.csv", as.is = T)
z

## Clean names
Pf3D7_genes_plasmoDB_w_old_ids <- Pf3D7_genes_plasmoDB %>% separate_wider_delim(`Previous ID(s)`, delim = ";", names = c("old_id1", "old_id2", "old_id3", "old_id4", "old_id5", "old_id6", "old_id7", "old_id8"), too_few = "align_start", too_many = "debug")

## Load test set
# We use the sorbitol replicates from Le Roch et al. (2003).
# Substitute other test set here
#x <- read.csv("winzeler_stage_expression.csv", as.is=T, sep = ';')

##Load the input data file
x <- read.csv("data/raw/Pf/andrade_lemiux_MLE/wgenes.csv", as.is=T, sep = ';')
x <- x[,c(1,2:13)]

################################################################################
###
### Pre-process
###


## Sync data
# Ensures that data sets contain the same genes and are in the same order
# Note that first column must be gene description titled "Name"
data <- sync_data(x, z)
x <- data[[1]]
z <- data[[2]]


## Store names
#genenames <- x$?..Name

#genenames <- x$Name
## Make data ordinal
# This is necessary for comparisons between one- and two-color arrays
# Returns a matrix and if `use.name = T', removes first "Name" column
x <- ordinal(x, use.name = T)
z <- ordinal(z, use.name = T)


## If same arry type, remove first column and make into a matrix
#x <- as.matrix(x[,-1])
#z <- as.matrix(z[,-1])


## Impute missing values
# For later interpretability, we impute the missing time points (23 and 29 HPI)
# for the Hb3 reference set. This method interpolates missing values using
# smoothing splines, though other approaches are possible.
z.na <- cbind(z[,1:22], rep(NA, nrow(z)), z[,23:27],
              rep(NA, nrow(z)), z[,28:46])
z <- t(apply(z.na, 1, smooth.missing))


## Load sigma_epsilon
# Estimate of background noise obtained by comparing Hb3 and 3D7 samples
sigma.epsilon <- 789.9056


################################################################################
###
### Get smooth reference data
###


## Check default log-likelihood parameters, stored in ll.par


## Smooth reference data set
# This smooths the time course for each gene in the reference set using a
# smoothing spline. There is also a loess implementation.
# `spar' is the parameter of smoothing
z.smooth <- smooth.ref(z, method = "spline", spar = 0.5)




## Estimate residual of smoothing
z.smooth.hourly <- z.smooth[,ll.par$hourly]
sigma.eta <- mean(sd(z - z.smooth.hourly, na.rm = T), na.rm=T)


## Calculate new sigma based on sigma_eta and sigma_epsilon
new.sigma <- sqrt(sigma.eta^2 + sigma.epsilon^2)


################################################################################
###
### Get Log-Likelihood, MLE, and Confidence Intervals
###


#-----------------------
# Main function options
#-----------------------
# B: number of iterations; 100 is default; 500 is better
# sample.rate: rate of subsampling; 0.50 is default
# alpha: for confidence intervals; alpha = 0.05 is default


ll <- compute.ll(x = x, z = z.smooth, sigma = new.sigma,
                 bootstrap = T,
                 B = 100,
                 sample.rate = 0.50)



# mle(ll) --> extracts MLEs
# loglik(ll) --> extracts log-likelihoods
# ci(ll) --> extracts confidence intervals


## To obtain only the log-likelihoods without bootstrapped confidence intervals
ll <- compute.ll(x = x, z = z.smooth, sigma = new.sigma,
                 bootstrap = F)


###############################################################################
###
### Plot Logliks and MLEs
###

## Plot the log-likelihood curves
plot.ll(ll)

## Plot the maximum likelihood estimates with 95% confidence intervals
plot.mle(ll)


