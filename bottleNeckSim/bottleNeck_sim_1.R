####### initialise packages and python environmnet
library(igraph)
library(reticulate)
use_condaenv("simupop-env", required = TRUE)

### set working dir
#setwd("~/bottleNeckSim")



## load simulation functions
source("simFunctions.R")


####initialise haplotypes and emperical counts
## this is the sampling haplotype freq 
## will sample using a Dirichlet dist and create initial mtDNA freq inputs
#hapFreq <-  c(c(2, 1, 1, 2, 2, 3, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1) + 0.05, rep(0.05,6)) #historic + additional 6
## new corrected frequency:
hapFreq <-c(c(rep(1, 7), rep(2, 6), rep(3, 2)) + 0.05, rep(0.05,6)) #historic + additional 6



runs.vec <- c(5000,7) ## number of initial mtDNA files and runs per initial mtDNA sequence.
## total number of runs will be:
# runs.vec[1] * runs.vec[2]
set.seed(123456)#sets an initial seed to generate the seed matrix
seedMatrix <- generate_sim_seeds(runs.vec = runs.vec) ## seed matrix, first col are initial mtDNA seeds, cols 2:end are pop seeds
run_sim(runs.vec = runs.vec, Ne1 = 10, seedMatx = F, cleanLogs = F)



