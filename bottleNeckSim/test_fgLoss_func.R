
hap.dat <- read.table('/Users/adamglucksman/Local PhD/bottleNeckSim/results/127109/127109_haploiso_75.txt', header = TRUE, sep = " ")

hap.dat %>%
  mutate(feeding_ground. = if_else(
    feeding_ground. == 0,
    sample(c(1, 2), size = n(), replace = TRUE)[feeding_ground. == 0],
    feeding_ground.
  ))

hap.dat %>%
  rowwise() %>%
  mutate(feeding_ground. = if_else(
    feeding_ground. == 0,
    sample(c(1, 2), 1),
    feeding_ground.
  )) %>%
  ungroup()


calc.hap.div_test('/Users/adamglucksman/Local PhD/bottleNeckSim/results/127109/127109_haploiso_end.txt')
calc.hap.div('/Users/adamglucksman/Local PhD/bottleNeckSim/results/127109/127109_haploiso_end.txt')

library(dplyr)

################################ 
#### calc hap diversity and S.D. 
################################

## calculates haplotype diversity
## input: filepath of generation summary output
## input files should be in the form of {SEED}_haploiso_{GEN}.txt, i.e. 587_haploiso_75.txt
## returns total num of haplotypes, diversity, and variance 

calc.hap.div_test <- function(hap_file) {
  
  # # Initialize results table
  # results <- data.frame(seed = character(),
  #                       hap_div = numeric(),
  #                       hap_var = numeric(),
  #                       stringsAsFactors = FALSE)
  
  
  # Read and process haplotype data
  hap.dat <- read.table(hap_file, header = TRUE, sep = " ")
  
  hap.dat <- hap.dat %>%
    filter(feeding_ground. != 0)
  
  hap.tab <- t(as.matrix(table(hap.dat$haplotype.)))
  
  n.samp <- sum(hap.tab)
  samp.prop <- hap.tab / n.samp
  samp.H <- ((n.samp) / (n.samp - 1)) * (1 - sum(samp.prop^2))
  v.H.samp <- (2 / (n.samp * (n.samp - 1))) * (2 * (n.samp - 2)) *
    ((sum(samp.prop^3) - sum(samp.prop^2)^2) + sum(samp.prop^2) - sum(samp.prop^2)^2)
  
  
  return(list(length(hap.tab), samp.H, v.H.samp))
  
}


################################################################ 
#### run test function 
################################################################

######################


FG_loss_test <- function(seeds) {
  
  ### initialise output table
  results <- data.frame(matrix(NA, nrow = length(seeds), ncol = 4))
  colnames(results) <- c("simSEED", "hapCur", "HapDiv", "HapVar")
  
  for (i in seq_along(seeds)) {
    subfolder <- file.path("results", seeds[i])
    hap_files <- list.files(subfolder, pattern = "_haploiso_end\\.txt$", full.names = TRUE)
    
    if (length(hap_files) == 0) {
      warning(paste("No haploiso file found for seed:", seeds[i]))
      next
    }
    
    hap_file <- hap_files[1]
    
    # Attempt to calculate haplo diversity with error handling
    haplo_diversity <- tryCatch({
      calc.hap.div_test(hap_file)
    }, error = function(e) {
      warning(paste("Error reading or processing file for seed", seeds[i], ":", e$message))
      return(NULL)
    })
    
    # If failed, skip
    if (is.null(haplo_diversity) || length(haplo_diversity) != 3) {
      warning(paste("Skipping seed", seeds[i], "- invalid haplo_diversity output"))
      next
    }
    
    # Add result to table
    results[i, ] <- c(
      simSEED = seeds[i],
      hapCur = haplo_diversity[[1]],
      HapDiv = haplo_diversity[[2]],
      HapVar = haplo_diversity[[3]]
    )
  }
  
  return(results)
}



####################################
####################################
####################################
####################################
####################################
####################################

## get seed vector

setwd('/Users/adamglucksman/Local PhD/PopSim')

#source("/Users/adamglucksman/Library/CloudStorage/Dropbox/Vic/PhD/PopSim/analysis/gen_ratios_function.R")

##### setup

# Load data
sim_data_orig <- readRDS('/Users/adamglucksman/Local PhD/PopSim/output_tables/[19-06 18:28:35]_sim_results_table.rds')

seeds <- sim_data_orig$simSEED
seeds <- seeds[!is.na(seeds)]

## run test 

FG_loss_test <- FG_loss_test(seeds)

mean(FG_loss_test$HapDiv, na.rm = T)

mean(sim_data_orig$HapDiv, na.rm = T)
