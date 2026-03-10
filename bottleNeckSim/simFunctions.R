##########################################
##########################################
######## Functions for bottleneck sim
##########################################
##########################################

## created by Adam Glucksman 
## June 2025


################################ 
#### generate initial haplotypes 
################################


## samples from the dirichlet distribution
## sets up initial mtDNA file
## prints frequencies into a text file to be read by simuPop
## input: Ne1 is effective pop size
## freq is a vector of haplotype frequencies 
## seed is any integer — this is used later by the seed_matx below 

starting.div<-function(Ne1, freq, seed) {
  s1<-igraph::sample_dirichlet(1,(Ne1*freq))
  
  # check if file directory exists 
  if (!file.exists("mtDNA/")){ #if the necessary file does not exist, create it
    dir.create("mtDNA/")
  } 

  write.table(s1, paste("mtDNA/",seed, ".mtDNA.txt", sep = ""),
              col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  
}

################################ 
#### calc hap diversity and S.D. 
################################

## calculates haplotype diversity
## input: filepath of generation summary output
## input files should be in the form of {SEED}_haploiso_{GEN}.txt, i.e. 587_haploiso_75.txt
## returns total num of haplotypes, diversity, and variance 

calc.hap.div <- function(hap_file) {
  
  # Initialize results table
  results <- data.frame(seed = character(),
                        hap_div = numeric(),
                        hap_var = numeric(),
                        stringsAsFactors = FALSE)
  
  
  # Read and process haplotype data
  hap.dat <- read.table(hap_file, header = TRUE, sep = " ")
  hap.tab <- t(as.matrix(table(hap.dat$haplotype.)))
  
  n.samp <- sum(hap.tab)
  samp.prop <- hap.tab / n.samp
  samp.H <- ((n.samp) / (n.samp - 1)) * (1 - sum(samp.prop^2))
  v.H.samp <- (2 / (n.samp * (n.samp - 1))) * (2 * (n.samp - 2)) *
    ((sum(samp.prop^3) - sum(samp.prop^2)^2) + sum(samp.prop^2) - sum(samp.prop^2)^2)
  
  
  return(list(length(hap.tab), samp.H, v.H.samp))
  
}

## note: for a simplified version that does not take input file,
## use a vector of haplotypes for hap.tab


################################ 
#### generate seeds
################################

## generates a matrix of seeds to use for simulations
## the first column is used to generate new initial mtDNA files 
## cols 2-j are seeds used to in the python script to initialise populations

generate_sim_seeds <- function(runs.vec) {
  
  #initialise matrix
  seed_matx <- matrix(NA, nrow = runs.vec[1], ncol = runs.vec[2] + 1)
  
  for (i in 1:nrow(seed_matx)) { #rows for mtDNA output seeds
    for (j in 1:ncol(seed_matx)) { #cols for pop_sim in whale.py 
      seed_matx[i,j] <- sample.int(2^24 - 1, 1)
    }
  }
  #seed_hex <- sprintf("0x%08x", sample.int(2^20 - 1, runs.vec[2]))
  
  return(seed_matx)
}


################################ 
#### read generation summary function
################################

## reads a file and turns it into a dataframe 

read_gen_sum <- function(file){
  sum_table <- read.table(file, header = FALSE, sep = "\t")
  df <- t(sum_table[75:82,])
  rownames(df) <- NULL
  colnames(df) <- NULL
  
  return(df)
}




################################################################ 
#### run simulation function 
################################################################

## this function runs the whale.py code and returns a dataframe of results
## input details are listed below
## results returns a dataframe consisting of: 
## mtSEED   col-vector of seeds to initialise mtDNA file
## simSEED  row-vector of seeds to initialise populations 
## hap75, hap76, hap77, hap78, hap79, hap80, hap81   Number of haplotypes for gen 75-81  
## div75, div76, div77, div78, div79, div80, div81   haplotype diversity for gen 75-81
## var75, var76, var77, var78, var79, var80, var81   diversity variance for gen 75-81
## hapCur   number of haplotypes for gen 82 (current generation) 
## HapDiv   haplotype diversity for gen 82 (current generation)
## HapVar   diversity variance for gen 82 (current generation) 
## 75,76,77,78,79,80,81,82   Generation sizes (number of individuals) for gen 75-81


run_sim <- function(runs.vec,
                    haploFreq = hapFreq, 
                    Ne1 = 10, 
                    seedMatx = TRUE, 
                    cleanLogs = TRUE) {
  
  
  ######### setup ########
  #get total runtime of the whole function
  tot.start <- proc.time()
  
  cat("Starting simulation at", format(Sys.time(), "%d-%m %H:%M:%S"), "\n")
  
  # check if file directory exists for output tables 
  if (!file.exists("output_tables/")){ #if the necessary file does not exist, create it
    dir.create("output_tables/")
  } 
  
  # check if file directory exists for logs 
  if (!file.exists("logs/")){ #if the necessary file does not exist, create it
    dir.create("logs/")
  } 
  
  #setup logs
  file_name <- "output_tables/"
  timeStamp <- format(Sys.time(), "[%d-%m %H:%M:%S]")
  filename <- paste0(file_name, timeStamp, "_sim_results_log.txt")
  
  #sink output to log 
  #promptly shut it off after loop runs
  sink(file = filename) 
  #on.exit(if (sink.number() > 0) sink(NULL), add = TRUE)
  
  cat("Starting simulation at", format(Sys.time(), "%d-%m %H:%M:%S"), "\n")
  
  #generate matrix of seeds 
  # leave turned on to generate a unique seedMatrix for each run
  # turn off for comparison tests:
  ## if turned off, generate matrix before running function
  if (seedMatx) {
    seedMatrix <- generate_sim_seeds(runs.vec = runs.vec)
    #print(seedMatrix) 
  } 
  
  ### print info about the simulation run
  cat("Ne1 set to", Ne1, "\n")
  cat("seed matrix dim:", dim(seedMatrix), "\n")
  
  
  ### initialise output table
  # Initialize results table
  generations_cols <- as.character(75:82)
  results <- data.frame(matrix(NA, nrow = nrow(seedMatrix) * (ncol(seedMatrix) - 1), ncol = 34))
  colnames(results) <- c("mtSEED", "simSEED", "hap75", "div75", "var75", "hap76", "div76", "var76",
                         "hap77", "div77", "var77","hap78", "div78", "var78","hap79", "div79", "var79",
                         "hap80", "div80", "var80","hap81", "div81", "var81", "hapCur",
                         "HapDiv", "HapVar", generations_cols)
  
  ######### run 1 x 1 seed matx ########
  
  if (nrow(seedMatrix) == 1 && ncol(seedMatrix) == 1) {
    startingDiv <- starting.div(Ne1=Ne1, freq=hapFreq, seed = seedMatrix[1,1])
    log_file <- sprintf("logs/run_%s.log", seedMatrix[i,j]) #store python output 
    system2("python3", args = "whale.py", env = c(paste0("SEED=", seedMatrix[i,j]),
                                                  paste0("mtDNA_SEED=", seedMatrix[i,1])), stdout = log_file, stderr = log_file)
    
    ## store results
    
    # check if file directory exists 
    if (!file.exists("results/")){ #if the necessary file does not exist, create it
      dir.create("results/")
    } 
    
    # Build pattern to find the matching haploiso output file
    subfolder <- file.path("results", seedMatrix[i,j])
    hap_file <- list.files(subfolder, pattern = "_haploiso_end\\.txt$", full.names = TRUE)
    
    #calc haplotype diversity
    haplo_diversity <- calc.hap.div(hap_file)
    
    ### get population sizes for each generation 
    gen_file <- list.files(subfolder, pattern = "gen_summary.txt", full.names = TRUE)
    generationSize <- read_gen_sum(gen_file)
    
    ### get hap diversity
    #initial_div_file <- list.files(subfolder, pattern = "_haploiso_75\\.txt$", full.names = TRUE)
    #initial_div <- calc.hap.div(initial_div_file)
    # Pre-allocate list to store diversity results for each generation
    gen_div <- list()
    
    # Get all haploiso files (make sure they're sorted correctly)
    div_files <- list.files(subfolder, pattern = "_haploiso_.*\\.txt$", full.names = TRUE)
    # Remove any files containing "None"
    #div_files <- div_files[!grepl("None", div_files)]
    div_files <- div_files[!grepl("end", div_files)]
    # Sort files by generation number
    div_files <- div_files[order(as.numeric(gsub("[^0-9]", "", div_files)))]
    
    
    for (i in seq_along(div_files)) {
      file <- div_files[i]
      result <- calc.hap.div(file)
      gen_div[[i]] <- result
    }
    
    
    # Add result to table
    results[1, colnames(results)] <- c(
      list(
        seedMatrix[1, 1],
        seedMatrix[1, 1],
        
        gen_div[[1]][[1]], #gen75
        gen_div[[1]][[2]], #gen75
        gen_div[[1]][[3]], #gen75
        
        gen_div[[2]][[1]], #gen76
        gen_div[[2]][[2]], #gen76
        gen_div[[2]][[3]], #gen76
        
        gen_div[[3]][[1]], #gen77
        gen_div[[3]][[2]], #gen77
        gen_div[[3]][[3]], #gen77
        
        gen_div[[4]][[1]], #gen78
        gen_div[[4]][[2]], #gen78
        gen_div[[4]][[3]], #gen78
        
        gen_div[[5]][[1]], #gen79
        gen_div[[5]][[2]], #gen79
        gen_div[[5]][[3]], #gen79
        
        gen_div[[6]][[1]], #gen80
        gen_div[[6]][[2]], #gen80
        gen_div[[6]][[3]], #gen80
        
        gen_div[[7]][[1]], #gen81
        gen_div[[7]][[2]], #gen81
        gen_div[[7]][[3]], #gen81
        
        #initial_div[[3]],
        #initial_div[[1]],
        #initial_div[[2]],
        haplo_diversity[[1]],#current
        haplo_diversity[[2]],#current
        haplo_diversity[[3]]#current
      ),
      #((initial_div[[1]] - haplo_diversity[[1]]) / initial_div[[1]])*100),
      
      setNames(as.list(generationSize[2, ]), generations_cols))
    
    
    ######### run n x m seed matx ########  
  } else {
    
    idx <- 1 #index to store results
    
    for (i in 1:nrow(seedMatrix)) { #initialise mtDNA 
      
      ## generate initial mtDNA file
      set.seed(seedMatrix[i,1])
      #generate Hap frequencies
      startingDiv <- starting.div(Ne1=10, freq=hapFreq, seed = seedMatrix[i,1])
      
      for (j in 2:ncol(seedMatrix)) { #step through each sim per mtDNA file
        
        ## run sim
        log_file <- sprintf("logs/run_%s.log", seedMatrix[i,j]) #store python output 
        #####
        #### tryCatch incase the sim malfunctions 
        ####
        tryCatch({system2("python3", args = "whale.py", env = c(paste0("SEED=", seedMatrix[i,j]),
                                                                paste0("mtDNA_SEED=", seedMatrix[i,1])), stdout = log_file, stderr = log_file)
          
          ## store results
          # Build pattern to find the matching haploiso output file
          subfolder <- file.path("results", seedMatrix[i,j])
          hap_file <- list.files(subfolder, pattern = "_haploiso_end\\.txt$", full.names = TRUE)
          
          #calc haplotype diversity
          haplo_diversity <- calc.hap.div(hap_file)
          
          ### get population sizes for each generation 
          gen_file <- list.files(subfolder, pattern = "gen_summary.txt", full.names = TRUE)
          generationSize <- read_gen_sum(gen_file)
          
          ### get hap diversity
          #initial_div_file <- list.files(subfolder, pattern = "_haploiso_75\\.txt$", full.names = TRUE)
          #initial_div <- calc.hap.div(initial_div_file)
          gen_div <- list()
          
          # Get all haploiso files (make sure they're sorted correctly)
          div_files <- list.files(subfolder, pattern = "_haploiso_.*\\.txt$", full.names = TRUE)
          # Remove any files containing "None"
          #div_files <- div_files[!grepl("None", div_files)]
          div_files <- div_files[!grepl("end", div_files)]
          # Sort files by generation number
          div_files <- div_files[order(as.numeric(gsub("[^0-9]", "", div_files)))]
          
          
          for (k in seq_along(div_files)) {
            file <- div_files[k]
            result <- calc.hap.div(file)
            gen_div[[k]] <- result
          }
          
          # Add result to table
          results[idx, colnames(results)] <- c(
            list(
              seedMatrix[i, 1],
              seedMatrix[i, j],
              
              gen_div[[1]][[1]], #gen75
              gen_div[[1]][[2]], #gen75
              gen_div[[1]][[3]], #gen75
              
              gen_div[[2]][[1]], #gen76
              gen_div[[2]][[2]], #gen76
              gen_div[[2]][[3]], #gen76
              
              gen_div[[3]][[1]], #gen77
              gen_div[[3]][[2]], #gen77
              gen_div[[3]][[3]], #gen77
              
              gen_div[[4]][[1]], #gen78
              gen_div[[4]][[2]], #gen78
              gen_div[[4]][[3]], #gen78
              
              gen_div[[5]][[1]], #gen79
              gen_div[[5]][[2]], #gen79
              gen_div[[5]][[3]], #gen79
              
              gen_div[[6]][[1]], #gen80
              gen_div[[6]][[2]], #gen80
              gen_div[[6]][[3]], #gen80
              
              gen_div[[7]][[1]], #gen81
              gen_div[[7]][[2]], #gen81
              gen_div[[7]][[3]], #gen81
              
              #initial_div[[3]],
              #initial_div[[1]],
              #initial_div[[2]],
              haplo_diversity[[1]],
              haplo_diversity[[2]],
              haplo_diversity[[3]]
            ),
            # ((initial_div[[1]] - haplo_diversity[[1]]) / initial_div[[1]])*100),
            setNames(as.list(generationSize[2, ]), generations_cols))
          
          #####
          #### for errors in the simulation: 
          ####
        }, error = function(e) {
          cat("Run failed for seed", seedMatrix[i,j], "– inserting NA row\n")
          
          # Return a default row with the same structure but with NAs
          # Add result to table
          results[idx, ] <- NA
        })
        
        #print(results[idx,]) #print each line
        cat("Run number:", idx, "\n")
        
        #update index 
        idx <- idx + 1
        
      } #end of the j loop
      
      
    } #end of the i loop 
  } #end of the else statement
  
  #round numeric values to 5dp
  results[ , sapply(results, is.numeric)] <- round(results[ , sapply(results, is.numeric)], 5)
  
  print(results)
  saveRDS(results, file = paste0(file_name, timeStamp, "_sim_results_table.rds"))
  
  tot.end <- proc.time()
  tot.time <- tot.end - tot.start
  cat("\ntotal elapsed time:", round(tot.time[3]/60,2), "minutes\n")
  
  closeAllConnections()
  
  cat("finished run, output file:", filename,"\n")
  cat("total elapsed time:", round(tot.time[3]/60,2), "minutes\n")
  
  if (cleanLogs) {
    ## [WARNING] Clear old files before starting
    unlink("results/*", recursive = TRUE)
    warning("results folder cleared of all data. to turn this off, set cleanLogs = FALSE")
    unlink("logs/*", recursive = TRUE)
    warning("logs folder cleared of all data. to turn this off, set cleanLogs = FALSE")
  }
  
  return(results)
  
} #end of function 



