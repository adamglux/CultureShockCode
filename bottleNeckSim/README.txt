Whale Population Bottleneck Read Me

----------------------------------
----------------------------------

Coordinated by Dr Emma Carroll (EC) 
Edited by Adam Glucksman (AG)

----------------------------------

This README file contains a brief description of file contents as well as a brief description of the simulation architecture

----------------------------------

Description of Contents

      bottleNeck_sim_1.R
      
This script file will set up the necessary environment and packages and then run functions to carry out the bottle neck population simulations. You may have to install packages and setup SimuPop separately. You will also have to set a working directory to ensure all necessary files are located.

To Run the simulations, you'll first set up a vector of seeds. This is best run using the `generate_sim_seeds()` function.

the function `run_sim()` will carry out the simulations. This does not need to be saved to a variable, as the programme will save an RDS data.frame separately, along with all the necessary output files to separate locations. 


      simFunctions.R
  
This script file contains all of the necessary functions to carry out the simulations in R, using the python script whale.py. 


      whale.py
      
A python script which uses the simuPop library to run a multi-generational population simulation. 


      snp_1.single.txt
      
Initial SNP values for generation 1 of the population. This file is necessary for initialising and carrying out simulations, but the purposes of this study were to examine the relationships of mtDNA, and therefore no manipulations of the snp_1.single.txt are carried out anywhere in these scripts or files.


      [logs/] 
      
If this sub-directory doesn't exist, it will be created when the `run_sum()` function is run in R. This folder contains logs of python output. These contain population outputs and run outputs for python script. 


      [mtDNA/]

If this sub-directory doesn't exist, it will be created when the `run_sum()` function is run in R. This folder contains initial mtDNA files for each population sim. 


      [output_tables/]

If this sub-directory doesn't exist, it will be created when the `run_sum()` function is run in R. This folder contains console output from R, as well as RDS data.frames from the simulations.  

      [results/]

If this sub-directory doesn't exist, it will be created when the `run_sum()` function is run in R. This folder contains sub-directories of each population's relevant files. The numbered prefixes are the population simulation seed, gathered from the seed matrix
      

----------------------------------

Bottleneck Simulations function architecture for `run_sim()`

1. Initial Setup

   * Record the start time of the function to track how long the simulation takes.
   * Print the current time to show when the simulation started.
   * Check whether the `output_tables/` and `logs/` directories exist. If not, create them.
   * Set up a log file to capture output from each simulation run.

2. Seed Matrix Generation

   * If the user has `seedMatx = TRUE`, create a matrix of random seeds using `generate_sim_seeds()`.
   * This matrix determines the number of simulations and their variation.
   * Print basic information about the simulation (e.g., `Ne1`, dimensions of the seed matrix).

3. Initialize the Results Table

   * A data frame is created to store summary statistics for each simulation.
   * Columns include diversity measures across generations and population sizes.

4. Single Simulation (1x1 Seed Matrix)

   * If the seed matrix has only one value, a single simulation is run.
   * The function sets up the initial state, runs the simulation using Python, and processes the output files.
   * Haplotype diversity is calculated for generations 75–81 and the final generation.
   * Results are stored in the results table.

5. Multiple Simulations (n x m Seed Matrix)

   * For each row of the seed matrix (representing one mtDNA background):

     * Initialize the mtDNA using the first seed in the row.
     * For each remaining seed in the row:

       * Run the simulation using `system2()` to call a Python script with the relevant environment variables.
       * If the simulation fails, record an `NA` row instead.
       * Extract haplotype diversity and population size across generations.
       * Append the results to the results table.

6. Finalization

   * Round numeric results to 5 decimal places.
   * Print the full results and save them as an `.rds` file.
   * Print total runtime in minutes.
   * Close any open file connections.
   * Optionally (if `cleanLogs = TRUE`), delete all files in the `results/` and `logs/` directories and display warnings.

7. Return

   * The function returns the full `results` data frame.
