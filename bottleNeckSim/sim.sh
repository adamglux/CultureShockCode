#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=200G
#SBATCH --partition=longrun
#SBATCH --time=10-05:00:00
#SBATCH -o /nfs/home/glucksad/bottleNeckSim/bottleNeckSim_%A_%a.out
#SBATCH -e /nfs/home/glucksad/bottleNeckSim/bottleNeckSim_%A_%a.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=adam.glucksman@vuw.ac.nz
#SBATCH --job-name=popSim


ml purge; ml GCC/11.2.0 OpenMPI/4.1.1 R/4.2.0 Miniconda3/23.9.0-0 ; export PIP_NO_CACHE_DIR=1; export PYTHONNOUSERSITE=1; source $(conda info --base)/etc/profile.d/conda.sh;

unset PYTHONPATH
conda activate simupop-env

Rscript /nfs/home/glucksad/bottleNeckSim/bottleNeck_sim_1.R 