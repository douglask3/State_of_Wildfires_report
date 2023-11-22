#!/bin/bash 
#SBATCH --partition=par-multi  #high #long
#SBATCH --ntasks-per-node=1            # Number of tasks per node
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --output=/work/scratch-pw2/mbarbosa/output_%j.txt         # Output file (%j expands to job ID)
#SBATCH --error=/work/scratch-pw2/mbarbosa/error_%j.txt           # Error file (%j expands to job ID)
#SBATCH --time=22:00:00 

conda activate pymc5_env
#conda activate maxent
python pymc_MaxEnt.py
