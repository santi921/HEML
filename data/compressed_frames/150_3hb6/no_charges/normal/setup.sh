#!/bin/bash
#SBATCH -p RM-shared
#SBATCH -t 08:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=santiagovargas921@gmail.com 
#SBATCH --mail-type=ALL

module load AI/anaconda3-tf2.2020.11
conda activate py27
setupturbomole.py -t

