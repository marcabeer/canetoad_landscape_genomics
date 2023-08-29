#!/bin/bash
#SBATCH --partition=storfer
#SBATCH --cpus-per-task=9
#SBATCH --job-name Bayenv_mat5
#SBATCH --time=6-23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marc.beer@wsu.edu



/data/storfer/cane_toad_2021/analyses/bayenv/cov_matrix/bayenv2 -i /data/storfer/cane_toad_2021/analyses/bayenv/cov_matrix/ct_bayenv_05022022.txt -p 59 -k 500000 -r 50105 > bayenv_matrix05 
