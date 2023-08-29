#!/bin/bash
#SBATCH --partition=storfer
#SBATCH --array=1-5723%13
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=15G
#SBATCH --job-name BayEnv2_r1
#SBATCH --time=6-23:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marc.beer@wsu.edu

snps=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' /data/storfer/cane_toad_2021/analyses/bayenv/gea_spare/split_file_name.txt)  

/data/storfer/cane_toad_2021/analyses/bayenv/gea_spare/bayenv2 -i /data/storfer/cane_toad_2021/analyses/bayenv/gea_spare/tempdir/$snps -m cov_mat_avg -e environfile.txt -p 59 -n 4 -k 750000 -t -c -r 10101 -o bayenv2-snpset-$SLURM_ARRAY_TASK_ID.out 
