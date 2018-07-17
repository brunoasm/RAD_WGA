#!/bin/bash
#SBATCH -n 16 #Number of cores 
#SBATCH -N 1 #All cores in same machine
#SBATCH -p interact,serial_requeue  #Partition to submit to 
#SBATCH --mem=20000
#SBATCH -J tabulate
#SBATCH -o tabulate.%A.out
#SBATCH -t 2-00:00:00

./tabulate_s7_stats.py -i sample_info_new.csv 
