#!/bin/bash
#SBATCH -p serial requeue
#SBATCH -t 0-00:05
#SBATCH -c 1
#SBATCH --mem=1000
#SBATCH -e lr-input.%j.err
#SBATCH -o lr-input.%j.out
export R LIBS USER=$HOME/apps/R 4.1.0:$R LIBS USER
module load R/4.1.0-fasrc01
Rscript lin-reg-input-script.R $1