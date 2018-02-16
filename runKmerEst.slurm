#!/bin/sh
#SBATCH --job-name=kmerEst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=50:00:00
#SBATCH --mem=200gb
#SBATCH --output=kmerEst.out
#SBATCH --error=kmerEst.err
module load compiler/gcc/7.1 
/usr/bin/time -v ./kmerEst ../../ERR309932/ERR309932.fastq 31 50000000 0 
