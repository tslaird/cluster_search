#!/bin/bash -l

# Name of the job
#SBATCH -J name    # name of the job
#SBATCH -o %x-%j.out           # stdout
#SBATCH -e %x-%j.err           # stderr
#SBATCH -p med                 # partition, or queue, to assign to
#SBATCH -N 1                   # number of nodes or computer
#SBATCH -n 1                   # one task for this node
#SBATCH -c 16                  # cores per task
#SBATCH -t 120:00:00           # maximum run time
#SBATCH --mail-type=ALL        # mail all reports
#SBATCH --mail-user=tslaird@ucdavis.edu  # email

# initialize conda
. ~/miniconda3/etc/profile.d/conda.sh

# write commands below
snakemake --cores 16 --latency-wait 4000 --use-conda
