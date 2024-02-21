#!/bin/bash

#SBATCH -c 8                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
# -N 1 means all cores will be on the same node)
#SBATCH -p short		# partition name
#SBATCH -t 0-03:00 		# hours:minutes runlimit after which job will be killed
#SBATCH --mem 150G 		# amount of memory requested
#SBATCH --job-name FNRMS 		# Job name
#SBATCH -o /n/scratch/users/s/sad167/RMS/logs/sbatch_out.out		# File to which standard out will be written
#SBATCH -e /n/scratch/users/s/sad167/RMS/logs/sbatch_err.err 		# File to which standard err will be written

set -e

# Load required modules
module load gcc/9.2.0 R/4.2.1 udunits/2.2.28 geos/3.10.2 python/3.8.12

# Run the scripts
Rscript /n/scratch/users/s/sad167/RMS/codes/11b_Heatmap_aggregate_O2.R