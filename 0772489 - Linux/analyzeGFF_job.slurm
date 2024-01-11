#!/bin/bash -l
#SBATCH --account=lp_edu_bioinformatics
#SBATCH --clusters=genius
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --output=output_slurm.log
#SBATCH --error=error_slurm.log

#Specify permission to execute bash script
chmod u+x analyze_GFF_features.sh

#Analyze features of chromosome 1
./analyze_GFF_features.sh 1

echo "Slurm job is launched. Check output file for more information."
