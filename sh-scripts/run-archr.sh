#!/bin/bash
#SBATCH -t 1-23:59:00
#SBATCH -J archr
#SBATCH -o archr.out.%j
#SBATCH -e archr.err.%j
#SBATCH --ntasks=1
#SBATCH --partition=normal
#SBATCH --cpus-per-task=15
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

# load conda environment
module load anaconda
source activate archr

Rscript /lustre/scratch/kiviaho/spatac/scripts/TonATAC-ArchR-preprocessing-nfr-peaks.R
