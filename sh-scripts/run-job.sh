#!/bin/bash
#SBATCH -t 23:59:00
#SBATCH -J preprocess
#SBATCH -o preprocess.out.%j
#SBATCH -e preprocess.err.%j
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

# load conda environment
module load anaconda
source activate spatac

python3 -u preprocess-data.py


