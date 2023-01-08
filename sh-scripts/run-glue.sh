#!/bin/bash
#SBATCH -t 6-23:59:00
#SBATCH -J train-glue-model
# #SBATCH --begin=now+3hour
#SBATCH -o train-out.%j
#SBATCH -e train-err.%j
#SBATCH --partition=gpu
#SBATCH --gres=gpu:teslav100:1
# #SBATCH --nodelist=nag16
#SBATCH --exclude=meg[10-12],nag[01-09]
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100G
#SBATCH --mail-type=END
#SBATCH --mail-user=antti.kiviaho@tuni.fi

module load CUDA/11.2
echo "CUDA loaded"

# load conda environment
module load anaconda
source activate spatac

python -u scripts/train-glue-model.py