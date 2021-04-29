#!/bin/sh
#SBATCH --exclusive \
#SBATCH --ntasks=28 \
#SBATCH --ntasks=1 \
#SBATCH --mem=247g \
#SBATCH --ntasks-per-core=28 \
#SBATCH --ntasks-per-core=1 \
#SBATCH --gres=lscratch:100 \
#SBATCH --time 24:00:00

## Sbatch for running Drew's normalization protocol for idat 
## Check nodes: freen |  Singularity container:  rbox4_v0.sif 

cd /data/MDATA
export TMPDIR=/lscratch/$SLURM_JOB_ID
export SINGULARITYENV_THREADS=26
gpfs_dirs="$(echo /gs* | tr ' ' ',')"
export SINGULARITY_BINDPATH="/vf,${gpfs_dirs},/spin1,/data,/lscratch,/scratch,/fdb"
module load singularity
singularity exec  NormRcode/rbox_v0.sif /usr/bin/Rscript NormRcode/normalization.R $1
