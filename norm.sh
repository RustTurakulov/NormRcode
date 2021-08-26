#!/bin/sh
#SBATCH --partition=largemem \
#SBATCH --mem=1507g \
#SBATCH --exclusive \
#SBATCH --ntasks=144 \
#SBATCH --ntasks-per-core=2 \
#SBATCH --gres=lscratch:700 \
#SBATCH --time 24:00:00

## Sbatch for running Drew's normalization protocol for idat 
## Check nodes: freen |  Singularity container:  rbox4_v0.sif 
echo "Parameters line parsed: $1"; 
umask 007
cd /data/MDATA
export TMPDIR=/lscratch/$SLURM_JOB_ID
export SINGULARITYENV_THREADS=140
gpfs_dirs="$(echo /gs* | tr ' ' ',')"
export SINGULARITY_BINDPATH="/vf,${gpfs_dirs},/spin1,/data,/lscratch,/fdb"
module load singularity
singularity exec NormRcode/rbox_v0.sif /usr/bin/Rscript NormRcode/pancancer_pipeline_NORM.R $1