#!/bin/sh
#SBATCH --partition=largemem \
#SBATCH --mem=1507g \
#SBATCH --ntasks=144 \
###SBATCH --mem=247g \
###SBATCH --ntasks=54 \
#SBATCH --exclusive \
#SBATCH --ntasks-per-core=2 \
#SBATCH --gres=lscratch:300 \
#SBATCH --time 24:00:00

cd /data/MDATA
export TMPDIR=/lscratch/$SLURM_JOB_ID
export SINGULARITYENV_THREADS=142
source /usr/local/current/singularity/app_conf/sing_binds
module load singularity 
singularity exec NormRcode/rbox4met_v1.5.sif /usr/bin/Rscript NormRcode/CNS13Krefset/normalization_cns.R 
#singularity exec  NormRcode/rbox_v0.sif /usr/bin/Rscript NormRcode/CNS13Krefset/normalization_cns_topup.R 
