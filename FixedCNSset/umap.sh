#!/bin/sh
#SBATCH --mem=247g \
#SBATCH --ntasks=56 \
#SBATCH --exclusive \
#SBATCH --ntasks-per-core=2 \
#SBATCH --gres=lscratch:300 \
#SBATCH --time 24:00:00

outdir="/data/MDATA/TRANSFER/SAMPLESHEETS/WeeklyCNS"
newsmpls="/data/MDATA/NormRcode/CNS13Krefset/newsamples.csv"

DIR=${1:-$outdir}               # Defaults to outdir
NEWSAMPLES=${2:-$newsmpls}      # Default value is newsmpls

echo $DIR 
echo $NEWSAMPLES

cd /data/MDATA
export TMPDIR=/lscratch/$SLURM_JOB_ID
export SINGULARITYENV_THREADS=54
source /usr/local/current/singularity/app_conf/sing_binds
module load singularity 
singularity exec  NormRcode/rbox4met_v1.2.sif /usr/bin/Rscript NormRcode/FixedCNSset/umap.R $DIR $NEWSAMPLES
