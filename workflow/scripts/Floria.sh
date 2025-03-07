#!/bin/bash
#SBATCH -J FLORIA
#SBATCH --output=sbatch-%x-%j.out # j:jobid, x: job name
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=15G

# Check if an argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 path_to_bam path_to_assemb"
    exit 1
fi

MODEL="r1041_e82_400bps_sup_v5.0.0"

BAM=$(realpath $1)
FA=$(realpath $2)
DIR=$(dirname $BAM)

cd $DIR

source /data/brbloemen/mambaforge/etc/profile.d/conda.sh
conda activate longshot

srun longshot --bam $BAM --ref $FA --out Ecol_variants.vcf

conda deactivate
conda activate floria

srun floria -b $BAM -v Ecol_variants.vcf -r $FA -o Floria -t 10
srun floria-strainer -b $BAM Floria/ -m split
