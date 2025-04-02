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

# This scirpt runs Floria on a BAM file and a reference genome.
# Dependencies:
# - longshot: https://github.com/pjedge/longshot 
# - floria: https://github.com/bluenote-1577/floria
# - floria-strainer: https://github.com/maxibor/floria-strainer
#  --> Here these have been installed through conda environments. The floria/floria-strainer env is found in workflow/envs/floria.yaml


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
