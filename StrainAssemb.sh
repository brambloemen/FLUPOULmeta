#!/bin/bash
#SBATCH -J StrainAssemb
#SBATCH --output=sbatch-%x-%j.out # j:jobid, x: job name
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G

# Check if an argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 path_to_ecolreads"
    exit 1
fi

DIR=$(dirname $1)
cd $DIR 
echo $DOR

ml load flye/2.9.3
srun flye --nano-hq $1 --out-dir Flye_F6D0_enrich_downsampled -t 20 --meta --no-alt-contigs --keep-haplotypes -i 0 --read-error 0.03