#!/bin/bash
#SBATCH -J StrainAssemb_FLORIA
#SBATCH --output=sbatch-%x-%j.out # j:jobid, x: job name
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G

# Check if an argument is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 path_to_ecolreads_bam"
    exit 1
fi

MODEL="r1041_e82_400bps_sup_v5.0.0"

BAM=$(realpath $1)
cd $(dirname $BAM)

ASSEMB_DIR=$(basename $BAM .bam)
mkdir $ASSEMB_DIR

ml load samtools
srun samtools fastq -@8 $BAM | pigz -p8 > $ASSEMB_DIR/$ASSEMB_DIR.fastq.gz
FQ=$(realpath $ASSEMB_DIR/$ASSEMB_DIR.fastq.gz)
cd $ASSEMB_DIR
pwd
source /data/brbloemen/mambaforge/etc/profile.d/conda.sh
conda deactivate

ml load flye
srun flye --nano-hq $FQ -i 2 --read-error 0.03 --out-dir Flye -t 8 --asm-coverage 50 --genome-size 5m

ml -flye
ml load samtools minimap2/2.28 htslib/1.21 medaka

mkdir Medaka
srun mini_align -i $FQ -r Flye/assembly.fasta -m -p Medaka/calls_to_draft.tmp -t 8
srun samtools sort Medaka/calls_to_draft.tmp.bam -@8 -o Medaka/calls_to_draft.bam
srun samtools index -@8 Medaka/calls_to_draft.bam

srun medaka inference Medaka/calls_to_draft.bam Medaka/consensus_probs.hdf --model $MODEL --batch 200 --threads 8

srun medaka sequence Medaka/consensus_probs.hdf Flye/assembly.fasta Medaka/consensus.fasta

ml load hmmer prodigal pplacer CheckM
srun checkm lineage_wf -t8 -x fasta Medaka/ CheckM --pplacer_threads 8 > CheckM_summary.txt