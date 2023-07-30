#!/bin/bash --login

#SBATCH --job-name=magic

#SBATCH --output=logs/gdv.magic.out.%J

#SBATCH --error=logs/gdv.magic.err.%J

#SBATCH --ntasks=16

#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpusmall

srun /bin/hostname

srun pwd

export TMPDIR=tmp
srun ls $TMPDIR

export MB=/impacs/gdv1/.config/bin/magicblast
srun ls $MB

srun $MB -help

srun $MB -limit_lookup F -no_unaligned -outfmt sam \
    -sra_batch sra -db db/xylo -out out/aln.sam \
    -word_size 12 -penalty -2 -num_threads 16
