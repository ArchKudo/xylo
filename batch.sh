#!/bin/bash --login

#SBATCH --job-name=xylo

#SBATCH --output=gdv.xylo.out.%J

#SBATCH --error=gdv.xylo.err.%J

#SBATCH --ntasks=16

#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpusmall

srun /bin/hostname

MB=~/.config/bin/magicblast

srun $MB -help

srun $MB -limit_lookup F -no_unaligned -outfmt sam \
    -sra_batch sra -db db/xylo -out out/aln.sam \
    -word_size 12 -penalty -2 -gapopen 5 -gapextend 2 -num_threads 16
