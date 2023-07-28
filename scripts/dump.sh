#!/bin/bash --login

#SBATCH --job-name=dump

#SBATCH --output=gdv.dump.out.%J

#SBATCH --error=gdv.dump.err.%J

#SBATCH --ntasks=16

# Doesn't work!
#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpusmall

#SBATCH --time=00:01:00

srun /bin/hostname

$FD=~/.config/sratools/bin/fasterq-dump

srun $FD -help | head

for sra in $(cat sra); do srun $FD $sra -o sras/$sra.fastq -p -e 16 -v -L debug --concatenate-reads; done
