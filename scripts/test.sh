#!/bin/bash --login

#SBATCH --job-name=bowtie
#SBATCH --output=logs/gdv.bowtie.out.%J
#SBATCH --error=logs/gdv.bowtie.err.%J

#SBATCH --mail-type=ALL
#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpubig
#SBATCH --ntasks=96

# Setup checks
# echo "No srun $(/bin/hostname)"
# srun echo "With srun $(/bin/hostname)"

function align {
    echo "Step 1: $1"
    echo "Step 2: $1"
    echo "Step 3: $1"
}

export -f align

while IFS= read -r run; do
    echo "Aligning for <> $run </>"
    align "$run"
done < ./runs
