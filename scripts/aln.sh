#!/bin/bash --login

#SBATCH --job-name=bowtie

#SBATCH --output=logs/gdv.bowtie.out.%J

#SBATCH --error=logs/gdv.bowtie.err.%J

#SBATCH --ntasks=16

#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpusmall

srun /bin/hostname

srun pwd

declare -a setup=(pairs tmp aln)
mkdir -p "${setup[@]}"
ls -R "${setup[@]}"

function align {
    echo "$0ing $1"
    fasterq-dump ${1} --outdir pairs/ --temp tmp/ \
    --bufsize 10MB --curcache 100MB --mem 4000MB --threads 16 \
    --progress --verbose --details --log-level debug && \
    bowtie2 --threads 16 --time -x db/xylo \
        -q --phred33 --local --very-sensitive-local --no-unal \
        -N 1 -L 12 --rfg 5,2 \
        -1 pairs/${1}_1.fastq -2 pairs/${1}_2.fastq -S aln/${1}.sam && \
    rm -f pairs/${1}_{1,2}.fastq && \
    echo "Aligned $1"
}

export -f align
export parallel=/impacs/gdv1/.local/bin/parallel

srun $parallel --joblog tmp/parallel -j 16 -a runs align
