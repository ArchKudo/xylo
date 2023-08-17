#!/bin/bash --login

#SBATCH --job-name=bowtie
#SBATCH --output=logs/gdv.bowtie.out.%J
#SBATCH --error=logs/gdv.bowtie.err.%J

#SBATCH --mail-type=ALL
#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpusmall
#SBATCH --ntasks=16

# Setup checks
srun /bin/hostname
srun echo "$PATH"
srun pwd
srun which fasterq-dump
srun which bowtie2
srun which parallel
export TMPDIR=tmp
srun echo $TMPDIR
srun ls $TMPDIR

# Create required setup directories
declare -a setup=(pairs tmp aln logs)
mkdir -p "${setup[@]}"
ls -R "${setup[@]}"

# get fastq && align && delete fastq
function align {
    echo "Starting alignment for $1"

    # Dump fastq files
    fasterq-dump "${1}" --outdir pairs/ --temp tmp/ \
    --bufsize 10MB --curcache 100MB --mem 4000MB --threads 16 \
    --progress --verbose --details --log-level debug

    # Check if fasterq-dump command succeeded
    if [ $? -eq 0 ]; then
        echo "Downloaded fastq files for $1"
    else
        echo "Download failed for $1"
        exit 1
    fi

    # Run bowtie2 on downloaded files
    if [ -f "pairs/$1_1.fastq" ] && [ -f "pairs/$1_2.fastq" ]; then
        bowtie2 --threads 16 --time -x db/xylo \
            -q --phred33 --local --very-sensitive-local --no-unal \
            -N 1 -L 12 --rfg 5,2 \
            -1 pairs/"$1"_1.fastq -2 pairs/"$1"_2.fastq -S aln/"$1".sam
    else
        echo "Incomplete pairs present for $1"
        exit 1
    fi

    # Check if bowtie2 command succeeded
    if [ $? -eq 0 ]; then
        echo "Bowtie2 align completed for $1"
        echo "Deleting fastq files"
        rm -f pairs/$1_{1,2}.fastq
    else
        echo "Bowtie2 align failed for $1"
        exit 1
    fi

    echo "Finished alignment for $1"

}

export -f align

srun parallel --joblog logs/parallel --tmpdir $TMPDIR \
    --compress --keep-order --line-buffer \
    --retries 3 --jobs 2 --arg-file runs align
