#!/bin/bash --login

#SBATCH --job-name=bowtie
#SBATCH --output=logs/gdv.bowtie.out.%J
#SBATCH --error=logs/gdv.bowtie.err.%J

#SBATCH --mail-type=ALL
#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpubig
#SBATCH --ntasks=96

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
    if fasterq-dump "${1}" --outdir pairs/ --temp tmp/ \
    --bufsize 10MB --curcache 100MB --mem 16384MB --threads 96 \
    --progress --verbose --details --log-level debug;
        then
            echo "Downloaded fastq files for $1"
        else
            echo "Download failed for $1"
            exit 1
    fi

    # Run bowtie2 on downloaded files
    if [ -f "pairs/$1_1.fastq" ] && [ -f "pairs/$1_2.fastq" ]; then
        if bowtie2 --threads 96 --time -x db/xylo \
            -q --phred33 --local --very-sensitive-local --no-unal \
            -N 1 -L 12 --rfg 5,2 \
            -1 "pairs/$1_1.fastq" -2 "pairs/$1_2.fastq" -S aln/"$1".sam;
            then
                echo "Bowtie2 align completed for $1"
                echo "Deleting fastq files"
                rm -f pairs/"$1"_{1,2}.fastq
            else
                echo "Bowtie2 align failed for $1"
                exit 1
            fi
    else
        echo "Incomplete pairs present for $1"
        exit 1
    fi

    echo "Finished alignment for $1"

}

export -f align

function iter {
    while IFS= read -r run; do
        align "$run"
    done < "$1"
}

export -f iter

srun iter ./runs
