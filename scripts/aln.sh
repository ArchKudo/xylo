#!/bin/bash --login

#SBATCH --job-name=bowtie
#SBATCH --output=logs/gdv.bowtie.out.%J
#SBATCH --error=logs/gdv.bowtie.err.%J

#SBATCH --mail-type=ALL
#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpubig
#SBATCH --ntasks=96

# Create required setup directories
declare -a setup=(data pairs aln logs tmp)
mkdir -p "${setup[@]}"
# Setup checks
/bin/hostname
echo "$PATH"
pwd
which prefetch
which fasterq-dump
which bowtie2
which parallel
export TMPDIR=tmp
echo $TMPDIR
ls $TMPDIR

function chkdsk {
    local avl
    avl=$(df -BG . | awk 'NR==2 {print $4}' | tr -d 'G')
    if [ "$avl" -lt 100 ]; then
        echo "Delaying execution of $1 as only ${avl}GB available"
        while [ "$avl" -lt 100 ]; do
            echo "Sleeping 10 minutes..."
            sleep "10m"
            avl=$(df -BG . | awk 'NR==2 {print $4}' | tr -d 'G')
        done
        echo "Resuming $1 as ${avl}GB now available"
    fi
}

export -f chkdsk

# get fastq && align && delete fastq
function align {
    # Check if enough disk space is avl before starting job
    chkdsk "$1"
    
    echo "Starting alignment for <>$1</>"
    
    # Download compressed sra files
    # No exists check required thanks to resume
    if prefetch "$1" --output-directory data/ \
    --verbose --progress \
    --resume yes -L debug;
    then
        echo "Downloaded fastq files for $1"
    else
        echo "Download failed for $1"
        exit 1
    fi
    
    # Dump fastq files from sra
    # Skip if file already exists
    if [ ! -f "pairs/$1_1.fastq" ] && [ ! -f "pairs/$1_2.fastq" ]; then
        # Check if space available for extracting more fastq files
        chkdsk "$1"
        if fasterq-dump "data/$1/$1.sra" --outdir pairs/ --temp tmp/ \
        --bufsize 10MB --curcache 100MB --mem 4000MB --threads 96 \
        --progress --verbose --details --log-level debug;
        then
            echo "Extracted fastq files for $1"
        else
            echo "Extraction failed for $1"
            exit 1
        fi
    else
        echo "Skipping extraction as file: $1_1.fastq & $1_2.fastq already exists"
    fi
    
    # Run bowtie2 on downloaded files
    # Check if pairs exist
    if [ -f "pairs/$1_1.fastq" ] && [ -f "pairs/$1_2.fastq" ]; then
        if bowtie2 --threads 96 --time -x db/xylo \
        -q --phred33 --local --very-sensitive-local --no-unal \
        -N 1 -L 12 --rfg 5,2 \
        -1 "pairs/$1_1.fastq" -2 "pairs/$1_2.fastq" -S aln/"$1".sam;
        then
            echo "Bowtie2 align completed for $1"
            echo "Deleting sra files"
            rm -rf data/"$1"
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

parallel --joblog "logs/parallel.$SLURM_JOB_ID" --tmpdir tmp/ \
--compress --keep-order --group \
--retries 3 --jobs 10 --arg-file ./reduced align
