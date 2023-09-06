#!/bin/bash --login

#SBATCH --job-name=bowtie
#SBATCH --output=logs/gdv.bowtie.out.%J
#SBATCH --error=logs/gdv.bowtie.err.%J

#SBATCH --mail-type=ALL
#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpubig
#SBATCH --ntasks=94

# Create required setup directories
declare -a setup=(data pairs aln logs tmp)
mkdir -p "${setup[@]}" #creates directories listed in setup array
# Setup checks
/bin/hostname         #adds line of text to logfile to state which slurm node sbatch was run on
echo "$PATH"         #adds all directories in PATH to logfile
pwd                   #adds current directory to logfile
which prefetch         #adds binary path for prefetch to logfile
which fasterq-dump     #adds binary path for faster-dump to logfile
which bowtie2         #adds binary path for bowtie2 to logfile
which parallel         #adds binary path for parallel to logfile
export TMPDIR=tmp     #creates tmp directory in home to avoid /tmp disc storage limit
echo $TMPDIR            #adding path of tmp to logfile
ls $TMPDIR            #adding contents of tmp directory to logfile

#function to check that available storage is greater than 100GB, else, pause the addition of more downloads
function chkdsk {
    local avl
    avl=$(df -BG . | awk 'NR==2 {print $4}' | tr -d 'G')
    if [ "$avl" -lt 100 ]; then
        echo "Delaying execution of $1 as only ${avl}GB available"
        while [ "$avl" -lt 100 ]; do
            echo "Sleeping 10 minutes..." 
            sleep "10m"                #pause duration if disk is full
            avl=$(df -BG . | awk 'NR==2 {print $4}' | tr -d 'G')
        done
        echo "Resuming $1 as ${avl}GB now available"
    fi
}

export -f chkdsk #makes chkdsk function to other parts of the code

# extract fastq from SRA file && align && delete fastq after
function align {
    # Check if enough disk space is available before starting job
    chkdsk "$1"
    
    echo "Starting alignment for <>$1</>" #added statement to logfile
    
    # Download compressed SRA files
    # No exists check required thanks to resume
    # if bash script or download fails, resume flag will reuse already downloaded SRA files if not already aligned
    if prefetch "$1" --output-directory data/ \
    --verbose --progress --max-size 50G \
    --resume yes -L debug;
    then
        echo "Downloaded fastq files for $1"
    else
        echo "Download failed for $1"
        exit 1
    fi
    
    # extract fastq files from SRA
    # Skip if file already exists
    if [ ! -f "pairs/$1_1.fastq" ] && [ ! -f "pairs/$1_2.fastq" ]; then
        # Check if space available for extracting more fastq files
        chkdsk "$1"
        if fasterq-dump "data/$1/$1.sra" --outdir pairs/ --temp tmp/ \
        --bufsize 50MB --curcache 500MB --mem 5000MB --threads 94 \ #RAM allocated to each core of 8 parallel jobs (40GB RAM total)
        --progress --verbose --details --log-level debug; #added to logfile
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
        if bowtie2 --threads 94 --time -x db/xylo \
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
--retries 3 --jobs 8 --arg-file ./runs align
