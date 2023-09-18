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

# Create directories listed in setup array
mkdir -p "${setup[@]}"

# Setup checks

# Add line of text to logfile to state which slurm node sbatch was run on
/bin/hostname

# Add all directories in PATH to logfile
echo "$PATH"

# Add current directory to logfile
pwd

# Add binary path for prefetch to logfile
which prefetch

# Add binary path for faster-dump to logfile
which fasterq-dump

# Add binary path for bowtie2 to logfile
which bowtie2

# Add binary path for parallel to logfile
which parallel

# Create tmp directory environmental variable
# Required for magicblast to avoid the 200G /tmp limit
export TMPDIR=tmp

# Add path of tmp to logfile
echo $TMPDIR

# Add contents of tmp directory to logfile
ls $TMPDIR

# Function to check that available storage is greater than 100GB,
# else, pause the addition of more downloads
function chkdsk {
    
    # Create a local variable
    local avl
    
    # Get the disk space where the current directory is mounted
    # Get the 4th column of second line
    # Trim the G at the end
    avl=$(df -BG . | awk 'NR==2 {print $4}' | tr -d 'G')
    if [ "$avl" -lt 100 ]; then
        echo "Delaying execution of $1 as only ${avl}GB available"
        
        # Continually check for disk storage
        while [ "$avl" -lt 100 ]; do
            echo "Sleeping 10 minutes..."
            
            # Pause duration for 10 minutes while disk is full
            sleep "10m"
            
            # Re-check for available storage
            avl=$(df -BG . | awk 'NR==2 {print $4}' | tr -d 'G')
        done
        
        echo "Resuming $1 as ${avl}GB now available"
    fi
}

# Explicityly make chkdsk function available to other parts of the code
export -f chkdsk

# Download SRA files,
# Extract fastq from SRA file
# Align them
# If successful delete all artifacts
function align {
    
    # Check if enough disk space is available before starting job
    chkdsk "$1"
    
    # Add statement to logfile
    echo "Starting alignment for <>$1</>"
    
    # Download compressed SRA files
    # If remaining bash script or download fails,
    # resume flag will reuse already downloaded SRA files if not already aligned
    if prefetch "$1" --output-directory data/ \
    --verbose --progress --max-size 50G \
    --resume yes -L debug; then
        echo "Downloaded fastq files for $1"
    else
        echo "Download failed for $1"
        exit 1
    fi
    
    # Extract fastq files from SRA
    # Skip if file already exists
    if [ ! -f "pairs/$1_1.fastq" ] && [ ! -f "pairs/$1_2.fastq" ]; then
        # Check if space available for extracting more fastq files
        chkdsk "$1"
        
        # Logs added to logfile
        # Arbitary memory requirements can be further optimized dynamically
        if fasterq-dump "data/$1/$1.sra" --outdir pairs/ --temp tmp/ \
        --bufsize 50MB --curcache 500MB --mem 5000MB --threads 94 \
        --progress --verbose --details --log-level debug; then
            echo "Extracted fastq files for $1"
        else
            echo "Extraction failed for $1"
            
            # Fast exit script on error for that run
            exit 1
        fi
    else
        echo "Skipping extraction as file: $1_1.fastq & $1_2.fastq already exists"
    fi
    
    # Run bowtie2 on downloaded files
    # Check if both pairs exist
    if [ -f "pairs/$1_1.fastq" ] && [ -f "pairs/$1_2.fastq" ]; then
        if bowtie2 --threads 94 --time -x db/xylo \
        -q --phred33 --local --very-sensitive-local --no-unal \
        -N 1 -L 12 --rfg 5,2 \
        -1 "pairs/$1_1.fastq" -2 "pairs/$1_2.fastq" -S aln/"$1".sam; then
            # Clean-up artifacts
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

# Explicitly make align function available to other parts of the code
export -f align

# Start 8 parallel align task for the each run accession in ./runs file seperated by newline
parallel --joblog "logs/parallel.$SLURM_JOB_ID" --tmpdir tmp/ \
--compress --keep-order --group \
--retries 3 --jobs 8 --arg-file ./runs align
