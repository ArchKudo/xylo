#!/bin/bash --login

#SBATCH --job-name=bowtie
#SBATCH --output=logs/gdv.bowtie.out.%J
#SBATCH --error=logs/gdv.bowtie.err.%J

#SBATCH --mail-type=ALL
#SBATCH --mail-user=gdv1@aber.ac.uk

#SBATCH --partition=cpusmall
#SBATCH --ntasks=16

# Staging

mkdir -p staging
cd staging || exit

# Java binary in slurm is too old please update Java before running the next command
# Use picard to merge all sam files into one
java -jar ~/workplace/downloads/picard.jar MergeSamFiles \
$(find ../aln -maxdepth 1 -type f -name "*.sam" -printf "-I %p ") \
-O "aln.sam"

# Binarize the sam file
samtools view -@ 16 --verbosity 100 -b aln.sam -o aln.bam

# Sort the bam file
samtools sort -m 512M --threads 16 --verbosity 100 -o sort.bam aln.bam

# Mark duplicates in the sorted bam file
java -jar ~/workplace/downloads/picard.jar MarkDuplicates \
INPUT=sort.bam \
OUTPUT=dedup.bam \
METRICS_FILE=dedup.metrics

# Create an index file required by gretel
samtools index dedup.bam

# Do pileup and variant calling
bcftools mpileup -Ou --threads 12 -f ../xylo.fasta dedup.bam | \
bcftools call --threads 12 -mv -Oz -o calls.vcf

# Reperform variant calling using snpper
# Get all gene contigs from the reference fasta file
# Doesn't work in slurm due to PyVCF currently broken on Python3.7+ need to run this locally
grep -o -E "^>[^ ]+" xylo.fasta | cut -c 2- | while read -r gene; do
    # Run snpper.py for each gene and save the output to a VCF file
    anaconda3-launch --env xylo3.7.16 python ../scripts/snpper.py -b dedup.bam -r "$gene" > "$gene.vcf"
    bgzip -k -@ 16 "$gene.vcf"
    tabix "$gene.vcf.gz"
    anaconda3-launch gretel dedup.bam "$gene.vcf.gz" "$gene" --master ../xylo.fasta
done
