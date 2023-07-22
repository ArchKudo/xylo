#! /bin/bash
# Get first 5000 reads
head -n 5000 ../data/ERR3201375.fastq >fivek.fastq
# Fetch ecoli reference genome
# https://www.ncbi.nlm.nih.gov/datasets/taxonomy/562/
# Create a bwa index
bwa index ecoli.fna
# Align the reads to ecoli
bwa mem -t 12 ecoli.fna fivek.fastq >aln.sam
# Convert to bam file
samtools view -b aln.sam >aln.bam

samtools view -@ 12 --verbosity 100 -b aln.sam -o aln.bam

# Sort the file
samtools sort -m 512M --threads 12 --verbosity 100 -o sort.bam aln.bam
# Dedup
java -jar ~/workplace/downloads/picard.jar MarkDuplicates \
    INPUT=sort.bam \
    OUTPUT=dedup.bam \
    METRICS_FILE=dedup.metrics
# Add fake read-groups
java -jar ~/workplace/downloads/picard.jar AddOrReplaceReadGroups \
    INPUT=dedup.bam \
    OUTPUT=rg.bam \
    SORT_ORDER=coordinate \
    RGLB=10678_0001 \
    RGPL=ILLUMINA \
    RGPU=ERR3201375 \
    RGSM=SAMEA5383982 \
    RGID=ERR3201375 \
    RGCN=PRJEB31266 \
    RGDS="Cattle rumen microbiome. breed:Aax; sex:M; age:589" \
    CREATE_INDEX=True
# Validate for errors
java -jar ~/workplace/downloads/picard.jar ValidateSamFile -I rg.bam -MODE SUMMARY
# Create dict & vcf
samtools faidx ecoli.fna
~/workplace/gatk/gatk CreateSequenceDictionary -R ecoli.fna
~/workplace/gatk/gatk HaplotypeCaller -R ecoli.fna -I rg.bam -O var.vcf
# bgzip + tabix
bgzip var.vcf
tabix var.vcf.gz
# Convert fastq to fasta
# https://bioinformaticsworkbook.org/dataWrangling/fastaq-manipulations/converting-fastq-format-to-fasta.html#gsc.tab=0
sed -n '1~4s/^@/>/p;2~4p' fivek.fastq >fivek.fasta
# Install hmmer3
cpanm Bio::SearchIO::hmmer3 --sudo
# Run prokka
prokka --cpus 0 fivek.fastq --prefix fk --outdir ann # not required
# Run gretel
# gretel rg.bam var.vcf.gz ecoli -s 3694124 -e 3695164 -o snp.fasta (Doesn't work)
# Do variant calling using bcftools
bcftools mpileup -Ou --threads 12 -f MF045434.1.fna ERR3201375.aln.bam | bcftools call --threads 12 -mv -Ob -o ERR3201375.calls.bcf
# Filter for mapped sequences
samtools view -F 4 -h ERR3201375.aln.bam
# Align using bowtie
bowtie2-build MF045434.1.fna MF045434.1.fna.idx
bowtie2 -p 12 -x MF045434.1.fna.bt2.idx -U ERR3201375.fastq -S ERR3201375.aln.bt2.sam --local

bowtie2 -p 12 -x MF045434.1.fna.idx -U ERR3201409.fastq -S ERR3201409.aln.sam --local

# Search using blastn
makeblastdb -in lambda_virus.fa -dbtype nucl -parse_seqids -out lambda_virus
blastn -query longreads.fq -db lambda_virus -out xl.sam -task blastn -outfmt 17
blastn -outfmt \
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out summary.txt -num_threads 4 -query longreads.fq -db lambda_virus

bcftools mpileup -Ou --threads 12 -f xylo.fna ERR3201375.aln.bam | bcftools call --threads 12 -mv -Ob -o ERR3201375.calls.bcf

java -jar ~/workplace/downloads/picard.jar AddOrReplaceReadGroups \
    INPUT=dedup.bam \
    OUTPUT=rg.bam \
    SORT_ORDER=coordinate \
    RGLB=10678_0001 \
    RGPL=ILLUMINA \
    RGPU=ERR3201409 \
    RGSM=SAMEA5383982 \
    RGID=ERR3201409 \
    RGCN=PRJEB31266 \
    RGDS="Cattle rumen microbiome. breed:Aax; sex:M; age:589" \
    CREATE_INDEX=True

paste <(samtools view -h aln.sam | awk '!/^@/ { print $1 }') <(awk 'NR%4==0 { print $0 }' gen.fastq) | awk '{ printf("%s\tZQ:Z:%s\n", $1, $2) }' | samtools view -bS - >aln_with_quality.bam

magicblast -out aln.sam -num_threads 12 -outfmt sam -no_aligned -db db/xylo

magicblast -query ERR3201409.fastq -infmt fastq -db MF045434.1 -out ERR3201409.aln.sam -num_threads

efetch -db nucleotide -format fasta \
    -id MK690391.1,MF045437.1,MF045436.1,MF045435.1,MF045434.1,AB540108.1,KU668557.1,KU668556.1,KU668555.1,KU170775.1,KU170773.1,KU170772.1,KU170771.1,KU170770.1,KU170769.1,KU170768.1,KU170774.1,KU170767.1,KR611316.1,KR611315.1,FJ804147.1,KC818626.1,JN631038.1,AF278715.1,JF523203.1,FJ763639.2,HM769331.1,L36993.1,DQ665826.1,AY854501.1 \
    >xylo.fasta

makeblastdb -parse_seqids -dbtype nucl -in xylo.fasta -out db/xylo -title "XYLOSE REDUCTASE"

magicblast -sra_batch sra -db db/xylo -out out/aln.sam -num_threads 12 -outfmt sam -no_aligned

fasterq-dump ../ERR3201409 -o ERR3201409.fastq -p -e 12 -v -L debug

magicblast -sra ERX3229027 -db db/xylo -out out/aln.sam \
    -num_threads 12 -outfmt sam \
    -word_size 12 -penalty -2 -gapopen 5 -gapextend 2 -limit_lookup F

bcftools mpileup -Ou --threads 12 -f xylo.fasta sort.bam | bcftools call --threads 12 -mv -Ob -o calls.bcf
