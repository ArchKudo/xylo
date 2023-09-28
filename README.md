# Code for MSc Thesis on ``Comparing and contrasting the protein embedding landscape with commercially important enzymes haplotypes recovered from meta-genomes." 

## Usage

### Aligning files

```sh
git clone git@github.com:ArchKudo/xylo.git

cd xylo/

sbatch scripts/aln.sh

# bash scripts/aln.sh # locally
```

### Recovering haplotypes


#### On Slurm

```sh
# pyVCF anaconda version is broken better to run this locally instead without python
anaconda3-launch conda create -n "xylo3.7.16" python=3.7.16
anaconda3-launch --env xylo3.7.16 pip install hanselx pysam PyVCF gretel
sbatch scripts/recover.sh
```

#### Local install

```sh







## File structure

- db/ - Contains the bowtie index for reference files, in git just for backup
- report/ - Source for the latex report

- scripts/aln.sh - The main workflow used for aligning reads in slurm
- scripts/snpper.py - Stolen from gretel-test repository with breaking deprecated api fix
- scripts/snippets.sh - Snippets for adhoc data analysis

- .tool-versions - Required python version
- requirements.dev - Required python developer dependencies

The below two files are required for the aln.sh workflow:
- runs - Contains the SRA accessions for the metagenomes reads to align to
- xylo.fasta - The pseudoreference fasta sequence
