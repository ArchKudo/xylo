# Code for MSc Thesis on ``Comparing and contrasting the protein embedding landscape with commercially important enzymes haplotypes recovered from meta-genomes." 

## Usage

```sh
git clone git@github.com:ArchKudo/xylo.git

cd xylo/

sbatch scripts/aln.sh
```

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
