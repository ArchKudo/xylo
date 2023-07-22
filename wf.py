from Bio import Entrez
from Bio import SeqIO

Entrez.email = "your_email@example.com"

search = Entrez.esearch(db="nucleotide", term="xylose reductase (xyl1) gene, complete cds", retmax=200)
records = Entrez.read(search)
ids = records["IdList"]
