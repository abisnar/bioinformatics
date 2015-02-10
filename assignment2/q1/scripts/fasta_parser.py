from Bio import SeqIO
import re

def read_fasta(fasta):
    handle = open(fasta)
    records = list(SeqIO.parse(handle,"fasta"))
    handle.close()
    return[(re.sub("^sp\|[\w\d]*\|", "", str(rec.name)), str(rec.seq)) for rec in records]

def get_alignment_seq(fasta):
	return [fasta[1] for f in fasta]
