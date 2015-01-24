from Bio import SeqIO

def read_fasta(fasta):
    handle = open(fasta)
    records = list(SeqIO.parse(handle,"fasta"))
    handle.close()
    return[(rec.name, str(rec.seq)) for rec in records]

def get_alignment_seq(fasta):
	return [fasta[1] for f in fasta]
