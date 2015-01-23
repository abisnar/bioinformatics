from Bio import SeqIO

def read_fasta(fasta):
    handle = open(fasta)
    records = list(SeqIO.parse(handle,"fasta"))
    handle.close()
    return[(rec.id, str(rec.seq)) for rec in records]
