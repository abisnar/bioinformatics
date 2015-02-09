def get_seqs_from_records(records):
    return [record[1] for record in records]

def get_seq_id_from_records(records):
    return[record[0] for record in records]