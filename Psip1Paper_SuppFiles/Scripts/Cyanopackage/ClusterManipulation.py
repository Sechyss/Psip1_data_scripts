import os
import pickle

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def representative_sequence_extraction_fna(directory):
    os.chdir(directory)
    for file in os.listdir():
        if file.endswith('.fna'):
            db = open(file.replace('.$', '.fasta'), 'w')
            for seq_record in SeqIO.parse(file, 'fasta'):
                cluster = str(file.replace('.faa', ''))
                locustag = str(cluster)
                ntseq = seq_record.seq
                aa = SeqRecord(ntseq, id=locustag, description='')
                SeqIO.write(aa, db, "fasta")
                db.close()
                break


def representative_sequence_extraction_faa(directory):
    os.chdir(directory)
    for file in os.listdir():
        if file.endswith('.fna'):
            db = open(file.replace('.$', '.fasta'), 'w')
            for seq_record in SeqIO.parse(file, 'fasta'):
                cluster = str(file.replace('.faa', ''))
                locustag = str(cluster)
                ntseq = seq_record.seq
                aaseq = ntseq.translate(table=11, cds=False, stop_symbol="", gap=None)
                aa = SeqRecord(aaseq, id=locustag, description='')
                SeqIO.write(aa, db, "fasta")
                db.close()
                break


def translate_cluster_locustag_id(locustag: str, dictionary):
    dictionary1 = open(dictionary, "rb")
    clusters_gh = pickle.load(dictionary1)
    del dictionary1

    cluster_id = str(clusters_gh[locustag])

    return cluster_id
