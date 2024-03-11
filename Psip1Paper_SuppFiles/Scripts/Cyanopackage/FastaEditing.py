from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def extract_common_sequences(referencefile, infile, outfile):
    align = SeqIO.parse(referencefile, 'fasta')
    fasta_seqs = SeqIO.parse(infile, 'fasta')
    extract_list = []
    for sequence in align:
        extract_list.append(sequence.id)

    for sequence2 in fasta_seqs:
        if sequence2.id in extract_list:
            aa = SeqRecord(sequence2.seq, id=sequence2.id, description='')
            with open(outfile, 'a') as DB:
                SeqIO.write(aa, DB, "fasta")


def fasta_len(fastafile):
    FastaFile = open(fastafile, 'r')
    SeqLen = []
    for rec in SeqIO.parse(FastaFile, 'fasta'):
        name = rec.description
        seq = rec.seq
        ContigLen = int(len(seq))
        print(ContigLen)
        SeqLen.append(ContigLen)
    TotalLength = sum(SeqLen)
    print(TotalLength)
    FastaFile.close()


def edit_NCBI_ids(fastafile, outfile):
    fasta_seqs = SeqIO.parse(fastafile, 'fasta')
    features_to_extract = []
    for sequence in fasta_seqs:
        id = str(sequence.id)+'_'+str(sequence.description).split('[')[1].replace(']','').replace(' ','_')
        seq = sequence.seq
        record = SeqRecord(seq, id=id, description="")
        features_to_extract.append(record)

    SeqIO.write(features_to_extract, outfile, "fasta")  # Write everything in one single file.
    return outfile
