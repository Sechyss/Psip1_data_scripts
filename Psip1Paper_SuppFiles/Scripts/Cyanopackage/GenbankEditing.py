import os
import random
import sys

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


class GenbankEditing:
    """""This file contains methods to manipulate genbank files for bacteria"""

    def __init__(self, genbank):
        self.genbank = genbank

    def protein_fasta_creation(self, outfile):
        gbank = SeqIO.parse(self.genbank, 'genbank')  # Open the genbank file
        for genome in gbank:
            for feature in genome.features:
                if feature.type == "CDS":  # Find CDS to collect the information
                    locustag = str(feature.qualifiers["locus_tag"][0]).replace(' ', '_')
                    ntseq = feature.extract(genome.seq)
                    aaseq = ntseq.translate(table=11, cds=False, stop_symbol="")
                    fasta_aa = SeqRecord(aaseq, locustag, description='')
                    SeqIO.write(fasta_aa, outfile, 'fasta')

        return outfile

    def extract_sequence(self, outfile):  # Extract the sequence of the genome simply
        gbank = SeqIO.parse(self.genbank, 'genbank')
        for genome in gbank:
            ntseq = genome.seq
            fasta_nt = SeqRecord(ntseq, id=str(self.genbank).replace('.gbk', '.fna'), description='')
            SeqIO.write(fasta_nt, outfile, "fasta")

        return outfile

    def extraction_feature_fasta(self, outfile, feature):  # Extract given feature name
        gbank = SeqIO.parse(self.genbank, 'genbank')
        features_to_extract = []
        for genome in gbank:
            for i in genome.features:
                # Collection of features in  list for each rRNA.
                if i.type == str(feature):
                    locustag = str(i.qualifiers["locus_tag"][0])
                    seq = feature.extract(genome.seq)
                    record = SeqRecord(seq, id=locustag, description="")
                    features_to_extract.append(record)

        SeqIO.write(features_to_extract, outfile, "fasta")  # Write everything in one single file.
        return outfile

    def genbank_revamp(self, outfile):
        genbank = SeqIO.parse(self.genbank, "genbank")
        locustag = []
        stored_random = random.randint(1, 9999)
        for genome in genbank:
            final_feature = []
            for feature in genome.features:
                if feature.type == "CDS":
                    try:
                        locustag = feature.qualifiers["locus_tag"][0]
                        final_feature.append(feature)
                    finally:
                        stored_random += 1
                        feature.qualifiers["locus_tag"] = str(locustag[:5]) + "_" + str(stored_random) + "m"
                        final_feature.append(feature)
                elif feature.type == 'rRNA':
                    try:
                        locustag = feature.qualifiers["locus_tag"][0]
                        final_feature.append(feature)
                    finally:
                        stored_random += 1
                        feature.qualifiers["locus_tag"] = str(locustag[:5]) + "_" + str(stored_random) + "m"
                        final_feature.append(feature)
                elif feature.type == 'tRNA':
                    try:
                        locustag = feature.qualifiers["locus_tag"][0]
                        final_feature.append(feature)
                    finally:
                        stored_random += 1
                        feature.qualifiers["locus_tag"] = str(locustag[:5]) + "_" + str(stored_random) + "m"
                        final_feature.append(feature)
                elif feature.type == "source":
                    final_feature.append(feature)
            genome.features = final_feature
            with open(outfile, 'a') as newgbk:
                SeqIO.write(genome, newgbk, "genbank")

        return outfile

    def gene_length_dict(self, outdictionary):
        genbank = SeqIO.parse(self.genbank, "genbank")
        for genome in genbank:
            for feature in genome.features:
                if feature.type == "CDS" or feature.type == "rRNA" or feature.type == "tRNA":
                    locustag = str(feature.qualifiers["locus_tag"][0])
                    length = len(str(feature.extract(genome.seq)))
                    outdictionary.update({locustag: length})

    def genome_extraction_from_list(self, inputlist, outfile):
        genbank = SeqIO.parse(self.genbank, "genbank")
        for genome in genbank:
            if genome.name in inputlist:
                aa = SeqRecord(genome.seq, id=genome.id, description='')
                SeqIO.write(aa, outfile, "fasta")

        return outfile


def referenceseqs_single_feature(indirectory, outdirectory, featuretype):
    os.chdir(indirectory)
    for file in tqdm(os.listdir()):
        if file.endswith('.gbk'):
            for genome in SeqIO.parse(file, 'genbank'):
                for feature in genome.features:
                    if feature.type == featuretype:
                        try:
                            locustag = str(feature.qualifiers["locus_tag"][0]).replace(' ', '')
                            ntseq = feature.extract(genome.seq)
                            aa = SeqRecord(ntseq, id=locustag, description='')
                            with open(outdirectory + locustag + '.fna', 'w') as DB:
                                SeqIO.write(aa, DB, "fasta")
                        except NameError:
                            sys.stderr.write('Error: Check your genomes contain all its features')
                            continue


def extract_genome_sequences(genbank_file, fasta_file):
    with open(genbank_file, 'r') as genbank, open(fasta_file, 'w') as fasta:
        for record in SeqIO.parse(genbank, 'genbank'):
            genome_name = record.name
            genome_sequence = str(record.seq)

            fasta.write(f'>{genome_name}\n{genome_sequence}\n')

    print(f'FASTA file "{fasta_file}" has been created.')

