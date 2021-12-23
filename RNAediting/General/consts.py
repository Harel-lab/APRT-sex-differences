from datetime import datetime

def final_words():
    final_words_message = "A star shines on the hour of our scripting"
    print(f"\n\n\n{datetime.now()}\t\t{final_words_message}\n", end="")


# RNA_TO_AA_DICT = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L',
#                   'CUG': 'L', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V',
#                   'GUA': 'V', 'GUG': 'V', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P',
#                   'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
#                   'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '*',
#                   'UAG': '*', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N',
#                   'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'UGU': 'C',
#                   'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
#                   'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
#                   'GGG': 'G'}


# HE_TRANSFORMATIONS = ("A2G", "A2C", "A2T", "T2C", "T2G", "G2C")

GFF_COLS = ["chrom", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
# BED6_COLS = ["chrom", "start", "end", "name", "score", "strand"] # todo refactor BED6_COLS to bed_6_cols

BED_6_COLS = ["Chrom", "Start", "End", "Name", "Score", "Strand"]
BED_3_COLS = BED_6_COLS[:3]


STRANDED_MISMATCH_TYPES = {f"{x}2{y}" for x in "ATCG" for y in "ATCG" if x != y}
UNSTRANDED_MISMATCH_TYPES = {"A2G", "A2T", "G2A", "C2A", "A2C", "G2C"}
HE_UNSTRANDED_MISMATCH_TYPES = UNSTRANDED_MISMATCH_TYPES
INDEX_UNSTRANDED_MISMATCH_TYPES = {"A2G", "A2T", "C2T", "C2A", "A2C", "C2G"}
