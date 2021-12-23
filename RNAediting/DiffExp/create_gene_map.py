"""
Create a tab-delimited file where where each line contains the name of a transcript and the gene to which it belongs,
separated by a tab.
This is a file-format accepted by Salmon's '--geneMap' flag.


Example for given parameters:

transcriptome = "Genome/RefseqAnnotations/GCF_001465895.1_Nfu_20140520_rna.fna"
gff = "Genome/RefseqAnnotations/GCF_001465895.1_Nfu_20140520_genomic.gff"
id_prefix = "rna-"
parent_prefix = "gene-"
out_path = "Differential_Expression/MappingBasedSalmon/tx2gene.tab"
"""

from pathlib import Path
import argparse
import inspect
import sys

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from EditingUtils.seq_reader import GenomeReader
from EDA.gff3 import read_gff

parser = argparse.ArgumentParser()
parser.add_argument("--transcriptome", required=True, help="Transcriptome fasta file")
parser.add_argument("--gff", required=True, help="Gff file")
parser.add_argument("--id_prefix", default="",
                    help="For example, if 'ID' in the attributes column of the gff file is 'ID=rna-XM_015955279.1', "
                         "then id_prefix should be 'rna-'")
parser.add_argument("--parent_prefix", default="",
                    help="For example, if 'Parent' in the attributes column of the gff file is "
                         "'Parent=gene-LOC107382895', then parent_prefix should be 'gene-'")
parser.add_argument("--rna_terms", default=['mRNA', "transcript", "lnc_RNA"], nargs="*",
                    help="Terms of features to use from the gff file, as they appear in the 3rd column of the gff file")
parser.add_argument("--out_path", required=True)
args = parser.parse_args()

transcriptome = GenomeReader.parse_genome(args.transcriptome)
transcripts_names = [key for key in transcriptome]

gff = read_gff(args.gff)
gff_rna = gff.loc[gff["type"].isin(args.rna_terms)]
assert len(gff_rna) == len(transcripts_names)

id_prefix = f"={args.id_prefix}"
parent_prefix = f"={args.parent_prefix}"
id_col = gff_rna["attributes"].str.split(";").str[0].str.split(id_prefix).str[1]
parent_col = gff_rna["attributes"].str.split(";").str[1].str.split(parent_prefix).str[1]

gff_rna.insert(len(gff_rna.columns), "ID", id_col)
gff_rna.insert(len(gff_rna.columns), "Parent", parent_col)

out_df = gff_rna.loc[:, ["ID", "Parent"]]
out_df.to_csv(args.out_path, sep="\t", header=False, index=False)
