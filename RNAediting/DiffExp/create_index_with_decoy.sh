#!/usr/bin/env bash

# create argparse
ARGPARSE_DESCRIPTION="Create salmon index. Notice it is best to keep all files in "     # this is optional
source /private/common/Software/BashLibs/argparse-bash/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('-g', '--genome', type=str, help='Path of genome reference file - chromosomes (fasta)',required=True)
parser.add_argument('-t', '--transcriptome', type=str, help='Path of transcriptome file - transcripts (fasta). Transcript ID should match the transcript-to-geneID file')
parser.add_argument('-ucg', '--uncompress_fasta_genome',type=str, help='uncompress genome fasta command',default="zcat")
parser.add_argument('-uct', '--uncompress_fasta_transcriptome',type=str, help='uncompress transcriptome fasta command',default="zcat")
parser.add_argument('-n', '--name_genome_index', type=str, help='Genome name, to create relevant directory and output index',required=True)
parser.add_argument('-sid', '--salmon_index_dir',type=str, help='Output directory', required=True)
EOF

# print all commands
set -x

# create directory
output_dir=${SALMON_INDEX_DIR}/${NAME_GENOME_INDEX}
echo -e $output_dir
[ ! -d $output_dir ] && mkdir -p $output_dir

# Salmon indexing requires the names of the genome targets, which is extractable by using the grep command
grep "^>" <($UNCOMPRESS_FASTA_GENOME $GENOME) | cut -d " " -f 1 | sed -e 's/>//g' > $output_dir/decoys.txt

# Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index.
# NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
cat <($UNCOMPRESS_FASTA_TRANSCRIPTOME $TRANSCRIPTOME) <($UNCOMPRESS_FASTA_GENOME $GENOME) > $output_dir/gentrome.fa.gz

# index
#/private/common/Software/salmon/salmon1.4.0/salmon-latest_linux_x86_64/bin/salmon index -t $output_dir/gentrome.fa.gz -d $output_dir/decoys.txt -i $output_dir --gencode -p 40
salmon index -t $output_dir/gentrome.fa.gz -d $output_dir/decoys.txt -i $output_dir --gencode -p 40

#grep "^>" <(gunzip -c GRCm38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 | sed -e 's/>//g' > decoys.txt
#cat gencode.vM23.transcripts.fa.gz GRCm38.primary_assembly.genome.fa.gz > gentrome.fa.gz
#salmon1.4.0/salmon-latest_linux_x86_64/bin/salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode -p 40

