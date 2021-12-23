#%%
from multiprocessing import Pool
from pathlib import Path
import argparse
import inspect
import sys

import pandas as pd

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from EDA.gff3 import read_gff, attributes_to_dict

parser = argparse.ArgumentParser()
# general args
parser.add_argument("--gff", required=True,
                    help="GFF3 file.")
parser.add_argument("--data_table", required=True,
                    help="A CSV file such as SRA run table. "
                         "Expected to have 'Run', 'Tissue' & 'Group' columns (all other are ignored)."
                         "See alternate columns labels below.")
parser.add_argument("--run_col", default="Run",
                    help="Label of samples' col in data_table")
parser.add_argument("--tissue_col", default="Tissue",
                    help="Label of tissues' col in data_table")
parser.add_argument("--group_col", default="Group",
                    help="Label of groups' col in data_table")
parser.add_argument("--data_table_sep", default=",",
                    help="data_table's separator.")
parser.add_argument("--processes", default=20, type=int)
parser.add_argument("--id_prefix", default="",
                    help="For example, if 'ID' in the attributes column of the gff file is 'ID=gene-LOC107382895', "
                         "then id_prefix should be 'gene-'.")
parser.add_argument("--out_dir", required=True)
# make_gene_exp_file args
parser.add_argument("--salmon_dir", required=True, type=Path,
                    help="Salmon's main dir, containing a separate sub-dir for each sample.")
# make_groups_file
parser.add_argument("--alignment_dir", required=True, type=Path,
                    help="Dir with BAM files.")
parser.add_argument("--alignment_postfix", default="Aligned.sortedByCoord.out.bam")
# make_refseq_file args
parser.add_argument("--common_as_refseq_name", action="store_true",
                    help="As the gff_file doesn't have a 'Name' attribute, use its ID for both cols")


# parse args
args = parser.parse_args()

gff = args.gff
data_table = args.data_table
run_col = args.run_col
tissue_col = args.tissue_col
group_col = args.group_col
data_table_sep = args.data_table_sep
processes = args.processes
id_prefix = args.id_prefix
out_dir = Path(args.out_dir)
salmon_dir = args.salmon_dir
gene_exp_table = Path(out_dir, "IndexExpression.bed")  # gene expression out file
alignment_dir = args.alignment_dir
alignment_postfix = args.alignment_postfix
groups_table = Path(out_dir, "IndexGroups.csv")  # groups out file
common_as_refseq_name = args.common_as_refseq_name
refseq_table = Path(out_dir, "IndexRefseq.bed")  # refseq out file

gff_df = read_gff(gff)
data_table_df = pd.read_csv(data_table, sep=data_table_sep, usecols=[run_col, tissue_col, group_col])
id_prefix = f"={id_prefix}"


#%%

# run make_gene_exp_file.py

samples_dirs = [d for d in salmon_dir.iterdir() if d.is_dir() and d.name not in ["logs", "flags"]]
samples_dfs = [pd.read_csv(Path(sample_dir, "quant.genes.sf"),
                           sep="\t",
                           usecols=["Name", "TPM"])
               for sample_dir in samples_dirs]

tissues = data_table_df[tissue_col].unique().tolist()
if len(tissues) >= 2:  # 2 tissues or more
    samples_names = [sample_sub_dir.name for sample_sub_dir in samples_dirs]
    tissues.sort()
    print("\nOrder of tissues:\n", tissues, "\n\n")
    tissues_averaged_dfs = []
    for tissue in tissues:
        tissue_dfs = []
        for sample_df, sample_name in zip(samples_dfs, samples_names):
            if (data_table_df.loc[data_table_df[run_col] == sample_name, tissue_col] == tissue).bool():
                tissue_dfs.append(sample_df)
        tissue_averaged_df = pd.concat(tissue_dfs).groupby("Name").mean().reset_index()
        tissue_averaged_df = tissue_averaged_df.rename(columns={"TPM": tissue})
        tissues_averaged_dfs.append(tissue_averaged_df)

    final_averaged_df = tissues_averaged_dfs[0]

    for tissue_averaged_df in tissues_averaged_dfs[1:]:
        tissue = tissue_averaged_df.columns[1]
        final_averaged_df.insert(len(final_averaged_df.columns), tissue, tissue_averaged_df.loc[:, tissue])

else:
    final_averaged_df = pd.concat(samples_dfs).groupby("Name").mean().reset_index()


gff_genes_df = gff_df.loc[gff_df["type"] == "gene"]
id_col = gff_genes_df["attributes"].str.split(";").str[0].str.split(id_prefix).str[1]
gff_genes_df.insert(len(gff_genes_df.columns), "ID", id_col)

genes_missing_from_salmon = gff_genes_df.loc[~gff_genes_df["ID"].isin(final_averaged_df["Name"])]
if 0 < len(genes_missing_from_salmon):
    print(f"WARNING: {len(genes_missing_from_salmon):,} out of the {len(gff_genes_df):,} genes in the gff file "
          f"are not present in the genes quantified by Salmon.\n"
          f"See these genes at $out_dir/gff_genes_missing_from_salmon.txt")

# with Path(out_dir, "gff_genes_missing_from_salmon.txt").open("w") as missing_genes_file:
#     for gene in genes_missing_from_salmon["Name"].unique():
genes_missing_from_salmon["ID"].\
    reset_index(drop=True).\
    to_csv(Path(out_dir, "gff_genes_missing_from_salmon.txt"), header=False, sep="\t")


genes_present_in_salmon = gff_genes_df.loc[gff_genes_df["ID"].isin(final_averaged_df["Name"])]
genes_present_in_salmon = genes_present_in_salmon.sort_values("ID").reset_index(drop=True)
final_gff_genes = genes_present_in_salmon.loc[:, ["chrom", "start", "end", "ID", "strand"]]

final_averaged_df = final_averaged_df.sort_values("Name").reset_index(drop=True)
final_averaged_df = final_averaged_df.loc[final_averaged_df["Name"].isin(gff_genes_df["ID"])]

gene_exp_df = final_gff_genes
gene_exp_df["start"] = gene_exp_df["start"] - 1
if len(final_averaged_df) == 2:  # only one tissue
    tpm_col = final_averaged_df.iloc[:, 1]
else:
    tpm_col = final_averaged_df.iloc[:, 1:].apply(lambda x: ",".join(x.astype(str)), axis=1)
gene_exp_df.insert(4, "TPMs", tpm_col)

gene_exp_df = gene_exp_df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)

gene_exp_df.to_csv(gene_exp_table, sep="\t", header=False, index=False)


#%%

# run make_groups_file.py

groups_table_df = data_table_df.loc[:, [run_col, group_col]]
groups_table_df = groups_table_df.rename(columns={run_col: "Sample", group_col: "Group"})
groups_table_df = groups_table_df.sort_values("Sample").reset_index(drop=True)

bams = [bam for bam in alignment_dir.glob("**/*.bam") if bam.name.endswith(alignment_postfix)]

def sort_bams(bam):
    try:
        return int(bam.name.rstrip(alignment_postfix))
    except ValueError:  # probably due to an "SRR" prefix of files downloaded from the SRA
        return bam.name.rstrip(alignment_postfix)

sorted_bams = sorted(bams, key=lambda x: sort_bams(x))
bams_paths = [f"{bam}*" for bam in sorted_bams]

groups_table_df.insert(1, "SamplePath", bams_paths)

groups_table_df.to_csv(groups_table, sep=",", index=False)


#%%

# run make_refseq_file.py

gff_df = gff_df.loc[:, ["chrom", "type", "start", "end", "strand", "attributes"]]
gff_df["start"] = gff_df["start"] - 1

def sub_gff_dfs_by_genes(gff_df):
    genes_indices = gff_df.loc[gff_df["type"] == "gene"].index
    sub_gff_dfs = []
    for i in range(len(genes_indices)):
        start = genes_indices[i]
        if i < len(genes_indices) - 1:
            end = genes_indices[i+1]
        else:
            end = len(gff_df)
        sub_gff_df = gff_df.iloc[start:end]
        sub_gff_dfs.append(sub_gff_df)
    return sub_gff_dfs

sub_gff_dfs = sub_gff_dfs_by_genes(gff_df)

gene_id_prefix = args.id_prefix

def extract_refseq_and_common_names(attributes):
    attributes = attributes_to_dict(attributes)
    refseq_name = attributes["ID"].lstrip(gene_id_prefix)
    common_name = refseq_name if common_as_refseq_name else attributes["Name"]
    return refseq_name, common_name

def sub_gff_df_to_sub_refseq_tuple(sub_gff_df):
    gene_line = sub_gff_df.loc[sub_gff_df["type"] == "gene"].squeeze()
    chrom = gene_line["chrom"]
    gene_start, gene_end = gene_line["start"], gene_line["end"]
    gene_attributes = sub_gff_df.loc[sub_gff_df["type"] == "gene", "attributes"].values[0]
    gene_refseq_name, gene_common_name = extract_refseq_and_common_names(gene_attributes)
    strand = gene_line["strand"]
    exons_df = sub_gff_df.loc[sub_gff_df["type"] == "exon"]
    exons_starts = exons_df["start"].astype(str).str.cat(sep=",")
    exons_ends = exons_df["end"].astype(str).str.cat(sep=",")
    sub_ref_tuple = (chrom, gene_start, gene_end, gene_refseq_name, gene_common_name, strand, exons_starts, exons_ends)
    return sub_ref_tuple

with Pool(processes=processes) as pool:
    sub_refseq_tuples = pool.map(func=sub_gff_df_to_sub_refseq_tuple, iterable=sub_gff_dfs)

refseq_cols = ["chrom", "start", "end", "refseq_name", "common_name", "strand", "exons_starts", "exons_ends"]
refseq_df = pd.DataFrame.from_records(sub_refseq_tuples, columns=refseq_cols)

# filter out genes from refseq_df that aren't present in gene_exp_df
refseq_df = refseq_df.loc[refseq_df["refseq_name"].isin(gene_exp_df["ID"])]

refseq_df = refseq_df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)

refseq_df.to_csv(refseq_table, sep="\t", header=False, index=False)
