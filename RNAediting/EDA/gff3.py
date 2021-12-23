from pathlib import Path
import inspect
import sys
import os

import pandas as pd

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from General.consts import GFF_COLS


def read_gff(gff_file_path):
    gff_df = pd.read_csv(gff_file_path,
                         sep="\t",
                         names=GFF_COLS,
                         comment="#")
    return gff_df


def attributes_to_dict(attributes):
    """
    Parse a string of a gff file into a dictionary.

    @param attributes: attributes of a single feature (one line) in a gff file, as they appear in the 9th column
    @type attributes: str
    @return: attributes
    @rtype: dict
    """
    attributes = [attribute.split("=") for attribute in attributes.split(";")]
    attributes = {attribute[0]: attribute[1] for attribute in attributes}  # {k: v for (k, v) in attributes}
    return attributes


def types_terms_in_gff3_file(gff3_file):
    """Return a set of the terms found in the 'type' column (the 3rd column) of a GFF3 file.

    A GFF3 is a tab-delimited 9-cols file for storing genomic features.
    Documentation:
    (1)  https://www.ensembl.org/info/website/upload/gff3.html
    (2)  http://gmod.org/wiki/GFF3"""
    terms = set()
    with open(gff3_file, "r") as file:
        for line in file:
            if line.startswith("#"):   # skip comment line
                continue
            components = line.split("\t")
            terms.add(components[2])
    return terms


def create_gff_file_of_one_type(out_dir, type):
    gff_file = os.path.join(out_dir, f"{type}.gff")
    file = open(gff_file, "w")
    return file


def divide_to_files_by_type(gff3_file, out_dir):
    """Divide a gff3_file to sub files by their type (3rd col in a gff3 file), and write them to out_dir."""
    types = types_terms_in_gff3_file(gff3_file)
    file_type_handles = {type: create_gff_file_of_one_type(out_dir, type) for type in types}
    with open(gff3_file, "r") as file:
        for line in file:
            if line.startswith("#"):   # skip comment line
                continue
            components = line.split("\t")
            type = components[2]
            file_handle = file_type_handles[type]
            try:
                file_handle.write(line)
            except Exception as e:
                print(e)
    for file_handle in file_type_handles.values():
        file_handle.close()


def sub_gff_dfs_by_genes(gff_df: pd.DataFrame) -> list[pd.DataFrame]:
    """Divide gff_df to multiple dataframes, where each sub df contains one gene alongside its exons, introns, etc."""
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
