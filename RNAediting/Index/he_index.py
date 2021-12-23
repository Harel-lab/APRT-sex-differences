from pathlib import Path
import subprocess
import argparse
import inspect
import sys

from pybedtools import BedTool

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from EditingUtils.summarize_utils import execute_notebook
from EditingUtils.seq_reader import GenomeReader
from General.consts import final_words
from Index.index import index


def main(*, root_dir, he_cluster_dir, bam_files_suffix, follow_links, stranded_he, stranded_index, paired_end,
         genome_file, main_out_dir, sample_threads, sample_strands_threads, distances, gff, data_table, run_col,
         tissue_col, main_group_col, data_table_sep, processes, index_path, salmon_dir, common_as_refseq_name,
         id_prefix):

    # 1 - set input paths

    root_dir = Path(root_dir).absolute()
    genome_file = Path(genome_file).absolute()
    gff = Path(gff).absolute()
    data_table = Path(data_table).absolute()
    salmon_dir = Path(salmon_dir).absolute()

    # 2 - set main output path

    main_out_dir = Path(main_out_dir).absolute()
    main_out_dir.mkdir(exist_ok=True)

    # 3 - prepare files needed for running the Index

    # 3.1 - src_dir stores files needed for the Index

    src_dir = Path(main_out_dir, "SrcFiles")
    src_dir.mkdir(exist_ok=True)

    # 3.2 - extend ue regions
    # (they will use as different inputs for the different runs of the Index)

    # find original A2G UE regions file from the best cluster, as found previously by Cluster Screening
    stranded_prefix = "stranded." if stranded_he else ""  # set prefix to find ue_regions_file in he_cluster_dir
    ue_regions_file = Path(he_cluster_dir, "UERegions", f"UEmergeSd20.{stranded_prefix}A2G.bed").absolute()
    # create the extended regions files, and return their paths
    extended_ue_dir = Path(src_dir, "ExtendedUERegions")
    extended_ue_dir.mkdir(exist_ok=True)
    extended_ue_regions_files = extend_regions_by_different_distances(extended_ue_dir, ue_regions_file, distances,
                                                                      stranded_he, genome_file)

    # 3.3 - create the gene_expression_file, refseq_file, and groups_file
    # (they will serve each Index run in the same manner)

    prepare_src_files(gff, data_table, run_col, tissue_col, main_group_col, data_table_sep, processes, src_dir,
                      salmon_dir, root_dir, bam_files_suffix, common_as_refseq_name, id_prefix)
    gene_expression_file = Path(src_dir, "IndexExpression.bed")
    refseq_file = Path(src_dir, "IndexRefseq.bed")
    groups_file = Path(src_dir, "IndexGroups.csv")

    # 4 - define mismatches for the Index command

    index_cmd_mismatches = "AllMismatches" if stranded_index else "UnstrandedMismatches"

    # 5 - run Index for each distance

    for extended_ue_regions_file, distance in zip(extended_ue_regions_files, distances):
        sub_out_dir = Path(main_out_dir, f"HEIndex{distance}")
        sub_out_dir.mkdir(exist_ok=True)
        index(root_dir, bam_files_suffix, follow_links, stranded_index, paired_end, gene_expression_file, refseq_file,
              groups_file, genome_file, extended_ue_regions_file, sub_out_dir, index_cmd_mismatches, sample_threads,
              sample_strands_threads, index_path)


def extend_regions_by_different_distances(extended_ue_dir, ue_regions_file, distances, stranded_he, genome_file):
    """
    Create subsets of extended UE regions.

    For each distance, extend each region with $distance bp both up- and downstream.
    Save the sets of extended regions in extended_ue_dir. Also return their paths in a list.

    @param extended_ue_dir: where to save the extended sets of UE regions
    @type extended_ue_dir: Path
    @param ue_regions_file: original set of regions (described by distance = 0)
    @type ue_regions_file: Path
    @param distances: bp to extend regions by
    @type distances: list[int]
    @param stranded_he: whether the Hyper Editing script was run as stranded
    @type stranded_he: bool
    @param genome_file: path to genome fasta file
    @type genome_file: str or Path
    @return: extended_ue_regions_files
    @rtype: list[Path]
    """
    genome = GenomeReader.parse_genome(genome_file)
    chromosomal_lengths = {chromosome: len(seq) for chromosome, seq in genome.items()}
    extended_ue_regions_files = []
    for distance in distances:
        if stranded_he:
            extended_ue_regions_file = Path(extended_ue_dir, f"UE.extended{distance}.merged.stranded.A2G.bed")
        else:
            extended_ue_regions_file = Path(extended_ue_dir, f"UE.extended{distance}.merged.A2G.bed")
        extend_regions_by_distance(ue_regions_file, distance, chromosomal_lengths, extended_ue_regions_file)
        extended_ue_regions_files.append(extended_ue_regions_file)
    return extended_ue_regions_files


def extend_regions_by_distance(ue_regions_file, distance, chromosomal_lengths, out_path):

    def extend(feature):
        start = feature.start
        extended_start = start - distance
        extended_start = extended_start if 0 <= extended_start else 0
        feature.start = extended_start

        chromosome_len = chromosomal_lengths[feature.chrom]
        end = feature.end
        extended_end = end + distance
        extended_end = extended_end if extended_end <= chromosome_len else chromosome_len
        feature.end = extended_end

        return feature

    bed = BedTool(ue_regions_file).\
        sort().\
        each(extend).\
        merge(d=20, c=[4, 5, 6], o=["count", "count", "distinct"]).\
        sort()

    # save only first 3 cols, as needed for the index
    BedTool([(x.chrom, x.start, x.end) for x in bed]).\
        saveas(out_path)


def prepare_src_files(gff, data_table, run_col, tissue_col, group_col, data_table_sep, processes, src_dir, salmon_dir,
                      root_dir, bam_files_suffix, common_as_refseq_name, id_prefix):
    cmd = f"python Index/prepare_index_files.py " \
          f"--gff {gff} " \
          f"--data_table {data_table} " \
          f"--run_col {run_col} " \
          f"--tissue_col {tissue_col} " \
          f"--group_col {group_col} " \
          f"--data_table_sep {data_table_sep} " \
          f"--processes {processes} " \
          f"--out_dir {src_dir} " \
          f"--salmon_dir {salmon_dir} " \
          f"--alignment_dir {root_dir} " \
          f"--alignment_postfix {bam_files_suffix}"
    if common_as_refseq_name:
        cmd += " --common_as_refseq_name "
    if id_prefix:
        cmd += f" --id_prefix {id_prefix}"
    subprocess.run(cmd, shell=True)


def summarize(*, out_dir,

              global_notebook_template, out_global_notebook_name,
              global_notebook_title, global_notebook_img, global_notebook_img_src, groups_order,

              pca_notebook_template, out_pca_notebook_name,
              he_cluster_dir, stranded_he, stranded_index, sub_groups, gff,
              ):
    out_dir = Path(out_dir).absolute()

    # execute_notebook(template_notebook, execution_out_dir, executed_notebook_name, **parameters)

    # 1 - global notebook
    global_notebook_template = Path(global_notebook_template).absolute()
    he_index_main_dir = str(out_dir)
    global_notebook_title = global_notebook_title if global_notebook_title else out_global_notebook_name
    execute_notebook(global_notebook_template, out_dir, out_global_notebook_name,
                     he_index_main_dir=he_index_main_dir, notebook_title=global_notebook_title,
                     notebook_img=global_notebook_img, notebook_img_src=global_notebook_img_src,
                     groups_order=groups_order)

    # 2 - pca notebook
    pca_notebook_template = Path(pca_notebook_template).absolute()
    gff = str(Path(gff).absolute())
    groups_file = str(Path(out_dir, "SrcFiles", "IndexGroups.csv"))
    stranded_prefix = "stranded." if stranded_he else ""  # set prefix to find ue_regions_file in he_cluster_dir
    ue_regions_file = str(Path(Path(he_cluster_dir).absolute(),
                               "UERegions",
                               f"UEmergeSd20.{stranded_prefix}A2G.bed").absolute())
    best_index_file_flag = [f for f in out_dir.iterdir() if f.name.startswith("BEST")]
    assert len(best_index_file_flag) == 1  # there should be one, and one only, best distance to extend UE regions with
    best_index_distance = best_index_file_flag[0].name.split("_")[1]   # f"BEST_{best_distance}_BP"
    best_index_dir = str(Path(out_dir, f"HEIndex{best_index_distance}"))
    execute_notebook(pca_notebook_template, out_dir, out_pca_notebook_name,
                     index_dir=best_index_dir, groups_file=groups_file, groups_cols=sub_groups, stranded=stranded_index,
                     gff_file=gff, ue_regions_file=ue_regions_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--he_cluster_dir",
                        required=True,
                        help="Best cluster after Cluster Screening")
    parser.add_argument("--distances",
                        default=[0, 30, 60, 90, 120],
                        help="One or more distances that define extension of UE regions both up- and downstream "
                             "(0 is the original file).")

    parser.add_argument("-d", "--root_dir",
                        required=True,
                        help="Input dir with BAM files.")
    parser.add_argument("-f", "--bam_files_suffix",
                        help="A suffix of the BAM files to run on (e.g. Sorted.By.Coord.bam). Should be the *full* "
                             "suffix.",
                        default=".Aligned.SortedByCoords.bam")
    parser.add_argument("--follow_links",
                        action="store_true")
    parser.add_argument("--stranded_he",
                        action="store_true",
                        help="Whether the Hyper Editing & Cluster Screening scripts were run as stranded.")
    parser.add_argument("--stranded_index",
                        action="store_true",
                        help="Whether the Index should run as stranded.")
    parser.add_argument("--paired_end",
                        action="store_true")
    parser.add_argument("--data_table",
                        required=True,
                        help="A CSV file such as SRA run table. "
                             "Expected to have 'Run', 'Tissue' & 'Group' columns (all other are ignored)."
                             "See alternate columns labels below.")
    parser.add_argument("--data_table_sep",
                        default=",",
                        help="data_table's separator.")
    parser.add_argument("--run_col",
                        default="Run",
                        help="Label of samples' col in data_table")
    parser.add_argument("--tissue_col",
                        default="Tissue",
                        help="Label of tissues' col in data_table")
    parser.add_argument("--main_group_col",
                        default="Group",
                        help="Label of groups' col in `data_table`. The values correspond to the order of "
                             "`sub_groups` (see next argument).")
    parser.add_argument("--sub_groups",
                        nargs="+",
                        default=["Age(weeks)", "Genotype", "Diet", "Sex"])
    parser.add_argument("--groups_order",
                        nargs="+",
                        default=["6.5-WT-Full-Male",
                                 "6.5-Het-Full-Male",
                                 "15-WT-Full-Male",
                                 "15-Het-Full-Male",

                                 "6.5-WT-Fasted-Male",
                                 "6.5-Het-Fasted-Male",
                                 "15-WT-Fasted-Male",
                                 "15-Het-Fasted-Male",

                                 "6.5-WT-Fasted-Female",
                                 "6.5-Het-Fasted-Female",
                                 "15-WT-Fasted-Female",
                                 "15-Het-Fasted-Female"],
                        help="Sorting order of `sub_groups`' values for some basic plots in executed notebooks.")

    parser.add_argument("--salmon_dir",
                        required=True,
                        help="Salmon's main dir, containing a separate sub-dir for each sample.")

    parser.add_argument("-gf", "--genome_file",
                        required=True)
    parser.add_argument("--gff",
                        required=True,
                        help="A GFF3 file.")
    parser.add_argument("--id_prefix", default="",
                        help="For example, if 'ID' in the attributes column of the gff file is 'ID=gene-LOC107382895', "
                             "then id_prefix should be 'gene-'.")
    parser.add_argument("--common_as_refseq_name", action="store_true",
                        help="If the gff_file doesn't have a 'Name' attribute, use its ID for both cols")

    parser.add_argument("-o", "--out_dir",
                        required=True,
                        help="The root directory for the cmpileups outputs and other temporary files.Outputs (will "
                             "create sub-directories per sample according to their original directories tree from the "
                             "root input directory to the BAM file). Files are deleted at the end of the run (except "
                             "for cmpileups if otherwise specified).")

    parser.add_argument("--index_path",
                        default="RNAEditingIndex1.1")
    parser.add_argument("--processes", default=20)
    parser.add_argument("-ts", "--sample_threads", default=10,
                        help="The number of samples to process in parallel.")
    parser.add_argument("-tsd", "--sample_strands_threads", default=50,
                        help="The maximal strand decisions per sample to process in parallel.")

    parser.add_argument("--global_notebook_template",
                        default="Notebooks/HEIndexAnalysis.ipynb")
    parser.add_argument("--out_global_notebook_name",
                        default="HEIndexAnalysis")
    parser.add_argument("--global_notebook_title",
                        help="Defaults to `out_global_notebook_name` if None.")
    parser.add_argument("--global_notebook_img",
                        default="https://upload.wikimedia.org/wikipedia/commons/2/22/Nothobranchius_furzeri_GRZ_thumb.jpg")
    parser.add_argument("--global_notebook_img_src",
                        default="https://de.wikipedia.org/wiki/Benutzer:Ugau")

    parser.add_argument("--pca_notebook_template",
                        default="Notebooks/PerRegionPerSamplePCA.ipynb")
    parser.add_argument("--out_pca_notebook_name",
                        default="PerRegionPerSamplePCA")

    parser.add_argument("--re_summarize_only",
                        help="Don't calculate the Index as it already exists - only run summarize again.",
                        action="store_true")

    args = parser.parse_args()

    if not args.re_summarize_only:
        main(root_dir=args.root_dir, he_cluster_dir=args.he_cluster_dir, bam_files_suffix=args.bam_files_suffix,
             follow_links=args.follow_links, stranded_he=args.stranded_he, stranded_index=args.stranded_index,
             paired_end=args.paired_end, genome_file=args.genome_file, main_out_dir=args.out_dir,
             sample_threads=args.sample_threads, sample_strands_threads=args.sample_strands_threads,
             distances=args.distances, gff=args.gff, data_table=args.data_table, run_col=args.run_col,
             tissue_col=args.tissue_col, main_group_col=args.main_group_col, data_table_sep=args.data_table_sep,
             processes=args.processes, salmon_dir=args.salmon_dir, common_as_refseq_name=args.common_as_refseq_name,
             id_prefix=args.id_prefix, index_path=args.index_path)


    summarize(out_dir=args.out_dir,

              global_notebook_template=args.global_notebook_template,
              out_global_notebook_name=args.out_global_notebook_name,
              global_notebook_title=args.global_notebook_title, global_notebook_img=args.global_notebook_img,
              global_notebook_img_src=args.global_notebook_img_src, groups_order=args.groups_order,

              pca_notebook_template=args.pca_notebook_template, out_pca_notebook_name=args.out_pca_notebook_name,

              he_cluster_dir=args.he_cluster_dir, stranded_he=args.stranded_he, stranded_index=args.stranded_index,
              sub_groups=args.sub_groups, gff=args.gff)

    final_words()
