from multiprocessing import Pool
from collections import Counter
from pathlib import Path
import argparse
import inspect
import sys

from pybedtools import BedTool

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from General.os_utils import copy_text_file
from ClusterScreening.summarize_cluster_screening import summarize_cluster_screening
from General.consts import STRANDED_MISMATCH_TYPES, UNSTRANDED_MISMATCH_TYPES, final_words


def screen(in_dir, out_dir, distances, unique_beds_prefix, stranded):
    """
    Screen editing sites.

    Filter out clusters with multiple mismatches and/or poor enrichment of a single mismatch.

    @param in_dir: dir of one HE run
    @type in_dir: str
    @param out_dir: where to write dirs with different distances screening
    @type out_dir: str
    @param distances: distances to use by the different screens
    @type distances: list[int]
    @param unique_beds_prefix: find unique beds with beds_prefix in in_dir
    @type unique_beds_prefix: str
    @param stranded: data is stranded
    @type stranded: bool
    """
    # prepare arguments
    in_dir = Path(in_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    # distances = [int(distance) for distance in distances.split("_")]
    # filter clusters
    for distance in distances:
        sub_out_dir = Path(out_dir, f"{distance}_BP_Clusters")
        sub_out_dir.mkdir(exist_ok=True)
        one_distance_cluster_screening(in_dir, sub_out_dir, distance, unique_beds_prefix, stranded)
    # copy original in_dir + bed files that start with unique_beds_prefix to out_dir - will be useful later for
    # intersections of the bed files there against certain regions, especially if the signal vs snr in that
    # original in_dir is better than any cluster screening
    in_dir_copied_to_out_dir = Path(out_dir, "Original")
    in_dir_copied_to_out_dir.mkdir(exist_ok=True)
    try:
        stranded_prefix = unique_beds_prefix.split(".")[1] + "."
    except:
        stranded_prefix = ""
    files_to_copy = []
    for file in in_dir.iterdir():
        copy_file = False
        if stranded:
            if file.name.startswith(unique_beds_prefix):
                copy_file = True
            elif file.name.startswith(stranded_prefix):
                copy_file = True
        else:
            if file.name.startswith(unique_beds_prefix) and "stranded" not in file.name:
                copy_file = True
            else:
                for mm_type in UNSTRANDED_MISMATCH_TYPES:
                    if file.name.startswith(mm_type):
                        copy_file = True
                        break
        if copy_file:
            files_to_copy.append(file)
    copied_files_paths = [Path(in_dir_copied_to_out_dir, file.name) for file in files_to_copy]
    for src_file, dest_file in zip(files_to_copy, copied_files_paths):
        copy_text_file(src_file, dest_file)


def one_distance_cluster_screening(in_dir, out_dir, distance=400, unique_beds_prefix="ESuniqS", stranded=False):
    """
    Exclude clusters with different MM, retain only enriched clusters, and write the new files to out_dir.

    First filter by using the unique sites bed files. Then use those filtered unique bed files in order to filter out
    the general bed files.

    @param in_dir: dir with original ES bed files
    @type in_dir: Path
    @param out_dir: where the clean & enriched bed files will be written to
    @type out_dir: Path
    @param distance: maximal cluster size
    @type distance: int
    @param unique_beds_prefix: prefix of wanted bed files
    @type unique_beds_prefix: str
    @param stranded: data is stranded
    @type stranded: bool
    """
    # define in & out paths
    bed_files_names = [file.name for file in in_dir.iterdir() if file.name.startswith(unique_beds_prefix)]
    bed_in_paths = [Path(in_dir, name) for name in bed_files_names]
    bed_out_paths = [Path(out_dir, name) for name in bed_files_names]

    # remove clusters with more than one type of MM
    exclude_multiple_substitutions_clusters(bed_in_paths, bed_out_paths, distance)

    # retain only clusters with at least two nearby MM of the same type
    exclude_poor_enriched_clusters(bed_out_paths, distance)

    # filter the general bed files according to the previous two filtering done on the unique files
    # filter_general_by_unique(in_dir, out_dir, unique_beds_prefix, stranded)
    filter_general_by_unique_one_dir(in_dir, out_dir, unique_beds_prefix, stranded)


def exclude_multiple_substitutions_clusters(bed_in_paths, bed_out_paths, distance):
    """
    Run exclude_multiple_substitutions_clusters_one_file on each unique bed file in bed_files.

    For each bed file, run all other files in bed_files as other_bed_files.
    Each bed file is tested against the original bed files, thus guaranteeing that no cluster is missing because it
    was removed in any previous exclusion (e.g. for files A, B and C, if cluster X was removed from file B when file B
    was compared to file C, file A should be also compared to the original B file such that X will be excluded from A
    too).

    @param bed_in_paths: original bed files
    @type bed_in_paths: list
    @param bed_out_paths: output paths for the bed files after the filtering
    @type bed_out_paths: list
    @param distance: minimal distance between two editing sites of different type (e.g. A2G and A2T) in order to
    consider them on separate clusters
    @type distance: int
    """
    starmap_inputs = []
    for x in range(len(bed_in_paths)):   # could be bed_out_paths as well
        bed_files = [BedTool(bed_in_path).sort() for bed_in_path in bed_in_paths]
        bed_file = bed_files[x]
        bed_file_out_path = bed_out_paths[x]
        other_bed_files = bed_files[:x] + bed_files[x+1:]
        starmap_inputs.append((bed_file, bed_file_out_path, other_bed_files, distance))
    with Pool(processes=20) as pool:
        pool.starmap(func=exclude_multiple_substitutions_clusters_one_file, iterable=starmap_inputs)


def exclude_multiple_substitutions_clusters_one_file(bed_file, bed_file_out_path, other_bed_files, distance):
    """
    Retain clusters of $distance bp in bed_file only if no other mutations exist in any of the other_bed_files.

    The window is used without respect to strand. This makes sense even if the data is stranded, as a mutation in one
    strand probably implies a mutation in the other strand.
    The header of the original bed_file is kept.
    The final bed file is saved according to bed_file_out_path.

    @param bed_file: bed file with certain mutation (e.g. A2G)
    @type bed_file: BedTool
    @param bed_file_out_path: path of final bed file after it's procession
    @type bed_file_out_path: Path
    @param other_bed_files: other bed files with other mutations (e.g. A2T, G2A)
    @type other_bed_files: list
    @param distance: minimal distance between two editing sites of different type (e.g. A2G and A2T) in order to
    consider them on separate clusters
    @type distance: int
    @return: bed_file
    @rtype: BedTool
    """
    for other_bed_file in other_bed_files:
        bed_file = bed_file.window(other_bed_file, header=True, v=True, w=distance)
    bed_file.saveas(bed_file_out_path)


def exclude_poor_enriched_clusters(bed_files, distance):
    """
    Run exclude_poor_enriched_clusters_one_file on each unique bed file in bed_files.

    @param bed_files: bed files
    @type bed_files: list
    @param distance: maximal distance between two editing sites in order to consider them on the same cluster
    @type distance: int
    """
    starmap_inputs = [(bed_file, distance) for bed_file in bed_files]
    with Pool(processes=20) as pool:
        pool.starmap(func=exclude_poor_enriched_clusters_one_file, iterable=starmap_inputs)


def exclude_poor_enriched_clusters_one_file(bed_file, distance):
    """
    Retain only editing sites that are enriched with other nearby (up to $distance bp) editing sites.

    The new file is saved according to its previous.

    @param bed_file: bed file, after multiple substitutions filtering
    @type bed_file: Path
    @param distance: maximal distance between two editing sites in order to consider them on the same cluster
    @type distance: int
    """
    # cluster editing sites with the same type that are within $distance bp of each other
    clustered_bed = BedTool(bed_file).sort().cluster(s=True, d=distance)
    # count how many times each cluster is shown - a cluster that is shown only once is a cluster with only one feature
    # (i.e. there's only one editing site in that cluster) and thus it will be excluded
    clusters = Counter()
    for feature in clustered_bed:
        clusters[feature[6]] += 1
    enriched_clusters = {cluster for cluster in clusters if 1 < clusters[cluster]}

    def bed_filter_func(bed_feature):
        return bed_feature[6] in enriched_clusters
    enriched_bed = clustered_bed.filter(bed_filter_func)
    # remove the cluster tag that was created by the clustering of bedtools
    bed_lines_without_cluster_tag = []
    for feature in enriched_bed:
        new_line = "\t".join([feature[x] for x in range(6)]) + "\n"
        bed_lines_without_cluster_tag.append(new_line)
    bed_lines_joined_to_str = "".join(line for line in bed_lines_without_cluster_tag)
    # save the final file, in place
    BedTool(bed_lines_joined_to_str, from_string=True).sort().saveas(bed_file)


def filter_general_by_unique_one_dir(original_beds_dir, cluster_dir, unique_beds_prefix="ESuniqS",
                                     stranded=False):
    """
    Use the unique filtered beds in cluster_dir in order to filter the general bed files in original_beds_dir.

    @param original_beds_dir: dir with original bed files
    @type original_beds_dir: Path
    @param cluster_dir: one dir of screened cluster
    @type cluster_dir: Path
    @param unique_beds_prefix: find unique beds according to unique_beds_prefix
    @type unique_beds_prefix: str
    @param stranded: data is stranded
    @type stranded: bool
    """
    mismatch_types = STRANDED_MISMATCH_TYPES if stranded else UNSTRANDED_MISMATCH_TYPES
    try:
        stranded_prefix = unique_beds_prefix.split(".")[1] + "."
    except:
        stranded_prefix = ""
    for mismatch_type in mismatch_types:
        general_unfiltered_bed = BedTool(Path(original_beds_dir, f"{stranded_prefix}{mismatch_type}.bed")).sort()
        unique_filtered_bed = BedTool(Path(cluster_dir, f"{unique_beds_prefix}.{mismatch_type}.bed")).sort()
        general_unfiltered_bed.\
            intersect(unique_filtered_bed, u=True, s=True).\
            saveas(Path(cluster_dir, f"{stranded_prefix}{mismatch_type}.bed"))


if __name__ == "__main__":
    # parse args
    # clustering args
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_dir",
                        required=True)
    parser.add_argument("--out_dir",
                        required=True)
    parser.add_argument("--distances",
                        default=[50, 100, 150, 200, 250, 300, 350, 400, 450, 500],
                        nargs="+",
                        help="One or more distances that define clusters of UE regions")
    # analysis args
    parser.add_argument("--genome_file")
    parser.add_argument("--repeats_file",
                        help="BAM/BED/GFF/VCF file")
    parser.add_argument("--cds_file",
                        help="BAM/BED/GFF/VCF file")
    parser.add_argument("--surrounding_bases",
                        type=int, default=1)
    parser.add_argument("--motif_col_position",
                        choices=(0, 1), default=0)
    parser.add_argument("--keep_temp_files",
                        action="store_true")
    parser.add_argument("--skip_header",
                        action="store_true")
    parser.add_argument("--fractions",
                        default=(0.7, 0.8, 0.9, 1.0), nargs="+",
                        help="Fractions of overlap between UE regions and CDS regions")
    # common args
    parser.add_argument("--beds_prefix",
                        default="ESuniqS",
                        help="prefix of bed files to be screened")
    parser.add_argument("--stranded",
                        action="store_true")
    parser.add_argument("--editing_type",
                        default="HE", choices=("HE", "regular"))
    parser.add_argument("--data_table",
                        help="Something like the SRA run table, required columns are ... ")  # TODO complete description
    parser.add_argument("--analyse_a2g_table",
                        help="HE analyse.HE_RUN.A2G table, before any other filtering such as cluster screening")
    parser.add_argument("--notebook_template_path",
                        default="Notebooks/ClusterScreening.ipynb",
                        help="Provide a different notebook template for the Cluster Screening report")
    parser.add_argument("--processes",
                        default=30,
                        type=int,
                        help="number of processes to run in parallel")
    # parse args
    args = parser.parse_args()
    # run
    screen(in_dir=args.in_dir, out_dir=args.out_dir, distances=args.distances, unique_beds_prefix=args.beds_prefix,
           stranded=args.stranded)
    # summarize
    summarize_cluster_screening(original_files_dir=args.in_dir, main_clusters_dir=args.out_dir,
                                unique_beds_prefix=args.beds_prefix, genome_file=args.genome_file,
                                repeats_file=args.repeats_file, cds_file=args.cds_file,
                                surrounding_bases=args.surrounding_bases,
                                motif_col_position=args.motif_col_position, keep_temp_files=args.keep_temp_files,
                                skip_header=args.skip_header, fractions=args.fractions,
                                editing_type=args.editing_type, distances=args.distances,
                                data_table=args.data_table, stranded=args.stranded,
                                analyse_a2g_table=args.analyse_a2g_table, processes=args.processes)
    # end
    final_words()
