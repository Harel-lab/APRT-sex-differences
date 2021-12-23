from collections import namedtuple
from multiprocessing import Pool
from pathlib import Path
import subprocess
import inspect
import sys

from pybedtools import BedTool
import papermill as pm

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from EditingUtils.pybedtools_utils import count_unique_sites, count_total_events
from EditingUtils.summarize_utils import intersect_editing_results, plot_intersections, summarize_intersections, \
    he_total_mmb
from EditingUtils.logo import multiple_logos
from HE.utils import ue_from_es
from General.consts import STRANDED_MISMATCH_TYPES, UNSTRANDED_MISMATCH_TYPES


OUT_NOTEBOOK_NAME = "ClusterScreeningSummary"


def signal_noise_ratio(signal_level, noise_level, ndigits=2):
    """Calc signal-to-noise ratio, rounded to ndigits precision after the decimal point."""
    snr = (signal_level / (signal_level + noise_level)) * 100
    snr = round(snr, ndigits)
    return snr


def find_second_frequent_mismatch(cluster_dir, unique_beds_prefix="ESuniqS"):
    """
    Define which type of mismatch is the 'noise', i.e. the second most common in some cluster_dir.

    Mainly consider total editing events. Break ties according to unique sites.

    @param cluster_dir: dir of some screening
    @type cluster_dir: Path
    @param unique_beds_prefix: unique bed files prefix (usually "ESuniqS" for unstranded data and "ESuniqS.stranded"
    for stranded)
    @type unique_beds_prefix: str
    @return: noise_type
    @rtype: str
    """
    all_other_beds = [BedTool(file) for file in cluster_dir.iterdir()
                      if file.name.startswith(unique_beds_prefix) and "A2G" not in file.name]
    noise_sites = 0
    noise_events = 0
    noise_type = None
    for other_bed in all_other_beds:
        other_sites = count_unique_sites(other_bed)
        other_events = count_total_events(other_bed)
        replace = False
        if other_events < noise_events:
            continue
        elif noise_events < other_events:
            replace = True
        else:  # other_events == noise_events
            if other_sites < noise_sites:
                replace = True
        if replace:
            noise_sites = other_sites
            noise_events = other_events
            noise_type = other_bed.fn.split(".")[-2]
    return noise_type


AllClustersStats = namedtuple("All_Clusters_Stats", ["events", "sites", "events_snrs", "sites_snrs", "distances"])


def optimize_signal_vs_snr(original_files_dir, main_clusters_dir, unique_beds_prefix="ESuniqS"):
    """
    Calc total events, unique sites signal and SNR of each, after sequential cluster screening.

    @param original_files_dir: folder with original bed files
    @type original_files_dir: Path
    @param main_clusters_dir: folder with bed files after cluster screening
    @type main_clusters_dir: Path
    @param unique_beds_prefix: unique bed files prefix (usually "ESuniqS" for unstranded data and "ESuniqS.stranded"
    for stranded)
    @type unique_beds_prefix: str
    @return: all_clusters_stats - containing events, sites, events_snrs, sites_snrs, and distances
    @rtype: AllClustersStats

    """
    cluster_dirs = [cluster_dir for cluster_dir in main_clusters_dir.iterdir()
                    if cluster_dir.is_dir() and "_BP_Clusters" in cluster_dir.name]
    distances = [int(cluster_dir.name.split("_")[0]) for cluster_dir in cluster_dirs]

    original_a2g_bed = BedTool(Path(original_files_dir, f"{unique_beds_prefix}.A2G.bed"))
    original_a2g_sites = count_unique_sites(original_a2g_bed)
    original_a2g_events = count_total_events(original_a2g_bed)

    original_noise_type = find_second_frequent_mismatch(original_files_dir)
    original_noise_bed = BedTool(Path(original_files_dir, f"{unique_beds_prefix}.{original_noise_type}.bed"))
    original_noise_sites = count_unique_sites(original_noise_bed)
    original_noise_events = count_total_events(original_noise_bed)

    sites = [original_a2g_sites]
    sites_snrs = [signal_noise_ratio(original_a2g_sites, original_noise_sites)]
    events = [original_a2g_events]
    events_snrs = [signal_noise_ratio(original_a2g_events, original_noise_events)]

    # for i, (distance, cluster_dir) in enumerate(zip(distances, cluster_dirs)):
    for distance, cluster_dir in zip(distances, cluster_dirs):

        noise_type = find_second_frequent_mismatch(cluster_dir)

        screened_a2g_bed = BedTool(Path(cluster_dir, f"{unique_beds_prefix}.A2G.bed"))
        screened_noise_bed = BedTool(Path(cluster_dir, f"{unique_beds_prefix}.{noise_type}.bed"))
        screened_a2g_sites = count_unique_sites(screened_a2g_bed)
        screened_a2g_events = count_total_events(screened_a2g_bed)
        screened_noise_sites = count_unique_sites(screened_noise_bed)
        screened_noise_events = count_total_events(screened_noise_bed)

        sites.append(screened_a2g_sites)
        events.append(screened_a2g_events)

        sites_snr = signal_noise_ratio(screened_a2g_sites, screened_noise_sites)
        events_snr = signal_noise_ratio(screened_a2g_events, screened_noise_events)
        sites_snrs.append(sites_snr)
        events_snrs.append(events_snr)

    distances = ["Original"] + distances
    all_clusters_stats = AllClustersStats(events, sites, events_snrs, sites_snrs, distances)
    return all_clusters_stats


def find_best_cluster(signals, snrs, distances):
    """
    Find best cluster considering signal vs. signal-to-noise ratio (SNR).

    @param signals: signals of clusters
    @type signals: list
    @param snrs: snrs of clusters
    @type snrs: list
    param distances: the distances used to filter the clusters
    @type distances: list
    @return: chosen_signal - signal of chosen cluster
    @rtype: float
    @return: chosen_snr - SNR of chosen cluster
    @rtype: float
    @return chosen_distance - distance used to create best cluster
    @rtype: int
    """
    chosen_signal = signals[0]
    chosen_snr = snrs[0]
    chosen_distance = distances[0]
    for signal, snr, distance in zip(signals[1:], snrs[1:], distances[1:]):
        if chosen_snr < snr:
            chosen_signal = signal
            chosen_snr = snr
            chosen_distance = distance
    return chosen_signal, chosen_snr, chosen_distance


def summarize_cluster_screening(*, original_files_dir, main_clusters_dir, unique_beds_prefix, genome_file, repeats_file,
                                cds_file, surrounding_bases, motif_col_position, keep_temp_files, skip_header,
                                fractions, distances, editing_type, data_table, analyse_a2g_table, stranded,
                                notebook_template_path="Notebooks/ClusterScreening.ipynb",
                                processes=1):
    """
    Summarize the cluster screening which is done by cluster_screening.py.

    The workflow goes like this:
    1) This function parameterizes and executes a jupyter notebook template with papermill.
    2) The notebook uses other functions defined in this module.
    3) This function converts the executed parameterized notebook.ipynb to an html file.

    @param original_files_dir: dir with original bed files
    @type original_files_dir: str
    @param main_clusters_dir: dir with sub dirs of clustered bed files
    @type main_clusters_dir: str
    @param unique_beds_prefix: find unique_beds with unique_beds_prefix in in_dir
    @type unique_beds_prefix: str
    @param genome_file: genome fasta file path
    @type genome_file: str
    @param repeats_file: repeats file path
    @type repeats_file: str
    @param cds_file: cds file path
    @type cds_file: str
    @param surrounding_bases: number of bases from each side of edited positions (e.g. 1 for ADAR, 3 for APOBEC)
    @type surrounding_bases: int
    @param motif_col_position: used for editing files of old editing_detection, should be removed
    @type motif_col_position: int
    @param keep_temp_files: a flag indicating whether to keep temp files
    @type keep_temp_files: bool
    param skip_header: a flag indicating the first line of bed_file if a comment (e.g. #chrom start end)
    @type skip_header: bool
    @param distances: distances used to define the different screens in main_clusters_dir, e.g. "50_100_150"
    @type distances: list[int]
    @param editing_type: editing type, either "HE" or "regular"
    @type editing_type: str
    @param data_table:
    @type data_table:
    @param analyse_a2g_table:
    @type analyse_a2g_table:
    @param processes: number of processes to run in parallel
    @type processes: int
    """

    # (1) summarize each cluster
    summarize_clusters(main_clusters_dir, distances, genome_file, repeats_file, cds_file, fractions,
                       unique_out_prefix=unique_beds_prefix, surrounding_bases=surrounding_bases,
                       motif_col_position=motif_col_position, keep_temp_files=keep_temp_files, skip_header=skip_header,
                       editing_type=editing_type, data_table=data_table, analyse_a2g_table=analyse_a2g_table,
                       processes=processes, stranded=stranded)

    # (2) choose best cluster and report it in a notebook
    notebook_output_path = f"{main_clusters_dir}/{OUT_NOTEBOOK_NAME}.ipynb"
    # (2.a) execute template notebook with papermill
    pm.execute_notebook(notebook_template_path,
                        notebook_output_path,
                        parameters=dict(original_files_dir=original_files_dir,
                                        main_clusters_dir=main_clusters_dir,
                                        unique_beds_prefix=unique_beds_prefix,
                                        genome_file=genome_file,
                                        repeats_file=repeats_file,
                                        cds_file=cds_file,
                                        stranded=stranded))
    # (2.b) convert executed notebook to html with nbconvert
    convert_cmd = f"jupyter nbconvert --to html {notebook_output_path}"
    subprocess.run(convert_cmd, shell=True)


def summarize_clusters(main_clusters_dir, distances, genome_file, repeats_file, cds_file, fractions,
                       unique_out_prefix="ESuniqS", surrounding_bases=1, motif_col_position=0, keep_temp_files=False,
                       skip_header=False, editing_type="HE", data_table="", analyse_a2g_table="", processes=1,
                       stranded=False):
    """
    Summarize each cluster in main_clusters_dir by running summarize_one_cluster on it.

    @param main_clusters_dir: dir with sub dirs of screened clusters
    @type main_clusters_dir: str
    @param distances: distances used to define the different screens in main_clusters_dir, e.g. "50_100_150"
    @type distances: list[int]
    @param genome_file: genome fasta file
    @type genome_file: str
    @param repeats_file: repeats file
    @type repeats_file: str
    @param cds_file: cds file
    @type cds_file: str
    @param unique_out_prefix: out prefix of unique bed files (in summarize_one_best_cluster.intersect_editing_results)
    @type unique_out_prefix: str
    @param surrounding_bases: number of bases from each side of edited positions (e.g. 1 for ADAR, 3 for APOBEC)
    @type surrounding_bases: int
    @param motif_col_position: used for editing files of old editing_detection, should be removed
    @type motif_col_position: int
    @param keep_temp_files: a flag indicating whether to keep temp files
    @type keep_temp_files: bool
    param skip_header: a flag indicating the first line of bed_file if a comment (e.g. #chrom start end)
    @type skip_header: bool
    @param editing_type: editing type, either "HE" or "regular"
    @type editing_type: str
    @param data_table:
    @type data_table:
    @param analyse_a2g_table:
    @type analyse_a2g_table:
    @param processes: number of processes to run in parallel
    @type processes: int
    """
    main_clusters_dir = Path(main_clusters_dir)
    genome_file = Path(genome_file)
    repeats_file = Path(repeats_file)
    cds_file = Path(cds_file)
    # distances = [int(distance) for distance in distances.split("_")]
    clusters_dirs = [Path(main_clusters_dir, f"{distance}_BP_Clusters") for distance in distances]
    original_dir = Path(main_clusters_dir, "Original")

    try:
        stranded_prefix = unique_out_prefix.split(".")[1] + "."
    except:
        stranded_prefix = ""

    summarize_es_of_one_cluster_inputs = [(cluster_dir, cluster_dir, genome_file, repeats_file, cds_file,
                                           unique_out_prefix, stranded_prefix, surrounding_bases, motif_col_position,
                                           keep_temp_files, skip_header, editing_type, data_table, analyse_a2g_table)
                                          for cluster_dir in clusters_dirs]
    summarize_es_of_one_cluster_inputs.append((original_dir, original_dir, genome_file, repeats_file, cds_file,
                                               unique_out_prefix, stranded_prefix, surrounding_bases,
                                               motif_col_position, keep_temp_files, skip_header, editing_type,
                                               data_table, analyse_a2g_table))
    with Pool(processes=processes) as pool:
        pool.starmap(func=summarize_es_of_one_cluster, iterable=summarize_es_of_one_cluster_inputs)

    mismatches = STRANDED_MISMATCH_TYPES if stranded else UNSTRANDED_MISMATCH_TYPES
    regions_dir_basename = "UERegions"
    sites_within_regions_dir_basename = "ESInCDSRegions"

    for cluster_dir in clusters_dirs + [original_dir]:
        summarize_ue_of_one_cluster(cluster_dir, Path(cluster_dir, regions_dir_basename), stranded_prefix, cds_file,
                                    repeats_file, processes, mismatches, fractions)

    summarize_cds_es_in_ue_of_one_cluster_inputs = [(cluster_dir, Path(cluster_dir, regions_dir_basename),
                                                     Path(cluster_dir, sites_within_regions_dir_basename),
                                                     unique_out_prefix, stranded_prefix, mismatches, fractions)
                                                    for cluster_dir in clusters_dirs + [original_dir]]
    with Pool(processes=processes) as pool:
        pool.starmap(func=summarize_es_in_cds_ue_of_one_cluster,
                     iterable=summarize_cds_es_in_ue_of_one_cluster_inputs)

    plotting_inputs = [(Path(cluster_dir, sites_within_regions_dir_basename), unique_out_prefix, fractions, genome_file,
                       surrounding_bases, motif_col_position, keep_temp_files, skip_header)
                      for cluster_dir in clusters_dirs + [original_dir]]
    with Pool(processes=processes) as pool:
        pool.starmap(func=plot_motifs_of_es_in_cds_ue_of_one_cluster, iterable=plotting_inputs)


def summarize_es_in_cds_ue_of_one_cluster(cluster_dir, regions_dir, sites_within_regions_dir, unique_out_prefix,
                                          stranded_prefix, mismatches, fractions):
    sites_within_regions_dir.mkdir(exist_ok=True)
    for mismatch in mismatches:
        all_editing_sites = BedTool(Path(cluster_dir, f"{stranded_prefix}{mismatch}.bed"))
        unique_editing_sites = BedTool(Path(cluster_dir, f"{unique_out_prefix}.{mismatch}.bed"))
        for fraction in fractions:
            # "pure" CDS UE regions
            cds_ue_regions = BedTool(Path(regions_dir, f"UEmergeSd20.{stranded_prefix}{mismatch}.CDS-f{fraction}.bed"))
            all_editing_sites.\
                intersect(cds_ue_regions, u=True, s=True).\
                sort().\
                saveas(Path(sites_within_regions_dir, f"{stranded_prefix}{mismatch}.CDS-f{fraction}.bed"))
            unique_editing_sites.\
                intersect(cds_ue_regions, u=True, s=True).\
                sort().\
                saveas(Path(sites_within_regions_dir, f"{unique_out_prefix}.{mismatch}.CDS-f{fraction}.bed"))
            # "pure" CDS UE regions, that are also not annotated by repeats
            cds_wo_repeats_ue_regions =\
                BedTool(Path(regions_dir, f"UEmergeSd20.{stranded_prefix}{mismatch}.CDS-f{fraction}.WoRepeats.bed"))
            all_editing_sites.\
                intersect(cds_wo_repeats_ue_regions, u=True, s=True).\
                sort().\
                saveas(Path(sites_within_regions_dir, f"{stranded_prefix}{mismatch}.CDS-f{fraction}.WoRepeats.bed"))
            unique_editing_sites. \
                intersect(cds_wo_repeats_ue_regions, u=True, s=True).\
                sort().\
                saveas(Path(sites_within_regions_dir, f"{unique_out_prefix}.{mismatch}.CDS-f{fraction}.WoRepeats.bed"))


def plot_motifs_of_es_in_cds_ue_of_one_cluster(sites_within_regions_dir, unique_out_prefix, fractions, genome_file,
                                               surrounding_bases, motif_col_position, keep_temp_files, skip_header,
                                               plots_file_type="svg"):
    # plot motif of A2G sites in "pure" CDS UE regions
    bed_files = [Path(sites_within_regions_dir, f"{unique_out_prefix}.A2G.CDS-f{fraction}.bed")
                 for fraction in fractions]
    titles = [f"-f {fraction}" for fraction in fractions]
    main_title = "ADAR motifs in UE regions with different overlap fractions to CDS regions"
    for plot_postfix in ["png", plots_file_type]:
        output_file = Path(sites_within_regions_dir, f"ADAR-motifs.A2G.CDS.{plot_postfix}")
        multiple_logos(bed_files, genome_file, titles=titles, output_file=output_file,
                       surrounding_bases=surrounding_bases, motif_col_position=motif_col_position,
                       keep_temp_files=keep_temp_files, skip_header=skip_header, main_title=main_title)
    # plot motif of A2G sites in "pure" CDS UE regions, that are also not annotated by repeats
    bed_files = [Path(sites_within_regions_dir, f"{unique_out_prefix}.A2G.CDS-f{fraction}.WoRepeats.bed")
                 for fraction in fractions]
    # titles = [f"-f {fraction} CDS UE regions, repeats excluded" for fraction in fractions]
    titles = [f"-f {fraction}" for fraction in fractions]
    main_title = "ADAR motifs in UE regions with different overlap fractions to CDS regions, repeats excluded"
    for plot_postfix in ["png", plots_file_type]:
        output_file = Path(sites_within_regions_dir, f"ADAR-motifs.A2G.CDS.WoRepeats.{plot_postfix}")
        multiple_logos(bed_files, genome_file, titles=titles, output_file=output_file,
                       surrounding_bases=surrounding_bases, motif_col_position=motif_col_position,
                       keep_temp_files=keep_temp_files, skip_header=skip_header, main_title=main_title)



def summarize_ue_of_one_cluster(cluster_dir, regions_dir, stranded_prefix, cds_file, repeats_file, processes,
                                mismatches, fractions):
    cds_bed = BedTool(cds_file).sort()
    repeats_bed = BedTool(repeats_file).sort()
    regions_dir.mkdir(exist_ok=True)
    for mismatch in mismatches:
        general_editing_file = Path(cluster_dir, f"{stranded_prefix}{mismatch}.bed")
        merged_regions_bed = ue_from_es(in_editing_sites=general_editing_file, processes=processes)
        for fraction in fractions:  # e.g. 0.7, 0.8, 0.9, 1.0
            merged_regions_bed.\
                intersect(cds_bed, sorted=True, s=True, u=True, f=float(fraction)).\
                sort().\
                saveas(Path(regions_dir, f"UEmergeSd20.{stranded_prefix}{mismatch}.CDS-f{fraction}.bed")).\
                intersect(repeats_bed, sorted=True, s=True, v=True).\
                sort().\
                saveas(Path(regions_dir, f"UEmergeSd20.{stranded_prefix}{mismatch}.CDS-f{fraction}.WoRepeats.bed"))
        merged_regions_bed.saveas(Path(regions_dir, f"UEmergeSd20.{stranded_prefix}{mismatch}.bed"))


def summarize_es_of_one_cluster(cluster_dir, out_dir, genome_file, repeats_file, cds_file, unique_out_prefix="ESuniqS",
                                stranded_prefix="", surrounding_bases=1, motif_col_position=0, keep_temp_files=False,
                                skip_header=False, editing_type="HE", data_table="", analyse_a2g_table="",
                                plots_file_type="svg"):
    """
    Summarize one cluster.

    @param cluster_dir: one cluster created by cluster_screening.py
    @type cluster_dir: Path
    @param out_dir: where to write output files to
    @type out_dir: Path
    @param genome_file: genome fasta file
    @type genome_file: Path
    @param repeats_file: repeats file
    @type repeats_file: Path
    @param cds_file: cds file
    @type cds_file: Path
    @param unique_out_prefix: out prefix of unique bed files (in summarize_one_best_cluster.intersect_editing_results)
    @type unique_out_prefix: str
    @param surrounding_bases: number of bases from each side of edited positions (e.g. 1 for ADAR, 3 for APOBEC)
    @type surrounding_bases: int
    @param motif_col_position: used for editing files of old editing_detection, should be removed
    @type motif_col_position: int
    @param keep_temp_files: a flag indicating whether to keep temp files
    @type keep_temp_files: bool
    param skip_header: a flag indicating the first line of bed_file if a comment (e.g. #chrom start end)
    @type skip_header: bool
    @param editing_type: editing type, either "HE" or "regular"
    @type editing_type: str
    @param data_table:
    @type data_table:
    @param analyse_a2g_table:
    @type analyse_a2g_table:
    """
    # (1) find editing in different regions
    unique_a2g_file = BedTool(Path(cluster_dir, f"{unique_out_prefix}.A2G.bed")).sort()
    general_a2g_file = BedTool(Path(cluster_dir, f"{stranded_prefix}A2G.bed")).sort()
    intersect_editing_results(unique_a2g_file, general_a2g_file, repeats_file, cds_file, out_dir, unique_out_prefix,
                              stranded_prefix)

    # (2) plot motifs of editing in the different regions
    plot_intersections(intersections_dir=out_dir,
                       out_dir=out_dir,
                       genome_file=genome_file,
                       unique_beds_prefix=unique_out_prefix,
                       surrounding_bases=surrounding_bases,
                       motif_col_position=motif_col_position,
                       keep_temp_files=keep_temp_files,
                       skip_header=skip_header,
                       plots_file_type=plots_file_type)

    # (3) summarize base statistics for the editing in the different regions
    if editing_type == "HE":
        total_mmb = he_total_mmb(data_table, analyse_a2g_table)
    else:
        raise NotImplementedError("Summarizing cluster is not implemented for regular editing - a way to provide"
                                  "number of mapped reads should be developed")
    summarize_intersections(intersections_dir=out_dir,
                            out_dir=out_dir,
                            total_mmb=total_mmb,
                            unique_beds_prefix=unique_out_prefix)
