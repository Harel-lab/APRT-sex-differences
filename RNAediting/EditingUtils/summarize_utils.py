from itertools import chain
from pathlib import Path
import subprocess
import inspect
import sys

from pybedtools import BedTool
import papermill as pm
import pandas as pd

sys.path.append(f"{Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent}")
from EditingUtils.pybedtools_utils import count_total_events, count_unique_sites
from EditingUtils.logo import multiple_logos


# intersected bed files base names
ALL_A2G = "A2G"  # actually not intersected - this is the original file
A2G_IN_REPEATS = "A2G_in_REPEATS"
A2G_IN_ALL_CDS = "A2G_in_all_CDS"
A2G_IN_CDS_W_REPEATS = "A2G_in_CDS_w_REPEATS"
A2G_IN_CDS_WO_REPEATS = "A2G_in_CDS_wo_REPEATS"
ALL_FILES = [ALL_A2G, A2G_IN_REPEATS, A2G_IN_ALL_CDS, A2G_IN_CDS_W_REPEATS, A2G_IN_CDS_WO_REPEATS]


def intersect_editing_results(unique_a2g_file, general_a2g_file, repeats_file, orfs_file, out_dir,
                              unique_out_prefix="ESuniqS", stranded_prefix=""):
    """
    Intersect A2G editing bed file with repeats and orfs files.

    Also create general A2G files for A2G_in_CDS_wo_REPEATS and A2G_in_REPEATS.

    @param unique_a2g_file: unique A2G bed file
    @type unique_a2g_file: Path
    @param general_a2g_file: general A2G bed file
    @type general_a2g_file: Path
    @param repeats_file: a file with repeats annotations
    @type repeats_file: Path
    @param orfs_file: a file with orfs annotations
    @type orfs_file: Path
    @param out_dir: write intersected bed files to out_dir
    @type out_dir: Path
    @param unique_out_prefix: out prefix of unique bed files
    @type unique_out_prefix: str
    @param stranded_prefix: prefix of stranded bed files, either "stranded." or "" if files aren't stranded, default=""
    @type stranded_prefix: str
    """
    # annotation files
    repeats_bed = BedTool(repeats_file).sort()
    orfs_bed = BedTool(orfs_file).sort()

    # all a2g bed files
    unique_a2g_bed = BedTool(unique_a2g_file).sort()
    general_a2g_bed = BedTool(general_a2g_file).sort()

    # editing sites in repetitive regions
    unique_a2g_in_repeats = unique_a2g_bed.\
        intersect(repeats_bed, u=True, s=True).\
        sort().\
        saveas(Path(out_dir, f"{unique_out_prefix}.{A2G_IN_REPEATS}.bed"))
    general_a2g_in_repeats = general_a2g_bed.\
        intersect(repeats_bed, u=True, s=True).\
        sort().\
        saveas(Path(out_dir, f"{stranded_prefix}{A2G_IN_REPEATS}.bed"))

    # editing sites in all cds regions
    unique_a2g_in_all_cds = unique_a2g_bed.\
        intersect(orfs_bed, u=True, s=True).\
        sort().\
        saveas(Path(out_dir, f"{unique_out_prefix}.{A2G_IN_ALL_CDS}.bed"))
    general_a2g_in_all_cds = general_a2g_bed.\
        intersect(orfs_bed, u=True, s=True).\
        sort().\
        saveas(Path(out_dir, f"{stranded_prefix}{A2G_IN_ALL_CDS}.bed"))

    # editing sites in cds regions, which are also annotated as repetitive regions
    unique_a2g_in_cds_w_repeats = unique_a2g_in_all_cds.\
        intersect(unique_a2g_in_repeats, u=True, s=True).\
        sort().\
        saveas(Path(out_dir, f"{unique_out_prefix}.{A2G_IN_CDS_W_REPEATS}.bed"))
    general_a2g_in_cds_w_repeats = general_a2g_in_all_cds.\
        intersect(general_a2g_in_repeats, u=True, s=True).\
        sort().\
        saveas(Path(out_dir, f"{stranded_prefix}{A2G_IN_CDS_W_REPEATS}.bed"))

    # editing sites in bona-fide cds regions, i.e. without the ambiguous editing sites that are found both in cds and
    # repeats
    if 0 < len(unique_a2g_in_cds_w_repeats):
        unique_a2g_in_all_cds.\
            intersect(unique_a2g_in_cds_w_repeats, v=True, s=True).\
            sort().\
            saveas(Path(out_dir, f"{unique_out_prefix}.{A2G_IN_CDS_WO_REPEATS}.bed"))
    else:
        unique_a2g_in_all_cds.saveas(Path(out_dir, f"{unique_out_prefix}.{A2G_IN_CDS_WO_REPEATS}.bed"))
    if 0 < len(general_a2g_in_cds_w_repeats):
        general_a2g_in_all_cds.\
            intersect(general_a2g_in_cds_w_repeats, v=True, s=True).\
            sort().\
            saveas(Path(out_dir, f"{stranded_prefix}{A2G_IN_CDS_WO_REPEATS}.bed"))
    else:
        general_a2g_in_all_cds.saveas(Path(out_dir, f"{stranded_prefix}{A2G_IN_CDS_WO_REPEATS}.bed"))


def plot_intersections(intersections_dir, out_dir, genome_file, unique_beds_prefix="ESuniqS", surrounding_bases=1,
                       motif_col_position=0, keep_temp_files=False, skip_header=False, plots_file_type="svg"):
    """
    Plot each intersection's motif (and the original A2G bed file's motif too).

    @param intersections_dir: dir with some intersected bed files (summarize_editing_utils.intersect_editing_results)
    @type intersections_dir: Path
    @param out_dir: where to write the plots to
    @type out_dir: Path
    @param genome_file: genome fasta file
    @type genome_file: Path
    @param unique_beds_prefix: use only unique bed files with unique_beds_prefix in intersections_dir
    @type unique_beds_prefix: str
    @param surrounding_bases: number of bases from each side of edited positions (e.g. 1 for ADAR, 3 for APOBEC)
    @type surrounding_bases: int
    @param motif_col_position: used for editing files of old editing_detection, should be removed
    @type motif_col_position: int
    @param keep_temp_files: a flag indicating whether to keep temp files
    @type keep_temp_files: bool
    @param skip_header: a flag indicating the first line of bed_file is a header line (e.g. #chrom start end)
    @type skip_header: bool
    """
    bed_files = [Path(intersections_dir, f"{unique_beds_prefix}.{bed_file_name}.bed") for bed_file_name in
                 ("A2G", "A2G_in_REPEATS", "A2G_in_CDS_wo_REPEATS")]
    titles = ["All A2G", "A2G in Repeats", "A2G in CDS (Repeats excluded)"]
    main_title = "ADAR motifs in different regions"
    for plot_postfix in ["png", plots_file_type]:
        output_file = Path(out_dir, f"ADAR_motifs.{plot_postfix}")
        multiple_logos(bed_files, genome_file, titles=titles, output_file=output_file, surrounding_bases=surrounding_bases,
                       motif_col_position=motif_col_position, keep_temp_files=keep_temp_files, skip_header=skip_header,
                       main_title=main_title)


def summarize_intersections(intersections_dir, out_dir, total_mmb, unique_beds_prefix="ESuniqS"):
    """
    Create a pandas table reporting all the numbers, including overlapping percentages, of the different intersections.

    @param intersections_dir: dir with intersected unique bed files
    @type intersections_dir: Path
    @param out_dir: write intersections stats to out_dir
    @type out_dir: Path
    @param unique_beds_prefix: prefix of intersected unique bed files
    @type unique_beds_prefix: str
    @param total_mmb: total millions of mapped bases
    @type total_mmb: float
    """
    # find intersected input files
    all_a2g = Path(intersections_dir, f"{unique_beds_prefix}.{ALL_A2G}.bed")
    a2g_in_repeats = Path(intersections_dir, f"{unique_beds_prefix}.{A2G_IN_REPEATS}.bed")
    a2g_in_all_cds = Path(intersections_dir, f"{unique_beds_prefix}.{A2G_IN_ALL_CDS}.bed")
    a2g_in_cds_w_repeats = Path(intersections_dir, f"{unique_beds_prefix}.{A2G_IN_CDS_W_REPEATS}.bed")
    a2g_in_cds_wo_repeats = Path(intersections_dir, f"{unique_beds_prefix}.{A2G_IN_CDS_WO_REPEATS}.bed")
    in_bed_files = [all_a2g, a2g_in_repeats, a2g_in_all_cds, a2g_in_cds_w_repeats, a2g_in_cds_wo_repeats]
    in_bed_files = [BedTool(bed_file) for bed_file in in_bed_files]

    # summarize stats in a pandas df
    unique_sites = [count_unique_sites(bed_file) for bed_file in in_bed_files]
    total_events = [count_total_events(bed_file) for bed_file in in_bed_files]
    all_a2g_sites = count_unique_sites(in_bed_files[0])
    all_a2g_events = count_total_events(in_bed_files[0])
    unique_sites_shared_with_all_a2g = [sites / all_a2g_sites * 100 for sites in unique_sites]
    total_events_shared_with_all_a2g = [events / all_a2g_events * 100 for events in total_events]

    # avg_events_per_site = [events / sites for events, sites in zip(total_events, unique_sites)]
    avg_events_per_site = []
    for events, sites in zip(total_events, unique_sites):
        if 0 < sites:
            avg_events_per_site.append(events / sites)
        else:
            avg_events_per_site.append(0)


    intersections_df = pd.DataFrame({"uniq_sites": unique_sites,
                                     "uniq_sites_overlap_A2G(%)": unique_sites_shared_with_all_a2g,
                                     "tot_events": total_events,
                                     "tot_events_overlap_A2G(%)": total_events_shared_with_all_a2g,
                                     "events_per_site": avg_events_per_site},
                                    index=ALL_FILES)
    intersections_df["uniq_sites_per_mmb"] = intersections_df["uniq_sites"] / total_mmb
    intersections_df["tot_events_per_mmb"] = intersections_df["tot_events"] / total_mmb

    # write df to out_dir
    intersections_df.to_csv(path_or_buf=f"{out_dir}/intersections.csv")


def he_total_mmb(run_table, analyse_a2g_table, paired_end=False):
    """
    Calc millions of mapped bases in an HE run.

    @param run_table: SRA run table
    @type run_table: Path or str
    @param analyse_a2g_table: HE analyse.HE_RUN.A2G table - before any other filtering such as cluster screening
    @type analyse_a2g_table: Path or str
    @return: stats
    @rtype: Pandas.DataFrame
    """
    run_table_df = pd.read_csv(run_table)
    run_table_df.sort_values(by="Run", inplace=True)
    analyse_a2g_names = [
        "sample", "src_path", "out_pre", "src_reads", "Map_reads", "%Map_reads", "unMap_reads", "%unMap_reads",
        "Edit_type",
        "UE_clusters", "%UE_of_tot", "Avg._clust_len", "Tot_ES", "%ES_of_tot", "Avg._ES", "merge_UE", "%mergeUE_of_tot",
        "Avg._Reg_len", "uniq_ES", "%uniqES_of_tot", "Avg._ES_2", "Norm:HEreads/mapReads(M)",
        "Norm:UniqSites/mapReads(M)",
        "(+)strand", "G-up", "G-down"]
    analyse_a2g_df = pd.read_csv(analyse_a2g_table, sep="\t", names=analyse_a2g_names, skiprows=5)
    total_lines = analyse_a2g_df.loc[analyse_a2g_df["sample"] == "Total"]
    first_total_line_index = total_lines.index[0]
    analyse_a2g_df = analyse_a2g_df.iloc[:first_total_line_index]
    analyse_a2g_df.insert(1, "Run", analyse_a2g_df["sample"].str.split("-").str[0])
    if paired_end:
        reads_lengths = [int(avg_spot_len) / 2 for avg_spot_len in run_table_df["AvgSpotLen"]]
    else:
        reads_lengths = [int(avg_spot_len) for avg_spot_len in run_table_df["AvgSpotLen"]]
    if len(run_table_df) < len(analyse_a2g_df):  # paired-end was run as single-end
        reads_lengths = list(chain.from_iterable((reads_len, reads_len) for reads_len in reads_lengths))
    mmr = [(int(mapped_reads) / 10 ** 6) for mapped_reads in analyse_a2g_df["Map_reads"]]
    mmb = [reads_len * million_map_reads for reads_len, million_map_reads in zip(reads_lengths, mmr)]
    total_mmb = sum(mmb)
    return total_mmb


def strand_bias_in_general_bed(general_bed_path, ndigits=2):
    general_bed_df = pd.read_csv(general_bed_path, sep="\t",
                                 names=["chromosome", "start", "end", "read_data", "triplet", "strand"])
    general_bed_df["plus_orientation_mapped_reads"] = \
        general_bed_df["read_data"].apply(lambda x: True if x[1] == "+" else False)

    try:
        plus_orientation_mapped_reads = len(general_bed_df.loc[general_bed_df["plus_orientation_mapped_reads"] == True])
        minus_orientation_mapped_reads = len(general_bed_df) - plus_orientation_mapped_reads
        plus_orientation_mapped_reads_percentage = plus_orientation_mapped_reads / len(general_bed_df) * 100
        minus_orientation_mapped_reads_percentage = minus_orientation_mapped_reads / len(general_bed_df) * 100
        plus_orientation_mapped_reads_percentage = round(plus_orientation_mapped_reads_percentage, ndigits)
        minus_orientation_mapped_reads_percentage = round(minus_orientation_mapped_reads_percentage, ndigits)
        return plus_orientation_mapped_reads_percentage, minus_orientation_mapped_reads_percentage
    except:
        error_msg = "len(general_bed_df) = " + str(len(general_bed_df)) +\
                    "\ngeneral_bed_path.parent = " + general_bed_path.parent + "   general_bed_path.name = " + \
                    general_bed_path.name
        raise(error_msg)


def execute_notebook(template_notebook, execution_out_dir, executed_notebook_name, **parameters):

    full_executed_notebook_ipynb_path = f"{execution_out_dir}/{executed_notebook_name}.ipynb"
    clean_executed_notebook_html_path = f"{execution_out_dir}/{executed_notebook_name}.no-input.no-prompt.html"

    print(f"{template_notebook = }")
    print(f"{execution_out_dir = }")
    print(f"{executed_notebook_name = }")
    print(f"{parameters = }")
    print(f"{full_executed_notebook_ipynb_path = }")
    print(f"{clean_executed_notebook_html_path = }")
    # return   # todo delete

    pm.execute_notebook(template_notebook,
                        full_executed_notebook_ipynb_path,
                        parameters=parameters)

    full_convert_cmd = f"jupyter nbconvert --to html {full_executed_notebook_ipynb_path} " \
                       f"--TagRemovePreprocessor.remove_input_tags item"
    subprocess.run(full_convert_cmd, shell=True)

    clean_convert_cmd = f"jupyter nbconvert --to html --stdout --no-input --no-prompt " \
                        f"{full_executed_notebook_ipynb_path} > {clean_executed_notebook_html_path}"
    subprocess.run(clean_convert_cmd, shell=True)
