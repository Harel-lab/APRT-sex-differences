from multiprocessing import Pool

from pybedtools import BedTool


def ue_from_es(in_editing_sites, out_regions=None, processes=1):
    """
    Process an 'ES.bed_files/$type.bed' into 'UE.bed_files/UEmergeSd20.$type.bed'.

    If out_regions is None, return ue_merged_df. Else write ue_merged_df to out_regions.

    @param in_editing_sites: either a path or a BedTool of an 'ES.bed_files/$type.bed' file
    @type in_editing_sites: Path or str or BedTool
    @param out_regions: write 'UE.bed_files/UEmergeSd20.$type.bed' to out_regions if out_regions != None
    @type out_regions: str or Path
    @param processes: use $processes in parallel
    @type processes: int
    @return: ue_merged_df
    @rtype: BedTool or None
    """
    es_bed = BedTool(in_editing_sites).sort()
    es_df = es_bed.to_dataframe(header=None,
                                names=["chrom", "start", "end", "read", "triplet", "strand"])
    es_df.loc[:, "read"] = es_df.loc[:, "read"].astype(str).str.split(";").str[0]
    reads = es_df["read"].unique()
    get_cluster_tuple_inputs = [(es_df, read) for read in reads]
    with Pool(processes=processes) as pool:
        clusters_tuples = pool.starmap(func=get_cluster_tuple, iterable=get_cluster_tuple_inputs)
    clusters_bed = BedTool(clusters_tuples)
    ue_merged_bed = clusters_bed.sort().merge(s=True, d=20, c=6, o="distinct")
    ue_merged_df = ue_merged_bed.to_dataframe()
    ue_merged_df.insert(3, "empty_1", 0)
    ue_merged_df.insert(3, "empty2", 0)
    if out_regions:
        ue_merged_df.to_csv(out_regions, sep="\t", index=False, header=False)
    else:
        ue_merged_bed = BedTool.from_dataframe(ue_merged_df)
        return ue_merged_bed


def get_cluster_tuple(es_df, read):
    read_df = es_df.loc[es_df["read"] == read]
    cluster_chrom = read_df["chrom"].values[0]
    cluster_start = min(read_df["start"].astype(int))
    cluster_end = max(read_df["end"].astype(int))
    cluster_strand = read_df["strand"].values[0]
    cluster_tuple = (cluster_chrom, cluster_start, cluster_end, 0, 0, cluster_strand)
    return cluster_tuple
