import pandas as pd
from pybedtools import BedTool, featurefuncs


def count_unique_sites(bedtool):
    """
    Count unique editing sites.

    @param bedtool: general RNA editing bed file
    @type bedtool: BedTool
    @return: unique_editing_sites
    @rtype: int
    """
    unique_editing_sites = len(bedtool)
    return unique_editing_sites


def count_total_events(bedtool):
    """
    Count total editing events.

    @param bedtool: unique RNA editing bed file
    @type bedtool: BedTool
    @return: total_editing_events
    @rtype: int
    """
    editing_events = [int(interval[4]) for interval in bedtool]  # 4 is the 5th column
    total_editing_events = sum(editing_events)
    return total_editing_events


def extract_column(bedtool, col_num, casting=None):
    """
    Extract column at index col_num from bedtool.

    @param bedtool: bed file
    @type bedtool: BedTool
    @param col_num: index of wanted column (0-based)
    @type col_num: int
    @param casting: cast features in column to casting if casting != None
    @type casting: TypeVar
    @return: column
    @rtype: list
    """
    column = [line[col_num] for line in bedtool]
    if casting:
        column = [casting(x) for x in column]
    return column


def intersect_dfs(df1, df2, intersect_kwargs=None, read_table_names=None, other_read_table_kwargs=None):
    """
    Intersect two pandas DataFrames using pybedtools.

    @param df1: file A
    @type df1: pandas.DataFrame
    @param df2 file B
    @type df2: pandas.DataFrame
    @param intersect_kwargs: kwargs passed to pybedtools.intersect
    @type intersect_kwargs: dict
    @param read_table_names: list of column names passed to pandas.read_table instead of the default ones given by
    pybedtools
    @type read_table_names: list[str]
    @param other_read_table_kwargs: kwargs passed to pandas.read_table other than `header` and `names`
    @type other_read_table_kwargs: dict
    @return: intersected_df
    @rtype: pandas.DataFrame
    """
    bed1 = BedTool.from_dataframe(df1)
    bed2 = BedTool.from_dataframe(df2)

    intersect_kwargs = {} if intersect_kwargs is None else intersect_kwargs

    intersected_bed = bed1.intersect(bed2, **intersect_kwargs)

    read_table_kwargs = {}
    if read_table_names is not None:
        read_table_kwargs["header"] = None
        read_table_kwargs["names"] = read_table_names
    if other_read_table_kwargs is not None:
        read_table_kwargs.update(other_read_table_kwargs)

    intersected_df = intersected_bed.to_dataframe(**read_table_kwargs)

    return intersected_df


def no_opposing_no_overlapping(bed):
    """
    Filter out either opposing features (on different strands) or overlapping features (on the same strand).

    @param bed: bed file
    @type bed: BedTool
    @return: no_opposing_no_overlapping_features
    @rtype: BedTool
    """
    # all_features = bed
    # # opposing_features - on opposing strands
    # opposing_features = bed.intersect(all_features, S=True, u=True)
    # no_opposing_features = bed.intersect(all_features, S=True, v=True)
    # assert len(all_features) == len(opposing_features) + len(no_opposing_features)
    no_opposing_features = get_non_opposing_features(bed).sort()
    merged_no_opposing_features = no_opposing_features.merge(s=True, c=[4, 5, 6], o=["distinct", "count", "distinct"])
    # overlapping features - on the same strand
    no_opposing_no_overlapping_features = merged_no_opposing_features.filter(lambda f: int(f[4]) <= 1).saveas()
    return no_opposing_no_overlapping_features


def gff_to_no_opposing_no_overlapping_bed6(gff, name_field=None):
    """
    Convert to bed6 file wile filtering out opposing- and overlapping-features (on the opposite/same strand).

    @param gff: gff file
    @type gff: str or BedTool
    @param name_field: see https://daler.github.io/pybedtools/pybedtools.featurefuncs.gff2bed.html#pybedtools.featurefuncs.gff2bed
    @type name_field: str or int
    @return: no_opposing_no_overlapping_features
    @rtype: BedTool
    """
    gff = BedTool(gff).sort()
    # opposing_features - on opposing strands
    no_opposing_features = get_non_opposing_features(gff)
    # convert the Intervals from gff Intervals to bed6 Intervals
    if name_field:
        no_opposing_features = BedTool([featurefuncs.gff2bed(f, name_field=name_field)
                                        for f in no_opposing_features])
    else:
        no_opposing_features = BedTool([featurefuncs.gff2bed(f) for f in no_opposing_features])
    no_opposing_features = no_opposing_features.saveas().sort()
    merged_no_opposing_features = no_opposing_features.merge(s=True, c=[4, 5, 6], o=["distinct", "count", "distinct"])
    # overlapping features - on the same strand
    no_opposing_no_overlapping_features = merged_no_opposing_features.filter(lambda f: int(f[4]) <= 1).saveas()
    return no_opposing_no_overlapping_features


def get_non_opposing_features(bedtool):
    all_features = bedtool
    # opposing_features - on opposing strands
    opposing_features = bedtool.intersect(all_features, S=True, u=True)
    no_opposing_features = bedtool.intersect(all_features, S=True, v=True)
    assert len(all_features) == len(opposing_features) + len(no_opposing_features)
    return no_opposing_features
