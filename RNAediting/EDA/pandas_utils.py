def reorder_df_by_wanted_cols(df, wanted_first_cols):
    """
    Reorder the columns in df such that columns in wanted_first_cols will be the first ones.

    The relative order of the other columns will remain unchanged.

    @param df: data frame
    @type df: pandas.DataFrame
    @param wanted_first_cols: columns to bring upfront (to the left)
    @type wanted_first_cols: list
    @return: df
    @rtype: pandas.DataFrame
    """
    wanted_first_cols_set = set(wanted_first_cols)
    all_cols_set = set(df.columns)
    other_cols_set = all_cols_set - wanted_first_cols_set
    other_cols = [col for col in df.columns if col in other_cols_set]
    new_cols = wanted_first_cols + other_cols
    df = df.reindex(new_cols, axis=1)
    return df


def sort_by_group_mean_val(df, group_col, signal_col, ascending=True):
    """
    Sort df by the mean value of each group.

    For example, for a given df such as

        group 	val
    0 	A 	    1
    1 	B 	    2
    2 	B    	2
    3 	A    	6

    the function returns

         group  val
    0    B      2
    1    B      2
    2    A      1
    3    A      6

    @param df: data frame
    @type df: pandas.DataFrame
    @param group_col: label of the column the holds the group of each row/instance/sample
    @type group_col: str
    @param signal_col: label of the column the average upon
    @type signal_col: str
    @param ascending: sort by ascending order of mean group value
    @type ascending: bool
    @return: df
    @rtype: pandas.DataFrame
    """
    # mean_df = df.groupby(group_col).mean(skipna=False).fillna(0).reset_index()
    mean_df = df.groupby(group_col).mean().reset_index()
    mean_col = \
        df.apply(lambda row: mean_df.loc[mean_df[group_col].values == row[group_col], signal_col].values[0], axis=1)
    df = df.assign(mean_col=mean_col)
    df = df.sort_values("mean_col", ascending=ascending).reset_index(drop=True)
    df = df.drop("mean_col", axis=1)
    return df


def sort_groups_by_median(df, group_col, signal_col, ascending=True):
    median_df = df.groupby(group_col).median().reset_index()
    median_col = \
        df.apply(lambda row: median_df.loc[median_df[group_col].values == row[group_col], signal_col].values[0], axis=1)
    df = df.assign(median_col=median_col)
    df = df.sort_values("median_col", ascending=ascending).reset_index(drop=True)
    df = df.drop("median_col", axis=1)
    return df


def sort_groups_by_mean(df, group_col, signal_col, ascending=True):
    """
    An alias for sort_by_group_mean_val
    """
    return sort_by_group_mean_val(df, group_col, signal_col, ascending=ascending)


def define_groups_order(df, group_col, signal_col, ascending=True, strategy="median"):
    if strategy == "median":
        df = sort_groups_by_median(df, group_col, signal_col, ascending=ascending)
    elif strategy == "mean":
        df = sort_groups_by_mean(df, group_col, signal_col, ascending=ascending)
    else:
        raise KeyError("choose either 'median' or 'mean' as a strategy to sort groups in df")
    groups_order = []
    for group in df[group_col]:
        if group not in groups_order:
            groups_order.append(group)
        else:
            # check if groups_order = [..., A, B, ...] but group = A
            if groups_order.index(group) != len(groups_order) - 1:
                raise Exception("df isn't sorted by discrete homogeneous groups")
    if len(groups_order) != len(df[group_col].unique()):
        raise Exception("len(groups_order) != len(df[group_col])")
    return groups_order
