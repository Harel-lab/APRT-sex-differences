import pandas as pd


def read_index_df(index_df_path):
    """
    Read and parse EditingIndex.csv table into a DataFrame.

    Remove unnecessary cols, and rows whose values are cols names (must be a bug in the Index).

    @param index_df_path: path of the EditingIndex.csv file
    @type index_df_path: Path
    @return: index_df
    @rtype: pandas.DataFrame
    """
    index_df = pd.read_csv(index_df_path,
                           usecols=lambda x: x not in ["StrandDecidingMethod", "SamplePath"])

    # remove rows whose values are the cols names
    rows_repeating_col_names = index_df.apply(lambda row: (row.values == index_df.columns).all(), axis=1)
    index_df = index_df.loc[~rows_repeating_col_names]
    # convert signal cols from str to float - they were originally inferred as str if some rows were cols
    # names as above-mentioned
    signal_col_names = index_df.iloc[:, 2:].columns
    index_df = index_df.astype({col: float for col in signal_col_names})

    return index_df


def signal_noise_ratio(signal_level, noise_level, ndigits=2):
    """Calc signal-to-noise ratio, rounded to ndigits precision after the decimal point."""
    snr = (signal_level / (signal_level + noise_level)) * 100
    snr = round(snr, ndigits)
    return snr


def find_highest_noise(index_df, signal_col="A2GEditingIndex"):
    # index_cols = [col for col in index_df.columns if "EditingIndex" in col]
    index_cols = [col for col in index_df.columns if col.endswith("EditingIndex")]
    # highest_noise = index_df.agg({col: "mean" for col in signal_col_names if col != "A2GEditingIndex"}).max()
    highest_noise = index_df.agg({col: "mean" for col in index_cols if col != signal_col}).max()
    return highest_noise
