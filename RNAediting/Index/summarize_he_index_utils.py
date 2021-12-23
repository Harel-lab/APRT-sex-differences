from pathlib import Path
import inspect
import sys
import re

import pandas as pd

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from Index.utils import find_highest_noise, signal_noise_ratio, read_index_df


def optimize_signal_vs_snr(he_index_main_out_dir):
    """
    Calc total events, unique sites signal and SNR of each, after sequential cluster screening.

    @param he_index_main_out_dir: folder with different Index runs
    @type he_index_main_out_dir: Path
    @return: signals_snrs_df, a table containing distances, A2G signals and SNRs
    @rtype: pandas.DataFrame
    """
    sub_out_dirs = [sub_out_dir for sub_out_dir in he_index_main_out_dir.iterdir()
                    if sub_out_dir.is_dir() and re.match(r"HEIndex[0-9]+", sub_out_dir.name)]

    distances = [int(sub_out_dir.name.lstrip("HEIndex")) for sub_out_dir in sub_out_dirs]

    signals = []
    highest_noises = []
    snrs = []
    actual_signals = []

    index_dfs = [read_index_df(Path(sub_out_dir, "EditingIndex.csv")) for sub_out_dir in sub_out_dirs]

    for index_df in index_dfs:

        a2g_signal = index_df["A2GEditingIndex"].mean()
        highest_noise = find_highest_noise(index_df)
        snr = signal_noise_ratio(a2g_signal, highest_noise)
        actual_signal = a2g_signal - highest_noise

        signals.append(a2g_signal)
        highest_noises.append(highest_noise)
        snrs.append(snr)
        actual_signals.append(actual_signal)

    signals_snrs_df = pd.DataFrame({"Distance": distances, "Signal": signals, "Noise": highest_noises,
                                    "SNR": snrs, "ActualSignal": actual_signals})

    return signals_snrs_df


def find_best_distance(signals, snrs, distances):
    """
    Find best distance to extend UE regions by, considering signal vs. signal-to-noise ratio (SNR).

    @param signals: signals of indices, whose UE regions were extended by different number of BP from each side
    @type signals: list
    @param snrs: snrs of those indices
    @type snrs: list
    param distances: the corresponding distances
    @type distances: list
    @return: best_signal - signal of best index
    @rtype: float
    @return: best_snr - SNR of chosen cluster
    @rtype: float
    @return best_distance - distance used to create best cluster
    @rtype: int
    """
    best_signal = signals[0]
    best_snr = snrs[0]
    best_distance = distances[0]
    for signal, snr, distance in zip(signals[1:], snrs[1:], distances[1:]):
        if best_snr < snr:
            best_signal = signal
            best_snr = snr
            best_distance = distance
    return best_signal, best_snr, best_distance
