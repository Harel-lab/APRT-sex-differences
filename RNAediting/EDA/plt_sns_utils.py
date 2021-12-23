from collections import Counter
from pathlib import Path
import inspect
import sys

from statsmodels.stats.multitest import fdrcorrection
from more_itertools import distinct_combinations
from scipy.stats import ttest_ind
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from EDA.pandas_utils import define_groups_order
from Statannot.statannot import add_stat_annotation


def legal_groups_in_category(df, category_col, group_min_repeats=3):
    groups = df[category_col].tolist()
    groups_repeats_counter = Counter(groups)  # num of repeats for each group in the experiment
    legal_groups = [group for group in groups_repeats_counter
                    if groups_repeats_counter[group] >= group_min_repeats]
    return legal_groups


def one_category_box_pairs(df, category_col, group_min_repeats=3):
    # category_vals = df[category_col].unique().tolist()
    # pairs = [(x, y) for i, x in enumerate(category_vals) for y in category_vals[i + 1:]]

    legal_groups = legal_groups_in_category(df, category_col, group_min_repeats=group_min_repeats)

    pairs = [(x, y)
             for i, x in enumerate(legal_groups)
             for y in legal_groups[i + 1:]]

    return pairs


def two_categories_box_pairs(df, category_1_col, category_2_col, group_min_repeats=3):
    # category_1_groups = df[category_1_col].unique().tolist()
    # category_2_groups = df[category_2_col].unique().tolist()

    category_1_groups = legal_groups_in_category(df, category_1_col, group_min_repeats=group_min_repeats)
    category_2_groups = legal_groups_in_category(df, category_2_col, group_min_repeats=group_min_repeats)

    all_vals = category_1_groups + category_2_groups
    positions = []
    for x, y in distinct_combinations(all_vals, r=2):
        position = (x, y)
        if x == y:
            continue
        if x in category_1_groups and y in category_1_groups:
            continue
        if x in category_2_groups and y in category_2_groups:
            continue
        positions.append(position)
    box_pairs = [(pos_a, pos_b) for pos_a, pos_b in distinct_combinations(positions, r=2)
                 if pos_a != pos_b]
    return box_pairs


def pairs_y_vals(df, pairs, x, y, hue=None, sample=None):
    y_vals = []
    for pair in pairs:
        a, b = pair
        if hue:
            a_x, a_hue, b_x, b_hue = *a, *b
            if sample:
                a_y_vals = (df.
                            loc[(df[x] == a_x) & (df[hue] == a_hue)].
                            groupby(sample).
                            mean().
                            loc[:, y].
                            to_numpy())
                b_y_vals = (df.
                            loc[(df[x] == b_x) & (df[hue] == b_hue)].
                            groupby(sample).
                            mean().
                            loc[:, y].
                            to_numpy())
            else:
                a_y_vals = (df.
                            loc[(df[x] == a_x) & (df[hue] == a_hue)].
                            loc[:, y].
                            to_numpy())
                b_y_vals = (df.
                            loc[(df[x] == b_x) & (df[hue] == b_hue)].
                            loc[:, y].
                            to_numpy())
        else:
            a_x, b_x = pair
            if sample:
                a_y_vals = (df.
                            loc[df[x] == a_x].
                            groupby(sample).
                            mean().
                            loc[:, y].
                            to_numpy())
                b_y_vals = (df.
                            loc[df[x] == b_x].
                            groupby(sample).
                            mean().
                            loc[:, y].
                            to_numpy())
            else:
                a_y_vals = df.loc[df[x] == a_x].loc[:, y].to_numpy()
                b_y_vals = df.loc[df[x] == b_x].loc[:, y].to_numpy()
        y_vals.append((a_y_vals, b_y_vals))
    return y_vals


def statannot_plot(df, x, y, pairs=None, hue=None, group_min_repeats=3, sample=None, x_order=None,
                   ascending_x_order=True, x_ordering_strategy="median", hue_order=None, ascending_hue_order=True,
                   hue_ordering_strategy="median", outer_test="t-test-ind", inner_test=None, main_title="",
                   text_format='star', dodge=True, base_plot="boxplot", points_plot="swarmplot", base_plot_color=None,
                   base_plot_pallete=None, points_plot_color="black", points_plot_pallete=None,
                   points_plot_alpha=0.6, ylabel=None, xlabel=None, top_ylim=None,
                   bottom_ylim=None, main_title_pad=20, xticklabels_rotation=None, xticklabels_fontsize=None,
                   xlables_pad=10, ylables_pad=10, barplot_capsize=.05, palette=None, points_plot_size=6,
                   plots_out_dir=None, dpi=600, scale_fig_width_inches_by=1.1, scale_fig_height_inches_by=1.1,
                   explicit_pvs_and_ns_theme=False, explicit_pvs_and_ns_theme_full_experience=False,
                   text_offset=1, line_offset=None, fontsize="medium", line_offset_to_box=None, line_height=0.02,
                   width=0.8, linewidth=None, legend_loc="lower right", legend_frameon=True, bbox_to_anchor=None):

    if explicit_pvs_and_ns_theme_full_experience:
        explicit_pvs_and_ns_theme = True
        text_offset = 1.5
        line_offset = 0.35
        fontsize = "x-small"
        line_offset_to_box = 0.15
        line_height = 0.07
        linewidth = 2
        legend_loc = "upper left"
        bbox_to_anchor = (1.03, 0.7)
        legend_frameon = False
        points_plot_size = 6

    if not x_order:
        x_order = define_groups_order(df=df, group_col=x, signal_col=y, ascending=ascending_x_order,
                                      strategy=x_ordering_strategy)
    if hue and not hue_order:
        hue_order = define_groups_order(df=df, group_col=hue, signal_col=y, ascending=ascending_hue_order,
                                        strategy=hue_ordering_strategy)
    if not pairs:
        if hue:
            pairs = two_categories_box_pairs(df, category_1_col=x, category_2_col=hue,
                                             group_min_repeats=group_min_repeats)
        else:
            pairs = one_category_box_pairs(df, category_col=x, group_min_repeats=group_min_repeats)

    if outer_test and inner_test:
        raise Exception("decide if to use an outer test or an inner one (that's included in statannot)")

    # make a grouped boxplot/barplot and save it in a variable
    if base_plot == "boxplot":
        ax = sns.boxplot(data=df, x=x, y=y, order=x_order,
                         hue=hue, hue_order=hue_order if hue else None,
                         color=base_plot_color if base_plot_color else None,
                         palette=base_plot_pallete if base_plot_pallete else None,
                         whis=np.inf, manage_ticks=True, width=width,
                         linewidth=linewidth if linewidth else None)
    elif base_plot == "barplot":
        ax = sns.barplot(data=df, x=x, y=y, order=x_order,
                         hue=hue, hue_order=hue_order if hue else None,
                         color=base_plot_color if base_plot_color else None,
                         palette=base_plot_pallete if base_plot_pallete else None,
                         capsize=barplot_capsize,
                         linewidth=linewidth if linewidth else None)
    else:
        raise Exception("choose either 'boxplot' or 'barplot' for base_plot")
    # make a grouped stripplot/swarmplot and save it in a variable
    if points_plot == "stripplot":
        ax = sns.stripplot(data=df, x=x, y=y, order=x_order,
                           hue=hue, hue_order=hue_order if hue else None,
                           palette=points_plot_pallete if points_plot_pallete else None,
                           color=points_plot_color, alpha=points_plot_alpha,
                           dodge=dodge, jitter=True, size=points_plot_size, ax=ax)
    elif points_plot == "swarmplot":
        ax = sns.swarmplot(data=df, x=x, y=y, order=x_order,
                           hue=hue, hue_order=hue_order if hue else None,
                           palette=points_plot_pallete if points_plot_pallete else None,
                           color=points_plot_color, alpha=points_plot_alpha,
                           dodge=dodge, size=points_plot_size, ax=ax)
    else:
        raise Warning("A proper name of function for plotting each sample wasn't supplied")

    # add stats
    if outer_test and len(pairs) > 0:
        y_vals = pairs_y_vals(df, pairs, x, y, hue, sample)
        p_vals = []
        for a_y_vals, b_y_vals in y_vals:
            if outer_test == "t-test-ind":
                _, p_val = ttest_ind(a_y_vals, b_y_vals)
            else:
                raise NotImplementedError("outer_test not implemented yet")
            p_vals.append(p_val)

        corrected_tests_rejected, corrected_p_values = fdrcorrection(p_vals)
        significant_box_pairs = []
        significant_p_vals = []
        for corrected_rejected, corrected_p_val, pair in zip(corrected_tests_rejected, corrected_p_values, pairs):
            if corrected_rejected:
                significant_box_pairs.append(pair)
                if explicit_pvs_and_ns_theme:
                    pvalue_thresholds = [[1e-4, "****"], [1e-3, "***"], [1e-2, "**"], [0.05, "*"]]
                    for threshold, stars in pvalue_thresholds:
                        if corrected_p_val <= threshold:
                            significant_p_vals.append(f"{stars}\n{corrected_p_val:.2e}")
                            break
                    else:
                        significant_p_vals.append(f"{corrected_p_val:.2e}")
                else:
                    significant_p_vals.append(corrected_p_val)
            else:
                if explicit_pvs_and_ns_theme:
                    significant_box_pairs.append(pair)
                    significant_p_vals.append("NS")

        if len(significant_box_pairs) == len(significant_p_vals) > 0:
            if explicit_pvs_and_ns_theme:
                ax, test_results = add_stat_annotation(ax, data=df, x=x, y=y, hue=hue, plot=base_plot,
                                                       box_pairs=significant_box_pairs,
                                                       order=x_order, hue_order=hue_order if hue else None,
                                                       perform_stat_test=False, text_annot_custom=significant_p_vals,
                                                       pvalues=[0 for _ in range(len(significant_p_vals))],
                                                       test_short_name=outer_test, text_format=text_format, verbose=2,
                                                       text_offset=text_offset, line_offset=line_offset,
                                                       fontsize=fontsize, line_offset_to_box=line_offset_to_box,
                                                       line_height=line_height)
            else:
                ax, test_results = add_stat_annotation(ax, data=df, x=x, y=y, hue=hue, plot=base_plot,
                                                       box_pairs=significant_box_pairs,
                                                       order=x_order, hue_order=hue_order if hue else None,
                                                       perform_stat_test=False, pvalues=significant_p_vals,
                                                       test_short_name=outer_test, text_format=text_format, verbose=2,
                                                       text_offset=text_offset, line_offset=line_offset,
                                                       fontsize=fontsize, line_offset_to_box=line_offset_to_box,
                                                       line_height=line_height)
    elif inner_test and len(pairs) > 0:
        ax, test_results = add_stat_annotation(ax, data=df, x=x, y=y, hue=hue, plot=base_plot,
                                               order=x_order, hue_order=hue_order if hue else None,
                                               box_pairs=pairs, test=inner_test, verbose=2, text_format=text_format,
                                               text_offset=text_offset, line_offset=line_offset, fontsize=fontsize,
                                               line_offset_to_box=line_offset_to_box, line_height=line_height)

    # set legend only for for hue
    if hue:
        if points_plot_pallete:  # make a custom legend
            handles = [mpatches.Patch(color=color, label=label) for color, label in zip(points_plot_pallete, hue_order)]
            if bbox_to_anchor:
                plt.legend(handles=handles,
                           bbox_to_anchor=bbox_to_anchor,
                           fontsize=fontsize,
                           frameon=legend_frameon,
                           loc=legend_loc,
                           title=hue)
            else:
                plt.legend(handles=handles,
                           fontsize=fontsize,
                           frameon=legend_frameon,
                           loc=legend_loc,
                           title=hue)
        else:  # use existing legend
            # get legend information from the plot object
            handles, labels = ax.get_legend_handles_labels()
            # specify just one legend
            if bbox_to_anchor:
                plt.legend(handles[0:2], labels[0:2],
                           bbox_to_anchor=bbox_to_anchor,
                           fontsize=fontsize,
                           frameon=legend_frameon,
                           loc='lower right',
                           title=hue)
            else:
                plt.legend(handles[0:2], labels[0:2],
                           fontsize=fontsize,
                           frameon=legend_frameon,
                           loc='lower right',
                           title=hue)

    # y axis title
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    # main title of the plot
    plt.title(main_title, fontdict={"fontsize": 14}, pad=main_title_pad)

    _, current_top_ylim = plt.ylim()
    if not bottom_ylim:
        bottom_ylim = 0
    if not top_ylim:
        top_ylim = current_top_ylim
        top_ylim += (top_ylim - bottom_ylim) * 0.1
    plt.ylim(bottom_ylim, top_ylim)

    # change xticklabels font size and/or rotate xticklables
    if xticklabels_fontsize and xticklabels_rotation:
        ax.set_xticklabels(fontdict={"fontsize": xticklabels_fontsize},
                           rotation=xticklabels_rotation,
                           labels=ax.get_xticklabels())
    elif xticklabels_fontsize:
        ax.set_xticklabels(fontdict={"fontsize": xticklabels_fontsize},
                           labels=ax.get_xticklabels())
    elif xticklabels_rotation:
        ax.set_xticklabels(rotation=xticklabels_rotation,
                           labels=ax.get_xticklabels())
    # pad x and ylabels
    ax.xaxis.labelpad, ax.yaxis.labelpad = xlables_pad, ylables_pad
    # remove left and upper borders
    sns.despine()

    # save plot
    if plots_out_dir is not None:

        fig = plt.gcf()
        width, height = fig.get_size_inches()
        width, height = width * scale_fig_width_inches_by, height * scale_fig_height_inches_by
        fig.set_size_inches(width, height)
        fig.tight_layout()

        plots_out_dir = Path(plots_out_dir)
        existing_plots = len(list(plots_out_dir.iterdir()))

        plt.savefig(Path(plots_out_dir,
                         f"{existing_plots + 1} {main_title}.svg"), dpi=dpi)

    print()

    plt.show()
