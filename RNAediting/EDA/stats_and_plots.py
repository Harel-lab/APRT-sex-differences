from statistics import stdev, mean, StatisticsError
from itertools import chain
from pathlib import Path
from math import ceil
import inspect
import sys

import numpy as np
import pandas as pd
from scipy.stats import ranksums, ttest_ind, ttest_rel
from statsmodels.stats.multitest import fdrcorrection
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from IPython.display import display_html

sys.path.append(str(Path(inspect.getfile(inspect.currentframe())).absolute().parent.parent))
from EDA.pandas_utils import sort_by_group_mean_val


def independent_t_test(group_a_df, group_b_df, signal_column, sample_column):
    group_a_means_df = group_a_df.loc[:, [sample_column, signal_column]].groupby(sample_column).mean()
    group_b_means_df = group_b_df.loc[:, [sample_column, signal_column]].groupby(sample_column).mean()
    group_a_means = group_a_means_df.loc[:, signal_column].to_numpy()
    group_b_means = group_b_means_df.loc[:, signal_column].to_numpy()
    _, p_val = ttest_ind(group_a_means, group_b_means)
    return p_val


def related_t_test(group_a_df, group_b_df, signal_column, sample_column):
    group_a_means_df = group_a_df.loc[:, [sample_column, signal_column]].groupby(sample_column).mean()
    group_b_means_df = group_b_df.loc[:, [sample_column, signal_column]].groupby(sample_column).mean()
    group_a_means = group_a_means_df.loc[:, signal_column].to_numpy()
    group_b_means = group_b_means_df.loc[:, signal_column].to_numpy()
    _, p_val = ttest_rel(group_a_means, group_b_means)
    return p_val


def wilcoxon_ranksums(group_a_df, group_b_df, signal_column, _):
    group_a_signals = group_a_df.loc[:, signal_column]
    group_b_signals = group_b_df.loc[:, signal_column]
    _, p_val = ranksums(group_a_signals, group_b_signals)
    return p_val


def display_two_pairwise_comparisons(groups_1_df, groups_2_df, group_column, signal_column, sample_column,
                                     test=independent_t_test, alpha=0.05, title_1="", title_2=""):
    styled_results_1_df, *_ = groups_pairwise_comparisons_2(groups_1_df, group_column, signal_column, sample_column,
                                                            test=test, alpha=alpha)
    styled_results_2_df, *_ = groups_pairwise_comparisons_2(groups_2_df, group_column, signal_column, sample_column,
                                                            test=test, alpha=alpha)

    styled_results_1_df = styled_results_1_df.set_table_attributes("style='display:inline'").set_caption(title_1)
    styled_results_2_df = styled_results_2_df.set_table_attributes("style='display:inline'").set_caption(title_2)

    space = "\xa0" * 10

    display_html(styled_results_1_df._repr_html_() + space + styled_results_2_df._repr_html_(), raw=True)


def groups_pairwise_comparisons_2(groups_df, group_column, signal_column, sample_column, test=independent_t_test,
                                  alpha=0.05):
    groups = groups_df[group_column].unique()
    groups_w_two_or_more_samples = []
    for group in groups:
        group_samples = groups_df.loc[groups_df[group_column] == group]
        num_of_samples = len(group_samples)
        if 2 <= num_of_samples:
            groups_w_two_or_more_samples.append(group)
    groups = groups_w_two_or_more_samples
    tests_df = pd.DataFrame({group: ["" for _ in range(len(groups))] for group in groups},
                            index=groups)

    for i, group_a in enumerate(groups):
        group_a_df = groups_df.loc[groups_df[group_column] == group_a]
        for group_b in groups[i+1:]:
            group_b_df = groups_df.loc[groups_df[group_column] == group_b]
            p_val = test(group_a_df, group_b_df, signal_column, sample_column)
            tests_df.loc[group_a, group_b] = p_val
    concatenated_tests_w_nan = np.concatenate(tests_df.to_numpy(), axis=None)
    concatenated_tests_wo_nan = [p_val for p_val in concatenated_tests_w_nan if p_val]  # filter out tests with nan

    a_1 = 1
    a_n = len(groups) - 1
    n = len(groups) - 1
    s = (n * (a_1 + a_n)) / 2
    if s != len(concatenated_tests_wo_nan):
        raise Exception("s != len(concatenated_tests_wo_nan)")

    corrected_rejected_tests, corrected_p_values = fdrcorrection(concatenated_tests_wo_nan, alpha=alpha)
    corrected_p_values_df = pd.DataFrame({group: ["" for _ in range(len(groups))] for group in groups},
                                         index=groups)
    corrected_rejected_tests_df = pd.DataFrame({group: ["" for _ in range(len(groups))] for group in groups},
                                               index=groups)
    corrected_rejected_tests_iter = iter(corrected_rejected_tests)
    corrected_p_values_iter = iter(corrected_p_values)
    for i, group_a in enumerate(groups):
        try:
            for group_b in groups[i+1:]:
                corrected_rejected_test = next(corrected_rejected_tests_iter)
                corrected_p_val = next(corrected_p_values_iter)
                corrected_p_values_df.loc[group_a, group_b] = corrected_p_val
                corrected_rejected_tests_df.loc[group_a, group_b] = corrected_rejected_test
        except StopIteration:
            break
    def significant_red(pval_rejection):
        color = 'red' if pval_rejection else 'black'
        return 'color: %s' % color
    styled_corrected_p_values_df = \
        corrected_p_values_df.style.apply(lambda x: corrected_rejected_tests_df.applymap(significant_red), axis=None)
    return styled_corrected_p_values_df, corrected_p_values_df, corrected_rejected_tests_df


def groups_pairwise_comparisons(groups_df, group_column, signal_column, alpha=0.05):
    groups = groups_df[group_column].unique()
    tests_df = pd.DataFrame({group: ["" for _ in range(len(groups))] for group in groups},
                            index=groups)
    for i, group_a in enumerate(groups):
        group_a_signals = groups_df.loc[groups_df[group_column] == group_a, signal_column]
        for group_b in groups[i+1:]:
            group_b_signals = groups_df.loc[groups_df[group_column] == group_b, signal_column]
            _, p_val = ranksums(group_a_signals, group_b_signals)
            tests_df.loc[group_a, group_b] = p_val
    concatenated_tests_w_nan = np.concatenate(tests_df.to_numpy(), axis=None)
    concatenated_tests_wo_nan = [p_val for p_val in concatenated_tests_w_nan if p_val]  # filter out tests with nan
    corrected_rejected_tests, corrected_p_values = fdrcorrection(concatenated_tests_wo_nan, alpha=alpha)
    corrected_p_values_df = pd.DataFrame({group: ["" for _ in range(len(groups))] for group in groups},
                                         index=groups)
    corrected_rejected_tests_df = pd.DataFrame({group: ["" for _ in range(len(groups))] for group in groups},
                                               index=groups)
    corrected_rejected_tests_iter = iter(corrected_rejected_tests)
    corrected_p_values_iter = iter(corrected_p_values)
    for i, group_a in enumerate(groups):
        try:
            for group_b in groups[i+1:]:
                corrected_rejected_test = next(corrected_rejected_tests_iter)
                corrected_p_val = next(corrected_p_values_iter)
                corrected_p_values_df.loc[group_a, group_b] = corrected_p_val
                corrected_rejected_tests_df.loc[group_a, group_b] = corrected_rejected_test
        except StopIteration:
            break
    def significant_red(pval_rejection):
        color = 'red' if pval_rejection else 'black'
        return 'color: %s' % color
    styled_corrected_p_values_df = \
        corrected_p_values_df.style.apply(lambda x: corrected_rejected_tests_df.applymap(significant_red), axis=None)
    return styled_corrected_p_values_df, corrected_p_values_df, corrected_rejected_tests_df


def all_a2g_and_cds_plots_and_tables(analyse_all_a2g_df, analyse_cds_a2g_df, all_a2g_df, cds_a2g_df,
                                     barplot_signal_column, pairwise_signal_column, group_column, main_title,
                                     barplot_y_title, subplot_titles, all_a2g_alpha=0.05, cds_a2g_alpha=0.05,
                                     width=1250, coloraxis_colorbar_title="", horizontal_spacing=0.04,
                                     vertical_spacing=0.11, table_font_size=10):
    groups = analyse_all_a2g_df[group_column].unique()

    # *** preparing bar plots of groups' mean signals ***

    # get mean signal of each group in each df (absent group will have a 0 value)
    all_a2g_mean_signals = pd.Series([analyse_all_a2g_df.loc[analyse_all_a2g_df[group_column] == group,
                                                             barplot_signal_column].mean(skipna=False)
                                      for group in groups])
    cds_a2g_mean_signals = pd.Series([analyse_cds_a2g_df.loc[analyse_cds_a2g_df[group_column] == group,
                                                             barplot_signal_column].mean(skipna=False)
                                      for group in groups])
    all_a2g_mean_signals = all_a2g_mean_signals.fillna(0)
    cds_a2g_mean_signals = cds_a2g_mean_signals.fillna(0)
    # couple mean signals and groups of each data set into a list of tuples, and sort that list by ascending order of
    # mean signals
    all_a2g_mean_signals_and_groups = [(signal, group) for signal, group in zip(all_a2g_mean_signals, groups)]
    cds_a2g_mean_signals_and_groups = [(signal, group) for signal, group in zip(cds_a2g_mean_signals, groups)]
    all_a2g_mean_signals_and_groups.sort(key=lambda x: x[0])
    cds_a2g_mean_signals_and_groups.sort(key=lambda x: x[0])
    # final x and y values for the bar plots
    final_general_signals = [signal for signal, _ in all_a2g_mean_signals_and_groups]
    final_general_groups = [group for _, group in all_a2g_mean_signals_and_groups]
    final_cds_signals = [signal for signal, _ in cds_a2g_mean_signals_and_groups]
    final_cds_groups = [group for _, group in cds_a2g_mean_signals_and_groups]

    # *** preparing pairwise comparisons tables ***

    # getting p-vals of tests after fdr correction
    _, all_p_vals_df, all_rejected_tests_df = groups_pairwise_comparisons(all_a2g_df, group_column,
                                                                          pairwise_signal_column, all_a2g_alpha)
    _, cds_p_vals_df, cds_rejected_tests_df = groups_pairwise_comparisons(cds_a2g_df, group_column,
                                                                          pairwise_signal_column, cds_a2g_alpha)
    # keeping only p-vals of rejected tests
    joined_p_vals_dfs = [all_p_vals_df, cds_p_vals_df]
    joined_rejected_tests_dfs = [all_rejected_tests_df, cds_rejected_tests_df]
    joined_ordered_groups = [final_general_groups, final_cds_groups]
    joined_final_dfs = [pd.DataFrame({group: ["" for _ in range(len(ordered_groups))] for group in ordered_groups},
                                     index=ordered_groups)
                        for ordered_groups in joined_ordered_groups]
    for p_vals_df, rejected_tests_df, ordered_groups, final_df in zip(joined_p_vals_dfs, joined_rejected_tests_dfs,
                                                                      joined_ordered_groups, joined_final_dfs):
        for i, group_a in enumerate(ordered_groups):
            for group_b in ordered_groups[i+1:]:
                try:
                    p_val = p_vals_df.loc[group_a, group_b]
                    test_rejected = rejected_tests_df.loc[group_a, group_b]
                    if test_rejected:
                        p_val = round(p_val, 3)
                        if p_val == 0:
                            p_val = "0"
                        else:
                            p_val = "." + str(p_val).split(".")[1]
                    else:
                        p_val = ""
                except KeyError:
                    p_val = ""
                final_df.loc[group_a, group_b] = p_val
    # convert joined_p_vals_dfs to lists of lists
    all_cols = [""] + [col for col in joined_final_dfs[0].columns]
    all_cells = [[col for col in joined_final_dfs[0].columns]]
    all_cells.extend(joined_final_dfs[0][col].tolist() for col in joined_final_dfs[0].columns)
    cds_cols = [""] + [col for col in joined_final_dfs[1].columns]
    cds_cells = [[col for col in joined_final_dfs[1].columns]]
    cds_cells.extend(joined_final_dfs[1][col].tolist() for col in joined_final_dfs[1].columns)

    #  *** finally, assemble the fig ***

    # create a fig object
    fig = make_subplots(
        rows=2, cols=2,
        shared_yaxes=True,
        horizontal_spacing=horizontal_spacing, vertical_spacing=vertical_spacing,
        specs=[[{"type": "bar"}, {"type": "bar"}], [{"type": "table"}, {"type": "table"}]],
        subplot_titles=subplot_titles
    )

    # create bar plots
    fig.add_trace(go.Bar(x=final_general_groups, y=final_general_signals,
                         marker=dict(color=final_general_signals, coloraxis="coloraxis")),
                  row=1, col=1)
    fig.add_trace(go.Bar(x=final_cds_groups, y=final_cds_signals,
                         marker=dict(color=final_cds_signals, coloraxis="coloraxis")),
                  row=1, col=2)

    # create pairwise comparisons tables
    letters = sum(len(col) for col in all_cols)
    letters += max(len(col) for col in all_cols)
    width_per_letter = int(width / letters) - 1
    all_column_width = [width_per_letter * len(col) for col in all_cols]
    all_column_width[0] = max(all_column_width) - 25
    cds_column_width = [width_per_letter * len(col) for col in cds_cols]
    cds_column_width[0] = max(cds_column_width) - 25
    fig.add_trace(go.Table(columnwidth=all_column_width,
                           header=dict(values=all_cols, font=dict(size=table_font_size), align="center", height=30),
                           cells=dict(values=all_cells, font=dict(size=table_font_size), align="center", height=30)),
                  row=2, col=1)
    fig.add_trace(go.Table(columnwidth=cds_column_width,
                           header=dict(values=cds_cols, font=dict(size=table_font_size), align="center", height=30),
                           cells=dict(values=cds_cells, font=dict(size=table_font_size), align="center", height=30)),
                  row=2, col=2)

    # update axes & layout
    fig.update_yaxes(title_text=barplot_y_title, row=1, col=1)
    fig.update_layout(title_text=main_title, title_x=0.47,
                      height=1200, width=width,
                      coloraxis=dict(colorscale='Bluered_r'),
                      coloraxis_colorbar=dict(len=0.4, lenmode="fraction", yanchor="top",
                                              y=1, ticks="", showticklabels=True,
                                              dtick=((int(max(final_general_signals) / 9) // 1000) * 1000),
                                              title=coloraxis_colorbar_title),
                      showlegend=False)

    return fig


def make_pretty(*, replacements=(("_", "<br>"), ("-", "<br>")), capitalize=True):
    if replacements or capitalize:
        if replacements:
            replacements = list(chain.from_iterable(replacements))
        def prettify_labels(labels):
            new_labels = []
            for label in labels:
                for x in range(0, len(replacements), 2):
                    old_exp, new_exp = replacements[x], replacements[x+1]
                    label = label.replace(old_exp, new_exp)
                if capitalize:
                    label = label.capitalize()
                new_labels.append(label)
            return new_labels
        return prettify_labels
    return None


def try_stdev(vals):
    try:
        return stdev(vals)
    except StatisticsError:
        return 0


def two_error_bar_plots(df_1, df_2, group_col, signal_col, *, width=800, height=500, horizontal_spacing=0.04,
                        vertical_spacing=0.13, coloraxis_colorbar_title="", main_title="", barplot_y_title="",
                        subplot_titles="", title_x=0.5, tickangle=0,
                        x_labels_replacements=(("_", "<br>"), ("-", "<br>")), capitalize_x_labels=True):

    df_1 = sort_by_group_mean_val(df_1, group_col, signal_col)
    df_2 = sort_by_group_mean_val(df_2, group_col, signal_col)

    groups_1 = df_1[group_col].unique()
    groups_2 = df_2[group_col].unique()

    prettify = make_pretty(replacements=x_labels_replacements, capitalize=capitalize_x_labels)
    x_vals_1 = prettify(groups_1) if prettify else groups_1
    x_vals_2 = prettify(groups_2) if prettify else groups_2

    # all_y_vals are nested lists:
    # each inner-list holds the signal values (as defined by "signal_col") of the samples of each group (as defined
    # by "group_col")
    all_y_vals_1 = [df_1.loc[df_1[group_col] == group, signal_col].to_list()
                    for group in groups_1]
    all_y_vals_2 = [df_2.loc[df_2[group_col] == group, signal_col].to_list()
                    for group in groups_2]

    sds_1 = [round(try_stdev(y_vals), 1) for y_vals in all_y_vals_1]
    sds_2 = [round(try_stdev(y_vals), 1) for y_vals in all_y_vals_2]

    mean_y_vals_1 = [mean(y_vals) for y_vals in all_y_vals_1]
    mean_y_vals_2 = [mean(y_vals) for y_vals in all_y_vals_2]

    fig = make_subplots(
        rows=1, cols=2,
        shared_yaxes=True,
        horizontal_spacing=horizontal_spacing, vertical_spacing=vertical_spacing,
        specs=[[{"type": "bar"}, {"type": "bar"}]],
        subplot_titles=subplot_titles)

    fig.add_trace(go.Bar(x=x_vals_1, y=mean_y_vals_1, error_y=dict(type='data', array=sds_1),
                         marker=dict(color=mean_y_vals_1, coloraxis="coloraxis")),
                  row=1, col=1)
    fig.add_trace(go.Bar(x=x_vals_2, y=mean_y_vals_2, error_y=dict(type='data', array=sds_2),
                         marker=dict(color=mean_y_vals_2, coloraxis="coloraxis")),
                  row=1, col=2)

    fig.update_yaxes(title_text=barplot_y_title, row=1, col=1)
    fig.update_xaxes(
        tickangle=tickangle,
        #         title_text = "Month",
        #         title_font = {"size": 20},
        title_standoff=35)
    fig.update_layout(title_text=main_title, title_x=title_x,
                      height=height, width=width,
                      coloraxis=dict(colorscale='Bluered_r'),
                      coloraxis_colorbar=dict(lenmode="fraction", yanchor="top",
                                              y=1, ticks="", showticklabels=True,
                                              title=coloraxis_colorbar_title),
                      showlegend=False)

    return fig


def scatter_error_bar_plot(df, group_col, signal_col, *, width=600, height=400, coloraxis_colorbar_title="",
                           main_title="", title_x=0.5, y_title="", tickangle=0, x_title="",
                           first_bar_coord=1, bar_width=1.5, bar_space=3, scatter_offset=-1.3, scatter_color="black",
                           scatter_symbol="circle", scatter_size=8,
                           x_labels_replacements=(("_", "<br>"), ("-", "<br>")), capitalize_x_labels=True):

    plot_df = df.copy()
    plot_df[group_col] = plot_df[group_col].astype(str)

    plot_df = sort_by_group_mean_val(plot_df, group_col, signal_col)
    groups = plot_df[group_col].unique()

    prettify = make_pretty(replacements=x_labels_replacements, capitalize=capitalize_x_labels)
    x_labels = prettify(groups) if prettify else groups

    # all_bar_y_vals is a nested lists:
    # each inner-list holds the signal values (as defined by "signal_col") of the samples of each group (as defined
    # by "group_col")
    all_bar_y_vals = [plot_df.loc[plot_df[group_col] == group, signal_col].to_list()
                      for group in groups]
    bar_x_coords = [x for x in range(first_bar_coord,
                                     int(len(all_bar_y_vals) * (bar_width + bar_space)),
                                     int(bar_width + bar_space))]
    bar_x_coords = bar_x_coords[:len(all_bar_y_vals)]
    bar_sds = [round(try_stdev(y_vals), 1) for y_vals in all_bar_y_vals]
    mean_bar_y_vals = [mean(y_vals) for y_vals in all_bar_y_vals]

    fig = go.Figure()

    fig.add_trace(go.Bar(x=bar_x_coords, y=mean_bar_y_vals, error_y=dict(type='data', array=bar_sds),
                         marker=dict(color=mean_bar_y_vals, coloraxis="coloraxis"),
                         width=[bar_width for _ in range(len(mean_bar_y_vals))]))

    scatter_x_coords = []
    scatter_y_vals = []
    for i, bar_x_coord in enumerate(bar_x_coords):
        for y_val in all_bar_y_vals[i]:
            scatter_x_coord = bar_x_coord + scatter_offset
            scatter_x_coords.append(scatter_x_coord)
            scatter_y_vals.append(y_val)

    fig.add_trace(go.Scatter(x=scatter_x_coords, y=scatter_y_vals,
                             mode='markers', marker_color=scatter_color,
                             marker=dict(line_width=1, symbol=scatter_symbol, size=scatter_size)))

    fig.update_yaxes(title_text=y_title)
    fig.update_xaxes(title_text=x_title, tickangle=tickangle)
    fig.update_layout(title_text=main_title, title_x=title_x,
                      height=height, width=width,
                      coloraxis=dict(colorscale='Bluered_r'),
                      coloraxis_colorbar=dict(lenmode="fraction", yanchor="top", y=1, ticks="",
                                              showticklabels=True, title=coloraxis_colorbar_title),
                      showlegend=False,
                      xaxis=dict(tickmode='array', tickvals=bar_x_coords, ticktext=x_labels))

    return fig


def multiple_scatter_error_bar_plots(dfs, group_cols, signal_col, *, rows=None, cols=None, max_cols=3, width=600,
                                     height=400, coloraxis_colorbar_title="", main_title="", title_x=0.5, y_title="",
                                     subplot_titles=(), tickangle=0, x_titles=(), x_titles_font_size=14,
                                     title_standoff=None,
                                     first_bar_coord=1, bar_width=1.5, bar_space=3, scatter_offset=-1.3,
                                     scatter_color="black", scatter_symbol="circle", scatter_size=6,
                                     horizontal_spacing=0.04, vertical_spacing=0.11,
                                     x_labels_replacements=(("_", "<br>"), ("-", "<br>")), capitalize_x_labels=True,
                                     coloraxis_colorscale="Bluered_r"):
    num_of_plots = len(dfs)
    if (rows and not cols) or (cols and not rows):
        raise Exception("rows and cols must be provided together or not at all")
    if not rows and not cols:
        cols = num_of_plots if num_of_plots <= max_cols else max_cols
        rows = ceil(num_of_plots / cols)

    # create a fig object
    fig = make_subplots(
        rows=rows, cols=cols,
        shared_yaxes=True,
        horizontal_spacing=horizontal_spacing, vertical_spacing=vertical_spacing,
        subplot_titles=subplot_titles,
        y_title=y_title
    )

    prettify = make_pretty(replacements=x_labels_replacements, capitalize=capitalize_x_labels)

    for plot_num in range(1, num_of_plots + 1):
        col = ((plot_num - 1) % cols) + 1
        row = ceil(plot_num / cols)
        group_col = group_cols[plot_num - 1]
        df = dfs[plot_num - 1].copy()
        df = sort_by_group_mean_val(df, group_col, signal_col)
        groups = df[group_col].unique()
        x_title = x_titles[plot_num - 1] if x_titles else ""
        plot_x_labels = prettify(groups) if prettify else groups

        all_bar_y_vals = [df.loc[df[group_col] == group, signal_col].to_list()
                          for group in groups]
        bar_x_coords = [x for x in range(first_bar_coord,
                                         int(len(all_bar_y_vals) * (bar_width + bar_space)),
                                         int(bar_width + bar_space))]
        bar_x_coords = bar_x_coords[:len(all_bar_y_vals)]
        bar_sds = [round(try_stdev(y_vals), 1) for y_vals in all_bar_y_vals]
        mean_bar_y_vals = [mean(y_vals) for y_vals in all_bar_y_vals]

        scatter_x_coords = []
        scatter_y_vals = []
        for i, bar_x_coord in enumerate(bar_x_coords):
            for y_val in all_bar_y_vals[i]:
                scatter_x_coord = bar_x_coord + scatter_offset
                scatter_x_coords.append(scatter_x_coord)
                scatter_y_vals.append(y_val)

        fig.add_trace(go.Bar(x=bar_x_coords, y=mean_bar_y_vals, error_y=dict(type='data', array=bar_sds),
                             marker=dict(color=mean_bar_y_vals, coloraxis="coloraxis"),
                             width=[bar_width for _ in range(len(mean_bar_y_vals))]),
                      row=row, col=col)
        fig.add_trace(go.Scatter(x=scatter_x_coords, y=scatter_y_vals,
                                 mode='markers', marker_color=scatter_color,
                                 marker=dict(line_width=1, symbol=scatter_symbol, size=scatter_size)),
                      row=row, col=col)
        fig.update_xaxes(title_text=x_title, tickangle=tickangle,
                         tickmode='array', tickvals=bar_x_coords, ticktext=plot_x_labels,
                         title_font=dict(size=x_titles_font_size),
                         row=row, col=col)
        if title_standoff:
            fig.update_xaxes(title_standoff=title_standoff)

    fig.update_layout(title_text=main_title, title_x=title_x,
                      height=height, width=width,
                      coloraxis=dict(colorscale=coloraxis_colorscale),
                      showlegend=False,
                      coloraxis_colorbar=dict(lenmode="fraction", yanchor="top",
                                              y=1, ticks="", showticklabels=True,
                                              title=coloraxis_colorbar_title)
                      )
    if 1 < rows:
        fig.update_layout(
            coloraxis_colorbar=dict(lenmode="fraction", yanchor="middle",
                                    len=(1 / rows) + (vertical_spacing / 2),
                                    # len=(1 / rows),
                                    y=0.5, ticks="", showticklabels=True,
                                    title=coloraxis_colorbar_title)
        )

    return fig



