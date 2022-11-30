#!/usr/bin/env python3
"""
File:           visualizer.py
Description:    This file contains multiple functions that can be used for
                visualising data.
"""

import os
import math
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.optimize import curve_fit
from prospr import Protein
from prospr.visualize import plot_protein

from .functions import probit_func, linear_func
from .statistics import get_extrema
from .point import Point


# Colors and lengths in order of lengths.
LENGTHS = [30, 25, 20, 15, 10]
COLORS_SORTED = ["#5875A4", "#CC8963", "#5F9E6E", "#B55D60", "#857AAB"]
DUTCH_GLORY = ["#C6676B", "#F7F6F3", "#557DBE", "#DB8554"]


def _style_hist_plot(
    job, lengths, column, fig, title, ylim, show=True, save=False
):
    """
    Style histogram plot to usual style. Show and save according to arguments.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param [int]    lengths:    Lengths that were plotted.
    :param str      column:     The column from the results to use as data.
    :param Figure   fig:        Figure to style.
    :param str      title:      Title to give the plot.
    :param float    ylim:       Maximum value that is plotted.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    ax = fig.gca()
    ax.set_yscale("log")

    # Set y-axis label based on selected column.
    if column == "placed":
        ax.set_ylabel("Number of recursions", fontsize=11)
        ax.set_ylim([1, 10 ** (math.ceil(math.log(ylim, 10)) + 0.5)])
    elif column == "checked":
        ax.set_ylabel("Number of solutions checked", fontsize=11)
        ax.set_ylim([1, 10 ** (math.floor(math.log(ylim, 10)) + 0.3)])

    # Set labels that are non length-dependent.
    ax.set_xlabel("Problem instance rank", fontsize=11)
    ax.tick_params(labelsize=10)

    # Set hatch size smaller.
    plt.rcParams["hatch.linewidth"] = 0.15

    # Fix legend and add to plot.
    h, l = ax.get_legend_handles_labels()
    new_h = h[-len(lengths) :]
    new_l = l[-len(lengths) :]
    new_h.append(h[-len(lengths) - 1])
    new_l.append(l[-len(lengths) - 1])

    # Add handle for fitted function if present.
    if l[0] != "Quartiles" and l[0][:6] != "Length":
        new_h.append(h[0])
        new_l.append(l[0])

    ax.legend(new_h, new_l, loc="upper left", fontsize=10)
    plt.tight_layout()

    # Save plot if specified, title is length dependent.
    if save:
        if len(lengths) == 1:
            fpath = f"single_{lengths[0]}"
        else:
            fpath = f"group_{lengths}"

        dirpath = f"figures/{job}"
        Path(dirpath).mkdir(parents=True, exist_ok=True)

        # Determine filename based on selected column, and presence of a fitted
        # function handle.
        fextension = ""
        if l[0] != "Quartiles" and l[0][:6] != "Length":
            fextension = "_fitted"

        if column == "placed":
            plt.savefig(
                f"{dirpath}/{fpath}_placed{fextension}.png",
                format="png",
                dpi=600,
            )
        elif column == "checked":
            plt.savefig(
                f"{dirpath}/{fpath}_checked{fextension}.png",
                format="png",
                dpi=600,
            )

        print(f"\n~ Saved figure in {dirpath}/")

    # Show plot if specified.
    if show:
        plt.show()


def plot_single(job, length, column="placed", show=True, save=False):
    """
    Plot a barplot with checked solutions per protein_id. This can be done
    with or without sorting by the number of checked solutions.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param int      length:     Length to use.
    :param str      column:     The column from the results to use as data.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    root = f"jobs/{job}/results"

    # Load source file into a DataFrame.
    filename = [
        f
        for f in os.listdir(root)
        if f.startswith("HP") and int(f.split("_")[1]) == length
    ][0]
    source_df = pd.read_csv(f"{root}/{filename}")

    # Print statistics on the loaded data.
    print("Values stats:")
    print(f"\tMin:             {source_df[column].min()}")
    print(f"\tMax:             {source_df[column].max()}")
    print(f"\tMean:            {source_df[column].mean()}")
    print(f"\tStd:             {source_df[column].std()}\n")

    # Compute values for the bar chart.
    ys = sorted(source_df[column].values, reverse=False)
    xs = np.array(range(1, len(ys) + 1))

    # Create figure, plot bar chart and plot fitted functions.
    fig = plt.figure(figsize=(7.5, 4.2), dpi=600)
    sns.set(style="whitegrid")
    plt.bar(
        xs,
        ys,
        width=1,
        lw=0.05,
        label="Computational costs",
        color=COLORS_SORTED[LENGTHS.index(length) - 1],
    )

    # Add 1st and last quartile lines.
    plt.axvline(len(xs) / 4, color="black", ls=":", lw=1.5, label="Quartiles")
    plt.axvline(len(xs) / 4 * 3, color="black", ls=":", lw=1.5)

    # Style, save, and show plot according to arguments.
    title = f"Computational costs per problem instance for length {length}"
    _style_hist_plot(
        job, [length], column, fig, title, max(ys), show=show, save=save
    )


def plot_group(job, lengths, column="placed", show=True, save=False):
    """
    Plot a barplot with checked solutions per protein_id for multiple lengths.
    This can be done with or without sorting by the number of checked
    solutions.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    root = f"jobs/{job}/results"

    # Load source file into a DataFrame.
    lengths.sort(reverse=True)
    filenames = [
        f
        for length in lengths
        for f in os.listdir(root)
        if f.startswith("HP") and int(f.split("_")[1]) == length
    ]
    source_dfs = [pd.read_csv(f"{root}/{f}") for f in filenames]

    fig = plt.figure(figsize=(7.5, 4.2), dpi=600)
    sns.set(style="whitegrid")
    max_ys = []

    # Plot each length, and store largest y value for ylim.
    for source_df, length, color in zip(
        source_dfs,
        lengths,
        COLORS_SORTED,
    ):
        # Print statistics on the loaded data.
        print(f"Values stats for length: {length}")
        print(f"\tMin:             {source_df[column].min()}")
        print(f"\tMax:             {source_df[column].max()}")
        print(f"\tMean:            {source_df[column].mean()}")
        print(f"\tStd:             {source_df[column].std()}\n")

        ys = sorted(source_df[column].values, reverse=False)
        max_ys.append(max(ys))
        xs = np.array(range(1, len(ys) + 1))
        plt.bar(
            xs,
            ys,
            width=1,
            linewidth=0.05,
            label=f"Length {length}",
            color=color,
        )

    # Add 1st and last quartile lines.
    plt.axvline(len(xs) / 4, color="black", ls=":", lw=1.5, label="Quartiles")
    plt.axvline(len(xs) / 4 * 3, color="black", ls=":", lw=1.5)

    # Style, save, and show plot according to arguments.
    title = "Computational costs per problem instance for all lengths"
    _style_hist_plot(
        job, lengths, column, fig, title, max(max_ys), show=show, save=save
    )


def plot_fit_single(job, length, column="placed", show=True, save=False):
    """
    Plot a barplot with checked solutions per protein_id. This can be done
    with or without sorting by the number of checked solutions.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param int      length:     Length to use.
    :param str      column:     The column from the results to use as data.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    root = f"jobs/{job}/results"

    # Load source file into a DataFrame.
    filename = [
        f
        for f in os.listdir(root)
        if f.startswith("HP") and int(f.split("_")[1]) == length
    ][0]
    source_df = pd.read_csv(f"{root}/{filename}")

    # Print statistics on the loaded data.
    print("Values stats:")
    print(f"\tMin:             {source_df[column].min()}")
    print(f"\tMax:             {source_df[column].max()}")
    print(f"\tMean:            {source_df[column].mean()}")
    print(f"\tStd:             {source_df[column].std()}\n")

    ys = sorted(source_df[column].values, reverse=False)
    xs = np.array(range(1, len(ys) + 1))

    # Fit data to the scaled probit function.
    # Scale x-data to [0,1] and log-scale y-data for fit.
    machine_precision = np.finfo(float).eps
    unit_xs = np.linspace(machine_precision, 1 - machine_precision, len(ys))
    log_ys = np.log(ys)
    p0 = [math.sqrt(2), log_ys[xs[len(xs) // 2]]]
    popt = curve_fit(probit_func, unit_xs, log_ys, p0=p0)[0]

    print(f"Optimal params:\n{popt}")

    # Generate Y-data from fit and rescale to original size.
    fit_xs = np.linspace(machine_precision, 1 - machine_precision, 5000)
    fit_ys_probit = [probit_func(x, *popt) for x in fit_xs]
    fit_ys = np.power(np.e, fit_ys_probit)
    fit_xs_scaled = np.linspace(1, len(ys), 5000)

    fig = plt.figure(figsize=(7.5, 4.2), dpi=600)
    sns.set(style="whitegrid")

    plt.bar(
        xs,
        ys,
        width=1,
        linewidth=0.05,
        label=f"Length {length}",
        color=COLORS_SORTED[LENGTHS.index(length) - 1],
    )
    plt.plot(fit_xs_scaled, fit_ys, color="orange", label="Probit")

    # Add 1st and last quartile lines.
    plt.axvline(len(xs) / 4, color="black", ls=":", lw=1.5, label="Quartiles")
    plt.axvline(len(xs) / 4 * 3, color="black", ls=":", lw=1.5)

    # Style, save, and show plot according to arguments.
    title = f"Probit fitted on the instance hardness for length {length}"
    _style_hist_plot(
        job, [length], column, fig, title, max(ys), show=show, save=save
    )


def plot_fit_group(job, lengths, column="placed", show=True, save=False):
    """
    Plot a barplot with checked solutions per protein_id. This can be done
    with or without sorting by the number of checked solutions.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    root = f"jobs/{job}/results"

    # Load source file into a DataFrame.
    lengths.sort(reverse=True)
    filenames = [
        f
        for length in lengths
        for f in os.listdir(root)
        if f.startswith("HP") and int(f.split("_")[1]) == length
    ]
    source_dfs = [pd.read_csv(f"{root}/{f}") for f in filenames]

    fig = plt.figure(figsize=(7.5, 4.2), dpi=600)
    sns.set(style="whitegrid")
    max_ys = []

    # Plot each length, and store largest y value for ylim.
    for source_df, length, color in zip(
        source_dfs,
        lengths,
        DUTCH_GLORY,
    ):
        # Print statistics on the loaded data.
        print(f"Values stats for length: {length}")
        print(f"\tMin:             {source_df[column].min()}")
        print(f"\tMax:             {source_df[column].max()}")
        print(f"\tMean:            {source_df[column].mean()}")
        print(f"\tStd:             {source_df[column].std()}\n")

        # Plot the regular histogram.
        ys = sorted(source_df[column].values, reverse=False)
        max_ys.append(max(ys))
        xs = np.array(range(1, len(ys) + 1))
        plt.bar(
            xs,
            ys,
            width=1,
            linewidth=0.05,
            label=f"Length {length}",
            color=color,
        )

        # Fit data to the scaled probit function.
        # Scale x-data to [0,1] and log-scale y-data for fit.
        machine_precision = np.finfo(float).eps
        unit_xs = np.linspace(
            machine_precision, 1 - machine_precision, len(ys)
        )
        log_ys = np.log(ys)
        p0 = [math.sqrt(2), log_ys[xs[len(xs) // 2]]]
        popt = curve_fit(probit_func, unit_xs, log_ys, p0=p0)[0]

        print(f"Optimal params:\n{popt}")

        # Generate Y-data from fit and rescale to original size.
        fit_xs = np.linspace(machine_precision, 1 - machine_precision, 5000)
        fit_ys_probit = [probit_func(x, *popt) for x in fit_xs]
        fit_ys = np.power(np.e, fit_ys_probit)
        fit_xs_scaled = np.linspace(1, len(ys), 5000)

        # plt.plot(fit_xs_scaled, fit_ys, color="orange", label="Probit")
        plt.plot(
            fit_xs_scaled, fit_ys, color="mediumaquamarine", label="Probit"
        )

    # Add 1st and last quartile lines.
    plt.axvline(len(xs) / 4, color="black", ls=":", lw=1.5, label="Quartiles")
    plt.axvline(len(xs) / 4 * 3, color="black", ls=":", lw=1.5)

    # Style, save, and show plot according to arguments.
    title = "Probit fitted on the instance hardness for all lengths"
    _style_hist_plot(
        job, lengths, column, fig, title, max(max_ys), show=show, save=save
    )


def plot_params(job, lengths, column="placed", show=True, save=False):
    """
    Fit a Probit to each length, store the best parameters and create a plot
    showing the development of the params across lengths.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    root = f"jobs/{job}/results"

    # Load source file into a DataFrame.
    lengths.sort()
    filenames = [
        f
        for length in lengths
        for f in os.listdir(root)
        if f.startswith("HP") and int(f.split("_")[1]) == length
    ]
    source_dfs = [pd.read_csv(f"{root}/{f}") for f in filenames]
    popts = pd.DataFrame(columns=["alpha", "beta"])

    # Plot each length, and store largest y value for ylim.
    for source_df, length in zip(
        source_dfs,
        lengths,
    ):
        ys = sorted(source_df[column].values, reverse=False)
        xs = np.array(range(1, len(ys) + 1))

        # Fit data to the scaled probit function.
        # Scale x-data to [0,1] and log-scale y-data for fit.
        machine_precision = np.finfo(float).eps
        unit_xs = np.linspace(
            machine_precision, 1 - machine_precision, len(ys)
        )
        log_ys = np.log10(ys)
        p0 = [math.sqrt(2), log_ys[xs[len(xs) // 2]]]
        popt = curve_fit(probit_func, unit_xs, log_ys, p0=p0)[0]

        # Print statistics on the loaded data.
        print(f"Values stats for length: {length}")
        print(f"\tMin:             {source_df[column].min()}")
        print(f"\tMax:             {source_df[column].max()}")
        print(f"\tMean:            {source_df[column].mean()}")
        print(f"\tStd:             {source_df[column].std()}")
        print(f"\tOptimal params:  {popt}\n")
        popts.loc[length] = {
            "alpha": popt[0],
            "beta": popt[1],
        }

    # Fit the fitted alpha parameter with a linear function and report MSE
    # against the popts.
    alpha_ys = [popt for popt in popts["alpha"]]
    p0 = [alpha_ys[1] / alpha_ys[0], alpha_ys[1] - alpha_ys[0]]
    popt = curve_fit(linear_func, lengths, alpha_ys, p0=p0)[0]
    lin_ys = [linear_func(x, *popt) for x in lengths]
    popts["alpha_fit"] = lin_ys
    mse = ((np.array(alpha_ys) - np.array(lin_ys)) ** 2).mean()
    print(f"Alpha popt:  {popt}")
    print(f"MSE:         {mse}\n")

    # Fit the fitted beta parameter with a linear function and report MSE
    # against the popts.
    beta_ys = [popt for popt in popts["beta"]]
    p0 = [beta_ys[1] / beta_ys[0], beta_ys[1] - beta_ys[0]]
    popt = curve_fit(linear_func, lengths, beta_ys, p0=p0)[0]
    lin_ys = [linear_func(x, *popt) for x in lengths]
    popts["beta_fit"] = lin_ys
    mse = ((np.array(beta_ys) - np.array(lin_ys)) ** 2).mean()
    print(f"Beta popt:   {popt}")
    print(f"MSE:         {mse}")

    # Plot the lines for each parameter.
    plt.figure(figsize=(4, 3.3), dpi=600)
    sns.set_style("whitegrid")
    sns.lineplot(
        popts[["alpha", "beta"]],
        lw=1.6,
        palette=["#4c72b0", "#dd8452"],
        markers=["X", "X"],
        markersize=7,
    )
    ax = sns.lineplot(
        popts[["alpha_fit", "beta_fit"]],
        lw=1.8,
        palette=["#55a868", "#c44e52"],
    )
    [line.set_linestyle("-") for line in ax.lines[:2]]
    [line.set_linestyle("--") for line in ax.lines[2:]]

    # Add text above the line points.
    for l in lengths:
        plt.text(
            l,
            popts.loc[l]["alpha"] + 0.25,
            f"{popts.loc[l]['alpha']:.2f}",
            ha="center",
            va="center",
            fontsize=9,
        )

        if l != lengths[-1]:
            plt.text(
                l,
                popts.loc[l]["beta"] + 0.3,
                f"{popts.loc[l]['beta']:.2f}",
                ha="center",
                va="center",
                fontsize=9,
            )
        else:
            plt.text(
                l,
                popts.loc[l]["beta"] + 0.2,
                f"{popts.loc[l]['beta']:.2f}",
                ha="center",
                va="center",
                fontsize=9,
            )

    # Style the plot.
    plt.xticks(popts.index.values.astype(int))
    plt.xlabel("Protein length", fontsize=10)
    plt.ylabel("Parameter value", fontsize=10)
    plt.tick_params(labelsize=9)
    plt.tick_params(axis="both", which="major", pad=-1)
    plt.ylim(0, 6.999)

    # Fix style of alpha and beta. Set labels for fit.
    h, l = ax.get_legend_handles_labels()
    h[0].set_linestyle("-")
    h[1].set_linestyle("-")
    l[2] = "0.035x + 0.357"
    l[3] = "0.269x - 0.127"
    plt.legend(h, l, loc="upper left", fontsize=9)
    plt.tight_layout()

    # Save figure if specified.
    if save:
        dirpath = f"figures/{job}"
        Path(dirpath).mkdir(parents=True, exist_ok=True)

        plt.savefig(
            f"{dirpath}/params_{lengths}.pdf",
            format="pdf",
            dpi=600,
        )
        print(f"\n~ Saved figure in {dirpath}/")

    # Show figure if specified.
    if show:
        plt.show()


def plot_dataset_boxplot(dataset, lengths, show=True, save=False):
    """
    Plot a boxplot with H-ratio values per specified group.
    :param func     dataset:    Function call to loading the data.
    :param [int]    lengths:    Lengths to use.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    # Create DataFrame with all sequences for each length.
    lengths.sort()
    aminos = dataset(lengths[0])
    aminos["length"] = lengths[0]

    for l in lengths[1:]:
        df = dataset(l)
        df["length"] = l
        aminos = pd.concat([aminos, df])

    # Compute the H-ratio for each sequence.
    aminos = aminos.set_index("length")
    aminos["H-ratio"] = aminos.sequence.str.count("H")
    aminos["H-ratio"] = aminos["H-ratio"].div(aminos.index.to_series(), axis=0)
    print(aminos)

    # Create boxplot for each length.
    sns.set(style="whitegrid")
    ax = sns.boxplot(x=aminos.index.to_series(), y="H-ratio", data=aminos)
    ax.set_xlabel("Protein length")
    ax.set_ylabel("H-ratio")
    ax.set_title("H-ratio of the generated sequences per protein length")

    # Save figure if specified.
    if save:
        dirpath = f"figures/{dataset.__name__[5:]}"
        Path(dirpath).mkdir(parents=True, exist_ok=True)

        plt.savefig(
            f"{dirpath}/boxplot_{lengths}.svg",
            format="svg",
            dpi=600,
        )
        print(f"\n~ Saved figure in {dirpath}/")

    # Show figure if specified.
    if show:
        plt.show()


def plot_easiest_hardest_small(
    job, dataset, lengths, column="placed", show=True, save=False
):
    """
    Plot the conformation of the easiest and hardest problem instances ordered
    by their protein length.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param func     dataset:    Function call to loading the data.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    easiest_hardest = get_extrema(job, dataset, lengths, column=column)
    lengths.sort(reverse=False)

    # Create base for the subplots.
    fig, axs = plt.subplots(
        len(lengths), 2, figsize=(4 * 2, 2.5 * len(lengths))
    )

    for length, (ax_left, ax_right) in zip(lengths, axs):
        easy_p = easiest_hardest[length]["easiest"][0]
        hard_p = easiest_hardest[length]["hardest"][0]

        plot_protein(easy_p, style="paper", ax=ax_left, show=False)
        plot_protein(hard_p, style="paper", ax=ax_right, show=False)

    # Show if specified.
    if show:
        plt.show()


def plot_easiest_hardest(
    job, dataset, lengths, column="placed", show=True, save=False
):
    """
    Plot the conformation of the easiest and hardest problem instances ordered
    by their protein length.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param func     dataset:    Function call to loading the data.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    nprot = 5
    easiest_hardest = get_extrema(
        job, dataset, lengths, column=column, n=nprot
    )
    lengths.sort(reverse=False)

    # Create base for the subplots.
    fig, axs = plt.subplots(
        len(lengths) * 2, nprot, figsize=(4 * nprot, 2.5 * len(lengths) * 2)
    )

    for ax_i, length in zip(range(0, len(lengths) * 2, 2), lengths):
        easiest_p = easiest_hardest[length]["easiest"]
        hardest_p = easiest_hardest[length]["hardest"]

        for p_idx in range(nprot):
            plot_protein(
                easiest_p[p_idx],
                style="paper",
                ax=axs[ax_i][p_idx],
                show=False,
            )
            plot_protein(
                hardest_p[p_idx],
                style="paper",
                ax=axs[ax_i + 1][p_idx],
                show=False,
            )

    plt.tight_layout()

    # Save if specified.
    if save:
        plt.savefig(
            f"figures/{job}/easiest_hardest_{lengths}_n{nprot}.png",
            format="png",
            dpi=600,
        )
        print(f"\n~ Saved figure in figures/{job}/")

    # Show if specified.
    if show:
        plt.show()


def plot_conformation(
    job, dataset, lengths, column="placed", show=True, save=False
):
    """
    Plot the conformation of the hardest problem instances ordered by their
    protein length.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param func     dataset:    Function call to loading the data.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param bool     show:       True if plot.show() needs to be called.
    :param bool     save:       True if figure needs to be saved.
    """
    nprot = 2
    easiest_hardest = get_extrema(
        job, dataset, lengths, column=column, n=nprot, add_placed=True
    )
    lengths.sort(reverse=False)

    fig = plt.figure(dpi=600, figsize=(10, 7))
    gs = fig.add_gridspec(
        len(lengths), nprot, hspace=0.07, wspace=0.04, width_ratios=[2, 4]
    )
    (ax1, ax2), (ax3, ax4) = gs.subplots()
    axs = [[ax1, ax2], [ax3, ax4]]

    # Dictionary for rotating an image 90 degrees.
    transform_90deg = {-2: -1, -1: 2, 1: -2, 2: 1}

    for i, length in enumerate(lengths):
        easiest_p = easiest_hardest[length]["easiest"]
        hardest_p = easiest_hardest[length]["hardest"]

        for p_idx in range(nprot):
            p_easy = easiest_p[0][0]
            p_hard = hardest_p[1][0]

            # Top right window has different axis.
            if i == 0 and p_idx == 0:
                plot_protein(
                    p_easy,
                    style="paper",
                    ax=axs[p_idx][i],
                    legend_style="inner",
                    show=False,
                )
                axs[p_idx][i].set_xticks([-2.5, -2.0, -1.0, 0.0, 0.5])
                axs[p_idx][i].set_yticks(
                    [-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
                )
            elif i == 0 and p_idx == 1:
                new_hash = [transform_90deg[d] for d in p_hard.hash_fold()]
                p_hard.set_hash(new_hash, track=False)
                plot_protein(
                    p_hard,
                    style="paper",
                    ax=axs[p_idx][i],
                    legend=False,
                    show=False,
                )
                axs[p_idx][i].set_xticks(
                    [
                        -4.0,
                        -3,
                        -2.0,
                        -1,
                        0.0,
                        1,
                        2.0,
                        3,
                        4.0,
                        5,
                        6.0,
                        7,
                        8.0,
                        9,
                        10.0,
                        11,
                        12.0,
                    ]
                )
                axs[p_idx][i].set_yticks(
                    [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
                )
            # Create one axis with legend.
            elif i == 1 and p_idx == 1:
                new_hash = [transform_90deg[d] for d in p_hard.hash_fold()]
                p_hard.set_hash(new_hash, track=False)
                plot_protein(
                    p_hard,
                    style="paper",
                    legend=False,
                    ax=axs[i][p_idx],
                    show=False,
                )
                axs[p_idx][i].set_xticks(
                    [
                        -2.0,
                        -1,
                        0.0,
                        1,
                        2.0,
                        3,
                        4.0,
                        5,
                        6.0,
                        7,
                        8.0,
                        9,
                        10.0,
                        11,
                        12.0,
                        13,
                        14.0,
                    ]
                )
                axs[p_idx][i].set_yticks([-1.5, -1.0, 0, 1.0, 2.0, 2.5])
            else:
                new_hash = [transform_90deg[d] for d in p_easy.hash_fold()]
                p_easy.set_hash(new_hash, track=False)
                plot_protein(
                    p_easy,
                    style="paper",
                    ax=axs[0][1],
                    legend=False,
                    show=False,
                )
                axs[0][1].set_xticks(
                    [
                        -4.0,
                        -3,
                        -2.0,
                        -1,
                        0.0,
                        1,
                        2.0,
                        3,
                        4.0,
                        5,
                        6.0,
                        7,
                        8.0,
                        9,
                        10.0,
                        11,
                        12.0,
                    ]
                )
                axs[0][1].set_yticks(
                    [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
                )

    txt_e10 = f"#recs:\n{easiest_hardest[10]['easiest'][0][1]:,}".replace(
        ",", " "
    )
    txt_h10 = f"#recs:\n{easiest_hardest[10]['hardest'][1][1]:,}".replace(
        ",", " "
    )
    txt_e25 = f"#recs:\n{easiest_hardest[25]['easiest'][0][1]:,}".replace(
        ",", " "
    )
    txt_h25 = f"#recs:\n{easiest_hardest[25]['hardest'][1][1]:,}".replace(
        ",", " "
    )

    # Add text with number of recursions and annotate first amino acid.
    bbox_style = {
        "facecolor": "grey",
        "edgecolor": "black",
        "alpha": 0.3,
    }

    axs[0][0].text(-0.29, 6.25, txt_e10, ha="center", va="center").set_bbox(
        bbox_style
    )
    axs[0][0].text(0.06, -0.50, "1", fontsize=8, fontweight="demibold")
    axs[0][1].text(8.1, 2.5, txt_h10, ha="center", va="center").set_bbox(
        bbox_style
    )
    axs[0][1].text(0.14, -0.36, "1", fontsize=8, fontweight="demibold")
    axs[1][0].text(-0.29, 2.52, txt_e25, ha="center", va="center").set_bbox(
        bbox_style
    )
    axs[1][0].text(0.06, -0.21, "1", fontsize=8, fontweight="demibold")
    axs[1][1].text(9.5, 1.5, txt_h25, ha="center", va="center").set_bbox(
        bbox_style
    )
    axs[1][1].text(0.12, -0.22, "1", fontsize=8, fontweight="demibold")

    # Style axis.
    for ax in fig.get_axes():
        ax.axis("on")
        ax.margins(0.12)
        ax.label_outer()
        ax.tick_params(
            left=False,
            right=False,
            bottom=False,
            top=False,
            labelleft=False,
            labelbottom=False,
        )
        ax.grid(
            which="both",
            axis="both",
            color="grey",
            linestyle="--",
            linewidth=0.5,
            alpha=0.3,
            zorder=5,
        )

    # Fix legend to be roughly in the middle of plot.
    handles, labels = ax1.get_legend_handles_labels()
    score_patch = Line2D(
        [],
        [],
        color="indianred",
        linestyle=":",
        alpha=0.9,
        label="Bond",
        lw=2,
    )
    handles.append(score_patch)
    labels.append(score_patch.get_label())
    ax1.legend(handles=handles, labels=labels, loc="center", prop={"size": 12})

    # Fix labels.
    fig.supxlabel("Protein length", fontsize=14)
    ax3.set_xlabel(f"{lengths[0]}", fontsize=12)
    ax4.set_xlabel(f"{lengths[1]}", fontsize=12)
    ax1.set_ylabel("Easiest", labelpad=10, fontsize=12)
    ax3.set_ylabel("Hardest", labelpad=10, fontsize=12)
    plt.subplots_adjust(bottom=0.08, top=0.93, left=0.07, right=0.93)

    # Save if specified.
    if save:
        plt.savefig(
            f"figures/{job}/conf_{lengths}_n{nprot}.pdf",
            format="pdf",
            dpi=600,
        )
        print(f"\n~ Saved figure in figures/{job}/")

    # Show if specified.
    if show:
        plt.show()


def _rescale(coordinate_list):
    """
    Rescale the coordinates and return the new coordinates and max
    x,y-dimensions.
    """
    max_x = coordinate_list[0].x
    max_y = coordinate_list[0].y
    min_x = coordinate_list[0].x
    min_y = coordinate_list[0].y

    for c in coordinate_list:
        if c.x > max_x:
            max_x = c.x
        if c.x < min_x:
            min_x = c.x
        if c.y > max_y:
            max_y = c.y
        if c.y < min_y:
            min_y = c.y

    for i in range(len(coordinate_list)):
        coordinate_list[i] = (
            coordinate_list[i].x - min_x,
            coordinate_list[i].y - min_y,
        )

    return coordinate_list, max_x, max_y


def _create_graph(coordinate_list, dictionary, max_x, max_y, bond_idxs):
    M = max_x + 20
    N = max_y + 20
    size = max(N, M)

    G = nx.grid_2d_graph(N, M)

    pos = {(i, j): np.array([i, j]) for i in range(N) for j in range(M)}
    elist = []

    # Create list of edges.
    for i in range(len(coordinate_list)):
        if i < len(coordinate_list) - 1:
            elist.append([coordinate_list[i], coordinate_list[i + 1]])

    # Create list of bond edges.
    blist = []
    for b_i, b_j in bond_idxs:
        coord = [coordinate_list[b_i], coordinate_list[b_j]]
        coord_flip = [coordinate_list[b_j], coordinate_list[b_i]]
        if (coord not in elist and coord_flip not in elist) and (
            coord not in blist and coord_flip not in blist
        ):
            blist.append(coord)

    # Remove weird nodes?
    for i in range(N):
        for j in range(M):
            if (i, j) not in coordinate_list:
                G.remove_node((i, j))

    # Color nodes.
    node_color = []
    for node in pos:
        for key in dictionary:
            if node == key:
                if dictionary[key] == "H":
                    node_color.append("k")
                else:
                    node_color.append("white")

    # Draw the network.
    nodes = nx.draw_networkx_nodes(
        G, pos, node_color=node_color, node_size=300 / (0.05 * size)
    )

    nodes.set_edgecolor("k")
    nx.draw_networkx_edges(G, pos, elist, edge_color="k", width=1.0)
    nx.draw_networkx_edges(
        G, pos, blist, edge_color="r", style="--", width=1.0
    )

    plt.gca().invert_yaxis()
    plt.box(False)
    plt.tight_layout()

    plt.savefig("figures/front_imgs/front_img_2.pdf", dpi=600, format="pdf")
    plt.show()


def plot_front_figure():
    """
    Create the sub-figure used as front figure in paper. Inspiration taken from
    a colleague's plotting code:
        https://github.com/yingying-m/protein-folding.
    """
    sequence = "HPHHPPPHHHPPPHHHHPH"
    moves = [
        -1,
        -2,
        -1,
        2,
        -1,
        -2,
        -2,
        1,
        -2,
        -2,
        1,
        1,
        2,
        -1,
        2,
        1,
        1,
        -2,
    ]
    posit = [0, 0]
    coord_list = [Point(x=0, y=0)]

    for m in moves:
        posit[abs(m) - 1] += m // abs(m)
        coord_list.append(Point(posit[0], posit[1]))

    # Fetch bonds formed between amino acids.
    p = Protein(sequence, model="HP")
    p.set_hash(moves)
    bond_idxs = p.get_bonds()

    # Get new dimensions from network.
    coordinate_list, max_x, max_y = _rescale(coord_list)

    seq_array = [char for char in sequence]
    dictionary = {
        coordinate_list[i]: seq_array[i] for i in range(len(coordinate_list))
    }

    # Create and show conformation using NetworkX.
    _create_graph(coordinate_list, dictionary, max_x, max_y, bond_idxs)
