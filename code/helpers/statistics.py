#!/usr/bin/env python3
"""
File:           statistics.py
Description:    This file contains multiple functions that can be used for
                computing statistics on acquired data.
"""

import os
import ast
import math
from statistics import mean
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
from scipy.optimize import curve_fit
from prospr import Protein

from .functions import probit_func, linear_func


def comp_dist_gof(job, lengths, column="placed"):
    """
    Perform a 2-sample Kolmogoriv-Smirnov test and compute MSE for the Probit
    fits for specific lengths of a jobs results.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
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

    # Compute MSE and KS-test results for each length.
    for source_df, length in zip(
        source_dfs,
        lengths,
    ):
        # Print statistics on the loaded data.
        print(f"Values stats for length: {length}")
        print(f"\tMin:             {source_df[column].min()}")
        print(f"\tMax:             {source_df[column].max()}")
        print(f"\tMean:            {source_df[column].mean()}")
        print(f"\tStd:             {source_df[column].std()}")

        # Plot the regular histogram.
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
        print(f"Optimal params:   {popt}")

        # Compute Kolmogorov-Smirnov goodness-of-fit.
        probit_ys = [probit_func(x, *popt) for x in unit_xs]

        mse = ((log_ys - probit_ys) ** 2).mean()
        print(f"MSE:              {mse}")

        ks_ret = ks_2samp(log_ys, probit_ys)
        print(ks_ret, end="\n\n\n")


def comp_params_gof(job, lengths, column="placed"):
    """
    Perform a 2-sample Kolmogoriv-Smirnov test and compute MSE for the Probit
    fits for specific lengths of a jobs results.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    """
    root = f"jobs/{job}/results"

    # Load source file into a DataFrame.
    lengths.sort(reverse=False)
    filenames = [
        f
        for length in lengths
        for f in os.listdir(root)
        if f.startswith("HP") and int(f.split("_")[1]) == length
    ]
    source_dfs = [pd.read_csv(f"{root}/{f}") for f in filenames]
    popts = []

    # Fit the data distributions for each length, storing the popts.
    for source_df, length in zip(
        source_dfs,
        lengths,
    ):
        # Plot the regular histogram.
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
        popts.append(popt)
        print(f"Length {length}:   {popt}")
    print("")

    # Fit the fitted alpha parameter with a linear function and report MSE
    # against the popts.
    alpha_ys = [popt[0] for popt in popts]
    p0 = [alpha_ys[1] / alpha_ys[0], alpha_ys[1] - alpha_ys[0]]
    popt = curve_fit(linear_func, lengths, alpha_ys, p0=p0)[0]
    lin_ys = [linear_func(x, *popt) for x in lengths]
    mse = ((np.array(alpha_ys) - np.array(lin_ys)) ** 2).mean()
    ks_ret = ks_2samp(alpha_ys, lin_ys)
    print(f"Alpha popt:  {popt}")
    print(f"MSE:         {mse}")
    print(f"KS:          {ks_ret}\n")

    # Fit the fitted beta parameter with a linear function and report MSE
    # against the popts.
    beta_ys = [popt[1] for popt in popts]
    p0 = [beta_ys[1] / beta_ys[0], beta_ys[1] - beta_ys[0]]
    popt = curve_fit(linear_func, lengths, beta_ys, p0=p0)[0]
    lin_ys = [linear_func(x, *popt) for x in lengths]
    mse = ((np.array(beta_ys) - np.array(lin_ys)) ** 2).mean()
    ks_ret = ks_2samp(beta_ys, lin_ys)
    print(f"Beta popt:   {popt}")
    print(f"MSE:         {mse}")
    print(f"KS:          {ks_ret}\n")


def get_extrema(job, dataset, lengths, column="placed", n=1, add_placed=False):
    """
    Get the extrema proteins for specific lengths of a jobs results.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param func     dataset:    Function call to loading the data.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param int      n:          Amount easiest or hardest to take per length.
    :param bool     add_placed: True if number of placed amino acids is needed.
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
    extrema = {}

    # Plot for each length.
    for source_df, length in zip(source_dfs, lengths):
        proteins = dataset(length)
        source_df = source_df.sort_values(column)
        easiest, hardest = source_df.iloc[:n], source_df.iloc[-n:]
        easiest_p = []
        hardest_p = []

        # Create lists of Protein objects ordered from easy to hard.
        for easy, hard in zip(
            easiest.itertuples(index=False), hardest.itertuples(index=False)
        ):
            # Create Protein objects for easiest and hardest instances.
            easy_p = Protein(
                proteins.iloc[easy.protein_id]["sequence"], model="HP"
            )
            easy_p.set_hash(ast.literal_eval(easy.hash))
            hard_p = Protein(
                proteins.iloc[hard.protein_id]["sequence"], model="HP"
            )
            hard_p.set_hash(ast.literal_eval(hard.hash))

            if add_placed:
                easiest_p.append([easy_p, easy.placed])
                hardest_p.append([hard_p, hard.placed])
            else:
                easiest_p.append([easy_p])
                hardest_p.append([hard_p])

        # Store Proteins.
        extrema[length] = {"easiest": easiest_p, "hardest": hardest_p}

    return extrema


def print_n_extrema(job, dataset, lengths, column="placed", n=5, side="easy"):
    """
    Get the extrema proteins for specific lengths of a jobs results.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param func     dataset:    Function call to loading the data.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param int      n:          Amount easiest or hardest to take per length.
    :param str      side:       Either "easy" or "hard" to pick end of results.
    """
    root = f"jobs/{job}/results"

    # Load source file into a DataFrame.
    lengths.sort(reverse=False)
    filenames = [
        f
        for length in lengths
        for f in os.listdir(root)
        if f.startswith("HP") and int(f.split("_")[1]) == length
    ]
    source_dfs = [pd.read_csv(f"{root}/{f}") for f in filenames]

    if side == "easy":
        print(f"Easiest {n} for lengths {lengths}.")
    else:
        print(f"Hardest {n} for lengths {lengths}.")

    # Plot for each length.
    for source_df, length in zip(source_dfs, lengths):
        proteins = dataset(length)
        source_df = source_df.sort_values(column)

        if side == "easy":
            subset = source_df.iloc[:n]
        else:
            subset = source_df.iloc[-n:]

        print(f"\nLength {length}:")
        print("=========")

        for instance in subset.itertuples(index=False):
            print(proteins.iloc[instance.protein_id]["sequence"])


def comp_stats_extrema(job, dataset, lengths, column="placed", e_n=5, h_n=5):
    """
    Get the extrema proteins for specific lengths of a jobs results.
    :param str      job:        Folder name of the job which results will be
                                used.
    :param func     dataset:    Function call to loading the data.
    :param [int]    lengths:    Lengths to use.
    :param str      column:     The column from the results to use as data.
    :param int      n:          Amount easiest or hardest to take per length.
    :param str      side:       Either "easy" or "hard" to pick end of results.
    """
    root = f"jobs/{job}/results"

    # Load source file into a DataFrame.
    lengths.sort(reverse=False)
    filenames = [
        f
        for length in lengths
        for f in os.listdir(root)
        if f.startswith("HP") and int(f.split("_")[1]) == length
    ]
    source_dfs = [pd.read_csv(f"{root}/{f}") for f in filenames]

    for source_df, length in zip(source_dfs, lengths):
        proteins = dataset(length)
        source_df = source_df.sort_values(column)
        easy_subset = source_df.iloc[:e_n]
        hard_subset = source_df.iloc[-h_n:]

        print(f"\nLength {length}:")
        print("=========")

        # Compute H-ratio in first half of sequence.
        print("Avg P-ratio second half:")
        p_ratio_end = []
        for instance in easy_subset.itertuples(index=False):
            seq = proteins.iloc[instance.protein_id]["sequence"]
            p_ratio_end.append(seq[length // 2 :].count("P") / (length // 2))
        print(f"\tEasiest:  {mean(p_ratio_end)}")

        p_ratio_end = []
        for instance in hard_subset.itertuples(index=False):
            seq = proteins.iloc[instance.protein_id]["sequence"]
            p_ratio_end.append(
                seq[len(seq) // 2 :].count("P") / len(seq[len(seq) // 2 :])
            )
        print(f"\tHardest:  {mean(p_ratio_end)}")

        # Compute what fraction of bonds are in the first sequence half.
        print("Avg fraction of bonds in first half:")
        bond_frac_begin = []
        for instance in easy_subset.itertuples(index=False):
            p = Protein(
                proteins.iloc[instance.protein_id]["sequence"], model="HP"
            )
            p.set_hash(ast.literal_eval(instance.hash))
            bonds = p.get_bonds()

            if len(bonds) != 0:
                begin_bonds = [
                    b
                    for b in bonds
                    if b[0] <= length // 2 and b[1] <= length // 2
                ]
                bond_frac_begin.append((len(begin_bonds)) / len(bonds))
        print(f"\tEasiest:  {mean(bond_frac_begin)}")

        bond_frac_begin = []
        for instance in hard_subset.itertuples(index=False):
            p = Protein(
                proteins.iloc[instance.protein_id]["sequence"], model="HP"
            )
            p.set_hash(ast.literal_eval(instance.hash))
            bonds = p.get_bonds()

            if len(bonds) != 0:
                begin_bonds = [
                    b
                    for b in bonds
                    if b[0] <= length // 2 and b[1] <= length // 2
                ]
                bond_frac_begin.append((len(begin_bonds)) / len(bonds))
        print(f"\tHardest:  {mean(bond_frac_begin)}")

        # Compute what fraction of bonds are formed by placements in the second
        # half.
        print("Avg fraction of bonds formed by placement last quarter:")
        bond_frac_end = []
        for instance in easy_subset.itertuples(index=False):
            p = Protein(
                proteins.iloc[instance.protein_id]["sequence"], model="HP"
            )
            p.set_hash(ast.literal_eval(instance.hash))
            bonds = p.get_bonds()

            if len(bonds) != 0:
                begin_bonds = [
                    b
                    for b in bonds
                    if b[0] > 3 * length // 4 or b[1] > 3 * length // 4
                ]
                bond_frac_end.append((len(begin_bonds)) / len(bonds))
        print(f"\tEasiest:  {mean(bond_frac_end)}")

        bond_frac_end = []
        for instance in hard_subset.itertuples(index=False):
            p = Protein(
                proteins.iloc[instance.protein_id]["sequence"], model="HP"
            )
            p.set_hash(ast.literal_eval(instance.hash))
            bonds = p.get_bonds()

            if len(bonds) != 0:
                begin_bonds = [
                    b
                    for b in bonds
                    if b[0] > 3 * length // 4 or b[1] > 3 * length // 4
                ]
                bond_frac_end.append((len(begin_bonds)) / len(bonds))
        print(f"\tHardest:  {mean(bond_frac_end)}")

        # Compute H-ratio in first half of sequence.
        print("Avg #recursions:")
        e_nrecs = [i.placed for i in easy_subset.itertuples(index=False)]
        print(f"\tEasiest:  min:  {min(e_nrecs)}")
        print(f"\tEasiest:  mean: {mean(e_nrecs)}")
        print(f"\tEasiest:  max:  {max(e_nrecs)}")

        m_nrecs = [
            i.placed for i in source_df.iloc[e_n:-h_n].itertuples(index=False)
        ]
        print(f"\tMiddle:   min:  {min(m_nrecs)}")
        print(f"\tMiddle:   mean: {mean(m_nrecs)}")
        print(f"\tMiddle:   max:  {max(m_nrecs)}")

        h_nrecs = [i.placed for i in hard_subset.itertuples(index=False)]
        print(f"\tHardest:  min:  {min(h_nrecs)}")
        print(f"\tHardest:  mean: {mean(h_nrecs)}")
        print(f"\tHardest:  max:  {max(h_nrecs)}")

        print(f"\tEasiest:  frac: {mean(e_nrecs) / mean(m_nrecs)}")
        print(f"\tMiddle:   frac: {mean(m_nrecs) / mean(m_nrecs)}")
        print(f"\tHardest:  frac: {mean(h_nrecs) / mean(m_nrecs)}")
