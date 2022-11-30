#!/usr/bin/env python3
"""
File:           dfs_bnb_sequential.py
Description:    This file uses the depth_first_bnb algorithm from Prospr to
                compute statistics for the specified lengths of proteins from
                the specified dataset.
                The computed statistics are:
                    - Time it took to compute exact minimum of the score.
                    - Minimum found score.
                    - Number of placements for an amino acid.
                    - Number of full conformations checked.
"""

from time import time
import csv
import os
from prospr import Protein, depth_first_bnb
from pathlib import Path


def dfs_bnb_serial(job, dataset, lengths, dim=2, exec_lisa=False, debug=False):
    """
    Executes a serial session of depth_first_bnb searches over the dataset.
    :param str      job:        Folder of the job that is running.
    :param func     dataset:    Function call to loading the data.
    :param int[]    lengths:    What lengths to use, e.g. [10, 20].
    :param int      dim:        Dimension in which the Proteins will be folded.
    :param bool     exec_lisa:  Flag set if run environment is the Lisa
                                cluster.
    :param bool     debug:      Flag set if debug info needs to be printed.
    """
    # Set paths according to run environment.
    if exec_lisa:
        results_path = f"{os.environ['LISA_DIR']}/{job}/results"
    else:
        results_path = f"jobs/{job}/results"

    # Create results directory, if not exist.
    Path(results_path).mkdir(parents=True, exist_ok=True)

    # Load specified proteins into Protein objects.
    proteins = [[] for _ in range(len(lengths))]
    for i, l in enumerate(lengths):
        aminos = dataset(l)

        for _, row in aminos.iterrows():
            proteins[i].append(Protein(row.sequence, dim, model="HP"))

    # Print debug info.
    if debug:
        print("Experiment Info:")
        print("\tLargest group:             ", max(map(len, proteins)))
        print("\tDataset:                   ", dataset.__name__)
        print("\tResults path:              ", results_path, end="\n\n")
        print("Progress:")
        print("=========")

    # Perform timed checks and store results per length.
    for i, l in enumerate(lengths):
        if debug:
            print(f"\t[{i + 1}/{len(lengths)}]  Starting computations..")

        # Skip group if there are already results.
        result_filepath = f"{results_path}/HP_{l}_bnb_dfs.csv"

        if os.path.exists(result_filepath):
            if debug:
                print(f"\t\tResults detected, skipping group {l}")
            continue

        # Open CSV to store header and results.
        with open(result_filepath, "w") as fp:
            csv_writer = csv.writer(fp)
            csv_writer.writerow(
                [
                    "protein_id",
                    "algorithm",
                    "time",
                    "score",
                    "checked",
                    "placed",
                    "hash",
                ]
            )

            # Time execution of the bnb_dfs algorithm.
            for j, p in enumerate(proteins[i]):
                start_time = time()
                depth_first_bnb(p)
                end_time = time()

                csv_writer.writerow(
                    [
                        j,
                        "bnb_dfs",
                        end_time - start_time,
                        p.score,
                        p.solutions_checked,
                        p.aminos_placed,
                        p.hash_fold(),
                    ]
                )
                p.reset()

                if debug:
                    print(
                        f"\t\tProtein [{j + 1}/{len(proteins[i])}] done in "
                        f"{end_time - start_time} \tseconds"
                    )
