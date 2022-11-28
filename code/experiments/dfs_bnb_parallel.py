#!/usr/bin/env python3
"""
File:           bnb_dfs_parallel.py
Description:    This file uses the bnb_dfs algorithm in parallel to compute
                statistics for the specified groups of the specified dataset.
                The computed statistics are:
                    - Time it took to compute exact minimum of the score.
                    - Minimum found score.
                    - Number of placements for an amino acid.
                    - Number of full conformations checked.
"""

import multiprocessing as mp
from time import time
import pandas as pd
import csv
import os
from pathlib import Path
from prospr import Protein, depth_first_bnb


def load_proteins(results_path, dataset, lengths, dim):
    """
    Detect complete or partially complete results, and only load the proteins
    that need to be finished.
    :param results_path:    Path to the results folder.
    :param dataset:         Function to load dataset data.
    :param lengths:         Lengths to process.
    :param dim:             Dimension in which the Proteins will be folded.
    :return:                List of Protein objects, grouped by length.
    """
    proteins = {}
    finished_ids = []

    # Print debug header for loading dataset if debug is on.
    print("Loading Dataset Info:")

    # Load protein strings into Protein objects from groups that are not
    # finished.
    for l in lengths:
        length_filepath = f"{results_path}/HP_{l}_bnb_dfs"

        if os.path.exists(f"{length_filepath}.csv"):
            # Skip group if all proteins are already finished.
            print(f"\tResults detected, skipping length {l}..")
            continue
        elif os.path.exists(f"{length_filepath}_TMP.csv"):
            # Only load not finished proteins if TMP file is present.
            print(f"\tPartial results detected for length {l}..")

            # Fetch IDs of the finished proteins.
            with open(f"{length_filepath}_TMP.csv", "r") as fp:
                df = pd.read_csv(fp)
                finished_ids = df["protein_id"].values
        else:
            # Load all proteins if group has never been computed.
            print(
                f"\tNo results detected for length {l}, loading all "
                f"proteins.."
            )

        # Only load the protein with length l which are not finished.
        df = dataset(l)
        proteins[l] = [
            (row.id, Protein(row.sequence, dim, model="HP"))
            for _, row in df.iterrows()
            if row.id not in finished_ids
        ]

    # Create whitespace between last debug row and upcoming debug info.
    print()

    return proteins


def print_debug_head(largest_group, proc_count, dataset, results_path):
    """
    Print debug info header and start of the progress debug info.
    :param int  largest_group:      Number of proteins in the largest group.
    :param int  proc_count:         Number of used processes.
    :param func  dataset:           Path to the used dataset.
    :param str  results_path:       Path to where the results will be stored.
    """
    print("Experiment Info:")
    print("\tLargest group:             ", largest_group)
    print("\tNumber of processes:       ", proc_count)
    print("\tData path:                 ", dataset.__name__)
    print("\tResults path:              ", results_path, end="\n\n")
    print("Progress:")
    print("=========")


def listener(group_filepath, size, q):
    """
    Listener that fetches all results from queue and then writes them in
    order of protein id.
    :param str      group_filepath:     Path where the results will be written
                                        to.
    :param int      size:               Number of proteins that produce
                                        results.
    :param Queue    q:                  Manager.Queue() object for receiving
                                        results.
    """
    results = []
    header = [
        "protein_id",
        "algorithm",
        "time",
        "score",
        "checked",
        "placed",
        "hash",
    ]
    add_header = not os.path.exists(f"{group_filepath}_TMP.csv")

    with open(f"{group_filepath}_TMP.csv", "a") as fp:
        csv_writer = csv.writer(fp)

        if add_header:
            csv_writer.writerow(header)

        # Fetch all results and store in the group's temporary file.
        while len(results) != size:
            message = q.get()
            results.append(message)
            csv_writer.writerow(message)

    # When all results have been gathered, save them ordered by protein id.
    # Load all proteins from TMP file, since previous results also need to be
    # stored.
    with open(f"{group_filepath}_TMP.csv", "r") as fp:
        fp.readline()
        reader = csv.reader(fp)
        results = [(int(row[0]), *row[1:]) for row in reader]
        results.sort()

    # Store all results
    with open(f"{group_filepath}.csv", "w") as fp:
        csv_writer = csv.writer(fp)
        csv_writer.writerow(header)

        for row in results:
            csv_writer.writerow(row)

    # Delete temporary file and end listener process.
    os.remove(f"{group_filepath}_TMP.csv")


def worker(protein, p_id, q, debug=False):
    """
    Worker that folds one protein and stores statistics in queue.
    :param Protein  protein:    Protein object to process.
    :param int      p_id:       ID of the protein to process.
    :param Queue    q:          Manager.Queue() object for sending results.
    :param bool     debug:      Flag set if debug info needs to be printed.
    """
    start_time = time()
    depth_first_bnb(protein)
    end_time = time()

    q.put(
        [
            p_id,
            "dfs_bnb",
            end_time - start_time,
            protein.score,
            protein.solutions_checked,
            protein.aminos_placed,
            protein.hash_fold(),
        ]
    )

    if debug:
        print(
            f"\t       Worker for {p_id} done in {end_time - start_time} \t"
            f"seconds"
        )


def dfs_bnb_parallel(
    job, dataset, lengths, max_procs, dim=2, exec_lisa=False, debug=False
):
    """
    Test to measure speed of the SurfSara supercomputer.
    :param str      job:        Folder of the job that is running.
    :param func     dataset:    Function call to loading the data.
    :param int[]    lengths:    What lengths to use, e.g. [10, 20].
    :param int      max_procs:  Maximum number of processes to create.
    :param bool     exec_lisa:  Flag set if run environment is the Lisa
                                cluster.
    :param bool     debug:      Flag set if debug info needs to be printed.

    NOTE:   This code does not work on a Windows machine! Create a __main__
            function and check output for Windows multiprocessing.
    """
    # Set paths according to run environment.
    if exec_lisa:
        results_path = f"{os.environ['LISA_DIR']}/{job}/results"
    else:
        results_path = f"jobs/{job}/results"

    # Create results directory, if not exist.
    Path(results_path).mkdir(parents=True, exist_ok=True)

    # Load proteins of all groups in a dictionary.
    proteins = load_proteins(results_path, dataset, lengths, dim)

    # Exit if there are no proteins that need to be processed.
    if not proteins.values():
        print_debug_head(1, 1, dataset, results_path)
        print("\tNo protein requires processing, shutting down..")
        return

    # Set number of processes equal the maximum number of proteins in a group,
    # or the given maximum.
    proc_count = max(map(len, proteins.values()))

    if proc_count > max_procs:
        proc_count = max_procs

    # Print debug info.
    if debug:
        print_debug_head(
            max(map(len, proteins.values())), proc_count, dataset, results_path
        )

    # Setup Manager Queue and the Process Pool.
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(proc_count)

    # Perform parallel execution for every group and store results
    # sequentially.
    for i, (group_id, group_proteins) in enumerate(proteins.items()):
        group_filepath = f"{results_path}/HP_{group_id}_bnb_dfs"

        # Initialise listener for writing statistics for current group.
        workers = [
            pool.apply_async(
                listener, (group_filepath, len(group_proteins), q)
            )
        ]

        if debug:
            print(f"\t[{i + 1}/{len(lengths)}]  Listener created")

        # Create a task per protein for the workers.
        for j, p in group_proteins:
            workers.append(pool.apply_async(worker, (p, j, q, debug)))

            if debug:
                print(f"\t[{i + 1}/{len(lengths)}]  Worker for {j} created")

        # Wait for all workers to finish the current group.
        [w.get() for w in workers]

        if debug:
            print(f"\t[{i + 1}/{len(lengths)}]  All workers done\n")

    if debug:
        print("\tAll work done, closing pool..")

    pool.close()
    pool.join()
