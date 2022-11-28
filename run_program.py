#!/usr/bin/env python3
"""
File:           run_program.py
Description:    This file parses user input and runs the requested script for
                folding proteins.
"""

# General imports.
import os
import argparse
from time import time
from arg_verifiers import job_type, dataset_type

# Import Prospr protein related functions.
from prospr import Protein, depth_first, depth_first_bnb
from prospr.datasets import load_vanEck250
from prospr.visualize import plot_protein

# Import experiments.
from code.experiments.dfs_bnb_serial import dfs_bnb_serial
from code.experiments.dfs_bnb_parallel import dfs_bnb_parallel


def parse_arguments(commands):
    """
    Create argument parser, parse given arguments and return results.
    :param  [str*]      commands:           List with possible commands to run.
    :return Namespace, ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description="Process input for execution."
    )
    parser.add_argument(
        "program",
        metavar="P",
        type=str,
        choices=commands,
        help="Program to execute. Options: %(choices)s",
    )
    parser.add_argument(
        "-j",
        "--job",
        type=job_type(),
        help="Job folder to use if required by selected program.",
    )
    parser.add_argument(
        "-d",
        "--debug",
        default=False,
        action="store_true",
        help="Use to print debug information during runtime",
    )
    parser.add_argument(
        "-dim",
        "--dimension",
        dest="dim",
        type=int,
        default=2,
        metavar="DIM",
        help="Dimension of the space for the protein to fold",
    )
    parser.add_argument(
        "-ds",
        "--dataset",
        default=load_vanEck250,
        type=dataset_type(),
        help="Dataset from prospr to use as input.",
    )
    parser.add_argument(
        "-le",
        "--lengths",
        nargs="+",
        dest="lengths",
        type=int,
        default=[10],
        help="Lengths of proteins to fold from dataset.",
    )
    parser.add_argument(
        "-v",
        "--visualise",
        default=False,
        action="store_true",
        help="Use to print visuals at end of execution",
    )
    parser.add_argument(
        "-l",
        "--lisa",
        default=False,
        action="store_true",
        help="Use if the execution is on the Lisa cluster",
    )

    # Parse args and return.
    return parser.parse_args()


def run_algorithm(args, algorithm):
    """
    Runs an algorithm with one protein and outputs results. Has the option
    to visualise the resulting fold.
    :param Namespace    args:           Parsed arguments.
    :param func         algorithm:      Function that runs the algorithm.
    """
    # Request needed information for creating the protein.
    print(
        "Requesting needed information.\nLeave blank for the default " "value."
    )
    print("==================================", end="\n\n")

    aminos = None

    while not aminos:
        aminos = input("Protein to use: ").strip()

    model = input("Model to use (default=('HP')): ").strip()
    model = model if model else "HP"

    # Create protein and run algorithm.
    protein = Protein(aminos, args.dim, model=model)
    start_time = time()
    algorithm(protein)
    end_time = time()

    # Print results.
    print("\nResults:\n========")
    print(f"\tTime:       {end_time - start_time}")
    print(f"\tScore:      {protein.score}")
    print(f"\tChecked:    {protein.solutions_checked}")
    print(f"\tHash:       {protein.hash_fold()}")

    # Output visualisation if specified as flag.
    if args.visualise:
        plot_protein(protein)

    return protein, protein.score


def run_dfs_bnb_serial(args):
    """
    Set path, run the time_lisa script and write duration.
    :param Namespace    args:   Results from the parsed arguments.
    """
    # Derive root path from run environment.
    if args.lisa:
        root = f"{os.environ['LISA_DIR']}/{args.job}"
    else:
        root = f"jobs/{args.job}"

    # Run test.
    start_time = time()
    dfs_bnb_serial(
        args.job,
        args.dataset,
        args.lengths,
        dim=args.dim,
        exec_lisa=args.lisa,
        debug=args.debug,
    )
    end_time = time()

    # Print results.
    print("\nResults:\n========")
    print(f"\tTime: {end_time - start_time}")

    # Write execution time to file if there is no timing of the groups yet.
    time_filepath = f"{root}/results/time_HP_{args.lengths}"

    if not os.path.exists(time_filepath):
        with open(time_filepath, "w") as fp:
            fp.write(str(end_time - start_time) + "\n")


def run_bnb_dfs_parallel(args):
    """
    Set path, run the time_lisa_parallel script and write duration.
    :param Namespace    args:   Results from the parsed arguments.
    """
    # Derive root path and max number of processes from run environment.
    if args.lisa:
        root = f"{os.environ['LISA_DIR']}/{args.job}"
        max_procs = 16
    else:
        root = f"jobs/{args.job}"
        max_procs = 4

    start_time = time()
    dfs_bnb_parallel(
        args.job,
        args.dataset,
        args.lengths,
        max_procs,
        dim=args.dim,
        exec_lisa=args.lisa,
        debug=args.debug,
    )
    end_time = time()

    # Print results.
    print("\nResults:\n========")
    print(f"\tTime: {end_time - start_time}")

    # Write execution time to file if there is no timing of the groups yet.
    time_filepath = f"{root}/results/time_HP_{args.groups}"

    if not os.path.exists(time_filepath):
        with open(time_filepath, "w") as fp:
            fp.write(str(end_time - start_time) + "\n")


if __name__ == "__main__":
    # Compose list of commands.
    algorithm_commands = ["run_dfs", "run_dfs_bnb"]

    # Setup list with commands that need a job specified.
    job_commands = ["dfs_bnb_serial", "dfs_bnb_parallel"]

    # Setup list with commands and parse arguments.
    command_lists = [
        algorithm_commands,
        job_commands,
    ]
    commands = [item for sublist in command_lists for item in sublist]
    args = parse_arguments(commands)

    # Print general debug info.
    print("\n\nDebug Info:")
    print("\tExecution on Lisa-Cluster:    ", args.lisa)
    print("\tDataset to use:               ", args.dataset)
    print("\tLengths to fold:              ", args.lengths)
    print("\tDimension of fold:            ", args.dim, end="\n\n")

    # Parse the program to run, which can be of the list above.
    # Parse an algorithm to run.
    if args.program == "run_dfs":
        run_algorithm(args, depth_first)
    elif args.program == "run_dfs_bnb":
        run_algorithm(args, depth_first_bnb)
    # Parse an experiment to run.
    elif args.program == "dfs_bnb_serial":
        run_dfs_bnb_serial(args)
    elif args.program == "dfs_bnb_parallel":
        run_dfs_bnb_serial(args)
    # Give error when nothing could be parsed from the input.
    else:
        print("\nInput program could not be parsed..")
