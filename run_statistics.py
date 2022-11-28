#!/usr/bin/env python3
"""
File:           run_visualizer.py
Description:    This file parses user input and runs the requested visualizer
                functions.
"""

import argparse
from prospr import load_vanEck1000

from arg_verifiers import job_type, dataset_type
from code.helpers.statistics import (
    get_extrema,
    comp_dist_gof,
    comp_params_gof,
    print_n_extrema,
    comp_stats_extrema,
)


def parse_arguments(commands):
    """
    Create argument parser, parse given arguments and return results.
    :param  [str*]      commands:           List with possible commands to run.
    :return Namespace, ArgumentParser
    """
    parser = argparse.ArgumentParser(
        description="Process input for creating visualizations."
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
        "-le",
        "--lengths",
        nargs="+",
        dest="lengths",
        type=int,
        default=[10],
        help="Lengths of proteins to fold from dataset.",
    )
    parser.add_argument(
        "-c",
        "--column",
        dest="column",
        type=str,
        default="placed",
        help="The column from the results to use as data.",
    )
    parser.add_argument(
        "-ds",
        "--dataset",
        default=load_vanEck1000,
        type=dataset_type(),
        help="Dataset from prospr to use as input.",
    )

    # Parse args and return.
    return parser.parse_args()


if __name__ == "__main__":
    # Compose list of commands and parse user input.
    commands = [
        "get_extrema",
        "comp_dist_gof",
        "comp_params_gof",
        "print_n_easiest",
        "print_n_hardest",
        "comp_stats_extrema",
    ]
    args = parse_arguments(commands)

    # Parse the program to run, which can be of the list above.
    # Plot a single length of a job's placed amino acids.
    if args.program == "get_extrema":
        get_extrema(args.job, args.dataset, args.lengths, column=args.column)
    # Execute the 2 sample Kolmogorov-Smirnov test.
    elif args.program == "comp_dist_gof":
        comp_dist_gof(args.job, args.lengths, column=args.column)
    # Execute the 2 sample Kolmogorov-Smirnov test.
    elif args.program == "comp_params_gof":
        comp_params_gof(args.job, args.lengths, column=args.column)
    # Print the n easiest proteins.
    elif args.program == "print_n_easiest":
        print_n_extrema(
            args.job,
            args.dataset,
            args.lengths,
            column=args.column,
            n=30,
            side="easy",
        )
    # Print the n hardest proteins.
    elif args.program == "print_n_hardest":
        print_n_extrema(
            args.job,
            args.dataset,
            args.lengths,
            column=args.column,
            n=30,
            side="hard",
        )
    # Print statistics on the extrema.
    elif args.program == "comp_stats_extrema":
        comp_stats_extrema(
            args.job,
            args.dataset,
            args.lengths,
            column=args.column,
            e_n=50,
            h_n=50,
        )
    # Give error when nothing could be parsed from the input.
    else:
        print("\nInput program could not be parsed..")
