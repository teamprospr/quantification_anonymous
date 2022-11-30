#!/usr/bin/env python3
"""
File:           run_visualizer.py
Description:    This file parses user input and runs the requested visualizer
                functions.
"""

import argparse

from arg_verifiers import job_type, dataset_type
from code.helpers.visualizer import (
    plot_single,
    plot_group,
    plot_fit_single,
    plot_fit_group,
    plot_params,
    plot_dataset_boxplot,
    plot_easiest_hardest,
    plot_conformation,
    plot_front_figure,
)
from prospr import load_vanEck1000


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
        "-ns",
        "--no_show",
        default=True,
        action="store_false",
        help="Set flag if you want plt.show() to be disabled.",
    )
    parser.add_argument(
        "-s",
        "--save",
        default=False,
        action="store_true",
        help="Set flag if generated figures should be saved.",
    )
    parser.add_argument(
        "-ds",
        "--dataset",
        default=load_vanEck1000,
        type=dataset_type(),
        help="Dataset from prospr to use as input.",
    )
    parser.add_argument(
        "-c",
        "--column",
        dest="column",
        type=str,
        default="placed",
        help="The column from the results to use as data.",
    )

    # Parse args and return.
    return parser.parse_args()


if __name__ == "__main__":
    # Compose list of commands and parse user input.
    commands = [
        "plot_single",
        "plot_group",
        "fit_single",
        "fit_group",
        "params",
        "dataset",
        "easiest_hardest",
        "plot_conf",
        "front_figure",
    ]
    args = parse_arguments(commands)

    # Parse the program to run, which can be of the list above.
    # Plot a single length of a job's results.
    if args.program == "plot_single":
        plot_single(
            args.job,
            args.lengths[0],
            column=args.column,
            show=args.no_show,
            save=args.save,
        )
    # Plot multiple lengths of a job's results.
    elif args.program == "plot_group":
        plot_group(
            args.job,
            args.lengths,
            column=args.column,
            show=args.no_show,
            save=args.save,
        )
    # Fit single length with a function.
    elif args.program == "fit_single":
        plot_fit_single(
            args.job, args.lengths[0], show=args.no_show, save=args.save
        )
    # Fit multiple lengths with a function.
    elif args.program == "fit_group":
        plot_fit_group(
            args.job, args.lengths, show=args.no_show, save=args.save
        )
    # Fit multiple lengths with a function and plot the flow of parameters.
    elif args.program == "params":
        plot_params(args.job, args.lengths, show=args.no_show, save=args.save)
    # Create boxplot from a dataset with HP-ratio per length.
    elif args.program == "dataset":
        plot_dataset_boxplot(
            args.dataset, args.lengths, show=args.no_show, save=args.save
        )
    # Plot the conformations of the easiest and hardest instances for the given
    # lengths.
    elif args.program == "easiest_hardest":
        plot_easiest_hardest(
            args.job,
            args.dataset,
            args.lengths,
            show=args.no_show,
            save=args.save,
        )
    # Plot the conformations of the easiest and hardest instances for the given
    # lengths.
    elif args.program == "plot_conf":
        plot_conformation(
            args.job,
            args.dataset,
            args.lengths,
            show=args.no_show,
            save=args.save,
        )
    # Plot the conformations of the easiest and hardest instances for the given
    # lengths.
    elif args.program == "front_figure":
        plot_front_figure()
    # Give error when nothing could be parsed from the input.
    else:
        print("\nInput program could not be parsed..")
