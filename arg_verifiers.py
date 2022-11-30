#!/usr/bin/env python3
"""
File:           arg_verifiers.py
Description:    This file contains verifiers for the argument parser in main.
"""

import os
import argparse
from code.helpers.datasets import load_proteins250, load_proteins1000


def job_type():
    """
    Returns function handle which acts as a job checker for argparse.
    Argument is valid if the job exists, and there is an entry script with the
    same name.
    """

    def job_checker(jobname):
        """
        New Type function for argparse, verifies if the job directory exists
        and if there is an entry script in the directory.
        """
        # Check if directory exists.
        if not os.path.isdir(f"jobs/{jobname}"):
            raise argparse.ArgumentTypeError(
                "Provided job directory must " "exist."
            )

        if not os.path.isfile(f"jobs/{jobname}/{jobname}.sh"):
            raise argparse.ArgumentTypeError(
                "Provided job directory must have"
                "a entry script with the same "
                "name."
            )

        return jobname

    return job_checker


def dataset_type():
    """
    Returns function handle which acts as a dataset checker for argparse.
    Argument is valid if the dataset name can be parsed into a dataset loader
    function from prospr.
    """

    def dataset_checker(dataset_name):
        """
        New Type function for argparse, verifies if the dataset exists within
        prospr and returns function pointer to the loader.
        """
        if dataset_name == "proteins250":
            return load_proteins250
        elif dataset_name == "proteins1000":
            return load_proteins1000

        # If nothing could be parsed, throw error.
        raise argparse.ArgumentTypeError(
            "Provided dataset could not be parsed into a prospr dataset."
        )

    return dataset_checker
