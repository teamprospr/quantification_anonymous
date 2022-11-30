#!/usr/bin/env python3
"""
File:           datasets.py
Description:    This file contains functions for loading datasets.
"""
import pandas as pd


def load_proteins250(length=10):
    """Returns a proteins250 dataset as a pandas dataframe."""
    return pd.read_csv(f"data/proteins250/HP{length}.csv")


def load_proteins1000(length=10):
    """Returns a proteins1000 dataset as a pandas dataframe."""
    return pd.read_csv(f"data/proteins1000/HP{length}.csv")
