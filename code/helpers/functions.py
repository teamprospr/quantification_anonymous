#!/usr/bin/env python3
"""
File:           functions.py
Description:    This file contains multiple functions that can be used for
                fitting.
"""

import numpy as np
from scipy.special import erfinv


def probit_func(x, alpha, beta):
    """
    Probit function with scaling parameters to fit data to.
    Note:   Log-scale y-values when fitting this function.
    :param  int      x:      Rank of the point in the figure (x-axis).
    :param  float    alpha:  Scaling factor (sqrt(2) in original probit)
    :param  float    beta:   Y-axis shift.
    """
    return alpha * erfinv(2 * x - 1) + beta


def probit_func_roofed(x, alpha, beta, gamma):
    """
    Probit function with scaling parameters to fit data to.
    Note:   Log-scale y-values when fitting this function.
    :param  int      x:      Rank of the point in the figure (x-axis).
    :param  float    alpha:  Scaling factor (sqrt(2) in original probit)
    :param  float    beta:   Y-axis shift.
    :param  float    gamma:  Value of the theoretical plateau.
    """
    return np.minimum(alpha * erfinv(2 * x - 1) + beta, gamma)


def logit_func(x, alpha, beta):
    """
    Logit function with scaling parameters to fit data to.
    Note:   Log-scale y-values when fitting this function.
    :param  int      x:      Rank of the point in the figure (x-axis).
    :param  float    alpha:  Scaling factor (sqrt(2) in original probit)
    :param  float    beta:   Y-axis shift.
    """
    return alpha * np.log(x / (1 - x)) + beta


def linear_func(x, alpha, beta):
    """
    Linear function with scaling parameters to fit data to.
    :param  int      x:      Rank of the point in the figure (x-axis).
    :param  float    alpha:  Slope factor.
    :param  float    beta:   Y-axis shift.
    """
    return alpha * x + beta
