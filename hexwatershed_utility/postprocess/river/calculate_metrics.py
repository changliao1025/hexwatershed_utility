"""
Hydrological model performance evaluation metrics.
"""

from typing import Union
import numpy as np


def _filter_valid(observed: np.ndarray, simulated: np.ndarray):
    """Filter out NaN, Inf, and masked entries."""
    obs = np.asarray(observed, dtype=float)
    sim = np.asarray(simulated, dtype=float)
    valid = np.isfinite(obs) & np.isfinite(sim)
    return obs[valid], sim[valid]


def calculate_nse(observed: Union[list, np.ndarray], simulated: Union[list, np.ndarray]) -> float:
    """
    Calculate the Nash–Sutcliffe model efficiency coefficient (NSE).

    Parameters
    ----------
    observed : array_like
        Observed values.
    simulated : array_like
        Simulated values.

    Returns
    -------
    float
        NSE value (-inf to 1.0, where 1.0 represents a perfect match).
    """
    obs, sim = _filter_valid(observed, simulated)
    if len(obs) == 0:
        return np.nan

    mean_obs = np.mean(obs)
    denominator = np.sum((obs - mean_obs) ** 2)
    if denominator == 0:
        return np.nan

    numerator = np.sum((obs - sim) ** 2)
    return float(1.0 - (numerator / denominator))


def calculate_kge(observed: Union[list, np.ndarray], simulated: Union[list, np.ndarray]) -> float:
    """
    Calculate the Kling-Gupta Efficiency (KGE, Gupta et al. 2009).

    Parameters
    ----------
    observed : array_like
        Observed values.
    simulated : array_like
        Simulated values.

    Returns
    -------
    float
        KGE value (-inf to 1.0, where 1.0 represents a perfect match).
    """
    obs, sim = _filter_valid(observed, simulated)
    if len(obs) < 2:
        return np.nan

    std_obs = np.std(obs)
    std_sim = np.std(sim)
    mean_obs = np.mean(obs)
    mean_sim = np.mean(sim)

    if std_obs == 0 or mean_obs == 0:
        return np.nan

    r = np.corrcoef(obs, sim)[0, 1]
    alpha = std_sim / std_obs
    beta = mean_sim / mean_obs

    return float(1.0 - np.sqrt((r - 1.0) ** 2 + (alpha - 1.0) ** 2 + (beta - 1.0) ** 2))


def calculate_rmse(observed: Union[list, np.ndarray], simulated: Union[list, np.ndarray]) -> float:
    """
    Calculate Root Mean Squared Error (RMSE).

    Parameters
    ----------
    observed : array_like
        Observed values.
    simulated : array_like
        Simulated values.

    Returns
    -------
    float
        RMSE value.
    """
    obs, sim = _filter_valid(observed, simulated)
    if len(obs) == 0:
        return np.nan
    return float(np.sqrt(np.mean((obs - sim) ** 2)))


def calculate_pbias(observed: Union[list, np.ndarray], simulated: Union[list, np.ndarray]) -> float:
    """
    Calculate Percent Bias (PBIAS).

    Parameters
    ----------
    observed : array_like
        Observed values.
    simulated : array_like
        Simulated values.

    Returns
    -------
    float
        PBIAS percentage (0.0 is optimal).
    """
    obs, sim = _filter_valid(observed, simulated)
    if len(obs) == 0:
        return np.nan
    sum_obs = np.sum(obs)
    if sum_obs == 0:
        return np.nan
    return float(100.0 * np.sum(sim - obs) / sum_obs)
