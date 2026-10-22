import numpy as np
def calculate_nse(observed, simulated):
    """
    Calculate the Nash–Sutcliffe model efficiency coefficient (NSE) with consideration of missing values.

    Parameters:
    observed (numpy array): Array of observed discharge values.
    simulated (numpy array): Array of simulated discharge values.

    Returns:
    float: NSE value.
    """
    observed = np.array(observed)
    simulated = np.array(simulated)

    # Filter out pairs where either observed or simulated is NaN
    valid_mask = ~np.isnan(observed) & ~np.isnan(simulated)
    observed = observed[valid_mask]
    simulated = simulated[valid_mask]

    # Mean of filtered observed data
    mean_observed = np.mean(observed)

    # Numerator: sum of squared differences between filtered observed and simulated values
    numerator = np.sum((observed - simulated) ** 2)

    # Denominator: sum of squared differences between filtered observed values and mean of filtered observed values
    denominator = np.sum((observed - mean_observed) ** 2)

    # NSE calculation
    nse = 1 - (numerator / denominator)

    return nse