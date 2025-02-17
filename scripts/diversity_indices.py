import pandas as pd
import numpy as np
import math

def calculate_shannon_index(proportions: list[float] | pd.Series) -> float:
    """
    Calculates the Shannon Diversity Index, a measure of species or land cover 
    diversity based on the proportions of land cover classes.

    Args:
        proportions (list[float]): A list of proportions (percentages) representing 
                                   different land cover classes.

    Returns:
        float: The Shannon Diversity Index value.
    """

    proportions = np.array(proportions, dtype=float)
    proportions /= 100

    sum = 0
    for proportion in proportions:
        if proportion > 0:
            sum += proportion * math.log(proportion)
    return -1 * sum


def calculate_simpson_index(proportions: list[float] | pd.Series) -> float:
    """
    Calculates the Simpson Diversity Index, which measures the probability that 
    two randomly selected individuals belong to the same category.

    Args:
        proportions (list[float]): A list of proportions (percentages) representing 
                                   different land cover classes.

    Returns:
        float: The Simpson Diversity Index value (1 - D, where D is the dominance index).
    """

    proportions = np.array(proportions, dtype=float)
    proportions /= 100
    
    simpson_index = 1 - (proportions ** 2).sum()

    return simpson_index

def proportion_developed(lu_classes: pd.DataFrame) -> float:
    pass

# only one agricultural class, "Cultivated Crops", excluding "Pasture/Hay"
    



