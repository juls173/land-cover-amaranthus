import pandas as pd
import numpy as np
import math

def calculate_shannon_index(proportions: np.ArrayLike) -> float:
    proportions = np.array(proportions) / 100

    sum = 0
    for proportion in proportions:
        if proportion > 0:
            sum += proportion * math.log(proportion)
    return -1 * sum

def calculate_simpson_index(proportions: np.ArrayLike) -> float:
    proportions = np.array(proportions) / 100
    simpson_index = 1 - (proportions ** 2).sum()

    return simpson_index

def proportion_developed(lu_classes: pd.DataFrame) -> float:
    pass

# only one agricultural class, "Cultivated Crops", excluding "Pasture/Hay"
    



