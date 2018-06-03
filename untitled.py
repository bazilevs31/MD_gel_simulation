import numpy as np
import pandas as pd
def in1d(a, b, tol=1e-3):
    """looks for similarities in two arrays with a certain tolerance

    Args:
        a (array): Description
        b (array): Description
        tol (float, optional): Description

    Returns:
        TYPE: Description
    """
    a = np.unique(a)
    intervals = np.empty(2*a.size, float)
    intervals[::2] = a - tol
    intervals[1::2] = a + tol
    overlaps = intervals[:-1] >= intervals[1:]
    overlaps[1:] = overlaps[1:] | overlaps[:-1]
    keep = np.concatenate((~overlaps, [True]))
    intervals = intervals[keep]
    return np.searchsorted(intervals, b, side='right') & 1 == 1

file1 = 'Cond_q1_st1000_ef1500_k1_ds10_sc1844_nm100_l50_fl01_fr02.csv'
file2 = 'Cond_q1_st1000_ef1500_k1_ds10_sc1844_nm100_l50_fl01_fr02_new.csv'
d1 = pd.read_csv(file1)
d2 = pd.read_csv(file2)

temp_array = np.linspace(0.1, 1.1, 11)

in1d(temp_array, d1['T'])
