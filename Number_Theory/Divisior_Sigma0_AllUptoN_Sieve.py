from numba import jit
import numpy as np

@jit(nopython=True)
def divcount_all_upto(N):
    divcount = np.full(int(N + 1), 2, dtype=np.int32)
    for i in range(2, int(N + 1) // 2):
        divcount[2 * i :: i] += 1
    return divcount

divcount_all_upto(1e7)
