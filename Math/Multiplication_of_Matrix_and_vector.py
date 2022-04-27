# Python program to multiply two matrices without numpy
import copy
import pprint as pp
import numpy as np
import os
import sys
import time

np.set_printoptions(
    edgeitems=30, linewidth=100000, formatter={"float": lambda x: "{0:0.2f}".format(x)}
)

A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
B = [25, 24, 23]


def matrix_vector_muliplication(mat, vec):
    ans = []
    for row in mat:
        ans.append(sum([x * y for x, y in zip(row, vec)]))
    return ans


print(matrix_vector_muliplication(A, B))
