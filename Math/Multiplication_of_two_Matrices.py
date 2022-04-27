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
B = [[25, 24, 23, 22], [21, 20, 19, 18], [17, 16, 15, 14]]


def matrix_muliplication(A, B):
    res = [
        [sum(a * b for a, b in zip(A_row, B_col)) for B_col in zip(*B)] for A_row in A
    ]
    return res


print(matrix_muliplication(A, B))
