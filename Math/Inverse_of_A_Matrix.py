# import random
# import operator as op
# from functools import reduce
import copy
import pprint as pp
import os
import sys
import time


sys.stdin = open("matrix.in", "r")
# sys.stdout = open("matrix_numpy.out", "w")
n = int(input())
mat = [[0] * (n) for i in range(n)]
for i in range(0, n):
    mat[i] = list(map(int, input().strip().split()))[:n]


def identity_matrix(n):
    m = [[0 for x in range(n)] for y in range(n)]
    for i in range(0, n):
        m[i][i] = 1
    return m


def invert_matrix(A, tol=None):
    """
    Returns the inverse of the passed in matrix.
        :param A: The matrix to be inversed

        :return: The inverse of the matrix A
    """
    # Section 1: Make sure A can be inverted.
    # check_squareness(A)
    # check_non_singular(A)

    # Section 2: Make copies of A & I, AM & IM, to use for row ops
    n = len(A)
    AM = copy.deepcopy(A)
    I = identity_matrix(n)
    IM = copy.deepcopy(I)

    # Section 3: Perform row operations
    indices = list(range(n))  # to allow flexible row referencing ***
    for fd in range(n):  # fd stands for focus diagonal
        fdScaler = 1.0 / AM[fd][fd]
        # FIRST: scale fd row with fd inverse.
        for j in range(n):  # Use j to indicate column looping.
            AM[fd][j] *= fdScaler
            IM[fd][j] *= fdScaler
        # SECOND: operate on all rows except fd row as follows:
        for i in indices[0:fd] + indices[fd + 1 :]:
            # *** skip row with fd in it.
            crScaler = AM[i][fd]  # cr stands for "current row".
            for j in range(n):
                # cr - crScaler * fdRow, but one element at a time.
                AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
                IM[i][j] = IM[i][j] - crScaler * IM[fd][j]

    # Section 4: Make sure IM is an inverse of A with specified tolerance
    # if check_matrix_equality(I, matrix_multiply(A, IM), tol):
    return IM
    # else:
    #     raise ArithmeticError("Matrix inverse out of tolerance.")


def print_2d(matrix, decimal_places=2):
    for i, matrix_row in enumerate(matrix):
        matrix[i] = [round(num, decimal_places) for num in matrix_row]
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = "\t".join("{{:{}}}".format(x) for x in lens)
    table = [fmt.format(*(row)) for row in s]
    print("\n".join(table))


inv_mat = invert_matrix(mat)
print_2d(inv_mat)
