# Last updated: 16-Jul-2018
# Udaya Maurya (udaya_cbscients@yahoo.com, telegram: https://content-available-to-author-only.me/udy11)
# Source: https://g...content-available-to-author-only...b.com/udy11, https://g...content-available-to-author-only...b.com/udy11
# Function to find non-negative solutions of Generalized Pell's
# Equation using Convergents of Continued Fractions

# Generalized Pell's Equation is defined as:
# x**2 - D * y**2 = N
# whose integer solutions in x and y are seeked for a positive
# non-square integer D and integer N

# ALL YOU NEED TO DO:
# Call pellgeq(D, N, z) with posiitve non-square integer D,
#   integer N and positive integer z for first z solutions
#   (sorted by ascending x)

# NOTE
# Don't forget to import math
# Only non-negative solutions are considered because if
#   (x, y) is a solution, then (-x, y), (x, -y) and
#   (-x, -y) are obviously also solutions

import math


def is_sqrn(n):
    """(int) -> logical
    To check if a non-negative integer n is a square number
    Sq(n) = n * n
    """
    if n % 4 == 0 or n % 8 == 1:
        r = round(math.sqrt(n))
        if r * r == n:
            return True
    return False


def cfsq(r):
    """(int) -> tuple
    Returns Continued Fractions of sqrt(r)
    Assumes input to be positive integer
    """
    ss = math.sqrt(r)
    n = int(ss)
    m1 = 1
    p1 = n
    m = p1
    p = r - m * m
    c1 = [n]
    w = []
    while True:
        if [m1, p1] in w:
            break
        else:
            w.append([m1, p1])
        s = (ss + m) / p
        n = int(s)
        m1 = p
        p1 = p * n - m
        m = p1
        p = (r - m * m) // m1
        c1.append(n)
    return (c1[0], c1[1:])


def pfp(x, y, D):
    """(int, int, int) -> logical
    Returns true if parameters satisfy positive Pell's Equation
    """
    return x * x - D * y * y == 1


def pfg(x, y, D, N):
    """(int, int, int, int) -> logical
    Returns true if parameters satisfy generalized Pell's Equation
    """
    return x * x - D * y * y == N


def pellfp(m):
    """Returns the fundamental solution of the positive Pell's Equation"""
    b = cfsq(m)
    n = len(b[1])
    x = [b[0], 0, 0]
    y = [1, 0, 0]
    if pfp(x[0], y[0], m):
        xfp = x[0]
        yfp = y[0]
    elif n > 0:
        x[1] = b[0] * b[1][0] + 1
        y[1] = b[1][0]
        if pfp(x[1], y[1], m):
            xfp = x[1]
            yfp = y[1]
        elif n > 1:
            k = 1
            while True:
                kk = k % n
                x[2] = b[1][kk] * x[1] + x[0]
                y[2] = b[1][kk] * y[1] + y[0]
                if pfp(x[2], y[2], m):
                    xfp = x[2]
                    yfp = y[2]
                    break
                else:
                    x[0] = x[1]
                    x[1] = x[2]
                    y[0] = y[1]
                    y[1] = y[2]
                    k += 1
    return (xfp, yfp)


def pellgeq(m, q, z):
    """(int, int, int) -> list of tuple
    Solves Generalized Pell's Equation x**2 - D * y**2 = N
    using Convergents of ContinuedbFractions. Returns
    list of tuple of first z solutions of the equation
    (sorted by ascending x)
    """
    if is_sqrn(m):
        if q == 0:
            return [(j * round(math.sqrt(m)), j) for j in range(z)]
        else:
            return
    if q == 0:
        return [(0, 0)]
    else:
        xf, yf = pellfp(m)
    sf = []
    sqm = math.sqrt(m)
    uf = xf + sqm * yf
    for i in range(1 + math.ceil(math.sqrt(abs(q) * uf))):
        for j in range(1 + math.ceil(math.sqrt(abs(q) * uf) / sqm)):
            if pfg(i, j, m, q):
                sf.append((i, j))
                if i > sqm * j:
                    sf.append((i, -j))
                else:
                    sf.append((-i, j))
    sf = sorted(set(sf))
    if len(sf) == 0:
        return
    sol = []
    sfc = [0 for i in range(len(sf))]
    ns = 0
    while ns < z:
        s = min(sf)
        k = sf.index(s)
        if (s[0] >= 0 and s[1] >= 0) and (ns == 0 or (ns > 0 and s != sol[-1])):
            sol.append(s)
            ns += 1
        sf[k] = (xf * sf[k][0] + m * yf * sf[k][1], xf * sf[k][1] + yf * sf[k][0])
        sfc[k] += 1
    return sol


D = 5
N = 1
sln = pellgeq(D, N, 20)
for s in sln:
    print(s, pfg(s[0], s[1], D, N))
