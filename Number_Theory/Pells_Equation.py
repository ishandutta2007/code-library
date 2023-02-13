import time
import math

start = time.time()


def solvePell(n):
    x = int(math.sqrt(n))
    if x * x == n:
        return -1, -1
    y, z, r = x, 1, x << 1
    e1, e2 = 1, 0
    f1, f2 = 0, 1
    while True:
        y = r * z - y
        z = (n - y * y) // z
        r = (x + y) // z

        e1, e2 = e2, e1 + e2 * r
        f1, f2 = f2, f1 + f2 * r

        a, b = f2 * x + e2, f2
        if a * a - n * b * b == 1:
            return a, b


for D in range(1, 1000):
    x, y = solvePell(D)
    print(x, y)
print(time.time() - start)
