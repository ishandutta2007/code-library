import time
import math

start = time.time()
from decimal import Decimal, getcontext, ROUND_FLOOR

# We use decimal to control decimal floating issue

getcontext().prec = 512


def solve_pells_equation(D):
    # Find the solution in the form of
    # x^2 - D*y^2 = 1
    # https://www.alpertron.com.ar/METHODS.HTM#Hyperb

    x = int(math.sqrt(D))
    if x * x == D:
        return -1, -1

    r = Decimal(4 * D).sqrt()
    t = Decimal(2)
    c = (r / t).to_integral_exact(rounding=ROUND_FLOOR)
    r, t = t, r - t * c
    xn_2 = 1
    xn_1 = c
    yn_2 = 0
    yn_1 = 1
    idx = 1
    while True:
        c = (r / t).to_integral_exact(rounding=ROUND_FLOOR)
        r, t = t, r - t * c
        xn, yn = c * xn_1 + xn_2, c * yn_1 + yn_2
        if xn**2 - D * yn**2 == 1:
            if xn > 0 and yn > 0:
                return xn, yn
            elif xn < 0 and yn < 0:
                return -xn, -yn
        idx += 1
        yn_2, yn_1 = yn_1, yn
        xn_2, xn_1 = xn_1, xn


for D in range(1, 1000):
    x, y = solve_pells_equation(D)
    print(x, y)
print(time.time() - start)
