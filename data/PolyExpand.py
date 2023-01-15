from sympy import Symbol, prod, expand, Poly
import time

p = 10**7 + 19
pp = p * p
start = time.time()
x = Symbol("x")
expanded_expr = 1
coeffs = []
for i in range(1, 3162 + 1):
    expanded_expr = expanded_expr * (x + i)
    expanded_expr = expand(expanded_expr)
    if i % 2 == 0:
        p = Poly(expanded_expr, x)
        coeffs = []
        for c in p.coeffs():
            if c >= pp:
                coeffs.append(c % pp)
            else:
                coeffs.append(c)

        coeffs.reverse()
        print(i, time.time() - start)
        print("=====")

with open("poly_coeffs_3162.txt", "w") as f:
    f.write("%s\n" % coeffs)
