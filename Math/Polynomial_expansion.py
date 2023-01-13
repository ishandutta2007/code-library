from sympy import Symbol, prod, expand, Poly

x = Symbol("x")
# primes = [5, 7, 11, 13, 17, 19, 23, 29]
primes = [19]# [10000019]
for p in primes:
    expr = prod([x + i for i in range(1, p)])
    expanded_expr = expand(expr)
    print(expanded_expr)

    for k in range(0, p):
        print("P(", k * p, ")=", expanded_expr.subs(x, k * p) % (p * p))

    poly_f = Poly(expanded_expr, x)
    coeffs_f = poly_f.coeffs()
    print(p, coeffs_f)
    coeffs_f = [x % (p * p) for x in coeffs_f]
    last = (coeffs_f[-1] + 1) / p
    print(p, coeffs_f)
    # coeffs_f = [x / (p) for x in coeffs_f]
    # print(p, coeffs_f, last)
    # coeffs_f.sort()
    # print(p, coeffs_f)
    print("====")

    # mod_t = tuple(x % (p*p) for x in expanded_expr.args)
    # print(mod_t)
