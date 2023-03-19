def _integer_nthroot_python(y, n):
    if y in (0, 1):
        return y, True
    if n == 1:
        return y, True
    if n == 2:
        # from mpmath.libmp import sqrtrem as mpmath_sqrtrem
        # x, rem = mpmath_sqrtrem(y)
        x = isqrt(y)
        return x, (y == (x * x))
    if n >= y.bit_length():
        return 1, False
    try:
        guess = int(y ** (1.0 / n) + 0.5)
    except OverflowError:
        exp = math.log(y, 2) / n
        if exp > 53:
            shift = int(exp - 53)
            guess = int(2.0 ** (exp - shift) + 1) << shift
        else:
            guess = int(2.0 ** exp)
    if guess > 2 ** 50:
        # Newton iteration
        xprev, x = -1, guess
        while 1:
            t = x ** (n - 1)
            xprev, x = x, ((n - 1) * x + y // t) // n
            if abs(x - xprev) < 2:
                break
    else:
        x = guess
    # Compensate
    t = x ** n
    while t < y:
        x += 1
        t = x ** n
    while t > y:
        x -= 1
        t = x ** n
    return int(x), t == y  # int converts long to int if possible
