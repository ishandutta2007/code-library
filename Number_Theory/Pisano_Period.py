def getPisanoPeriod(m):
    f_0 = 1
    f_1 = 1
    pisano_period = 1
    while f_1 != 1 or f_0 != 0:
        f_0, f_1 = f_1, (f_0 + f_1) % m
        pisano_period += 1

    return pisano_period
