def sum_of_divisors(prime_factorization):
    # sum_of_divisors = (p1^0 + p1^1 + ... + p1^a1) * (p2^0 + p2^1 + ... + p2^a2) * ... * (pk^0 + pk^1 + ... + pk^ak)
    _sum = 1

    # Iterate over each prime factor and its exponent in the prime factorization
    for prime, exponent in prime_factorization.items():
        # Calculate the sum of powers of the prime factor using the geometric series formula
        sum_of_powers = (prime ** (exponent + 1) - 1) // (prime - 1)

        # Multiply the sum by the sum of powers
        _sum = _sum * sum_of_powers % MOD

    return _sum


def sum_of_divisors(prime_factorization, mod_prime):
    # sum_of_divisors = (p1^0 + p1^1 + ... + p1^a1) * (p2^0 + p2^1 + ... + p2^a2) * ... * (pk^0 + pk^1 + ... + pk^ak)
    _sum = 1

    # Iterate over each prime factor and its exponent in the prime factorization
    for prime, exponent in prime_factorization.items():
        # Calculate the sum of powers of the prime factor using the geometric series formula
        sum_of_powers = (prime ** (exponent + 1) - 1) * pow(prime - 1, -1, mod_prime)

        # Multiply the sum by the sum of powers
        _sum = _sum * sum_of_powers % MOD

    return _sum
