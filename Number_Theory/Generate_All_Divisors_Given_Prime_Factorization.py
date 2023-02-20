
def generate_divisors(prime_factorization):
    divisors = [1]
    for prime in prime_factorization.keys():
        new_divisors = []
        for divisor in divisors:
            for i in range(1, prime_factorization[prime] + 1):
                new_divisors.append((divisor * prime ** i)%MOD)
        divisors.extend(new_divisors)
    return divisors
