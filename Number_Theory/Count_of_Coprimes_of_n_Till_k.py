import sympy
import bisect
import random
import math

# Supports till 10^7 (to be precise till 2*3*5*7*11*13*17*19)

# Example scenario
# n = 100
# m = 40
# d = 4
# we are to find sum of all places under(inclusive) 40 where 100 has gcd(100,i)=4 ie
# 4+8+12+16+24+28+32+36, (but not 20 and 40)
# basically 4*(sum of coprimes of 100/4 under 40/4 ie 1+2+3+4+6+7+8+9,(but not 5 and 10))


def count_of_coprimes_of_n_till_k(n, k):
    primes = list(sympy.factorint(n, multiple=True))
    primes = list(set(primes))
    primes.sort()
    ans = k
    for i0 in range(0, len(primes)):
        p1 = primes[i0]
        pp1 = p1
        # we have to remove the series of p under k ie p+2p+3p....[k/p]p
        t1 = k // pp1
        # print("pp1=", pp1, "k=", k, "terms=", t1)
        ans -= t1
        for i1 in range(i0 + 1, len(primes)):
            p2 = primes[i1]
            pp2 = pp1 * p2
            # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
            t2 = k // pp2
            # print("pp2=", pp2, "k=", k, "terms=", t2)
            ans += t2
            for i2 in range(i1 + 1, len(primes)):
                p3 = primes[i2]
                pp3 = pp2 * p3
                # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                t3 = k // pp3
                # print("pp3=", pp3, "k=", k, "terms=", t3)
                ans -= t3
                for i3 in range(i2 + 1, len(primes)):
                    p4 = primes[i3]
                    pp4 = pp3 * p4
                    # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                    t4 = k // pp4
                    # print("pp4=", pp4, "k=", k, "terms=", t4)
                    ans += t4
                    for i4 in range(i3 + 1, len(primes)):
                        p5 = primes[i4]
                        pp5 = pp4 * p5
                        # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                        t5 = k // pp5
                        # print("pp5=", pp5, "k=", k, "terms=", t5)
                        ans -= t5
                        for i5 in range(i4 + 1, len(primes)):
                            p6 = primes[i5]
                            pp6 = pp5 * p6
                            # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                            t6 = k // pp6
                            # print("pp6=", pp6, "k=", k, "terms=", t6)
                            ans += t6
                            for i6 in range(i5 + 1, len(primes)):
                                p7 = primes[i6]
                                pp7 = pp6 * p7
                                # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                                t7 = k // pp7
                                # print("pp7=", pp7, "k=", k, "terms=", t7)
                                ans -= t7
                                for i7 in range(i6 + 1, len(primes)):
                                    p8 = primes[i7]
                                    pp8 = pp7 * p8
                                    # we have to add back the series of pp2 under k ie pp2+2pp2+3pp2....[k/pp2]p
                                    t8 = k // pp8
                                    # print("pp8=", pp8, "k=", k, "terms=", t8)
                                    ans += t8
    return ans


def count_of_coprimes_of_n_till_k_brute(n, k):
    ans = 0
    for i in range(1, k + 1):
        if math.gcd(i, n) == 1:
            ans += 1
    return ans


# print(count_of_coprimes_of_n_till_k(25, 10))

algo_ans = count_of_coprimes_of_n_till_k(50, 20)
print(algo_ans)
assert algo_ans == 8  # 1 , 3 , 7 , 9 , 11 , 13 , 17 , 19
# (not 2,4,6,8,10,12,14,16,18,20 and not 5,10,15,20 note 10, 20 have been removed twice)


algo_ans = count_of_coprimes_of_n_till_k(50, 30)
print(algo_ans)
assert algo_ans == 12  # 1 , 3 , 7 , 9 , 11 , 13 , 17 , 19 , 21 , 23 , 27 , 29
# (not 2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 and not 5,10,15,20,25,30 note 10, 20,30 have been removed twice)


algo_ans = count_of_coprimes_of_n_till_k(60, 30)
print(algo_ans)
assert algo_ans == 8  # 1 , 7 , 11 , 13 , 17 , 19 , 23 , 29
# not 2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 and not 3,6,9,12,15,18,21,24,27,30 and not 5,10,15,20,25,30
# note 6,12,18,24,30 have been removed twice and 10,20,30 has been removed twice and 15,30 removed twice
# note 30 has been removed twice then added back thrice


algo_ans = count_of_coprimes_of_n_till_k(210, 20)
print(algo_ans)
assert algo_ans == 5  # 1 , 11 , 13 , 17 , 19


algo_ans = count_of_coprimes_of_n_till_k(210 * 11, 20)
print(algo_ans)
assert algo_ans == 4  # 1 , 13 , 17 , 19


algo_ans = count_of_coprimes_of_n_till_k(210 * 11 * 13, 20)
print(algo_ans)
assert algo_ans == 3  # 1 , 17 , 19


print(count_of_coprimes_of_n_till_k(2 * 2 * 2 * 3 * 3 * 3 * 7 * 7 * 17, 10820))
print(count_of_coprimes_of_n_till_k_brute(2 * 2 * 2 * 3 * 3 * 3 * 7 * 7 * 17, 10820))


for i in range(10 ** 6, 10 ** 7):
    r = random.randint(1, i)
    algo_ans = count_of_coprimes_of_n_till_k(i, r)
    brute_ans = count_of_coprimes_of_n_till_k_brute(i, r)
    print(i, r, algo_ans, brute_ans)
    assert algo_ans == brute_ans
