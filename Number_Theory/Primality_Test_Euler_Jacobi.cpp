/* Euler (aka Euler-Jacobi) pseudoprime:  a^((n-1)/2) = (a|n) mod n */
int is_euler_pseudoprime(UV const n, UV a) {
  if (n < 5)
    return (n == 2 || n == 3);
  if (!(n & 1))
    return 0;
  if (a < 2)
    croak("Base %" UVuf " is invalid", a);
  if (a > 2) {
    if (a >= n) {
      a %= n;
      if (a <= 1)
        return (a == 1);
      if (a == n - 1)
        return !(a & 1);
    }
    if ((n % a) == 0)
      return 0;
  }
  {
#if USE_MONTMATH
    const uint64_t npi = mont_inverse(n), mont1 = mont_get1(n);
    const uint64_t monta = mont_geta(a, n);
    UV ap = mont_powmod(monta, (n - 1) >> 1, n);
    if (ap != mont1 && ap != n - mont1)
      return 0;
    if (a == 2) {
      uint32_t nmod8 = n & 0x7;
      return (nmod8 == 1 || nmod8 == 7) ? (ap == mont1) : (ap == n - mont1);
    } else {
      return (kronecker_uu(a, n) >= 0) ? (ap == mont1) : (ap == n - mont1);
    }
#else
    UV ap = powmod(a, (n - 1) >> 1, n);
    if (ap != 1 && ap != n - 1)
      return 0;
    if (a == 2) {
      uint32_t nmod8 = n & 0x7;
      return (nmod8 == 1 || nmod8 == 7) ? (ap == 1) : (ap == n - 1);
    } else {
      return (kronecker_uu(a, n) >= 0) ? (ap == 1) : (ap == n - 1);
    }
#endif
  }
}

// https://github.com/danaj/Math-Prime-Util/blob/7cf12584381aea0878b7b8d438a51b6c2f3a40b6/primality.c
