
/* Fermat pseudoprime */
int is_pseudoprime(UV const n, UV a) {
  if (n < 4)
    return (n == 2 || n == 3);
  if (!(n & 1) && !(a & 1))
    return 0;
  if (a < 2)
    croak("Base %" UVuf " is invalid", a);
  if (a >= n) {
    a %= n;
    if (a <= 1)
      return (a == 1);
    if (a == n - 1)
      return !(a & 1);
  }

#if USE_MONTMATH
  if (n & 1) { /* The Montgomery code only works for odd n */
    const uint64_t npi = mont_inverse(n), mont1 = mont_get1(n);
    const uint64_t monta = (a == 2) ? mont_get2(n) : mont_geta(a, n);
    return mont_powmod(monta, n - 1, n) == mont1;
  }
#endif
  return powmod(a, n - 1, n) == 1; /* a^(n-1) = 1 mod n */
}
