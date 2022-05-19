
inline int64 mulmod(int64 a, int64 b, int64 mod) {
  int64 res = (a * ((long double)b / (long double)mod));
  res = a * b - res * mod;
  if (res >= mod)
    res -= mod;
  if (res < 0)
    res += mod;
  return res;
}

inline bool witness(uint64 n, uint64 s, uint64 d, uint64 a) {
  uint64 x = power(a, d, n);
  uint64 y;

  do {
    y = mulmod(x, x, n);
    if (y == 1 && x != 1 && x != n - 1)
      return false;
    x = y;
    --s;
  } while (s);
  if (y != 1)
    return false;
  return true;
}

inline bool is_prime_mr(int64 n) {
  if (((!(n & 1)) && n != 2) || (n < 2) || (n % 3 == 0 && n != 3))
    return false;
  if (n <= 3)
    return true;

  uint64 d = n >> 1;
  uint64 s = 1;
  while (!(d & 1)) {
    d >>= 1;
    ++s;
  }

  if (n < 1373653)
    return witness(n, s, d, 2) && witness(n, s, d, 3);
  if (n < 9080191)
    return witness(n, s, d, 31) && witness(n, s, d, 73);
  if (n < 4759123141LL)
    return witness(n, s, d, 2) && witness(n, s, d, 7) && witness(n, s, d, 61);
  if (n < 1122004669633LL)
    return witness(n, s, d, 2) && witness(n, s, d, 13) &&
           witness(n, s, d, 23) && witness(n, s, d, 1662803);
  if (n < 2152302898747LL)
    return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) &&
           witness(n, s, d, 7) && witness(n, s, d, 11);
  if (n < 3474749660383LL)
    return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) &&
           witness(n, s, d, 7) && witness(n, s, d, 11) && witness(n, s, d, 13);
  return witness(n, s, d, 2) && witness(n, s, d, 3) && witness(n, s, d, 5) &&
         witness(n, s, d, 7) && witness(n, s, d, 11) && witness(n, s, d, 13) &&
         witness(n, s, d, 17);
}