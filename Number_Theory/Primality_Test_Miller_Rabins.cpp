// O(it * (logn)^3), it = number of rounds performed
inline bool miller_rabin(ll n) {
  if (n <= 2 || (n & 1 ^ 1))
    return (n == 2);
  if (n < P)
    return spf[n] == n;
  ll c, d, s = 0, r = n - 1;
  for (; !(r & 1); r >>= 1, s++) {
  }
  // each iteration is a round
  for (int i = 0; primes[i] < n && primes[i] < 32; i++) {
    c = pow_mod(primes[i], r, n);
    for (int j = 0; j < s; j++) {
      d = mul_mod(c, c, n);
      if (d == 1 && c != 1 && c != (n - 1))
        return false;
      c = d;
    }
    if (c != 1)
      return false;
  }
  return true;
}