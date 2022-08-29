ll fastmod(ll x, ll y = mod - 2) {
  ll res = 1;
  x %= mod;
  while (y > 0) {
    if (y & 1)
      res = res * x % mod;
    x = x * x % mod;
    y /= 2;
  }
  return res;
}