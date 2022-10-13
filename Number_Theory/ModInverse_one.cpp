//O(log mod)
ll mod_inv(x, mod) { return mod - (mod / x) * mod_inv(mod % x, mod) % mod; }

//O(log mod_prime)
ll mod_inv(x, mod_prime) { return bigpow(x, mod_prime - 2, mod_prime); }

//O(log mod)
ll mod_inv(x, mod) { return bigpow(x, etf(mod) - 1, mod); }

mint inv() const {
  if (prime) {
    assert(_v);
    return pow(umod() - 2);
  } else {
    auto eg = internal::inv_gcd(_v, m);
    assert(eg.first == 1);
    return eg.second;
  }
}