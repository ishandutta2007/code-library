#include <bits/stdc++.h>

using namespace std;
using namespace atcoder;

// const int N = 1e6 + 3, mod = 1e6 + 3;

template <const int32_t MOD> struct modint {
  int32_t value;
  modint() = default;
  modint(int32_t value_) : value(value_) {}
  inline modint<MOD> operator+(modint<MOD> other) const {
    int32_t c = this->value + other.value;
    return modint<MOD>(c >= MOD ? c - MOD : c);
  }
  inline modint<MOD> operator-(modint<MOD> other) const {
    int32_t c = this->value - other.value;
    return modint<MOD>(c < 0 ? c + MOD : c);
  }
  inline modint<MOD> operator*(modint<MOD> other) const {
    int32_t c = (int64_t)this->value * other.value % MOD;
    return modint<MOD>(c < 0 ? c + MOD : c);
  }
  inline modint<MOD> &operator+=(modint<MOD> other) {
    this->value += other.value;
    if (this->value >= MOD)
      this->value -= MOD;
    return *this;
  }
  inline modint<MOD> &operator-=(modint<MOD> other) {
    this->value -= other.value;
    if (this->value < 0)
      this->value += MOD;
    return *this;
  }
  inline modint<MOD> &operator*=(modint<MOD> other) {
    this->value = (int64_t)this->value * other.value % MOD;
    if (this->value < 0)
      this->value += MOD;
    return *this;
  }
  inline modint<MOD> operator-() const {
    return modint<MOD>(this->value ? MOD - this->value : 0);
  }
  modint<MOD> pow(uint64_t k) const {
    modint<MOD> x = *this, y = 1;
    for (; k; k >>= 1) {
      if (k & 1)
        y *= x;
      x *= x;
    }
    return y;
  }
  modint<MOD> inv() const { return pow(MOD - 2); } // MOD must be a prime
  inline modint<MOD> operator/(modint<MOD> other) const {
    return *this * other.inv();
  }
  inline modint<MOD> operator/=(modint<MOD> other) {
    return *this *= other.inv();
  }
  inline bool operator==(modint<MOD> other) const {
    return value == other.value;
  }
  inline bool operator!=(modint<MOD> other) const {
    return value != other.value;
  }
  inline bool operator<(modint<MOD> other) const { return value < other.value; }
  inline bool operator>(modint<MOD> other) const { return value > other.value; }
};
template <int32_t MOD> modint<MOD> operator*(int32_t value, modint<MOD> n) {
  return modint<MOD>(value) * n;
}
template <int32_t MOD> modint<MOD> operator*(int64_t value, modint<MOD> n) {
  return modint<MOD>(value % MOD) * n;
}
template <int32_t MOD> istream &operator>>(istream &in, modint<MOD> &n) {
  return in >> n.value;
}
template <int32_t MOD> ostream &operator<<(ostream &out, modint<MOD> n) {
  return out << n.value;
}

using ll = long long;
using mint = modint;

const int mod = 200003;

mint fact[mod];
mint ifact[mod];

// binom(n, k) mod 200003
mint binom(ll n, ll k) {
  if (k < 0 or k > n)
    return 0;
  mint res = 1;
  while (n) {
    ll n0 = n % mod;
    ll k0 = k % mod;
    if (n0 < k0)
      return 0;
    res *= fact[n0] * ifact[k0] * ifact[n0 - k0];
    n /= mod;
    k /= mod;
  }
  return res;
}

int main() {
  mint::set_mod(mod);
  fact[0] = 1;
  for (int i = 1; i < mod; i++)
    fact[i] = fact[i - 1] * i;
  ifact[mod - 1] = fact[mod - 1].inv();
  for (int i = mod - 1; i >= 1; i--)
    ifact[i - 1] = ifact[i] * i;
  cout << binom(100, 7) << endl;
}
