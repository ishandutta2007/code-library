#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using i64 = long long;
using i128 = __int128;
using u64 = unsigned long long;
using u128 = __uint128_t;
using f64 = double;
using f128 = __float128;

#pragma once
#include <iostream>

template <int MOD> class Modular {
  // using ModInt = Modular<1'000'000'007>;
public:
  Modular(long long v = 0) {
    value = v % MOD;
    if (value < 0)
      value += MOD;
  }

  Modular(long long a, long long b) : value(0) {
    *this += a;
    *this /= b;
  }

  Modular &operator+=(Modular const &b) {
    value += b.value;
    if (value >= MOD)
      value -= MOD;
    return *this;
  }

  Modular &operator-=(Modular const &b) {
    value -= b.value;
    if (value < 0)
      value += MOD;
    return *this;
  }

  Modular &operator*=(Modular const &b) {
    value = (long long)value * b.value % MOD;
    return *this;
  }

  friend Modular power(Modular a, long long e) {
    Modular res = 1;
    while (e) {
      if (e & 1)
        res *= a;
      a *= a;
      e >>= 1;
    }
    return res;
  }

  friend Modular inverse(Modular a) { return power(a, MOD - 2); }

  Modular &operator/=(Modular const &b) { return *this *= inverse(b); }

  friend Modular operator+(Modular a, Modular const b) { return a += b; }

  friend Modular operator-(Modular a, Modular const b) { return a -= b; }

  friend Modular operator-(Modular const a) { return 0 - a; }

  friend Modular operator*(Modular a, Modular const b) { return a *= b; }

  friend Modular operator/(Modular a, Modular const b) { return a /= b; }

  friend std::ostream &operator<<(std::ostream &os, Modular const &a) {
    return os << a.value;
  }

  friend std::istream &operator>>(std::istream &is, Modular &a) {
    is >> a.value;
    a.value %= MOD;
    if (a.value < 0)
      a.value += MOD;
    return is;
  }

  friend bool operator==(Modular const &a, Modular const &b) {
    return a.value == b.value;
  }

  friend bool operator!=(Modular const &a, Modular const &b) {
    return a.value != b.value;
  }

  friend Modular &operator++(Modular &a, int) { return a += 1; }

  friend Modular operator++(Modular const &a, int) { return Modular(a)++; }

  friend Modular &operator--(Modular &a, int) { return a -= 1; }

  friend Modular operator--(Modular const &a, int) { return Modular(a)--; }

  int value;
  static const int MOD_value = MOD;
};

int main(int argc, char **argv) {
  // multiplication
  Modular<13> a = 5, m = 12;
  a *= m; // a = 519994069
  cout << a << endl;

  // Vector Multiplication
  vector<Modular<13>> p(10);
  vector<Modular<13>> q(10);
  for (int i = 0; i < p.size(); i++) {
    p[i] = i;
    q[i] = i * i;
    q[i] *= p[i];
    cout << q[i] << " ";
  }

  // using mint = dynamic_mint;
  // mint::set_mod(17);
  // vector<mint>a(10);
  // vector<mint>b(10);
  // vector<mint>c(10);
  // a.push_back((mint)13);
  // b.push_back((mint)(7));
  // for(int i = 0;i < 10; i++){
  // c[i]=a[i]+b[i];
  //	  cout<<a[i]<<b[i]<<c[i];
  //}
  return 0;
}