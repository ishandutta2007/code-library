#include <cassert>
#include <cmath>
#include <cstdio>

#include <algorithm>
#include <vector>

#include <array>
#include <initializer_list>
#include <tuple>

using i32 = int32_t;
using i64 = int64_t;
using u8 = uint8_t;
using u32 = uint32_t;
using u64 = uint64_t;
using u128 = __uint128_t;

template <typename T> struct fast_div;

template <> struct fast_div<int> {
  constexpr fast_div() : m(), s(), x() {}
  constexpr fast_div(int n)
      : m(n), s(std::__lg(n - 1)), x(((u128(1) << (s + 64)) + n - 1) / n) {}
  constexpr friend u64 operator/(u64 n, const fast_div &d) {
    return (u128(n) * d.x >> d.s) >> 64;
  }
  constexpr friend int operator%(u64 n, const fast_div &d) {
    return n - n / d * d.m;
  }
  constexpr std::pair<i64, int> divmod(u64 n) const {
    u64 q = n / *this;
    return {q, n - q * m};
  }
  int m;
  int s;
  u64 x;
};

template <typename T> struct promote_integral;
template <> struct promote_integral<i32> { using type = i64; };

template <typename T, typename M>
inline constexpr T mod_inv_pe(T a, T p, T e, const M &pe) {
  using TT = typename promote_integral<T>::type;
  T ret = pow_mod(a, p - 2, pe);
  for (T i = 1; i < e; i <<= 1)
    ret = TT(ret) * (2 + pe - TT(ret) * a % pe) % pe;
  return ret;
}

std::vector<std::pair<int, int>> factors(int n) {
  std::vector<std::pair<int, int>> ret;
  for (int i = 2; i64(i) * i <= n; ++i) {
    if (n % i == 0) {
      int e = 1;
      n /= i;
      while (n % i == 0)
        n /= i, ++e;
      ret.emplace_back(i, e);
    }
  }
  if (n > 1)
    ret.emplace_back(n, 1);
  return ret;
}

template <typename T, typename M>
inline constexpr T pow_mod(T a, T e, const M &m) {
  using TT = typename promote_integral<T>::type;
  T ret = 1;
  for (; e; e >>= 1, a = TT(a) * a % m) {
    if (e & 1)
      ret = TT(ret) * a % m;
  }
  return ret;
}

template <typename T> struct cint {
private:
  using TT = typename promote_integral<T>::type;

public:
  cint() : d() {}
  constexpr explicit cint(T n) : d(n) {}
  friend constexpr TT operator/(TT lhs, const cint &rhs) { return lhs / rhs.d; }
  friend constexpr T operator/(T lhs, const cint &rhs) { return lhs / rhs.d; }
  friend constexpr T operator/(const cint &lhs, const cint &rhs) {
    return lhs.m / rhs.d;
  }
  friend constexpr TT operator%(TT lhs, const cint &rhs) { return lhs % rhs.d; }
  friend constexpr T operator%(T lhs, const cint &rhs) { return lhs % rhs.d; }
  friend constexpr std::pair<TT, T> divmod(TT n, const cint &d) {
    return d.d.divmod(n);
  }
  constexpr operator T() const { return d.m; }
  fast_div<T> d;
};
using c32 = cint<int>;

struct binomial_mod_pe {
  binomial_mod_pe() {}
  binomial_mod_pe(int p, int e, size_t pre_size = 1 << 14)
      : p(p), e(e), ppows(e + 1, 1) {
    for (int i = 1; i <= e; ++i)
      ppows[i] = ppows[i - 1] * p;
    this->pe = c32(ppows[e]), this->p1 = c32(p - 1);
    int ppow = pe;
    while (size_t(ppow) * p <= pre_size)
      ppow *= p;
    this->pp = c32(ppow), this->q = c32(ppow / pe * p);
    facts.resize(ppow, 1);
    for (int qq = 1, l = pp / p; qq < q; qq *= p, l /= p) {
      for (int i = 0; i < l; ++i)
        for (int j = i * p + 1; j < i * p + p; ++j) {
          facts[j * qq] = j;
        }
    }
    for (int i = 1; i < ppow; ++i)
      facts[i] = i64(facts[i - 1]) * facts[i] % pe;
    mask = (facts[ppow - 1] == pe - 1);
    ds.resize(q, 0);
    for (int i = 0; i < ppow / pe; ++i) {
      for (int j = 0, s = ds[i]; j < p; ++j)
        ds[i * p + j] = s + j;
    }
  }

  int operator()(i64 n, i64 m) const {
    i64 r = n - m;
    int numer = 1, denom = 1, x = 0, s = 0;
    if (e > 1) {
      while (n > 0) {
        int n1, nr1;
        std::tie(n1, nr1) = divmod(n, pp);
        int m1, mr1;
        std::tie(m1, mr1) = divmod(m, pp);
        int r1, rr1;
        std::tie(r1, rr1) = divmod(r, pp);
        numer = i64(numer) * facts[nr1] % pp;
        denom = i64(denom) * facts[mr1] % pp * facts[rr1] % pp;
        s += n1 - m1 - r1;
        int nr2;
        std::tie(n, nr2) = divmod(n, q);
        int mr2;
        std::tie(m, mr2) = divmod(m, q);
        int rr2;
        std::tie(r, rr2) = divmod(r, q);
        x += ds[mr2] + ds[rr2] - ds[nr2];
      }
    } else {
      while (n > 0) {
        int nr;
        std::tie(n, nr) = divmod(n, pp);
        int mr;
        std::tie(m, mr) = divmod(m, pp);
        int rr;
        std::tie(r, rr) = divmod(r, pp);
        s += n - m - r;
        numer = i64(numer) * facts[nr] % pp;
        denom = i64(denom) * facts[mr] % pp * facts[rr] % pp;
        x += ds[mr] + ds[rr] - ds[nr];
      }
    }
    if (x >= e * (p - 1))
      return 0;
    if (p > 2)
      x /= p1;
    int ret = i64(numer) * mod_inv_pe(denom, p, e - x, pe) % pe * ppows[x] % pe;
    if ((s & mask) && ret > 0)
      ret = pe - ret;
    return ret;
  }

  int p, e, mask;
  c32 pe, q, pp, p1;
  std::vector<int> ppows, facts, ds;
};

struct binomial_mod {
  binomial_mod(int mod, size_t pre_size = 1 << 14)
      : mod(mod), f(factors(mod)), binom_pp(f.size()), iprods(f.size()) {
    if (f.size() == 1)
      pre_size = 1 << 20;
    int prod = 1;
    for (size_t i = 0; i < f.size(); ++i) {
      int p, e;
      std::tie(p, e) = f[i];
      binom_pp[i] = binomial_mod_pe(p, e, pre_size);
      iprods[i] = mod_inv_pe(prod, p, e, binom_pp[i].pe);
      prod *= binom_pp[i].pe;
    }
  }

  int operator()(i64 n, i64 m) const {
    assert(m >= 0);
    int s = 0;
    if (n < 0)
      s = (m & 1), n = -n + m - 1;
    if (n < m)
      return 0;
    if (m == 0)
      return (mod > 1);
    int ret = 0, prod = 1;
    for (size_t i = 0; i < binom_pp.size(); ++i) {
      const auto &d = binom_pp[i].pe;
      ret += i64(binom_pp[i](n, m) + d - ret % d) * iprods[i] % d * prod;
      prod *= d;
    }
    return (s && ret > 0) ? mod - ret : ret;
  }

  int mod;
  std::vector<std::pair<int, int>> f;
  std::vector<binomial_mod_pe> binom_pp;
  std::vector<int> iprods;
};

#define getchar getchar_unlocked
#define putchar putchar_unlocked

i64 get_i64() {
  int c;
  i64 n;
  while ((c = getchar()) < '0')
    ;
  n = c - '0';
  while ((c = getchar()) >= '0')
    n = n * 10 + (c - '0');
  return n;
}

void put_int(int n) {
  if (n < 0)
    putchar('-'), n = -n;
  int res[11], i = 0;
  do {
    res[i++] = n % 10, n /= 10;
  } while (n);
  while (i)
    putchar(res[--i] + '0');
  putchar('\n');
}

void solve() {
  int T = get_i64(), m = get_i64();
  const auto binom = binomial_mod(m);
  for (int i = 0; i < T; ++i) {
    i64 n = get_i64();
    i64 m = get_i64();
    put_int(binom(n, m));
  }
}

int main() {
  solve();
  return 0;
}

// https://judge.yosupo.jp/submission/46090
