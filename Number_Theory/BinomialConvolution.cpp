#include <cstdio>
#include <cassert>
#include <cmath>

#include <algorithm>
#include <vector>

#include <random>

// credits min_25
using namespace std;

using i64 = int64_t;
using i128 = __int128_t;
using u32 = uint32_t;
using u64 = uint64_t;
using u128 = __uint128_t;

constexpr u64 mul_inv(u64 n) {
  u64 x = n;
  for (int i = 0; i < 5; ++i)
    x *= u64(2) - n * x;
  return x;
}

template <typename T,
          std::enable_if_t<std::is_integral<T>::value, nullptr_t> = nullptr>
constexpr T mod_inv(T a, T mod) {
  using S = typename std::make_signed<T>::type;
  S b = mod, s = 1, t = 0, old_a = a;
  while (b > 0) {
    S u = s - t * (a / b);
    s = t;
    t = u;
    S c = a % b;
    a = b;
    b = c;
  }
  if (a > 1) {
    fprintf(stderr, "ModInverseError: gcd(%ld, %ld) = %ld\n", old_a, mod, S(a));
    exit(1);
  };
  return s < 0 ? s + mod : s;
}

template <u64 mod, u64 prim> struct ntt {
  static_assert(!(mod >> 60), "mod should be < 2^{60}");
  static constexpr int lv = __builtin_ctzll(mod - 1);
  static_assert(lv >= 3, "(mod - 1) should be divisible by 8");

  static constexpr u64 modulus() { return mod; }

  struct LMod {
    static constexpr u64 r = -u128(mod) % mod;
    static constexpr u64 v = mul_inv(mod);
    constexpr LMod() : n(0) {}
    constexpr LMod(u64 n) : n(init(n)) {}
    constexpr LMod(const LMod &rhs) : n(rhs.n) {}
    constexpr u64 init(u64 n) const { return reduce(u128(n) * r); }
    constexpr u64 reduce(u128 x) const {
      return u64(x >> 64) + mod - u64((u128(u64(x) * v) * mod) >> 64);
    }
    constexpr u64 get() const { return reduce(n) % mod; }
    constexpr LMod pow(u64 e) const {
      LMod x = LMod(1);
      for (LMod b = *this; e; e >>= 1, b *= b)
        if (e & 1)
          x *= b;
      return x;
    }
    constexpr LMod inv() const { return pow(mod - 2); }
    constexpr LMod &operator*=(const LMod &rhs) {
      n = reduce(u128(n) * rhs.n);
      return *this;
    }
    constexpr LMod operator*(const LMod &rhs) const {
      return LMod(*this) *= rhs;
    }
    constexpr LMod &operator+=(const LMod &rhs) {
      n += rhs.n;
      return *this;
    }
    constexpr LMod operator+(const LMod &rhs) const {
      return LMod(*this) += rhs;
    }
    constexpr LMod &operator-=(const LMod &rhs) {
      n += (mod << 2) - rhs.n;
      return *this;
    }
    constexpr LMod operator-(const LMod &rhs) const {
      return LMod(*this) -= rhs;
    }
    u64 n;
  };

  constexpr ntt() {
    // < (8 - 4 * sqrt(3)) * 2^{60} < 1.0718 * 2^{60}
    rs[lv] = LMod(prim).pow((mod - 1) >> lv);
    for (int i = lv - 1; i >= 0; --i)
      rs[i] = rs[i + 1] * rs[i + 1];
    irs[lv] = LMod(rs[lv]).inv();
    for (int i = lv - 1; i >= 0; --i)
      irs[i] = irs[i + 1] * irs[i + 1];
    itwopows[lv] = LMod(u64(1) << lv).inv();
    LMod two = LMod(2);
    for (int i = lv - 1; i >= 0; --i)
      itwopows[i] = itwopows[i + 1] * two;
    dw[0] = rs[3];
    for (int i = 1; i < lv - 2; ++i)
      dw[i] = dw[i - 1] * irs[i + 1] * rs[i + 3];
    dw[lv - 2] = dw[lv - 3] * irs[lv - 1];
    idw[0] = irs[3];
    for (int i = 1; i < lv - 2; ++i)
      idw[i] = idw[i - 1] * rs[i + 1] * irs[i + 3];
    idw[lv - 2] = idw[lv - 3] * rs[lv - 1];
  }

  void trans(int l, vector<LMod> &A) const {
    const size_t n = size_t(1) << l, nh = (n >> 1);
    const LMod one = rs[0], imag = rs[2];
    if (l & 1)
      for (size_t i = 0; i < nh; ++i) {
        LMod a = A[i], b = A[i + nh];
        A[i] = a + b;
        A[i + nh] = a - b; // < 2.1436 * 2^{60}, < 5.0718 * 2^{60}
      }
    for (int e = l & ~1; e >= 2; e -= 2) {
      const size_t m = size_t(1) << e, m4 = m >> 2;
      LMod w2 = one;
      for (size_t i = 0; i < n; i += m) {
        const LMod w1 = w2 * w2, w3 = w1 * w2; // < 1.0718 * 2^{60}
        // assume A[*] < 9.65 * 2^{60}
        for (size_t j = i; j < i + m4; ++j) {
          LMod a0 = A[j + m4 * 0] * one,
               a1 = A[j + m4 * 1] * w2; // < 1.65 * 2^{60}
          LMod a2 = A[j + m4 * 2] * w1,
               a3 = A[j + m4 * 3] * w3;        // < 1.65 * 2^{60}
          LMod t02p = a0 + a2, t13p = a1 + a3; // < 3.3 * 2^{60}, < 3.3 * 2^{60}
          LMod t02m = a0 - a2,
               t13m = (a1 - a3) * imag; // < 5.65 * 2^{60}, < 1.38 * 2^{60}
          A[j + m4 * 0] = t02p + t13p;
          A[j + m4 * 1] = t02p - t13p; // < 6.6 * 2^{60}, < 7.3 * 2^{60}
          A[j + m4 * 2] = t02m + t13m;
          A[j + m4 * 3] = t02m - t13m; // < 7.03 * 2^{60}, < 9.65 * 2^{60}
        }
        w2 *= dw[__builtin_ctzll(~(i >> e))];
      }
    }
  }

  void itrans(int l, vector<LMod> &A) const {
    const size_t n = size_t(1) << l, nh = (n >> 1);
    const LMod one = rs[0], imag = irs[2];
    for (int e = 2; e <= l; e += 2) {
      const size_t m = size_t(1) << e, m4 = m >> 2;
      LMod w2 = one;
      for (size_t i = 0; i < n; i += m) {
        const LMod w1 = w2 * w2, w3 = w1 * w2; // < 1.0718 * 2^{60}
        // assume A[*] < 1.65 * 2^{60}
        for (size_t j = i; j < i + m4; ++j) {
          LMod a0 = A[j + m4 * 0], a1 = A[j + m4 * 1];
          LMod a2 = A[j + m4 * 2], a3 = A[j + m4 * 3];
          LMod t01p = a0 + a1, t23p = a2 + a3; // 3.3 * 2^{60}, 3.3 * 2^{60}
          LMod t01m = a0 - a1,
               t23m = (a2 - a3) * imag; // 5.65 * 2^{60}, 1.38 * 2^{60}
          A[j + m4 * 0] = (t01p + t23p) * one;
          A[j + m4 * 2] = (t01p - t23p) * w1; // 1.45 * 2^{60}, 1.49 * 2^{60}
          A[j + m4 * 1] = (t01m + t23m) * w2;
          A[j + m4 * 3] = (t01m - t23m) * w3; // 1.48 * 2^{60}, 1.65 * 2^{60}
        }
        w2 *= idw[__builtin_ctzll(~(i >> e))];
      }
    }
    if (l & 1)
      for (size_t i = 0; i < nh; ++i) {
        LMod a = A[i], b = A[i + nh];
        A[i] = a + b;
        A[i + nh] = a - b; // < 3.3 * 2^{60}, 5.65 * 2^{60}
      }
  }

  template <typename T>
  vector<u64> convolve(const vector<T> &f, const vector<T> &g) const {
    // assume f[i], g[i] < 2^{60}
    const size_t s = f.size() + g.size() - 1;
    const int l = __lg(2 * s - 1);
    assert(l <= lv);
    const size_t sz = size_t(1) << l;
    vector<LMod> A(sz);
    for (size_t i = 0; i < f.size(); ++i)
      A[i] = LMod(f[i]); // < 1.0718 * 2^{60}
    trans(l, A);         // < 9.65 * 2^{60}
    const LMod inv = itwopows[l];
    if (&f == &g) {
      for (size_t i = 0; i < sz; ++i)
        (A[i] *= A[i]) *= inv; // < 1.46 * 2^{60}
    } else {
      vector<LMod> B(sz);
      for (size_t i = 0; i < g.size(); ++i)
        B[i] = LMod(g[i]);
      trans(l, B);
      for (size_t i = 0; i < sz; ++i)
        (A[i] *= B[i]) *= inv; // < 1.46 * 2^{60}
    }
    itrans(l, A); // < 5.65 * 2^{60}
    vector<u64> ret(s);
    for (size_t i = 0; i < s; ++i)
      ret[i] = A[i].get();
    return ret;
  }

  template <typename T>
  vector<u64> convolve_binom(const vector<T> &f, const vector<T> &g, int p,
                             const std::vector<int> &vs) const {
    // assume f[i], g[i] < 2^{60}
    const size_t s = f.size() + g.size() - 1;
    const int l = __lg(2 * s - 1);
    assert(l <= lv);
    const size_t sz = size_t(1) << l;
    vector<LMod> A(sz);
    const LMod lp = LMod(p), ilp = lp.inv(), one = rs[0];
    int e = 0;
    LMod ippow = one;
    for (size_t i = 0; i < f.size(); ++i) {
      for (; e < vs[i]; ++e)
        ippow *= ilp;
      A[i] = LMod(f[i]) * ippow; // < 1.0718 * 2^{60}
    }
    trans(l, A); // < 9.65 * 2^{60}
    const LMod inv = itwopows[l];
    if (&f == &g) {
      for (size_t i = 0; i < sz; ++i)
        (A[i] *= A[i]) *= inv; // < 1.46 * 2^{60}
    } else {
      vector<LMod> B(sz);
      e = 0;
      ippow = one;
      for (size_t i = 0; i < g.size(); ++i) {
        for (; e < vs[i]; ++e)
          ippow *= ilp;
        B[i] = LMod(g[i]) * ippow;
      }
      trans(l, B);
      for (size_t i = 0; i < sz; ++i)
        (A[i] *= B[i]) *= inv; // < 1.46 * 2^{60}
    }
    itrans(l, A); // < 5.65 * 2^{60}
    e = 0;
    LMod ppow = one;
    vector<u64> ret(s);
    for (size_t i = 0; i < s; ++i) {
      for (; e < vs[i]; ++e)
        ppow *= lp;
      ret[i] = (A[i] * ppow).get();
    }
    return ret;
  }

  LMod dw[lv - 1], idw[lv - 1], rs[lv + 1], irs[lv + 1], itwopows[lv + 1];
};

constexpr auto ntt64_1 = ntt<1148418935771627521, 19>();
constexpr auto ntt64_2 = ntt<1151514404601200641, 19>();

vector<pair<int, int>> factors(int n) {
  vector<pair<int, int>> ret;
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

#define getchar getchar_unlocked
#define putchar putchar_unlocked

int get_int() {
  int c, n;
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

struct ExactDiv {
  ExactDiv() {}
  constexpr ExactDiv(u64 n)
      : mod(n), shift(__builtin_ctzll(n)), mask((u64(1) << shift) - 1),
        inv(mul_inv(n >> shift)) {}
  constexpr friend u64 operator/(u64 n, ExactDiv d) {
    return (n >> d.shift) * d.inv;
  };
  constexpr bool divide(u64 n) const {
    if (n & mask)
      return false;
    return !__builtin_mul_overflow_p((n >> shift) * inv, mod, u64(0));
  }
  u64 mod;
  size_t shift;
  u64 mask, inv;
};

struct FastDiv {
  FastDiv() {}
  constexpr FastDiv(u64 n)
      : m(n), s((n == 1) ? 0 : 64 + __lg(n - 1)),
        x(((u128(1) << s) + n - 1) / n) {}
  constexpr friend u64 operator/(u64 n, FastDiv d) {
    return u128(n) * d.x >> d.s;
  }
  constexpr friend u64 operator%(u64 n, FastDiv d) { return n - n / d * d.m; }
  u64 m, s, x;
};

struct BinConvPre {
  BinConvPre(int n, int p, int pe) : facts(n + 1), ifacts(n + 1), vs(n + 1) {
    const ExactDiv ediv(p);
    const FastDiv fdiv(pe);
    facts[0] = 1;
    vs[0] = 0;
    for (int i = 1; i <= n; ++i) {
      int j = i, e = 0;
      while (ediv.divide(j))
        j = j / ediv, e += 1;
      facts[i] = u64(facts[i - 1]) * j % fdiv;
      vs[i] = vs[i - 1] + e;
    }
    int f = ifacts[n] = mod_inv(facts[n], pe);
    for (int i = n - 1; i >= 0; --i) {
      int j = i + 1;
      while (ediv.divide(j))
        j = j / ediv;
      ifacts[i] = f = u64(f) * j % fdiv;
    }
  }
  std::vector<int> facts;
  std::vector<int> ifacts;
  std::vector<int> vs;
};

int main() {
  int n, m, mod;
  while (~scanf("%d %d %d", &n, &m, &mod)) {
    vector<int> a(n + 1);
    for (int i = 0; i <= n; ++i)
      a[i] = get_int();
    vector<int> b(m + 1);
    for (int i = 0; i <= m; ++i)
      b[i] = get_int();

    static constexpr u64 mod1 = ntt64_1.modulus();
    static constexpr u64 mod2 = ntt64_2.modulus();
    static constexpr u64 inv = mod_inv(mod1, mod2);
    static_assert(mod1 < mod2, "");

    int s = n + m + 1;
    vector<int> ans(s);
    int prod = 1;

    for (auto &&pp : factors(mod)) {
      int p = pp.first, e = pp.second;
      int pe = p;
      for (int k = 1; k < e; ++k)
        pe *= p;
      auto pre = BinConvPre(s, p, pe);
      auto fpe = FastDiv(pe);
      vector<int> af(a);
      for (int i = 0; i <= n; ++i)
        af[i] = u64(af[i]) * pre.ifacts[i] % fpe;
      vector<int> bf(b);
      for (int i = 0; i <= m; ++i)
        bf[i] = u64(bf[i]) * pre.ifacts[i] % fpe;
      auto res1 = ntt64_1.convolve_binom(af, bf, p, pre.vs);
      auto res2 = ntt64_2.convolve_binom(af, bf, p, pre.vs);
      const int m1 = mod1 % fpe;
      const int iprod = mod_inv(prod, pe);
      for (int i = 0; i < s; ++i) {
        // crt1
        u64 r1 = res1[i], r2 = res2[i];
        int k = u64(u128(r2 + mod2 - r1) * inv % mod2) % fpe;
        int r = u64((r1 + u64(m1) * k) % fpe) * pre.facts[i] % fpe;
        // crt2
        int k2 = u64(r + pe - ans[i] % fpe) * iprod % fpe;
        ans[i] += k2 * prod;
      }
      prod *= pe;
    }
    printf("%d", ans[0]);
    for (int i = 1; i < s; ++i) {
      printf(" %d", ans[i]);
    }
    puts("");
  }
}
// https://loj.ac/p/174