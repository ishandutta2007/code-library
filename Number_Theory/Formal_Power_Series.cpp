#line 1 "main.cpp"
#define PROBLEM "https://judge.yosupo.jp/problem/pow_of_formal_power_series"

#line 1 "lib/mod/modint.hpp"

#include <cassert>
#include <numeric>
#include <type_traits>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#line 1 "lib/ac-library/atcoder/internal_math.hpp"

#include <utility>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace atcoder {

namespace internal {

// @param m `1 <= m`
// @return x mod m
constexpr long long safe_mod(long long x, long long m) {
  x %= m;
  if (x < 0)
    x += m;
  return x;
}

// Fast modular multiplication by barrett reduction
// Reference: https://en.wikipedia.org/wiki/Barrett_reduction
// NOTE: reconsider after Ice Lake
struct barrett {
  unsigned int _m;
  unsigned long long im;

  // @param m `1 <= m < 2^31`
  explicit barrett(unsigned int m)
      : _m(m), im((unsigned long long)(-1) / m + 1) {}

  // @return m
  unsigned int umod() const { return _m; }

  // @param a `0 <= a < m`
  // @param b `0 <= b < m`
  // @return `a * b % m`
  unsigned int mul(unsigned int a, unsigned int b) const {
    // [1] m = 1
    // a = b = im = 0, so okay

    // [2] m >= 2
    // im = ceil(2^64 / m)
    // -> im * m = 2^64 + r (0 <= r < m)
    // let z = a*b = c*m + d (0 <= c, d < m)
    // a*b * im = (c*m + d) * im = c*(im*m) + d*im = c*2^64 + c*r + d*im
    // c*r + d*im < m * m + m * im < m * m + 2^64 + m <= 2^64 + m * (m + 1) <
    // 2^64 * 2
    // ((ab * im) >> 64) == c or c + 1
    unsigned long long z = a;
    z *= b;
#ifdef _MSC_VER
    unsigned long long x;
    _umul128(z, im, &x);
#else
    unsigned long long x =
        (unsigned long long)(((unsigned __int128)(z)*im) >> 64);
#endif
    unsigned int v = (unsigned int)(z - x * _m);
    if (_m <= v)
      v += _m;
    return v;
  }
};

// @param n `0 <= n`
// @param m `1 <= m`
// @return `(x ** n) % m`
constexpr long long pow_mod_constexpr(long long x, long long n, int m) {
  if (m == 1)
    return 0;
  unsigned int _m = (unsigned int)(m);
  unsigned long long r = 1;
  unsigned long long y = safe_mod(x, m);
  while (n) {
    if (n & 1)
      r = (r * y) % _m;
    y = (y * y) % _m;
    n >>= 1;
  }
  return r;
}

// Reference:
// M. Forisek and J. Jancina,
// Fast Primality Testing for Integers That Fit into a Machine Word
// @param n `0 <= n`
constexpr bool is_prime_constexpr(int n) {
  if (n <= 1)
    return false;
  if (n == 2 || n == 7 || n == 61)
    return true;
  if (n % 2 == 0)
    return false;
  long long d = n - 1;
  while (d % 2 == 0)
    d /= 2;
  constexpr long long bases[3] = {2, 7, 61};
  for (long long a : bases) {
    long long t = d;
    long long y = pow_mod_constexpr(a, t, n);
    while (t != n - 1 && y != 1 && y != n - 1) {
      y = y * y % n;
      t <<= 1;
    }
    if (y != n - 1 && t % 2 == 0) {
      return false;
    }
  }
  return true;
}
template <int n> constexpr bool is_prime = is_prime_constexpr(n);

// @param b `1 <= b`
// @return pair(g, x) s.t. g = gcd(a, b), xa = g (mod b), 0 <= x < b/g
constexpr std::pair<long long, long long> inv_gcd(long long a, long long b) {
  a = safe_mod(a, b);
  if (a == 0)
    return {b, 0};

  // Contracts:
  // [1] s - m0 * a = 0 (mod b)
  // [2] t - m1 * a = 0 (mod b)
  // [3] s * |m1| + t * |m0| <= b
  long long s = b, t = a;
  long long m0 = 0, m1 = 1;

  while (t) {
    long long u = s / t;
    s -= t * u;
    m0 -= m1 * u; // |m1 * u| <= |m1| * s <= b

    // [3]:
    // (s - t * u) * |m1| + t * |m0 - m1 * u|
    // <= s * |m1| - t * u * |m1| + t * (|m0| + |m1| * u)
    // = s * |m1| + t * |m0| <= b

    auto tmp = s;
    s = t;
    t = tmp;
    tmp = m0;
    m0 = m1;
    m1 = tmp;
  }
  // by [3]: |m0| <= b/g
  // by g != b: |m0| < b/g
  if (m0 < 0)
    m0 += b / s;
  return {s, m0};
}

// Compile time primitive root
// @param m must be prime
// @return primitive root (and minimum in now)
constexpr int primitive_root_constexpr(int m) {
  if (m == 2)
    return 1;
  if (m == 167772161)
    return 3;
  if (m == 469762049)
    return 3;
  if (m == 754974721)
    return 11;
  if (m == 998244353)
    return 3;
  int divs[20] = {};
  divs[0] = 2;
  int cnt = 1;
  int x = (m - 1) / 2;
  while (x % 2 == 0)
    x /= 2;
  for (int i = 3; (long long)(i)*i <= x; i += 2) {
    if (x % i == 0) {
      divs[cnt++] = i;
      while (x % i == 0) {
        x /= i;
      }
    }
  }
  if (x > 1) {
    divs[cnt++] = x;
  }
  for (int g = 2;; g++) {
    bool ok = true;
    for (int i = 0; i < cnt; i++) {
      if (pow_mod_constexpr(g, (m - 1) / divs[i], m) == 1) {
        ok = false;
        break;
      }
    }
    if (ok)
      return g;
  }
}
template <int m> constexpr int primitive_root = primitive_root_constexpr(m);

// @param n `n < 2^32`
// @param m `1 <= m < 2^32`
// @return sum_{i=0}^{n-1} floor((ai + b) / m) (mod 2^64)
unsigned long long floor_sum_unsigned(unsigned long long n,
                                      unsigned long long m,
                                      unsigned long long a,
                                      unsigned long long b) {
  unsigned long long ans = 0;
  while (true) {
    if (a >= m) {
      ans += n * (n - 1) / 2 * (a / m);
      a %= m;
    }
    if (b >= m) {
      ans += n * (b / m);
      b %= m;
    }

    unsigned long long y_max = a * n + b;
    if (y_max < m)
      break;
    // y_max < m * (n + 1)
    // floor(y_max / m) <= n
    n = (unsigned long long)(y_max / m);
    b = (unsigned long long)(y_max % m);
    std::swap(m, a);
  }
  return ans;
}

} // namespace internal

} // namespace atcoder

#line 1 "lib/ac-library/atcoder/internal_type_traits.hpp"

#line 7 "lib/ac-library/atcoder/internal_type_traits.hpp"

namespace atcoder {

namespace internal {

#ifndef _MSC_VER
template <class T>
using is_signed_int128 =
    typename std::conditional<std::is_same<T, __int128_t>::value ||
                                  std::is_same<T, __int128>::value,
                              std::true_type, std::false_type>::type;

template <class T>
using is_unsigned_int128 =
    typename std::conditional<std::is_same<T, __uint128_t>::value ||
                                  std::is_same<T, unsigned __int128>::value,
                              std::true_type, std::false_type>::type;

template <class T>
using make_unsigned_int128 =
    typename std::conditional<std::is_same<T, __int128_t>::value, __uint128_t,
                              unsigned __int128>;

template <class T>
using is_integral =
    typename std::conditional<std::is_integral<T>::value ||
                                  is_signed_int128<T>::value ||
                                  is_unsigned_int128<T>::value,
                              std::true_type, std::false_type>::type;

template <class T>
using is_signed_int =
    typename std::conditional<(is_integral<T>::value &&
                               std::is_signed<T>::value) ||
                                  is_signed_int128<T>::value,
                              std::true_type, std::false_type>::type;

template <class T>
using is_unsigned_int =
    typename std::conditional<(is_integral<T>::value &&
                               std::is_unsigned<T>::value) ||
                                  is_unsigned_int128<T>::value,
                              std::true_type, std::false_type>::type;

template <class T>
using to_unsigned = typename std::conditional<
    is_signed_int128<T>::value, make_unsigned_int128<T>,
    typename std::conditional<std::is_signed<T>::value, std::make_unsigned<T>,
                              std::common_type<T>>::type>::type;

#else

template <class T> using is_integral = typename std::is_integral<T>;

template <class T>
using is_signed_int =
    typename std::conditional<is_integral<T>::value && std::is_signed<T>::value,
                              std::true_type, std::false_type>::type;

template <class T>
using is_unsigned_int =
    typename std::conditional<is_integral<T>::value &&
                                  std::is_unsigned<T>::value,
                              std::true_type, std::false_type>::type;

template <class T>
using to_unsigned =
    typename std::conditional<is_signed_int<T>::value, std::make_unsigned<T>,
                              std::common_type<T>>::type;

#endif

template <class T>
using is_signed_int_t = std::enable_if_t<is_signed_int<T>::value>;

template <class T>
using is_unsigned_int_t = std::enable_if_t<is_unsigned_int<T>::value>;

template <class T> using to_unsigned_t = typename to_unsigned<T>::type;

} // namespace internal

} // namespace atcoder

#line 14 "lib/mod/modint.hpp"

namespace atcoder {

namespace internal {

struct modint_base {};
struct static_modint_base : modint_base {};

template <class T> using is_modint = std::is_base_of<modint_base, T>;
template <class T> using is_modint_t = std::enable_if_t<is_modint<T>::value>;

} // namespace internal

template <int m, std::enable_if_t<(1 <= m)> * = nullptr>
struct static_modint : internal::static_modint_base {
  using mint = static_modint;

public:
  static constexpr int mod() { return m; }
  static constexpr mint raw(int v) {
    mint x;
    x._v = v;
    return x;
  }

  constexpr static_modint() : _v(0) {}
  template <class T, internal::is_signed_int_t<T> * = nullptr>
  constexpr static_modint(T v)
      : _v() {
    long long x = (long long)(v % (long long)(umod()));
    if (x < 0)
      x += umod();
    _v = (unsigned int)(x);
  }
  template <class T, internal::is_unsigned_int_t<T> * = nullptr>
  constexpr static_modint(T v)
      : _v() {
    _v = (unsigned int)(v % umod());
  }

  constexpr unsigned int val() const { return _v; }

  constexpr mint &operator++() {
    _v++;
    if (_v == umod())
      _v = 0;
    return *this;
  }
  constexpr mint &operator--() {
    if (_v == 0)
      _v = umod();
    _v--;
    return *this;
  }
  constexpr mint operator++(int) {
    mint result = *this;
    ++*this;
    return result;
  }
  constexpr mint operator--(int) {
    mint result = *this;
    --*this;
    return result;
  }

  constexpr mint &operator+=(const mint &rhs) {
    _v += rhs._v;
    if (_v >= umod())
      _v -= umod();
    return *this;
  }
  constexpr mint &operator-=(const mint &rhs) {
    _v -= rhs._v;
    if (_v >= umod())
      _v += umod();
    return *this;
  }
  constexpr mint &operator*=(const mint &rhs) {
    unsigned long long z = _v;
    z *= rhs._v;
    _v = (unsigned int)(z % umod());
    return *this;
  }
  constexpr mint &operator/=(const mint &rhs) {
    return *this = *this * rhs.inv();
  }

  constexpr mint operator+() const { return *this; }
  constexpr mint operator-() const { return mint() - *this; }

  constexpr mint pow(long long n) const {
    assert(0 <= n);
    mint x = *this, r = 1;
    while (n) {
      if (n & 1)
        r *= x;
      x *= x;
      n >>= 1;
    }
    return r;
  }
  constexpr mint inv() const {
    if (prime) {
      assert(_v);
      return pow(umod() - 2);
    } else {
      auto eg = internal::inv_gcd(_v, m);
      assert(eg.first == 1);
      return eg.second;
    }
  }

  constexpr friend mint operator+(const mint &lhs, const mint &rhs) {
    return mint(lhs) += rhs;
  }
  constexpr friend mint operator-(const mint &lhs, const mint &rhs) {
    return mint(lhs) -= rhs;
  }
  constexpr friend mint operator*(const mint &lhs, const mint &rhs) {
    return mint(lhs) *= rhs;
  }
  constexpr friend mint operator/(const mint &lhs, const mint &rhs) {
    return mint(lhs) /= rhs;
  }
  constexpr friend bool operator==(const mint &lhs, const mint &rhs) {
    return lhs._v == rhs._v;
  }
  constexpr friend bool operator!=(const mint &lhs, const mint &rhs) {
    return lhs._v != rhs._v;
  }

private:
  unsigned int _v;
  static constexpr unsigned int umod() { return m; }
  static constexpr bool prime = internal::is_prime<m>;
};

template <int id> struct dynamic_modint : internal::modint_base {
  using mint = dynamic_modint;

public:
  static int mod() { return (int)(bt.umod()); }
  static void set_mod(int m) {
    assert(1 <= m);
    bt = internal::barrett(m);
  }
  static mint raw(int v) {
    mint x;
    x._v = v;
    return x;
  }

  dynamic_modint() : _v(0) {}
  template <class T, internal::is_signed_int_t<T> * = nullptr>
  dynamic_modint(T v) {
    long long x = (long long)(v % (long long)(mod()));
    if (x < 0)
      x += mod();
    _v = (unsigned int)(x);
  }
  template <class T, internal::is_unsigned_int_t<T> * = nullptr>
  dynamic_modint(T v) {
    _v = (unsigned int)(v % mod());
  }

  unsigned int val() const { return _v; }

  mint &operator++() {
    _v++;
    if (_v == umod())
      _v = 0;
    return *this;
  }
  mint &operator--() {
    if (_v == 0)
      _v = umod();
    _v--;
    return *this;
  }
  mint operator++(int) {
    mint result = *this;
    ++*this;
    return result;
  }
  mint operator--(int) {
    mint result = *this;
    --*this;
    return result;
  }

  mint &operator+=(const mint &rhs) {
    _v += rhs._v;
    if (_v >= umod())
      _v -= umod();
    return *this;
  }
  mint &operator-=(const mint &rhs) {
    _v += mod() - rhs._v;
    if (_v >= umod())
      _v -= umod();
    return *this;
  }
  mint &operator*=(const mint &rhs) {
    _v = bt.mul(_v, rhs._v);
    return *this;
  }
  mint &operator/=(const mint &rhs) { return *this = *this * rhs.inv(); }

  mint operator+() const { return *this; }
  mint operator-() const { return mint() - *this; }

  mint pow(long long n) const {
    assert(0 <= n);
    mint x = *this, r = 1;
    while (n) {
      if (n & 1)
        r *= x;
      x *= x;
      n >>= 1;
    }
    return r;
  }
  mint inv() const {
    auto eg = internal::inv_gcd(_v, mod());
    assert(eg.first == 1);
    return eg.second;
  }

  friend mint operator+(const mint &lhs, const mint &rhs) {
    return mint(lhs) += rhs;
  }
  friend mint operator-(const mint &lhs, const mint &rhs) {
    return mint(lhs) -= rhs;
  }
  friend mint operator*(const mint &lhs, const mint &rhs) {
    return mint(lhs) *= rhs;
  }
  friend mint operator/(const mint &lhs, const mint &rhs) {
    return mint(lhs) /= rhs;
  }
  friend bool operator==(const mint &lhs, const mint &rhs) {
    return lhs._v == rhs._v;
  }
  friend bool operator!=(const mint &lhs, const mint &rhs) {
    return lhs._v != rhs._v;
  }

private:
  unsigned int _v;
  static internal::barrett bt;
  static unsigned int umod() { return bt.umod(); }
};
template <int id> internal::barrett dynamic_modint<id>::bt(998244353);

using modint998244353 = static_modint<998244353>;
using modint1000000007 = static_modint<1000000007>;
using modint = dynamic_modint<-1>;

namespace internal {

template <class T>
using is_static_modint = std::is_base_of<internal::static_modint_base, T>;

template <class T>
using is_static_modint_t = std::enable_if_t<is_static_modint<T>::value>;

template <class> struct is_dynamic_modint : public std::false_type {};
template <int id>
struct is_dynamic_modint<dynamic_modint<id>> : public std::true_type {};

template <class T>
using is_dynamic_modint_t = std::enable_if_t<is_dynamic_modint<T>::value>;

} // namespace internal

} // namespace atcoder

#line 2 "lib/prelude.hpp"
#ifndef LOCAL
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
#endif
#include <bits/stdc++.h>
using namespace std;
using ll = long long;
#define rep2(i, m, n) for (auto i = (m); i < (n); i++)
#define rep(i, n) rep2(i, 0, n)
#define repr2(i, m, n) for (auto i = (n); i-- > (m);)
#define repr(i, n) repr2(i, 0, n)
#define all(x) begin(x), end(x)
template <class T> auto ndvec(size_t n, T &&x) {
  return vector(n, forward<T>(x));
}
template <class... Ts> auto ndvec(size_t n, Ts &&... xs) {
  return vector(n, ndvec(forward<Ts>(xs)...));
}
#line 3 "lib/io.hpp"

template <size_t BUF_SIZE = 1 << 26> class stdin_reader {
public:
  stdin_reader() { buf[fread(buf, 1, sizeof(buf), stdin)] = 0; }

  template <class T> enable_if_t<is_integral_v<T>> read(T &x) {
    skip();
    [[maybe_unused]] bool neg = false;
    if
      constexpr(is_signed_v<T>)neg = *p == '-' ? (p++, true) : false;
    x = 0;
    while (*p > ' ')
      x = x * 10 + (*p++ & 0x0F);
    if
      constexpr(is_signed_v<T>)x = neg ? -x : x;
  }
  template <class T> void_t<decltype(&T::val)> read(T &x) {
    x = T((unsigned)(*this));
  }
  void read(char *q) {
    skip();
    char *p0 = p;
    while (*p > ' ')
      p++;
    copy(p0, p, q);
  }
  template <size_t N> void read(char(&s)[N]) { read(s); }
  void read(string &s) {
    skip();
    char *p0 = p;
    while (*p > ' ')
      p++;
    s.assign(p0, p);
  }
  template <class T, class U> void read(pair<T, U> &x) {
    read(x.first), read(x.second);
  }
  template <class... Ts> void read(tuple<Ts...> &x) {
    read_tuple(x, make_index_sequence<sizeof...(Ts)>{});
  }
  template <class T, size_t N> void read(T(&a)[N]) {
    for (auto &e : a)
      read(e);
  }

  template <class T> operator T() {
    T x;
    return read(x), x;
  }
  template <class... Ts> void operator()(Ts &... xs) { (read(xs), ...); }
  int operator--() { return (int)*this - 1; }
  template <class T> void vec(vector<T> &v, int n) {
    v.resize(n);
    for (auto &e : v)
      read(e);
  }
  template <class T> vector<T> vec(int n) {
    vector<T> v(n);
    return vec(v, n), v;
  }

private:
  char buf[BUF_SIZE], *p = buf;
  void skip() {
    while (*p <= ' ')
      p++;
  }
  template <class T, size_t... Is>
  void read_tuple(T &x, index_sequence<Is...>) {
    (*this)(get<Is>(x)...);
  }
};

template <size_t BUF_SIZE = 1 << 26> class stdout_writer {
public:
  ~stdout_writer() { flush(); }
  void flush() { fwrite(buf, 1, p - buf, stdout), p = buf; }
  void write_char(char c) { *p++ = c; }
  void write(char c) { write_char(c); }
  template <class T> enable_if_t<is_integral_v<T>> write(T x) {
    if (!x)
      return write_char('0');
    if
      constexpr(is_signed_v<T>)if (x < 0) write_char('-'), x = -x;
    static char tmp[16];
    char *q = end(tmp);
    while (x >= 10000)
      memcpy(q -= 4, four_digits.data + x % 10000 * 4, 4), x /= 10000;
    if (x < 10)
      write_char('0' + x);
    else if (x < 100)
      write_char('0' + (uint8_t)x / 10), write_char('0' + (uint8_t)x % 10);
    else if (x < 1000)
      memcpy(p, four_digits.data + x * 4 + 1, 3), p += 3;
    else
      memcpy(p, four_digits.data + x * 4, 4), p += 4;
    memcpy(p, q, end(tmp) - q), p += end(tmp) - q;
  }
  template <class T> void_t<decltype(&T::val)> write(T x) { write(x.val()); }
  void write(double x) {
    ll integer = x;
    write(integer), write_char('.'), write((int)((x - integer) * 1000000000));
  }
  void write(const char *s) {
    while (*s)
      *p++ = *s++;
  }
  void write(const string &s) { memcpy(p, s.c_str(), s.size()), p += s.size(); }
  template <class T, class U> void write(const pair<T, U> &x) {
    write(x.first), write_char(' '), write(x.second);
  }
  template <class... Ts> void write(const tuple<Ts...> &x) {
    write_tuple(x, make_index_sequence<sizeof...(Ts)>{});
  }
  template <class... Ts> void write(const Ts &... xs) {
    ((write(xs), write_char(' ')), ...), --p;
  }
  template <class... Ts> void writeln(const Ts &... xs) {
    write(xs...), write_char('\n');
  }

  template <class... Ts> void operator()(const Ts &... xs) { writeln(xs...); }
  template <class It> void iter(It first, It last, char sep = ' ') {
    if (first == last)
      write_char('\n');
    else {
      while (first != last)
        write(*first++), write_char(sep);
      p[-1] = '\n';
    }
  }

#define INSTANT(s)                                                             \
  void s() { writeln("s"); }
  INSTANT(No)
  INSTANT(NO)
  INSTANT(Aoki)
  INSTANT(possible)
  INSTANT(Possible)
  INSTANT(POSSIBLE)
  INSTANT(impossible)
  INSTANT(Impossible)
  INSTANT(IMPOSSIBLE)
#undef INSTANT
  void Yes(bool b = true) { writeln(b ? "Yes" : "No"); }
  void YES(bool b = true) { writeln(b ? "YES" : "NO"); }
  void Takahashi(bool b = true) { writeln(b ? "Takahashi" : "Aoki"); }

private:
  char buf[BUF_SIZE], *p = buf;
  template <class T, size_t... Is>
  void write_tuple(const T &x, index_sequence<Is...>) {
    ((write(get<Is>(x)), write_char(' ')), ...), --p;
  }
  struct four_digits {
    char data[40000];
    constexpr four_digits() : data() {
      for (int i = 0; i < 10000; i++)
        for (int n = i, j = 4; j--;)
          data[i * 4 + j] = n % 10 + '0', n /= 10;
    }
  } static constexpr four_digits{};
};

static stdin_reader<> in;
static stdout_writer<> out;
#line 3 "lib/mod/inv.hpp"

template <class T> T inverse(int n) {
  static vector<T> inv = {T(0), T(1)};
  if (inv.size() <= n) {
    int l = inv.size();
    inv.resize(n + 1);
    rep2(i, l, n + 1) inv[i] = -inv[T::mod() % i] * (T::mod() / i);
  }
  return inv[n];
}
#line 3 "lib/bit/ctz.hpp"

template <class T> __attribute__((target("bmi"))) int ctz(T x) {
  if (!x)
    return sizeof(T) * 8;
  if
    constexpr(sizeof(T) <= sizeof(unsigned)) {
      return __builtin_ctz((unsigned)x);
    }
  else if
    constexpr(sizeof(T) <= sizeof(unsigned long long)) {
      return __builtin_ctzll((unsigned long long)x);
    }
  else if
    constexpr(sizeof(T) <= sizeof(unsigned long long) * 2) {
      unsigned long long y = x;
      return y ? ctz(y)
               : sizeof(y) * 8 + ctz((unsigned long long)(x >> sizeof(y) * 8));
    }
}
#line 3 "lib/util/seed.hpp"

auto seed() {
#if defined(LOCAL) || defined(FIX_SEED)
  return 314169265258979;
#endif
  return chrono::steady_clock::now().time_since_epoch().count();
}
#line 3 "lib/util/rand.hpp"

uint32_t rand32() {
  static uint32_t x = seed();
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  return x;
}
#line 4 "lib/mod/sqrt.hpp"

template <class T> optional<T> mod_sqrt(T a) {
  // Tonelli-Shanks
  if (T::mod() <= 2)
    return a;
  if (a == T(0))
    return T(0);
  if (a.pow((T::mod() - 1) / 2) == -1)
    return nullopt;
  int s = ctz(T::mod() - 1);
  int q = (T::mod() - 1) >> s;
  T x = a.pow((q + 1) / 2);
  T b = rand32();
  while (b.pow((T::mod() - 1) / 2) != -1)
    b = rand32();
  b = b.pow(q);
  T ia = a.inv();
  s -= 2;
  for (T e = ia * x * x; e != 1; b *= b, s--) {
    if (e.pow(1 << s) != 1)
      x *= b, e = ia * x * x;
  }
  return x;
}
#line 3 "lib/arith/sat.hpp"

template <class T> T sat_add(T a, T b) {
  T res;
  return __builtin_add_overflow(a, b, &res)
             ? (a < 0 ? numeric_limits<T>::min() : numeric_limits<T>::max())
             : res;
}
template <class T> T sat_sub(T a, T b) {
  T res;
  return __builtin_sub_overflow(a, b, &res)
             ? (a < 0 ? numeric_limits<T>::min() : numeric_limits<T>::max())
             : res;
}
template <class T> T sat_mul(T a, T b) {
  T res;
  return __builtin_mul_overflow(a, b, &res)
             ? ((a < 0) == (b < 0) ? numeric_limits<T>::max()
                                   : numeric_limits<T>::min())
             : res;
}
#line 3 "lib/types.hpp"

template <class It> using val_t = typename iterator_traits<It>::value_type;
#line 4 "lib/ps/fft.hpp"

#line 1 "lib/ac-library/atcoder/internal_bit.hpp"

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace atcoder {

namespace internal {

// @param n `0 <= n`
// @return minimum non-negative `x` s.t. `n <= 2**x`
int ceil_pow2(int n) {
  int x = 0;
  while ((1U << x) < (unsigned int)(n))
    x++;
  return x;
}

// @param n `1 <= n`
// @return minimum non-negative `x` s.t. `(n & (1 << x)) != 0`
constexpr int bsf_constexpr(unsigned int n) {
  int x = 0;
  while (!(n & (1 << x)))
    x++;
  return x;
}

// @param n `1 <= n`
// @return minimum non-negative `x` s.t. `(n & (1 << x)) != 0`
int bsf(unsigned int n) {
#ifdef _MSC_VER
  unsigned long index;
  _BitScanForward(&index, n);
  return index;
#else
  return __builtin_ctz(n);
#endif
}

} // namespace internal

} // namespace atcoder

#line 8 "lib/ps/fft.hpp"

#ifndef ATCODER_CONVOLUTION_HPP
#define ATCODER_CONVOLUTION_HPP

namespace atcoder {

namespace internal {

template <class mint, int g = internal::primitive_root<mint::mod()>,
          internal::is_static_modint_t<mint> * = nullptr>
struct fft_info {
  static constexpr int rank2 = bsf_constexpr(mint::mod() - 1);
  mint root[rank2 + 1];  // root[i]^(2^i) == 1
  mint iroot[rank2 + 1]; // root[i] * iroot[i] == 1

  mint rate2[std::max(0, rank2 - 2 + 1)];
  mint irate2[std::max(0, rank2 - 2 + 1)];

  mint rate3[std::max(0, rank2 - 3 + 1)];
  mint irate3[std::max(0, rank2 - 3 + 1)];

  constexpr fft_info() : root(), iroot(), rate2(), irate2(), rate3(), irate3() {
    root[rank2] = mint(g).pow((mint::mod() - 1) >> rank2);
    iroot[rank2] = root[rank2].inv();
    for (int i = rank2 - 1; i >= 0; i--) {
      root[i] = root[i + 1] * root[i + 1];
      iroot[i] = iroot[i + 1] * iroot[i + 1];
    }

    {
      mint prod = 1, iprod = 1;
      for (int i = 0; i <= rank2 - 2; i++) {
        rate2[i] = root[i + 2] * prod;
        irate2[i] = iroot[i + 2] * iprod;
        prod *= iroot[i + 2];
        iprod *= root[i + 2];
      }
    }
    {
      mint prod = 1, iprod = 1;
      for (int i = 0; i <= rank2 - 3; i++) {
        rate3[i] = root[i + 3] * prod;
        irate3[i] = iroot[i + 3] * iprod;
        prod *= iroot[i + 3];
        iprod *= root[i + 3];
      }
    }
  }
};

template <class It, internal::is_static_modint_t<val_t<It>> * = nullptr>
void butterfly(It a, It a_last) {
  using mint = val_t<It>;
  int n = a_last - a;
  int h = internal::ceil_pow2(n);

  static constexpr fft_info<mint> info;

  int len = 0; // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
  while (len < h) {
    if (h - len == 1) {
      int p = 1 << (h - len - 1);
      mint rot = 1;
      for (int s = 0; s < (1 << len); s++) {
        int offset = s << (h - len);
        for (int i = 0; i < p; i++) {
          auto l = a[i + offset];
          auto r = a[i + offset + p] * rot;
          a[i + offset] = l + r;
          a[i + offset + p] = l - r;
        }
        if (s + 1 != (1 << len))
          rot *= info.rate2[bsf(~(unsigned int)(s))];
      }
      len++;
    } else {
      // 4-base
      int p = 1 << (h - len - 2);
      mint rot = 1, imag = info.root[2];
      for (int s = 0; s < (1 << len); s++) {
        mint rot2 = rot * rot;
        mint rot3 = rot2 * rot;
        int offset = s << (h - len);
        for (int i = 0; i < p; i++) {
          auto mod2 = 1ULL * mint::mod() * mint::mod();
          auto a0 = 1ULL * a[i + offset].val();
          auto a1 = 1ULL * a[i + offset + p].val() * rot.val();
          auto a2 = 1ULL * a[i + offset + 2 * p].val() * rot2.val();
          auto a3 = 1ULL * a[i + offset + 3 * p].val() * rot3.val();
          auto a1na3imag = 1ULL * mint(a1 + mod2 - a3).val() * imag.val();
          auto na2 = mod2 - a2;
          a[i + offset] = a0 + a2 + a1 + a3;
          a[i + offset + 1 * p] = a0 + a2 + (2 * mod2 - (a1 + a3));
          a[i + offset + 2 * p] = a0 + na2 + a1na3imag;
          a[i + offset + 3 * p] = a0 + na2 + (mod2 - a1na3imag);
        }
        if (s + 1 != (1 << len))
          rot *= info.rate3[bsf(~(unsigned int)(s))];
      }
      len += 2;
    }
  }
}

template <class mint> void butterfly(std::vector<mint> &a) {
  butterfly(a.begin(), a.end());
}

template <class It, internal::is_static_modint_t<val_t<It>> * = nullptr>
void butterfly_inv(It a, It a_last) {
  using mint = val_t<It>;
  int n = a_last - a;
  int h = internal::ceil_pow2(n);

  static constexpr fft_info<mint> info;

  int len = h; // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
  while (len) {
    if (len == 1) {
      int p = 1 << (h - len);
      mint irot = 1;
      for (int s = 0; s < (1 << (len - 1)); s++) {
        int offset = s << (h - len + 1);
        for (int i = 0; i < p; i++) {
          auto l = a[i + offset];
          auto r = a[i + offset + p];
          a[i + offset] = l + r;
          a[i + offset + p] =
              (unsigned long long)(mint::mod() + l.val() - r.val()) *
              irot.val();
          ;
        }
        if (s + 1 != (1 << (len - 1)))
          irot *= info.irate2[bsf(~(unsigned int)(s))];
      }
      len--;
    } else {
      // 4-base
      int p = 1 << (h - len);
      mint irot = 1, iimag = info.iroot[2];
      for (int s = 0; s < (1 << (len - 2)); s++) {
        mint irot2 = irot * irot;
        mint irot3 = irot2 * irot;
        int offset = s << (h - len + 2);
        for (int i = 0; i < p; i++) {
          auto a0 = 1ULL * a[i + offset + 0 * p].val();
          auto a1 = 1ULL * a[i + offset + 1 * p].val();
          auto a2 = 1ULL * a[i + offset + 2 * p].val();
          auto a3 = 1ULL * a[i + offset + 3 * p].val();

          auto a2na3iimag =
              1ULL * mint((mint::mod() + a2 - a3) * iimag.val()).val();

          a[i + offset] = a0 + a1 + a2 + a3;
          a[i + offset + 1 * p] =
              (a0 + (mint::mod() - a1) + a2na3iimag) * irot.val();
          a[i + offset + 2 * p] =
              (a0 + a1 + (mint::mod() - a2) + (mint::mod() - a3)) * irot2.val();
          a[i + offset + 3 * p] =
              (a0 + (mint::mod() - a1) + (mint::mod() - a2na3iimag)) *
              irot3.val();
        }
        if (s + 1 != (1 << (len - 2)))
          irot *= info.irate3[bsf(~(unsigned int)(s))];
      }
      len -= 2;
    }
  }
}

template <class mint> void butterfly_inv(vector<mint> &a) {
  butterfly_inv(a.begin(), a.end());
}

template <class mint, internal::is_static_modint_t<mint> * = nullptr>
std::vector<mint> convolution_naive(const std::vector<mint> &a,
                                    const std::vector<mint> &b) {
  int n = int(a.size()), m = int(b.size());
  std::vector<mint> ans(n + m - 1);
  if (n < m) {
    for (int j = 0; j < m; j++) {
      for (int i = 0; i < n; i++) {
        ans[i + j] += a[i] * b[j];
      }
    }
  } else {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        ans[i + j] += a[i] * b[j];
      }
    }
  }
  return ans;
}

template <class mint, internal::is_static_modint_t<mint> * = nullptr>
std::vector<mint> convolution_fft(std::vector<mint> a, std::vector<mint> b) {
  int n = int(a.size()), m = int(b.size());
  int z = 1 << internal::ceil_pow2(n + m - 1);
  a.resize(z);
  internal::butterfly(a);
  b.resize(z);
  internal::butterfly(b);
  for (int i = 0; i < z; i++) {
    a[i] *= b[i];
  }
  internal::butterfly_inv(a);
  a.resize(n + m - 1);
  mint iz = mint(z).inv();
  for (int i = 0; i < n + m - 1; i++)
    a[i] *= iz;
  return a;
}

} // namespace internal

template <class mint, internal::is_static_modint_t<mint> * = nullptr>
std::vector<mint> convolution(std::vector<mint> &&a, std::vector<mint> &&b) {
  int n = int(a.size()), m = int(b.size());
  if (!n || !m)
    return {};
  if (std::min(n, m) <= 60)
    return convolution_naive(a, b);
  return internal::convolution_fft(a, b);
}

template <class mint, internal::is_static_modint_t<mint> * = nullptr>
std::vector<mint> convolution(const std::vector<mint> &a,
                              const std::vector<mint> &b) {
  int n = int(a.size()), m = int(b.size());
  if (!n || !m)
    return {};
  if (std::min(n, m) <= 60)
    return convolution_naive(a, b);
  return internal::convolution_fft(a, b);
}

template <unsigned int mod = 998244353, class T,
          std::enable_if_t<internal::is_integral<T>::value> * = nullptr>
std::vector<T> convolution(const std::vector<T> &a, const std::vector<T> &b) {
  int n = int(a.size()), m = int(b.size());
  if (!n || !m)
    return {};

  using mint = static_modint<mod>;
  std::vector<mint> a2(n), b2(m);
  for (int i = 0; i < n; i++) {
    a2[i] = mint(a[i]);
  }
  for (int i = 0; i < m; i++) {
    b2[i] = mint(b[i]);
  }
  auto c2 = convolution(move(a2), move(b2));
  std::vector<T> c(n + m - 1);
  for (int i = 0; i < n + m - 1; i++) {
    c[i] = c2[i].val();
  }
  return c;
}

std::vector<long long> convolution_ll(const std::vector<long long> &a,
                                      const std::vector<long long> &b) {
  int n = int(a.size()), m = int(b.size());
  if (!n || !m)
    return {};

  static constexpr unsigned long long MOD1 = 754974721; // 2^24
  static constexpr unsigned long long MOD2 = 167772161; // 2^25
  static constexpr unsigned long long MOD3 = 469762049; // 2^26
  static constexpr unsigned long long M2M3 = MOD2 * MOD3;
  static constexpr unsigned long long M1M3 = MOD1 * MOD3;
  static constexpr unsigned long long M1M2 = MOD1 * MOD2;
  static constexpr unsigned long long M1M2M3 = MOD1 * MOD2 * MOD3;

  static constexpr unsigned long long i1 =
      internal::inv_gcd(MOD2 * MOD3, MOD1).second;
  static constexpr unsigned long long i2 =
      internal::inv_gcd(MOD1 * MOD3, MOD2).second;
  static constexpr unsigned long long i3 =
      internal::inv_gcd(MOD1 * MOD2, MOD3).second;

  auto c1 = convolution<MOD1>(a, b);
  auto c2 = convolution<MOD2>(a, b);
  auto c3 = convolution<MOD3>(a, b);

  std::vector<long long> c(n + m - 1);
  for (int i = 0; i < n + m - 1; i++) {
    unsigned long long x = 0;
    x += (c1[i] * i1) % MOD1 * M2M3;
    x += (c2[i] * i2) % MOD2 * M1M3;
    x += (c3[i] * i3) % MOD3 * M1M2;
    // B = 2^63, -B <= x, r(real value) < B
    // (x, x - M, x - 2M, or x - 3M) = r (mod 2B)
    // r = c1[i] (mod MOD1)
    // focus on MOD1
    // r = x, x - M', x - 2M', x - 3M' (M' = M % 2^64) (mod 2B)
    // r = x,
    //     x - M' + (0 or 2B),
    //     x - 2M' + (0, 2B or 4B),
    //     x - 3M' + (0, 2B, 4B or 6B) (without mod!)
    // (r - x) = 0, (0)
    //           - M' + (0 or 2B), (1)
    //           -2M' + (0 or 2B or 4B), (2)
    //           -3M' + (0 or 2B or 4B or 6B) (3) (mod MOD1)
    // we checked that
    //   ((1) mod MOD1) mod 5 = 2
    //   ((2) mod MOD1) mod 5 = 3
    //   ((3) mod MOD1) mod 5 = 4
    long long diff =
        c1[i] - internal::safe_mod((long long)(x), (long long)(MOD1));
    if (diff < 0)
      diff += MOD1;
    static constexpr unsigned long long offset[5] = {0, 0, M1M2M3, 2 * M1M2M3,
                                                     3 * M1M2M3};
    x -= offset[diff % 5];
    c[i] = x;
  }

  return c;
}

} // namespace atcoder

#endif // ATCODER_CONVOLUTION_HPP

template <class It> void fft(It a, It a_last) {
  atcoder::internal::butterfly(a, a_last);
}

template <class T> void fft(vector<T> &a, int n = -1) {
  if (n != -1)
    a.resize(n);
  fft(all(a));
}

template <class It> void ifft(It a, It a_last) {
  atcoder::internal::butterfly_inv(a, a_last);
}

template <class T> void ifft(vector<T> &a) { ifft(all(a)); }

template <class T> void double_fft(vector<T> &a) {
  static constexpr atcoder::internal::fft_info<T> info{};
  int m = a.size();
  a.reserve(m * 2), a.insert(a.end(), all(a));
  ifft(a.begin() + m, a.end());
  T z = T(m).inv();
  T w = info.root[ctz(m * 2)];
  rep2(i, m, m * 2) a[i] *= z, z *= w;
  fft(a.begin() + m, a.end());
}
#line 4 "lib/ps/exp.hpp"

template <class T> vector<T> exp(const vector<T> &p, int deg = -1) {
  if (p.empty())
    return vector<T>{T(1)};
  assert(p[0] == 0);
  if (deg == -1)
    deg = p.size();
  int z = 1 << atcoder::internal::ceil_pow2(deg);
  vector<T> f = {1}, ffft = {1}, g = {1}, gfft = {1};
  vector<T> deriv(z);
  rep(i, p.size() - 1) deriv[i] = p[i + 1] * (i + 1);
  vector<T> h(1), q(1);
  for (auto p : {&f, &ffft, &g, &gfft, &h, &q})
    p->reserve(z);
  const T i2 = T(2).inv();
  T im = 1, imm = 1, im2 = i2;
  for (int m = 1; m < z; m <<= 1, im = im2, imm = im * im, im2 *= i2) {
    ffft.assign(all(f)), fft(ffft, m * 2);
    rep(i, m) h[i] = ffft[i] * gfft[i];
    ifft(h), fill(h.begin(), h.begin() + (m + 1) / 2, 0), fft(h);
    rep(i, m) h[i] *= -gfft[i];
    ifft(h);
    rep2(i, (m + 1) / 2, m) g.push_back(h[i] * imm);
    gfft.assign(all(g)), fft(gfft, m * 2);
    h.assign(deriv.begin(), deriv.begin() + m - 1), fft(h, m);
    rep(i, m) q[i] = ffft[i] * h[i];
    ifft(q);
    h.resize(m * 2);
    rep(i, m - 1) h[i + 1] = f[i + 1] * (i + 1);
    rep(i, m) h[i + 1] -= q[i] * im;
    h[0] = -q[m - 1] * im;
    fft(h);
    q.resize(m * 2);
    rep(i, m * 2) q[i] = gfft[i] * h[i];
    ifft(q);
    h.assign(p.begin() + m, p.begin() + min((int)p.size(), m * 2));
    repr(i, m) h[i] -= q[i] * im2 * inverse<T>(i + m);
    fft(h, m * 2);
    rep(i, m * 2) q[i] = ffft[i] * h[i];
    ifft(q);
    f.resize(m * 2);
    rep(i, m) f[i + m] = q[i] * im2;
  }
  f.resize(deg);
  return f;
}
#line 3 "lib/ps/inv.hpp"

template <class T> vector<T> inv(const vector<T> &f, int deg = -1) {
  assert(f[0] != 0);
  if (deg == -1)
    deg = f.size();
  int z = 1 << atcoder::internal::ceil_pow2(deg);
  vector<T> g = {1 / f[0]}, gfft(1);
  vector<T> h;
  g.reserve(z), gfft.reserve(z), h.reserve(z);
  const T i4 = T(4).inv();
  T imm4 = i4;
  for (int m = 1; m < z; m <<= 1, imm4 *= i4) {
    h.assign(f.begin(), f.begin() + min((int)f.size(), m * 2));
    copy(all(g), gfft.begin());
    fft(h, m * 2), fft(gfft, m * 2);
    rep(i, m * 2) h[i] *= gfft[i];
    ifft(h), fill(h.begin(), h.begin() + m, 0), fft(h);
    rep(i, m * 2) h[i] *= -gfft[i];
    ifft(h);
    rep2(i, m, m * 2) g.push_back(h[i] * imm4);
  }
  g.resize(deg);
  return g;
}
#line 3 "lib/ps/bostan_mori.hpp"

template <class T> T bostan_mori(vector<T> a, vector<T> b, ll n) {
  int k = max(2, 1 + atcoder::internal::ceil_pow2(max(a.size(), b.size())));
  int d = 1 << k, half = d / 2;
  a.resize(d), b.resize(d);
  fft(a), fft(b);

  const T w =
      T(atcoder::internal::primitive_root<T::mod()>).pow((T::mod() - 1) >> k);
  vector<T> z(half);
  z[0] = T(d).inv() * 2;
  rep(i, half - 1) z[i + 1] = z[i] * w;
  vector<T> iw(k - 1), y(half); // bit-reversed
  iw[k - 2] = w.inv();
  repr(i, k - 2) iw[i] = iw[i + 1] * iw[i + 1];
  y[0] = 1;
  rep(t, k - 1) rep(i, 1 << t) y[i | 1 << t] = y[i] * iw[t];

  auto dbl = [&](vector<T> &v) {
    copy(v.begin(), v.begin() + half, v.begin() + half);
    ifft(v.begin() + half, v.end());
    rep(i, half) v[i + half] *= z[i];
    fft(v.begin() + half, v.end());
  };

  T c = 1;
  while (true) {
    rep(i, d) a[i] *= b[i ^ 1];
    if (n & 1)
      rep(i, half) a[i] = (a[i << 1] - a[i << 1 | 1]) * y[i];
    else
      rep(i, half) a[i] = a[i << 1] + a[i << 1 | 1];
    if (n == 1) {
      return accumulate(a.begin(), a.begin() + half, T(0)) / (c * d);
    }
    dbl(a);
    rep(i, half) b[i] = b[i << 1] * b[i << 1 | 1];
    dbl(b);
    n >>= 1, c *= 2;
  }
}
#line 9 "lib/ps/fps.hpp"

template <class T, enable_if_t<is_base_of_v<atcoder::internal::modint_base, T>>
                       * = nullptr>
class formal_power_series : public vector<T> {
private:
  using fps = formal_power_series<T>;

public:
  using sparse = vector<pair<int, T>>;

  using vector<T>::vector;
  formal_power_series(vector<T> v) : vector<T>(move(v)) {}
  formal_power_series(sparse p) : vector<T>() {
    resize(p.back().first + 1);
    for (auto[k, c] : p)
      (*this)[k] += c;
  }
  static fps one() { return fps({T(1)}); }
  static fps zero() { return fps{}; }

  int size() const { return vector<T>::size(); }
  const vector<T> &as_vec() const { return (const vector<T> &)*this; }
  vector<T> &as_vec() { return (vector<T> &)*this; }
  vector<T> into_vec() { return move(as_vec()); }
  T eval(T x) const {
    T pow(1), ans(0);
    for (auto e : *this)
      ans += e * pow, pow *= x;
    return ans;
  }
  void trunc() {
    while (!this->empty() && this->back() == T(0))
      this->pop_back();
  }
  fps &mul(T c) {
    for (auto &e : *this)
      e *= c;
    return *this;
  }
  fps add(const fps &v) && {
    this->resize(max(size(), v.size()));
    rep(i, v.size())(*this)[i] += v[i];
    return move(*this);
  }
  fps add(const fps &v) const & { return fps(*this).add(v); }
  fps sub(const fps &v) && {
    this->resize(max(size(), v.size()));
    rep(i, v.size())(*this)[i] -= v[i];
    return move(*this);
  }
  fps sub(const fps &v) const & { return fps(*this).sub(v); }
  fps conv(fps v, int deg = -1) && {
    if (~deg)
      this->resize(min(size(), deg)), v.resize(min(v.size(), deg));
    auto f = convolution(into_vec(), v.into_vec());
    if (~deg)
      f.resize(deg);
    return f;
  }
  fps conv(fps v, int deg = -1) const & {
    return fps(*this).conv(move(v), deg);
  }
  fps diff() && {
    rep(i, size() - 1)(*this)[i] = (*this)[i + 1] * (i + 1);
    this->pop_back();
    return move(*this);
  }
  fps diff() const & { return fps(*this).diff(); }
  fps integr() && {
    this->push_back(0);
    repr(i, size() - 1)(*this)[i + 1] = (*this)[i] * inverse<T>(i + 1);
    (*this)[0] = 0;
    return move(*this);
  }
  fps integr() const & { return fps(*this).integr(); }
  fps inv(int deg = -1) const & { return fps(::inv(as_vec(), deg)); }
  fps div(const fps &v, int deg = -1) && {
    return move(*this).conv(v.inv(deg), deg);
  }
  fps div(const fps &v, int deg = -1) const & { return fps(*this).div(v, deg); }
  fps log(int deg = -1) && {
    return inv(deg - 1).conv(move(*this).diff(), deg - 1).integr();
  }
  fps log(int deg = -1) const & { return fps(*this).log(deg); }
  fps exp(int deg = -1) const & { return ::exp(as_vec(), deg); }
  fps pow(ll k, int deg = -1) && {
    if (deg == -1)
      deg = size();
    int z = -1;
    rep(i, size()) if ((*this)[i] != 0) {
      z = i;
      break;
    }
    if (z == -1 || sat_mul<ll>(z, k) > deg) {
      fps res(deg, 0);
      res[0] = k == 0;
      return res;
    }
    ll rest = deg - z * k;
    this->erase(this->begin(), this->begin() + z);
    T c = (*this)[0].pow(k);
    fps f = move(*this).log(rest).mul(k).exp(rest);
    for (auto &e : f)
      e *= c;
    f.resize(deg);
    copy_backward(f.begin(), f.begin() + rest, f.end());
    fill(f.begin(), f.begin() + z * k, 0);
    return f;
  }
  fps pow(ll k, int deg = -1) const & { return fps(*this).pow(k, deg); }
  fps square(int deg = -1) && {
    if (deg == -1)
      deg = size() * 2 - 1;
    int n = 1 << atcoder::internal::ceil_pow2(deg);
    fft(as_vec(), n * 2);
    for (auto &e : *this)
      e *= e;
    ifft(as_vec());
    auto in2 = inverse<T>(n * 2);
    this->resize(deg);
    for (auto &e : *this)
      e *= in2;
    return move(*this);
  }
  fps square(int deg = -1) const & { return fps(*this).square(deg); }
  T div_at(fps f, ll n) && { return bostan_mori(into_vec(), f.into_vec(), n); }
  T div_at(fps f, ll n) const & { return fps(*this).div_at(move(f), n); }
  optional<fps> sqrt(int deg = -1) && {
    if (deg == -1)
      deg = size();
    this->resize(deg);
    if (this->empty())
      return move(*this);
    if ((*this)[0] == 0) {
      int b = 0;
      while (b < size() && (*this)[b] == 0)
        b++;
      if (b == size())
        return move(*this);
      if (b % 2 != 0)
        return nullopt;
      this->erase(this->begin(), this->begin() + b);
      auto ans = move(*this).sqrt(deg - b / 2);
      if (ans)
        ans->insert(ans->begin(), b / 2, T(0));
      return ans;
    }
    auto x = mod_sqrt((*this)[0]);
    if (!x)
      return nullopt;
    fps f = {*x};
    int z = 1 << atcoder::internal::ceil_pow2(deg);
    f.reserve(z);
    const T i2 = inverse<T>(2);
    for (int m = 1; m < z; m *= 2) {
      fps h(this->begin(), this->begin() + min(m * 2, size()));
      fps hf = move(h).div(f, m * 2);
      f = move(f).add(hf).mul(i2);
    }
    f.resize(deg);
    return f;
  }

  fps div_poly(fps g) && {
    int d = size() - g.size() + 1;
    if (d <= 0)
      return zero();
    reverse(all(*this));
    reverse(all(g));
    fps q = move(*this).div(move(g), d);
    reverse(all(q));
    return q;
  }
  fps div_poly(fps g) const & { return fps(*this).div_poly(move(g)); }
  pair<fps, fps> div_rem_poly(fps g) && {
    int d = g.size() - 1;
    fps q = div_poly(g);
    fps r = move(*this).sub(move(g).conv(q, d));
    r.resize(d);
    r.trunc();
    return pair(move(q), move(r));
  }
  pair<fps, fps> div_rem_poly(fps g) const & {
    return fps(*this).div_rem_poly(move(g));
  }

  fps conv(sparse v) && {
    if (v.empty())
      return zero();
    if (v.front().first == 0)
      v.front().second -= T(1);
    repr(i, size()) for (auto[k, c] : v) {
      if (k > i)
        break;
      (*this)[i] += (*this)[i - k] * c;
    }
    return move(*this);
  }
  fps conv(sparse v) const & { return fps(*this).conv(move(v)); }
  fps div(sparse v) && {
    auto[k0, r] = v.front();
    assert(k0 == 0 && r != T(0));
    T ir = r.inv();
    v.pop_front();
    rep(i, size()) {
      for (auto[k, c] : v) {
        if (k > i)
          break;
        (*this)[i] -= (*this)[i - k] * c;
      }
      (*this)[i] *= ir;
    }
    return move(*this);
  }
  fps div(sparse v) const & { return fps(*this).div(move(v)); }

  template <class It> static fps prod(It a, It a_last) {
    if (a == a_last)
      return one();
    vector<fps> vec(a, a_last);
    vec.reserve(distance(a, a_last) * 2);
    for (int i = 0; i + 1 < vec.size(); i += 2)
      vec.push_back(move(vec[i]).conv(move(vec[i + 1])));
    return vec.back();
  }
};
#line 6 "main.cpp"
using mint = atcoder::modint998244353;
using fps = formal_power_series<mint>;

int main() {
  int n = in;
  ll m = in;
  fps f(in.vec<mint>(n));
  f = move(f).pow(m);
  out.iter(all(f));
}
// https://judge.yosupo.jp/submission/103275