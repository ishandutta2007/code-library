// #pragma once
#include <bits/stdc++.h>
using namespace std;

#define each(x, v) for (auto &x : v)
#define rep(i, N) for (long long i = 0; i < (long long)(N); i++)

#include <chrono>

struct Timer {
  chrono::high_resolution_clock::time_point st;

  Timer() { reset(); }

  void reset() { st = chrono::high_resolution_clock::now(); }

  chrono::milliseconds::rep elapsed() {
    auto ed = chrono::high_resolution_clock::now();
    return chrono::duration_cast<chrono::milliseconds>(ed - st).count();
  }
};

template <uint32_t mod>

struct LazyMontgomeryModInt {
  using mint = LazyMontgomeryModInt;
  using i32 = int32_t;
  using u32 = uint32_t;
  using u64 = uint64_t;

  static constexpr u32 get_r() {
    u32 ret = mod;
    for (i32 i = 0; i < 4; ++i)
      ret *= 2 - mod * ret;
    return ret;
  }

  static constexpr u32 r = get_r();
  static constexpr u32 n2 = -u64(mod) % mod;
  static_assert(r * mod == 1, "invalid, r * mod != 1");
  static_assert(mod < (1 << 30), "invalid, mod >= 2 ^ 30");
  static_assert((mod & 1) == 1, "invalid, mod % 2 == 0");

  u32 a;

  constexpr LazyMontgomeryModInt() : a(0) {}
  constexpr LazyMontgomeryModInt(const int64_t &b)
      : a(reduce(u64(b % mod + mod) * n2)){};

  static constexpr u32 reduce(const u64 &b) {
    return (b + u64(u32(b) * u32(-r)) * mod) >> 32;
  }

  constexpr mint &operator+=(const mint &b) {
    if (i32(a += b.a - 2 * mod) < 0)
      a += 2 * mod;
    return *this;
  }

  constexpr mint &operator-=(const mint &b) {
    if (i32(a -= b.a) < 0)
      a += 2 * mod;
    return *this;
  }

  constexpr mint &operator*=(const mint &b) {
    a = reduce(u64(a) * b.a);
    return *this;
  }

  constexpr mint &operator/=(const mint &b) {
    *this *= b.inverse();
    return *this;
  }

  constexpr mint operator+(const mint &b) const { return mint(*this) += b; }
  constexpr mint operator-(const mint &b) const { return mint(*this) -= b; }
  constexpr mint operator*(const mint &b) const { return mint(*this) *= b; }
  constexpr mint operator/(const mint &b) const { return mint(*this) /= b; }
  constexpr bool operator==(const mint &b) const {
    return (a >= mod ? a - mod : a) == (b.a >= mod ? b.a - mod : b.a);
  }
  constexpr bool operator!=(const mint &b) const {
    return (a >= mod ? a - mod : a) != (b.a >= mod ? b.a - mod : b.a);
  }
  constexpr mint operator-() const { return mint() - mint(*this); }

  constexpr mint pow(u64 n) const {
    mint ret(1), mul(*this);
    while (n > 0) {
      if (n & 1)
        ret *= mul;
      mul *= mul;
      n >>= 1;
    }
    return ret;
  }

  constexpr mint inverse() const { return pow(mod - 2); }

  friend ostream &operator<<(ostream &os, const mint &b) {
    return os << b.get();
  }

  friend istream &operator>>(istream &is, mint &b) {
    int64_t t;
    is >> t;
    b = LazyMontgomeryModInt<mod>(t);
    return (is);
  }

  constexpr u32 get() const {
    u32 ret = reduce(a);
    return ret >= mod ? ret - mod : ret;
  }

  static constexpr u32 get_mod() { return mod; }
};

using mint = LazyMontgomeryModInt<998244353>;

template <typename T> struct Schoenhage_Strassen_radix2 {
  T *buf = nullptr;

  void rot(T *dst, T *src, int s, int d) {
    assert(0 <= d and d < 2 * s);
    bool f = s <= d;
    if (s <= d)
      d -= s;
    int i = 0;
    if (f) {
      for (; i < s - d; i++)
        dst[i + d] = -src[i];
      for (; i < s; i++)
        dst[i + d - s] = src[i];
    } else {
      for (; i < s - d; i++)
        dst[i + d] = src[i];
      for (; i < s; i++)
        dst[i + d - s] = -src[i];
    }
  }

  void in_add(T *dst, T *src, int s) {
    for (int i = 0; i < s; i++)
      dst[i] += src[i];
  }
  void in_sub(T *dst, T *src, int s) {
    for (int i = 0; i < s; i++)
      dst[i] -= src[i];
  }

  void cp(T *dst, T *src, int s) { memcpy(dst, src, s * sizeof(T)); }
  void reset(T *dst, int s) { fill(dst, dst + s, T()); }

  // R[x] / (1 + x^(2^m)) 上の長さ2^LのFFT
  void fft(T *a, int l, int m) {
    if (l == 0)
      return;
    int L = 1 << l, M = 1 << m;
    assert(M * 2 >= L);
    assert(buf != nullptr);

    vector<int> dw(l - 1);
    for (int i = 0; i < l - 1; i++) {
      dw[i] = (1 << (l - 2 - i)) + (1 << (l - 1 - i)) - (1 << (l - 1));
      if (dw[i] < 0)
        dw[i] += L;
      if (L == M)
        dw[i] *= 2;
      if (2 * L == M)
        dw[i] *= 4;
    }

    for (int d = L; d >>= 1;) {
      int w = 0;
      for (int s = 0, k = 0;;) {
        for (int i = s, j = s + d; i < s + d; ++i, ++j) {
          T *ai = a + i * M, *aj = a + j * M;
          rot(buf, aj, M, w);
          cp(aj, ai, M);
          in_add(ai, buf, M);
          in_sub(aj, buf, M);
        }
        if ((s += 2 * d) >= L)
          break;
        w += dw[__builtin_ctz(++k)];
        if (w >= 2 * M)
          w -= 2 * M;
      }
    }
  }

  // R[x] / (1 + x^(2^m)) 上の長さ2^LのIFFT
  void ifft(T *a, int l, int m) {
    if (l == 0)
      return;
    int L = 1 << l, M = 1 << m;
    assert(M * 2 >= L);
    assert(buf != nullptr);

    vector<int> dw(l - 1);
    for (int i = 0; i < l - 1; i++) {
      dw[i] = (1 << (l - 2 - i)) + (1 << (l - 1 - i)) - (1 << (l - 1));
      if (dw[i] < 0)
        dw[i] += L;
      if (L == M)
        dw[i] *= 2;
      if (2 * L == M)
        dw[i] *= 4;
    }

    for (int d = 1; d < L; d *= 2) {
      int w = 0;
      for (int s = 0, k = 0;;) {
        for (int i = s, j = s + d; i < s + d; ++i, ++j) {
          T *ai = a + i * M, *aj = a + j * M;
          cp(buf, ai, M);
          in_add(ai, aj, M);
          in_sub(buf, aj, M);
          rot(aj, buf, M, w);
        }
        if ((s += 2 * d) >= L)
          break;
        w -= dw[__builtin_ctz(++k)];
        if (w < 0)
          w += 2 * M;
      }
    }
  }

  // a <- ab / (x^(2^n)+1)
  int naive_mul(T *a, T *b, int n) {
    int N = 1 << n;
    reset(buf, N);
    for (int i = 0; i < N; i++) {
      int j = 0;
      for (; j < N - i; j++)
        buf[i + j] += a[i] * b[j];
      for (; j < N; j++)
        buf[i + j - N] -= a[i] * b[j];
    }
    cp(a, buf, N);
    return 0;
  }

  // a <- ab / (x^(2^n)+1)
  int inplace_mul(T *a, T *b, int n) {
    if (n <= 5) {
      naive_mul(a, b, n);
      return 0;
    }

    int l = (n + 1) / 2;
    int m = n / 2;
    int L = 1 << l, M = 1 << m, N = 1 << n;
    int cnt = 0;

    // R[x] (x^(2^(m+1))-1) R[y] (y^(2^l)-1)
    vector<T> A(N * 2), B(N * 2);
    reset(buf + M, M);
    for (int i = 0, s = 0, ds = 2 * M / L; i < L; i++) {
      // y -> x^{2m/r} y
      cp(buf, a + i * M, M);
      rot(A.data() + i * M * 2, buf, 2 * M, s);
      cp(buf, b + i * M, M);
      rot(B.data() + i * M * 2, buf, 2 * M, s);
      s += ds;
      if (s >= 4 * M)
        s -= 4 * M;
    }

    fft(A.data(), l, m + 1);
    fft(B.data(), l, m + 1);
    for (int i = 0; i < L; i++) {
      cnt = inplace_mul(A.data() + i * M * 2, B.data() + i * M * 2, m + 1);
    }
    ifft(A.data(), l, m + 1);

    for (int i = 0, s = 0, ds = 2 * M / L; i < L; i++) {
      // y -> x^{-2m/r} y
      cp(buf, A.data() + i * M * 2, 2 * M);
      rot(A.data() + i * M * 2, buf, 2 * M, s);
      s -= ds;
      if (s < 0)
        s += 4 * M;
    }

    for (int i = 0; i < L; i++) {
      for (int j = 0; j < M; j++) {
        a[i * M + j] = A[i * M * 2 + j];
        if (i == 0) {
          a[i * M + j] -= A[(L - 1) * M * 2 + M + j];
        } else {
          a[i * M + j] += A[(i - 1) * M * 2 + M + j];
        }
      }
    }
    return cnt + l;
  }

  // a <- ab / (x^(2^n)-1)
  int inplace_mul2(T *a, T *b, int n) {
    if (n <= 5) {
      naive_mul(a, b, n);
      return 0;
    }

    int l = (n + 1) / 2;
    int m = n / 2;
    int L = 1 << l, M = 1 << m, N = 1 << n;
    int cnt = 0;

    // R[x] (x^(2^(m+1))-1) R[y] (y^(2^l)-1)
    vector<T> A(N * 2), B(N * 2);
    for (int i = 0; i < L; i++) {
      cp(A.data() + i * M * 2, a + i * M, M);
      cp(B.data() + i * M * 2, b + i * M, M);
    }
    fft(A.data(), l, m + 1);
    fft(B.data(), l, m + 1);
    for (int i = 0; i < L; i++) {
      cnt = inplace_mul(A.data() + i * M * 2, B.data() + i * M * 2, m + 1);
    }
    ifft(A.data(), l, m + 1);
    for (int i = 0; i < L; i++) {
      for (int j = 0; j < M; j++) {
        a[i * M + j] = A[i * M * 2 + j];
        a[i * M + j] += A[(i ? i - 1 : L - 1) * M * 2 + M + j];
      }
    }
    return cnt + l;
  }

  pair<vector<T>, int> multiply(const vector<T> &a, const vector<T> &b) {
    int L = a.size() + b.size() - 1;
    int M = 1, m = 0;
    while (M < L)
      M *= 2, m++;
    buf = new T[M];
    vector<T> s(M), t(M);
    for (int i = 0; i < (int)a.size(); i++)
      s[i] = a[i];
    for (int i = 0; i < (int)b.size(); i++)
      t[i] = b[i];
    int cnt = inplace_mul2(s.data(), t.data(), m);
    vector<T> u(L);
    for (int i = 0; i < L; i++)
      u[i] = s[i];
    delete[] buf;
    return make_pair(u, cnt);
  }
};

Schoenhage_Strassen_radix2<mint> ss;
void calc() {
  int N = 100000 * 5;
  vector<mint> a(N), b(N);
  int e = 1;
  each(x, a) x = e, ++e;
  each(x, b) x = e, ++e;
  Timer timer;
  auto c = ss.multiply(a, b);
  cerr << c.first.size() << " " << timer.elapsed() << endl;
}
void test() {
  for (int n = 1, m = 1; n < 100; n += 2, m += 2) {
    vector<mint> a(n), b(m);
    int e = 1;
    each(x, a) x = e, e += 1;
    // each(x, b) x = e, e += 1;
    b = a;
    auto[c, d] = ss.multiply(a, b);
    mint inv = mint(1 << d).inverse();
    each(x, c) x *= inv;
    vector<mint> C(n + m - 1);
    rep(i, n) rep(j, m) C[i + j] += a[i] * b[j];
    assert(C == c);
  }
}

void in() {}
template <typename T, class... U> void in(T &t, U &... u) {
  cin >> t;
  in(u...);
}

void out() { cout << "\n"; }
template <typename T, class... U, char sep = ' '>
void out(const T &t, const U &... u) {
  cout << t;
  if (sizeof...(u))
    cout << sep;
  out(u...);
}

#define ini(...)                                                               \
  int __VA_ARGS__;                                                             \
  in(__VA_ARGS__)

int main() {
  // test();
  // calc();

  ini(N, M);
  vector<mint> a(N), b(M);
  in(a, b);
  auto[c, d] = ss.multiply(a, b);
  mint inv = mint(1 << d).inverse();
  each(x, c) x *= inv;

  out(c);
}
