#define NDEBUG
#include <bits/stdc++.h>
#include <immintrin.h>

#include <algorithm>

#include <cassert>

#include <chrono>

#include <cstdint>

#include <cstdio>

#include <cstring>

#include <iostream>

#include <random>

#include <type_traits>

#include <utility>

#include <vector>

using namespace std;
#include <type_traits>
using namespace std;
namespace fastio {
static constexpr int SZ = 1 << 17;
char inbuf[SZ], outbuf[SZ];
int in_left = 0, in_right = 0, out_right = 0;
struct Pre {
  char num[40000];
  constexpr Pre() : num() {
    for (int i = 0; i < 10000; i++) {
      int n = i;
      for (int j = 3; j >= 0; j--) {
        num[i * 4 + j] = n % 10 + '0';
        n /= 10;
      }
    }
  }
} constexpr pre;
inline void load() {
  int len = in_right - in_left;
  memmove(inbuf, inbuf + in_left, len);
  in_right = len + fread(inbuf + len, 1, SZ - len, stdin);
  in_left = 0;
}
inline void flush() {
  fwrite(outbuf, 1, out_right, stdout);
  out_right = 0;
}
inline void skip_space() {
  if (in_left + 32 > in_right)
    load();
  while (inbuf[in_left] <= ' ')
    in_left++;
}
inline void rd(char &c) {
  if (in_left + 32 > in_right)
    load();
  c = inbuf[in_left++];
}
template <typename T> inline void rd(T &x) {
  if (in_left + 32 > in_right)
    load();
  char c;
  do
    c = inbuf[in_left++];
  while (c < '-');
  [[maybe_unused]] bool minus = false;
  if
    constexpr(is_signed<T>::value == true) {
      if (c == '-')
        minus = true, c = inbuf[in_left++];
    }
  x = 0;
  while (c >= '0') {
    x = x * 10 + (c & 15);
    c = inbuf[in_left++];
  }
  if
    constexpr(is_signed<T>::value == true) {
      if (minus)
        x = -x;
    }
}
inline void rd() {}
template <typename Head, typename... Tail>
inline void rd(Head &head, Tail &... tail) {
  rd(head);
  rd(tail...);
}
inline void wt(char c) {
  if (out_right > SZ - 32)
    flush();
  outbuf[out_right++] = c;
}
inline void wt(bool b) {
  if (out_right > SZ - 32)
    flush();
  outbuf[out_right++] = b ? '1' : '0';
}
template <typename T> inline void wt(T x) {
  if (out_right > SZ - 32)
    flush();
  if (!x) {
    outbuf[out_right++] = '0';
    return;
  }
  if
    constexpr(is_signed<T>::value == true) {
      if (x < 0)
        outbuf[out_right++] = '-', x = -x;
    }
  int i = 12;
  char buf[16];
  while (x >= 10000) {
    memcpy(buf + i, pre.num + (x % 10000) * 4, 4);
    x /= 10000;
    i -= 4;
  }
  if (x < 100) {
    if (x < 10) {
      outbuf[out_right] = '0' + x;
      ++out_right;
    } else {
      uint32_t q = (uint32_t(x) * 205) >> 11;
      uint32_t r = uint32_t(x) - q * 10;
      outbuf[out_right] = '0' + q;
      outbuf[out_right + 1] = '0' + r;
      out_right += 2;
    }
  } else {
    if (x < 1000) {
      memcpy(outbuf + out_right, pre.num + (x << 2) + 1, 3);
      out_right += 3;
    } else {
      memcpy(outbuf + out_right, pre.num + (x << 2), 4);
      out_right += 4;
    }
  }
  memcpy(outbuf + out_right, buf + i + 4, 12 - i);
  out_right += 12 - i;
}
inline void wt() {}
template <typename Head, typename... Tail>
inline void wt(Head &&head, Tail &&... tail) {
  wt(head);
  wt(forward<Tail>(tail)...);
}
template <typename... Args> inline void wtn(Args &&... x) {
  wt(forward<Args>(x)...);
  wt('\n');
}
struct Dummy {
  Dummy() { atexit(flush); }
} dummy;
} // namespace fastio
using fastio::rd;
using fastio::skip_space;
using fastio::wt;
using fastio::wtn;
//
template <uint32_t mod> struct LazyMontgomeryModInt {
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
//
#pragma GCC optimize("O3,unroll-loops")
#pragma GCC target("avx2")
using m256 = __m256i;
struct alignas(32) mmint {
  m256 x;
  static m256 R, mR, M0, M1, M2, N2;
  mmint() : x() {}
  inline mmint(const m256 &_x) : x(_x) {}
  inline mmint(int a) : x(_mm256_set1_epi32(a)) {}
  inline mmint(int a0, int a1, int a2, int a3, int a4, int a5, int a6, int a7)
      : x(_mm256_set_epi32(a7, a6, a5, a4, a3, a2, a1, a0)) {}
  inline operator m256 &() { return x; }
  inline operator const m256 &() const { return x; }
  inline int &operator[](int i) { return *(reinterpret_cast<int *>(&x) + i); }
  inline const int &operator[](int i) const {
    return *(reinterpret_cast<const int *>(&x) + i);
  }
  friend ostream &operator<<(ostream &os, const mmint &m) {
    unsigned r = R[0], mod = M1[0];
    auto reduce1 = [&](const uint64_t &b) {
      unsigned res = (b + uint64_t(unsigned(b) * unsigned(-r)) * mod) >> 32;
      return res >= mod ? res - mod : res;
    };
    for (int i = 0; i < 8; i++) {
      os << reduce1(m[i]) << (i == 7 ? "" : " ");
    }
    return os;
  }
  template <typename mint> static void set_mod() {
    R = _mm256_set1_epi32(mint::r);
    M0 = _mm256_setzero_si256();
    M1 = _mm256_set1_epi32(mint::get_mod());
    M2 = _mm256_set1_epi32(mint::get_mod() * 2);
    N2 = _mm256_set1_epi32(mint::n2);
    mR = _mm256_sub_epi32(M0, R);
  }
  static inline mmint reduce(const mmint &prod02, const mmint &prod13) {
    m256 unpalo = _mm256_unpacklo_epi32(prod02, prod13);
    m256 unpahi = _mm256_unpackhi_epi32(prod02, prod13);
    m256 prodlo = _mm256_unpacklo_epi64(unpalo, unpahi);
    m256 prodhi = _mm256_unpackhi_epi64(unpalo, unpahi);
    m256 hiplm1 = _mm256_add_epi32(prodhi, M1);
    m256 lomulr = _mm256_mullo_epi32(prodlo, R);
    m256 a13 = _mm256_shuffle_epi32(lomulr, 0xF5);
    m256 prod02_ = _mm256_mul_epu32(lomulr, M1);
    m256 prod13_ = _mm256_mul_epu32(a13, M1);
    m256 unpalo_ = _mm256_unpacklo_epi32(prod02_, prod13_);
    m256 unpahi_ = _mm256_unpackhi_epi32(prod02_, prod13_);
    m256 prod = _mm256_unpackhi_epi64(unpalo_, unpahi_);
    return _mm256_sub_epi32(hiplm1, prod);
  }
  static inline mmint itom(const mmint &A) { return A * N2; }
  static inline mmint mtoi(const mmint &A) {
    m256 lomulr = _mm256_mullo_epi32(A, mR);
    m256 a13 = _mm256_shuffle_epi32(lomulr, 0xF5);
    m256 prod02_ = _mm256_mul_epu32(lomulr, M1);
    m256 prod13_ = _mm256_mul_epu32(a13, M1);
    m256 unpalo_ = _mm256_unpacklo_epi32(prod02_, prod13_);
    m256 unpahi_ = _mm256_unpackhi_epi32(prod02_, prod13_);
    m256 prod = _mm256_unpackhi_epi64(unpalo_, unpahi_);
    m256 cmpA = _mm256_cmpgt_epi32(A, M0);
    m256 addA = _mm256_and_si256(cmpA, mmint(1));
    return _mm256_add_epi32(prod, addA);
  }
  friend inline mmint operator+(const mmint &A, const mmint &B) {
    m256 apb = _mm256_add_epi32(A, B);
    m256 ret = _mm256_sub_epi32(apb, M2);
    m256 cmp = _mm256_cmpgt_epi32(M0, ret);
    m256 add = _mm256_and_si256(cmp, M2);
    return _mm256_add_epi32(add, ret);
  }
  friend inline mmint operator-(const mmint &A, const mmint &B) {
    m256 ret = _mm256_sub_epi32(A, B);
    m256 cmp = _mm256_cmpgt_epi32(M0, ret);
    m256 add = _mm256_and_si256(cmp, M2);
    return _mm256_add_epi32(add, ret);
  }
  friend inline mmint operator*(const mmint &A, const mmint &B) {
    m256 a13 = _mm256_shuffle_epi32(A, 0xF5);
    m256 b13 = _mm256_shuffle_epi32(B, 0xF5);
    m256 prod02 = _mm256_mul_epu32(A, B);
    m256 prod13 = _mm256_mul_epu32(a13, b13);
    return reduce(prod02, prod13);
  }
  inline mmint &operator+=(const mmint &A) { return (*this) = (*this) + A; }
  inline mmint &operator-=(const mmint &A) { return (*this) = (*this) - A; }
  inline mmint &operator*=(const mmint &A) { return (*this) = (*this) * A; }
  bool operator==(const mmint &A) {
    m256 sub = _mm256_sub_epi32(x, A.x);
    return _mm256_testz_si256(sub, sub) == 1;
  }
  bool operator!=(const mmint &A) { return !((*this) == A); }
};
__attribute__((aligned(32))) __m256i mmint::R, mmint::mR;
__attribute__((aligned(32))) __m256i mmint::M0, mmint::M1, mmint::M2, mmint::N2;
/**
 * @brief vectorize modint
 */
constexpr int B = 1 << 7;
// B*Bの正方行列を高速に乗算するライブラリ。
// B*B行列a,bを タテB行 ヨコB/8行の行列と見なす.
// s : 正順に配置。すなわちa_{i,k}をs[i * (B / 8) + k]に配置する。
// t : 逆順に配置。すなわちb_{k,j}をt[j * B + k]に配置する。
// u : 正順に配置。すなわちc_{i,j}をu[i * (B / 8) + j]に配置する。
namespace block_mat_prod_impl {
constexpr int B8 = B / 8;
void mul_simd(mmint *__restrict__ s, mmint *__restrict__ t,
              mmint *__restrict__ u) {
  // normalize
  for (int i = 0; i < B * B8; i++) {
    const m256 cmpS = _mm256_cmpgt_epi32(s[i], mmint::M1);
    const m256 cmpT = _mm256_cmpgt_epi32(t[i], mmint::M1);
    const m256 difS = _mm256_and_si256(cmpS, mmint::M1);
    const m256 difT = _mm256_and_si256(cmpT, mmint::M1);
    s[i] = _mm256_sub_epi32(s[i], difS);
    t[i] = _mm256_sub_epi32(t[i], difT);
  }
  // mul
  mmint th1, th2, zero = _mm256_setzero_si256();
  th1[1] = th1[3] = th1[5] = th1[7] = mmint::M1[0];
  th2[1] = th2[3] = th2[5] = th2[7] = mmint::M2[0];
#define INIT_X(x, y)                                                           \
  m256 prod02##x##y = _mm256_setzero_si256();                                  \
  m256 prod13##x##y = _mm256_setzero_si256()
#define INIT_Y(j, k, l, y)                                                     \
  m256 T##y = t[(j + y) * B + k + l];                                          \
  const m256 T13##y = _mm256_shuffle_epi32(T##y, 0xF5);
#define PROD(x, y)                                                             \
  m256 S##x##y = _mm256_set1_epi32(s[(i + x) * B8 + k / 8][l]);                \
  const m256 ST02##x##y = _mm256_mul_epu32(S##x##y, T##y);                     \
  const m256 ST13##x##y = _mm256_mul_epu32(S##x##y, T13##y);                   \
  prod02##x##y = _mm256_add_epi64(prod02##x##y, ST02##x##y);                   \
  prod13##x##y = _mm256_add_epi64(prod13##x##y, ST13##x##y)
#define COMP(x, y)                                                             \
  m256 cmp02##x##y = _mm256_cmpgt_epi64(zero, prod02##x##y);                   \
  m256 cmp13##x##y = _mm256_cmpgt_epi64(zero, prod13##x##y);                   \
  m256 dif02##x##y = _mm256_and_si256(cmp02##x##y, th2);                       \
  m256 dif13##x##y = _mm256_and_si256(cmp13##x##y, th2);                       \
  prod02##x##y = _mm256_sub_epi64(prod02##x##y, dif02##x##y);                  \
  prod13##x##y = _mm256_sub_epi64(prod13##x##y, dif13##x##y)
#define REDUCE(x, y)                                                           \
  for (int _ = 0; _ < 2; _++) {                                                \
    m256 cmp02 = _mm256_cmpgt_epi64(prod02##x##y, th1);                        \
    m256 cmp13 = _mm256_cmpgt_epi64(prod13##x##y, th1);                        \
    m256 dif02 = _mm256_and_si256(cmp02, th1);                                 \
    m256 dif13 = _mm256_and_si256(cmp13, th1);                                 \
    prod02##x##y = _mm256_sub_epi64(prod02##x##y, dif02);                      \
    prod13##x##y = _mm256_sub_epi64(prod13##x##y, dif13);                      \
  }                                                                            \
  u[(i + x) * B8 + j + y] = mmint::reduce(prod02##x##y, prod13##x##y)
  for (int i = 0; i < B; i += 4) {
    for (int j = 0; j < B8; j += 2) {
      INIT_X(0, 0);
      INIT_X(1, 0);
      INIT_X(0, 1);
      INIT_X(1, 1);
      INIT_X(2, 0);
      INIT_X(2, 1);
      INIT_X(3, 0);
      INIT_X(3, 1);
      for (int k = 0; k < B; k += 8) {
        for (int l = 0; l < 8; l++) {
          INIT_Y(j, k, l, 0);
          INIT_Y(j, k, l, 1);
          PROD(0, 0);
          PROD(0, 1);
          PROD(1, 0);
          PROD(1, 1);
          PROD(2, 0);
          PROD(2, 1);
          PROD(3, 0);
          PROD(3, 1);
        }
        COMP(0, 0);
        COMP(0, 1);
        COMP(1, 0);
        COMP(1, 1);
        COMP(2, 0);
        COMP(2, 1);
        COMP(3, 0);
        COMP(3, 1);
      }
      REDUCE(0, 0);
      REDUCE(0, 1);
      REDUCE(1, 0);
      REDUCE(1, 1);
      REDUCE(2, 0);
      REDUCE(2, 1);
      REDUCE(3, 0);
      REDUCE(3, 1);
    }
  }
}
#undef INIT
#undef PROD
#undef COMP
#undef REDUCE
} // namespace block_mat_prod_impl
using block_mat_prod_impl::mul_simd;
void strassen(int N, mmint *__restrict__ s, mmint *__restrict__ t,
              mmint *__restrict__ u) {
  memset((void *)u, 0, N * N * sizeof(int));
  if (N == B) {
    mul_simd(s, t, u);
    return;
  }
  mmint *ps = s + N * N / 8;
  mmint *pt = t + N * N / 8;
  mmint *pu = u + N * N / 8;
  int nx = N * N / 32;
  int o11 = nx * 0, o12 = nx * 1, o21 = nx * 2, o22 = nx * 3;
  // P1
  for (int i = 0; i < nx; i++)
    ps[i] = s[o11 + i] + s[o22 + i];
  for (int i = 0; i < nx; i++)
    pt[i] = t[o11 + i] + t[o22 + i];
  strassen(N / 2, ps, pt, pu);
  for (int i = 0; i < nx; i++)
    u[o11 + i] = pu[i], u[o22 + i] = pu[i];
  // P2
  for (int i = 0; i < nx; i++)
    ps[i] = s[o21 + i] + s[o22 + i];
  for (int i = 0; i < nx; i++)
    pt[i] = t[o11 + i];
  strassen(N / 2, ps, pt, pu);
  for (int i = 0; i < nx; i++)
    u[o21 + i] = pu[i], u[o22 + i] -= pu[i];
  // P3
  for (int i = 0; i < nx; i++)
    ps[i] = s[o11 + i];
  for (int i = 0; i < nx; i++)
    pt[i] = t[o12 + i] - t[o22 + i];
  strassen(N / 2, ps, pt, pu);
  for (int i = 0; i < nx; i++)
    u[o12 + i] = pu[i], u[o22 + i] += pu[i];
  // P4
  for (int i = 0; i < nx; i++)
    ps[i] = s[o22 + i];
  for (int i = 0; i < nx; i++)
    pt[i] = t[o21 + i] - t[o11 + i];
  strassen(N / 2, ps, pt, pu);
  for (int i = 0; i < nx; i++)
    u[o11 + i] += pu[i], u[o21 + i] += pu[i];
  // P5
  for (int i = 0; i < nx; i++)
    ps[i] = s[o11 + i] + s[o12 + i];
  for (int i = 0; i < nx; i++)
    pt[i] = t[o22 + i];
  strassen(N / 2, ps, pt, pu);
  for (int i = 0; i < nx; i++)
    u[o11 + i] -= pu[i], u[o12 + i] += pu[i];
  // P6
  for (int i = 0; i < nx; i++)
    ps[i] = s[o21 + i] - s[o11 + i];
  for (int i = 0; i < nx; i++)
    pt[i] = t[o11 + i] + t[o12 + i];
  strassen(N / 2, ps, pt, pu);
  for (int i = 0; i < nx; i++)
    u[o22 + i] += pu[i];
  // P7
  for (int i = 0; i < nx; i++)
    ps[i] = s[o12 + i] - s[o22 + i];
  for (int i = 0; i < nx; i++)
    pt[i] = t[o21 + i] + t[o22 + i];
  strassen(N / 2, ps, pt, pu);
  for (int i = 0; i < nx; i++)
    u[o11 + i] += pu[i];
}
constexpr int S = 1024;
constexpr int S8 = S / 8;
void place_s(int N, int a, int b, mmint *__restrict__ dst,
             mmint *__restrict__ src) {
  if (N == B) {
    int B8 = B / 8;
    for (int i = 0; i < B; i++) {
      memcpy(dst + i * B8, src + (a + i) * S8 + b / 8, B8 * sizeof(mmint));
    }
    return;
  }
  int nx = N * N / 32, M = N / 2;
  place_s(M, a + 0, b + 0, dst + nx * 0, src);
  place_s(M, a + 0, b + M, dst + nx * 1, src);
  place_s(M, a + M, b + 0, dst + nx * 2, src);
  place_s(M, a + M, b + M, dst + nx * 3, src);
}
void place_t(int N, int a, int b, mmint *__restrict__ dst,
             mmint *__restrict__ src) {
  if (N == B) {
    int B8 = B / 8;
    // t : 逆順に配置。すなわちb_{k,j}をt[j * B + k]に配置する。
    for (int k = 0; k < B; k++) {
      for (int j = 0; j < B8; j++) {
        dst[j * B + k] = src[(a + k) * S8 + j + b / 8];
      }
    }
    return;
  }
  int nx = N * N / 32, M = N / 2;
  place_t(M, a + 0, b + 0, dst + nx * 0, src);
  place_t(M, a + 0, b + M, dst + nx * 1, src);
  place_t(M, a + M, b + 0, dst + nx * 2, src);
  place_t(M, a + M, b + M, dst + nx * 3, src);
}
void place_rev(int N, int a, int b, mmint *__restrict__ dst,
               mmint *__restrict__ src) {
  if (N == B) {
    int w = B / 8;
    for (int i = 0; i < B; i++) {
      memcpy(src + (a + i) * S8 + b / 8, dst + i * w, w * sizeof(mmint));
    }
    return;
  }
  int nx = N * N / 32, M = N / 2;
  place_rev(M, a + 0, b + 0, dst + nx * 0, src);
  place_rev(M, a + 0, b + M, dst + nx * 1, src);
  place_rev(M, a + M, b + 0, dst + nx * 2, src);
  place_rev(M, a + M, b + M, dst + nx * 3, src);
}
mmint s[S * S8 * 3 / 2], t[S * S8 * 3 / 2], u[S * S8 * 3 / 2];
void prod(unsigned int *__restrict__ a, unsigned int *__restrict__ b,
          unsigned int *__restrict__ c) {
  place_s(S, 0, 0, s, reinterpret_cast<mmint *>(a));
  place_t(S, 0, 0, t, reinterpret_cast<mmint *>(b));
  for (int i = 0; i < S * S8; i++)
    s[i] = mmint::itom(s[i]);
  for (int i = 0; i < S * S8; i++)
    t[i] = mmint::itom(t[i]);
  strassen(S, s, t, u);
  for (int i = 0; i < S * S8; i++)
    u[i] = mmint::mtoi(u[i]);
  place_rev(S, 0, 0, u, reinterpret_cast<mmint *>(c));
}
void test() {
  unsigned int *a = reinterpret_cast<unsigned int *>(t);
  unsigned int *b = reinterpret_cast<unsigned int *>(u);
  unsigned int *c = reinterpret_cast<unsigned int *>(s);
  for (int i = 0; i < S; i++) {
    for (int j = 0; j < S; j++) {
      b[i * S + j] = a[i * S + j] = i + j;
    }
  }
  for (int loop = 0; loop < 100; loop++)
    prod(a, b, c);
}

int main(int argc, char **argv) {
  freopen("Matrix_Multiplication.in", "r", stdin);
  // freopen("Matrix_Multiplication.out", "w", stdout);
  using mint = LazyMontgomeryModInt<998244353>;
  mmint::set_mod<mint>();
  // #ifdef PROFILER
  //   test();
  //   return 0;
  // #endif
  int N, M, K;
  cin >> N >> M >> K;
  // rd(N, M, K);
  unsigned int *a = reinterpret_cast<unsigned int *>(t);
  unsigned int *b = reinterpret_cast<unsigned int *>(u);
  unsigned int *c = reinterpret_cast<unsigned int *>(s);
  unsigned int x;
  for (int i = 0; i < N; i++) {
    for (int k = 0; k < M; k++) {
      // rd(x);
      cin >> x;
      // x = i + k;
      a[i * S + k] = x;
    }
  }
  for (int k = 0; k < M; k++) {
    for (int j = 0; j < K; j++) {
      // rd(x);
      cin >> x;
      // x = k - j;
      b[k * S + j] = x;
    }
  }
  prod(a, b, c);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < K; j++) {
      x = c[i * S + j];
      // wt(x);
      cout << x << " ";
      wt(j == K - 1 ? '\n' : ' ');
    }
    cout << endl;
  }
}
