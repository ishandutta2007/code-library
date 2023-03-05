#include <bits/stdc++.h>
using u32 = unsigned;
using u64 = unsigned long long;
using u128 = __uint128_t;

struct istream {
  static const u32 size = 1 << 22;
  char buf[size], *vin;
  inline istream() {
    fread(buf, 1, size, stdin);
    vin = buf - 1;
  }
  inline istream &operator>>(u32 &x) {
    for (x = *++vin & 15; isdigit(*++vin);)
      x = x * 10 + (*vin & 15);
    return *this;
  }
} cin;
struct ostream {
  static const u32 size = 1 << 23;
  char buf[size], *vout;
  unsigned map[10000];
  inline ostream() {
    for (u32 i = 0; i < 10000; ++i) {
      map[i] = (i % 10 + 48) << 24 | (i / 10 % 10 + 48) << 16 |
               (i / 100 % 10 + 48) << 8 | (i / 1000 + 48);
    }
    vout = buf + size;
  }
  inline ~ostream() { fwrite(vout, 1, buf + size - vout, stdout); }
  inline ostream &operator<<(u32 x) {
    if (x) {
      for (; x >= 1000; x /= 10000)
        *--(unsigned *&)vout = map[x % 10000];
      while (x)
        *--vout = x % 10 + 48, x /= 10;
    } else {
      *--vout = 48;
    }
    return *this;
  }
  inline ostream &operator<<(char x) {
    *--vout = x;
    return *this;
  }
} cout;

namespace MPE {
constexpr u32 N = 1 << 18 | 1;
u32 mod;
u32 mod2 = mod * 2;

inline u64 getdiv(u64 x) {
  u64 base = (-1ull) / mod;
  return u64((u128)x * base >> 64);
}

inline u32 getmod(u64 x) { return (u32)x - mod * u32(getdiv(x)); }

struct multi_integer {
  u32 val, ival;
  inline multi_integer() {}
  inline explicit multi_integer(u32 v) {
    val = v;
    ival = ((u64)v << 32) / mod;
  }
  __attribute((always_inline)) u32 operator*(u32 x) const {
    return val * x - u32((u64)x * ival >> 32) * mod;
  }
} wn[N << 1], iwn[N << 1];
u32 lim, shift;

u32 dfta[N], rev[N], w[N];

inline u32 get_index(u32 x, u32 lim) {
  static u32 bak[50];
  int i = 1, lg = std::__lg(lim);
  *bak = x;
  for (int i = 1; i < lg; ++i)
    bak[i] = (u64)bak[i - 1] * bak[i - 1] % mod;
  int res = 0;
  for (int i = lg - 1; i >= 0; --i)
    if ((u64)bak[i] * w[res << i & lim - 1] % mod != 1)
      res += 1 << lg - 1 - i;
  return res ? lim - res : 0;
}

__attribute((always_inline)) u32 norm1(u32 x) { return x >= mod ? x - mod : x; }
static u32 map[4] = {0, mod, mod2, mod + mod2};
inline u32 norm2(u32 x) { return x - map[x >> 30]; }
inline u32 norm2_lazy(u32 x) { return x - map[x >> 30]; }
inline u32 norm2_ex(u32 x) { return x - map[x >> 30]; }
inline u32 pow(u32 coeff, u32 point, u32 ans = 1) {
  for (; point; point >>= 1, coeff = (u64)coeff * coeff % mod)
    if (point & 1)
      ans = (u64)ans * coeff % mod;
  return ans;
}

inline u32 get(u32 x) { return ((u64)x << 32) / mod; }
inline void base_init(u32 len) {
  u32 N = 1;
  for (; N < len;)
    N <<= 1;
  const u32 mid = N >> 1;
  const u32 w = pow(3, mod / N);
  const u32 iw = pow((mod + 1) / 3, mod / N);
  wn[mid] = multi_integer(1);
  iwn[mid] = multi_integer(1);
  for (int i = 1; i < mid; ++i) {
    wn[mid + i] = multi_integer((u64)wn[mid + i - 1].val * w % mod);
    iwn[mid + i] = multi_integer((u64)iwn[mid + i - 1].val * iw % mod);
  }
  for (int i = mid - 1; i >= 0; --i) {
    wn[i] = wn[i << 1];
    iwn[i] = iwn[i << 1];
  }
}
inline void init(u32 len) {
  lim = len;
  shift = std::__lg(lim);
}
inline void safe_init(u32 len) {
  lim = 1, shift = 0;
  for (; lim < len;)
    lim <<= 1, ++shift;
}
inline u32 multi(u32 w, u32 idx) { return wn[idx] * w; }

inline void dft(u32 *coeff) {
#define trans(coeff, point, idx)                                               \
  {                                                                            \
    const u32 A = norm2(coeff + point);                                        \
    point = wn[idx] * (coeff + mod2 - point), coeff = A;                       \
  }
#define trans2(coeff, point)                                                   \
  {                                                                            \
    const u32 A = norm2(coeff + point);                                        \
    point = norm2(coeff + mod2 - point), coeff = A;                            \
  }

  for (int mid = lim >> 1; mid != 4; mid >>= 1) {
    for (int j = 0; j != lim; j += mid + mid) {
      for (int k = 0; k != mid; k += 4) {
        const u32 A0 = wn[mid + k + 0] *
                       (coeff[j + k + 0] + mod2 - coeff[mid + j + k + 0]);
        const u32 A1 = wn[mid + k + 1] *
                       (coeff[j + k + 1] + mod2 - coeff[mid + j + k + 1]);
        const u32 A2 = wn[mid + k + 2] *
                       (coeff[j + k + 2] + mod2 - coeff[mid + j + k + 2]);
        const u32 A3 = wn[mid + k + 3] *
                       (coeff[j + k + 3] + mod2 - coeff[mid + j + k + 3]);
        coeff[j + k + 0] = norm2(coeff[j + k + 0] + coeff[mid + j + k + 0]);
        coeff[j + k + 1] = norm2(coeff[j + k + 1] + coeff[mid + j + k + 1]);
        coeff[j + k + 2] = norm2(coeff[j + k + 2] + coeff[mid + j + k + 2]);
        coeff[j + k + 3] = norm2(coeff[j + k + 3] + coeff[mid + j + k + 3]);
        coeff[mid + j + k + 0] = A0;
        coeff[mid + j + k + 1] = A1;
        coeff[mid + j + k + 2] = A2;
        coeff[mid + j + k + 3] = A3;
      }
    }
  }
  for (int j = 0; j != lim; j += 8) {
    trans2(coeff[j + 0], coeff[j + 4]);
    trans(coeff[j + 1], coeff[j + 5], 5);
    trans(coeff[j + 2], coeff[j + 6], 6);
    trans(coeff[j + 3], coeff[j + 7], 7);

    trans2(coeff[j + 0], coeff[j + 2]);
    trans(coeff[j + 1], coeff[j + 3], 3);
    trans2(coeff[j + 4], coeff[j + 6]);
    trans(coeff[j + 5], coeff[j + 7], 3);

    trans2(coeff[j + 0], coeff[j + 1]);
    trans2(coeff[j + 2], coeff[j + 3]);
    trans2(coeff[j + 4], coeff[j + 5]);
    trans2(coeff[j + 6], coeff[j + 7]);
  }
#undef trans
#undef trans2
}

inline void dft_last(u32 *coeff) {
  int mid = lim >> 1;

  {
    int k = 0, limit = mid / 2 + 1;
    for (; k + 1 != limit; k += 4) {
      coeff[mid + k + 0] = wn[mid + k + 0] * coeff[k + 0];
      coeff[mid + k + 1] = wn[mid + k + 1] * coeff[k + 1];
      coeff[mid + k + 2] = wn[mid + k + 2] * coeff[k + 2];
      coeff[mid + k + 3] = wn[mid + k + 3] * coeff[k + 3];
    }
    coeff[mid + k + 0] = wn[mid + k + 0] * coeff[k + 0];
  }

  mid >>= 1;

  const u32 A0 = wn[mid + 0] * (coeff[0] + mod2 - coeff[mid + 0]);
  coeff[0] = norm2(coeff[0] + coeff[mid + 0]), coeff[mid + 0] = A0;

  {
    int k = 1;
    for (; k + 3 != mid; k += 4) {
      coeff[mid + k + 0] = wn[mid + k + 0] * coeff[k + 0];
      coeff[mid + k + 1] = wn[mid + k + 1] * coeff[k + 1];
      coeff[mid + k + 2] = wn[mid + k + 2] * coeff[k + 2];
      coeff[mid + k + 3] = wn[mid + k + 3] * coeff[k + 3];
    }
    for (; k != mid; k += 1)
      coeff[mid + k + 0] = wn[mid + k + 0] * coeff[k + 0];
  }

  coeff += mid + mid;

  const u32 A1 = wn[mid + 0] * (coeff[0] + mod2 - coeff[mid + 0]);
  coeff[0] = norm2(coeff[0] + coeff[mid + 0]), coeff[mid + 0] = A1;

  {
    int k = 1;
    for (; k + 3 != mid; k += 4) {
      coeff[mid + k + 0] = wn[mid + k + 0] * coeff[k + 0];
      coeff[mid + k + 1] = wn[mid + k + 1] * coeff[k + 1];
      coeff[mid + k + 2] = wn[mid + k + 2] * coeff[k + 2];
      coeff[mid + k + 3] = wn[mid + k + 3] * coeff[k + 3];
    }
    for (; k != mid; k += 1)
      coeff[mid + k + 0] = wn[mid + k + 0] * coeff[k + 0];
  }

  init(lim >> 2);
  dft(coeff - lim);
  dft(coeff);
  dft(coeff + lim);
}

inline u32 div_lim(u32 x) { return (x + u64(-x & (lim - 1)) * mod) >> shift; }
inline void base_idft(u32 *coeff) {
#define trans(coeff, point, idx)                                               \
  {                                                                            \
    u32 A = norm2_lazy(coeff), B = iwn[idx] * point;                           \
    coeff = A + B;                                                             \
    point = A + mod2 - B;                                                      \
  }
#define trans2(coeff, point)                                                   \
  {                                                                            \
    const u32 A = norm2(coeff), B = norm2(point);                              \
    coeff = A + B;                                                             \
    point = A + mod2 - B;                                                      \
  }

  for (int j = 0; j != lim; j += 8) {
    trans2(coeff[j + 0], coeff[j + 1]);
    trans2(coeff[j + 2], coeff[j + 3]);
    trans2(coeff[j + 4], coeff[j + 5]);
    trans2(coeff[j + 6], coeff[j + 7]);

    trans2(coeff[j + 0], coeff[j + 2]);
    trans(coeff[j + 1], coeff[j + 3], 3);
    trans2(coeff[j + 4], coeff[j + 6]);
    trans(coeff[j + 5], coeff[j + 7], 3);

    trans2(coeff[j + 0], coeff[j + 4]);
    trans(coeff[j + 1], coeff[j + 5], 5);
    trans(coeff[j + 2], coeff[j + 6], 6);
    trans(coeff[j + 3], coeff[j + 7], 7);
  }
  for (int mid = 8; mid != lim; mid <<= 1) {
    for (int j = 0; j != lim; j += mid + mid) {
      for (int k = 0; k != mid; k += 4) {
        const u32 A0 = norm2_lazy(coeff[j + k + 0]),
                  B0 = iwn[mid + k + 0] * coeff[mid + j + k + 0];
        const u32 A1 = norm2_lazy(coeff[j + k + 1]),
                  B1 = iwn[mid + k + 1] * coeff[mid + j + k + 1];
        const u32 A2 = norm2_lazy(coeff[j + k + 2]),
                  B2 = iwn[mid + k + 2] * coeff[mid + j + k + 2];
        const u32 A3 = norm2_lazy(coeff[j + k + 3]),
                  B3 = iwn[mid + k + 3] * coeff[mid + j + k + 3];
        coeff[mid + j + k + 0] = A0 + mod2 - B0;
        coeff[mid + j + k + 1] = A1 + mod2 - B1;
        coeff[mid + j + k + 2] = A2 + mod2 - B2;
        coeff[mid + j + k + 3] = A3 + mod2 - B3;
        coeff[j + k + 0] = A0 + B0;
        coeff[j + k + 1] = A1 + B1;
        coeff[j + k + 2] = A2 + B2;
        coeff[j + k + 3] = A3 + B3;
      }
    }
  }
#undef trans
#undef trans2
}
inline void idft_last_copy(u32 *coeff, u32 *res) {

#define trans(coeff, point, idx)                                               \
  {                                                                            \
    u32 A = norm2_lazy(coeff), B = iwn[idx] * point;                           \
    coeff = A + B;                                                             \
    point = A + mod2 - B;                                                      \
  }
#define trans2(coeff, point)                                                   \
  {                                                                            \
    const u32 A = norm2(coeff), B = norm2(point);                              \
    coeff = A + B;                                                             \
    point = A + mod2 - B;                                                      \
  }
  for (int j = 0; j != lim; j += 8) {
    trans2(coeff[j + 0], coeff[j + 1]);
    trans2(coeff[j + 2], coeff[j + 3]);
    trans2(coeff[j + 4], coeff[j + 5]);
    trans2(coeff[j + 6], coeff[j + 7]);

    trans2(coeff[j + 0], coeff[j + 2]);
    trans(coeff[j + 1], coeff[j + 3], 3);
    trans2(coeff[j + 4], coeff[j + 6]);
    trans(coeff[j + 5], coeff[j + 7], 3);

    trans2(coeff[j + 0], coeff[j + 4]);
    trans(coeff[j + 1], coeff[j + 5], 5);
    trans(coeff[j + 2], coeff[j + 6], 6);
    trans(coeff[j + 3], coeff[j + 7], 7);
  }
  for (int mid = 8; mid < lim >> 2; mid <<= 1) {
    for (int j = 0; j != lim; j += mid + mid) {
      for (int k = 0; k != mid; k += 4) {
        const u32 A0 = norm2_lazy(coeff[j + k + 0]),
                  B0 = iwn[mid + k + 0] * coeff[mid + j + k + 0];
        const u32 A1 = norm2_lazy(coeff[j + k + 1]),
                  B1 = iwn[mid + k + 1] * coeff[mid + j + k + 1];
        const u32 A2 = norm2_lazy(coeff[j + k + 2]),
                  B2 = iwn[mid + k + 2] * coeff[mid + j + k + 2];
        const u32 A3 = norm2_lazy(coeff[j + k + 3]),
                  B3 = iwn[mid + k + 3] * coeff[mid + j + k + 3];
        coeff[mid + j + k + 0] = A0 + mod2 - B0;
        coeff[mid + j + k + 1] = A1 + mod2 - B1;
        coeff[mid + j + k + 2] = A2 + mod2 - B2;
        coeff[mid + j + k + 3] = A3 + mod2 - B3;
        coeff[j + k + 0] = A0 + B0;
        coeff[j + k + 1] = A1 + B1;
        coeff[j + k + 2] = A2 + B2;
        coeff[j + k + 3] = A3 + B3;
      }
    }
  }
  int mid = lim >> 2;
  for (int j = 0; j != lim; j += mid + mid) {
    for (int k = 0; k != mid; k += 4) {
      const u32 A0 = norm2_lazy(coeff[j + k + 0]),
                B0 = iwn[mid + k + 0] * coeff[mid + j + k + 0];
      const u32 A1 = norm2_lazy(coeff[j + k + 1]),
                B1 = iwn[mid + k + 1] * coeff[mid + j + k + 1];
      const u32 A2 = norm2_lazy(coeff[j + k + 2]),
                B2 = iwn[mid + k + 2] * coeff[mid + j + k + 2];
      const u32 A3 = norm2_lazy(coeff[j + k + 3]),
                B3 = iwn[mid + k + 3] * coeff[mid + j + k + 3];
      coeff[mid + j + k + 0] = A0 + mod2 - B0;
      coeff[mid + j + k + 1] = A1 + mod2 - B1;
      coeff[mid + j + k + 2] = A2 + mod2 - B2;
      coeff[mid + j + k + 3] = A3 + mod2 - B3;
    }
  }
  res -= mid;
  mid <<= 1;
  for (int k = mid >> 1; k != mid; k += 4) {
    const u32 A0 = norm2_lazy(coeff[k + 0]),
              B0 = iwn[mid + k + 0] * coeff[mid + k + 0];
    const u32 A1 = norm2_lazy(coeff[k + 1]),
              B1 = iwn[mid + k + 1] * coeff[mid + k + 1];
    const u32 A2 = norm2_lazy(coeff[k + 2]),
              B2 = iwn[mid + k + 2] * coeff[mid + k + 2];
    const u32 A3 = norm2_lazy(coeff[k + 3]),
              B3 = iwn[mid + k + 3] * coeff[mid + k + 3];
    res[k + 0] = norm2_ex(A0 + mod2 - B0);
    res[k + 1] = norm2_ex(A1 + mod2 - B1);
    res[k + 2] = norm2_ex(A2 + mod2 - B2);
    res[k + 3] = norm2_ex(A3 + mod2 - B3);
  }
#undef trans
#undef trans2
}

inline void idft(u32 *coeff) {
  base_idft(coeff);
  for (int i = 0; i != lim; ++i)
    coeff[i] = div_lim(coeff[i]);
}

inline void fill(u32 *coeff, const u32 *point, u32 len) {
  memcpy(coeff, point, len << 2);
  memset(coeff + len, 0, (lim - len) << 2);
}
inline void sub(u32 &x) { x = x ? x - 1 : mod - 1; }
static const u32 brute_limit = 32;
u32 coeff[N << 1], point[N << 1], b2[N], b4[N], c[N << 1], d[N << 1];
std::vector<u32> ans(N, 0);
u32 o[9][N];
u32 val[9][N << 2];
u32 sgt[9][N << 2];

u32 power_b[N], ex_b[N];
u32 dft_val[N];

inline void prod(u32 *coeff, const u32 *point, const u32 *c) {
  for (int i = 0; i < lim; i += 4) {
    coeff[i + 0] = getmod((u64)point[i + 0] * c[i + 0]);
    coeff[i + 1] = getmod((u64)point[i + 1] * c[i + 1]);
    coeff[i + 2] = getmod((u64)point[i + 2] * c[i + 2]);
    coeff[i + 3] = getmod((u64)point[i + 3] * c[i + 3]);
  }
}
inline void solve(u32 m, u32 dep, u32 l, u32 r, bool good) {

  if (l >= m) {
    std::fill(val[dep] + l * 4, val[dep] + r * 4, 1);
    return;
  }

  u32 n = r - l;

  static u32 t[N];
  if (n < brute_limit) {
    u32 *x = o[dep] + l;
    for (int i = l; i < r; i += 4) {
      const u32 v0 = point[i + 0] ? mod - point[i + 0] : 0;
      const u32 v1 = point[i + 1] ? mod - point[i + 1] : 0;
      const u32 v2 = point[i + 2] ? mod - point[i + 2] : 0;
      const u32 v3 = point[i + 3] ? mod - point[i + 3] : 0;

      const u32 v01 = (u64)v0 * v1 % mod;
      const u32 v23 = (u64)v2 * v3 % mod;

      const u32 a4 = (u64)v01 * v23 % mod;
      const u32 a3 = ((u64)v01 * (v2 + v3) + (u64)v23 * (v0 + v1)) % mod;
      const u32 a2 = (v01 + v23 + (u64)(v0 + v1) * (v2 + v3)) % mod;
      const u32 a1 = (v0 + v1 + v2 + v3) % mod;

      for (int j = i - l + 3; j > 3; --j) {
        x[j] = ((u64)x[j - 4] * a4 + (u64)x[j - 3] * a3 + (u64)x[j - 2] * a2 +
                (u64)x[j - 1] * a1 + x[j]) %
               mod;
      }
      x[3] =
          (a4 + (u64)x[0] * a3 + (u64)x[1] * a2 + (u64)x[2] * a1 + x[3]) % mod;
      x[2] = (a3 + (u64)x[0] * a2 + (u64)x[1] * a1 + x[2]) % mod;
      x[1] = (a2 + (u64)x[0] * a1 + x[1]) % mod;
      x[0] = norm1(x[0] + a1);
    }

  } else {
    u32 mid0 = (l * 3 + r * 1) >> 2;
    u32 mid1 = (l * 2 + r * 2) >> 2;
    u32 mid2 = (l * 1 + r * 3) >> 2;
    solve(m, dep + 1, l, mid0, 1);
    solve(m, dep + 1, mid0, mid1, 1);
    solve(m, dep + 1, mid1, mid2, 1);
    solve(m, dep + 1, mid2, r, 1);
    init(n);
    u32 *da = val[dep + 1] + l * 4;
    u32 *db = val[dep + 1] + mid0 * 4;
    u32 *dc = val[dep + 1] + mid1 * 4;
    u32 *dd = val[dep + 1] + mid2 * 4;
    if (mid0 >= m) {
      memcpy(dft_val, da, lim << 2);
      memcpy(o[dep] + l, o[dep + 1] + l, lim);
    } else {
      prod(sgt[dep] + l + l, da, db);
      prod(sgt[dep] + mid1 + mid1, dc, dd);
      prod(dft_val, sgt[dep] + l + l, sgt[dep] + mid1 + mid1);
      memcpy(t, dft_val, lim << 2);
      idft(t);
      memcpy(o[dep] + l, t + 1, (lim - 1) << 2);
      sub(o[dep][l + n - 1] = t[0]);
    }
  }

  if (good) {
    init(n << 2);
    memcpy(val[dep] + l * 4 + 1, o[dep] + l, n << 2);
    val[dep][l * 4] = 1;
    if (n >= brute_limit) {
      dft_last(val[dep] + l * 4);
      memcpy(val[dep] + l * 4, dft_val, n << 2);
    } else {
      dft(val[dep] + l * 4);
    }
  }
}
inline void getans(u32 m, u32 dep, u32 l, u32 r, u32 *coeff, u32 *cur,
                   bool good, u32 inv) {

  if (l >= m)
    return;

  u32 n = r - l;
  if (n < brute_limit) {
    static u32 g[N], B[N];
    memcpy(B + 1, o[dep] + l, (n - 1) << 2);
    *B = 1;
    for (int i = 0; i < n; ++i) {
      coeff[i] = norm1(coeff[i]);
      u64 sum = 0;
      int j = 0;
      for (; j + 3 <= i; j += 4) {
        sum += (u64)coeff[j + 0] * B[i - j - 0] +
               (u64)coeff[j + 1] * B[i - j - 1] +
               (u64)coeff[j + 2] * B[i - j - 2] +
               (u64)coeff[j + 3] * B[i - j - 3];
      }
      for (; j <= i; ++j)
        sum += (u64)coeff[j] * B[i - j];
      g[i] = sum % mod;
    }
    for (int i = l; i < r; ++i) {
      u32 &x = ans[i];
      const u32 b1 = point[i];
      const u32 b2 = MPE::b2[i];
      const u32 b3 = (u64)b2 * b1 % mod;
      const u32 b4 = MPE::b4[i];
      x = ((u64)x * b4 + (u64)g[0 + 0] * b3 + (u64)g[0 + 1] * b2 +
           (u64)g[0 + 2] * b1 + g[0 + 3]) %
          mod;
      x = ((u64)x * b4 + (u64)g[4 + 0] * b3 + (u64)g[4 + 1] * b2 +
           (u64)g[4 + 2] * b1 + g[4 + 3]) %
          mod;
      if (n != 8) {
        x = ((u64)x * b4 + (u64)g[8 + 0] * b3 + (u64)g[8 + 1] * b2 +
             (u64)g[8 + 2] * b1 + g[8 + 3]) %
            mod;
        x = ((u64)x * b4 + (u64)g[12 + 0] * b3 + (u64)g[12 + 1] * b2 +
             (u64)g[12 + 2] * b1 + g[12 + 3]) %
            mod;
      }
      x = (u64)x * inv % mod;
    }
    return;
  }
  static u32 c[N];
  u32 mid0 = (l * 3 + r * 1) >> 2;

  if (mid0 >= m) {
    for (int i = 0; i != n / 4; ++i)
      coeff[n + i] = coeff[n / 4 * 3 + i];
    getans(m, dep + 1, l, mid0, coeff + n, cur + n, 1, inv);
    return;
  }

  u32 mid1 = (l * 2 + r * 2) >> 2;
  u32 mid2 = (l * 1 + r * 3) >> 2;

  u32 *da = val[dep + 1] + l * 4;
  u32 *db = val[dep + 1] + mid0 * 4;
  u32 *dc = val[dep + 1] + mid1 * 4;
  u32 *dd = val[dep + 1] + mid2 * 4;
  u32 *dab = sgt[dep] + l * 2;
  u32 *dcd = sgt[dep] + mid1 * 2;

  init(n);
  dft(coeff);

  prod(cur, coeff, dcd);
  prod(c, cur, db);

  inv = div_lim(inv);
  idft_last_copy(c, coeff + n);
  getans(m, dep + 1, l, mid0, coeff + n, cur + n, 1, inv);

  if (mid0 >= m)
    return;
  init(n);
  prod(c, cur, da);
  idft_last_copy(c, coeff + n);
  getans(m, dep + 1, mid0, mid1, coeff + n, cur + n, 1, inv);

  if (mid1 >= m)
    return;
  init(n);
  prod(cur, coeff, dab);
  prod(c, cur, dd);
  idft_last_copy(c, coeff + n);
  getans(m, dep + 1, mid1, mid2, coeff + n, cur + n, 1, inv);

  if (mid2 >= m)
    return;
  init(n);
  prod(c, cur, dc);
  idft_last_copy(c, coeff + n);
  getans(m, dep + 1, mid2, r, coeff + n, cur + n, 1, inv);
}

inline void init_inv(const u32 *coeff, u32 *res, int n) {
  static u32 pre[N];
  *pre = *coeff;
  for (int i = 1; i != n; ++i) {
    pre[i] = (u64)pre[i - 1] * coeff[i] % mod;
  }
  u32 &all_inv = res[0] = pow(pre[n - 1], mod - 2);
  for (int i = n - 1; i; --i) {
    res[i] = (u64)all_inv * pre[i - 1] % mod;
    all_inv = (u64)all_inv * coeff[i] % mod;
  }
}

inline void naive(u32 n, u32 m) {
  for (u32 i = m - 1; ~i; --i) {
    u32 x = point[i], &ret = ans[i] = 0;
    for (int i = n - 1; ~i; --i) {
      ret = ((u64)ret * x + coeff[i]) % mod;
    }
  }
}

std::vector<u32> multi_point_evaluation(u32 n, u32 m, u32 M) {
  if ((u64)n * m < 500000) {
    naive(n, m);
  } else {
    init(M);
    fill(dfta, coeff, n);
    dft(dfta);
    for (int i = 1; i < lim; ++i) {
      rev[i] = rev[i >> 1] >> 1 | i % 2 * lim / 2;
    }
    const int GG = pow(3, mod / lim);
    for (int i = 0, multi = 1; i < lim; ++i) {
      coeff[rev[lim - i]] = (u64)dfta[rev[i]] * multi % mod;
      w[i] = multi;
      multi = (u64)multi * GG % mod;
    }
    *coeff = *dfta;

    u32 res = 0;
    for (int i = 0; i < m; ++i) {
      power_b[i] = point[i];
    }
    static u32 og[N];
    init(M), lim = m;

    prod(power_b, power_b, power_b);
    memcpy(b2, power_b, m << 2);
    prod(power_b, power_b, power_b);
    memcpy(b4, power_b, m << 2);

    for (int i = 2; i < shift; ++i) {
      prod(power_b, power_b, power_b);
    }
    init(M);
    for (int i = m; i < lim; ++i) {
      power_b[i] = 1;
    }
    for (int i = 0; i < lim; ++i) {
      power_b[i] = norm1(norm2(power_b[i] + mod - 1));
      if (power_b[i] == 0) {
        og[i] = point[i];
        point[i] = 3;
      }
    }

    solve(m, 0, 0, M, 0);
    init(M);
    init_inv(dft_val, d, M);
    prod(coeff, d, coeff);
    idft(coeff);
    for (u32 i = 0; i < M; ++i) {
      coeff[i] = norm2(coeff[i]);
    }
    static u32 cur[N << 1];
    getans(m, 0, 0, M, coeff, cur, 0, 1);
    for (int i = 0; i < m; ++i) {
      if (power_b[i]) {
        ans[i] = (u64)ans[i] * (mod - power_b[i]) % mod;
      } else {
        ans[i] = norm1(dfta[rev[get_index(og[i], M)]]);
      }
    }
  }
  return ans;
}
} // namespace MPE

int main() {
  u32 n, m, M;
  cin >> n >> m;

  for (M = 1; M < n || M < m;) {
    M <<= 1;
  }
  const int max_ask = 1 << 18;
  if (M > max_ask) {
    M = max_ask;
  }

  MPE::mod = 998244353;
  MPE::base_init(M);

  for (u32 i = 0; i < n; ++i) {
    cin >> MPE::coeff[i];
  }

  for (u32 i = 0; i < m; ++i) {
    cin >> MPE::point[i];
  }

  auto ans = MPE::multi_point_evaluation(n, m, M);
  for (u32 i = m - 1; ~i; --i) {
    cout << '\n' << ans[i];
  }
}

// https://judge.yosupo.jp/submission/16397
