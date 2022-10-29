#include <bits/stdc++.h>
#define lg2 std::__lg
#define ctz __builtin_ctz

typedef unsigned short u16;
typedef unsigned int u32;
typedef unsigned long long u64;

struct istream {
  static const int size = 1 << 26;
  static const u32 b = 0x30303030;
  char buf[size], *vin;
  istream() : vin(buf - 1) { fread(buf, 1, size, stdin); }
  inline istream &operator>>(u64 &x) {
    x = *++vin & 15, ++vin;
    u32 *&idx = (u32 *&)vin;
    for (; (*idx & b) == b; ++idx) {
      *idx ^= b;
      *idx = (*idx >> 8 & 0x00ff00ff) + (*idx & 0x00ff00ff) * 10;
      x = x * 10000 + (*idx & 65535) * 100 + (*idx >> 16);
    }
    for (; isdigit(*vin); ++vin)
      x = x * 10 + (*vin & 15);
    return *this;
  }
} cin;

struct ostream {
  static const int size = 1 << 27;
  char buf[size], *vout;
  u32 map[10000];
  ostream() {
    for (int i = 0; i < 10000; ++i)
      map[i] = i % 10 + 48 << 24 | i / 10 % 10 + 48 << 16 |
               i / 100 % 10 + 48 << 8 | i / 1000 + 48;
    vout = buf + size;
  }
  ~ostream() { fwrite(vout, 1, buf + size - vout, stdout); }
  inline ostream &operator<<(u64 x) {
    for (; x >= 10000; x /= 10000)
      *--(u32 *&)vout = map[x % 10000];
    do
      *--vout = x % 10 + 48;
    while (x /= 10);
    return *this;
  }
  inline ostream &operator<<(char x) { return *--vout = x, *this; }
  inline ostream &operator<<(const char *s) {
    const int L = strlen(s);
    return memcpy(vout -= L, s, L), *this;
  }
} cout;

namespace nimbers {
constexpr u32 n2f[16] = {0x0001u, 0x8ff5u, 0x6cbfu, 0xa5beu, 0x761au, 0x8238u,
                         0x4f08u, 0x95acu, 0xf340u, 0x1336u, 0x7d5eu, 0x86e7u,
                         0x3a47u, 0xe796u, 0xb7c3u, 0x0008u},
              f2n[16] = {0x0001u, 0x2827u, 0x392bu, 0x8000u, 0x20fdu, 0x4d1du,
                         0xde4au, 0x0a17u, 0x3464u, 0xe3a9u, 0x6d8du, 0x34bcu,
                         0xa921u, 0xa173u, 0x0ebcu, 0x0e69u};
inline u32 nimber2field(u32 x) {
  u32 y = 0;
  for (; x; x &= x - 1)
    y ^= n2f[ctz(x)];
  return y;
}
inline u32 field2nimber(u32 x) {
  u32 y = 0;
  for (; x; x &= x - 1)
    y ^= f2n[ctz(x)];
  return y;
}
inline u32 __builtin_double(u32 x) {
  return x << 1 ^ (x < 0x8000u ? 0 : 0x1681fu);
}

int ln[65536];
u16 buf_[262200], *exp = buf_ + 131082, *Hexp = exp + 3, *H2exp = exp + 6;

inline void init() {
  int i;
  *ln = -65541;
  for (*exp = i = 1; i < 65535; ++i)
    exp[i] = __builtin_double(exp[i - 1]);
  for (i = 1; i < 65535; ++i)
    exp[i] = field2nimber(exp[i]), ln[exp[i]] = i;
  memcpy(exp + 65535, exp, 131070);
  memcpy(exp + 131070, exp, 10);
}

inline u16 product(u16 A, u16 B) { return exp[ln[A] + ln[B]]; }
inline u16 H(u16 A) { return Hexp[ln[A]]; }
inline u16 H2(u16 A) { return H2exp[ln[A]]; }
inline u16 Hproduct(u16 A, u16 B) { return Hexp[ln[A] + ln[B]]; }

inline u32 product(u32 A, u32 B) {
  u16 a = A & 65535, b = B & 65535, c = A >> 16, d = B >> 16, e = product(a, b);
  return u32(product(u16(a ^ c), u16(b ^ d)) ^ e) << 16 | (Hproduct(c, d) ^ e);
}

inline u32 H(u32 A) {
  u16 a = A & 65535, b = A >> 16;
  return H(u16(a ^ b)) << 16 | H2(b);
}

inline u64 product(u64 A, u64 B) {
  u32 a = A & UINT_MAX, b = B & UINT_MAX, c = A >> 32, d = B >> 32,
      e = product(a, b);
  return u64(product(a ^ c, b ^ d) ^ e) << 32 | (H(product(c, d)) ^ e);
}
}

u64 ans[1000001];

int main() {
  int i, T;
  u64 a, b;
  cin >> a, T = a, nimbers::init();
  for (i = 0; i < T; ++i)
    cin >> a >> b, ans[i] = nimbers::product(a, b);
  for (i = T - 1; i >= 0; --i)
    cout << '\n' << ans[i];
  return 0;
}
// https://judge.yosupo.jp/problem/nim_product_64
