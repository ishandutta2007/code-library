#pragma GCC optimize("O3")
// #pragma GCC target ("avx")
#pragma GCC target("sse4.2") // SPOJ, codechef, CodeIQ

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#define _rep(_1, _2, _3, _4, name, ...) name
#define rep2(i, n) rep3(i, 0, n)
#define rep3(i, a, b) rep4(i, a, b, 1)
#define rep4(i, a, b, c) for (int i = int(a); i < int(b); i += int(c))
#define rep(...) _rep(__VA_ARGS__, rep4, rep3, rep2, _)(__VA_ARGS__)

using namespace std;

using i64 = long long;
using u8 = unsigned char;
using u32 = unsigned;
using u64 = unsigned long long;
using f80 = long double;

const u32 MOD = 786433;
const u32 G = 10;

const u32 N_MAX = 250000;
u32 X[N_MAX + 10];
u32 A[MOD - 1];
u32 B[MOD - 1];

u32 ind[MOD];
u32 pows[MOD];

u32 pow_mod(u32 a, u32 e) {
  u32 ret = 1;
  for (; e; e >>= 1, a = u64(a) * a % MOD) {
    if (e & 1)
      ret = u64(ret) * a % MOD;
  }
  return ret;
}

u32 add_mod(u32 a, u32 b) {
  a += b;
  return a >= MOD ? a - MOD : a;
}
u32 sub_mod(u32 a, u32 b) {
  a -= b;
  return int(a) < 0 ? a + MOD : a;
}

void ntt(u32 N) {
  const u32 n = MOD - 1;
  const u32 logn = 31 - __builtin_clz(n / 3);
  u32 omega = pow_mod(G, n / 3);
  u32 omega2 = u64(omega) * omega % MOD;
  rep(i, N, n) A[i] = 0;
  rep(i, n) {
    int k = i;
    int r = 0;
    rep(j, logn) r = r * 2 + (k & 1), k >>= 1;
    r = r * 3 + k;
    B[r] = A[i];
  }
  rep(i, 0, n, 3) {
    u32 a = B[i + 0];
    u32 b = B[i + 1];
    u32 c = B[i + 2];
    B[i + 0] = (a + b + c) % MOD;
    B[i + 1] = (a + u64(omega) * b + u64(omega2) * c) % MOD;
    B[i + 2] = (a + u64(omega2) * b + u64(omega) * c) % MOD;
  }
  rep(lg, 1, logn + 1) {
    u32 m = 3 << lg;
    u32 mh = m >> 1;
    u32 dw = pow_mod(G, n / m);
    u32 w = 1;
    rep(i, mh) {
      rep(j, i, n, m) {
        u32 a = B[j];
        u32 b = u64(B[j + mh]) * w % MOD;
        B[j] = add_mod(a, b);
        B[j + mh] = sub_mod(a, b);
      }
      w = u64(w) * dw % MOD;
    }
  }
}

u32 eval(u32 N, u32 x) {
  u32 ret = 0;
  rep(i, N) ret = (u64(ret) * x + A[N - 1 - i]) % MOD;
  return ret;
}

void solve() {
  int N;
  while (~scanf("%d", &N)) {
    rep(i, N + 1) scanf("%u", &A[i]);
    int Q;
    scanf("%d", &Q);
    rep(i, Q) scanf("%u", &X[i]);
    u32 pw = 1;
    rep(i, MOD - 1) {
      ind[pw] = i;
      pows[i] = pw;
      pw = u64(pw) * G % MOD;
    }
    ntt(N + 1);
    rep(i, Q) {
      u32 ans = A[0];
      if (X[i])
        ans = B[ind[X[i]]];
      printf("%u\n", ans);
    }
  }
}

int main() {
  clock_t beg = clock();
  solve();
  clock_t end = clock();
  fprintf(stderr, "%.3f sec\n", double(end - beg) / CLOCKS_PER_SEC);
  return 0;
}