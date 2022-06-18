#include <bits/stdc++.h>
// credits nurbijoy
using namespace std;

using f64 = double;
using f128 = __float128;
using i64 = long long;
using i128 = __int128;
using u64 = unsigned long long;
using u128 = __uint128_t;

template <class T> void read(T &x) {
  char c;
  for (c = getchar(); c < '0' || c > '9'; c = getchar())
    ;
  for (x = 0; c <= '9' && c >= '0'; c = getchar())
    x = x * 10 + (c & 15);
}
template <class I> string to_str(I x) {
  if (!x)
    return "0";
  if (x < 0)
    return "-" + to_str(-x);
  string res = "";
  while (x)
    res += x % 10 + '0', x /= 10;
  reverse(res.begin(), res.end());
  return res;
}
ostream &operator<<(ostream &cout, i128 x) { return (cout << to_str(x)); }
ostream &operator<<(ostream &cout, u128 x) { return (cout << to_str(x)); }

namespace Solver {

u64 pow_64(u64 a, u64 t, u64 mod) {
  u64 r = 1;
  for (; t; t >>= 1, a = u128(a) * a % mod)
    if (t & 1)
      r = u128(r) * a % mod;
  return r;
}

bool check_prime(u64 n) {
  static const int jp[] = {2, 3, 5, 7, 11, 13, 17, 19};
  if (n == 1)
    return false;
  for (int p : jp)
    if (n % p == 0)
      return n == p;
  u64 r = n - 1, x, y;
  int e = 0;
  while (~r & 1)
    r >>= 1, ++e;
  for (int p : jp) {
    x = pow_64(p, r, n);
    for (int t = 0; t < e && x > 1; ++t) {
      y = (i128)x * x % n;
      if (y == 1 && x != n - 1)
        return false;
      x = y;
    }
    if (x != 1)
      return false;
  }
  return true;
}
u64 next_prime(u64 n) {
  n += ~n & 1;
  while (!check_prime(n += 2))
    ;
  return n;
}

u128 ctz(u128 x) {
  return u64(x) ? __builtin_ctzll(x) : __builtin_ctzll(x >> 64) + 64;
}
u128 gcd(u128 a, u128 b) {
  if (!a || !b)
    return a | b;
  int shift = ctz(a | b);
  for (b >>= ctz(b); a; a -= b)
    if ((a >>= ctz(a)) < b)
      swap(a, b);
  return b << shift;
}

struct u256 {
  u128 lo, hi;
  u256() {}
  u256(u128 lo, u128 hi) : lo(lo), hi(hi) {}
  static u256 mul128(u128 a, u128 b) {
    u64 a_hi = a >> 64, a_lo = u64(a);
    u64 b_hi = b >> 64, b_lo = u64(b);
    u128 p01 = u128(a_lo) * b_lo;
    u128 p12 = u128(a_hi) * b_lo + u64(p01 >> 64);
    u64 t_hi = p12 >> 64, t_lo = p12;
    p12 = u128(a_lo) * b_hi + t_lo;
    u128 p23 = u128(a_hi) * b_hi + u64(p12 >> 64) + t_hi;
    return u256(u64(p01) | (p12 << 64), p23);
  }
};

struct Mont {
  u128 mod, inv, r2;
  Mont(u128 n) : mod(n) {
    assert(n & 1);
    inv = n;
    for (int i = 0; i < 6; ++i)
      inv *= 2 - n * inv;
    r2 = -n % n;
    for (int i = 0; i < 4; ++i)
      if ((r2 <<= 1) >= mod)
        r2 -= mod;
    for (int i = 0; i < 5; ++i)
      r2 = mul(r2, r2);
  }
  u128 reduce(u256 x) const {
    u128 y = x.hi - u256::mul128(x.lo * inv, mod).hi;
    return i128(y) < 0 ? y + mod : y;
  }
  u128 reduce(u128 x) const { return reduce(u256(x, 0)); }
  u128 init(u128 n) const { return reduce(u256::mul128(n, r2)); }
  u128 mul(u128 a, u128 b) const { return reduce(u256::mul128(a, b)); }
} mont(1);

i128 N;
i128 add(i128 a, i128 b) { return (a += b) <= N ? a : a - N; }
i128 sub(i128 a, i128 b) { return (a -= b) < 0 ? a + N : a; }
i128 mul(i128 a, i128 b) { return mont.mul(a, b); }
void exgcd(i128 a, i128 b, i128 &x, i128 &y) {
  if (b == 0)
    return (void)(x = 1, y = 0);
  exgcd(b, a % b, y, x), y -= a / b * x;
}
i128 inv128(i128 t) {
  i128 q = N, x, y;
  exgcd(t, q, x, y);
  return (x % N + N) % N;
}

const int L = 250;
int prime[L], pcnt, root[L];
u128 _mont_prime[L];
f64 logp[L];

struct word {
  bitset<L> mask;
  u128 lhs, rhs;
  int bit;
  word() { mask = lhs = rhs = 0, bit = L; }
  void merge(word &x) {
    lhs = mul(lhs, x.lhs), rhs = mul(rhs, x.rhs);
    bitset<L> cons = mask & x.mask;
    for (int pos = cons._Find_first(); pos < L; pos = cons._Find_next(pos))
      rhs = mul(rhs, _mont_prime[pos]);
    mask ^= x.mask, bit = mask._Find_first();
  }
};
vector<word> smooth;
unordered_map<i64, word> data;
u128 insert(word &o) {
  int x;
  for (x = 0; x < smooth.size(); ++x) {
    if (smooth[x].bit > o.bit)
      break;
    if (smooth[x].bit == o.bit)
      o.merge(smooth[x]);
  }
  if (o.bit == L) {
    u128 g = gcd(o.lhs + o.rhs, N);
    if (g != 1 && g != N) {
      return g;
      // cout<<"g="<<g<<endl;
      // if (g > N / g)
      //     g = N / g;
      // cout << g << ' ' << N / g << endl;
      // cerr << clock() << endl;
      // exit(0);
    }
  } else
    smooth.insert(smooth.begin() + x, o);
  return 0;
}
u128 try_insert(i128 x, i128 y) {
  while (x < 0)
    x += N;
  word ins;
  ins.lhs = x, ins.rhs = 1;
  for (int k = 1; k <= pcnt; ++k) {
    int cnt = 0;
    while (y % prime[k] == 0)
      y /= prime[k], ++cnt;
    if (cnt)
      for (ins.mask[k] = cnt & 1, cnt >>= 1; cnt; --cnt)
        ins.rhs *= prime[k];
  }
  ins.lhs = mont.init(ins.lhs % N);
  ins.rhs = mont.init(ins.rhs % N);
  ins.bit = ins.mask._Find_first();
  if (y != 1) {
    if (data.count(y)) {
      ins.merge(data[y]);
      ins.rhs = mul(ins.rhs, mont.init(add(y, N)));
      y = 1;
    } else
      data[y] = ins;
  }
  if (y == 1) {
    // cout<<"ins={"<<ins.lhs<<" "<<ins.rhs<<" "<<ins.bit<<  endl;
    return insert(ins);
  }
  return 0;
}

void init() {
  int i, j;
  int B = pow(log(1. * N), 2) * .6;
  vector<char> mark((B >> 1) + 5);
  for (i = 3; i * i <= B; i += 2)
    if (!mark[i >> 1])
      for (j = i * i; j <= B; j += i << 1)
        mark[j >> 1] = true;
  auto append = [&](int p) {
    prime[++pcnt] = p, _mont_prime[pcnt] = mont.init(p), logp[pcnt] = log(p);
  };
  for (append(2), i = 3; i <= B; i += 2)
    if (!mark[i >> 1])
      if (pow_64(N % i, i >> 1, i) == 1)
        append(i);
  for (i = 1; i <= pcnt; ++i)
    for (j = N % prime[i]; root[i] * root[i] % prime[i] != j; ++root[i])
      ;
  // cout<<"pcnt="<<pcnt<<endl;
  // for(int i=0;i<pcnt;i++)cout<<root[i]<<" ";cout<<endl;
  // for(int i=0;i<pcnt;i++)cout<<prime[i]<<" ";cout<<endl;
  // for(int i=0;i<pcnt;i++)cout<<_mont_prime[i]<<" ";cout<<endl;
  // for(int i=0;i<pcnt;i++)cout<<logp[i]<<" ";cout<<endl;
}
u128 main(i128 n) {
  N = n;
  mont = Mont(N);
  init();
  int M = pcnt * 50;
  i64 D = sqrt(sqrt(2. * N) / M);
  f64 bound = log(M * sqrt(.5 * N)) * 0.75;
  // cout<<"M="<<M<<endl;
  // cout<<"D="<<D<<endl;
  // cout<<"bound="<<bound<<endl;
  while (true) {
    do {
      D = next_prime(D);
      // cout<<"D(in)="<<D<<endl;
    } while ((D & 3) == 1);
    // cout<<"D(exit)="<<D<<endl;
    u64 y0 = pow_64(u64(N % D), (D + 1) >> 2, D);
    if ((i128(y0) * y0 - N) % D)
      continue;
    u64 y1 = u64((N - y0 * y0) / D % D + D) % D;
    y1 = u128(y1) * (D / 2 + 1) % D * pow_64(y0, D - 2, D) % D;
    i128 A = i128(D) * D;
    i128 B = u128(y1) * D + y0;
    i128 C = (B * B - N) / A;
    i128 E = mul(mont.init(B), inv128(D));
    vector<f64> pos(M), neg(M);
    for (int x = 1; x <= pcnt; ++x) {
      int p = prime[x], s = A % p, a = A % p, inv_a = 1;
      if (!a)
        continue;
      while (a > 1)
        inv_a = (p - inv_a) * (p / a) % p, a = p % a;
      int u = (p - B % p) * inv_a % p, v = root[x] * inv_a % p;
      int r1 = (u + v) % p, r2 = (u - v + p) % p;
      for (int su = 0; su < M; su += p) {
        if (su + r1 < M)
          pos[su + r1] += logp[x];
        if (su + r2 < M)
          pos[su + r2] += logp[x];
        if (su > 0)
          neg[su - r1] += logp[x], neg[su - r2] += logp[x];
      }
    }
    for (int x = 0; x < M; ++x) {
      u128 prime1;
      if (pos[x] > bound)
        prime1 = try_insert(E + D * x, A * x * x + 2 * B * x + C);
      if (neg[x] > bound)
        prime1 = try_insert(E - D * x, A * x * x - 2 * B * x + C);
      if (prime1)
        return prime1;
    }
  }
}
} // namespace Solver

int main() {
  i128 n;
  read(n);
  u128 prime1 = Solver::main(n);
  cout << "prime1=" << prime1 << " prime2=" << n / prime1 << endl;
}