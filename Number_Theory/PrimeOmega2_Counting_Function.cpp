#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using int64 = long long;
using uint64 = long long;
#define MOD 1000000087
// https://oeis.org/A036352
// https://oeis.org/A066265
// https://oeis.org/A001358

vector<int> primes;

inline uint64 mul_mod(uint64 a, uint64 b, uint64 mod) {
  assert(0 <= a && a < mod);
  assert(0 <= b && b < mod);
  if (mod < int(1e9))
    return a * b % mod;
  uint64 k = (uint64)((long double)a * b / mod);
  uint64 res = a * b - k * mod;
  if ((int64)res < 0)
    res += mod;
  return res;
}

int64 pow_mod(int64 a, int64 n, int64 m) {
  int64 res = 1;
  for (a %= m; n; n >>= 1) {
    if (n & 1)
      res = mul_mod(res, a, m);
    a = mul_mod(a, a, m);
  }
  return res;
}

const int MAXGLOBAL_ARRAYSIZE = 22000000;
vector<bool> is_prime_bool;
static int prime_cnt[MAXGLOBAL_ARRAYSIZE];

// credit: min_25
// takes 0.6s for n = 1e9
vector<int> sieve(const int N, const int Q = 17, const int L = 1 << 15) {
  static const int rs[] = {1, 7, 11, 13, 17, 19, 23, 29};
  struct P {
    P(int p) : p(p) {}
    int p;
    int pos[8];
  };
  auto approx_prime_count = [](const int N) -> int {
    return N > 60184 ? N / (log(N) - 1.1) : max(1., N / (log(N) - 1.11)) + 1;
  };

  const int v = sqrt(N), vv = sqrt(v);
  vector<bool> isp(v + 1, true);
  for (int i = 2; i <= vv; ++i)
    if (isp[i]) {
      for (int j = i * i; j <= v; j += i)
        isp[j] = false;
    }

  const int rsize = approx_prime_count(N + 30);
  vector<int> primes = {2, 3, 5};
  int psize = 3;
  primes.resize(rsize);

  vector<P> sprimes;
  size_t pbeg = 0;
  int prod = 1;
  for (int p = 7; p <= v; ++p) {
    if (!isp[p])
      continue;
    if (p <= Q)
      prod *= p, ++pbeg, primes[psize++] = p;
    auto pp = P(p);
    for (int t = 0; t < 8; ++t) {
      int j = (p <= Q) ? p : p * p;
      while (j % 30 != rs[t])
        j += p << 1;
      pp.pos[t] = j / 30;
    }
    sprimes.push_back(pp);
  }

  vector<unsigned char> pre(prod, 0xFF);
  for (size_t pi = 0; pi < pbeg; ++pi) {
    auto pp = sprimes[pi];
    const int p = pp.p;
    for (int t = 0; t < 8; ++t) {
      const unsigned char m = ~(1 << t);
      for (int i = pp.pos[t]; i < prod; i += p)
        pre[i] &= m;
    }
  }

  const int block_size = (L + prod - 1) / prod * prod;
  vector<unsigned char> block(block_size);
  unsigned char *pblock = block.data();
  const int M = (N + 29) / 30;

  for (int beg = 0; beg < M; beg += block_size, pblock -= block_size) {
    int end = min(M, beg + block_size);
    for (int i = beg; i < end; i += prod) {
      copy(pre.begin(), pre.end(), pblock + i);
    }
    if (beg == 0)
      pblock[0] &= 0xFE;
    for (size_t pi = pbeg; pi < sprimes.size(); ++pi) {
      auto &pp = sprimes[pi];
      const int p = pp.p;
      for (int t = 0; t < 8; ++t) {
        int i = pp.pos[t];
        const unsigned char m = ~(1 << t);
        for (; i < end; i += p)
          pblock[i] &= m;
        pp.pos[t] = i;
      }
    }
    for (int i = beg; i < end; ++i) {
      for (int m = pblock[i]; m > 0; m &= m - 1) {
        primes[psize++] = i * 30 + rs[__builtin_ctz(m)];
      }
    }
  }
  assert(psize <= rsize);
  while (psize > 0 && primes[psize - 1] > N)
    --psize;
  primes.resize(psize);
  int j = 0, maxj = min(MAXGLOBAL_ARRAYSIZE, N);
  prime_cnt[j] = 0;
  for (int i = 0; i < primes.size() and j < maxj; i++) {
    while (j != primes[i]) {
      is_prime_bool.push_back(0);
      if (j > 0)
        prime_cnt[j] = prime_cnt[j - 1];
      j++;
    }
    is_prime_bool.push_back(1);
    prime_cnt[j] = prime_cnt[j - 1] + 1;
    j++;
  }
  while (j < maxj) {
    prime_cnt[j] = prime_cnt[j - 1];
    j++;
  }
  return primes;
}

using ll = long long;

int isqrt(ll n) { return sqrtl(n); }

ll prime_pi(const ll N) {
  if (N <= 1)
    return 0;
  if (N == 2)
    return 1;
  const int v = isqrt(N);
  int s = (v + 1) / 2;
  vector<int> smalls(s);
  for (int i = 1; i < s; i++)
    smalls[i] = i;
  vector<int> roughs(s);
  for (int i = 0; i < s; i++)
    roughs[i] = 2 * i + 1;
  vector<ll> larges(s);
  for (int i = 0; i < s; i++)
    larges[i] = (N / (2 * i + 1) - 1) / 2;
  vector<bool> skip(v + 1);
  const auto divide = [](ll n, ll d) -> int { return (double)n / d; };
  const auto half = [](int n) -> int { return (n - 1) >> 1; };
  int pc = 0;
  for (int p = 3; p <= v; p += 2)
    if (!skip[p]) {
      int q = p * p;
      if ((ll)q * q > N)
        break;
      skip[p] = true;
      for (int i = q; i <= v; i += 2 * p)
        skip[i] = true;
      int ns = 0;
      for (int k = 0; k < s; k++) {
        int i = roughs[k];
        if (skip[i])
          continue;
        ll d = (ll)i * p;
        larges[ns] = larges[k] -
                     (d <= v ? larges[smalls[d >> 1] - pc]
                             : smalls[half(divide(N, d))]) +
                     pc;
        roughs[ns++] = i;
      }
      s = ns;
      for (int i = half(v), j = ((v / p) - 1) | 1; j >= p; j -= 2) {
        int c = smalls[j >> 1] - pc;
        for (int e = (j * p) >> 1; i >= e; i--)
          smalls[i] -= c;
      }
      pc++;
    }
  larges[0] += (ll)(s + 2 * (pc - 1)) * (s - 1) / 2;
  for (int k = 1; k < s; k++)
    larges[0] -= larges[k];
  for (int l = 1; l < s; l++) {
    ll q = roughs[l];
    ll M = N / q;
    int e = smalls[half(M / q)] - pc;
    if (e < l + 1)
      break;
    ll t = 0;
    for (int k = l + 1; k <= e; k++)
      t += smalls[half(divide(M, roughs[k]))];
    larges[0] += t - (ll)(e - l) * (pc + l - 1);
  }
  return larges[0] + 1;
}

int myprimepi(int n) {
  if (n < MAXGLOBAL_ARRAYSIZE)
    return prime_cnt[n];
  else
    return prime_pi(n);
}

void prime_omega2(int n) {
  int sqn = 1000;
  primes = sieve(n / sqn);
  cout << "sieve time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC
       << endl;

  int s = 0;
  for (int p : primes) {
    int upto = min(p, (int)n / p);
    s += myprimepi(upto);
  }

  int pcnt_nbysqn = prime_pi(n / sqn);
  s -= prime_pi(sqn) * pcnt_nbysqn;
  for (int p : primes) {
    if (p > sqn)
      break;
    s += (prime_pi(n / p));
  }
}

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  auto start_time = clock();
  int n = (int)1e8;
  int s = prime_omega2(n);
  cout << "W2(" << n << ") = " << s << endl;
  cout << "time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << endl;
  return 0;
}

// Projecteuler 187
