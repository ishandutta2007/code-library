// Generates all primes from 3.1622766601683*10^13 to 3.1622776601683*10^13
// range gap = 10000000
// pi(31622766601683) = 1052369844709
// pi(31622776601683) = 1052370166553
// total primes in range = 321844

#include <bits/stdc++.h>
using namespace std;
using ll = long long;

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
  return primes;
}

vector<ll> sieve_range(vector<int> primes, ll range_start, ll range_endin) {
  int sqrt_range_start = (int)ceil(sqrt(range_start));
  int sqrt_range_endin = (int)sqrtl(range_endin);

  vector<ll> primes_range;
  vector<bool> isp2(range_endin - range_start + 1, true);

  for (int i = 0; primes[i] <= sqrt_range_endin /* and i < primes.size()*/;
       i++) {
    ll p = primes[i];
    ll min_multiple =
        (p < sqrt_range_start)
            ? ((ll)(range_start / p + ((range_start % p == 0) ? 0 : 1)) * p)
            : p * p;
    ll max_multiple = (ll)(range_endin / p) * p;
    for (ll p2 = min_multiple; p2 <= max_multiple; p2 += p) {
      isp2[p2 - range_start] = false;
    }
  }
  for (int i = 0; i < isp2.size(); i++) {
    if (isp2[i])
      primes_range.push_back((ll)i + range_start);
  }
  return primes_range;
}

int main(int argc, char *argv[]) {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  auto start_time = clock();

  ll range_start = 31622766601683;
  ll range_endin = range_start + 10000000; // 31622776601683

  auto primes = sieve((int)sqrtl(range_endin) + 10000);
  cout << "time(sieve): " << (double)(clock() - start_time) / CLOCKS_PER_SEC
       << "s" << endl;

  start_time = clock();
  vector<ll> primes_range = sieve_range(primes, range_start, range_endin);
  cout << "primes_range.size() = " << primes_range.size() << endl;
  cout << "smallest prime in range = " << primes_range[0] << endl;
  cout << "largest prime in range = " << primes_range[primes_range.size() - 1]
       << endl;
  cout << "time(sieve_range): "
       << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s" << endl;
}
