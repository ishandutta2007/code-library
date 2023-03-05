#include <bits/stdc++.h>
using namespace std;
using ll = long long;
//counts squarefree numbers upto 1e15 in 0.5 sec
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

int main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  auto start_time = clock();
  auto p = sieve(4e7);
  //cout << p.size() << '\n';
  cout << "time till prime sieve: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s"
       << endl;

const int M = 33554432;
const ll M2 = 1125899906842624ll;

  ll S = M2;

  long pp[9];
  for(int j0 = 0;j0<2063689;j0++) {pp[0] = p[j0];                          S -= M2/pp[0]/pp[0];
  for(int j1 = 0;j1<j0;j1++) {pp[1] = pp[0] * p[j1]; if(pp[1] >= M) break; S += M2/pp[1]/pp[1]; 
  for(int j2 = 0;j2<j1;j2++) {pp[2] = pp[1] * p[j2]; if(pp[2] >= M) break; S -= M2/pp[2]/pp[2];
  for(int j3 = 0;j3<j2;j3++) {pp[3] = pp[2] * p[j3]; if(pp[3] >= M) break; S += M2/pp[3]/pp[3];
  for(int j4 = 0;j4<j3;j4++) {pp[4] = pp[3] * p[j4]; if(pp[4] >= M) break; S -= M2/pp[4]/pp[4];
  for(int j5 = 0;j5<j4;j5++) {pp[5] = pp[4] * p[j5]; if(pp[5] >= M) break; S += M2/pp[5]/pp[5];
  for(int j6 = 0;j6<j5;j6++) {pp[6] = pp[5] * p[j6]; if(pp[6] >= M) break; S -= M2/pp[6]/pp[6];
  for(int j7 = 0;j7<j6;j7++) {pp[7] = pp[6] * p[j7]; if(pp[7] >= M) break; S += M2/pp[7]/pp[7];
  }}}}}}}}
  
  
  cout << S << " in " << (double)(clock()-start_time)/CLOCKS_PER_SEC << " seconds" << endl;

  return 0;
}
// https://judge.yosupo.jp/problem/enumerate_primes
