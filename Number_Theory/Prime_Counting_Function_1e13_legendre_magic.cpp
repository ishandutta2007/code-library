#include <bits/stdc++.h>
using namespace std;
// n = 1e+10
// result: 455052511
// time: 0.008363s

// n = 1e+11
// result: 4118054813
// time: 0.04879s

// n = 1e+12
// result: 37607912018
// time: 0.174681s

// n = 1e+13
// result: 346065536839
// time: 0.8755s

// n = 1e+14
// result: 3204941750802
// time: 4.41496s

//  Time Complexity Calculation:
//  the loop for p is upto n^{1/4} and for each p we see the k is n^{1/2}/2, n^{1/2}/2 - primepi(n^{1/2}), n^{1/2}/4, n^{1/2}5,......times
//  \sqrt{\sqrt{n}}$ + \sqrt{\sqrt{n/(2)}}$ + ...... + \sqrt{\sqrt{n/(n-2)}}$ + \sqrt{\sqrt{n/(n-1)}}$ = 
//  \int\limits_{1}^{n^{1/4}} \sqrt {n/t} dt
//  O(n^{1/2}.(t)^{1/2})
//  O(n^{1/2}.(n^{1/4})^{1/2})
//  O(n^{1/2 + 1/8})
//  O(n^{5/8})
//  Note that number of terms is not ncessarily n^{1/4} but primepi(n^{1/4}), so the component n^{1/8} is actually very relaxed upper bound

//  More accurate estimate might be this https://math.stackexchange.com/questions/4668989/how-many-iterations-does-it-take-to-converge-to-1-or-lesser

// $\int\limits_{u=1}^{N^{1/4}}u(N^{1/2}/u-(N^{1/2}/u)log (N^{1/2}/u))du$

// $N^{1/2}\int\limits_{u=1}^{N^{1/4}}du -N^{1/2}\int\limits_{u=1}^{N^{1/4}} (0.5log N - log u)/u \ du$

// $N^{1/2}\int\limits_{u=1}^{N^{1/4}}du -N^{1/2}(0.5log N)\int\limits_{u=1}^{N^{1/4}} du/u +
// N^{1/2}\int\limits_{u=1}^{N^{1/4}} (log u)/u \ du$


// $N^{3/4} - (1/8)N^{1/2}(log N)^2 +(2 + (1/4)log N)/N^{1/4} $

//  O(n^{3/4}/(log(n))^2)

// https://math.stackexchange.com/questions/4669811/what-will-be-the-sum-for-n1-4-iterations
// https://ideone.com/eKenZC

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

int main() {
  //  freopen("input.txt", "r", stdin);
  //  freopen("output.txt", "w", stdout);

  // ll n;
  // scanf("%lld", &n);
  // n = 1234567891234;
  // printf("%lld\n", prime_pi(n));
  // n = 123456789123;
  // printf("%lld\n", prime_pi(n));
  // n = 12345678912;
  // printf("%lld\n", prime_pi(n));
  // n = 1234567891;
  // printf("%lld\n", prime_pi(n));
  ios_base::sync_with_stdio(false);
  cin.tie(NULL);
  cout.tie(NULL);

  vector<ll> ns = {(ll)1e10, (ll)1e11, (ll)1e12, (ll)1e13};
  for (ll n : ns) {
    auto start_time = clock();
    cout << "n = " << (double)n << endl;
    ll res = prime_pi(n);
    cout << "result: " << res << endl;
    cout << "time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s"
         << endl;
    cout << endl;
  }

  return 0;
}
