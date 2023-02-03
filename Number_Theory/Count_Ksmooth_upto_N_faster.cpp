#include <bits/stdc++.h>
// This is faster than other becuaase it reverse engineers prime counting
// function of legendre
// and constructs own counting function
// Credits: ironman353 and min_25
// n = 1e+10
// time: 0.3s

using namespace std;

typedef long long int lli;
typedef long double ld;

lli count_k_smooth(lli K, lli N) {
  vector<int> s_cnt(K + 1);
  vector<lli> l_cnt(K + 1);
  for (int i = 1; i <= K; i++) {
    s_cnt[i] = i - 1;
    l_cnt[i] = N / i - 1;
  }
  for (lli p = 2; p <= K; p++) {
    if (s_cnt[p] == s_cnt[p - 1]) {
      continue;
    }
    int p_cnt = s_cnt[p - 1];
    lli q = p * p;
    lli end = min(K, N / q);
    for (lli i = 1; i <= end; i++) {
      lli d = i * p;
      if (d <= v) {
        l_cnt[i] -= l_cnt[d] - p_cnt;
      } else {
        l_cnt[i] -= s_cnt[N / d] - p_cnt;
      }
    }
    for (lli i = K; i >= q; i--) {
      s_cnt[i] -= s_cnt[i / p] - p_cnt;
    }
  }
  lli ans = N;
  for (int i = 1; i <= K; i++) {
    ans -= (l_cnt[i] - s_cnt[i - 1]);
  }
  return ans;
}

int main() {
  lli k = lli(1e5);
  lli n = lli(1e10);
  lli ans = count_k_smooth(k, n);
  cout << ans << "\n";
  cout << "Time elapsed : " << (1.0 * clock() / CLOCKS_PER_SEC) << " seconds\n";
  return 0;
}
