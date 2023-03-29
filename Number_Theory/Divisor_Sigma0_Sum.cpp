// https://oeis.org/A006218
// a[6] = 6/1+6/2+6/3+6/4+6/5+6/6 = 14
// a[5] = 5/1+5/2+5/3+5/4+5/5     = 10
// a[4] = 4/1+4/2+4/3+4/4         =  8
// a[3] = 3/1+3/2+3/3             =  5

// a[6] = d(1)+d(2)+d(3)+d(4)+d(5)+d(6) = 1+2+2+3+2+4 = 14
// a[5] = d(1)+d(2)+d(3)+d(4)+d(5)      = 1+2+2+3+2   = 10
// a[4] = d(1)+d(2)+d(3)+d(4)           = 1+2+2+3     =  8
// a[3] = d(1)+d(2)+d(3)                = 1+2+2       =  5

//Works till ~1e16 in 1 sec

#include <bits/stdc++.h>
using namespace std;

using int64 = unsigned long long;
using pii = std::pair<int64, int64>;

const int N = 2.16e7, M = 1e4, S = 316228;

std::vector<pii> sqf[S];
int64 sf[S], sg[S], nn, m;
int omega[S], mu[N], sigma[N], e[N];
int ps[M];

void prepare(int n) {
  mu[1] = sigma[1] = 1;
  for (int i = 2, m = 0; i <= n; ++i) {
    if (!sigma[i]) {
      sigma[i] = e[i] = 2;
      mu[i] = -1;
      if ((int64)i * i <= n)
        ps[m++] = i;
    }
    for (int j = 0, u = n / i; j < m && ps[j] <= u; ++j) {
      int p = ps[j], v = p * i;
      if (i % p == 0) {
        mu[v] = 0;
        e[v] = e[i] + 1;
        sigma[v] = sigma[i] / e[i] * e[v];
      } else {
        e[v] = 2;
        mu[v] = -mu[i];
        sigma[v] = sigma[i] * 2;
      }
    }
  }
  std::vector<int64> p2(20, 1), p3(20, 1);
  for (int i = 1; i < 20; ++i) {
    p2[i] = p2[i - 1] * 2;
    p3[i] = p3[i - 1] * 3;
  }
  for (int i = 1; i < S; ++i)
    if (sigma[i] == 2) {
      for (int j = i; j < S; j += i)
        omega[j]++;
    }
  for (int i = 1; i < S; ++i)
    if (mu[i]) {
      for (int j = i, k = 1; j < S; j += i, ++k)
        if (mu[j]) {
          int e3 = omega[k], e2 = omega[j] - e3;
          sqf[i].emplace_back(j, p2[e2] * p3[e3]);
        }
    }
  for (int i = 1; i <= n; ++i) {
    e[i] = e[i - 1] + sigma[i] * mu[i] * mu[i];
    sigma[i] += sigma[i - 1];
  }
}

int64 sum_sigma(int64 n) {
  int64 v = sqrt(n), ret = 0;
  for (int i = 1; i <= v; ++i)
    ret += n / i;
  ret = ret * 2 - (int64)v * v;
  return ret;
}

int64 sum_sigma_mu(int64 n) {
  int64 ret = 0, ud = cbrt(n);
  for (int64 d = 1; d <= ud; ++d)
    if (mu[d]) {
      const auto &s = sqf[d];
      int64 nd = n / d, v = sqrt(nd), sum = 0;
      size_t ui =
          std::upper_bound(s.begin(), s.end(), pii(v + 1, 0)) - s.begin();
      for (size_t i = 0; i < ui; ++i) {
        int64 a = s[i].first, b = s[i].second;
        int64 u = nd / (a * a);
        sum += mu[a] * b * (u <= m ? sigma[u] : sg[nn / u]);
      }
      ret += mu[d] * sum;
    }
  return ret;
}

int main() {
  int T;
  scanf("%d", &T);
  std::vector<int64> tests(T);
  for (int i = 0; i < T; ++i) {
    scanf("%lld", &tests[i]);
  }
  m = *std::max_element(tests.begin(), tests.end());
  m = cbrt(m);
  m = m * m;
  m = std::max<int64>(m, 10000);
  prepare(m);
  for (int i = 1; i <= 10; i++)
    cout << sum_sigma(i) << ", ";
  cout << endl;

  return 0;
}
