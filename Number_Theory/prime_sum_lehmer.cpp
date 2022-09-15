#include <bits/stdc++.h>
using namespace std;
// Meissel–Lehmer Algorithm
// Fast Calculation of Prime Counting Fucntion, or sum of F(p) where F is a
// multiplicative function
// O(n^(2/3 + eps)) time complexity and O(n^(1/3 + eps)) space complexity for
// all eps > 0. (Correct me if I'm wrong)
// Credit: chemthan
using T = long long;

// struct meissel_lehmer{
const int maxx = 1e2 + 5, maxy = 1e5 + 5, maxn = 1e7 + 5;
vector<int> lpf, prime, cn;
vector<T> sum;
vector<vector<T>> f;
T F(long long x) { return 1; }
T sum_F(long long x) { return x; }

pair<deque<int>, vector<int>> linearSieve(int n) {
  deque<int> p_list;
  vector<int> lpf(n + 1);
  for (int d = 2; d < n + 1; ++d) {
    if (!lpf[d]) {
      lpf[d] = d;
      p_list.emplace_back(d);
    }
    for (const auto &p : p_list) {
      if (p * d > n || p > lpf[d]) {
        break;
      }
      lpf[p * d] = p;
    }
  }
  return make_pair(p_list, lpf);
}

// meissel_lehmer(): sum(maxn), f(maxx, vector<T>(maxy)), cn(maxn){
void init : sum(maxn), f(maxx, vector<T>(maxy)), cn(maxn) {
  pair<deque<int>, vector<int>> tup = linearSieve(maxn - 1);
  // // tie(lpf, prime)
  deque<int> lpf = tup.first;
  vector<int> prime = tup.second;
  for (int i = 2, cnt = 0; i < maxn; ++i) {
    cout << lpf[i] << " ";
    sum[i] = sum[i - 1];
    // if(lpf[i] == i) sum[i] += F(i), ++ cnt;
    // 	cn[i] = cnt;
  }
  // for(int i = 0; i < maxx; ++ i) for(int j = 0; j < maxy; ++ j){
  // 	f[i][j] = i ? f[i - 1][j] - f[i - 1][j / prime[i - 1]] * F(prime[i - 1])
  // : sum_F(j);
  // }
}
T legendre_sum(long long m, int n) {
  if (!n)
    return sum_F(m);
  if (m <= prime[n - 1])
    return F(1);
  if (m < maxy && n < maxx)
    return f[n][m];
  return legendre_sum(m, n - 1) -
         legendre_sum(m / prime[n - 1], n - 1) * F(prime[n - 1]);
}
T pi(long long m) {
  if (m <= maxn)
    return sum[m];
  int x = sqrt(m + 0.9), y = cbrt(m + 0.9), a = cn[y];
  T res = legendre_sum(m, a) - F(1) + sum[y];
  for (int i = a; prime[i] <= x; ++i)
    res -= (pi(m / prime[i]) - pi(prime[i] - 1)) * F(prime[i]);
  return res;
}
// };

int main() {
  init(10);
  // T x=
  // meissel_lehmer ml;
  // meissel_lehmer.
  // pi(100);
  // cout<<x<<endl;
  return 0;
}