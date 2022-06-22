// https://csacademy.com/contest/archive/task/and-closure/statement/
//######################################################################## short
// clean for -        AND/OR/XOR

#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define sz (1 << 20)
#define maxN (sz + 11)
#define mod 1000000007
int AND = 0, OR = 1, XOR = 2;
ll n, x, cnt = 0;
vector<ll> a(maxN, 0);
ll poww(ll a, ll b) {
  if (b <= 0)
    return 1;
  ll c = poww(a, b / 2);
  c *= c;
  c %= mod;
  return b % 2 ? c * a % mod : c;
}
void FST(bool inv, int flag) {
  for (int step = 1; 2 * step <= sz; step *= 2) {
    for (int i = 0; i < sz; i += 2 * step) {
      for (int j = 0; j < step; j++) {
        ll &u = a[j + i], &v = a[j + step + i];
        if (flag == AND)
          tie(u, v) = inv ? make_pair((v - u + mod) % mod, u)
                          : make_pair(v, (u + v) % mod); // AND
        if (flag == OR)
          tie(u, v) = inv ? make_pair(v, (u - v + mod) % mod)
                          : make_pair((u + v) % mod, v); // OR
        if (flag == XOR)
          tie(u, v) = make_pair((u + v) % mod, (u - v + mod) % mod);
      }
    }
  }
  if (inv && flag == XOR) {
    for (int i = 0; i < sz; i++) {
      a[i] /= sz; // may be mod inv?
    }
  }
}
int main() {
  cin >> n;
  for (int i = 1; i <= n; i++) {
    cin >> x;
    a[x]++;
  }
  a[0]++;
  FST(0, AND);
  for (int i = 0; i < sz; i++)
    a[i] = poww(a[i], n);
  FST(1, AND);
  for (int i = 0; i < sz; i++)
    if (a[i] != 0)
      cnt++;
  cout << cnt;
  return 0;
}
