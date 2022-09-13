#include <bits/stdc++.h>
using namespace std;

using ll = long long;
// In 1975

namespace PollardRho {
mt19937 rnd(chrono::steady_clock::now().time_since_epoch().count());
const int P = 1e6 + 9;
ll seq[P];
int primes[P], spf[P];
inline ll add_mod(ll x, ll y, ll m) { return (x += y) < m ? x : x - m; }
inline ll mul_mod(ll x, ll y, ll m) {
  ll res = __int128(x) * y % m;
  return res;
  // ll res = x * y - (ll)((long double)x * y / m + 0.5) * m;
  // return res < 0 ? res + m : res;
}
inline ll pow_mod(ll x, ll n, ll m) {
  ll res = 1 % m;
  for (; n; n >>= 1) {
    if (n & 1)
      res = mul_mod(res, x, m);
    x = mul_mod(x, x, m);
  }
  return res;
}
// O(it * (logn)^3), it = number of rounds performed
inline bool miller_rabin(ll n) {
  if (n <= 2 || (n & 1 ^ 1))
    return (n == 2);
  if (n < P)
    return spf[n] == n;
  ll c, d, s = 0, r = n - 1;
  for (; !(r & 1); r >>= 1, s++) {
  }
  // each iteration is a round
  for (int i = 0; primes[i] < n && primes[i] < 32; i++) {
    c = pow_mod(primes[i], r, n);
    for (int j = 0; j < s; j++) {
      d = mul_mod(c, c, n);
      if (d == 1 && c != 1 && c != (n - 1))
        return false;
      c = d;
    }
    if (c != 1)
      return false;
  }
  return true;
}
void init() {
  int cnt = 0;
  for (int i = 2; i < P; i++) {
    if (!spf[i])
      primes[cnt++] = spf[i] = i;
    for (int j = 0, k; (k = i * primes[j]) < P; j++) {
      spf[k] = primes[j];
      if (spf[i] == spf[k])
        break;
    }
  }
}
// returns O(n^(1/4))
ll pollard_rho(ll n) {
  while (1) {
    ll x = rnd() % n, y = x, c = rnd() % n, u = 1, v, t = 0;
    ll *px = seq, *py = seq;
    while (1) {
      *py++ = y = add_mod(mul_mod(y, y, n), c, n);
      *py++ = y = add_mod(mul_mod(y, y, n), c, n);
      if ((x = *px++) == y)
        break;
      v = u;
      u = mul_mod(u, abs(y - x), n);
      if (!u)
        return __gcd(v, n);
      if (++t == 32) {
        t = 0;
        if ((u = __gcd(u, n)) > 1 && u < n)
          return u;
      }
    }
    if (t && (u = __gcd(u, n)) > 1 && u < n)
      return u;
  }
}
vector<ll> factorize(ll n) {
  if (n == 1)
    return vector<ll>();
  if (miller_rabin(n))
    return vector<ll>{n};
  vector<ll> v, w;
  while (n > 1 && n < P) {
    v.push_back(spf[n]);
    n /= spf[n];
  }
  if (n >= P) {
    ll x = pollard_rho(n);
    v = factorize(x);
    w = factorize(n / x);
    v.insert(v.end(), w.begin(), w.end());
  }
  sort(v.begin(), v.end());
  return v;
}
} // namespace PollardRho ends

namespace FactorHelper {
auto flat_to_unique_flat(vector<ll> v) {
  v.erase(unique(v.begin(), v.end()), v.end());
  return v;
}

auto flat_to_map_format(vector<ll> v) {
  map<ll, int> mp;
  int n = v.size();
  for (int i = 0; i < n; i++) {
    int count = 1;
    while (i < n - 1 && v[i] == v[i + 1]) {
      count++;
      i++;
    }
    // cout << v[i] << ":" << count << ", ";
    mp[v[i]] = count;
  }
  return mp;
}

auto flat_to_two_arrays_format(vector<ll> v) {
  vector<ll> prime_divisors;
  vector<int> multiplicity;
  int n = v.size();
  for (int i = 0; i < n; i++) {
    int count = 1;
    while (i < n - 1 && v[i] == v[i + 1]) {
      count++;
      i++;
    }
    prime_divisors.push_back(v[i]);
    multiplicity.push_back(count);
  }
  return make_pair(prime_divisors, multiplicity);
}

vector<ll> generate_all_divisors_from_two_arrays_format(
    vector<ll> prime_divisors, vector<int> multiplicity, int current_divisor,
    long current_result) {
  vector<ll> all_facs;
  if (current_divisor == prime_divisors.size()) {
    all_facs.push_back(current_result);
    // cout<<current_result<<endl;
  } else {
    for (int i = 0; i <= multiplicity[current_divisor]; ++i) {
      vector<ll> sub_facs = generate_all_divisors_from_two_arrays_format(
          prime_divisors, multiplicity, current_divisor + 1, current_result);
      for (int j = 0; j < sub_facs.size(); j++) {
        all_facs.push_back(sub_facs[j]);
      }
      current_result *= prime_divisors[current_divisor];
    }
  }
  sort(all_facs.begin(), all_facs.end());
  return all_facs;
}

vector<ll> get_divisors(ll n) {
  auto flat_format = PollardRho::factorize(n);
  auto two_arrays = flat_to_two_arrays_format(flat_format);
  auto all_f = generate_all_divisors_from_two_arrays_format(
      two_arrays.first, two_arrays.second, 0, 1);
  return all_f;
}
} // namespace FactorHelper ends

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  PollardRho::init();
  int t;
  cin >> t;
  while (t--) {
    ll n;
    cin >> n;

    auto flat_format = PollardRho::factorize(n);
    cout << "All flat_format factors:" << endl;
    for (int j = 0; j < flat_format.size(); j++)
      cout << flat_format[j] << " ";
    cout << endl;

    auto unique_flat_format = FactorHelper::flat_to_unique_flat(flat_format);
    cout << "All unique_flat_format factors:" << endl;
    for (int j = 0; j < unique_flat_format.size(); j++)
      cout << unique_flat_format[j] << " ";
    cout << endl;

    map<ll, int> map_format = FactorHelper::flat_to_map_format(flat_format);
    cout << "All map_format factors:" << endl;
    for (map<ll, int>::iterator it = map_format.begin(); it != map_format.end();
         it++)
      cout << it->first << " : " << it->second << endl;

    auto two_arrays = FactorHelper::flat_to_two_arrays_format(flat_format);
    cout << "All two_arrays factors:" << endl;
    cout << "First array(Prime Factors):" << endl;
    for (int j = 0; j < two_arrays.first.size(); j++)
      cout << two_arrays.first[j] << "  ";
    cout << endl;
    cout << "Second array(Counts):" << endl;
    for (int j = 0; j < two_arrays.second.size(); j++)
      cout << two_arrays.second[j] << "  ";
    cout << endl;

    auto all_d = FactorHelper::generate_all_divisors_from_two_arrays_format(
        two_arrays.first, two_arrays.second, 0, 1);
    cout << "generate_all_divisors_from_two_arrays_format:" << endl;
    for (int j = 0; j < all_d.size(); j++)
      cout << all_d[j] << " ";
    cout << endl;

    auto all_d2 = FactorHelper::get_divisors(n);
    cout << "get_divisors(directly):" << endl;
    for (int j = 0; j < all_d2.size(); j++)
      cout << all_d2[j] << " ";
  }
  return 0;
}
// https://judge.yosupo.jp/problem/factorize
