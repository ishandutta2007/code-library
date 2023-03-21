#include <bits/stdc++.h>
using namespace std;

// Without any precomputation it can answer 20 queries in 1 sec
// but if n is small (n<10^7) and we are allowed to prcompute all mobius then it can answer 10^4 queries in ~1 sec
// TODO : This is not performing as fast as python as my divisors function is slower than sympy's by a factor of 100

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
} // namespace PollardRho

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
} // namespace FactorHelper


// Sieves till 1e8 in 2 sec

const ll m = (int)(1e8) + 1;
bool flag[m];
ll p[m], u[m], i, j, k, n, tot;
void mobius_sieve() {
  ll i, j, k;
  // memset(flag,0,sizeof(flag));
  tot = 0;
  for (i = 2, u[1] = 1; i < m; i++) {
    if (!flag[i])
      p[++tot] = i, u[i] = -1;
    for (j = 1; j <= tot; j++) {
      k = i * p[j];
      if (k >= m)
        break;
      flag[k] = 1;
      if (i % p[j] == 0) {
        u[k] = 0;
        break;
      }
      u[k] = -u[i];
    }
  }
}


inline ll T(int n){
    return n * (n + 1) / 2;  
}


ll sum_of_coprimes_of_n_till_k(int n, int k){
    auto divs = FactorHelper::get_divisors(n);
    ll ans = 0;
    for (ll d : divs){
        ans += u[d] * d * T(k/d);
    }
    return ans;
}

int main() {
  auto start_time = clock();
  mobius_sieve();
  cout << "Time till prime and mobius combined sieve: "
       << (1.0 * (clock() - start_time) / CLOCKS_PER_SEC) << "s" << endl;
  PollardRho::init();
  cout << sum_of_coprimes_of_n_till_k(2 * 3 * 5 * 7 * 11 * 13 * 17 * 19, 2 * 3 * 5 * 7 * 11 * 13) << endl;
  cout << "Time till prime and mobius combined sieve: "
       << (1.0 * (clock() - start_time) / CLOCKS_PER_SEC) << "s" << endl;
  return 0;
}
