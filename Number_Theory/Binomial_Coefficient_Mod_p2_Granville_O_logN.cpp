#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef unsigned long long ull;
typedef pair<ll, int> ii;
typedef vector<int> vi;
typedef vector<ll> vll;
typedef vector<ii> vii;

typedef vector<vii> wgraph;
typedef vector<vi> graph;
typedef vector<vll> vvll;

#define rep(i, n) for (int i = 0; i < n; i++)
#define repx(i, a, b) for (int i = a; i < (int)b; i++)
#define invrep(i, a, b) for (int i = b; i-- > (int)a;)

#define umap unordered_map
#define eb emplace_back

#define debugx(x) cerr << #x << ": " << x << endl
#define debugv(v)                                                              \
  cerr << #v << ": ";                                                          \
  for (auto e : v)                                                             \
    cerr << e << " ";                                                          \
  cerr << endl
#define debugm(m)                                                              \
  cerr << #m << ":\n";                                                         \
  for (auto &v : m) {                                                          \
    for (auto e : v)                                                           \
      cerr << e << " ";                                                        \
    cerr << "\n";                                                              \
  }                                                                            \
  cerr << endl
#define debugmp(m)                                                             \
  cerr << #m << ":\n";                                                         \
  for (auto &v : m) {                                                          \
    cerr << v.first << ": ";                                                   \
    for (auto e : v.second)                                                    \
      cerr << e << " ";                                                        \
    cerr << "\n";                                                              \
  }                                                                            \
  cerr << endl;

struct timer {
  clock_t t;
  timer() { t = clock(); }
  ~timer() { printf("runtime %.6f secs\n", getTime()); }
  double getTime() {
    return ((double)clock() - (double)t) / (double)CLOCKS_PER_SEC;
  }
};

template <typename _Ty1, typename _Ty2>
std::ostream &operator<<(std::ostream &_os, const std::pair<_Ty1, _Ty2> &_p) {
  _os << '(' << _p.first << ',' << _p.second << ')';
  return _os;
}

vi sieve(int n) {
  vi primes;

  vector<bool> is_prime(n + 1, true);
  int limit = (int)floor(sqrt(n));
  repx(i, 2, limit + 1) if (is_prime[i]) for (int j = i * i; j <= n; j += i)
      is_prime[j] = false;

  repx(i, 2, n + 1) if (is_prime[i]) primes.eb(i);

  return primes;
}

ll inline mod(ll x, ll m) { return ((x %= m) < 0) ? x + m : x; }
ll inline mul(ll x, ll y, ll m) { return (x * y) % m; }
ll inline add(ll x, ll y, ll m) { return (x + y) % m; }

ll mulmod(ull a, ull b, ull c) {
  ull x = 0, y = a % c;
  while (b) {
    if (b & 1)
      x = (x + y) % c;
    y = (y << 1) % c;
    b >>= 1;
  }
  return x % c;
}

ll fastPow(ll x, ll n, ll MOD) {
  ll ret = 1;
  while (n) {
    if (n & 1)
      ret = mulmod(ret, x, MOD);
    x = mulmod(x, x, MOD);
    n >>= 1;
  }
  return ret;
}

// extended euclidean algorithm
// finds g, x, y such that
//    a * x + b * y = g = GCD(a,b)
ll gcdext(ll a, ll b, ll &x, ll &y) {
  ll r2, x2, y2, r1, x1, y1, r0, x0, y0, q;
  r2 = a, x2 = 1, y2 = 0;
  r1 = b, x1 = 0, y1 = 1;
  while (r1) {
    q = r2 / r1;
    r0 = r2 % r1;
    x0 = x2 - q * x1;
    y0 = y2 - q * y1;
    r2 = r1, x2 = x1, y2 = y1;
    r1 = r0, x1 = x0, y1 = y0;
  }
  ll g = r2;
  x = x2, y = y2;
  if (g < 0)
    g = -g, x = -x, y = -y; // make sure g > 0
  // for debugging (in case you think you might have bugs)
  // assert (g == a * x + b * y);
  // assert (g == __gcd(abs(a),abs(b)));
  return g;
}
pair<ll, ll> CRT(ll r1, ll m1, ll r2, ll m2) {
  ll g, x, y;
  g = gcdext(m1, m2, x, y);
  if ((r1 - r2) % g != 0)
    return {-1, -1}; // no solution
  ll z = m2 / g;
  ll lcm = m1 * z;
  ll sol =
      add(mod(r1, lcm), m1 * mul(mod(x, z), mod((r2 - r1) / g, z), z), lcm);
  // for debugging (in case you think you might have bugs)
  // assert (0 <= sol and sol < lcm);
  // assert (sol % m1 == r1 % m1);
  // assert (sol % m2 == r2 % m2);
  return {sol, lcm}; // solution + lcm(m1,m2)
}

ll num = 142857;
umap<ll, int> MOD;         // MOD[P]=V_p(mod)
umap<ll, vector<ll>> FMOD; // n! mod p if MOD[p]=1 else the product of all i mod
                           // P^MOD[P], where 1<=i<=n and (i,p)=1
umap<ll, vector<ll>> invFMOD; // the inverse of FMOD[n] in the corresponding MOD

void preCompute() {
  // Factor mod->MOD
  vi primes = sieve(num);
  ll m = num;
  for (auto p : primes) {
    if (p * p > m)
      break;
    while (m % p == 0) {
      MOD[p]++;
      if ((m /= p) == 1)
        goto next;
    }
  }
  if (m > 1)
    MOD[m] = 1;
next:
  // Compute FMOD and invFMOD
  for (auto p : MOD) {
    int m = pow(p.first, p.second); // p^V_p(n)
    FMOD[p.first].assign(m, 1);
    invFMOD[p.first].assign(m, 1);
    repx(i, 2, FMOD[p.first].size()) {
      if (i % p.first == 0 and p.second > 1)
        FMOD[p.first][i] = FMOD[p.first][i - 1];
      else
        FMOD[p.first][i] = mul(FMOD[p.first][i - 1], i, FMOD[p.first].size());

      // Compute using Euler's theorem i.e. a^phi(m)=1 mod m with (a.m)=1
      invFMOD[p.first][i] =
          fastPow(FMOD[p.first][i], m / p.first * (p.first - 1) - 1, m);
    }
  }
}

// Compute nCr using Granville's theorem (prime powers)
// Auxiliary functions

// V_p(n!) using Legendre's theorem
int V(ll n, int p) {
  int e = 0;
  while ((n /= p) > 0)
    e += n;
  return e;
}

//
ll f(ll n, ll p) {
  ll m = pow(p, MOD[p]);
  int e = n / m;
  return mul(fastPow(FMOD[p][m - 1], e, m), FMOD[p][n % m], m);
}
ll F(ll n, ll p) {
  ll m = pow(p, MOD[p]);
  ll ans = 1;
  ll temp = 1;
  while (temp <= n) {
    ans = mul(ans, f(n / temp, p), m);
    temp *= p;
  }
  return ans;
}
// Granville theorem
ll granville(ll n, ll r, int p) {
  int e = V(n, p) - V(n - r, p) - V(r, p);
  ll m = pow(p, MOD[p]);
  if (e >= MOD[p])
    return 0;
  ll ans = fastPow(p, e, m);
  ans = mul(ans, F(n, p), m);
  ans = mul(ans, fastPow(F(r, p), pow(p, MOD[p] - 1) * (p - 1) - 1, m), m);
  ans = mul(ans, fastPow(F(n - r, p), pow(p, MOD[p] - 1) * (p - 1) - 1, m), m);
  return ans;
}

// Compute nCr using Lucas theorem (primes)
ll lucas(ll n, ll r, int p) {
  // Trivial cases
  if (r > n or r < 0)
    return 0;
  if (r == 0 or n == r)
    return 1;
  if (r == 1 or r == n - 1)
    return n % p;
  // Base case
  if (n < p and r < p) {
    ll ans = mul(invFMOD[p][r], invFMOD[p][n - r], p); // 1/(r!(n-r)!) mod p
    ans = mul(ans, FMOD[p][n], p);                     // n!/(r!(n-r!)) mod p
    return ans;
  }
  ll ans = lucas(n / p, r / p, p);           // Recursion
  ans = mul(ans, lucas(n % p, r % p, p), p); // False recursion
  return ans;
}

// Given the prime decomposition of mod;
ll nCr(ll n, ll r) {
  // Trivial cases
  if (n < r or r < 0)
    return 0;
  if (r == 0 or r == n)
    return 1;
  if (r == 1 or r == n - 1)
    return (n % num);
  // Non-trivial cases
  ll ans = 0;
  ll mod = 1;
  for (auto p : MOD) {
    ll temp = pow(p.first, p.second);
    if (p.second > 1) {
      ans = CRT(ans, mod, granville(n, r, p.first), temp).first;
    } else {
      ans = CRT(ans, mod, lucas(n, r, p.first), temp).first;
    }
    mod *= temp;
  }
  return ans;
}

int main() {
  ios::sync_with_stdio(0);
  cin.tie(0);
  cout.tie(0);
  timer t;
  preCompute();

  int t;
  cin >> t;
  rep(_, t) {
    ll n, r;
    cin >> n >> r;
    cout << nCr(n, r) << '\n';
  }
}

// https://github.com/N9199/Programacion-Competitiva/blob/c02fec8a5571a4f1ed0650fccefc6cda9c734a00/Taller%20de%20Programaci%C3%B3n/2019-2/Contest%2011/A.cpp
// https://web.archive.org/web/20170202003812/http://www.dms.umontreal.ca/~andrew/PDF/BinCoeff.pdf
