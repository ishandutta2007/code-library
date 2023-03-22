// https://math.stackexchange.com/questions/4661467/calculate-sum-k-1n-k-cdot-varphik?noredirect=1&lq=1
// O(n^{3/4}) (for sum k * u(k)) + O(n^{1/2}) (for sum k * phi(k) using O(1) access to sum k * u(k))

// # i*phi(i)

#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using i128 = __int128_t;
#define MOD 999999017


std::ostream &operator<<(std::ostream &dest, __int128_t value) {
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[128];
    char *d = std::end(buffer);
    do {
      --d;
      *d = "0123456789"[tmp % 10];
      tmp /= 10;
    } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}


namespace mu_siv{
const ll m = 100000001;
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
}


namespace i_mu{
const ll MAX_N = (ll)(1e11);
const int MAX_NSQ = (int)(4e5);
ll NN;
int N_SQ;


// unordered_map<ll,ll>mp;
static i128 snu[MAX_NSQ];
static i128 snu_large_inv[MAX_NSQ];


i128 T(ll n) {
    return (i128)n * (n + 1) / 2;
}

i128 S(ll n) {
    return (i128)n * (n + 1) * (2 * n + 1) / 6;
}

i128 sum_k_mu(ll x) {
    if (x == 1) return snu[x] = 1;
    if (x <= N_SQ) {
        if (snu[x] != 0) return snu[x];
    } else {
        if (snu_large_inv[NN/x] !=0) return snu_large_inv[NN/x];
        // if (mp.find(x) != mp.end()) return mp[x];
    }
    int x_sr = sqrtl(x);
    i128 ans = 1;

    for (int d = 2; d <= x_sr; d++) {
        ans = (ans - (i128)(d) * sum_k_mu(x/d) );
    }

    int st = x/x_sr - 1;
    for (int m = st; m >= 1; m--) {
        ans = (ans - (i128)(T(x/m) - T(x/(m+1))) * sum_k_mu(m));
    }

    if (x <= N_SQ) {
        return snu[x] = ans;
    } else {
        return snu_large_inv[NN/x] = ans;
        // return mp[x] = ans;
    }
}
}

ll i_phi_i_good(ll n) {
    ll n_sr = sqrtl(n);
    ll s = 0;
    for(ll d = 1; d <= n_sr ; d++) {
        s = (s + mu_siv::u[d] * d * (S(n / d) % MOD) % MOD + MOD) % MOD;
    }
    // cout << "s after 1st half = "<< s << endl;

    for(ll v = 1 ; v <= n / (n_sr + 1) ; v++ ) {
        ll l = n / (v + 1);
        NN = l;
        N_SQ = sqrtl(l);
        memset(snu_large_inv, 0, sizeof(snu_large_inv));
        ll s_l = i_mu::sum_k_mu(l) % MOD;

        ll r = n / v;
        NN = r;
        N_SQ = sqrtl(r);
        memset(snu_large_inv, 0, sizeof(snu_large_inv));
        ll s_r = i_mu::sum_k_mu(r) % MOD;

        ll sum_k_mu_range = (s_r - s_l + MOD) % MOD;
        // cout << "S(" << v << ") =" << S(v) << "*" << "sum_k_mu(" << r << ")=" << s_r << "-" << "sum_k_mu(" << l << ")=" << s_l <<endl;
        s = (s + (S(v) % MOD) * sum_k_mu_range % MOD + MOD) % MOD;
    }
    return s;
}

int32_t main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    ms::mobius_sieve();
    // for (int i = 1 ; i <= 30; i++) {
    //     memset(snu_large_inv,0,sizeof(snu_large_inv));
    //     NN = i;
    //     N_SQ = sqrtl(i);
    //     cout << i << " : " << sum_k_mu(i) << endl;
    // }
    // vector<ll> ns = {(ll)1e4, (ll)1e5, (ll)1e6, (ll)1e7, (ll)1e8, /*(ll)1e9, (ll)1e10,*/ 99999999019};
    vector<ll> ns = { (ll)99999999019};
    for (ll n : ns) {
        auto start_time = clock();
        i128 ret = i_phi_i_good(n);
        cout << "10^" << 1.0*log(n)/log(10) << " : " << ret << endl;
        cout << "t: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << endl;
    }
    return 0;
}

