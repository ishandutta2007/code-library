// https://math.stackexchange.com/questions/4661944/calculate-sum-k-1n-k-cdot-muk/4662128#4662128
// \sum k*phi(k)
// O(n^{3/4}) Time and O(n^{1/2}) Space
// It is to be noted recursive version is of same speed as iterative if not a tad faster as the recusion depth is 1 so doesn't make any difference

#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using i128 = __int128_t;
#define MOD 1000000007

const ll MAX_N = (ll)(1e11);
const int MAX_NSQ = (int)(4e5);
ll NN;
int N_SQ;

i128 T(ll n) {
    return (i128)n * (n + 1) / 2;
}

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


// unordered_map<ll,ll>mp;
static i128 snu[MAX_NSQ];
static i128 snu_large_inv[MAX_NSQ];

i128 f(ll x) {
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
        ans = (ans - (i128)(d) * f(x/d) );
    }

    int st = x/x_sr - 1;
    for (int m = st; m >= 1; m--) {
        ans = (ans - (i128)(T(x/m) - T(x/(m+1))) * f(m));
    }

    if (x <= N_SQ) {
        return snu[x] = ans;
    } else {
        return snu_large_inv[NN/x] = ans;
        // return mp[x] = ans;
    }
}

int32_t main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    for (int i = 1 ; i <= 30; i++) {
        memset(snu_large_inv,0,sizeof(snu_large_inv));
        NN = i;
        N_SQ = sqrtl(i);
        cout << i << " : " << f(i) << endl;
    }
    vector<ll> ns = {(ll)1e4, (ll)1e5, (ll)1e6, (ll)1e7, (ll)1e8, (ll)1e9, (ll)1e10, (ll)1e11};
    // vector<ll> ns = { (ll)1e11};
    for (ll n : ns) {
        auto start_time = clock();
        memset(snu_large_inv,0,sizeof(snu_large_inv));
        NN = n;
        N_SQ = sqrtl(n);
        i128 ret = f(n);
        cout << "10^" << 1.0*log(n)/log(10) << " : " << ret << endl;
        cout << "t: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << endl;
    }
    return 0;
}
