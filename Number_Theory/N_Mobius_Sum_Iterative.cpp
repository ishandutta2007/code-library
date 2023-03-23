// https://math.stackexchange.com/questions/4661944/calculate-sum-k-1n-k-cdot-muk/4662128#4662128
// https://oeis.org/A068340
// \sum k*phi(k)
// The iterative version supports till 10^9 to support beyond me might need to change to int128(check recursive one)
// O(n^{3/4}) Time and O(n^{1/2}) Space

#include <bits/stdc++.h>
using namespace std;
using ll = long long;

const ll MAX_N = (ll)(1e11);
const int MAX_NSQ = (int)(1e6);

ll T(ll n) {
    return n * (n + 1) / 2;
}

// unordered_map <ll,ll> mp;
static ll snu[MAX_NSQ];
static ll snu_large_inv[MAX_NSQ];

ll sum_k_mu(ll n) {
    if (n == 1) return 1;
    double n_sr_d = sqrt(n);
    ll n_sr = sqrtl(n);
    vector <ll> smalls, larges;

    ll final_st = n_sr;
    if (final_st * final_st == n) final_st--;
    for (ll m = 1; m <= final_st; m++) {
        smalls.push_back(m);
    }

    for (ll d = n / (final_st + 1); d >= 2; d--) {
        larges.push_back(n / d);
    }
    larges.push_back(n);

    // cout<<"smalls:";for(int i = 0;i < smalls.size(); i++) cout << smalls[i] << " "; cout << endl;
    // cout<<"larges:";for(int i = 0;i < larges.size(); i++) cout << larges[i] << " "; cout << endl;
    // cout << "smalls[smalls.size() - 1]=" << smalls[smalls.size() - 1] << endl;
    // cout << "larges[0]=" << larges[0] << endl;
// auto small_time = clock();

    // assert(smalls[smalls.size() - 1] * smalls[smalls.size() - 1] < n);
    for(ll i = 0; i < smalls.size(); i++) {
        ll x = smalls[i];
        ll x_sr = sqrtl(x);
        snu[x] = 1;

        ll st = x_sr;
        if (st * st == x) st--;
        for (ll m = 1; m <= st; m++) {
            snu[x] -= (T(x / m) - T(x / (m + 1))) * snu[m];
        }

        for (ll d = x / (st + 1); d >= 2; d--) {
            snu[x] -= d * snu[x / d];
        }
        // cout << "donex=" << x << "snu[x]=" << snu[x] << endl;
    }
// cout << "st: " << (double)(clock() - small_time) / CLOCKS_PER_SEC << endl;
// auto large_time = clock();
// ll ltfh=0,ltsh=0;

    // assert(larges[0] * larges[0] >= n);
    for(int i = 0; i < larges.size(); i++) {
        ll x = larges[i];
        ll x_sr = sqrtl(x);
        // mp[x] = 1;
        snu_large_inv[n / x] = 1;
// auto large_timefh = clock();

        ll st = x_sr ;
        if (st * st == x) st--;
        for (ll m = 1; m <= st; m++) {
            // mp[x] -= (T(x / m) - T(x / (m + 1))) * snu[m];
            snu_large_inv[n / x] -= (T(x / m) - T(x / (m + 1))) * snu[m];
        }
// ltfh+=(clock()-large_timefh);
// auto large_timesh = clock();
        for (ll d = x / (st + 1); d >= 2; d--) {
            if ((x / d) < n_sr_d)
                // mp[x] -= d * snu[x / d];
                snu_large_inv[n / x] -= d * snu[x / d];
            else
                // mp[x] -= d * mp[x / d];
                snu_large_inv[n / x] -= d * snu_large_inv[n / (x / d)];
        }
// ltsh+=(clock()-large_timesh);
    }
        // cout << "lt: " << (double)(clock() - large_time) / CLOCKS_PER_SEC << endl;
        // cout << "ltfh: " << (double)(ltfh) / CLOCKS_PER_SEC << endl;
        // cout << "ltsh: " << (double)(ltsh) / CLOCKS_PER_SEC << endl;

    // return mp[n];
    return snu_large_inv[1];
}

int32_t main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    for (int i = 2; i <= 30; i++) {
        ll ret = sum_k_mu(i);
        cout << i << " : " << ret << endl;
    }

    vector<ll> ns = {(ll)1e4, (ll)1e5, (ll)1e6, (ll)1e7, (ll)1e8, (ll)1e9, (ll)1e10, (ll)1e11};
    for (ll n : ns) {
        auto start_time = clock();
        ll ret = sum_k_mu(n);
        cout << "10^" << 1.0*log(n)/log(10) << " : " << ret << endl;
        cout << "t: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << endl;
    }
    return 0;
}
