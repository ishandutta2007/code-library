#include <bits/stdc++.h>
using namespace std;
using ll = long long;

// https://oeis.org/A002088
// rev_large_phi_sum[n/num] = sum of all etf till num

const ll MOD = 998244353;
const int NMAX = 1e6 + 9;
int phi[NMAX];
ll phi_sum[NMAX];
ll rev_large_phi_sum[NMAX];

void small_totient_sum(int N) {
    for (int i = 1; i < N; i++)
    phi[i] = i;
    for (int i = 2; i < N; i++) {
    if (phi[i] == i) {
          for (int j = i; j < N; j += i)
            phi[j] -= phi[j] / i;
        }
    }
    for (ll i = 0; i <= N; i++) {
        phi_sum[i] = phi_sum[i - 1] + phi[i];
    }
}

ll C2(ll n) {
    n = n % MOD;
    return n * (n + 1) % MOD * ((1 + MOD) / 2) % MOD;
}

ll etf_sum(ll num) {
    if (num == 1) {
        return 1;
    }
    if (num <= sqrtl(n)) {
        return phi_sum[num]%MOD;
    }
    if (rev_large_phi_sum[n/num] == 0) {
        ll sum = C2(num);
        for (ll i = 2; i <= num;) {
            ll j = num / (num / i) + 1;
            sum -= ((j - i) % MOD) * etf_sum(num / i);
            i = j;
        }
        rev_large_phi_sum[n/num] = (MOD + sum % MOD) % MOD;
    }
    return rev_large_phi_sum[n/num];
}


int main() {
    ll n = 1e11;
    ll n_sq = (ll)(sqrtl(n));
    small_totient_sum(n_sq + 1);
    cout << etf_sum(1234567890);
}
