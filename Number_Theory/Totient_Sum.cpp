#include <bits/stdc++.h>
using namespace std;
using ll = long long;

// https://oeis.org/A002088
// rev_large_phi_sum[n/num] = sum of all etf till num

const ll MOD = 998244353;
const int N_SMALL_MAX = 1e6 + 9;
int phi[N_SMALL_MAX];
ll phi_sum[N_SMALL_MAX];
ll rev_large_phi_sum[N_SMALL_MAX];

ll n_sq;
ll reduced_N_SMALL_MAX;

void small_totient_sum(int N_sq) {
    for (int i = 1; i < N_sq; i++)
    phi[i] = i;
    for (int i = 2; i < N_sq; i++) {
    if (phi[i] == i) {
          for (int j = i; j < N_sq; j += i)
            phi[j] -= phi[j] / i;
        }
    }
    for (ll i = 0; i <= N_sq; i++) {
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
    if (num <= n_sq) {
        return phi_sum[num]%MOD;
    }
    if (rev_large_phi_sum[reduced_N_SMALL_MAX/num] == 0) {
        ll sum = C2(num);
        for (ll i = 2; i <= num;) {
            ll j = num / (num / i) + 1;
            sum -= ((j - i) % MOD) * etf_sum(num / i);
            i = j;
        }
        rev_large_phi_sum[reduced_N_SMALL_MAX/num] = (MOD + sum % MOD) % MOD;
    }
    return rev_large_phi_sum[reduced_N_SMALL_MAX/num];
}


int main() {
    ll n = 1e11;
    reduced_N_SMALL_MAX = n;
    n_sq = (ll)(sqrtl(n));
    small_totient_sum(n_sq + 1);
    cout << etf_sum(1234567890);
}