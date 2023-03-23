#include <bits/stdc++.h>
using namespace std;
// https://oeis.org/A011755

// 1, 3, 9, 17, 37, 49, 91

// F(7) = 7*8*15/6 - (2.F(3) + 3.F(2) + 4.F(1) + 5.F(1) + 6.F(1) + 7.F(1))
//      = 140-(2.9 + 3.3 + 4.1 + 5.1 + 6.1 + 7.1)
//      = 140-(18 + 9 + 4 + 5 + 6 + 7)
//      = 140-49
//      = 91

// F(8) = (2.F(4) + 3.F(2) + 4.F(2) + 5.F(1) + 6.F(1) + 7.F(1) + 8.F(1))
//      = 2.17 + 3.3 + 4.3 + 5.1 + 6.1 + 7.1 + 8.1
//      = 34 + 9 + 12 + 5 + 6 + 7 + 8
//      = 81

// O(n^{3/4})

// Credits Anton Lunyov
// F(n) = \sum i^2 - \sum_2<=k<=n k.F(n/k)
// No Mobius required

typedef long long LL;

const int MOD = 999999017 * 2;
const int N = 20000000;
const int K = 100000;

int a[N];
int phi[N];
int sphi[N];
int f[K];

int T2(LL n) {
    LL a[3] = { n, n + 1, 2 * n + 1 };
    for (int i = 0; i < 3; i++)
        if (a[i] % 2 == 0) {
            a[i] /= 2;
            break;
        }
    for (int i = 0; i < 3; i++)
        if (a[i] % 3 == 0) {
            a[i] /= 3;
            break;
        }
    int res = 1;
    for (int i = 0; i < 3; i++)
        res = a[i] % MOD * res % MOD;
    return res;
}

int T(LL n) {
    if (n % 2)
        return n % MOD * ((n + 1) / 2 % MOD) % MOD;
    return (n / 2 % MOD) * ((n + 1) % MOD) % MOD;
}

int F(LL n, LL k) {
    LL m = n / k;
    if (m < N)
        return sphi[m];
    int& res = f[k];
    if (res >= 0)
        return res;
    int q = sqrt(m + 0.5);
    int dd = m / (q + 1);
    res = 0;
    for (int d = 2; d <= dd; d++)
        res = (res + LL(d) % MOD * F(n, k * d)) % MOD;
    for (int d = 1; d <= q; d++)
        res = (res + LL(sphi[d]) * (T(m / d) - T(m / (d + 1)))) % MOD;
    if (res < 0)
        res += MOD;
    res = T2(m) - res;
    if (res < 0)
        res += MOD;
    return res;
}

void small_totients_sieve() {
    for (int i = 1; i < N; i++)
        phi[i] = i;
    for (int i = 2; i < N; i++)
        if (!a[i])
            for (int j = i; j < N; j += i) {
                a[j] = 1;
                phi[j] = phi[j] / i * (i - 1);
            }
    for (int i = 1; i < N; i++)
        sphi[i] = (sphi[i - 1] + LL(i) * phi[i]) % MOD;
}

int32_t main() {
    small_totients_sieve();
    auto start_time = clock();
    int n = 100;
    int q = sqrt(n + 0.5);
    for (int i = 1; i <= q; i++)
        cout << i << ":" << sphi[i] << endl;

    for (int i = 1; i < q; i++)
        cout << n / i << ":" << F(n, i) << endl;

    cout << "t: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << endl;
    return 0;
}
