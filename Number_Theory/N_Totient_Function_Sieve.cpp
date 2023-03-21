// Use of spf[] (smallest prime factor) and primes[] 
// has given a 3x speed boost by only increasing memory from O(N) to O(3N)

#include <bits/stdc++.h>
using namespace std;
using ll = long long;

const int N = 1e8 + 9;
int32_t lpf[N + 1];
int32_t primes[N + 1];
int64_t n_phi[N + 1];

// n_phi[] -> n*phi(n)

void totients_faster() {
  n_phi[1] = 1;
  int32_t nprimes = 0;
  for (int32_t k = 2; k <= N; ++k) {
    if (lpf[k] == 0) {
      n_phi[k] = (ll)k*(ll)(k - 1);
      lpf[k] = k;
      primes[nprimes++] = k;
    }
    for (ll j = 0;
         primes[j] < lpf[k] && (ll)primes[j] * (ll)k <= (ll)N;
         ++j) {
      lpf[primes[j] * k] = primes[j];
      n_phi[primes[j] * k] = (ll)(primes[j])*(ll)(primes[j] - 1) * n_phi[k];
    }

    if ((ll)lpf[k] * (ll)k <= (ll)N) {
      lpf[lpf[k] * k] = lpf[k];
      n_phi[lpf[k] * k] = (ll)(lpf[k])*(ll)lpf[k] * n_phi[k];
    }
  }
}


int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  totients_faster();
  for (int32_t d = 1; d <= 10; ++d) cout << d << ":" << n_phi[d] << endl;
  return 0;
}
