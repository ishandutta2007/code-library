// Use of spf[] (smallest prime factor) and primes[] 
// has given a 3x speed boost by only increasing memory from O(N) to O(3N)

// Time complexity O(N.loglog min(primepi(lpf[i]),primepi(N/i))) = O(N.loglog (N/log N))
// Space complexity O(N)

#include <bits/stdc++.h>
using namespace std;
using ll = long long;

const int N = 1e8 + 9;
int32_t lpf[N + 1];
int32_t primes[N + 1];
int32_t phi[N + 1];

void totients_faster() {
  phi[1] = 1;
  int32_t nprimes = 0;
  for (int32_t k = 2; k <= N; ++k) {
    if (lpf[k] == 0) {
      phi[k] = k - 1;
      lpf[k] = k;
      primes[nprimes++] = k;
    }

    for (ll j = 0; primes[j] < lpf[k] && (ll)primes[j] * (ll)k <= (ll)N; ++j) {
      lpf[primes[j] * k] = primes[j];
      phi[primes[j] * k] = (primes[j] - 1) * phi[k];
    }

    if ((ll)lpf[k] * (ll)k <= (ll)N) {
      lpf[lpf[k] * k] = lpf[k];
      phi[lpf[k] * k] = lpf[k] * phi[k];
    }
  }
}


int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  totients_faster();
  //for (int32_t d = 1; d <= 10; ++d) cout << d << ":" << phi[d] << endl;
  return 0;
}
