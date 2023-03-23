#include <bits/stdc++.h>
using namespace std;
// Time complexity O(N.loglog N)
// Space complexity O(N)

const int N = 1e8 + 9;
int phi[N];
void totients_sieve() {
  for (int i = 1; i < N; i++)
    phi[i] = i;
  for (int i = 2; i < N; i++) {
    if (phi[i] == i) {
      for (int j = i; j < N; j += i)
        phi[j] -= phi[j] / i;
    }
  }
}

void totients_sieve_faster() {
  for (int i = 1; i < N; i++)
    phi[i] = i;
  for (int j = 2; j < N; j += 2)
    phi[j] -= phi[j] >> 1;
  for (int i = 3; i < N; i += 2) {
    if (phi[i] == i) {
      for (int j = i; j < N; j += i)
        phi[j] -= phi[j] / i;
    }
  }
}

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  totients_sieve_faster();
  return 0;
}