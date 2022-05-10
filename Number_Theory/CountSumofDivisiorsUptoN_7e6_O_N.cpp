#include <bits/stdc++.h>
using namespace std;
#define MAX 7000004
#define ll unsigned long long

uint divcntsum[MAX + 1] = {0};

// This precomputes for all N=1 to MAX
inline void initDivSum() {
  for (uint i = 1; i <= MAX; i++) {
    for (int j = i; j <= MAX; j += i) {
      divcntsum[j] += 1;
    }
  }
  for (int i = 2; i <= MAX; i++) {
    divcntsum[i] = divcntsum[i - 1] + divcntsum[i];
  }
}
int main() {
  initDivSum();
  for (uint i = 1; i <= 13; i++)
    cout << divcntsum[i] << " " << endl;
  return 0;
}
