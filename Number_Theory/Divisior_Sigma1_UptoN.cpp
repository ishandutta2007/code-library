#include <bits/stdc++.h>
using namespace std;
#define MAX 7000004
#define ll unsigned long long

uint divcnt[MAX + 1] = {0};

inline void initDivSum() {
  for (uint i = 1; i <= MAX; i++) {
    for (int j = i; j <= MAX; j += i) {
      divcnt[j] += 1;
    }
  }
}
int main() {
  initDivSum();
  for (uint i = 1; i <= 13; i++)
    cout << divcnt[i] << " " << endl;
  return 0;
}