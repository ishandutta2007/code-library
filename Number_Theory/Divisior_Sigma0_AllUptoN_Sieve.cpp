#include <bits/stdc++.h>
using namespace std;

vector<int> divcount_all_upto(int N) {
  int sqtN = (int)sqrt(N) + 1;

  int isprime[sqtN], d, n, e, p, k, sum;
  vector <int> A(N + 1, 0);

  for (n = 1; n <= N; n++)
    A[n] = 0;
  for (n = 2; n < sqtN; n++)
    isprime[n] = 1;

  for (d = 2; d < sqtN; d++) {
    if (isprime[d]) {
      for (n = d * d; n < sqtN + 1; n += d)
        isprime[n] = 0;
      for (n = d * d; n <= N; n += d)
        A[n] = d;
    }
  }

  A[1] = 1;
  for (n = 2; n <= N; n++) {
    if (A[n] == 0)
      A[n] = 2;
    else {
      p = A[n], k = n / p, e = 2;
      while (k % p == 0)
        k /= p, e++;
      A[n] = A[k] * e;
    }
  }
  return A;
}

int main() {
  double dtime = clock();
  int N = 10000000;
  vector<int> dc = divcount_all_upto(N);
  for(int i = 1 ;i <= 10; i++) cout << i << ":" << dc[i] << endl;
  printf("time=%.3lf sec.\n", (double)(clock() - dtime) / CLOCKS_PER_SEC);
  return 0;
}

