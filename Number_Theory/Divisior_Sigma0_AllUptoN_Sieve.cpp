#include <bits/stdc++.h>

void divcount_all_upto(int N) {
  int sqtN = (int)sqrt(N) + 1;

  int *A, isprime[sqtN], d, n, e, p, k, sum;
  A = (int *)(malloc)((N + 1) * sizeof(int));

  for (n = 1; n <= N; n++)
    A[n] = 0;
  for (n = 2; n < sqtN; n++)
    isprime[n] = 1;

  for (d = 2; d < sqtN; d++) {
    if (isprime[d]) {
      for (n = d * d; n < 3163; n += d)
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
}

int main() {
  double dtime = clock();
  divcount_all_upto(10000000)
      l printf("time=%.3lf sec.\n", (double)(clock() - dtime) / CLOCKS_PER_SEC);
  return 0;
}
