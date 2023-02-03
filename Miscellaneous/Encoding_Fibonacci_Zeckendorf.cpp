// fibonacci encoding using Zeckendorf's theorem
#include <bits/stdc++.h>
using namespace std;
#define N 30
int fib[N];

int zeckendorf(int n) {
  fib[0] = 1;
  fib[1] = 2;
  int i;
  for (i = 2; fib[i - 1] <= n; i++)
    fib[i] = fib[i - 1] + fib[i - 2];
  return (i - 2);
}

void reprFibonacci(int n) {
  while (n > 0) {
    int f = fib[zeckendorf(n)];
    cout << f << " ";
    n = n - f;
  }
}

char *encodeFibonacci(int n) {
  int index = zeckendorf(n);
  char *codeword = (char *)malloc(sizeof(char) * (index + 3));
  int i = index;
  while (n) {
    codeword[i] = '1';
    n = n - fib[i];
    i = i - 1;
    while (i >= 0 && fib[i] > n) {
      codeword[i] = '0';
      i = i - 1;
    }
  }
  codeword[index + 1] = '1';
  codeword[index + 2] = '\0';
  return codeword;
}

int main() {
  int n;
  n = 30;
  cout << "Non-neighbouring Fibonacci Representation of " << n << " is :";
  reprFibonacci(n);
  cout << endl;
  n = 143;
  char *repr = encodeFibonacci(n);
  cout << "Fibonacci code for " << n << " is " << repr << endl;
  free(repr);
  return 0;
}