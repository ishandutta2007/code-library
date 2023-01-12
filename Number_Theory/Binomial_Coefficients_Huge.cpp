#include <bits/stdc++.h>

using namespace std;

struct huge_unsigned : public vector<uint32_t> {
  uint64_t carry;

  huge_unsigned() : carry(0) {}

  huge_unsigned &operator*=(uint32_t i) {
    register uint64_t p, q, x = i;

    for (register size_t k = 0, l = size(); k < l;
         at(k++) = q, carry = (q >> 32)) {
      p = at(k), q = p * x;

      if (carry != 0)
        q += carry;
    }

    if (carry != 0)
      push_back(carry), carry = 0;

    return *this;
  }

  huge_unsigned &operator/=(uint32_t j) {
    register uint64_t p, q, x = j;

    for (register size_t k = size(); k > 0; at(k) = q, carry = p - q * x)
      p = (carry << 32) | at(--k), q = p / x;

    while (size() > 0 && back() == 0)
      pop_back();

    return *this;
  }

  void write() {
    for (register size_t k = size() - 1; k > 0; k--)
      printf("%u ", at(k));

    printf("%u", at(0));
  }
};

struct binomial_coefficient : public huge_unsigned {
  binomial_coefficient(uint32_t n, uint32_t k) {
    if (k > n) {
      push_back(0);
      return;
    }

    register uint32_t p = n - k, q = min(k, p);

    if (q == 0) {
      push_back(1);
      return;
    }

    push_back(n--);

    for (register uint32_t i = 2; i <= q; *this /= i++)
      *this *= n--;
  }
};

int main() // A sample test program
{
  uint32_t n, k;
  scanf("%u %u", &n, &k);

  clock_t start_time = clock();

  binomial_coefficient c(n, k);

  clock_t end_time = clock();

  double elapsed_time = (end_time - start_time) * 1000.0 / CLOCKS_PER_SEC;

  printf("Computed c( %u, %u ) = ", n, k), c.write(), putchar('\n'),
      printf("Using %u huge unsigned digits(s)\n", c.size()),
      printf("Elapsed time = %f msec.\n", elapsed_time);
}
// https://codeforces.com/blog/entry/55311
