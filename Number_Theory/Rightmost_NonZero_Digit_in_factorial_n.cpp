// https://open.kattis.com/problems/loworderzeros
//
// The algorithm that obtains the rightmost non-zero digit of n!
#include <bits/extc++.h>

using namespace std;

typedef long long ll;

constexpr array<ll, 4> A = {1, 3, 4, 2};
constexpr array<ll, 5> B = {1, 6, 2, 8, 4};

auto f = [](ll n) {
  if (n == 1) {
    return (ll)1;
  }
  ll m = 1, c = 0, i;
  while (n > 1) {
    if (2 == (i = n % 5)) {
      m <<= 1;
    } else if (4 == i) {
      m <<= 2;
    }
    n /= 5;
    c += n;
  }
  if (c & 1) {
    return B[(A[c % 4] * 4 * m) % 5];
  }
  return B[(A[c % 4] * m) % 5];
};

int main() {

  ios_base::sync_with_stdio(false);
  cin.tie(nullptr);
  ll n;
  while (true) {
    cin >> n;
    if (n == 0) {
      break;
    }
    cout << f(n) << '\n';
  }

  return 0;
}