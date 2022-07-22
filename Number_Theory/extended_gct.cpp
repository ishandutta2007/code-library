//credits: razimantv

#include <cassert>
#include <iostream>
#include <tuple>

template <typename T>
std::tuple<T, T, T> extended_gcd(T a, T b) {
  T sa{1}, sb{0}, ta{0}, tb{1};
  while (b != 0) {
    T q = a / b;

    auto euclid = [q](T& a, T& b) {
      T temp = a;
      a = b;
      b = temp - q * b;
    };
    euclid(a, b);
    euclid(sa, sb);
    euclid(ta, tb);
  }
  return {a, sa, ta};
}

template<typename T>
T modular_inverse(T a, T mod) {
  auto [g, sa, sb] = extended_gcd(a, mod);
  assert(g == 1);
  return sa < 0 ? mod + sa : sa;
}

int main() {
  {
    auto a{314159}, b{265358};
    auto [g, sa, sb]{extended_gcd(a, b)};
    assert(sa * a + sb * b == g);
    std::cout << g << " " << sa << " " << sb << "\n";
  }
  {
    auto a{31415926535ll}, b{8979323846ll};
    auto [g, sa, sb]{extended_gcd(a, b)};
    assert(sa * a + sb * b == g);
    std::cout << g << " " << sa << " " << sb << "\n";
  }
  {
    std::cout << modular_inverse(9, 13) << "\n";
  }
  return 0
