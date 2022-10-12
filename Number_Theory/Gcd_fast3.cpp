
namespace fast_gcd {
using u64 = uint64_t;
using u32 = uint32_t;

__attribute__((target("bmi"))) constexpr u64 binary_gcd(u64 a, u64 b) {
  if (a == 0 || b == 0)
    return a + b;
  int n = __builtin_ctzll(a);
  int m = __builtin_ctzll(b);
  a >>= n;
  b >>= m;
  while (a != b) {
    int m = __builtin_ctzll(a - b);
    bool f = a > b;
    u64 c = f ? a : b;
    b = f ? b : a;
    a = (c - b) >> m;
  }
  return a << min(n, m);
}
} // namespace fast_gcd

using fast_gcd::binary_gcd;
