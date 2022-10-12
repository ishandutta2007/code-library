#include <bits/stdc++.h>
using namespace std;

#ifndef POPCNT_HPP
#define POPCNT_HPP

#include <stdint.h>

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

#if defined(__GNUC__) || __has_builtin(__builtin_popcountll)

inline uint64_t popcnt64(uint64_t x) { return __builtin_popcountll(x); }

#elif defined(_MSC_VER) && defined(_WIN64)

#include <nmmintrin.h>

inline uint64_t popcnt64(uint64_t x) { return _mm_popcnt_u64(x); }

#elif defined(_MSC_VER) && defined(_WIN32)

#include <nmmintrin.h>

inline uint64_t popcnt64(uint64_t x) {
  return _mm_popcnt_u32((uint32_t)x) + _mm_popcnt_u32((uint32_t)(x >> 32));
}

#else

inline uint64_t popcnt64(uint64_t x) {
  uint64_t m1 = 0x5555555555555555ull;
  uint64_t m2 = 0x3333333333333333ull;
  uint64_t m4 = 0x0F0F0F0F0F0F0F0Full;
  uint64_t h01 = 0x0101010101010101ull;

  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;

  return (x * h01) >> 56;
}

#endif

#endif

// Maximum cache line size of current CPUs
#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 512
#endif

/// Some functions in primesieve use a large number of variables
/// at the same time. If such functions are inlined then
/// performance drops because not all variables fit into registers
/// which causes register spilling. We annotate such functions
/// with NOINLINE in order to avoid these issues.
///
#if __has_attribute(noinline)
#define NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define NOINLINE __declspec(noinline)
#else
#define NOINLINE
#endif

template <typename T> void clear(std::unique_ptr<T> &ptr) {
  ptr.reset(nullptr);
}

enum {
  BIT0 = 0xfe, // 11111110
  BIT1 = 0xfd, // 11111101
  BIT2 = 0xfb, // 11111011
  BIT3 = 0xf7, // 11110111
  BIT4 = 0xef, // 11101111
  BIT5 = 0xdf, // 11011111
  BIT6 = 0xbf, // 10111111
  BIT7 = 0x7f  // 01111111
};

#if __cplusplus >= 202002L && __has_cpp_attribute(unlikely)
#define if_unlikely(x)                                                         \
  if (x)                                                                       \
  [[unlikely]]
#elif defined(__GNUC__) || __has_builtin(__builtin_expect)
#define if_unlikely(x) if (__builtin_expect(!!(x), 0))
#else
#define if_unlikely(x) if (x)
#endif

#if __cplusplus >= 201703L && __has_cpp_attribute(fallthrough)
#define FALLTHROUGH [[fallthrough]]
#elif __has_attribute(fallthrough)
#define FALLTHROUGH __attribute__((fallthrough))
#else
#define FALLTHROUGH
#endif

#if defined(__GNUC__) || __has_builtin(__builtin_unreachable)
#define UNREACHABLE __builtin_unreachable()
#elif defined(_MSC_VER)
#define UNREACHABLE __assume(0)
#else
#define UNREACHABLE
#endif

/// Update the current sieving prime's multipleIndex
/// and wheelIndex after sieving has finished.
///

// #define CHECK_FINISHED_SMALL(p, wheelIndex)                                       \
//   if_unlikely(p >= sieveEnd) {                                                 \
//     multipleIndex = (uint64_t)(p - sieveEnd);                                  \
//     prime.set(multipleIndex, wheelIndex);                                      \
//     break;                                                                     \
//   }

#define CHECK_FINISHED_MEDIUM(wheelIndex)                                      \
  if_unlikely(p >= sieveEnd) {                                                 \
    multipleIndex = (uint64_t)(p - sieveEnd);                                  \
    if (Bucket::isFull(buckets_[wheelIndex]))                                  \
      memoryPool_.addBucket(buckets_[wheelIndex]);                             \
    buckets_[wheelIndex]++->set(sievingPrime, multipleIndex, wheelIndex);      \
    break;                                                                     \
  }

/// unset bits < start
const array<uint8_t, 37> unsetSmaller = {
    0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xfe, 0xfe,
    0xfe, 0xfe, 0xfc, 0xfc, 0xf8, 0xf8, 0xf8, 0xf8, 0xf0, 0xf0,
    0xe0, 0xe0, 0xe0, 0xe0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0,
    0x80, 0x80, 0x00, 0x00, 0x00, 0x00, 0x00};

/// unset bits > stop
const array<uint8_t, 37> unsetLarger = {
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01, 0x01, 0x01,
    0x01, 0x03, 0x03, 0x07, 0x07, 0x07, 0x07, 0x0f, 0x0f, 0x1f,
    0x1f, 0x1f, 0x1f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x7f,
    0x7f, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};

namespace primesum {

const array<int, 10> small_primes_ = {0, 2, 3, 5, 7, 11, 13, 17, 19, 23};
// Small primes >= 7
// const array<uint64_t, 5> primes = { 7, 11, 13, 17, 19 };

// Prime products of primes >= 7
const array<uint64_t, 5> primeProducts = {210, 2310, 30030, 510510, 9699690};

using int128_t = __int128_t;
using uint128_t = __uint128_t;

/// Get the time in seconds
double get_time() {
  auto now = chrono::steady_clock::now();
  auto time = now.time_since_epoch();
  auto micro = chrono::duration_cast<chrono::microseconds>(time);
  return (double)micro.count() / 1e6;
}

/// Portable namespace, includes functions which (unlike the versions
/// form the C++ standard library) work with the int128_t and
/// uint128_t types (2014).
///
namespace prt {

template <typename T> struct numeric_limits {
  static constexpr T max() { return std::numeric_limits<T>::max(); }
};

template <> struct numeric_limits<int128_t> {
  static constexpr int128_t min() { return ((int128_t)1) << 127; }
  static constexpr int128_t max() { return ~min(); }
};

template <> struct numeric_limits<uint128_t> {
  static constexpr uint128_t min() { return 0; }
  static constexpr uint128_t max() { return ~min(); }
};

template <typename T> struct is_integral {
  enum {
    value = std::is_integral<T>::value || std::is_same<T, int128_t>::value ||
            std::is_same<T, uint128_t>::value
  };
};

template <typename T> struct is_signed {
  enum { value = std::is_signed<T>::value || std::is_same<T, int128_t>::value };
};

template <typename T> struct is_unsigned {
  enum {
    value = std::is_unsigned<T>::value || std::is_same<T, uint128_t>::value
  };
};

} // namespace prt

class int256_t {
public:
  int256_t() : low(0), high(0) {}

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t(T x)
      : low(x), high((x < 0) ? -1 : 0) {}

  bool operator==(const int256_t &other) const {
    return low == other.low && high == other.high;
  }

  bool operator!=(const int256_t &other) const { return !(*this == other); }

  bool operator<(const int256_t &other) const {
    int128_t hi1 = high;
    int128_t hi2 = other.high;

    return hi1 < hi2 || (hi1 == hi2 && low < other.low);
  }

  bool operator<=(const int256_t &other) const { return !(other < *this); }

  bool operator>(const int256_t &other) const { return other < *this; }

  bool operator>=(const int256_t &other) const { return !(*this < other); }

  int256_t operator+() const { return *this; }

  int256_t operator-() const {
    int256_t res = ~*this;
    ++res;
    return res;
  }

  int256_t operator~() const { return int256_t(~low, ~high); }

  int256_t &operator--() {
    if (low == 0)
      high--;

    low--;
    return *this;
  }

  int256_t &operator++() {
    low++;
    if (low == 0)
      high++;

    return *this;
  }

  int256_t operator--(int) {
    int256_t res = *this;
    --*this;
    return res;
  }

  int256_t operator++(int) {
    int256_t res = *this;
    ++*this;
    return res;
  }

  int256_t operator+(const int256_t &other) const {
    int256_t res;

    res.high = high + other.high;
    res.low = low + other.low;

    if (res.low < other.low)
      res.high++;

    return res;
  }

  int256_t operator-(const int256_t &other) const {
    int256_t res;

    res.high = high - other.high;
    res.low = low - other.low;

    if (res.low > low)
      res.high--;

    return res;
  }

  int256_t operator*(const int256_t &other) const {
    auto max64 = prt::numeric_limits<std::uint64_t>::max();

    if (low <= max64 && other.low <= max64 &&
        ((high == 0 || ~high == 0) && (other.high == 0 || ~other.high == 0))) {
      return int256_t(low * other.low, (high == other.high) ? 0 : -1);
    } else {
      auto al = low & max64;
      auto ah = low >> 64;
      auto bl = other.low & max64;
      auto bh = other.low >> 64;

      auto x = al * bl;
      auto y = al * bh;
      auto z = ah * bl;
      auto w = ah * bh;

      return int256_t(y << 64, y >> 64) + int256_t(z << 64, z >> 64) +
             int256_t(x, w + low * other.high + high * other.low);
    }
  }

  int256_t operator/(const int256_t &other) const {
    return div256(other).first;
  }

  int256_t operator%(const int256_t &other) const {
    return div256(other).second;
  }

  int256_t operator&(const int256_t &other) const {
    return int256_t(low & other.low, high & other.high);
  }

  int256_t operator|(const int256_t &other) const {
    return int256_t(low | other.low, high | other.high);
  }

  int256_t operator^(const int256_t &other) const {
    return int256_t(low ^ other.low, high ^ other.high);
  }

  int256_t operator<<(std::size_t bits) const {
    if (bits >= 128)
      return int256_t(0, low << (bits - 128));
    else
      return int256_t(low << bits, (high << bits) | (low >> (128 - bits)));
  }

  int256_t operator>>(std::size_t bits) const {
    if (bits >= 128)
      return int256_t(high >> (bits - 128), 0);
    else
      return int256_t((low >> bits) | (high << (128 - bits)), high >> bits);
  }

  int256_t &operator+=(const int256_t &other) {
    *this = *this + other;
    return *this;
  }

  int256_t &operator-=(const int256_t &other) {
    *this = *this - other;
    return *this;
  }

  int256_t &operator*=(const int256_t &other) {
    *this = *this * other;
    return *this;
  }

  int256_t &operator/=(const int256_t &other) {
    *this = *this / other;
    return *this;
  }

  int256_t &operator%=(const int256_t &other) {
    *this = *this % other;
    return *this;
  }

  int256_t &operator&=(const int256_t &other) {
    *this = *this & other;
    return *this;
  }

  int256_t &operator|=(const int256_t &other) {
    *this = *this | other;
    return *this;
  }

  int256_t &operator^=(const int256_t &other) {
    *this = *this ^ other;
    return *this;
  }

  int256_t &operator<<=(std::size_t bits) {
    *this = *this << bits;
    return *this;
  }

  int256_t &operator>>=(std::size_t bits) {
    *this = *this >> bits;
    return *this;
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  bool operator==(T x) const {
    return *this == int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  bool operator!=(T x) const {
    return *this != int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  bool operator<(T x) const {
    return *this < int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  bool operator<=(T x) const {
    return *this <= int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  bool operator>(T x) const {
    return *this > int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  bool operator>=(T x) const {
    return *this >= int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator+(T x) const {
    return *this + int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator-(T x) const {
    return *this - int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator*(T x) const {
    return *this * int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator/(T x) const {
    return *this / int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator%(T x) const {
    return *this % int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator&(T x) const {
    return *this & int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator|(T x) const {
    return *this | int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator^(T x) const {
    return *this ^ int256_t(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator<<(T x) const {
    return *this << static_cast<std::size_t>(x);
  }

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t operator>>(T x) const {
    return *this >> static_cast<std::size_t>(x);
  }

  operator std::int8_t() const {
    return (*this < 0)
               ? -static_cast<std::int8_t>(
                     (low - 1) ^ prt::numeric_limits<uint128_t>::max())
               : static_cast<std::int8_t>(low);
  }

  operator std::int16_t() const {
    return (*this < 0)
               ? -static_cast<std::int16_t>(
                     (low - 1) ^ prt::numeric_limits<uint128_t>::max())
               : static_cast<std::int16_t>(low);
  }

  operator std::int32_t() const {
    return (*this < 0)
               ? -static_cast<std::int32_t>(
                     (low - 1) ^ prt::numeric_limits<uint128_t>::max())
               : static_cast<std::int32_t>(low);
  }

  operator std::int64_t() const {
    return (*this < 0)
               ? -static_cast<std::int64_t>(
                     (low - 1) ^ prt::numeric_limits<uint128_t>::max())
               : static_cast<std::int64_t>(low);
  }

  operator int128_t() const {
    return (*this < 0)
               ? -static_cast<int128_t>((low - 1) ^
                                        prt::numeric_limits<uint128_t>::max())
               : static_cast<int128_t>(low);
  }

  operator std::uint8_t() const { return static_cast<std::uint8_t>(low); }

  operator std::uint16_t() const { return static_cast<std::uint16_t>(low); }

  operator std::uint32_t() const { return static_cast<std::uint32_t>(low); }

  operator std::uint64_t() const { return static_cast<std::uint64_t>(low); }

  operator uint128_t() const { return static_cast<uint128_t>(low); }

  friend std::ostream &operator<<(std::ostream &out, int256_t n);

private:
  uint128_t low;
  uint128_t high;

  int256_t(uint128_t low, uint128_t high) : low(low), high(high) {}

  static int256_t min_value() {
    return int256_t(0, prt::numeric_limits<int128_t>::min());
  }

  static int256_t max_value() {
    return int256_t(prt::numeric_limits<uint128_t>::max(),
                    prt::numeric_limits<int128_t>::max());
  }

  bool get_bit(std::size_t bit) const {
    if (bit >= 128)
      return (high >> (bit - 128)) & 1;
    else
      return (low >> bit) & 1;
  }

  void set_bit(std::size_t bit, bool value) {
    if (bit >= 256)
      return;

    if (bit >= 128) {
      uint128_t mask = static_cast<uint128_t>(1) << (bit - 128);

      if (value)
        high |= mask;
      else
        high &= ~mask;
    } else {
      uint128_t mask = static_cast<uint128_t>(1) << bit;

      if (value)
        low |= mask;
      else
        low &= ~mask;
    }
  }

  int find_most_significant_bit() const {
    auto x = *this;
    int pos = 0;

    while (x != 0) {
      pos++;
      x >>= 1;
    }

    return pos;
  }

  /// Unsigned division with remainder
  std::pair<int256_t, int256_t> udiv256(const int256_t &other) const {
    int256_t zero = 0;
    int256_t one = 1;

    if (other == 0) {
      assert(other != 0);
      std::abort();
      return {zero, zero};
    } else if (other == 1) {
      return {*this, zero};
    } else if (*this == other) {
      return {one, zero};
    } else if (*this == 0 || (*this != min_value() && *this < other)) {
      return {zero, *this};
    } else if (high == 0 && other.high == 0) {
      return {int256_t(low / other.low, 0), int256_t(low % other.low, 0)};
    } else {
      int256_t quotient = 0;
      int256_t remainder = 0;

      for (int i = find_most_significant_bit(); i >= 0 && i <= 256; i--) {
        remainder <<= 1;
        remainder.set_bit(0, get_bit(i));

        if (remainder >= other) {
          remainder -= other;
          quotient.set_bit(i, true);
        }
      }
      return {quotient, remainder};
    }
  }

  /// Signed division with remainder
  std::pair<int256_t, int256_t> div256(const int256_t &other) const {
    if (*this < 0) {
      auto x = -*this;

      if (other < 0) {
        auto res = x.udiv256(-other);
        return {res.first, -res.second};
      } else {
        auto res = x.udiv256(other);
        return {-res.first, -res.second};
      }
    } else {
      if (other < 0) {
        auto res = udiv256(-other);
        return {-res.first, res.second};
      } else
        return udiv256(other);
    }
  }
};

inline std::ostream &operator<<(std::ostream &stream, int256_t n) {
  std::string str;

  if (n < 0) {
    stream << "-";
    n = -n;
  }

  while (n > 0) {
    str += '0' + std::int8_t(n % 10);
    n /= 10;
  }

  if (str.empty())
    str = "0";

  stream << std::string(str.rbegin(), str.rend());

  return stream;
}

template <typename T> struct next_larger_type {
  typedef typename std::conditional<
      std::is_same<T, int64_t>::value, int128_t,
      typename std::conditional<
          std::is_same<T, uint64_t>::value, uint128_t,
          typename std::conditional<
              std::is_same<T, int128_t>::value, int256_t,
              typename std::conditional<std::is_same<T, uint128_t>::value,
                                        int256_t, T>::type>::type>::type>::type
      type;
};

inline std::ostream &operator<<(std::ostream &stream, uint128_t n) {
  std::string str;

  while (n > 0) {
    str += '0' + n % 10;
    n /= 10;
  }
  if (str.empty())
    str = "0";

  stream << std::string(str.rbegin(), str.rend());
  return stream;
}

inline std::ostream &operator<<(std::ostream &stream, int128_t n) {
  if (n < 0) {
    stream << "-";
    n = -n;
  }
  stream << (uint128_t)n;
  return stream;
}

template <typename T> struct fastdiv {
  typedef typename std::conditional<
      sizeof(T) / 2 <= sizeof(uint32_t), uint32_t,
      typename std::conditional<sizeof(T) / 2 <= sizeof(uint64_t), uint64_t,
                                T>::type>::type type;
};

template <typename X, typename Y>
typename std::enable_if<(sizeof(X) == sizeof(Y)), X>::type fast_div(X x, Y y) {
  static_assert(prt::is_integral<X>::value && prt::is_integral<Y>::value,
                "fast_div(x, y): types must be integral");

  using fastdiv_t = typename fastdiv<X>::type;

  if (x <= std::numeric_limits<fastdiv_t>::max() &&
      y <= std::numeric_limits<fastdiv_t>::max()) {
    return (fastdiv_t)x / (fastdiv_t)y;
  }

  return x / y;
}

template <typename X, typename Y>
typename std::enable_if<(sizeof(X) > sizeof(Y)), X>::type fast_div(X x, Y y) {
  static_assert(prt::is_integral<X>::value && prt::is_integral<Y>::value,
                "fast_div(x, y): types must be integral");

  using fastdiv_t = typename fastdiv<X>::type;

  if (x <= std::numeric_limits<fastdiv_t>::max())
    return (fastdiv_t)x / (fastdiv_t)y;

  return x / y;
}

/// Returns 2^64-1 if (x + y) > 2^64-1
inline uint64_t checkedAdd(uint64_t x, uint64_t y) {
  if (x >= std::numeric_limits<uint64_t>::max() - y)
    return std::numeric_limits<uint64_t>::max();
  else
    return x + y;
}

template <typename T> constexpr T numberOfBits(T) {
  return (T)std::numeric_limits<typename std::make_unsigned<T>::type>::digits;
}

template <typename T> inline T ilog2(T x) {
  T log2 = 0;
  T bits = numberOfBits(x);

  for (T i = bits / 2; i > 0; i /= 2) {
    T one = 1;
    if (x >= (one << i)) {
      x >>= i;
      log2 += i;
    }
  }

  return log2;
}

template <typename T> inline T ipow(T x, int n) {
  T r = 1;
  for (int i = 0; i < n; i++)
    r *= x;

  return r;
}

template <typename T> constexpr bool isPow2(T x) {
  return x != 0 && (x & (x - 1)) == 0;
}

/// Integer nth root
template <int N, typename T> inline T iroot(T x) {
  if (N == 0)
    return 0;

  T r = (T)std::pow((double)x, 1.0 / N);

  // fix root too large
  for (; r > 0; r--) {
    if (ipow(r, N - 1) <= x / r)
      break;
  }

  // fix root too small
  while (ipow(r + 1, N - 1) <= x / (r + 1))
    r += 1;

  return r;
}

inline int64_t isquare(int64_t x) { return x * x; }

template <typename A, typename B> inline A ceil_div(A a, B b) {
  assert(b > 0);
  return (A)((a + b - 1) / b);
}

// isqrt code starts
#if __cplusplus >= 201402L

/// C++14 compile time square root using binary search
template <typename T> constexpr T sqrt_helper(T x, T lo, T hi) {
  if (lo == hi)
    return lo;

  const T mid = (lo + hi + 1) / 2;

  if (x / mid < mid)
    return sqrt_helper<T>(x, lo, mid - 1);
  else
    return sqrt_helper(x, mid, hi);
}

template <typename T> constexpr T ct_sqrt(T x) {
  return sqrt_helper<T>(x, 0, x / 2 + 1);
}

#elif __cplusplus >= 201103L

#define MID ((lo + hi + 1) / 2)

/// C++11 compile time square root using binary search
template <typename T> constexpr T sqrt_helper(T x, T lo, T hi) {
  return lo == hi ? lo : ((x / MID < MID) ? sqrt_helper<T>(x, lo, MID - 1)
                                          : sqrt_helper<T>(x, MID, hi));
}

template <typename T> constexpr T ct_sqrt(T x) {
  return sqrt_helper<T>(x, 0, x / 2 + 1);
}

#endif

template <typename T> inline T isqrt(T x) {
  T r = (T)std::sqrt((double)x);

#if __cplusplus >= 201103L
  static const T sqrt_max = ct_sqrt(prt::numeric_limits<T>::max());
  r = std::min(r, sqrt_max);
#endif

  while (r * r > x)
    r--;
  while (x - r * r > r * 2)
    r++;

  return r;
}
// isqrt code ends

class PhiTiny {
public:
  PhiTiny();

  template <typename T> T phi(T x, int64_t a) const {
    assert(a <= max_a());

    T pp = prime_products[a];
    return (x / pp) * totients[a] + phi_[a][x % pp];
  }

  static int64_t get_c(int64_t y) {
    assert(y >= 0);

    if (y >= primes.back())
      return max_a();
    else
      return pi[y];
  }

  static int64_t max_a() { return primes.size() - 1; }

private:
  std::array<std::vector<int16_t>, 7> phi_;
  static const std::array<int, 7> primes;
  static const std::array<int, 7> prime_products;
  static const std::array<int, 7> totients;
  static const std::array<int, 13> pi;
};

extern const PhiTiny phiTiny;

inline bool is_phi_tiny(int64_t a) { return a <= PhiTiny::max_a(); }

template <typename T> T phi_tiny(T x, int64_t a) {
  if (x <= std::numeric_limits<uint32_t>::max())
    return phiTiny.phi((uint32_t)x, a);
  else
    return phiTiny.phi(x, a);
}

class PiTable {
public:
  PiTable(uint64_t max);

  /// Get number of primes <= n
  int64_t operator[](uint64_t n) const {
    assert(n <= max_);
    uint64_t bitmask = 0xffffffffffffffffull >> (63 - n % 64);
    return pi_[n / 64].prime_count + popcnt64(pi_[n / 64].bits & bitmask);
  }

  int64_t size() const { return max_ + 1; }

private:
  struct PiData {
    uint64_t prime_count = 0;
    uint64_t bits = 0;
  };

  std::vector<PiData> pi_;
  uint64_t max_;
};

/// The WheelInit data structure is used to calculate the
/// first multiple >= start of each sieving prime
///
struct WheelInit {
  uint8_t nextMultipleFactor;
  uint8_t wheelIndex;
};

extern const WheelInit wheel30Init[30];
extern const WheelInit wheel210Init[210];

/// The WheelElement data structure is used to skip multiples
/// of small primes using wheel factorization
///
struct WheelElement {
  /// Bitmask used to unset the bit corresponding to the current
  /// multiple of a SievingPrime object
  uint8_t unsetBit;
  /// Factor used to calculate the next multiple of a sieving prime
  /// that is not divisible by any of the wheel factors
  uint8_t nextMultipleFactor;
  /// Overflow needed to correct the next multiple index
  /// (due to sievingPrime = prime / 30)
  uint8_t correct;
  /// Used to calculate the next wheel index:
  /// wheelIndex += next;
  int8_t next;
};

extern const WheelElement wheel30[8 * 8];
extern const WheelElement wheel210[48 * 8];

/// The abstract Wheel class is used skip multiples of small
/// primes in the sieve of Eratosthenes. The EratSmall,
/// EratMedium and EratBig classes are derived from Wheel.
///
template <int MODULO, int SIZE, const WheelElement *WHEEL,
          const WheelInit *INIT>
class Wheel {
public:
  /// Add a new sieving prime to the sieving algorithm.
  /// Calculate the first multiple > segmentLow of prime and
  /// the position within the sieve array of that multiple
  /// and its wheel index. When done store the sieving prime.
  ///
  void addSievingPrime(uint64_t prime, uint64_t segmentLow) {
    assert(segmentLow % 30 == 0);

    // This hack is required because in primesieve the 8
    // bits of each byte (of the sieve array) correspond to
    // the offsets { 7, 11, 13, 17, 19, 23, 29, 31 }.
    // So we are looking for: multiples > segmentLow + 6.
    segmentLow += 6;

    // calculate the first multiple (of prime) > segmentLow
    uint64_t quotient = (segmentLow / prime) + 1;
    quotient = std::max(prime, quotient);
    uint64_t multiple = prime * quotient;
    // prime not needed for sieving
    if (multiple > stop_ || multiple < segmentLow)
      return;

    // calculate the next multiple of prime that is not
    // divisible by any of the wheel's factors
    uint64_t nextMultipleFactor = INIT[quotient % MODULO].nextMultipleFactor;
    uint64_t nextMultiple = prime * nextMultipleFactor;
    if (nextMultiple > stop_ - multiple)
      return;

    nextMultiple += multiple - segmentLow;
    uint64_t multipleIndex = nextMultiple / 30;
    uint64_t wheelIndex =
        wheelOffsets_[prime % 30] + INIT[quotient % MODULO].wheelIndex;
    storeSievingPrime(prime, multipleIndex, wheelIndex);
  }

protected:
  uint64_t stop_ = 0;
  virtual ~Wheel() = default;
  virtual void storeSievingPrime(uint64_t, uint64_t, uint64_t) = 0;

  static uint64_t getMaxFactor() { return WHEEL[0].nextMultipleFactor; }

  /// Cross-off the current multiple of sievingPrime
  /// and calculate its next multiple
  ///
  static void unsetBit(uint8_t *sieve, uint64_t sievingPrime,
                       uint64_t *multipleIndex, uint64_t *wheelIndex) {
    sieve[*multipleIndex] &= WHEEL[*wheelIndex].unsetBit;
    *multipleIndex += WHEEL[*wheelIndex].nextMultipleFactor * sievingPrime;
    *multipleIndex += WHEEL[*wheelIndex].correct;
    *wheelIndex += WHEEL[*wheelIndex].next;
  }

private:
  static const uint64_t wheelOffsets_[30];
};

template <int MODULO, int SIZE, const WheelElement *WHEEL,
          const WheelInit *INIT>
const uint64_t Wheel<MODULO, SIZE, WHEEL, INIT>::wheelOffsets_[30] = {
    0, SIZE * 7, 0, 0,        0, 0, 0, SIZE * 0, 0, 0,
    0, SIZE * 1, 0, SIZE * 2, 0, 0, 0, SIZE * 3, 0, SIZE * 4,
    0, 0,        0, SIZE * 5, 0, 0, 0, 0,        0, SIZE * 6};

/// 3rd wheel, skips multiples of 2, 3 and 5
using Wheel30_t = Wheel<30, 8, wheel30, wheel30Init>;

/// 4th wheel, skips multiples of 2, 3, 5 and 7
using Wheel210_t = Wheel<210, 48, wheel210, wheel210Init>;

/// Each SievingPrime object contains a sieving prime and the
/// position of its next multiple inside the sieve array i.e.
/// multipleIndex and a wheelIndex. In order to reduce the memory
/// usage the multipleIndex and wheelIndex are packed into a
/// single 32-bit variable.
///
class SievingPrime {
public:
  enum {
    MAX_MULTIPLEINDEX = (1 << 23) - 1,
    MAX_WHEELINDEX = (1 << (32 - 23)) - 1
  };

  SievingPrime() = default;

  SievingPrime(uint64_t sievingPrime, uint64_t multipleIndex,
               uint64_t wheelIndex) {
    set(sievingPrime, multipleIndex, wheelIndex);
  }

  void set(uint64_t multipleIndex, uint64_t wheelIndex) {
    assert(multipleIndex <= MAX_MULTIPLEINDEX);
    assert(wheelIndex <= MAX_WHEELINDEX);
    indexes_ = (uint32_t)(multipleIndex | (wheelIndex << 23));
  }

  void set(uint64_t sievingPrime, uint64_t multipleIndex, uint64_t wheelIndex) {
    assert(multipleIndex <= MAX_MULTIPLEINDEX);
    assert(wheelIndex <= MAX_WHEELINDEX);
    indexes_ = (uint32_t)(multipleIndex | (wheelIndex << 23));
    sievingPrime_ = (uint32_t)sievingPrime;
  }

  uint64_t getSievingPrime() const { return sievingPrime_; }

  uint64_t getMultipleIndex() const { return indexes_ & MAX_MULTIPLEINDEX; }

  uint64_t getWheelIndex() const { return indexes_ >> 23; }

private:
  /// multipleIndex = 23 least significant bits of indexes_
  /// wheelIndex = 9 most significant bits of indexes_
  uint32_t indexes_;
  uint32_t sievingPrime_;
}; // SievingPrime ends

/// primeCountApprox(x) >= pi(x)
inline std::size_t primeCountApprox(uint64_t start, uint64_t stop) {
  if (start > stop)
    return 0;
  if (stop <= 10)
    return 4;

  // pi(x) <= x / (log(x) - 1.1) + 5, for x >= 4
  double x = (double)stop;
  double logx = std::log(x);
  double div = logx - 1.1;
  double pix = (stop - start) / div + 5;

  return (std::size_t)pix;
}

inline std::size_t primeCountApprox(uint64_t stop) {
  return primeCountApprox(0, stop);
}

/// EratSmall is an implementation of the segmented sieve
/// of Eratosthenes optimized for small sieving primes that
/// have many multiples per segment.
///
class EratSmall : public Wheel30_t {
public:
  static uint64_t getL1CacheSize(uint64_t);
  void init(uint64_t, uint64_t, uint64_t);
  void crossOff(uint8_t *, uint64_t);
  bool enabled() const { return enabled_; }
  // bool CHECK_FINISHED_SMALL(p, int);
  // bool CHECK_FINISHED_SMALL(uint8_t *, uint8_t *, uint8_t *, uint8_t *,
  // uint64_t);
  bool CHECK_FINISHED_SMALL(uint8_t *, uint8_t *, uint8_t *, uint64_t *,
                            uint64_t);

private:
  uint64_t maxPrime_ = 0;
  uint64_t l1CacheSize_ = 0;
  std::vector<SievingPrime> primes_;
  bool enabled_ = false;
  void storeSievingPrime(uint64_t, uint64_t, uint64_t);
  NOINLINE void crossOff(uint8_t *, uint8_t *);
};

/// @stop:        Upper bound for sieving
/// @l1CacheSize: CPU L1 cache size
/// @maxPrime:    Sieving primes <= maxPrime
///
void EratSmall::init(uint64_t stop, uint64_t l1CacheSize, uint64_t maxPrime) {
  assert(maxPrime <= l1CacheSize * 3);
  assert(l1CacheSize <= SievingPrime::MAX_MULTIPLEINDEX + 1);

  enabled_ = true;
  stop_ = stop;
  maxPrime_ = maxPrime;
  l1CacheSize_ = l1CacheSize;

  size_t count = primeCountApprox(maxPrime);
  primes_.reserve(count);
}

/// Add a new sieving prime to EratSmall
void EratSmall::storeSievingPrime(uint64_t prime, uint64_t multipleIndex,
                                  uint64_t wheelIndex) {
  assert(prime <= maxPrime_);
  uint64_t sievingPrime = prime / 30;
  primes_.emplace_back(sievingPrime, multipleIndex, wheelIndex);
}

/// Both EratMedium and EratBig run fastest using a sieve size
/// that matches the CPU's L2 cache size (or slightly less).
/// However, proportionally EratSmall does a lot more memory
/// writes than both EratMedium and EratBig and hence EratSmall
/// runs fastest using a smaller sieve size that matches the
/// CPU's L1 cache size.
///
/// @sieveSize:   CPU L2 cache size / 2
/// @l1CacheSize: CPU L1 cache size
///
void EratSmall::crossOff(uint8_t *sieve, uint64_t sieveSize) {
  for (uint64_t i = 0; i < sieveSize; i += l1CacheSize_) {
    uint64_t end = i + l1CacheSize_;
    end = std::min(end, sieveSize);
    crossOff(&sieve[i], &sieve[end]);
  }
}

/// Segmented sieve of Eratosthenes with wheel factorization
/// optimized for small sieving primes that have many multiples
/// per segment. This algorithm uses a hardcoded modulo 30
/// wheel that skips multiples of 2, 3 and 5.
///
bool EratSmall::CHECK_FINISHED_SMALL(uint8_t *p, uint8_t *sieveEnd,
                                     uint8_t *prime, uint64_t *multipleIndex,
                                     uint64_t wheelIndex) {
  if_unlikely(p >= sieveEnd) {
    multipleIndex = (uint64_t)(p - sieveEnd);
    prime.set(multipleIndex, wheelIndex);
    return true;
  }
  return false;
}

void EratSmall::crossOff(uint8_t *sieve, uint8_t *sieveEnd) {
  for (auto &prime : primes_) {
    uint64_t sievingPrime = prime.getSievingPrime();
    uint64_t multipleIndex = prime.getMultipleIndex();
    uint64_t wheelIndex = prime.getWheelIndex();
    uint64_t maxLoopDist = sievingPrime * 28 + 27;
    uint8_t *loopEnd = std::max(sieveEnd, sieve + maxLoopDist) - maxLoopDist;

    // Pointer to the byte containing the first multiple of
    // sievingPrime within the current segment.
    uint8_t *p = &sieve[multipleIndex];

    switch (wheelIndex) {
      // sievingPrime % 30 == 7
      for (;;) {
      case 0: // Each iteration removes the next 8
        // multiples of the sievingPrime.
        for (; p < loopEnd; p += sievingPrime * 30 + 7) {
          p[sievingPrime * 0 + 0] &= BIT0;
          p[sievingPrime * 6 + 1] &= BIT4;
          p[sievingPrime * 10 + 2] &= BIT3;
          p[sievingPrime * 12 + 2] &= BIT7;
          p[sievingPrime * 16 + 3] &= BIT6;
          p[sievingPrime * 18 + 4] &= BIT2;
          p[sievingPrime * 22 + 5] &= BIT1;
          p[sievingPrime * 28 + 6] &= BIT5;
        }

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 0) {
            break;
          }
        *p &= BIT0;
        p += sievingPrime * 6 + 1;
        FALLTHROUGH;
      case 1:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 1) {
            break;
          }
        *p &= BIT4;
        p += sievingPrime * 4 + 1;
        FALLTHROUGH;
      case 2:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 2) {
            break;
          }
        *p &= BIT3;
        p += sievingPrime * 2 + 0;
        FALLTHROUGH;
      case 3:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 3) {
            break;
          }
        *p &= BIT7;
        p += sievingPrime * 4 + 1;
        FALLTHROUGH;
      case 4:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 4) {
            break;
          }
        *p &= BIT6;
        p += sievingPrime * 2 + 1;
        FALLTHROUGH;
      case 5:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 5) {
            break;
          }
        *p &= BIT2;
        p += sievingPrime * 4 + 1;
        FALLTHROUGH;
      case 6:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 6) {
            break;
          }
        *p &= BIT1;
        p += sievingPrime * 6 + 1;
        FALLTHROUGH;
      case 7:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 7) {
            break;
          }
        *p &= BIT5;
        p += sievingPrime * 2 + 1;
      }
      break;

      // sievingPrime % 30 == 11
      for (;;) {
      case 8:
        for (; p < loopEnd; p += sievingPrime * 30 + 11) {
          p[sievingPrime * 0 + 0] &= BIT1;
          p[sievingPrime * 6 + 2] &= BIT3;
          p[sievingPrime * 10 + 3] &= BIT7;
          p[sievingPrime * 12 + 4] &= BIT5;
          p[sievingPrime * 16 + 6] &= BIT0;
          p[sievingPrime * 18 + 6] &= BIT6;
          p[sievingPrime * 22 + 8] &= BIT2;
          p[sievingPrime * 28 + 10] &= BIT4;
        }

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 8) {
            break;
          }
        *p &= BIT1;
        p += sievingPrime * 6 + 2;
        FALLTHROUGH;
      case 9:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 9) {
            break;
          }
        *p &= BIT3;
        p += sievingPrime * 4 + 1;
        FALLTHROUGH;
      case 10:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 10) {
            break;
          }
        *p &= BIT7;
        p += sievingPrime * 2 + 1;
        FALLTHROUGH;
      case 11:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 11) {
            break;
          }
        *p &= BIT5;
        p += sievingPrime * 4 + 2;
        FALLTHROUGH;
      case 12:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 12) {
            break;
          }
        *p &= BIT0;
        p += sievingPrime * 2 + 0;
        FALLTHROUGH;
      case 13:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 13) {
            break;
          }
        *p &= BIT6;
        p += sievingPrime * 4 + 2;
        FALLTHROUGH;
      case 14:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 14) {
            break;
          }
        *p &= BIT2;
        p += sievingPrime * 6 + 2;
        FALLTHROUGH;
      case 15:

        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 15) {
            break;
          }
        *p &= BIT4;
        p += sievingPrime * 2 + 1;
      }
      break;

      // sievingPrime % 30 == 13
      for (;;) {
      case 16:
        for (; p < loopEnd; p += sievingPrime * 30 + 13) {
          p[sievingPrime * 0 + 0] &= BIT2;
          p[sievingPrime * 6 + 2] &= BIT7;
          p[sievingPrime * 10 + 4] &= BIT5;
          p[sievingPrime * 12 + 5] &= BIT4;
          p[sievingPrime * 16 + 7] &= BIT1;
          p[sievingPrime * 18 + 8] &= BIT0;
          p[sievingPrime * 22 + 9] &= BIT6;
          p[sievingPrime * 28 + 12] &= BIT3;
        }
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 16) {
            break;
          }
        *p &= BIT2;
        p += sievingPrime * 6 + 2;
        FALLTHROUGH;
      case 17:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 17) {
            break;
          }
        *p &= BIT7;
        p += sievingPrime * 4 + 2;
        FALLTHROUGH;
      case 18:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 18) {
            break;
          }
        *p &= BIT5;
        p += sievingPrime * 2 + 1;
        FALLTHROUGH;
      case 19:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 19) {
            break;
          }
        *p &= BIT4;
        p += sievingPrime * 4 + 2;
        FALLTHROUGH;
      case 20:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 20) {
            break;
          }
        *p &= BIT1;
        p += sievingPrime * 2 + 1;
        FALLTHROUGH;
      case 21:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 21) {
            break;
          }
        *p &= BIT0;
        p += sievingPrime * 4 + 1;
        FALLTHROUGH;
      case 22:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 22) {
            break;
          }
        *p &= BIT6;
        p += sievingPrime * 6 + 3;
        FALLTHROUGH;
      case 23:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 23) {
            break;
          }
        *p &= BIT3;
        p += sievingPrime * 2 + 1;
      }
      break;

      // sievingPrime % 30 == 17
      for (;;) {
      case 24:
        for (; p < loopEnd; p += sievingPrime * 30 + 17) {
          p[sievingPrime * 0 + 0] &= BIT3;
          p[sievingPrime * 6 + 3] &= BIT6;
          p[sievingPrime * 10 + 6] &= BIT0;
          p[sievingPrime * 12 + 7] &= BIT1;
          p[sievingPrime * 16 + 9] &= BIT4;
          p[sievingPrime * 18 + 10] &= BIT5;
          p[sievingPrime * 22 + 12] &= BIT7;
          p[sievingPrime * 28 + 16] &= BIT2;
        }
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 24) {
            break;
          }
        *p &= BIT3;
        p += sievingPrime * 6 + 3;
        FALLTHROUGH;
      case 25:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 25) {
            break;
          }
        *p &= BIT6;
        p += sievingPrime * 4 + 3;
        FALLTHROUGH;
      case 26:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 26) {
            break;
          }
        *p &= BIT0;
        p += sievingPrime * 2 + 1;
        FALLTHROUGH;
      case 27:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 27) {
            break;
          }
        *p &= BIT1;
        p += sievingPrime * 4 + 2;
        FALLTHROUGH;
      case 28:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 28) {
            break;
          }
        *p &= BIT4;
        p += sievingPrime * 2 + 1;
        FALLTHROUGH;
      case 29:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 29) {
            break;
          }
        *p &= BIT5;
        p += sievingPrime * 4 + 2;
        FALLTHROUGH;
      case 30:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 30) {
            break;
          }
        *p &= BIT7;
        p += sievingPrime * 6 + 4;
        FALLTHROUGH;
      case 31:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 31) {
            break;
          }
        *p &= BIT2;
        p += sievingPrime * 2 + 1;
      }
      break;

      // sievingPrime % 30 == 19
      for (;;) {
      case 32:
        for (; p < loopEnd; p += sievingPrime * 30 + 19) {
          p[sievingPrime * 0 + 0] &= BIT4;
          p[sievingPrime * 6 + 4] &= BIT2;
          p[sievingPrime * 10 + 6] &= BIT6;
          p[sievingPrime * 12 + 8] &= BIT0;
          p[sievingPrime * 16 + 10] &= BIT5;
          p[sievingPrime * 18 + 11] &= BIT7;
          p[sievingPrime * 22 + 14] &= BIT3;
          p[sievingPrime * 28 + 18] &= BIT1;
        }
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 32) {
            break;
          }
        *p &= BIT4;
        p += sievingPrime * 6 + 4;
        FALLTHROUGH;
      case 33:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 33) {
            break;
          }
        *p &= BIT2;
        p += sievingPrime * 4 + 2;
        FALLTHROUGH;
      case 34:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 34) {
            break;
          }
        *p &= BIT6;
        p += sievingPrime * 2 + 2;
        FALLTHROUGH;
      case 35:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 35) {
            break;
          }
        *p &= BIT0;
        p += sievingPrime * 4 + 2;
        FALLTHROUGH;
      case 36:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 36) {
            break;
          }
        *p &= BIT5;
        p += sievingPrime * 2 + 1;
        FALLTHROUGH;
      case 37:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 37) {
            break;
          }
        *p &= BIT7;
        p += sievingPrime * 4 + 3;
        FALLTHROUGH;
      case 38:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 38) {
            break;
          }
        *p &= BIT3;
        p += sievingPrime * 6 + 4;
        FALLTHROUGH;
      case 39:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 39) {
            break;
          }
        *p &= BIT1;
        p += sievingPrime * 2 + 1;
      }
      break;

      // sievingPrime % 30 == 23
      for (;;) {
      case 40:
        for (; p < loopEnd; p += sievingPrime * 30 + 23) {
          p[sievingPrime * 0 + 0] &= BIT5;
          p[sievingPrime * 6 + 5] &= BIT1;
          p[sievingPrime * 10 + 8] &= BIT2;
          p[sievingPrime * 12 + 9] &= BIT6;
          p[sievingPrime * 16 + 12] &= BIT7;
          p[sievingPrime * 18 + 14] &= BIT3;
          p[sievingPrime * 22 + 17] &= BIT4;
          p[sievingPrime * 28 + 22] &= BIT0;
        }
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 40) {
            break;
          }
        *p &= BIT5;
        p += sievingPrime * 6 + 5;
        FALLTHROUGH;
      case 41:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 41) {
            break;
          }
        *p &= BIT1;
        p += sievingPrime * 4 + 3;
        FALLTHROUGH;
      case 42:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 42) {
            break;
          }
        *p &= BIT2;
        p += sievingPrime * 2 + 1;
        FALLTHROUGH;
      case 43:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 43) {
            break;
          }
        *p &= BIT6;
        p += sievingPrime * 4 + 3;
        FALLTHROUGH;
      case 44:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 44) {
            break;
          }
        *p &= BIT7;
        p += sievingPrime * 2 + 2;
        FALLTHROUGH;
      case 45:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 45) {
            break;
          }
        *p &= BIT3;
        p += sievingPrime * 4 + 3;
        FALLTHROUGH;
      case 46:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 46) {
            break;
          }
        *p &= BIT4;
        p += sievingPrime * 6 + 5;
        FALLTHROUGH;
      case 47:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 47) {
            break;
          }
        *p &= BIT0;
        p += sievingPrime * 2 + 1;
      }
      break;

      // sievingPrime % 30 == 29
      for (;;) {
      case 48:
        for (; p < loopEnd; p += sievingPrime * 30 + 29) {
          p[sievingPrime * 0 + 0] &= BIT6;
          p[sievingPrime * 6 + 6] &= BIT5;
          p[sievingPrime * 10 + 10] &= BIT4;
          p[sievingPrime * 12 + 12] &= BIT3;
          p[sievingPrime * 16 + 16] &= BIT2;
          p[sievingPrime * 18 + 18] &= BIT1;
          p[sievingPrime * 22 + 22] &= BIT0;
          p[sievingPrime * 28 + 27] &= BIT7;
        }
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 48) {
            break;
          }
        *p &= BIT6;
        p += sievingPrime * 6 + 6;
        FALLTHROUGH;
      case 49:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 49) {
            break;
          }
        *p &= BIT5;
        p += sievingPrime * 4 + 4;
        FALLTHROUGH;
      case 50:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 50) {
            break;
          }
        *p &= BIT4;
        p += sievingPrime * 2 + 2;
        FALLTHROUGH;
      case 51:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 51) {
            break;
          }
        *p &= BIT3;
        p += sievingPrime * 4 + 4;
        FALLTHROUGH;
      case 52:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 52) {
            break;
          }
        *p &= BIT2;
        p += sievingPrime * 2 + 2;
        FALLTHROUGH;
      case 53:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 53) {
            break;
          }
        *p &= BIT1;
        p += sievingPrime * 4 + 4;
        FALLTHROUGH;
      case 54:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 54) {
            break;
          }
        *p &= BIT0;
        p += sievingPrime * 6 + 5;
        FALLTHROUGH;
      case 55:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 55) {
            break;
          }
        *p &= BIT7;
        p += sievingPrime * 2 + 2;
      }
      break;

      // sievingPrime % 30 == 1
      for (;;) {
      case 56:
        for (; p < loopEnd; p += sievingPrime * 30 + 1) {
          p[sievingPrime * 0 + 0] &= BIT7;
          p[sievingPrime * 6 + 1] &= BIT0;
          p[sievingPrime * 10 + 1] &= BIT1;
          p[sievingPrime * 12 + 1] &= BIT2;
          p[sievingPrime * 16 + 1] &= BIT3;
          p[sievingPrime * 18 + 1] &= BIT4;
          p[sievingPrime * 22 + 1] &= BIT5;
          p[sievingPrime * 28 + 1] &= BIT6;
        }
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 56) {
            break;
          }
        *p &= BIT7;
        p += sievingPrime * 6 + 1;
        FALLTHROUGH;
      case 57:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 57) {
            break;
          }
        *p &= BIT0;
        p += sievingPrime * 4 + 0;
        FALLTHROUGH;
      case 58:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 58) {
            break;
          }
        *p &= BIT1;
        p += sievingPrime * 2 + 0;
        FALLTHROUGH;
      case 59:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 59) {
            break;
          }
        *p &= BIT2;
        p += sievingPrime * 4 + 0;
        FALLTHROUGH;
      case 60:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 60) {
            break;
          }
        *p &= BIT3;
        p += sievingPrime * 2 + 0;
        FALLTHROUGH;
      case 61:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 61) {
            break;
          }
        *p &= BIT4;
        p += sievingPrime * 4 + 0;
        FALLTHROUGH;
      case 62:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 62) {
            break;
          }
        *p &= BIT5;
        p += sievingPrime * 6 + 0;
        FALLTHROUGH;
      case 63:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 63) {
            break;
          }
        *p &= BIT6;
        p += sievingPrime * 2 + 0;
      }
      break;

    default:
      UNREACHABLE;
    }
  }
}

class PreSieve {
public:
  void init(uint64_t start, uint64_t stop) {
    // The pre-sieve buffer should be at least 100
    // times smaller than the sieving distance
    // in order to reduce initialization overhead.
    uint64_t dist = stop - start;
    uint64_t threshold = max(dist, isqrt(stop)) / 100;
    auto last = primeProducts.end() - 1;
    auto iter = lower_bound(primeProducts.begin(), last, threshold);
    auto i = distance(primeProducts.begin(), iter);
    const array<uint64_t, 5> primes = {7, 11, 13, 17, 19};

    if (primes.at(i) > maxPrime_)
      initBuffer(primes[i], primeProducts[i]);
  }
  uint64_t getMaxPrime() const { return maxPrime_; }
  void copy(uint8_t *sieve, uint64_t sieveSize, uint64_t segmentLow) const {
    // Find segmentLow index
    uint64_t remainder = segmentLow % primeProduct_;
    uint64_t i = remainder / 30;
    uint64_t sizeLeft = size_ - i;
    auto buffer = buffer_.data();

    if (sieveSize <= sizeLeft)
      copy_n(&buffer[i], sieveSize, sieve);
    else {
      // Copy the last remaining bytes of buffer
      // to the beginning of the sieve array
      copy_n(&buffer[i], sizeLeft, sieve);

      // Restart copying at the beginning of buffer
      for (i = sizeLeft; i + size_ < sieveSize; i += size_)
        copy_n(buffer, size_, &sieve[i]);

      // Copy the last remaining bytes
      copy_n(buffer, sieveSize - i, &sieve[i]);
    }
  };

private:
  uint64_t maxPrime_ = 0;
  uint64_t primeProduct_ = 0;
  uint64_t size_ = 0;
  const array<uint64_t, 5> primes = {7, 11, 13, 17, 19};
  std::vector<uint8_t> buffer_;
  void initBuffer(uint64_t maxPrime, uint64_t primeProduct) {
    maxPrime_ = maxPrime;
    primeProduct_ = primeProduct;
    size_ = primeProduct_ / 30;

    buffer_.clear();
    buffer_.resize(size_, 0xff);

    EratSmall eratSmall;
    uint64_t stop = primeProduct_ * 2;
    eratSmall.init(stop, size_, maxPrime_);

    for (uint64_t prime : primes)
      if (prime <= maxPrime_)
        eratSmall.addSievingPrime(prime, primeProduct_);

    auto buffer = buffer_.data();
    eratSmall.crossOff(buffer, size_);
  }
};

namespace config {

enum {
  /// Number of sieving primes per Bucket in EratSmall, EratMedium
  /// and EratBig objects, affects performance by about 3%.
  /// @pre BUCKET_BYTES must be a power of 2.
  ///
  /// - For x86-64 CPUs after  2010 use 8192
  /// - For x86-64 CPUs before 2010 use 4096
  /// - For PowerPC G4 CPUs    2003 use 2048
  ///
  BUCKET_BYTES = 1 << 13,

  /// The MemoryPool allocates at most MAX_ALLOC_BYTES of new
  /// memory when it runs out of buckets.
  ///
  MAX_ALLOC_BYTES = (1 << 20) * 16,

  /// iterator::prev_prime() caches at least MIN_CACHE_ITERATOR
  /// bytes of primes. Larger is usually faster but also
  /// requires more memory.
  ///
  MIN_CACHE_ITERATOR = (1 << 20) * 8,

  /// iterator::prev_prime() maximum cache size in bytes, used
  /// if pi(sqrt(n)) * 8 bytes > MAX_CACHE_ITERATOR.
  ///
  MAX_CACHE_ITERATOR = (1 << 20) * 1024
};

/// Sieving primes <= (sieveSize in bytes * FACTOR_ERATSMALL)
/// are processed in EratSmall. The ideal value for
/// FACTOR_ERATSMALL has been determined experimentally by
/// running benchmarks near 10^10.
/// @pre FACTOR_ERATSMALL >= 0 && <= 3
///
const double FACTOR_ERATSMALL = 0.175;

/// Sieving primes > (sieveSize in bytes * FACTOR_ERATSMALL)
/// and <= (sieveSize in bytes * FACTOR_ERATMEDIUM)
/// are processed in EratMedium. The ideal value for
/// FACTOR_ERATMEDIUM has been determined experimentally by
/// running benchmarks near 10^14.
///
/// @pre FACTOR_ERATMEDIUM >= 0 && <= 9
/// FACTOR_ERATMEDIUM * max(sieveSize) / 30 * 6 + 6 <= max(multipleIndex)
/// FACTOR_ERATMEDIUM * 2^22 / 30 * 6 + 6 <= 2^23 - 1
/// FACTOR_ERATMEDIUM <= ((2^23 - 7) * 30) / (2^22 * 6)
/// FACTOR_ERATMEDIUM <= 9.999991655
///
const double FACTOR_ERATMEDIUM = 5.0;

/// Each thread sieves at least a distance of MIN_THREAD_DISTANCE
/// in order to reduce the initialization overhead.
/// @pre MIN_THREAD_DISTANCE >= 100
///
const uint64_t MIN_THREAD_DISTANCE = (uint64_t)1e7;

} // namespace config

/// The Bucket data structure is used to store sieving primes.
/// @see http://www.ieeta.pt/~tos/software/prime_sieve.html
/// The Bucket class is designed as a singly linked list, once
/// there is no more space in the current Bucket a new Bucket
/// is allocated.
///
class Bucket {
public:
  SievingPrime *begin() { return &sievingPrimes_[0]; }
  SievingPrime *end() { return end_; }
  Bucket *next() { return next_; }
  void setNext(Bucket *next) { next_ = next; }
  void setEnd(SievingPrime *end) { end_ = end; }
  void reset() { end_ = begin(); }

  /// Get the sieving prime's bucket.
  /// For performance reasons we don't keep an array with all
  /// buckets. Instead we find the sieving prime's bucket by
  /// doing pointer arithmetic using the sieving prime's address.
  /// Since all buckets are aligned by sizeof(Bucket) we
  /// calculate the next address that is smaller than the sieving
  /// prime's address and that is aligned by sizeof(Bucket).
  /// That's the address of the sieving prime's bucket.
  ///
  static Bucket *get(SievingPrime *sievingPrime) {
    assert(sievingPrime != nullptr);
    std::size_t address = (std::size_t)sievingPrime;
    // We need to adjust the address
    // in case the bucket is full.
    address -= 1;
    address -= address % sizeof(Bucket);
    return (Bucket *)address;
  }

  /// Returns true if the bucket is full with sieving primes
  /// (or if there is no bucket i.e. sievingPrime == nullptr).
  /// Each bucket's memory address is aligned by sizeof(Bucket)
  /// (which is a power of 2) in the MemoryPool. This allows
  /// us to quickly check if the bucket is full using the next
  /// sieving prime's address % sizeof(Bucket).
  ///
  static bool isFull(SievingPrime *sievingPrime) {
    std::size_t address = (std::size_t)sievingPrime;
    return address % sizeof(Bucket) == 0;
  }

private:
  enum {
    SIEVING_PRIMES_OFFSET = sizeof(SievingPrime *) + sizeof(Bucket *),
    SIEVING_PRIMES_SIZE =
        (config::BUCKET_BYTES - SIEVING_PRIMES_OFFSET) / sizeof(SievingPrime)
  };

  SievingPrime *end_;
  Bucket *next_;
  SievingPrime sievingPrimes_[SIEVING_PRIMES_SIZE];
};

class MemoryPool {
public:
  NOINLINE void addBucket(SievingPrime *&sievingPrime);
  void freeBucket(Bucket *bucket);

private:
  void allocateBuckets();
  void initBuckets(void *memory, std::size_t bytes);
  void increaseAllocCount();
  /// List of empty buckets
  Bucket *stock_ = nullptr;
  /// Number of buckets to allocate
  std::size_t count_ = 64;
  /// Pointers of allocated buckets
  std::vector<std::unique_ptr<char[]>> memory_;
};

void MemoryPool::addBucket(SievingPrime *&sievingPrime) {
  if (!stock_)
    allocateBuckets();

  Bucket *bucket = stock_;
  stock_ = stock_->next();
  bucket->setNext(nullptr);

  // In case we add a bucket to the front of a
  // non empty bucket list we need to set the
  // next pointer of the new bucket to the bucket
  // that was previously at the front of the list.
  if (sievingPrime) {
    Bucket *old = Bucket::get(sievingPrime);
    old->setEnd(sievingPrime);
    bucket->setNext(old);
  }

  sievingPrime = bucket->begin();
}

void MemoryPool::freeBucket(Bucket *bucket) {
  bucket->reset();
  bucket->setNext(stock_);
  stock_ = bucket;
}

/// primesieve throws a primesieve_error exception
/// if an error occurs e.g. prime > 2^64.
///
class primesieve_error : public std::runtime_error {
public:
  primesieve_error(const std::string &msg) : std::runtime_error(msg) {}
};

void MemoryPool::allocateBuckets() {
  if (memory_.empty())
    memory_.reserve(128);

  // Allocate a large chunk of memory
  size_t bytes = count_ * sizeof(Bucket);
  char *memory = new char[bytes];
  memory_.emplace_back(memory);
  void *ptr = memory;

  // Align pointer address to sizeof(Bucket)
  if (!std::align(sizeof(Bucket), sizeof(Bucket), ptr, bytes))
    throw primesieve_error("MemoryPool: failed to align memory!");

  initBuckets(ptr, bytes);
  increaseAllocCount();
}

void MemoryPool::initBuckets(void *memory, size_t bytes) {
  Bucket *buckets = (Bucket *)memory;
  count_ = bytes / sizeof(Bucket);
  size_t i = 0;

  if ((size_t)buckets % sizeof(Bucket) != 0)
    throw primesieve_error("MemoryPool: failed to align memory!");

  if (count_ < 10)
    throw primesieve_error("MemoryPool: insufficient buckets allocated!");

  for (; i + 1 < count_; i++) {
    buckets[i].reset();
    buckets[i].setNext(&buckets[i + 1]);
  }

  buckets[i].reset();
  buckets[i].setNext(nullptr);
  stock_ = buckets;
}

void MemoryPool::increaseAllocCount() {
  count_ += count_ / 8;
  size_t maxCount = config::MAX_ALLOC_BYTES / sizeof(Bucket);
  count_ = std::min(count_, maxCount);
}

/// EratBig is an implementation of the segmented sieve of
/// Eratosthenes optimized for big sieving primes that have
/// very few multiples per segment.
///
class EratBig : public Wheel210_t {
public:
  void init(uint64_t, uint64_t, uint64_t);
  NOINLINE void crossOff(uint8_t *);
  bool enabled() const { return enabled_; }

private:
  uint64_t maxPrime_ = 0;
  uint64_t log2SieveSize_ = 0;
  uint64_t moduloSieveSize_ = 0;
  std::vector<SievingPrime *> buckets_;
  MemoryPool memoryPool_;
  bool enabled_ = false;
  void storeSievingPrime(uint64_t, uint64_t, uint64_t);
  NOINLINE void crossOff(uint8_t *, Bucket *);
};

/// @stop:      Upper bound for sieving
/// @sieveSize: Sieve size in bytes
/// @maxPrime:  Sieving primes <= maxPrime
///
void EratBig::init(uint64_t stop, uint64_t sieveSize, uint64_t maxPrime) {
  // '>> log2SieveSize' requires power of 2 sieveSize
  assert(isPow2(sieveSize));
  assert(sieveSize <= SievingPrime::MAX_MULTIPLEINDEX + 1);

  enabled_ = true;
  stop_ = stop;
  maxPrime_ = maxPrime;
  log2SieveSize_ = ilog2(sieveSize);
  moduloSieveSize_ = sieveSize - 1;

  uint64_t maxSievingPrime = maxPrime_ / 30;
  uint64_t maxNextMultiple = maxSievingPrime * getMaxFactor() + getMaxFactor();
  uint64_t maxMultipleIndex = sieveSize - 1 + maxNextMultiple;
  uint64_t maxSegmentCount = maxMultipleIndex >> log2SieveSize_;
  uint64_t size = maxSegmentCount + 1;

  buckets_.resize(size);
}

/// Add a new sieving prime
void EratBig::storeSievingPrime(uint64_t prime, uint64_t multipleIndex,
                                uint64_t wheelIndex) {
  assert(prime <= maxPrime_);
  uint64_t sievingPrime = prime / 30;
  uint64_t segment = multipleIndex >> log2SieveSize_;
  multipleIndex &= moduloSieveSize_;

  if (Bucket::isFull(buckets_[segment]))
    memoryPool_.addBucket(buckets_[segment]);

  buckets_[segment]++->set(sievingPrime, multipleIndex, wheelIndex);
}

/// Iterate over the buckets related to the current segment
/// and for each bucket execute crossOff() to remove
/// the multiples of its sieving primes.
///
void EratBig::crossOff(uint8_t *sieve) {
  while (buckets_[0]) {
    Bucket *bucket = Bucket::get(buckets_[0]);
    bucket->setEnd(buckets_[0]);
    buckets_[0] = nullptr;

    while (bucket) {
      crossOff(sieve, bucket);
      Bucket *processed = bucket;
      bucket = bucket->next();
      memoryPool_.freeBucket(processed);
    }
  }

  // Move the bucket related to the next segment to
  // the 1st position so that it will be used when
  // sieving the next segment.
  std::rotate(buckets_.begin(), buckets_.begin() + 1, buckets_.end());
}

/// Removes the next multiple of each sieving prime from the
/// sieve array. After the next multiple of a sieving prime
/// has been removed we calculate its next multiple and
/// determine in which segment that multiple will occur. Then
/// we move the sieving prime to the bucket list related to
/// the previously computed segment.
///
void EratBig::crossOff(uint8_t *sieve, Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  auto buckets = buckets_.data();
  uint64_t moduloSieveSize = moduloSieveSize_;
  uint64_t log2SieveSize = log2SieveSize_;

  for (; prime != end; prime++) {
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint64_t wheelIndex = prime->getWheelIndex();
    uint64_t sievingPrime = prime->getSievingPrime();

    unsetBit(sieve, sievingPrime, &multipleIndex, &wheelIndex);
    uint64_t segment = multipleIndex >> log2SieveSize;
    multipleIndex &= moduloSieveSize;

    if (Bucket::isFull(buckets[segment]))
      memoryPool_.addBucket(buckets[segment]);

    buckets[segment]++->set(sievingPrime, multipleIndex, wheelIndex);
  }
}

/// EratMedium is an implementation of the segmented sieve of
/// Eratosthenes optimized for medium sieving primes
/// that have a few multiples per segment.
///
class EratMedium : public Wheel30_t {
public:
  void init(uint64_t, uint64_t, uint64_t);
  bool enabled() const { return enabled_; }
  NOINLINE void crossOff(uint8_t *, uint64_t);

private:
  bool enabled_ = false;
  uint64_t maxPrime_ = 0;
  MemoryPool memoryPool_;
  std::array<SievingPrime *, 64> buckets_;
  void storeSievingPrime(uint64_t, uint64_t, uint64_t);
  NOINLINE void crossOff_7(uint8_t *, uint8_t *, Bucket *);
  NOINLINE void crossOff_11(uint8_t *, uint8_t *, Bucket *);
  NOINLINE void crossOff_13(uint8_t *, uint8_t *, Bucket *);
  NOINLINE void crossOff_17(uint8_t *, uint8_t *, Bucket *);
  NOINLINE void crossOff_19(uint8_t *, uint8_t *, Bucket *);
  NOINLINE void crossOff_23(uint8_t *, uint8_t *, Bucket *);
  NOINLINE void crossOff_29(uint8_t *, uint8_t *, Bucket *);
  NOINLINE void crossOff_31(uint8_t *, uint8_t *, Bucket *);
};

/// @stop:      Upper bound for sieving
/// @sieveSize: Sieve size in bytes
/// @maxPrime:  Sieving primes <= maxPrime
///
void EratMedium::init(uint64_t stop, uint64_t sieveSize, uint64_t maxPrime) {
  assert(maxPrime <= sieveSize * 9);
  assert(sieveSize * 2 <= SievingPrime::MAX_MULTIPLEINDEX + 1);
  enabled_ = true;
  stop_ = stop;
  maxPrime_ = maxPrime;
  buckets_.fill(nullptr);
}

/// Add a new sieving prime to EratMedium
void EratMedium::storeSievingPrime(uint64_t prime, uint64_t multipleIndex,
                                   uint64_t wheelIndex) {
  assert(prime <= maxPrime_);
  uint64_t sievingPrime = prime / 30;

  if (Bucket::isFull(buckets_[wheelIndex]))
    memoryPool_.addBucket(buckets_[wheelIndex]);

  buckets_[wheelIndex]++->set(sievingPrime, multipleIndex, wheelIndex);
}

void EratMedium::crossOff(uint8_t *sieve, uint64_t sieveSize) {
  // Make a copy of buckets, then reset it
  auto buckets = buckets_;
  buckets_.fill(nullptr);
  uint8_t *sieveEnd = sieve + sieveSize;

  // Iterate over the 64 bucket lists.
  // The 1st list contains sieving primes with wheelIndex = 0.
  // The 2nd list contains sieving primes with wheelIndex = 1.
  // The 3rd list contains sieving primes with wheelIndex = 2.
  // ...
  for (uint64_t i = 0; i < 64; i++) {
    if (!buckets[i])
      continue;

    Bucket *bucket = Bucket::get(buckets[i]);
    bucket->setEnd(buckets[i]);
    uint64_t wheelIndex = i;

    // Iterate over the current bucket list.
    // For each bucket cross off the
    // multiples of its sieving primes.
    while (bucket) {
      switch (wheelIndex / 8) {
      case 0:
        crossOff_7(sieve, sieveEnd, bucket);
        break;
      case 1:
        crossOff_11(sieve, sieveEnd, bucket);
        break;
      case 2:
        crossOff_13(sieve, sieveEnd, bucket);
        break;
      case 3:
        crossOff_17(sieve, sieveEnd, bucket);
        break;
      case 4:
        crossOff_19(sieve, sieveEnd, bucket);
        break;
      case 5:
        crossOff_23(sieve, sieveEnd, bucket);
        break;
      case 6:
        crossOff_29(sieve, sieveEnd, bucket);
        break;
      case 7:
        crossOff_31(sieve, sieveEnd, bucket);
        break;
      default:
        UNREACHABLE;
      }

      Bucket *processed = bucket;
      bucket = bucket->next();
      memoryPool_.freeBucket(processed);
    }
  }
}

/// For sieving primes of type n % 30 == 7
void EratMedium::crossOff_7(uint8_t *sieve, uint8_t *sieveEnd, Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  uint64_t wheelIndex = prime->getWheelIndex();

  for (; prime != end; prime++) {
    uint64_t sievingPrime = prime->getSievingPrime();
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint8_t *p = sieve + multipleIndex;
    uint64_t dist0 = sievingPrime * 6 + 1;
    uint64_t dist1 = sievingPrime * 4 + 1;
    uint64_t dist2 = sievingPrime * 2 + 0;
    uint64_t dist4 = sievingPrime * 2 + 1;

    switch (wheelIndex) {
    default:
      UNREACHABLE;

      for (;;) {
      case 0:
        CHECK_FINISHED_MEDIUM(0);
        *p &= BIT0;
        p += dist0;
        FALLTHROUGH;
      case 1:
        CHECK_FINISHED_MEDIUM(1);
        *p &= BIT4;
        p += dist1;
        FALLTHROUGH;
      case 2:
        CHECK_FINISHED_MEDIUM(2);
        *p &= BIT3;
        p += dist2;
        FALLTHROUGH;
      case 3:
        CHECK_FINISHED_MEDIUM(3);
        *p &= BIT7;
        p += dist1;
        FALLTHROUGH;
      case 4:
        CHECK_FINISHED_MEDIUM(4);
        *p &= BIT6;
        p += dist4;
        FALLTHROUGH;
      case 5:
        CHECK_FINISHED_MEDIUM(5);
        *p &= BIT2;
        p += dist1;
        FALLTHROUGH;
      case 6:
        CHECK_FINISHED_MEDIUM(6);
        *p &= BIT1;
        p += dist0;
        FALLTHROUGH;
      case 7:
        CHECK_FINISHED_MEDIUM(7);
        *p &= BIT5;
        p += dist4;
      }
    }
  }
}

/// For sieving primes of type n % 30 == 11
void EratMedium::crossOff_11(uint8_t *sieve, uint8_t *sieveEnd,
                             Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  uint64_t wheelIndex = prime->getWheelIndex();

  for (; prime != end; prime++) {
    uint64_t sievingPrime = prime->getSievingPrime();
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint8_t *p = sieve + multipleIndex;
    uint64_t dist0 = sievingPrime * 6 + 2;
    uint64_t dist1 = sievingPrime * 4 + 1;
    uint64_t dist2 = sievingPrime * 2 + 1;
    uint64_t dist3 = sievingPrime * 4 + 2;
    uint64_t dist4 = sievingPrime * 2 + 0;

    switch (wheelIndex) {
    default:
      UNREACHABLE;

      for (;;) {
      case 8:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 8) {
            break;
          }
        *p &= BIT1;
        p += dist0;
        FALLTHROUGH;
      case 9:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 9) {
            break;
          }
        *p &= BIT3;
        p += dist1;
        FALLTHROUGH;
      case 10:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 10) {
            break;
          }
        *p &= BIT7;
        p += dist2;
        FALLTHROUGH;
      case 11:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 11) {
            break;
          }
        *p &= BIT5;
        p += dist3;
        FALLTHROUGH;
      case 12:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 12) {
            break;
          }
        *p &= BIT0;
        p += dist4;
        FALLTHROUGH;
      case 13:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 13) {
            break;
          }
        *p &= BIT6;
        p += dist3;
        FALLTHROUGH;
      case 14:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 14) {
            break;
          }
        *p &= BIT2;
        p += dist0;
        FALLTHROUGH;
      case 15:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 15) {
            break;
          }
        *p &= BIT4;
        p += dist2;
      }
    }
  }
}

/// For sieving primes of type n % 30 == 13
void EratMedium::crossOff_13(uint8_t *sieve, uint8_t *sieveEnd,
                             Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  uint64_t wheelIndex = prime->getWheelIndex();

  for (; prime != end; prime++) {
    uint64_t sievingPrime = prime->getSievingPrime();
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint8_t *p = sieve + multipleIndex;
    uint64_t dist0 = sievingPrime * 6 + 2;
    uint64_t dist1 = sievingPrime * 4 + 2;
    uint64_t dist2 = sievingPrime * 2 + 1;
    uint64_t dist5 = sievingPrime * 4 + 1;
    uint64_t dist6 = sievingPrime * 6 + 3;

    switch (wheelIndex) {
    default:
      UNREACHABLE;

      for (;;) {
      case 16:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 16) {
            break;
          }
        *p &= BIT2;
        p += dist0;
        FALLTHROUGH;
      case 17:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 17) {
            break;
          }
        *p &= BIT7;
        p += dist1;
        FALLTHROUGH;
      case 18:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 18) {
            break;
          }
        *p &= BIT5;
        p += dist2;
        FALLTHROUGH;
      case 19:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 19) {
            break;
          }
        *p &= BIT4;
        p += dist1;
        FALLTHROUGH;
      case 20:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 20) {
            break;
          }
        *p &= BIT1;
        p += dist2;
        FALLTHROUGH;
      case 21:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 21) {
            break;
          }
        *p &= BIT0;
        p += dist5;
        FALLTHROUGH;
      case 22:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 22) {
            break;
          }
        *p &= BIT6;
        p += dist6;
        FALLTHROUGH;
      case 23:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 23) {
            break;
          }
        *p &= BIT3;
        p += dist2;
      }
    }
  }
}

/// For sieving primes of type n % 30 == 17
void EratMedium::crossOff_17(uint8_t *sieve, uint8_t *sieveEnd,
                             Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  uint64_t wheelIndex = prime->getWheelIndex();

  for (; prime != end; prime++) {
    uint64_t sievingPrime = prime->getSievingPrime();
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint8_t *p = sieve + multipleIndex;
    uint64_t dist0 = sievingPrime * 6 + 3;
    uint64_t dist1 = sievingPrime * 4 + 3;
    uint64_t dist2 = sievingPrime * 2 + 1;
    uint64_t dist3 = sievingPrime * 4 + 2;
    uint64_t dist6 = sievingPrime * 6 + 4;

    switch (wheelIndex) {
    default:
      UNREACHABLE;

      for (;;) {
      case 24:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 24) {
            break;
          }
        *p &= BIT3;
        p += dist0;
        FALLTHROUGH;
      case 25:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 25) {
            break;
          }
        *p &= BIT6;
        p += dist1;
        FALLTHROUGH;
      case 26:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 26) {
            break;
          }
        *p &= BIT0;
        p += dist2;
        FALLTHROUGH;
      case 27:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 27) {
            break;
          }
        *p &= BIT1;
        p += dist3;
        FALLTHROUGH;
      case 28:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 28) {
            break;
          }
        *p &= BIT4;
        p += dist2;
        FALLTHROUGH;
      case 29:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 29) {
            break;
          }
        *p &= BIT5;
        p += dist3;
        FALLTHROUGH;
      case 30:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 30) {
            break;
          }
        *p &= BIT7;
        p += dist6;
        FALLTHROUGH;
      case 31:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 31) {
            break;
          }
        *p &= BIT2;
        p += dist2;
      }
    }
  }
}

/// For sieving primes of type n % 30 == 19
void EratMedium::crossOff_19(uint8_t *sieve, uint8_t *sieveEnd,
                             Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  uint64_t wheelIndex = prime->getWheelIndex();

  for (; prime != end; prime++) {
    uint64_t sievingPrime = prime->getSievingPrime();
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint8_t *p = sieve + multipleIndex;
    uint64_t dist0 = sievingPrime * 6 + 4;
    uint64_t dist1 = sievingPrime * 4 + 2;
    uint64_t dist2 = sievingPrime * 2 + 2;
    uint64_t dist4 = sievingPrime * 2 + 1;
    uint64_t dist5 = sievingPrime * 4 + 3;

    switch (wheelIndex) {
    default:
      UNREACHABLE;

      for (;;) {
      case 32:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 32) {
            break;
          }
        *p &= BIT4;
        p += dist0;
        FALLTHROUGH;
      case 33:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 33) {
            break;
          }
        *p &= BIT2;
        p += dist1;
        FALLTHROUGH;
      case 34:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 34) {
            break;
          }
        *p &= BIT6;
        p += dist2;
        FALLTHROUGH;
      case 35:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 35) {
            break;
          }
        *p &= BIT0;
        p += dist1;
        FALLTHROUGH;
      case 36:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 36) {
            break;
          }
        *p &= BIT5;
        p += dist4;
        FALLTHROUGH;
      case 37:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 37) {
            break;
          }
        *p &= BIT7;
        p += dist5;
        FALLTHROUGH;
      case 38:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 38) {
            break;
          }
        *p &= BIT3;
        p += dist0;
        FALLTHROUGH;
      case 39:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 39) {
            break;
          }
        *p &= BIT1;
        p += dist4;
      }
    }
  }
}

/// For sieving primes of type n % 30 == 23
void EratMedium::crossOff_23(uint8_t *sieve, uint8_t *sieveEnd,
                             Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  uint64_t wheelIndex = prime->getWheelIndex();

  for (; prime != end; prime++) {
    uint64_t sievingPrime = prime->getSievingPrime();
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint8_t *p = sieve + multipleIndex;
    uint64_t dist0 = sievingPrime * 6 + 5;
    uint64_t dist1 = sievingPrime * 4 + 3;
    uint64_t dist2 = sievingPrime * 2 + 1;
    uint64_t dist4 = sievingPrime * 2 + 2;

    switch (wheelIndex) {
    default:
      UNREACHABLE;

      for (;;) {
      case 40:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 40) {
            break;
          }
        *p &= BIT5;
        p += dist0;
        FALLTHROUGH;
      case 41:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 41) {
            break;
          }
        *p &= BIT1;
        p += dist1;
        FALLTHROUGH;
      case 42:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 42) {
            break;
          }
        *p &= BIT2;
        p += dist2;
        FALLTHROUGH;
      case 43:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 43) {
            break;
          }
        *p &= BIT6;
        p += dist1;
        FALLTHROUGH;
      case 44:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 44) {
            break;
          }
        *p &= BIT7;
        p += dist4;
        FALLTHROUGH;
      case 45:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 45) {
            break;
          }
        *p &= BIT3;
        p += dist1;
        FALLTHROUGH;
      case 46:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 46) {
            break;
          }
        *p &= BIT4;
        p += dist0;
        FALLTHROUGH;
      case 47:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 47) {
            break;
          }
        *p &= BIT0;
        p += dist2;
      }
    }
  }
}

/// For sieving primes of type n % 30 == 29
void EratMedium::crossOff_29(uint8_t *sieve, uint8_t *sieveEnd,
                             Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  uint64_t wheelIndex = prime->getWheelIndex();

  for (; prime != end; prime++) {
    uint64_t sievingPrime = prime->getSievingPrime();
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint8_t *p = sieve + multipleIndex;
    uint64_t dist0 = sievingPrime * 6 + 6;
    uint64_t dist1 = sievingPrime * 4 + 4;
    uint64_t dist2 = sievingPrime * 2 + 2;
    uint64_t dist6 = sievingPrime * 6 + 5;

    switch (wheelIndex) {
    default:
      UNREACHABLE;

      for (;;) {
      case 48:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 48) {
            break;
          }
        *p &= BIT6;
        p += dist0;
        FALLTHROUGH;
      case 49:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 49) {
            break;
          }
        *p &= BIT5;
        p += dist1;
        FALLTHROUGH;
      case 50:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 50) {
            break;
          }
        *p &= BIT4;
        p += dist2;
        FALLTHROUGH;
      case 51:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 51) {
            break;
          }
        *p &= BIT3;
        p += dist1;
        FALLTHROUGH;
      case 52:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 52) {
            break;
          }
        *p &= BIT2;
        p += dist2;
        FALLTHROUGH;
      case 53:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 53) {
            break;
          }
        *p &= BIT1;
        p += dist1;
        FALLTHROUGH;
      case 54:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 54) {
            break;
          }
        *p &= BIT0;
        p += dist6;
        FALLTHROUGH;
      case 55:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 55) {
            break;
          }
        *p &= BIT7;
        p += dist2;
      }
    }
  }
}

/// For sieving primes of type n % 30 == 1
void EratMedium::crossOff_31(uint8_t *sieve, uint8_t *sieveEnd,
                             Bucket *bucket) {
  SievingPrime *prime = bucket->begin();
  SievingPrime *end = bucket->end();
  uint64_t wheelIndex = prime->getWheelIndex();

  for (; prime != end; prime++) {
    uint64_t sievingPrime = prime->getSievingPrime();
    uint64_t multipleIndex = prime->getMultipleIndex();
    uint8_t *p = sieve + multipleIndex;
    uint64_t dist0 = sievingPrime * 6 + 1;
    uint64_t dist1 = sievingPrime * 4 + 0;
    uint64_t dist2 = sievingPrime * 2 + 0;
    uint64_t dist6 = sievingPrime * 6 + 0;

    switch (wheelIndex) {
    default:
      UNREACHABLE;

      for (;;) {
      case 56:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 56) {
            break;
          }
        *p &= BIT7;
        p += dist0;
        FALLTHROUGH;
      case 57:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 57) {
            break;
          }
        *p &= BIT0;
        p += dist1;
        FALLTHROUGH;
      case 58:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 58) {
            break;
          }
        *p &= BIT1;
        p += dist2;
        FALLTHROUGH;
      case 59:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 59) {
            break;
          }
        *p &= BIT2;
        p += dist1;
        FALLTHROUGH;
      case 60:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 60) {
            break;
          }
        *p &= BIT3;
        p += dist2;
        FALLTHROUGH;
      case 61:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 61) {
            break;
          }
        *p &= BIT4;
        p += dist1;
        FALLTHROUGH;
      case 62:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 62) {
            break;
          }
        *p &= BIT5;
        p += dist6;
        FALLTHROUGH;
      case 63:
        if
          CHECK_FINISHED_SMALL(p, sieveEnd, &prime, &multipleIndex, 63) {
            break;
          }
        *p &= BIT6;
        p += dist2;
      }
    }
  }
}

class Erat {
public:
  uint64_t getSieveSize() const;
  uint64_t getStop() const;

protected:
  /// Sieve primes >= start_
  uint64_t start_ = 0;
  /// Sieve primes <= stop_
  uint64_t stop_ = 0;
  /// Size of sieve_ in bytes (power of 2)
  uint64_t sieveSize_ = 0;
  /// Lower bound of the current segment
  uint64_t segmentLow_ = ~0ull;
  /// Upper bound of the current segment
  uint64_t segmentHigh_ = 0;
  /// Sieve of Eratosthenes array
  uint8_t *sieve_ = nullptr;
  Erat();
  Erat(uint64_t, uint64_t);
  void init(uint64_t, uint64_t, uint64_t, PreSieve &);
  void addSievingPrime(uint64_t);
  NOINLINE void sieveSegment();
  bool hasNextSegment() const;
  static uint64_t nextPrime(uint64_t, uint64_t);

private:
  uint64_t maxPreSieve_ = 0;
  uint64_t maxEratSmall_ = 0;
  uint64_t maxEratMedium_ = 0;
  std::unique_ptr<uint8_t[]> deleter_;
  PreSieve *preSieve_ = nullptr;
  EratSmall eratSmall_;
  EratBig eratBig_;
  EratMedium eratMedium_;
  static uint64_t byteRemainder(uint64_t);
  uint64_t getL1CacheSize() const;
  void initSieve(uint64_t);
  void initErat();
  void preSieve();
  void crossOff();
  void sieveLastSegment();
};

/// Convert 1st set bit into prime
inline uint64_t Erat::nextPrime(uint64_t bits, uint64_t low) {
#if defined(ctz64)
  // Find first set 1 bit
  auto bitIndex = ctz64(bits);
  uint64_t bitValue = bitValues[bitIndex];
  uint64_t prime = low + bitValue;
  return prime;
#else
  // Fallback if CTZ instruction is not avilable
  uint64_t debruijn = 0x3F08A4C6ACB9DBDull;
  uint64_t hash = ((bits ^ (bits - 1)) * debruijn) >> 58;
  uint64_t bitValue = bruijnBitValues[hash];
  uint64_t prime = low + bitValue;
  return prime;
#endif
}

inline void Erat::addSievingPrime(uint64_t prime) {
  if (prime > maxEratMedium_)
    eratBig_.addSievingPrime(prime, segmentLow_);
  else if (prime > maxEratSmall_)
    eratMedium_.addSievingPrime(prime, segmentLow_);
  else /* (prime > maxPreSieve) */
    eratSmall_.addSievingPrime(prime, segmentLow_);
}

inline uint64_t Erat::getStop() const { return stop_; }

/// Sieve size in KiB
inline uint64_t Erat::getSieveSize() const { return sieveSize_ >> 10; }

Erat::Erat() = default;

Erat::Erat(uint64_t start, uint64_t stop) : start_(start), stop_(stop) {}

/// @start:     Sieve primes >= start
/// @stop:      Sieve primes <= stop
/// @sieveSize: Sieve size in KiB
/// @preSieve:  Pre-sieve small primes
///
void Erat::init(uint64_t start, uint64_t stop, uint64_t sieveSize,
                PreSieve &preSieve) {
  if (start > stop)
    return;

  assert(start >= 7);
  start_ = start;
  stop_ = stop;
  preSieve_ = &preSieve;
  preSieve_->init(start, stop);
  maxPreSieve_ = preSieve_->getMaxPrime();
  initSieve(sieveSize);

  // The 8 bits of each byte of the sieve array correspond to
  // the offsets { 7, 11, 13, 17, 19, 23, 29, 31 }. If we
  // would set dist = sieveSize * 30 we would not include the
  // last bit of the last byte which corresponds to the offset
  // 31. For this reason we set dist = sieveSize * 30 + 6.
  uint64_t rem = byteRemainder(start);
  uint64_t dist = sieveSize_ * 30 + 6;
  segmentLow_ = start_ - rem;
  segmentHigh_ = checkedAdd(segmentLow_, dist);
  segmentHigh_ = min(segmentHigh_, stop);

  initErat();
}

void Erat::initSieve(uint64_t sieveSize) {
  sieveSize_ = floorPow2(sieveSize);
  sieveSize_ = inBetween(8, sieveSize_, 4096);
  sieveSize_ *= 1024;

  sieve_ = new uint8_t[sieveSize_];
  deleter_.reset(sieve_);
}

void Erat::initErat() {
  uint64_t sqrtStop = isqrt(stop_);
  uint64_t l1CacheSize = getL1CacheSize();

  maxEratSmall_ = (uint64_t)(l1CacheSize * config::FACTOR_ERATSMALL);
  maxEratMedium_ = (uint64_t)(sieveSize_ * config::FACTOR_ERATMEDIUM);

  if (sqrtStop > maxPreSieve_)
    eratSmall_.init(stop_, l1CacheSize, maxEratSmall_);
  if (sqrtStop > maxEratSmall_)
    eratMedium_.init(stop_, sieveSize_, maxEratMedium_);
  if (sqrtStop > maxEratMedium_)
    eratBig_.init(stop_, sieveSize_, sqrtStop);
}

/// EratMedium and EratBig usually run fastest using a sieve
/// size that matches the CPUs L2 cache size. EratSmall
/// however runs fastest using a sieve size that matches the
/// CPUs L1 cache size. Hence we use a smaller sieve size
/// (L1 cache size) in EratSmall and a larger sieve size (L2
/// cache size) in both EratMedium and EratBig.
///
uint64_t Erat::getL1CacheSize() const {
  if (!cpuInfo.hasL1Cache())
    return sieveSize_;

  uint64_t size = cpuInfo.l1CacheSize();
  uint64_t minSize = 8 << 10;
  uint64_t maxSize = 4096 << 10;

  size = std::min(size, sieveSize_);
  size = inBetween(minSize, size, maxSize);

  return size;
}

bool Erat::hasNextSegment() const { return segmentLow_ < stop_; }

uint64_t Erat::byteRemainder(uint64_t n) {
  n %= 30;
  if (n <= 6)
    n += 30;
  return n;
}

/// Pre-sieve multiples of small primes e.g. <= 19
/// to speed up the sieve of Eratosthenes
///
void Erat::preSieve() {
  preSieve_->copy(sieve_, sieveSize_, segmentLow_);

  // unset bits < start
  if (segmentLow_ <= start_) {
    if (start_ <= maxPreSieve_)
      sieve_[0] = 0xff;
    uint64_t rem = byteRemainder(start_);
    sieve_[0] &= unsetSmaller[rem];
  }
}

void Erat::crossOff() {
  if (eratSmall_.enabled())
    eratSmall_.crossOff(sieve_, sieveSize_);
  if (eratMedium_.enabled())
    eratMedium_.crossOff(sieve_, sieveSize_);
  if (eratBig_.enabled())
    eratBig_.crossOff(sieve_);
}

void Erat::sieveSegment() {
  if (segmentHigh_ == stop_)
    sieveLastSegment();
  else {
    preSieve();
    crossOff();

    uint64_t dist = sieveSize_ * 30;
    segmentLow_ = checkedAdd(segmentLow_, dist);
    segmentHigh_ = checkedAdd(segmentHigh_, dist);
    segmentHigh_ = min(segmentHigh_, stop_);
  }
}

void Erat::sieveLastSegment() {
  uint64_t rem = byteRemainder(stop_);
  uint64_t dist = (stop_ - rem) - segmentLow_;
  sieveSize_ = dist / 30 + 1;

  preSieve();
  crossOff();

  // unset bits > stop
  sieve_[sieveSize_ - 1] &= unsetLarger[rem];

  // unset bytes > stop
  uint64_t bytes = sieveSize_ % 8;
  bytes = (8 - bytes) % 8;
  fill_n(&sieve_[sieveSize_], bytes, (uint8_t)0);

  segmentLow_ = stop_;
}

class PrimeGenerator : public Erat {
public:
  PrimeGenerator(uint64_t start, uint64_t stop);
  void fill(std::vector<uint64_t> &);
  void fill(std::vector<uint64_t> &primes, std::size_t *size);
  static uint64_t maxCachedPrime();

private:
  uint64_t low_ = 0;
  uint64_t sieveIdx_ = ~0ull;
  uint64_t prime_ = 0;
  PreSieve preSieve_;
  SievingPrimes sievingPrimes_;
  bool isInit_ = false;
  std::size_t getStartIdx() const;
  std::size_t getStopIdx() const;
  void initErat();
  void init(std::vector<uint64_t> &);
  void init(std::vector<uint64_t> &, std::size_t *);
  bool sieveSegment(std::vector<uint64_t> &);
  bool sieveSegment(std::vector<uint64_t> &, std::size_t *);
  void sieveSegment();
};

/// First 128 primes
const array<uint64_t, 128> smallPrimes = {
    2,   3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,
    53,  59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
    199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379,
    383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571,
    577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719};

/// Number of primes <= n
const array<uint8_t, 720> primePi = {
    0,   0,   1,   2,   2,   3,   3,   4,   4,   4,   4,   5,   5,   6,   6,
    6,   6,   7,   7,   8,   8,   8,   8,   9,   9,   9,   9,   9,   9,   10,
    10,  11,  11,  11,  11,  11,  11,  12,  12,  12,  12,  13,  13,  14,  14,
    14,  14,  15,  15,  15,  15,  15,  15,  16,  16,  16,  16,  16,  16,  17,
    17,  18,  18,  18,  18,  18,  18,  19,  19,  19,  19,  20,  20,  21,  21,
    21,  21,  21,  21,  22,  22,  22,  22,  23,  23,  23,  23,  23,  23,  24,
    24,  24,  24,  24,  24,  24,  24,  25,  25,  25,  25,  26,  26,  27,  27,
    27,  27,  28,  28,  29,  29,  29,  29,  30,  30,  30,  30,  30,  30,  30,
    30,  30,  30,  30,  30,  30,  30,  31,  31,  31,  31,  32,  32,  32,  32,
    32,  32,  33,  33,  34,  34,  34,  34,  34,  34,  34,  34,  34,  34,  35,
    35,  36,  36,  36,  36,  36,  36,  37,  37,  37,  37,  37,  37,  38,  38,
    38,  38,  39,  39,  39,  39,  39,  39,  40,  40,  40,  40,  40,  40,  41,
    41,  42,  42,  42,  42,  42,  42,  42,  42,  42,  42,  43,  43,  44,  44,
    44,  44,  45,  45,  46,  46,  46,  46,  46,  46,  46,  46,  46,  46,  46,
    46,  47,  47,  47,  47,  47,  47,  47,  47,  47,  47,  47,  47,  48,  48,
    48,  48,  49,  49,  50,  50,  50,  50,  51,  51,  51,  51,  51,  51,  52,
    52,  53,  53,  53,  53,  53,  53,  53,  53,  53,  53,  54,  54,  54,  54,
    54,  54,  55,  55,  55,  55,  55,  55,  56,  56,  56,  56,  56,  56,  57,
    57,  58,  58,  58,  58,  58,  58,  59,  59,  59,  59,  60,  60,  61,  61,
    61,  61,  61,  61,  61,  61,  61,  61,  62,  62,  62,  62,  62,  62,  62,
    62,  62,  62,  62,  62,  62,  62,  63,  63,  63,  63,  64,  64,  65,  65,
    65,  65,  66,  66,  66,  66,  66,  66,  66,  66,  66,  66,  66,  66,  66,
    66,  67,  67,  67,  67,  67,  67,  68,  68,  68,  68,  68,  68,  68,  68,
    68,  68,  69,  69,  70,  70,  70,  70,  71,  71,  71,  71,  71,  71,  72,
    72,  72,  72,  72,  72,  72,  72,  73,  73,  73,  73,  73,  73,  74,  74,
    74,  74,  74,  74,  75,  75,  75,  75,  76,  76,  76,  76,  76,  76,  77,
    77,  77,  77,  77,  77,  77,  77,  78,  78,  78,  78,  79,  79,  79,  79,
    79,  79,  79,  79,  80,  80,  80,  80,  80,  80,  80,  80,  80,  80,  81,
    81,  82,  82,  82,  82,  82,  82,  82,  82,  82,  82,  83,  83,  84,  84,
    84,  84,  84,  84,  85,  85,  85,  85,  86,  86,  86,  86,  86,  86,  87,
    87,  87,  87,  87,  87,  87,  87,  88,  88,  88,  88,  89,  89,  90,  90,
    90,  90,  91,  91,  91,  91,  91,  91,  91,  91,  91,  91,  91,  91,  92,
    92,  92,  92,  92,  92,  92,  92,  93,  93,  93,  93,  94,  94,  94,  94,
    94,  94,  94,  94,  95,  95,  95,  95,  96,  96,  96,  96,  96,  96,  97,
    97,  97,  97,  97,  97,  97,  97,  97,  97,  97,  97,  98,  98,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,  99,  99,  99,  99,  99,  99,  99,
    99,  100, 100, 100, 100, 100, 100, 101, 101, 101, 101, 101, 101, 101, 101,
    101, 101, 102, 102, 102, 102, 102, 102, 103, 103, 103, 103, 103, 103, 104,
    104, 105, 105, 105, 105, 105, 105, 106, 106, 106, 106, 106, 106, 106, 106,
    106, 106, 107, 107, 107, 107, 107, 107, 108, 108, 108, 108, 108, 108, 109,
    109, 110, 110, 110, 110, 110, 110, 111, 111, 111, 111, 111, 111, 112, 112,
    112, 112, 113, 113, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114,
    114, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 116, 116, 117, 117,
    117, 117, 118, 118, 118, 118, 118, 118, 119, 119, 119, 119, 119, 119, 120,
    120, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 122, 122,
    122, 122, 123, 123, 123, 123, 123, 123, 124, 124, 124, 124, 124, 124, 124,
    124, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 126, 126, 126, 126,
    126, 126, 126, 126, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 128};

} // namespace

namespace primesieve {

PrimeGenerator::PrimeGenerator(uint64_t start, uint64_t stop)
    : Erat(start, stop) {}

/// Used by iterator::prev_prime()
void PrimeGenerator::init(vector<uint64_t> &primes) {
  size_t size = primeCountApprox(start_, stop_);
  primes.reserve(size);

  if (start_ <= maxCachedPrime()) {
    size_t a = getStartIdx();
    size_t b = getStopIdx();

    primes.insert(primes.end(), smallPrimes.begin() + a,
                  smallPrimes.begin() + b);
  }

  initErat();
}

/// Used by iterator::next_prime()
void PrimeGenerator::init(vector<uint64_t> &primes, size_t *size) {
  if (start_ <= maxCachedPrime()) {
    size_t a = getStartIdx();
    size_t b = getStopIdx();

    *size = b - a;
    assert(*size <= primes.size());

    copy(smallPrimes.begin() + a, smallPrimes.begin() + b, primes.begin());
  }

  initErat();
}

void PrimeGenerator::initErat() {
  uint64_t startErat = maxCachedPrime() + 1;
  startErat = max(startErat, start_);
  isInit_ = true;

  if (startErat <= stop_) {
    int sieveSize = get_sieve_size();
    Erat::init(startErat, stop_, sieveSize, preSieve_);
    sievingPrimes_.init(this, preSieve_);
  }
}

uint64_t PrimeGenerator::maxCachedPrime() { return smallPrimes.back(); }

size_t PrimeGenerator::getStartIdx() const {
  size_t startIdx = 0;

  if (start_ > 1)
    startIdx = primePi[start_ - 1];

  return startIdx;
}

size_t PrimeGenerator::getStopIdx() const {
  size_t stopIdx = 0;

  if (stop_ < maxCachedPrime())
    stopIdx = primePi[stop_];
  else
    stopIdx = smallPrimes.size();

  return stopIdx;
}

void PrimeGenerator::sieveSegment() {
  uint64_t sqrtHigh = isqrt(segmentHigh_);

  sieveIdx_ = 0;
  low_ = segmentLow_;

  if (!prime_)
    prime_ = sievingPrimes_.next();

  while (prime_ <= sqrtHigh) {
    addSievingPrime(prime_);
    prime_ = sievingPrimes_.next();
  }

  Erat::sieveSegment();
}

/// Used by iterator::prev_prime()
bool PrimeGenerator::sieveSegment(vector<uint64_t> &primes) {
  if (!isInit_)
    init(primes);

  if (hasNextSegment()) {
    sieveSegment();
    return true;
  }

  return false;
}

/// Used by iterator::next_prime()
bool PrimeGenerator::sieveSegment(vector<uint64_t> &primes, size_t *size) {
  *size = 0;

  if (!isInit_) {
    init(primes, size);
    if (*size > 0)
      return false;
  }

  if (hasNextSegment()) {
    sieveSegment();
    return true;
  }

  // primesieve only supports primes < 2^64. In case the next
  // prime would be > 2^64 we simply return UINT64_MAX.
  if (stop_ >= numeric_limits<uint64_t>::max()) {
    primes[0] = ~0ull;
    *size = 1;
  }

  return false;
}

/// This method is used by iterator::prev_prime().
/// This method stores all primes inside [a, b] into the primes
/// vector. (b - a) is about sqrt(stop) so the memory usage is
/// quite large. Also after primesieve::iterator has iterated
/// over the primes inside [a, b] we need to generate new
/// primes which incurs an initialization overhead of O(sqrt(n)).
///
void PrimeGenerator::fill(vector<uint64_t> &primes) {
  while (sieveSegment(primes)) {
    while (sieveIdx_ < sieveSize_) {
      uint64_t bits = littleendian_cast<uint64_t>(&sieve_[sieveIdx_]);

      for (; bits != 0; bits &= bits - 1)
        primes.push_back(nextPrime(bits, low_));

      low_ += 8 * 30;
      sieveIdx_ += 8;
    }
  }
}

/// This method is used by iterator::next_prime().
/// This method stores only the next few primes (~ 200) in the
/// primes vector. Also for iterator::next_prime() there is no
/// recurring initialization overhead (unlike prev_prime()) for
/// this reason iterator::next_prime() runs up to 2x faster
/// than iterator::prev_prime().
///
void PrimeGenerator::fill(vector<uint64_t> &primes, size_t *size) {
  do {
    if (sieveIdx_ >= sieveSize_)
      if (!sieveSegment(primes, size))
        return;

    // Use local variables to prevent the compiler from
    // writing temporary results to memory.
    size_t i = 0;
    size_t maxSize = primes.size();
    assert(maxSize >= 64);
    uint64_t low = low_;
    uint8_t *sieve = sieve_;
    uint64_t sieveIdx = sieveIdx_;
    uint64_t sieveSize = sieveSize_;

    // Fill the buffer with at least (maxSize - 64) primes.
    // Each loop iteration can generate up to 64 primes
    // so we have to stop generating primes once there is
    // not enough space for 64 more primes.
    do {
      uint64_t bits = littleendian_cast<uint64_t>(&sieve[sieveIdx]);

      for (; bits != 0; bits &= bits - 1)
        primes[i++] = nextPrime(bits, low);

      low += 8 * 30;
      sieveIdx += 8;
    } while (i <= maxSize - 64 && sieveIdx < sieveSize);

    low_ = low;
    sieveIdx_ = sieveIdx;
    *size = i;
  } while (*size == 0);
}

uint64_t get_max_stop();

/// primesieve::iterator allows to easily iterate over primes both
/// forwards and backwards. Generating the first prime has a
/// complexity of O(r log log r) operations with r = n^0.5, after that
/// any additional prime is generated in amortized O(log n log log n)
/// operations. The memory usage is PrimePi(n^0.5) * 8 bytes.
///
class iterator {
public:
  iterator(uint64_t start = 0, uint64_t stop_hint = get_max_stop());
  iterator(const iterator &) = delete;
  iterator &operator=(const iterator &) = delete;
  iterator(iterator &&) noexcept;
  iterator &operator=(iterator &&) noexcept;

  ~iterator();
  void skipto(uint64_t start, uint64_t stop_hint = get_max_stop());
  uint64_t next_prime() {
    if (i_++ == last_idx_)
      generate_next_primes();
    return primes_[i_];
  }
  uint64_t prev_prime() {
    if (i_-- == 0)
      generate_prev_primes();
    return primes_[i_];
  }

private:
  std::size_t i_;
  std::size_t last_idx_;
  std::vector<uint64_t> primes_;
  uint64_t start_;
  uint64_t stop_;
  uint64_t stop_hint_;
  uint64_t dist_;
  std::unique_ptr<PrimeGenerator> primeGenerator_;
  void generate_next_primes();
  void generate_prev_primes();
}; // iterator class ends

class IteratorHelper {
public:
  // static void next(uint64_t* start,
  //                  uint64_t* stop,
  //                  uint64_t stopHint,
  //                  uint64_t* dist);

  // static void prev(uint64_t* start,
  //                  uint64_t* stop,
  //                  uint64_t stopHint,
  //                  uint64_t* dist);

  void next(uint64_t *start, uint64_t *stop, uint64_t stopHint,
            uint64_t *dist) {
    *start = checkedAdd(*stop, 1);
    uint64_t maxCachedPrime = PrimeGenerator::maxCachedPrime();

    if (*start < maxCachedPrime) {
      // When the stop number <= maxCachedPrime
      // primesieve::iterator uses the primes
      // cache instead of sieving and does not
      // even initialize Erat::init()
      *stop = maxCachedPrime;
      *dist = *stop - *start;
    } else {
      *dist = getNextDist(*start, *dist);
      *stop = checkedAdd(*start, *dist);

      if (useStopHint(*start, stopHint))
        *stop = checkedAdd(stopHint, maxPrimeGap(stopHint));
    }
  }

  void prev(uint64_t *start, uint64_t *stop, uint64_t stopHint,
            uint64_t *dist) {
    *stop = checkedSub(*start, 1);
    *dist = getPrevDist(*stop, *dist);
    *start = checkedSub(*stop, *dist);

    if (useStopHint(*start, *stop, stopHint))
      *start = checkedSub(stopHint, maxPrimeGap(stopHint));
  }

}; // iterator class ends

// iterator functions
iterator::~iterator() = default;

iterator::iterator(iterator &&) noexcept = default;

iterator &iterator::operator=(iterator &&) noexcept = default;

iterator::iterator(uint64_t start, uint64_t stop_hint) {
  skipto(start, stop_hint);
}

void iterator::skipto(uint64_t start, uint64_t stop_hint) {
  start_ = start;
  stop_ = start;
  stop_hint_ = stop_hint;
  i_ = 0;
  last_idx_ = 0;
  dist_ = 0;
  clear(primeGenerator_);
  primes_.clear();
}

void iterator::generate_next_primes() {
  while (true) {
    if (!primeGenerator_) {
      IteratorHelper::next(&start_, &stop_, stop_hint_, &dist_);
      auto p = new PrimeGenerator(start_, stop_);
      primeGenerator_.reset(p);
      primes_.resize(256);
    }

    primeGenerator_->fill(primes_, &last_idx_);

    // There are 3 different cases here:
    // 1) The primes array contains a few primes (<= 256).
    //    In this case we return the primes to the user.
    // 2) The primes array is empty because the next
    //    prime > stop. In this case we reset the
    //    primeGenerator object, increase the start & stop
    //    numbers and sieve the next segment.
    // 3) The next prime > 2^64. In this case the primes
    //    array contains an error code (UINT64_MAX) which
    //    is returned to the user.
    if (last_idx_ == 0)
      clear(primeGenerator_);
    else
      break;
  }

  i_ = 0;
  last_idx_--;
}

void iterator::generate_prev_primes() {
  if (primeGenerator_)
    start_ = primes_.front();

  primes_.clear();

  while (primes_.empty()) {
    IteratorHelper::prev(&start_, &stop_, stop_hint_, &dist_);
    if (start_ <= 2)
      primes_.push_back(0);
    auto p = new PrimeGenerator(start_, stop_);
    primeGenerator_.reset(p);
    primeGenerator_->fill(primes_);
    clear(primeGenerator_);
  }

  last_idx_ = primes_.size() - 1;
  i_ = last_idx_;
}
// iterator functions ends

// auto iter;
vector<int32_t> generate_n_primes(uint64_t n) {
  uint64_t start = 0;
  vector<int32_t> primes;
  if (n == 0)
    return primes;
  if (start > 0)
    start--;

  std::size_t size = primes.size() + (std::size_t)n;
  primes.reserve(size);
  using V = typename vector<int32_t>::value_type;

  double x = (double)start;
  x = std::max<double>(10.0, x);
  uint64_t logx = (uint64_t)std::log(x);
  uint64_t dist = n * (logx + 1);
  uint64_t stop = start + dist;

  iterator it(start, stop);
  int32_t prime = it.next_prime();
  for (; n > 0; n--, prime = it.next_prime())
    primes.push_back((V)prime);

  if (~prime == 0)
    throw primesieve_error("cannot generate primes > 2^64");
  return primes;
}

/// primeCountApprox(x) >= pi(x)
inline std::size_t prime_count_approx(uint64_t start, uint64_t stop) {
  if (start > stop)
    return 0;
  if (stop <= 10)
    return 4;

  // pi(x) <= x / (log(x) - 1.1) + 5, for x >= 4
  double x = (double)stop;
  double logx = std::log(x);
  double div = logx - 1.1;
  double pix = (stop - start) / div + 5;

  return (std::size_t)pix;
}

vector<int32_t> generate_primes(uint64_t stop) {
  uint64_t start = 0;
  vector<int32_t> primes;
  if (start > 0)
    start--;
  if (~stop == 0)
    stop--;

  if (start < stop) {
    using V = typename vector<int32_t>::value_type;
    std::size_t size = primes.size() + prime_count_approx(start, stop);
    primes.reserve(size);

    iterator it(start, stop);
    uint64_t prime = it.next_prime();
    for (; prime <= stop; prime = it.next_prime())
      primes.push_back((V)prime);
  }
  return primes;
}

/// Cache phi(x, a) results if a < MAX_A
const int MAX_A = 100;

class PhiCache {
public:
  PhiCache(vector<int32_t> &primes, PiTable &pi) : primes_(primes), pi_(pi) {}

  /// Calculate phi(x, a) using the recursive formula:
  /// phi(x, a) = phi(x, a - 1) - phi(x / primes_[a], a - 1)
  ///
  template <int SIGN> int64_t phi(int64_t x, int64_t a) {
    if (x <= primes_[a])
      return SIGN;
    else if (is_phi_tiny(a))
      return phi_tiny(x, a) * SIGN;
    else if (is_pix(x, a))
      return (pi_[x] - a + 1) * SIGN;
    else if (is_cached(x, a))
      return cache_[a][x] * SIGN;

    int64_t sqrtx = isqrt(x);
    int64_t pi_sqrtx = a;
    int64_t c = PhiTiny::get_c(sqrtx);
    int64_t sum = 0;

    if (sqrtx < pi_.size())
      pi_sqrtx = min(pi_[sqrtx], a);

    // Move out of the loop the calculations where phi(x2, i) = 1
    // phi(x, a) = 1 if primes_[a] >= x
    // x2 = x / primes_[i + 1]
    // phi(x2, i) = 1 if primes_[i] >= x / primes_[i + 1]
    // phi(x2, i) = 1 if primes_[i] >= sqrt(x)
    // phi(x2, i) = 1 if i >= pi(sqrt(x))
    // \sum_{i = pi(sqrt(x))}^{a - 1} phi(x2, i) = a - pi(sqrt(x))
    //
    sum += (pi_sqrtx - a) * SIGN;
    sum += phi_tiny(x, c) * SIGN;

    for (int64_t i = c; i < pi_sqrtx; i++) {
      int64_t x2 = fast_div(x, primes_[i + 1]);

      if (is_pix(x2, i))
        sum += (pi_[x2] - i + 1) * -SIGN;
      else
        sum += phi<-SIGN>(x2, i);
    }

    update_cache(x, a, sum);

    return sum;
  }

private:
  using T = uint16_t;
  array<vector<T>, MAX_A> cache_;
  vector<int32_t> &primes_;
  PiTable &pi_;

  void update_cache(uint64_t x, uint64_t a, int64_t sum) {
    if (a < cache_.size() && x <= numeric_limits<T>::max()) {
      if (x >= cache_[a].size())
        cache_[a].resize(x + 1, 0);

      cache_[a][x] = (T)abs(sum);
    }
  }

  bool is_pix(int64_t x, int64_t a) const {
    return x < pi_.size() && x < isquare(primes_[a + 1]);
  }

  bool is_cached(uint64_t x, uint64_t a) const {
    return a < cache_.size() && x < cache_[a].size() && cache_[a][x];
  }
};

/// Partial sieve function (a.k.a. Legendre-sum).
/// phi(x, a) counts the numbers <= x that are not divisible
/// by any of the first a primes.
///
int64_t phi(int64_t x, int64_t a, int threads) {
  if (x < 1)
    return 0;
  if (a > x)
    return 1;
  if (a < 1)
    return x;

  int64_t sum = 0;

  if (is_phi_tiny(a))
    sum = phi_tiny(x, a);
  else {
    auto primes = generate_n_primes(a);

    if (primes[a] >= x)
      sum = 1;
    else {
      // use large pi(x) lookup table for speed
      int64_t sqrtx = isqrt(x);
      PiTable pi(max(sqrtx, (int64_t)primes[a]));
      PhiCache cache(primes, pi);

      int64_t c = PhiTiny::get_c(sqrtx);
      int64_t pi_sqrtx = min(pi[sqrtx], a);
      int64_t thread_threshold = ipow(10ll, 10);
      threads = 1; // ideal_num_threads(threads, x, thread_threshold);

      sum = phi_tiny(x, c) - a + pi_sqrtx;

#pragma omp parallel for num_threads(threads)                                  \
    schedule(dynamic, 16) firstprivate(cache) reduction(+ : sum)
      for (int64_t i = c; i < pi_sqrtx; i++)
        sum += cache.phi<-1>(x / primes[i + 1], i);
    }
  }

  return sum;
}

/// Count the number of primes <= x using Legendre's formula.
/// Run time: O(x)
/// Memory usage: O(x^(1/2))
///
int64_t pi_legendre(int64_t x, int threads) {
  if (x < 2)
    return 0;

  int64_t a = pi_legendre(isqrt(x), 1);
  int64_t sum = phi(x, a, threads) + a - 1;

  return sum;
}

/// Generate a vector with Möbius function values.
/// This implementation is based on code by Rick Sladkey:
/// https://mathoverflow.net/q/99545
///
vector<int32_t> generate_moebius(int64_t max) {
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  vector<int32_t> mu(size, 1);

  for (int64_t i = 2; i <= sqrt; i++) {
    if (mu[i] == 1) {
      for (int64_t j = i; j < size; j += i)
        mu[j] *= (int32_t)-i;
      for (int64_t j = i * i; j < size; j += i * i)
        mu[j] = 0;
    }
  }

  for (int64_t i = 2; i < size; i++) {
    if (mu[i] == i)
      mu[i] = 1;
    else if (mu[i] == -i)
      mu[i] = -1;
    else if (mu[i] < 0)
      mu[i] = 1;
    else if (mu[i] > 0)
      mu[i] = -1;
  }

  return mu;
}

/// Generate a vector with the least prime factors
/// of the integers <= max
///
vector<int32_t> generate_lpf(int64_t max) {
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  vector<int32_t> lpf(size, 1);

  for (int64_t i = 2; i <= sqrt; i++)
    if (lpf[i] == 1)
      for (int64_t j = i * i; j < size; j += i)
        if (lpf[j] == 1)
          lpf[j] = (int32_t)i;

  for (int64_t i = 2; i < size; i++)
    if (lpf[i] == 1)
      lpf[i] = (int32_t)i;

  // phi(x / 1, c) contributes to the sum in the
  // Lagarias-Miller-Odlyzko prime counting algorithm,
  // thus set lpf[1] = MAX (normally lpf[1] = 1)
  if (lpf.size() > 1)
    lpf[1] = numeric_limits<int32_t>::max();

  return lpf;
}

template <int SIGN, typename Primes>
int128_t phi_sum128(int64_t x, int64_t a, Primes &&primes) {
  int128_t sum = 0;

  for (; a > 0; a--) {
    if (x <= primes[a])
      return sum + SIGN;

    int64_t x2 = fast_div(x, primes[a]);
    sum += phi_sum128<-SIGN>(x2, a - 1, primes) * primes[a];
  }

  int128_t n = x;
  int128_t fx = (n * (n + 1)) >> 1;
  sum += fx * SIGN;

  return sum;
}

template <int SIGN, typename Primes>
int256_t phi_sum256(int128_t x, int64_t a, Primes &&primes) {
  int256_t sum = 0;

  for (; a > 0; a--) {
    if (x <= primes[a])
      return sum + SIGN;

    int128_t x2 = fast_div(x, primes[a]);
    int256_t phi_sum;

    if (x2 <= numeric_limits<int64_t>::max())
      phi_sum = phi_sum128<-SIGN>((int64_t)x2, a - 1, primes);
    else
      phi_sum = phi_sum256<-SIGN>(x2, a - 1, primes);

    sum += phi_sum * primes[a];
  }

  int256_t n = x;
  int256_t fx = (n * (n + 1)) >> 1;
  sum += fx * SIGN;

  return sum;
}

int256_t phi_sum(int128_t x, int64_t a) {
  if (x <= numeric_limits<int64_t>::max())
    return phi_sum((int64_t)x, a);

  if (a < 10)
    return phi_sum256<1>(x, a, small_primes_);
  else
    return phi_sum256<1>(x, a, generate_n_primes(a));
}

int128_t prime_sum_tiny(int64_t x) {
  int64_t prime = 0;
  int128_t sum = 0;
  // iterator iter(0, x);

  // while ((prime = iter.next_prime()) <= x)
  //   sum += prime;

  return sum;
}

/// The aligned_vector class aligns each of its
/// elements on a new cache line in order to avoid
/// false sharing (cache trashing) when multiple
/// threads write to adjacent elements
///
template <typename T> class aligned_vector {
public:
  aligned_vector(std::size_t size) : vect_(size) {}
  std::size_t size() const { return vect_.size(); }
  T &operator[](std::size_t pos) { return vect_[pos].val[0]; }

private:
  struct align_t {
    T val[CACHE_LINE_SIZE / sizeof(T)];
  };
  std::vector<align_t> vect_;
};

/// Calculate the thread sieving distance. The idea is to
/// gradually increase the thread_distance in order to
/// keep all CPU cores busy.

template <typename T, typename X>
T P2_OpenMP_thread(X x, int64_t y, int64_t z, int64_t thread_distance,
                   int64_t thread_num, int64_t low, T &prime_sum, T &correct) {
  prime_sum = 0;
  correct = 0;
  low += thread_distance * thread_num;
  z = min(low + thread_distance, z);

  int64_t sqrtx = isqrt(x);
  int64_t start = (int64_t)max((int64_t)(x / z), (int64_t)y);
  int64_t stop = (int64_t)min((int64_t)(x / low), (int64_t)sqrtx);
  int64_t x_div_prime = 0;

  iterator rit(stop + 1, start);
  iterator it(low - 1, z);

  int64_t next = it.next_prime();
  int64_t prime = rit.prev_prime();
  T P2_thread = 0;

  while (prime > start && (x_div_prime = (int64_t)(x / prime)) < z) {
    // Sum the primes <= x / prime
    while (next <= x_div_prime) {
      if (next > y && next <= sqrtx) {
        P2_thread -= prime_sum * next;
        correct -= next;
      }

      prime_sum += next;
      next = it.next_prime();
    }

    P2_thread += prime_sum * prime;
    correct += prime;
    prime = rit.prev_prime();
  }

  // Sum the primes < z
  while (next < z) {
    if (next > y && next <= sqrtx) {
      P2_thread -= prime_sum * next;
      correct -= next;
    }

    prime_sum += next;
    next = it.next_prime();
  }

  return P2_thread;
}

/// P2(x, y) sums the numbers <= x that have exactly 2 prime
/// factors each exceeding the a-th prime, a = pi(y).
/// Space complexity: O((x / y)^(1/2)).
///
template <typename T>
typename next_larger_type<T>::type P2_OpenMP_master(T x, int64_t y,
                                                    int threads) {
  if (x < 4)
    return 0;

  int64_t a = pi_legendre(y, threads);
  int64_t b = pi_legendre(isqrt(x), threads);

  if (a >= b)
    return 0;

  int64_t low = 2;
  int64_t z = (int64_t)(x / max(y, (int64_t)1));
  int64_t min_distance = 1 << 23;
  int64_t thread_distance = min_distance;

  using res_t = typename next_larger_type<T>::type;

  aligned_vector<res_t> prime_sums(threads);
  aligned_vector<res_t> correct(threads);

  res_t p2 = 0;
  res_t prime_sum = prime_sum_tiny(y);

  while (low < z) {
    int64_t segments = ceil_div(z - low, thread_distance);
    threads = 1; // in_between(1, threads, segments);
    double time = get_time();

#pragma omp parallel for num_threads(threads) reduction(+ : p2)
    for (int i = 0; i < threads; i++)
      p2 += P2_OpenMP_thread(x, y, z, thread_distance, i, low, prime_sums[i],
                             correct[i]);

    for (int i = 0; i < threads; i++) {
      p2 += prime_sum * correct[i];
      prime_sum += prime_sums[i];
    }

    // low += thread_distance; //* threads;
    // balanceLoad(&thread_distance, low, z, threads, time);

    // if (is_print()) {
    // double percent = get_percent(low, z);
    // cout << "\rStatus: " << fixed << setprecision(get_status_precision(x))
    //      << percent << '%' << flush;
    // }
  }

  return p2;
}

int256_t P2(int128_t x, int64_t y, int threads) {
  cout << "" << endl;
  cout << "=== P2(x, y) ===" << endl;
  cout << "Computation of the 2nd partial sieve function" << endl;
  // print(x, y, threads);

  double time = get_time();
  int256_t p2;

  // uses less memory
  if (x <= numeric_limits<int64_t>::max())
    p2 = P2_OpenMP_master((int64_t)x, y, threads);
  else
    p2 = P2_OpenMP_master(x, y, threads);

  // print("P2", p2, time);
  return p2;
}

/// Calculate the number of primes below x using the
/// Lagarias-Miller-Odlyzko algorithm.
/// Run time: O(x^(2/3)) operations, O(x^(1/3)) space.
///
int256_t pi_lmo1(int128_t x) {
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x);
  int64_t pi_y = pi_legendre(y, 1);
  int64_t c = PhiTiny::get_c(y);
  int256_t S1 = 0;
  int256_t S2 = 0;

  vector<int32_t> primes = generate_primes(y);
  vector<int32_t> lpf = generate_lpf(y);
  vector<int32_t> mu = generate_moebius(y);

  // Calculate the contribution of the ordinary leaves
  for (int64_t n = 1; n <= y; n++)
    if (lpf[n] > primes[c])
      S1 += phi_sum(x / n, c) * (mu[n] * n);

  // Calculate the contribution of the special leaves
  for (int64_t b = c + 1; b < pi_y; b++)
    for (int128_t m = (y / primes[b]) + 1; m <= y; m++)
      if (lpf[m] > primes[b])
        S2 -= phi_sum(x / (primes[b] * m), b - 1) * (mu[m] * m * primes[b]);

  int256_t phi = S1 + S2;
  int256_t p2 = P2(x, y, 1);
  int256_t sum = phi + prime_sum_tiny(y) - 1 - p2;

  return sum;
}

} // namespace

int main() {
  auto start_time = clock();
  primesum::int128_t n = 1e15;
  cout << "n = " << (double)n << endl;
  primesum::int256_t res = primesum::pi_lmo1(n);
  cout << "result: " << res << endl;
  cout << "time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s"
       << endl;
  return 0;
}