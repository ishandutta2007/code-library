// #include <primesum.hpp>
// #include <primesum-internal.hpp>
// #include <generate.hpp>
// #include <PhiTiny.hpp>

#include <bits/stdc++.h>
using namespace std;
///
/// @file  macros.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MACROS_HPP
#define MACROS_HPP

#ifndef __has_attribute
#define __has_attribute(x) 0
#endif

#ifndef __has_builtin
#define __has_builtin(x) 0
#endif

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(x) 0
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

#if __cplusplus >= 202002L && __has_cpp_attribute(unlikely)
#define if_unlikely(x) if (x) [[unlikely]]
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

/// Use [[maybe_unused]] from C++17 once widely supported
#if defined(NDEBUG)
#define MAYBE_UNUSED(x)
#else
#define MAYBE_UNUSED(x) x
#endif

#endif

#ifndef IMATH_HPP
#define IMATH_HPP

#ifndef ISQRT_HPP
#define ISQRT_HPP
///
/// @file   int128_t.hpp
/// @brief  Support for int128_t, uint128_t types.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef INT128_T_HPP
#define INT128_T_HPP

#include <limits>
#include <stdint.h>
#include <type_traits>

/// The __int128_t type (GCC/Clang) is not well supported by
/// the C++ standard library (in 2016) so we have to define
/// some functions ourselves. We also define typedefs so we
/// can use int128_t instead of __int128_t. Once this is done
/// int128_t can be used like a regular integer type.
///
#if !defined(INT128_MAX) && defined(__SIZEOF_INT128__)

#include <ostream>
#include <string>

namespace primesum {

using int128_t = __int128_t;
using uint128_t = __uint128_t;

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

} // namespace primesum

#endif

namespace primesum {

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
} // namespace primesum

#endif

// namespace primesum {
// namespace prt {

// template <typename T> struct numeric_limits {
//   static constexpr T max() { return std::numeric_limits<T>::max(); }
// };

// template <> struct numeric_limits<int128_t> {
//   static constexpr int128_t min() { return ((int128_t)1) << 127; }
//   static constexpr int128_t max() { return ~min(); }
// };

// template <> struct numeric_limits<uint128_t> {
//   static constexpr uint128_t min() { return 0; }
//   static constexpr uint128_t max() { return ~min(); }
// };

// template <typename T> struct is_integral {
//   enum {
//     value = std::is_integral<T>::value || std::is_same<T, int128_t>::value ||
//             std::is_same<T, uint128_t>::value
//   };
// };

// template <typename T> struct is_signed {
//   enum { value = std::is_signed<T>::value || std::is_same<T, int128_t>::value
//   };
// };

// template <typename T> struct is_unsigned {
//   enum {
//     value = std::is_unsigned<T>::value || std::is_same<T, uint128_t>::value
//   };
// };

// } // namespace prt
// } // namespace primesum

namespace primesum {

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
  return lo == hi ? lo
                  : ((x / MID < MID) ? sqrt_helper<T>(x, lo, MID - 1)
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

} // namespace primesum

#endif

namespace primesum {

inline int64_t isquare(int64_t x) { return x * x; }

template <typename A, typename B> inline A ceil_div(A a, B b) {
  assert(b > 0);
  return (A)((a + b - 1) / b);
}

template <typename T> inline T number_of_bits(T) {
  return (T)std::numeric_limits<T>::digits;
}

template <typename T> inline T next_power_of_2(T x) {
  if (x == 0)
    return 1;

  x--;
  for (T i = 1; i < number_of_bits(x); i += i)
    x |= (x >> i);

  return ++x;
}

template <typename T> inline T prev_power_of_2(T x) {
  for (T i = 1; i < number_of_bits(x); i += i)
    x |= (x >> i);

  return x - (x >> 1);
}

template <typename T> inline int ilog(T x) { return (int)std::log((double)x); }

template <typename T> inline T ipow(T x, int n) {
  T r = 1;
  for (int i = 0; i < n; i++)
    r *= x;

  return r;
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

template <typename T1, typename T2, typename T3>
inline T2 in_between(T1 min, T2 x, T3 max) {
  if (x < min)
    return (T2)min;
  if (x > max)
    return (T2)max;

  return x;
}

template <typename T1, typename T2>
inline T2 pi_bsearch(const std::vector<T1> &primes, T2 x) {
  assert(primes.size() < 2 || primes[1] == 2);
  auto start = primes.begin() + 1;
  return (T2)(std::upper_bound(start, primes.end(), x) - start);
}

} // namespace primesum

#endif

#ifndef INT256_T_HPP
#define INT256_T_HPP

#ifndef INT128_T_HPP
#define INT128_T_HPP

#if !defined(INT128_MAX) && defined(__SIZEOF_INT128__)

namespace primesum {

using int128_t = __int128_t;
using uint128_t = __uint128_t;

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

} // namespace primesum

#endif

#endif

namespace primesum {

class int256_t {
public:
  int256_t() : low(0), high(0) {}

  template <typename T, typename = typename std::enable_if<
                            prt::is_integral<T>::value>::type>
  int256_t(T x) : low(x), high((x < 0) ? -1 : 0) {}

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
    return (*this < 0) ? -static_cast<std::int8_t>(
                             (low - 1) ^ prt::numeric_limits<uint128_t>::max())
                       : static_cast<std::int8_t>(low);
  }

  operator std::int16_t() const {
    return (*this < 0) ? -static_cast<std::int16_t>(
                             (low - 1) ^ prt::numeric_limits<uint128_t>::max())
                       : static_cast<std::int16_t>(low);
  }

  operator std::int32_t() const {
    return (*this < 0) ? -static_cast<std::int32_t>(
                             (low - 1) ^ prt::numeric_limits<uint128_t>::max())
                       : static_cast<std::int32_t>(low);
  }

  operator std::int64_t() const {
    return (*this < 0) ? -static_cast<std::int64_t>(
                             (low - 1) ^ prt::numeric_limits<uint128_t>::max())
                       : static_cast<std::int64_t>(low);
  }

  operator int128_t() const {
    return (*this < 0) ? -static_cast<int128_t>(
                             (low - 1) ^ prt::numeric_limits<uint128_t>::max())
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

} // namespace primesum

#ifdef _OPENMP

namespace primesum {

/// Requires OpenMP >= 4.0
#pragma omp declare reduction(+ : int256_t : omp_out += omp_in)

} // namespace primesum

#endif

#endif
// int256_t.hpp ends

#ifndef PITABLE_HPP
#define PITABLE_HPP

///
/// @file  popcnt.hpp
/// @brief Functions to count the number of 1 bits inside
///        a 64-bit word.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

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

namespace primesum {

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

} // namespace primesum

#endif

#ifndef PHITINY_HPP
#define PHITINY_HPP

namespace primesum {

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

} // namespace primesum

#endif

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

// #include <int128_t.hpp>

namespace primesum {

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

} // namespace primesum

#endif

///
/// @file   primesieve.hpp
/// @brief  primesieve C++ API. primesieve is a library for fast
///         prime number generation, in case an error occurs a
///         primesieve::primesieve_error exception (derived form
///         std::runtime_error) is thrown.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License.
///

#ifndef PRIMESIEVE_HPP
#define PRIMESIEVE_HPP

#define PRIMESIEVE_VERSION "7.6"
#define PRIMESIEVE_VERSION_MAJOR 7
#define PRIMESIEVE_VERSION_MINOR 6

// #include <primesieve/iterator.hpp>
// #include <primesieve/primesieve_error.hpp>
// #include <primesieve/StorePrimes.hpp>

/// Contains primesieve's C++ functions and classes.
namespace primesieve {

/// Store the primes <= stop in the primes vector.
template <typename T>
inline void generate_primes(uint64_t stop, std::vector<T> *primes) {
  if (primes)
    store_primes(0, stop, *primes);
}

/// Store the primes within the interval [start, stop]
/// in the primes vector.
///
template <typename T>
inline void generate_primes(uint64_t start, uint64_t stop,
                            std::vector<T> *primes) {
  if (primes)
    store_primes(start, stop, *primes);
}

/// Store the first n primes in the primes vector.
template <typename T>
inline void generate_n_primes(uint64_t n, std::vector<T> *primes) {
  if (primes)
    store_n_primes(n, 0, *primes);
}

/// Store the first n primes >= start in the primes vector.
template <typename T>
inline void generate_n_primes(uint64_t n, uint64_t start,
                              std::vector<T> *primes) {
  if (primes)
    store_n_primes(n, start, *primes);
}

/// Find the nth prime.
/// By default all CPU cores are used, use
/// primesieve::set_num_threads(int threads) to change the
/// number of threads.
///
/// Note that each call to nth_prime(n, start) incurs an
/// initialization overhead of O(sqrt(start)) even if n is tiny.
/// Hence it is not a good idea to use nth_prime() repeatedly in a
/// loop to get the next (or previous) prime. For this use case it
/// is better to use a primesieve::iterator which needs to be
/// initialized only once.
///
/// @param n  if n = 0 finds the 1st prime >= start, <br/>
///           if n > 0 finds the nth prime > start, <br/>
///           if n < 0 finds the nth prime < start (backwards).
///
uint64_t nth_prime(int64_t n, uint64_t start = 0);

/// Count the primes within the interval [start, stop].
/// By default all CPU cores are used, use
/// primesieve::set_num_threads(int threads) to change the
/// number of threads.
///
/// Note that each call to count_primes() incurs an initialization
/// overhead of O(sqrt(stop)) even if the interval [start, stop]
/// is tiny. Hence if you have written an algorithm that makes
/// many calls to count_primes() it may be preferable to use
/// a primesieve::iterator which needs to be initialized only once.
///
uint64_t count_primes(uint64_t start, uint64_t stop);

/// Count the twin primes within the interval [start, stop].
/// By default all CPU cores are used, use
/// primesieve::set_num_threads(int threads) to change the
/// number of threads.
///
uint64_t count_twins(uint64_t start, uint64_t stop);

/// Count the prime triplets within the interval [start, stop].
/// By default all CPU cores are used, use
/// primesieve::set_num_threads(int threads) to change the
/// number of threads.
///
uint64_t count_triplets(uint64_t start, uint64_t stop);

/// Count the prime quadruplets within the interval [start, stop].
/// By default all CPU cores are used, use
/// primesieve::set_num_threads(int threads) to change the
/// number of threads.
///
uint64_t count_quadruplets(uint64_t start, uint64_t stop);

/// Count the prime quintuplets within the interval [start, stop].
/// By default all CPU cores are used, use
/// primesieve::set_num_threads(int threads) to change the
/// number of threads.
///
uint64_t count_quintuplets(uint64_t start, uint64_t stop);

/// Count the prime sextuplets within the interval [start, stop].
/// By default all CPU cores are used, use
/// primesieve::set_num_threads(int threads) to change the
/// number of threads.
///
uint64_t count_sextuplets(uint64_t start, uint64_t stop);

/// Print the primes within the interval [start, stop]
/// to the standard output.
///
void print_primes(uint64_t start, uint64_t stop);

/// Print the twin primes within the interval [start, stop]
/// to the standard output.
///
void print_twins(uint64_t start, uint64_t stop);

/// Print the prime triplets within the interval [start, stop]
/// to the standard output.
///
void print_triplets(uint64_t start, uint64_t stop);

/// Print the prime quadruplets within the interval [start, stop]
/// to the standard output.
///
void print_quadruplets(uint64_t start, uint64_t stop);

/// Print the prime quintuplets within the interval [start, stop]
/// to the standard output.
///
void print_quintuplets(uint64_t start, uint64_t stop);

/// Print the prime sextuplets within the interval [start, stop]
/// to the standard output.
///
void print_sextuplets(uint64_t start, uint64_t stop);

/// Returns the largest valid stop number for primesieve.
/// @return 2^64-1 (UINT64_MAX).
///
uint64_t get_max_stop();

/// Get the current set sieve size in KiB.
int get_sieve_size();

/// Get the current set number of threads.
int get_num_threads();

/// Set the sieve size in KiB (kibibyte).
/// The best sieving performance is achieved with a sieve size
/// of your CPU's L1 or L2 cache size (per core).
/// @pre sieve_size >= 8 && <= 4096.
///
void set_sieve_size(int sieve_size);

/// Set the number of threads for use in
/// primesieve::count_*() and primesieve::nth_prime().
/// By default all CPU cores are used.
///
void set_num_threads(int num_threads);

/// Get the primesieve version number, in the form “i.j”.
std::string primesieve_version();
} // namespace primesieve

#endif

///
/// @file  PreSieve.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRESIEVE_HPP
#define PRESIEVE_HPP
namespace primesieve {

/// PreSieve objects are used to pre-sieve multiples of small primes
/// e.g. <= 19 to speed up the sieve of Eratosthenes. The idea is to
/// allocate an array (buffer_) and remove the multiples of small
/// primes from it at initialization. Then whilst sieving, the
/// buffer_ array is copied to the sieve array at the beginning of
/// each new segment to pre-sieve the multiples of small
/// primes <= maxPrime_. Pre-sieving speeds up my sieve of Eratosthenes
/// implementation by about 20 percent when sieving < 10^10.
///
/// <b> Memory Usage </b>
///
/// - PreSieve objects use: primeProduct(maxPrime_) / 30 bytes of memory
/// - PreSieve multiples of primes <=  7 uses    7    bytes
/// - PreSieve multiples of primes <= 11 uses   77    bytes
/// - PreSieve multiples of primes <= 13 uses 1001    bytes
/// - PreSieve multiples of primes <= 17 uses   17.02 kilobytes
/// - PreSieve multiples of primes <= 19 uses  323.32 kilobytes
/// - PreSieve multiples of primes <= 23 uses    7.44 megabytes
///
class PreSieve {
public:
  void init(uint64_t, uint64_t);
  uint64_t getMaxPrime() const { return maxPrime_; }
  void copy(uint8_t *, uint64_t, uint64_t) const;

private:
  uint64_t maxPrime_ = 0;
  uint64_t primeProduct_ = 0;
  uint64_t size_ = 0;
  std::vector<uint8_t> buffer_;
  void initBuffer(uint64_t, uint64_t);
};

} // namespace primesieve

#endif

#ifndef PRIMESIEVE_CLASS_HPP
#define PRIMESIEVE_CLASS_HPP

namespace primesieve {

using counts_t = std::array<uint64_t, 6>;
class ParallelSieve;

enum {
  COUNT_PRIMES = 1 << 0,
  COUNT_TWINS = 1 << 1,
  COUNT_TRIPLETS = 1 << 2,
  COUNT_QUADRUPLETS = 1 << 3,
  COUNT_QUINTUPLETS = 1 << 4,
  COUNT_SEXTUPLETS = 1 << 5,
  PRINT_PRIMES = 1 << 6,
  PRINT_TWINS = 1 << 7,
  PRINT_TRIPLETS = 1 << 8,
  PRINT_QUADRUPLETS = 1 << 9,
  PRINT_QUINTUPLETS = 1 << 10,
  PRINT_SEXTUPLETS = 1 << 11,
  PRINT_STATUS = 1 << 12
};

class PrimeSieve {
public:
  PrimeSieve();
  PrimeSieve(ParallelSieve *);
  virtual ~PrimeSieve();
  // Getters
  uint64_t getStart() const;
  uint64_t getStop() const;
  uint64_t getDistance() const;
  int getSieveSize() const;
  double getSeconds() const;
  PreSieve &getPreSieve();
  // Setters
  void setStart(uint64_t);
  void setStop(uint64_t);
  void updateStatus(uint64_t);
  void setSieveSize(int);
  void setFlags(int);
  void addFlags(int);
  // Bool is*
  bool isCount(int) const;
  bool isCountPrimes() const;
  bool isCountkTuplets() const;
  bool isPrint() const;
  bool isPrint(int) const;
  bool isPrintPrimes() const;
  bool isPrintkTuplets() const;
  bool isFlag(int) const;
  bool isFlag(int, int) const;
  bool isStatus() const;
  // Sieve
  virtual void sieve();
  void sieve(uint64_t, uint64_t);
  void sieve(uint64_t, uint64_t, int);
  // nth prime
  uint64_t nthPrime(uint64_t);
  uint64_t nthPrime(int64_t, uint64_t);
  // Count
  counts_t &getCounts();
  uint64_t getCount(int) const;
  uint64_t countPrimes(uint64_t, uint64_t);

protected:
  /// Sieve primes >= start_
  uint64_t start_ = 0;
  /// Sieve primes <= stop_
  uint64_t stop_ = 0;
  /// Time elapsed of sieve()
  double seconds_ = 0;
  /// Sieving status in percent
  double percent_ = 0;
  /// Prime number and prime k-tuplet counts
  counts_t counts_;
  void reset();
  void setStatus(double);

private:
  uint64_t sievedDistance_ = 0;
  uint64_t updateDistance_ = 0;
  /// Default flags
  int flags_ = COUNT_PRIMES;
  /// Sieve size in KiB
  int sieveSize_ = 0;
  /// Status updates must be synchronized by main thread
  ParallelSieve *parent_ = nullptr;
  PreSieve preSieve_;
  void processSmallPrimes();
  static void printStatus(double, double);
};

} // namespace primesieve

#endif

///
/// @file   ParallelSieve.hpp
/// @brief  The ParallelSieve class provides an easy API for
///         multi-threaded prime sieving.
///
/// Copyright (C) 2019 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PARALLELSIEVE_HPP
#define PARALLELSIEVE_HPP

// #include "PrimeSieve.hpp"

namespace primesieve {

class ParallelSieve : public PrimeSieve {
public:
  using PrimeSieve::sieve;

  ParallelSieve();
  static int getMaxThreads();
  int getNumThreads() const;
  int idealNumThreads() const;
  void setNumThreads(int numThreads);
  bool tryUpdateStatus(uint64_t);
  virtual void sieve();

private:
  std::mutex mutex_;
  int numThreads_ = 0;
  uint64_t getThreadDistance(int) const;
  uint64_t align(uint64_t) const;
};

} // namespace primesieve

#endif

///
/// @file   pmath.hpp
/// @brief  Auxiliary math functions for primesieve.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PMATH_HPP
#define PMATH_HPP

namespace {

template <typename X, typename Y> inline X ceilDiv(X x, Y y) {
  return (X)((x + y - 1) / y);
}

template <typename T> constexpr bool isPow2(T x) {
  return x != 0 && (x & (x - 1)) == 0;
}

template <typename T> constexpr T numberOfBits(T) {
  return (T)std::numeric_limits<typename std::make_unsigned<T>::type>::digits;
}

template <typename T> inline T floorPow2(T x) {
  for (T i = 1; i < numberOfBits(x); i += i)
    x |= (x >> i);

  return x - (x >> 1);
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

#if __cplusplus >= 201402L

/// C++14 compile time square root using binary search
template <typename T> constexpr T ctSqrt(T x, T lo, T hi) {
  if (lo == hi)
    return lo;

  const T mid = (lo + hi + 1) / 2;

  if (x / mid < mid)
    return ctSqrt<T>(x, lo, mid - 1);
  else
    return ctSqrt(x, mid, hi);
}

#else

#define MID ((lo + hi + 1) / 2)

/// C++11 compile time square root using binary search
template <typename T> constexpr T ctSqrt(T x, T lo, T hi) {
  return lo == hi ? lo
                  : ((x / MID < MID) ? ctSqrt<T>(x, lo, MID - 1)
                                     : ctSqrt<T>(x, MID, hi));
}

#endif

template <typename T> constexpr T ctSqrt(T x) {
  return ctSqrt<T>(x, 0, x / 2 + 1);
}

template <typename T> inline T isqrt(T x) {
  T r = (T)std::sqrt((double)x);

  constexpr T maxSqrt = ctSqrt(std::numeric_limits<T>::max());
  r = std::min(r, maxSqrt);

  while (r * r > x)
    r--;
  while (x - r * r > r * 2)
    r++;

  return r;
}

/// Returns 2^64-1 if (x + y) > 2^64-1
inline uint64_t checkedAdd(uint64_t x, uint64_t y) {
  if (x >= std::numeric_limits<uint64_t>::max() - y)
    return std::numeric_limits<uint64_t>::max();
  else
    return x + y;
}

/// Returns 0 if (x - y) < 0
inline uint64_t checkedSub(uint64_t x, uint64_t y) {
  if (x > y)
    return x - y;
  else
    return 0;
}

template <typename A, typename B, typename C>
inline B inBetween(A min, B x, C max) {
  using T = typename std::common_type<A, B, C>::type;

  if ((T)x < (T)min)
    return (B)min;
  if ((T)x > (T)max)
    return (B)max;

  return x;
}

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

/// Approximation of the maximum prime gap near n
template <typename T> inline T maxPrimeGap(T n) {
  double x = (double)n;
  x = std::max(8.0, x);
  double logx = std::log(x);
  double prime_gap = logx * logx;

  return (T)prime_gap;
}

} // namespace

#endif

///
/// @file   Wheel.hpp
/// @brief  Wheel factorization is used to skip multiles of
///         small primes in the sieve of Eratosthenes.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef WHEEL_HPP
#define WHEEL_HPP

namespace primesieve {

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

} // namespace primesieve

#endif

///
/// @file   config.hpp
/// @brief  primesieve compile time constants.
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef CONFIG_HPP
#define CONFIG_HPP

namespace {
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
} // namespace

#endif

#ifndef BUCKET_HPP
#define BUCKET_HPP

// #include "pmath.hpp"

namespace primesieve {

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
};

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

static_assert(isPow2(sizeof(Bucket)), "sizeof(Bucket) must be a power of 2!");

} // namespace primesieve

#endif

#ifndef ERATSMALL_HPP
#define ERATSMALL_HPP

// #include "macros.hpp"

namespace primesieve {

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

private:
  uint64_t maxPrime_ = 0;
  uint64_t l1CacheSize_ = 0;
  std::vector<SievingPrime> primes_;
  bool enabled_ = false;
  void storeSievingPrime(uint64_t, uint64_t, uint64_t);
  NOINLINE void crossOff(uint8_t *, uint8_t *);
};

} // namespace primesieve

#endif

///
/// @file  MemoryPool.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef MEMORYPOOL_HPP
#define MEMORYPOOL_HPP

// #include "Bucket.hpp"

namespace primesieve {

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

} // namespace primesieve

#endif

#ifndef ERATMEDIUM_HPP
#define ERATMEDIUM_HPP

// #include "Bucket.hpp"
// #include "macros.hpp"
// #include "Wheel.hpp"

namespace primesieve {

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

} // namespace primesieve

#endif

///
/// @file  EratBig.hpp
///
/// Copyright (C) 2020 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef ERATBIG_HPP
#define ERATBIG_HPP

// #include "Bucket.hpp"
// #include "macros.hpp"
// #include "MemoryPool.hpp"
// #include "Wheel.hpp"

namespace primesieve {

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

} // namespace primesieve

#endif

#ifndef ERAT_HPP
#define ERAT_HPP

// #include "forward.hpp"

/// In order convert 1 bits of the sieve array into primes we
/// need to quickly calculate the index of the first set bit.
/// This CPU instruction is usually named CTZ (Count trailing
/// zeros). Unfortunately only x86 and x64 CPUs currently have
/// an instruction for that, other CPU architectures like
/// ARM64 or PPC64 have to emulate the CTZ instruction using
/// multiple other instructions.
///
/// On x64 CPUs there are actually 2 instructions to count the
/// number of trailing zeros: BSF and TZCNT. BSF is an old
/// instruction whereas TZCNT is much more recent (Bit
/// Manipulation Instruction Set 1). Since I expect BSF to be
/// slow on future x64 CPUs (because it is legacy) we only use
/// __builtin_ctzll() if we can guarantee that TZCNT will be
/// generated.
///
/// There is also a quick, pure integer algorithm known for
/// quickly computing the index of the 1st set bit. This
/// algorithm is named the "De Bruijn bitscan".
/// https://www.chessprogramming.org/BitScan
///
/// Because of this situation, we only use __builtin_ctzll()
/// or std::countr_zero() when we know that the user's CPU
/// architecture can quickly compute CTZ, either using a single
/// instruction or emulated using very few instructions. For
/// all other CPU architectures we fallback to the "De Bruijn
/// bitscan" algorithm.

#if !defined(__has_builtin)
#define __has_builtin(x) 0
#endif

#if !defined(__has_include)
#define __has_include(x) 0
#endif

#if __cplusplus >= 202002L &&                                                  \
    __has_include(<bit>)&&(defined(__BMI__) /* TZCNT (x64) */ ||               \
                           defined(__aarch64__) /* CTZ = RBIT + CLZ */ ||      \
                           defined(_M_ARM64) /* CTZ = RBIT + CLZ */)

#include <bit>
#define ctz64(x) std::countr_zero(x)

#elif __has_builtin(__builtin_ctzll) &&                                        \
    (defined(__BMI__) /* TZCNT (x64) */ ||                                     \
     defined(__aarch64__) /* CTZ = RBIT + CLZ */ ||                            \
     defined(_M_ARM64) /* CTZ = RBIT + CLZ */)

#define ctz64(x) __builtin_ctzll(x)

#endif

namespace primesieve {

class PreSieve;

/// The abstract Erat class sieves primes using the segmented sieve
/// of Eratosthenes. It uses a bit array for sieving, the bit array
/// uses 8 flags for 30 numbers. Erat uses 3 different sieve of
/// Eratosthenes algorithms optimized for small, medium and big
/// sieving primes to cross-off multiples.
///
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

} // namespace primesieve

#endif

#ifndef PRINTPRIMES_HPP
#define PRINTPRIMES_HPP

// #include "macros.hpp"
// #include "PrimeSieve.hpp"

namespace primesieve {

class Store;

/// After a segment has been sieved PrintPrimes is
/// used to reconstruct primes and prime k-tuplets from
/// 1 bits of the sieve array
///
class PrintPrimes : public Erat {
public:
  PrintPrimes(PrimeSieve &);
  NOINLINE void sieve();

private:
  uint64_t low_ = 0;
  /// Count lookup tables for prime k-tuplets
  std::vector<uint8_t> kCounts_[6];
  counts_t &counts_;
  /// Reference to the associated PrimeSieve object
  PrimeSieve &ps_;
  void initCounts();
  void print();
  void countPrimes();
  void countkTuplets();
  void printPrimes() const;
  void printkTuplets() const;
};

} // namespace primesieve

#endif

// #include <primesieve/forward.hpp>
// #include <primesieve/PreSieve.hpp>

namespace {

struct SmallPrime {
  uint64_t first;
  uint64_t last;
  int index;
  string str;
};

const array<SmallPrime, 8> smallPrimes{{{2, 2, 0, "2"},
                                        {3, 3, 0, "3"},
                                        {5, 5, 0, "5"},
                                        {3, 5, 1, "(3, 5)"},
                                        {5, 7, 1, "(5, 7)"},
                                        {5, 11, 2, "(5, 7, 11)"},
                                        {5, 13, 3, "(5, 7, 11, 13)"},
                                        {5, 17, 4, "(5, 7, 11, 13, 17)"}}};

} // namespace

namespace primesieve {

PrimeSieve::PrimeSieve() {
  int sieveSize = get_sieve_size();
  setSieveSize(sieveSize);
}

/// Used for multi-threading
PrimeSieve::PrimeSieve(ParallelSieve *parent)
    : flags_(parent->flags_), sieveSize_(parent->sieveSize_), parent_(parent) {}

PrimeSieve::~PrimeSieve() = default;

void PrimeSieve::reset() {
  counts_.fill(0);
  percent_ = -1.0;
  seconds_ = 0.0;
  sievedDistance_ = 0;
}

bool PrimeSieve::isFlag(int flag) const { return (flags_ & flag) == flag; }

bool PrimeSieve::isFlag(int first, int last) const {
  int mask = (last * 2) - first;
  return (flags_ & mask) != 0;
}

bool PrimeSieve::isCountPrimes() const { return isFlag(COUNT_PRIMES); }

bool PrimeSieve::isPrintPrimes() const { return isFlag(PRINT_PRIMES); }

bool PrimeSieve::isPrint() const {
  return isFlag(PRINT_PRIMES, PRINT_SEXTUPLETS);
}

bool PrimeSieve::isCountkTuplets() const {
  return isFlag(COUNT_TWINS, COUNT_SEXTUPLETS);
}

bool PrimeSieve::isPrintkTuplets() const {
  return isFlag(PRINT_TWINS, PRINT_SEXTUPLETS);
}

bool PrimeSieve::isStatus() const { return isFlag(PRINT_STATUS); }

bool PrimeSieve::isCount(int i) const { return isFlag(COUNT_PRIMES << i); }

bool PrimeSieve::isPrint(int i) const { return isFlag(PRINT_PRIMES << i); }

uint64_t PrimeSieve::getStart() const { return start_; }

uint64_t PrimeSieve::getStop() const { return stop_; }

uint64_t PrimeSieve::getDistance() const {
  if (start_ <= stop_)
    return stop_ - start_;
  else
    return 0;
}

uint64_t PrimeSieve::getCount(int i) const { return counts_.at(i); }

counts_t &PrimeSieve::getCounts() { return counts_; }

int PrimeSieve::getSieveSize() const { return sieveSize_; }

double PrimeSieve::getSeconds() const { return seconds_; }

PreSieve &PrimeSieve::getPreSieve() { return preSieve_; }

void PrimeSieve::setFlags(int flags) { flags_ = flags; }

void PrimeSieve::addFlags(int flags) { flags_ |= flags; }

void PrimeSieve::setStart(uint64_t start) { start_ = start; }

void PrimeSieve::setStop(uint64_t stop) { stop_ = stop; }

/// Set the size of the sieve array in KiB (kibibyte)
void PrimeSieve::setSieveSize(int sieveSize) {
  sieveSize_ = inBetween(8, sieveSize, 4096);
  sieveSize_ = floorPow2(sieveSize_);
}

void PrimeSieve::setStatus(double percent) {
  if (!parent_) {
    auto old = percent_;
    percent_ = percent;
    if (isFlag(PRINT_STATUS))
      printStatus(old, percent_);
  }
}

void PrimeSieve::updateStatus(uint64_t dist) {
  if (parent_) {
    // This is a worker thread, so we need
    // to send the update status request
    // to the parent object which handles
    // thread synchronization.
    updateDistance_ += dist;
    if (parent_->tryUpdateStatus(updateDistance_))
      updateDistance_ = 0;
  } else {
    sievedDistance_ += dist;
    double percent = 100;
    if (getDistance() > 0)
      percent = sievedDistance_ * 100.0 / getDistance();
    auto old = percent_;
    percent_ = min(percent, 100.0);
    if (isFlag(PRINT_STATUS))
      printStatus(old, percent_);
  }
}

void PrimeSieve::printStatus(double old, double current) {
  int percent = (int)current;
  if (percent > (int)old) {
    cout << '\r' << percent << '%' << flush;
    if (percent == 100)
      cout << '\n';
  }
}

/// Process small primes <= 5 and small k-tuplets <= 17
void PrimeSieve::processSmallPrimes() {
  for (auto &p : smallPrimes) {
    if (p.first >= start_ && p.last <= stop_) {
      if (isCount(p.index))
        counts_[p.index]++;
      if (isPrint(p.index))
        cout << p.str << '\n';
    }
  }
}

uint64_t PrimeSieve::countPrimes(uint64_t start, uint64_t stop) {
  sieve(start, stop, COUNT_PRIMES);
  return getCount(0);
}

void PrimeSieve::sieve(uint64_t start, uint64_t stop) {
  setStart(start);
  setStop(stop);
  sieve();
}

void PrimeSieve::sieve(uint64_t start, uint64_t stop, int flags) {
  setStart(start);
  setStop(stop);
  setFlags(flags);
  sieve();
}

/// Sieve the primes and prime k-tuplets (twin primes,
/// prime triplets, ...) in [start, stop].
///
void PrimeSieve::sieve() {
  reset();

  if (start_ > stop_)
    return;

  setStatus(0);
  auto t1 = chrono::system_clock::now();

  if (start_ <= 5)
    processSmallPrimes();

  if (stop_ >= 7) {
    PrintPrimes printPrimes(*this);
    printPrimes.sieve();
  }

  auto t2 = chrono::system_clock::now();
  chrono::duration<double> seconds = t2 - t1;
  seconds_ = seconds.count();
  setStatus(100);
}

} // namespace primesieve

// #include <isqrt.hpp>
namespace primesum {

/// Generate a vector with the primes <= max.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
vector<int32_t> generate_primes(int64_t max) {
  vector<int32_t> primes = {0};
  primesieve::generate_primes(max, &primes);
  return primes;
}

/// Generate a vector with the first n primes.
/// The primes vector uses 1-indexing i.e. primes[1] = 2.
///
vector<int32_t> generate_n_primes(int64_t n) {
  vector<int32_t> primes = {0};
  primesieve::generate_n_primes(n, &primes);
  return primes;
}

/// Generate a vector with the prime counts <= max
/// using the sieve of Eratosthenes
///
vector<int32_t> generate_pi(int64_t max) {
  int64_t sqrt = isqrt(max);
  int64_t size = max + 1;
  vector<char> sieve(size, 1);

  for (int64_t i = 2; i <= sqrt; i++)
    if (sieve[i])
      for (int64_t j = i * i; j < size; j += i)
        sieve[j] = 0;

  vector<int32_t> pi(size, 0);
  int32_t pix = 0;

  for (int64_t i = 2; i < size; i++) {
    pix += sieve[i];
    pi[i] = pix;
  }

  return pi;
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

} // namespace primesum

// #include <primesum-internal.hpp>
// #include <imath.hpp>
// #include <min_max.hpp>
using namespace primesum;

namespace {

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

} // namespace

namespace primesum {
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
      PiTable pi(max(sqrtx, primes[a]));
      PhiCache cache(primes, pi);

      int64_t c = PhiTiny::get_c(sqrtx);
      int64_t pi_sqrtx = min(pi[sqrtx], a);
      int64_t thread_threshold = ipow(10ll, 10);
      threads = ideal_num_threads(threads, x, thread_threshold);

      sum = phi_tiny(x, c) - a + pi_sqrtx;

#pragma omp parallel for num_threads(threads)                                  \
    schedule(dynamic, 16) firstprivate(cache) reduction(+ : sum)
      for (int64_t i = c; i < pi_sqrtx; i++)
        sum += cache.phi<-1>(x / primes[i + 1], i);
    }
  }

  return sum;
}

} // namespace primesum

// #include <primesum-internal.hpp>
// #include <imath.hpp>

namespace primesum {

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

} // namespace primesum
namespace primesum {
int256_t pi_lmo1(int128_t x) {
  if (x < 2)
    return 0;

  int64_t y = iroot<3>(x);
  int64_t pi_y = pi_legendre(y, 1);
  // int64_t c = PhiTiny::get_c(y);
  // int256_t S1 = 0;
  // int256_t S2 = 0;

  // vector<int32_t> primes = generate_primes(y);
  // vector<int32_t> lpf = generate_lpf(y);
  // vector<int32_t> mu = generate_moebius(y);

  //// Calculate the contribution of the ordinary leaves
  // for (int64_t n = 1; n <= y; n++)
  //  if (lpf[n] > primes[c])
  //    S1 += phi_sum(x / n, c) * (mu[n] * n);

  //// Calculate the contribution of the special leaves
  // for (int64_t b = c + 1; b < pi_y; b++)
  //  for (int128_t m = (y / primes[b]) + 1; m <= y; m++)
  //    if (lpf[m] > primes[b])
  //      S2 -= phi_sum(x / (primes[b] * m), b - 1) * (mu[m] * m * primes[b]);

  // int256_t phi = S1 + S2;
  // int256_t p2 = P2(x, y, 1);
  // int256_t sum = phi + prime_sum_tiny(y) - 1 - p2;

  return sum;
}

} // namespace primesum
