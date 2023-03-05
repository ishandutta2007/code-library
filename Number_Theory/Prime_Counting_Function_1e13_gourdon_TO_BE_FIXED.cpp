///
/// @file  pi_gourdon.cpp
/// @brief Implementation of Xavier Gourdon's prime counting
///        algorithm. Xavier Gourdon's algorithm is an improved
///        version of the Deleglise-Rivat algorithm, according to my
///        benchmarks Gourdon's algorithm runs up to 2x faster than
///        the Deleglise-Rivat algorithm.
///
///        Xavier Gourdon formula:
///        pi(x) = A - B + C + D + Phi0 + Sigma
///
/// Copyright (C) 2021 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///
// https://github.com/kimwalisch/primecount

#include <bits/stdc++.h>
using namespace std;
using ll = long long;
typedef __int128_t i128;
typedef __int128_t maxint_t;
#define if_unlikely(x) if (x)

inline uint64_t popcnt64(uint64_t x) {
#if __cplusplus >= 201703L
  if constexpr (sizeof(int) >= sizeof(uint64_t))
    return (uint64_t)__builtin_popcount(x);
  else if constexpr (sizeof(long) >= sizeof(uint64_t))
    return (uint64_t)__builtin_popcountl(x);
  else if constexpr (sizeof(long long) >= sizeof(uint64_t))
    return (uint64_t)__builtin_popcountll(x);
#else
  return (uint64_t)__builtin_popcountll(x);
#endif
}

// Tuning factor used in Xavier Gourdon's algorithm
double alpha_y_ = -1;

// Tuning factor used in Xavier Gourdon's algorithm
double alpha_z_ = -1;

template <typename T> inline T ipow(T x, int n) {
  T r = 1;
  for (int i = 0; i < n; i++)
    r *= x;

  return r;
}

/// Integer nth root
template <int N, typename T> inline T iroot(T x) {
  T r;

  if (N == 3)
    r = (T)std::cbrt((double)x);
  else
    r = (T)std::pow((double)x, 1.0 / N);

  // fix root too large
  for (; r > 0; r--)
    if (ipow(r, N - 1) <= x / r)
      break;

  // fix root too small
  while (ipow(r + 1, N - 1) <= x / (r + 1))
    r += 1;

  return r;
}

template <typename T1, typename T2, typename T3>
inline T2 in_between(T1 min, T2 x, T3 max) {
  if (x < min || max < min)
    return (T2)min;
  if (x > max)
    return (T2)max;

  return x;
}

template <typename T> T D_approx(T x, T sigma, T phi0, T ac, T b) {
  T pix = Li(x);
  T d_approx = pix - (ac - b + phi0 + sigma);
  d_approx = std::max(d_approx, (T)0);
  return d_approx;
}

double truncate3(double n) { return (int64_t)(n * 1000) / 1000.0; }

template <typename A, typename B> inline A ceil_div(A a, B b) {
  return (A)((a + b - 1) / b);
}

template <typename T> T isqrt(T x) {
  T s = (T)std::sqrt((double)x);

  // By using constexpr for the sqrt_max variable type it
  // is guaranteed that ct_sqrt() is evaluated at compile
  // time. Compilation will fail if the compiler fails to
  // evaluate ct_sqrt() at compile time. This is great,
  // ct_sqrt() must be evaluated at compile time otherwise
  // the runtime complexity of isqrt(x) would deteriorate
  // from O(1) to O(log2(x)).
  //
  // If sqrt_max were declared without constexpr then the
  // compiler would be free to compute ct_sqrt() either at
  // compile time or at run time e.g. GCC-11 computes
  // ct_sqrt(MAX_INT128) at compile time whereas Clang-12
  // computes ct_sqrt(MAX_INT128) at run time even at -O2.
  //
  // C++20 fixed this annoying issue by adding consteval
  // to C++. Hence if the compiler supports C++20 ct_sqrt()
  // is defined as consteval instead of constexpr. Hence
  // using C++20 ct_sqrt() will be evaluated at compile
  // time in all cases i.e. even if sqrt_max were declared
  // without constexpr.
  //
  constexpr T sqrt_max = ct_sqrt(std::numeric_limits<T>::max());

  // For 128-bit integers we use uint64_t as the
  // result type. For all other types we use the
  // same result type as the input type.
  using R = typename std::conditional<sizeof(T) / 2 == sizeof(uint64_t),
                                      uint64_t, T>::type;
  R r = (R)std::min(s, sqrt_max);

  // In my tests the first corrections were needed above
  // 10^22 where the results were off by 1. Above 10^32 the
  // first results occurred that were off by > 1. Since
  // primecount only supports numbers up to 10^31 this is
  // not an issue for us.
  if (r * (T)r > x) {
    do {
      r--;
    } while (r * (T)r > x);
  }
  // Same as (r + 1)^2 < x but overflow safe
  else if ((T)(r * 2) < x - r * (T)r) {
    do {
      r++;
    } while ((T)(r * 2) < x - r * (T)r);
  }

  return r;
}

std::ostream &operator<<(std::ostream &dest, __int128_t value) {
  std::ostream::sentry s(dest);
  if (s) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[128];
    char *d = std::end(buffer);
    do {
      --d;
      *d = "0123456789"[tmp % 10];
      tmp /= 10;
    } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (dest.rdbuf()->sputn(d, len) != len) {
      dest.setstate(std::ios_base::badbit);
    }
  }
  return dest;
}

int64_t get_x_star_gourdon(maxint_t x, int64_t y) {
  // For some unknown reason it is necessary
  // to round up (x / y^2). Without rounding up
  // there are many miscalculations below 2000
  // in my implementation.
  y = max(y, (int64_t)1);
  maxint_t yy = (maxint_t)y * y;
  maxint_t x_div_yy = ceil_div(x, yy);

  int64_t x_star = (int64_t)max(iroot<4>(x), x_div_yy);
  int64_t sqrt_xy = (int64_t)isqrt(x / y);

  // x_star <= y
  // x_star <= (x / y)^(1/2)
  // The bounds above are missing in Xavier Gourdon's
  // paper. Without these bounds many of the 7 Sigma
  // formulas (Σ0 - Σ6) return incorrect results for
  // numbers below 10^6.
  x_star = min(x_star, y);
  x_star = min(x_star, sqrt_xy);
  x_star = max(x_star, (int64_t)1);

  return x_star;
}

/// Tuning factor used in Xavier Gourdon's algorithm.
double get_alpha_y(maxint_t x, int64_t y) {
  // y = x13 * alpha_y, thus alpha = y / x13
  double x13 = (double)iroot<3>(x);
  double alpha_y = (double)y / x13;
  double max_alpha_y = (double)y;
  int64_t verify_y = (int64_t)(x13 * alpha_y);

  // Prevent x^(1/3) * alpha_y = 23.99999...
  if (verify_y < y)
    alpha_y = std::nextafter(alpha_y, max_alpha_y);

  return alpha_y;
}

/// Tuning factor used in Xavier Gourdon's algorithm.
double get_alpha_z(int64_t y, int64_t z) {
  // z = y * alpha_z, thus alpha_z = z / y
  double alpha_z = (double)z / y;
  double max_alpha_z = (double)z;
  int64_t verify_z = (int64_t)(y * alpha_z);

  // Prevent y * alpha_z = 23.99999...
  if (verify_z < z)
    alpha_z = std::nextafter(alpha_z, max_alpha_z);

  return alpha_z;
}

void print_threads(int threads) {
  cout << "threads = " << threads << std::endl;
}

/// Used by pi_gourdon(x)
void print_gourdon(maxint_t x, int64_t y, int64_t z, int64_t k, int threads) {
  cout << "x = " << x << endl;
  cout << "y = " << y << endl;
  cout << "z = " << z << endl;
  cout << "k = " << k << endl;
  cout << "x_star = " << get_x_star_gourdon(x, y) << endl;
  cout << "alpha_y = " << std::fixed << std::setprecision(3)
       << get_alpha_y(x, y) << endl;
  cout << "alpha_z = " << std::fixed << std::setprecision(3)
       << get_alpha_z(y, z) << endl;
  print_threads(threads);
}

// #include <gourdon.hpp>
// #include <primecount.hpp>
// #include <primecount-internal.hpp>
// #include <imath.hpp>
// #include <PhiTiny.hpp>
// #include <print.hpp>
// #include <to_string.hpp>
/// In Xavier Gourdon's algorithm there are 2 alpha tuning
/// factors. The alpha_y tuning factor should grow like
/// O(log(x)^3) and the alpha_z tuning factor is a small
/// constant. Both alpha_y and alpha_z should be determined
/// experimentally by running benchmarks.
/// @see doc/alpha-factor-gourdon.pdf
///
/// y = x^(1/3) * alpha_y, with alpha_y >= 1.
/// z = y * alpha_z, with alpha_z >= 1.
/// alpha_y * alpha_z <= x^(1/6)
///
std::pair<double, double> get_alpha_gourdon(maxint_t x) {
  double alpha_y = alpha_y_;
  double alpha_z = alpha_z_;
  double x16 = (double)iroot<6>(x);
  double logx = std::log((double)x);
  double alpha_yz;

  // For x <= 10^11 our default formula does not
  // generate good alpha values. Hence we use
  // another formula optimized for small values.
  if (x <= 1e11) {
    double a = 0.078173;
    double b = 1;
    alpha_yz = a * logx + b;
  } else {
    double a = 0.00526934;
    double b = -0.495545;
    double c = 16.5791;
    double d = -183.836;
    double logx2 = logx * logx;
    double logx3 = logx * logx * logx;
    alpha_yz = a * logx3 + b * logx2 + c * logx + d;
  }

  // Use default alpha_z
  if (alpha_z < 1) {
    // y = x^(1/3) * alpha_y
    // z = y * alpha_z
    //
    // alpha_y should grow like O(log(x)^3) just like in the
    // Deleglise-Rivat algorithm whereas alpha_z is a small tuning
    // factor usually within [1, 4]. In my opinion the algorithm is
    // theoretically most efficient (i.e. uses the fewest number of
    // instructions) if (y == z), hence if alpha_z = 1. Because when
    // setting y to a value smaller than z this will decrease the
    // number of sparse easy leaves (which can be computed more
    // efficiently than other types of leaves) and increase the
    // number of other types of leaves.
    //
    // By setting alpha_z to a value > 1 this will cause y to be set
    // to a value < z which will generally improve the cache
    // efficiency of the algorithm but as a drawback also increase
    // the number of instructions used by the algorithm. The C1
    // algorithm (in AC.cpp) has severe scaling issues above 10^23
    // as it is not segmented and requires frequent thread
    // synchronization. The larger alpha_z, the less work there will
    // be in the C1 algorithm. Hence for computations >= 10^23 using
    // an alpha_z > 1 will likely improve performance.
    alpha_z = 2;

    // alpha_z should be significantly smaller than alpha_y
    alpha_z = in_between(1, alpha_yz / 5, alpha_z);
  }

  // Use default alpha_y
  if (alpha_y < 1)
    alpha_y = alpha_yz / alpha_z;

  // Preserve 3 digits after decimal point
  alpha_y = in_between(1, alpha_y, x16);
  alpha_y = truncate3(alpha_y);
  alpha_z = truncate3(alpha_z);

  // Ensure alpha_y * alpha_z <= x^(1/6)
  alpha_y = in_between(1, alpha_y, x16);
  double max_alpha_z = max(1.0, x16 / alpha_y);
  alpha_z = in_between(1, alpha_z, max_alpha_z);

  return std::make_pair(alpha_y, alpha_z);
}

namespace primecount {

template <typename T> class pod_vector {
public:
  static_assert(std::is_trivially_destructible<T>::value,
                "pod_vector<T> only supports types with trivial destructors!");

  using value_type = T;
  pod_vector() noexcept = default;

  pod_vector(std::size_t size) { resize(size); }

  ~pod_vector() { delete[] array_; }

  /// Free all memory, the pod_vector
  /// can be reused afterwards.
  void free() noexcept {
    delete[] array_;
    array_ = nullptr;
    end_ = nullptr;
    capacity_ = nullptr;
  }

  /// Reset the pod_vector, but do not free its
  /// memory. Same as std::vector.clear().
  void clear() noexcept { end_ = array_; }

  /// Copying is slow, we prevent it
  pod_vector(const pod_vector &) = delete;
  pod_vector &operator=(const pod_vector &) = delete;

  /// Move constructor
  pod_vector(pod_vector &&other) noexcept { swap(other); }

  /// Move assignment operator
  pod_vector &operator=(pod_vector &&other) noexcept {
    if (this != &other)
      swap(other);

    return *this;
  }

  /// Better assembly than: std::swap(vect1, vect2)
  void swap(pod_vector &other) noexcept {
    T *tmp_array = array_;
    T *tmp_end = end_;
    T *tmp_capacity = capacity_;

    array_ = other.array_;
    end_ = other.end_;
    capacity_ = other.capacity_;

    other.array_ = tmp_array;
    other.end_ = tmp_end;
    other.capacity_ = tmp_capacity;
  }

  bool empty() const noexcept { return array_ == end_; }

  T &operator[](std::size_t pos) noexcept {
    // ASSERT(pos < size());
    return array_[pos];
  }

  const T &operator[](std::size_t pos) const noexcept {
    // ASSERT(pos < size());
    return array_[pos];
  }

  T *data() noexcept { return array_; }

  const T *data() const noexcept { return array_; }

  std::size_t size() const noexcept {
    // ASSERT(end_ >= array_);
    return (std::size_t)(end_ - array_);
  }

  std::size_t capacity() const noexcept {
    // ASSERT(capacity_ >= array_);
    return (std::size_t)(capacity_ - array_);
  }

  T *begin() noexcept { return array_; }

  const T *begin() const noexcept { return array_; }

  T *end() noexcept { return end_; }

  const T *end() const noexcept { return end_; }

  T &front() noexcept {
    // ASSERT(!empty());
    return *array_;
  }

  const T &front() const noexcept {
    // ASSERT(!empty());
    return *array_;
  }

  T &back() noexcept {
    // ASSERT(!empty());
    return *(end_ - 1);
  }

  const T &back() const noexcept {
    // ASSERT(!empty());
    return *(end_ - 1);
  }

  void push_back(const T &value) {
    if_unlikely(end_ == capacity_)
        reserve_unchecked(std::max((std::size_t)1, capacity() * 2));
    *end_++ = value;
  }

  void push_back(T &&value) {
    if_unlikely(end_ == capacity_)
        reserve_unchecked(std::max((std::size_t)1, capacity() * 2));
    *end_++ = value;
  }

  template <class... Args> void emplace_back(Args &&...args) {
    if_unlikely(end_ == capacity_)
        reserve_unchecked(std::max((std::size_t)1, capacity() * 2));
    *end_++ = T(std::forward<Args>(args)...);
  }

  template <class InputIt>
  void insert(T *const pos, InputIt first, InputIt last) {
    // We only support appending to the vector
    // ASSERT(pos == end_);
    (void)pos;

    if (first < last) {
      std::size_t old_size = size();
      std::size_t new_size = old_size + (std::size_t)(last - first);
      reserve(new_size);
      end_ = array_ + new_size;
      std::copy(first, last, &array_[old_size]);
    }
  }

  void reserve(std::size_t n) {
    if (n > capacity())
      reserve_unchecked(n);
  }

  /// Resize without default initializing memory.
  /// If the pod_vector is not empty the current content
  /// will be copied into the new array.
  ///
  void resize(std::size_t n) {
    if (n > capacity())
      reserve_unchecked(n);
    else if (!std::is_trivial<T>::value && n > size()) {
      // This will only be used for classes
      // and structs with constructors.
      // ASSERT(n <= capacity());
      std::fill(end_, array_ + n, T());
    }

    end_ = array_ + n;
  }

private:
  T *array_ = nullptr;
  T *end_ = nullptr;
  T *capacity_ = nullptr;

  void reserve_unchecked(std::size_t n) {
    // ASSERT(n > capacity());
    std::size_t new_capacity = get_new_capacity<T>(n);
    std::size_t old_size = size();
    // ASSERT(new_capacity >= n);
    // ASSERT(new_capacity > old_size);

    // This default initializes memory of classes and
    // structs with constructors. But it does not default
    // initialize memory for POD types like int, long.
    T *old = array_;
    array_ = new T[new_capacity];
    end_ = array_ + old_size;
    capacity_ = array_ + new_capacity;

    if (old) {
      std::copy_n(old, old_size, array_);
      delete[] old;
    }
  }

  template <typename U>
  typename std::enable_if<std::is_trivial<U>::value, std::size_t>::type
  get_new_capacity(std::size_t size) {
    // ASSERT(size > 0);
    // GCC & Clang's std::vector grow the capacity by at least
    // 2x for every call to resize() with n > capacity(). We
    // grow by at least 1.5x as we tend to accurately calculate
    // the amount of memory we need upfront.
    std::size_t new_capacity = (std::size_t)(capacity() * 1.5);
    constexpr std::size_t min_alignment = sizeof(long) * 2;
    constexpr std::size_t min_capacity = min_alignment / sizeof(U);
    return std::max({min_capacity, size, new_capacity});
  }

  template <typename U>
  typename std::enable_if<!std::is_trivial<U>::value, std::size_t>::type
  get_new_capacity(std::size_t size) {
    // ASSERT(size > 0);
    // GCC & Clang's std::vector grow the capacity by at least
    // 2x for every call to resize() with n > capacity(). We
    // grow by at least 1.5x as we tend to accurately calculate
    // the amount of memory we need upfront.
    std::size_t new_capacity = (std::size_t)(capacity() * 1.5);
    return std::max(size, new_capacity);
  }
};

template <typename T, std::size_t N> class pod_array {
public:
  using value_type = T;
  T array_[N];

  T &operator[](std::size_t pos) noexcept {
    // ASSERT(pos < size());
    return array_[pos];
  }

  const T &operator[](std::size_t pos) const noexcept {
    // ASSERT(pos < size());
    return array_[pos];
  }

  void fill(const T &value) { std::fill_n(begin(), size(), value); }

  T *data() noexcept { return array_; }

  const T *data() const noexcept { return array_; }

  T *begin() noexcept { return array_; }

  const T *begin() const noexcept { return array_; }

  T *end() noexcept { return array_ + N; }

  const T *end() const noexcept { return array_ + N; }

  T &back() noexcept {
    // ASSERT(N > 0);
    return array_[N - 1];
  }

  const T &back() const noexcept {
    // ASSERT(N > 0);
    return array_[N - 1];
  }

  constexpr std::size_t size() const noexcept { return N; }
};
class BitSieve240 {
protected:
  static const pod_array<uint64_t, 6> pi_tiny_;
  static const pod_array<uint64_t, 240> set_bit_;
  static const pod_array<uint64_t, 240> unset_bit_;
  static const pod_array<uint64_t, 240> unset_larger_;
};

/// pi(x) for x < 6
const pod_array<uint64_t, 6> BitSieve240::pi_tiny_ = {0, 0, 1, 2, 2, 3};

/// Bitmasks needed to set a specific bit in the sieve array
const pod_array<uint64_t, 240> BitSieve240::set_bit_ = {
    0ull, 1ull << 0,  0ull, 0ull,       0ull, 0ull,       0ull, 1ull << 1,
    0ull, 0ull,       0ull, 1ull << 2,  0ull, 1ull << 3,  0ull, 0ull,
    0ull, 1ull << 4,  0ull, 1ull << 5,  0ull, 0ull,       0ull, 1ull << 6,
    0ull, 0ull,       0ull, 0ull,       0ull, 1ull << 7,  0ull, 1ull << 8,
    0ull, 0ull,       0ull, 0ull,       0ull, 1ull << 9,  0ull, 0ull,
    0ull, 1ull << 10, 0ull, 1ull << 11, 0ull, 0ull,       0ull, 1ull << 12,
    0ull, 1ull << 13, 0ull, 0ull,       0ull, 1ull << 14, 0ull, 0ull,
    0ull, 0ull,       0ull, 1ull << 15, 0ull, 1ull << 16, 0ull, 0ull,
    0ull, 0ull,       0ull, 1ull << 17, 0ull, 0ull,       0ull, 1ull << 18,
    0ull, 1ull << 19, 0ull, 0ull,       0ull, 1ull << 20, 0ull, 1ull << 21,
    0ull, 0ull,       0ull, 1ull << 22, 0ull, 0ull,       0ull, 0ull,
    0ull, 1ull << 23, 0ull, 1ull << 24, 0ull, 0ull,       0ull, 0ull,
    0ull, 1ull << 25, 0ull, 0ull,       0ull, 1ull << 26, 0ull, 1ull << 27,
    0ull, 0ull,       0ull, 1ull << 28, 0ull, 1ull << 29, 0ull, 0ull,
    0ull, 1ull << 30, 0ull, 0ull,       0ull, 0ull,       0ull, 1ull << 31,
    0ull, 1ull << 32, 0ull, 0ull,       0ull, 0ull,       0ull, 1ull << 33,
    0ull, 0ull,       0ull, 1ull << 34, 0ull, 1ull << 35, 0ull, 0ull,
    0ull, 1ull << 36, 0ull, 1ull << 37, 0ull, 0ull,       0ull, 1ull << 38,
    0ull, 0ull,       0ull, 0ull,       0ull, 1ull << 39, 0ull, 1ull << 40,
    0ull, 0ull,       0ull, 0ull,       0ull, 1ull << 41, 0ull, 0ull,
    0ull, 1ull << 42, 0ull, 1ull << 43, 0ull, 0ull,       0ull, 1ull << 44,
    0ull, 1ull << 45, 0ull, 0ull,       0ull, 1ull << 46, 0ull, 0ull,
    0ull, 0ull,       0ull, 1ull << 47, 0ull, 1ull << 48, 0ull, 0ull,
    0ull, 0ull,       0ull, 1ull << 49, 0ull, 0ull,       0ull, 1ull << 50,
    0ull, 1ull << 51, 0ull, 0ull,       0ull, 1ull << 52, 0ull, 1ull << 53,
    0ull, 0ull,       0ull, 1ull << 54, 0ull, 0ull,       0ull, 0ull,
    0ull, 1ull << 55, 0ull, 1ull << 56, 0ull, 0ull,       0ull, 0ull,
    0ull, 1ull << 57, 0ull, 0ull,       0ull, 1ull << 58, 0ull, 1ull << 59,
    0ull, 0ull,       0ull, 1ull << 60, 0ull, 1ull << 61, 0ull, 0ull,
    0ull, 1ull << 62, 0ull, 0ull,       0ull, 0ull,       0ull, 1ull << 63};

/// Bitmasks needed to unset a specific bit in the sieve array
const pod_array<uint64_t, 240> BitSieve240::unset_bit_ = {
    ~0ull, ~(1ull << 0),  ~0ull, ~0ull, ~0ull, ~0ull,
    ~0ull, ~(1ull << 1),  ~0ull, ~0ull, ~0ull, ~(1ull << 2),
    ~0ull, ~(1ull << 3),  ~0ull, ~0ull, ~0ull, ~(1ull << 4),
    ~0ull, ~(1ull << 5),  ~0ull, ~0ull, ~0ull, ~(1ull << 6),
    ~0ull, ~0ull,         ~0ull, ~0ull, ~0ull, ~(1ull << 7),
    ~0ull, ~(1ull << 8),  ~0ull, ~0ull, ~0ull, ~0ull,
    ~0ull, ~(1ull << 9),  ~0ull, ~0ull, ~0ull, ~(1ull << 10),
    ~0ull, ~(1ull << 11), ~0ull, ~0ull, ~0ull, ~(1ull << 12),
    ~0ull, ~(1ull << 13), ~0ull, ~0ull, ~0ull, ~(1ull << 14),
    ~0ull, ~0ull,         ~0ull, ~0ull, ~0ull, ~(1ull << 15),
    ~0ull, ~(1ull << 16), ~0ull, ~0ull, ~0ull, ~0ull,
    ~0ull, ~(1ull << 17), ~0ull, ~0ull, ~0ull, ~(1ull << 18),
    ~0ull, ~(1ull << 19), ~0ull, ~0ull, ~0ull, ~(1ull << 20),
    ~0ull, ~(1ull << 21), ~0ull, ~0ull, ~0ull, ~(1ull << 22),
    ~0ull, ~0ull,         ~0ull, ~0ull, ~0ull, ~(1ull << 23),
    ~0ull, ~(1ull << 24), ~0ull, ~0ull, ~0ull, ~0ull,
    ~0ull, ~(1ull << 25), ~0ull, ~0ull, ~0ull, ~(1ull << 26),
    ~0ull, ~(1ull << 27), ~0ull, ~0ull, ~0ull, ~(1ull << 28),
    ~0ull, ~(1ull << 29), ~0ull, ~0ull, ~0ull, ~(1ull << 30),
    ~0ull, ~0ull,         ~0ull, ~0ull, ~0ull, ~(1ull << 31),
    ~0ull, ~(1ull << 32), ~0ull, ~0ull, ~0ull, ~0ull,
    ~0ull, ~(1ull << 33), ~0ull, ~0ull, ~0ull, ~(1ull << 34),
    ~0ull, ~(1ull << 35), ~0ull, ~0ull, ~0ull, ~(1ull << 36),
    ~0ull, ~(1ull << 37), ~0ull, ~0ull, ~0ull, ~(1ull << 38),
    ~0ull, ~0ull,         ~0ull, ~0ull, ~0ull, ~(1ull << 39),
    ~0ull, ~(1ull << 40), ~0ull, ~0ull, ~0ull, ~0ull,
    ~0ull, ~(1ull << 41), ~0ull, ~0ull, ~0ull, ~(1ull << 42),
    ~0ull, ~(1ull << 43), ~0ull, ~0ull, ~0ull, ~(1ull << 44),
    ~0ull, ~(1ull << 45), ~0ull, ~0ull, ~0ull, ~(1ull << 46),
    ~0ull, ~0ull,         ~0ull, ~0ull, ~0ull, ~(1ull << 47),
    ~0ull, ~(1ull << 48), ~0ull, ~0ull, ~0ull, ~0ull,
    ~0ull, ~(1ull << 49), ~0ull, ~0ull, ~0ull, ~(1ull << 50),
    ~0ull, ~(1ull << 51), ~0ull, ~0ull, ~0ull, ~(1ull << 52),
    ~0ull, ~(1ull << 53), ~0ull, ~0ull, ~0ull, ~(1ull << 54),
    ~0ull, ~0ull,         ~0ull, ~0ull, ~0ull, ~(1ull << 55),
    ~0ull, ~(1ull << 56), ~0ull, ~0ull, ~0ull, ~0ull,
    ~0ull, ~(1ull << 57), ~0ull, ~0ull, ~0ull, ~(1ull << 58),
    ~0ull, ~(1ull << 59), ~0ull, ~0ull, ~0ull, ~(1ull << 60),
    ~0ull, ~(1ull << 61), ~0ull, ~0ull, ~0ull, ~(1ull << 62),
    ~0ull, ~0ull,         ~0ull, ~0ull, ~0ull, ~(1ull << 63)};

constexpr int right_shift(int n) {
  return (n % 30 >= 29)   ? 56 - (n / 30 * 8)
         : (n % 30 >= 23) ? 57 - (n / 30 * 8)
         : (n % 30 >= 19) ? 58 - (n / 30 * 8)
         : (n % 30 >= 17) ? 59 - (n / 30 * 8)
         : (n % 30 >= 13) ? 60 - (n / 30 * 8)
         : (n % 30 >= 11) ? 61 - (n / 30 * 8)
         : (n % 30 >= 7)  ? 62 - (n / 30 * 8)
         : (n % 30 >= 1)  ? 63 - (n / 30 * 8)
                          : 64 - (n / 30 * 8);
}

/// Bitmask to unset bits >= n
constexpr uint64_t unset_l(int n) {
  return (n == 0) ? 0 : ~0ull >> right_shift(n);
}

/// Unset bits > stop
const pod_array<uint64_t, 240> BitSieve240::unset_larger_ = {
    unset_l(0),   unset_l(1),   unset_l(2),   unset_l(3),   unset_l(4),
    unset_l(5),   unset_l(6),   unset_l(7),   unset_l(8),   unset_l(9),
    unset_l(10),  unset_l(11),  unset_l(12),  unset_l(13),  unset_l(14),
    unset_l(15),  unset_l(16),  unset_l(17),  unset_l(18),  unset_l(19),
    unset_l(20),  unset_l(21),  unset_l(22),  unset_l(23),  unset_l(24),
    unset_l(25),  unset_l(26),  unset_l(27),  unset_l(28),  unset_l(29),
    unset_l(30),  unset_l(31),  unset_l(32),  unset_l(33),  unset_l(34),
    unset_l(35),  unset_l(36),  unset_l(37),  unset_l(38),  unset_l(39),
    unset_l(40),  unset_l(41),  unset_l(42),  unset_l(43),  unset_l(44),
    unset_l(45),  unset_l(46),  unset_l(47),  unset_l(48),  unset_l(49),
    unset_l(50),  unset_l(51),  unset_l(52),  unset_l(53),  unset_l(54),
    unset_l(55),  unset_l(56),  unset_l(57),  unset_l(58),  unset_l(59),
    unset_l(60),  unset_l(61),  unset_l(62),  unset_l(63),  unset_l(64),
    unset_l(65),  unset_l(66),  unset_l(67),  unset_l(68),  unset_l(69),
    unset_l(70),  unset_l(71),  unset_l(72),  unset_l(73),  unset_l(74),
    unset_l(75),  unset_l(76),  unset_l(77),  unset_l(78),  unset_l(79),
    unset_l(80),  unset_l(81),  unset_l(82),  unset_l(83),  unset_l(84),
    unset_l(85),  unset_l(86),  unset_l(87),  unset_l(88),  unset_l(89),
    unset_l(90),  unset_l(91),  unset_l(92),  unset_l(93),  unset_l(94),
    unset_l(95),  unset_l(96),  unset_l(97),  unset_l(98),  unset_l(99),
    unset_l(100), unset_l(101), unset_l(102), unset_l(103), unset_l(104),
    unset_l(105), unset_l(106), unset_l(107), unset_l(108), unset_l(109),
    unset_l(110), unset_l(111), unset_l(112), unset_l(113), unset_l(114),
    unset_l(115), unset_l(116), unset_l(117), unset_l(118), unset_l(119),
    unset_l(120), unset_l(121), unset_l(122), unset_l(123), unset_l(124),
    unset_l(125), unset_l(126), unset_l(127), unset_l(128), unset_l(129),
    unset_l(130), unset_l(131), unset_l(132), unset_l(133), unset_l(134),
    unset_l(135), unset_l(136), unset_l(137), unset_l(138), unset_l(139),
    unset_l(140), unset_l(141), unset_l(142), unset_l(143), unset_l(144),
    unset_l(145), unset_l(146), unset_l(147), unset_l(148), unset_l(149),
    unset_l(150), unset_l(151), unset_l(152), unset_l(153), unset_l(154),
    unset_l(155), unset_l(156), unset_l(157), unset_l(158), unset_l(159),
    unset_l(160), unset_l(161), unset_l(162), unset_l(163), unset_l(164),
    unset_l(165), unset_l(166), unset_l(167), unset_l(168), unset_l(169),
    unset_l(170), unset_l(171), unset_l(172), unset_l(173), unset_l(174),
    unset_l(175), unset_l(176), unset_l(177), unset_l(178), unset_l(179),
    unset_l(180), unset_l(181), unset_l(182), unset_l(183), unset_l(184),
    unset_l(185), unset_l(186), unset_l(187), unset_l(188), unset_l(189),
    unset_l(190), unset_l(191), unset_l(192), unset_l(193), unset_l(194),
    unset_l(195), unset_l(196), unset_l(197), unset_l(198), unset_l(199),
    unset_l(200), unset_l(201), unset_l(202), unset_l(203), unset_l(204),
    unset_l(205), unset_l(206), unset_l(207), unset_l(208), unset_l(209),
    unset_l(210), unset_l(211), unset_l(212), unset_l(213), unset_l(214),
    unset_l(215), unset_l(216), unset_l(217), unset_l(218), unset_l(219),
    unset_l(220), unset_l(221), unset_l(222), unset_l(223), unset_l(224),
    unset_l(225), unset_l(226), unset_l(227), unset_l(228), unset_l(229),
    unset_l(230), unset_l(231), unset_l(232), unset_l(233), unset_l(234),
    unset_l(235), unset_l(236), unset_l(237), unset_l(238), unset_l(239)};

class PhiTiny : public BitSieve240 {
public:
  PhiTiny();

  /// Uses at most one level of phi(x, a) recursion
  /// to ensure that the runtime is O(1).
  template <typename T> T phi_recursive(T x, uint64_t a) const {
    // Unsigned integer division is usually
    // faster than signed integer division,
    // especially for int128_t.
    using UT = typename std::make_unsigned<T>::type;

    if (a < max_a())
      return phi((UT)x, a);
    else {
      // ASSERT(a == 8);
      // This code path will be executed most of the time.
      // In phi7(x) the variable a has been hardcoded to 7
      // which makes it run slightly faster than phi(x, a).
      // phi(x, 8) = phi(x, 7) - phi(x / prime[8], 7)
      return phi7((UT)x) - phi7((UT)x / 19);
    }
  }

  template <typename T> T phi(T x, uint64_t a) const {
    auto pp = prime_products[a];
    auto remainder = (uint64_t)(x % pp);
    T xpp = x / pp;
    T sum = xpp * totients[a];

    // For prime[a] <= 5 our phi(x % pp, a) lookup table
    // is a simple two dimensional array.
    if (a < phi_.size())
      sum += phi_[a][remainder];
    else {
      // For prime[a] > 5 we use a compressed phi(x % pp, a)
      // lookup table. Each bit of the sieve array corresponds
      // to an integer that is not divisible by 2, 3 and 5.
      // Hence the 8 bits of each byte correspond to the offsets
      // [ 1, 7, 11, 13, 17, 19, 23, 29 ].
      uint64_t count = sieve_[a][remainder / 240].count;
      uint64_t bits = sieve_[a][remainder / 240].bits;
      uint64_t bitmask = unset_larger_[remainder % 240];
      sum += (T)(count + popcnt64(bits & bitmask));
    }

    return sum;
  }

  /// In phi7(x) the variable a has been hardcoded to 7.
  /// phi7(x) uses division by a constant instead of regular
  /// integer division and hence phi7(x) is expected to run
  /// faster than the phi(x, a) implementation above.
  ///
  template <typename T> T phi7(T x) const {
    constexpr uint32_t a = 7;
    constexpr uint32_t pp = 510510;
    constexpr uint32_t totient = 92160;
    auto remainder = (uint64_t)(x % pp);
    T xpp = x / pp;
    T sum = xpp * totient;

    // For prime[a] > 5 we use a compressed phi(x % pp, a)
    // lookup table. Each bit of the sieve array corresponds
    // to an integer that is not divisible by 2, 3 and 5.
    // Hence the 8 bits of each byte correspond to the offsets
    // [ 1, 7, 11, 13, 17, 19, 23, 29 ].
    // ASSERT(sieve_.size() - 1 == a);
    uint64_t count = sieve_[a][remainder / 240].count;
    uint64_t bits = sieve_[a][remainder / 240].bits;
    uint64_t bitmask = unset_larger_[remainder % 240];
    sum += (T)(count + popcnt64(bits & bitmask));

    return sum;
  }

  static uint64_t get_c(uint64_t y) {
    if (y < pi.size())
      return pi[y];
    else
      return max_a();
  }

  /// In Xavier Gourdon's algorithm the small
  /// constant is named k instead of c.
  /// k <= PrimePi[min(x_star, sqrt(x / y))]
  ///
  template <typename T> static uint64_t get_k(T x) {
    return get_c(iroot<4>(x));
  }

  static constexpr uint64_t max_a() { return primes.size(); }

private:
  static const pod_array<uint32_t, 8> primes;
  static const pod_array<uint32_t, 8> prime_products;
  static const pod_array<uint32_t, 8> totients;
  static const pod_array<uint8_t, 20> pi;

/// Packing sieve_t increases the cache's capacity by 25%
/// which improves performance by up to 10%.
#pragma pack(push, 1)
  struct sieve_t {
    uint32_t count;
    uint64_t bits;
  };

#pragma pack(pop)

  /// sieve[a] contains only numbers that are not divisible
  /// by any of the the first a primes. sieve[a][i].count
  /// contains the count of numbers < i * 240 that are not
  /// divisible by any of the first a primes.
  pod_array<pod_vector<sieve_t>, 8> sieve_;
  pod_array<pod_vector<uint8_t>, 4> phi_;
};

const pod_array<uint32_t, 8> PhiTiny::primes = {0, 2, 3, 5, 7, 11, 13, 17};

// prime_products[n] = \prod_{i=1}^{n} primes[i]
const pod_array<uint32_t, 8> PhiTiny::prime_products = {
    1, 2, 6, 30, 210, 2310, 30030, 510510};

// totients[n] = \prod_{i=1}^{n} (primes[i] - 1)
const pod_array<uint32_t, 8> PhiTiny::totients = {1,  1,   2,    8,
                                                  48, 480, 5760, 92160};

// Number of primes <= next_prime(primes.back())
const pod_array<uint8_t, 20> PhiTiny::pi = {0, 0, 1, 2, 2, 3, 3, 4, 4, 4,
                                            4, 5, 5, 6, 6, 6, 6, 7, 7, 8};

// Singleton
const PhiTiny phiTiny;

PhiTiny::PhiTiny() {
  // The pi[x] lookup table must contain the number
  // of primes <= next_prime(primes.back()).
  // ASSERT(pi.back() == primes.size());
  // ASSERT(phi_.size() - 1 == (uint64_t)pi[5]);
  // ASSERT(sieve_.size() == primes.size());
  static_assert(prime_products.size() == primes.size(),
                "Invalid prime_products size!");
  static_assert(totients.size() == primes.size(), "Invalid totients size!");

  // a = 0
  phi_[0].resize(1);
  phi_[0][0] = 0;

  for (uint64_t a = 1; a < sieve_.size(); a++) {
    // For prime[a] <= 5 our phi(x % pp, a) lookup table
    // is a simple two dimensional array.
    if (a < phi_.size()) {
      uint64_t pp = prime_products[a];
      phi_[a].resize(pp);
      phi_[a][0] = 0;

      for (uint64_t x = 1; x < pp; x++) {
        uint64_t phi_xa = phi(x, a - 1) - phi(x / primes[a], a - 1);
        // ASSERT(phi_xa <= std::numeric_limits<uint8_t>::max());
        phi_[a][x] = (uint8_t)phi_xa;
      }
    } else {
      // For prime[a] > 5 we use a compressed phi(x % pp, a)
      // lookup table. Each bit of the sieve array corresponds
      // to an integer that is not divisible by 2, 3 and 5.
      // Hence the 8 bits of each byte correspond to the offsets
      // [ 1, 7, 11, 13, 17, 19, 23, 29 ].
      uint64_t pp = prime_products[a];
      uint64_t size = ceil_div(pp, 240);
      sieve_[a].resize(size);
      std::fill_n(sieve_[a].begin(), size, sieve_t{0, ~0ull});

      for (uint64_t i = pi[7]; i <= a; i++)
        for (uint64_t n = primes[i]; n < pp; n += primes[i] * 2)
          sieve_[a][n / 240].bits &= unset_bit_[n % 240];

      // Fill an array with the cumulative 1 bit counts.
      // sieve[i][j] contains the count of numbers < j * 240 that
      // are not divisible by any of the first i primes.
      uint64_t count = 0;
      for (auto &sieve : sieve_[a]) {
        sieve.count = (uint32_t)count;
        count += popcnt64(sieve.bits);
      }
    }
  }
}

extern const PhiTiny phiTiny;

inline bool is_phi_tiny(uint64_t a) { return a <= PhiTiny::max_a(); }

// template <typename T>
// typename std::enable_if<(sizeof(T) == sizeof(typename
// make_smaller<T>::type)),
//                         T>::type
// phi_tiny(T x, uint64_t a) {
//   return phiTiny.phi_recursive(x, a);
// }

// template <typename T>
// typename std::enable_if<(sizeof(T) > sizeof(typename make_smaller<T>::type)),
//                         T>::type
// phi_tiny(T x, uint64_t a) {
//   using smaller_t = typename make_smaller<T>::type;

//   // If possible use smaller integer type
//   // to speed up integer division.
//   if (x <= std::numeric_limits<smaller_t>::max())
//     return phiTiny.phi_recursive((smaller_t)x, a);
//   else
//     return phiTiny.phi_recursive(x, a);
// }
bool print_ = false;
bool print_variables_ = false;

bool is_print_variables() { return print_variables_; }
/// Only enabled for partial formulas
void print_gourdon_vars(maxint_t x, int64_t y, int threads) {
  if (is_print_variables()) {
    std::cout << "x = " << x << std::endl;
    std::cout << "y = " << y << std::endl;
    std::cout << "alpha_y = " << std::fixed << std::setprecision(3)
              << get_alpha_y(x, y) << std::endl;
    print_threads(threads);
    std::cout << std::endl;
  }
}

template <typename A, typename B, typename C> inline A max3(A a, B b, C c) {
  // static_assert(is_comparable_3<A, B, C>::value,
  //               "max3(A, B, C): Cannot compare types A, B and C");

  return std::max(a, (A)std::max(b, (B)c));
}

class PiTable : public BitSieve240 {
public:
  PiTable(uint64_t max_x, int threads);

  uint64_t size() const { return max_x_ + 1; }

  static int64_t max_cached() { return pi_cache_.size() * 240 - 1; }

  /// Get number of primes <= x
  int64_t operator[](uint64_t x) const {
    // ASSERT(x <= max_x_);

    if_unlikely(x < pi_tiny_.size()) return pi_tiny_[x];

    uint64_t count = pi_[x / 240].count;
    uint64_t bits = pi_[x / 240].bits;
    uint64_t bitmask = unset_larger_[x % 240];
    return count + popcnt64(bits & bitmask);
  }

  /// Get number of primes <= x
  static int64_t pi_cache(uint64_t x) {
    if_unlikely(x < pi_tiny_.size()) return pi_tiny_[x];

    uint64_t count = pi_cache_[x / 240].count;
    uint64_t bits = pi_cache_[x / 240].bits;
    uint64_t bitmask = unset_larger_[x % 240];
    return count + popcnt64(bits & bitmask);
  }

private:
  struct pi_t {
    uint64_t count;
    uint64_t bits;
  };

  void init(uint64_t limit, uint64_t cache_limit, int threads);
  void init_bits(uint64_t low, uint64_t high, uint64_t thread_num);
  void init_count(uint64_t low, uint64_t high, uint64_t thread_num);
  static const pod_array<pi_t, 64> pi_cache_;
  pod_vector<pi_t> pi_;
  pod_vector<uint64_t> counts_;
  uint64_t max_x_;
};

template <typename T> T Sigma0(T x, T a, int threads) {
  T pi_sqrtx = pi_noprint(isqrt(x), threads);
  return a - 1 + (pi_sqrtx * (pi_sqrtx - 1)) / 2 - (a * (a - 1)) / 2;
}

template <typename T> T Sigma1(T a, T b) { return (a - b) * (a - b - 1) / 2; }

template <typename T> T Sigma2(T a, T b, T c, T d) {
  return a * (b - c - (c * (c - 3)) / 2 + (d * (d - 3)) / 2);
}

template <typename T> T Sigma3(T b, T d) {
  return (b * (b - 1) * (2 * b - 1)) / 6 - b - (d * (d - 1) * (2 * d - 1)) / 6 +
         d;
}

/// Memory usage: O(x^(3/8))
template <typename T>
T Sigma456(T x, int64_t y, int64_t a, int64_t x_star, const PiTable &pi) {
  T sigma4 = 0;
  T sigma5 = 0;
  T sigma6 = 0;

  int64_t x13 = iroot<3>(x);
  int64_t sqrt_xy = isqrt(x / y);
  primesieve::iterator it(x_star, x13);
  int64_t prime = it.next_prime();

  // Sigma4: x_star < prime <= sqrt(x / y)
  // Sigma5: sqrt(x / y) < prime <= x^(1/3)
  // Sigma6: x_star < prime <= x^(1/3)
  for (; prime <= x13; prime = it.next_prime()) {
    if (prime <= sqrt_xy)
      sigma4 += pi[x / (prime * (T)y)];
    else
      sigma5 += pi[x / (prime * (T)prime)];

    // Note that in Xavier Gourdon's paper the actual
    // formula for Σ6 is: sum += pi(x^(1/2) / prime^(1/2))^2.
    // However when implemented this way using integers
    // the formula returns incorrect results.
    // Hence the formula must be implemented as below:
    int64_t sqrt_xp = isqrt(x / prime);
    int64_t pi_sqrt_xp = pi[sqrt_xp];
    sigma6 += pi_sqrt_xp * (T)pi_sqrt_xp;
  }

  sigma4 *= a;
  sigma6 = -sigma6;

  return sigma4 + sigma5 + sigma6;
}

int64_t Sigma(int64_t x, int64_t y, int threads, bool is_print) {
  if (is_print) {
    cout << " " << endl;
    cout << "=== Sigma(x, y) ===" << endl;
    print_gourdon_vars(x, y, threads);
  }

  // double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_pix_sigma4 = x / (x_star * y);
  int64_t max_pix_sigma5 = y;
  int64_t max_pix_sigma6 = isqrt(x / x_star);
  int64_t max_pix = max3(max_pix_sigma4, max_pix_sigma5, max_pix_sigma6);
  PiTable pi(max_pix, threads);

  int64_t a = pi[y];
  int64_t b = pi[iroot<3>(x)];
  int64_t c = pi[isqrt(x / y)];
  int64_t d = pi[x_star];

  int64_t sum = Sigma0(x, a, threads) + Sigma1(a, b) + Sigma2(a, b, c, d) +
                Sigma3(b, d) + Sigma456(x, y, a, x_star, pi);

  if (is_print)
    // print("Sigma", sum, time);
    cout << "Sigma" << sum << endl;

  return sum;
}

int64_t Phi0(int64_t x, int64_t y, int64_t z, int64_t k, int threads,
             bool is_print) {
  if (is_print) {
    cout << " " << endl;
    cout << "=== Phi0(x, y) ===" << endl;
    print_gourdon_vars(x, y, z, k, threads);
  }

  // double time = get_time();
  int64_t phi0 = Phi0_OpenMP(x, y, z, k, threads);

  if (is_print)
    // print("Phi0", phi0, time);
    cout < "Phi0" << phi0 << endl;

  return phi0;
}

int64_t AC(int64_t x, int64_t y, int64_t z, int64_t k, int threads,
           bool is_print) {
  if (is_print) {
    cout << " " << endl;
    cout << "=== AC(x, y) ===" << endl;
    print_gourdon_vars(x, y, z, k, threads);
  }

  // double time = get_time();
  int64_t x_star = get_x_star_gourdon(x, y);
  int64_t max_c_prime = y;
  int64_t max_a_prime = (int64_t)isqrt(x / x_star);
  int64_t max_prime = max(max_a_prime, max_c_prime);
  auto primes = generate_primes<uint32_t>(max_prime);

  int64_t sum = AC_OpenMP((uint64_t)x, y, z, k, x_star, max_a_prime, primes,
                          threads, is_print);

  if (is_print)
    // print("A + C", sum, time);
    cout << "A + C" << sum << endl;

  return sum;
}

/// Calculate the number of primes below x using
/// Xavier Gourdon's algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int64_t pi_gourdon_64(int64_t x, int threads, bool is_print) {
  if (x < 2)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t)1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t)1);

  if (is_print) {
    cout << "" << endl;
    cout << "=== pi_gourdon_64(x) ===" << endl;
    cout << "pi(x) = A - B + C + D + Phi0 + Sigma" << endl;
    print_gourdon(x, y, z, k, threads);
  }

  int64_t sigma = Sigma(x, y, threads, is_print);
  int64_t phi0 = Phi0(x, y, z, k, threads, is_print);
  int64_t b = B(x, y, threads, is_print);
  int64_t ac = AC(x, y, z, k, threads, is_print);
  int64_t d_approx = D_approx(x, sigma, phi0, ac, b);
  int64_t d = D(x, y, z, k, d_approx, threads, is_print);
  int64_t sum = ac - b + d + phi0 + sigma;

  return sum;
}

#if defined(HAVE_INT128_T)

/// Calculate the number of primes below x using
/// Xavier Gourdon's algorithm.
/// Run time: O(x^(2/3) / (log x)^2)
/// Memory usage: O(x^(1/3) * (log x)^3)
///
int128_t pi_gourdon_128(int128_t x, int threads, bool is_print) {
  if (x < 2)
    return 0;

  auto alpha = get_alpha_gourdon(x);
  double alpha_y = alpha.first;
  double alpha_z = alpha.second;
  maxint_t limit = get_max_x(alpha_y);

  if (x > limit)
    throw primecount_error("pi(x): x must be <= " + to_string(limit));

  int64_t x13 = iroot<3>(x);
  int64_t sqrtx = isqrt(x);
  int64_t y = (int64_t)(x13 * alpha_y);

  // x^(1/3) < y < x^(1/2)
  y = std::max(y, x13 + 1);
  y = std::min(y, sqrtx - 1);
  y = std::max(y, (int64_t)1);

  int64_t k = PhiTiny::get_k(x);
  int64_t z = (int64_t)(y * alpha_z);

  // y <= z < x^(1/2)
  z = std::max(z, y);
  z = std::min(z, sqrtx - 1);
  z = std::max(z, (int64_t)1);

  if (is_print) {
    cout << "" << endl;
    cout << "=== pi_gourdon_128(x) ===" << endl;
    cout << "pi(x) = A - B + C + D + Phi0 + Sigma" << endl;
    print_gourdon(x, y, z, k, threads);
  }

  int128_t sigma = Sigma(x, y, threads, is_print);
  int128_t phi0 = Phi0(x, y, z, k, threads, is_print);
  int128_t b = B(x, y, threads, is_print);
  int128_t ac = AC(x, y, z, k, threads, is_print);
  int128_t d_approx = D_approx(x, sigma, phi0, ac, b);
  int128_t d = D(x, y, z, k, d_approx, threads, is_print);
  int128_t sum = ac - b + d + phi0 + sigma;

  return sum;
}

#endif

} // namespace primecount

int main() {
  //  freopen("input.txt", "r", stdin);
  //  freopen("output.txt", "w", stdout);

  // ll n;
  // scanf("%lld", &n);
  // n = 1234567891234;
  // printf("%lld\n", prime_pi(n));
  // n = 123456789123;
  // printf("%lld\n", prime_pi(n));
  // n = 12345678912;
  // printf("%lld\n", prime_pi(n));
  // n = 1234567891;
  // printf("%lld\n", prime_pi(n));

  vector<ll> ns = {(ll)1e10, (ll)1e11, (ll)1e12, (ll)1e13};
  for (ll n : ns) {
    auto start_time = clock();
    cout << "n = " << (double)n << endl;
    ll res = primecount::pi_gourdon_64(n);
    cout << "result: " << res << endl;
    cout << "time: " << (double)(clock() - start_time) / CLOCKS_PER_SEC << "s"
         << endl;
    cout << endl;
  }

  return 0;
}
