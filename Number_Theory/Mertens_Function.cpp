// #pragma once

// #include "common/base.h"

#include <cassert>
#include <cstddef>
#include <cstdint>

// #include "common/factorization/table/mertens.h"

// #include "common/base.h"
// #include "common/factorization/table/mobius.h"

// #pragma once

// #include "common/base.h"
// #include "common/factorization/table/primes.h"
// #pragma once

// #include "common/base.h"

#include <vector>

namespace factorization {
namespace table {
class Primes {
protected:
  uint64_t table_size, squared_table_size;
  std::vector<uint64_t> primes, squared_primes;
  std::vector<unsigned> table;

public:
  explicit Primes(uint64_t size) {
    table_size = size;
    squared_table_size = table_size * table_size;
    table.resize(table_size + 1, 0);
    table[0] = table[1] = 1;
    primes.push_back(2);
    for (uint64_t i = 2; i <= table_size; i += 2)
      table[i] = 2;
    for (uint64_t i = 3; i <= table_size; i += 2) {
      if (!table[i]) {
        primes.push_back(i);
        table[i] = unsigned(i);
        for (uint64_t j = i * i; j <= table_size; j += 2 * i) {
          if (table[j] == 0)
            table[j] = unsigned(i);
        }
      }
    }
    squared_primes.reserve(primes.size());
    for (uint64_t p : primes)
      squared_primes.push_back(p * p);
  }

  const std::vector<uint64_t> &GetPrimes() const { return primes; }

  const std::vector<uint64_t> &GetSquaredPrimes() const {
    return squared_primes;
  }

  uint64_t GetTableSize() const { return table_size; }
  uint64_t GetSquaredTableSize() const { return squared_table_size; }
  const std::vector<unsigned> &GetTable() const { return table; }

  unsigned Get(uint64_t n) const { return table[n]; }
  unsigned operator()(uint64_t n) const { return Get(n); }

  bool IsPrime(uint64_t n) const { return (n > 1) && (table[n] == n); }
};
} // namespace table
} // namespace factorization

using PrimesTable = factorization::table::Primes;

#include <vector>

namespace factorization {
namespace table {
class Mobius : public Primes {
protected:
  std::vector<int> mobius;

public:
  explicit Mobius(uint64_t size) : Primes(size) {
    mobius.resize(Primes::table_size + 1);
    mobius[0] = 0;
    mobius[1] = 1;
    for (uint64_t i = 2; i <= Primes::table_size; ++i) {
      unsigned p = Primes::table[i];
      mobius[i] = (Primes::table[i / p] == p ? 0 : -1 * mobius[i / p]);
    }
  }

  int Get(uint64_t n) const { return mobius[n]; }
  int operator()(uint64_t n) const { return Get(n); }
};
} // namespace table
} // namespace factorization

#include <vector>

namespace factorization {
namespace table {
class Mertens : public Mobius {
public:
  using TBase = Mobius;

protected:
  std::vector<int> mertens;

public:
  explicit Mertens(uint64_t size) : TBase(size) {
    mertens.resize(TBase::table_size + 1);
    mertens[0] = 0;
    for (uint64_t i = 1; i <= TBase::table_size; ++i)
      mertens[i] = mertens[i - 1] + TBase::mobius[i];
  }

  int GetMobius(uint64_t n) const { return TBase::Get(n); }

  int Get(uint64_t n) const { return mertens[n]; }
  int operator()(uint64_t n) const { return Get(n); }
};
} // namespace table
} // namespace factorization

// #include "common/numeric/utils/usqrt.h"
// #pragma once

#include <cmath>
#include <limits>

namespace numeric {
namespace hidden {
template <class T> constexpr T CTUSqrtHelper(T x, T l, T h) {
  if (l == h)
    return l;
  T m = (l + h + 1) / 2;
  if (x / m < m)
    return CTUSqrtHelper<T>(x, l, m - 1);
  else
    return CTUSqrtHelper<T>(x, m, h);
}
} // namespace hidden
} // namespace numeric

template <class T> constexpr T CTUSqrt(T x) {
  return numeric::hidden::CTUSqrtHelper<T>(x, 0, x / 4 + 1);
}

template <class T> inline T USqrt(T x) {
  T r = T(sqrt(double(x)));
  T sqrt_max = CTUSqrt<T>(std::numeric_limits<T>::max());
  for (; (r < sqrt_max) && (r * r < x);)
    ++r;
  for (; (r > sqrt_max) || (r * r > x);)
    --r;
  return r;
}

// Memory: O(U)
// Time Build: O(U * log(log(U)))
// Time Get: O(X / (U^1/2))
// Optimal U ~ X^(2/3)
namespace factorization {
class Mertens {
protected:
  uint64_t u;
  table::Mertens mertens;

  // y <= U^2
  // Time: O(sqrt(y))
  int64_t S(uint64_t y) const {
    uint64_t v = USqrt(y), k = y / (v + 1);
    int64_t r = 1 + int64_t(k) * mertens(v);
    for (uint64_t n = y / u + 1; n <= k; ++n)
      r -= mertens(y / n);
    for (uint64_t n = 1; n <= v; ++n)
      r -= int64_t(y / n) * mertens.GetMobius(n);
    return r;
  }

public:
  explicit Mertens(uint64_t _u) : u(_u), mertens(u) {}

  // x <= U
  int GetMobius(uint64_t x) const { return mertens.GetMobius(x); }

  // x <= U^2
  int64_t GetMertens(uint64_t x) const {
    if (x <= u)
      return mertens(x);
    int64_t r = 0;
    for (uint64_t n = 1; n <= x / u; ++n) {
      int m = GetMobius(n);
      if (m)
        r += m * S(x / n);
    }
    return r;
  }

  int64_t operator()(uint64_t x) const { return GetMertens(x); }
};
} // namespace factorization

// #include "common/factorization/mertens.h"
// #include "common/vector/write.h"

#include <iostream>
#include <vector>

int main() {
  factorization::Mertens mertens(1ull << 20);
  // std::vector<int> expected{1,   0,     -1,  -2,   -1,    -4,   -1,    -2,
  //                           -1,  -4,    -4,  7,    -19,   22,   -32,   26,
  //                           14,  -20,   24,  -125, 257,   -362, 228,   -10,
  //                           211, -1042, 329, 330,  -1703, 6222, -10374};
  std::vector<int> output;
  for (unsigned i = 0; i <= 30; ++i) {
    // output.push_back(
    std::cout << mertens(1ull << i) << " ";
    // );
  }

  // if (output != expected) {
  //   std::cout << "Expected:" << std::endl;
  //   nvector::Write(expected);
  //   std::cout << "Output:" << std::endl;
  //   nvector::Write(output);
  //   // return false;
  // }
  // return true;
}

// https://github.com/Loks-/competitions/blob/83e4b2c57258821c78c96fbb41964ff82b1c5d1b/common/factorization/mertens.h
