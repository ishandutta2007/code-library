#include <bits/stdc++.h>

using namespace std;

#ifdef LOCAL
#include "algo/debug.h"
#else
#define debug(...) 42
#endif

template <typename T> T inverse(T a, T m) {
  T u = 0, v = 1;
  while (a != 0) {
    T t = m / a;
    m -= t * a;
    swap(a, m);
    u -= t * v;
    swap(u, v);
  }
  assert(m == 1);
  return u;
}

template <typename T> class Modular {
public:
  using Type = typename decay<decltype(T::value)>::type;

  constexpr Modular() : value() {}
  template <typename U> Modular(const U &x) { value = normalize(x); }

  template <typename U> static Type normalize(const U &x) {
    Type v;
    if (-mod() <= x && x < mod())
      v = static_cast<Type>(x);
    else
      v = static_cast<Type>(x % mod());
    if (v < 0)
      v += mod();
    return v;
  }

  const Type &operator()() const { return value; }
  template <typename U> explicit operator U() const {
    return static_cast<U>(value);
  }
  constexpr static Type mod() { return T::value; }

  Modular &operator+=(const Modular &other) {
    if ((value += other.value) >= mod())
      value -= mod();
    return *this;
  }
  Modular &operator-=(const Modular &other) {
    if ((value -= other.value) < 0)
      value += mod();
    return *this;
  }
  template <typename U> Modular &operator+=(const U &other) {
    return *this += Modular(other);
  }
  template <typename U> Modular &operator-=(const U &other) {
    return *this -= Modular(other);
  }
  Modular &operator++() { return *this += 1; }
  Modular &operator--() { return *this -= 1; }
  Modular operator++(int) {
    Modular result(*this);
    *this += 1;
    return result;
  }
  Modular operator--(int) {
    Modular result(*this);
    *this -= 1;
    return result;
  }
  Modular operator-() const { return Modular(-value); }

  template <typename U = T>
  typename enable_if<is_same<typename Modular<U>::Type, int>::value,
                     Modular>::type &
  operator*=(const Modular &rhs) {
#ifdef _WIN32
    uint64_t x = static_cast<int64_t>(value) * static_cast<int64_t>(rhs.value);
    uint32_t xh = static_cast<uint32_t>(x >> 32), xl = static_cast<uint32_t>(x),
             d, m;
    asm("divl %4; \n\t" : "=a"(d), "=d"(m) : "d"(xh), "a"(xl), "r"(mod()));
    value = m;
#else
    value = normalize(static_cast<int64_t>(value) *
                      static_cast<int64_t>(rhs.value));
#endif
    return *this;
  }
  template <typename U = T>
  typename enable_if<is_same<typename Modular<U>::Type, long long>::value,
                     Modular>::type &
  operator*=(const Modular &rhs) {
    long long q = static_cast<long long>(static_cast<long double>(value) *
                                         rhs.value / mod());
    value = normalize(value * rhs.value - q * mod());
    return *this;
  }
  template <typename U = T>
  typename enable_if<!is_integral<typename Modular<U>::Type>::value,
                     Modular>::type &
  operator*=(const Modular &rhs) {
    value = normalize(value * rhs.value);
    return *this;
  }

  Modular &operator/=(const Modular &other) {
    return *this *= Modular(inverse(other.value, mod()));
  }

  friend const Type &abs(const Modular &x) { return x.value; }

  template <typename U>
  friend bool operator==(const Modular<U> &lhs, const Modular<U> &rhs);

  template <typename U>
  friend bool operator<(const Modular<U> &lhs, const Modular<U> &rhs);

  template <typename V, typename U>
  friend V &operator>>(V &stream, Modular<U> &number);

private:
  Type value;
};

template <typename T>
bool operator==(const Modular<T> &lhs, const Modular<T> &rhs) {
  return lhs.value == rhs.value;
}
template <typename T, typename U>
bool operator==(const Modular<T> &lhs, U rhs) {
  return lhs == Modular<T>(rhs);
}
template <typename T, typename U>
bool operator==(U lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) == rhs;
}

template <typename T>
bool operator!=(const Modular<T> &lhs, const Modular<T> &rhs) {
  return !(lhs == rhs);
}
template <typename T, typename U>
bool operator!=(const Modular<T> &lhs, U rhs) {
  return !(lhs == rhs);
}
template <typename T, typename U>
bool operator!=(U lhs, const Modular<T> &rhs) {
  return !(lhs == rhs);
}

template <typename T>
bool operator<(const Modular<T> &lhs, const Modular<T> &rhs) {
  return lhs.value < rhs.value;
}

template <typename T>
Modular<T> operator+(const Modular<T> &lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) += rhs;
}
template <typename T, typename U>
Modular<T> operator+(const Modular<T> &lhs, U rhs) {
  return Modular<T>(lhs) += rhs;
}
template <typename T, typename U>
Modular<T> operator+(U lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) += rhs;
}

template <typename T>
Modular<T> operator-(const Modular<T> &lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) -= rhs;
}
template <typename T, typename U>
Modular<T> operator-(const Modular<T> &lhs, U rhs) {
  return Modular<T>(lhs) -= rhs;
}
template <typename T, typename U>
Modular<T> operator-(U lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) -= rhs;
}

template <typename T>
Modular<T> operator*(const Modular<T> &lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) *= rhs;
}
template <typename T, typename U>
Modular<T> operator*(const Modular<T> &lhs, U rhs) {
  return Modular<T>(lhs) *= rhs;
}
template <typename T, typename U>
Modular<T> operator*(U lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) *= rhs;
}

template <typename T>
Modular<T> operator/(const Modular<T> &lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) /= rhs;
}
template <typename T, typename U>
Modular<T> operator/(const Modular<T> &lhs, U rhs) {
  return Modular<T>(lhs) /= rhs;
}
template <typename T, typename U>
Modular<T> operator/(U lhs, const Modular<T> &rhs) {
  return Modular<T>(lhs) /= rhs;
}

template <typename T, typename U>
Modular<T> power(const Modular<T> &a, const U &b) {
  assert(b >= 0);
  Modular<T> x = a, res = 1;
  U p = b;
  while (p > 0) {
    if (p & 1)
      res *= x;
    x *= x;
    p >>= 1;
  }
  return res;
}

template <typename T> bool IsZero(const Modular<T> &number) {
  return number() == 0;
}

template <typename T> string to_string(const Modular<T> &number) {
  return to_string(number());
}

// U == std::ostream? but done this way because of fastoutput
template <typename U, typename T>
U &operator<<(U &stream, const Modular<T> &number) {
  return stream << number();
}

// U == std::istream? but done this way because of fastinput
template <typename U, typename T> U &operator>>(U &stream, Modular<T> &number) {
  typename common_type<typename Modular<T>::Type, long long>::type x;
  stream >> x;
  number.value = Modular<T>::normalize(x);
  return stream;
}

/*
using ModType = int;

struct VarMod { static ModType value; };
ModType VarMod::value;
ModType& md = VarMod::value;
using Mint = Modular<VarMod>;
*/

constexpr int md = 998244353;
using Mint = Modular<std::integral_constant<decay<decltype(md)>::type, md>>;

vector<Mint> fact(1, 1);
vector<Mint> inv_fact(1, 1);

Mint C(int n, int k) {
  if (k < 0 || k > n) {
    return 0;
  }
  while ((int)fact.size() < n + 1) {
    fact.push_back(fact.back() * (int)fact.size());
    inv_fact.push_back(1 / fact.back());
  }
  return fact[n] * inv_fact[k] * inv_fact[n - k];
}

const double eps = 1e-9;

bool IsZero(double v) { return abs(v) < 1e-9; }

enum GAUSS_MODE { DEGREE, ABS };

template <typename T>
void GaussianElimination(vector<vector<T>> &a, int limit,
                         GAUSS_MODE mode = DEGREE) {
  if (a.empty() || a[0].empty()) {
    return;
  }
  int h = static_cast<int>(a.size());
  int w = static_cast<int>(a[0].size());
  for (int i = 0; i < h; i++) {
    assert(w == static_cast<int>(a[i].size()));
  }
  assert(limit <= w);
  vector<int> deg(h);
  for (int i = 0; i < h; i++) {
    for (int j = 0; j < w; j++) {
      deg[i] += !IsZero(a[i][j]);
    }
  }
  int r = 0;
  for (int c = 0; c < limit; c++) {
    int id = -1;
    for (int i = r; i < h; i++) {
      if (!IsZero(a[i][c]) &&
          (id == -1 || (mode == DEGREE && deg[i] < deg[id]) ||
           (mode == ABS && abs(a[id][c]) < abs(a[i][c])))) {
        id = i;
      }
    }
    if (id == -1) {
      continue;
    }
    if (id > r) {
      swap(a[r], a[id]);
      swap(deg[r], deg[id]);
      for (int j = c; j < w; j++) {
        a[id][j] = -a[id][j];
      }
    }
    vector<int> nonzero;
    for (int j = c; j < w; j++) {
      if (!IsZero(a[r][j])) {
        nonzero.push_back(j);
      }
    }
    T inv_a = 1 / a[r][c];
    for (int i = r + 1; i < h; i++) {
      if (IsZero(a[i][c])) {
        continue;
      }
      T coeff = -a[i][c] * inv_a;
      for (int j : nonzero) {
        if (!IsZero(a[i][j]))
          deg[i]--;
        a[i][j] += coeff * a[r][j];
        if (!IsZero(a[i][j]))
          deg[i]++;
      }
    }
    ++r;
  }
  for (r = h - 1; r >= 0; r--) {
    for (int c = 0; c < limit; c++) {
      if (!IsZero(a[r][c])) {
        T inv_a = 1 / a[r][c];
        for (int i = r - 1; i >= 0; i--) {
          if (IsZero(a[i][c])) {
            continue;
          }
          T coeff = -a[i][c] * inv_a;
          for (int j = c; j < w; j++) {
            a[i][j] += coeff * a[r][j];
          }
        }
        break;
      }
    }
  }
}

template <typename T> T Determinant(vector<vector<T>> /*&*/ a) {
  if (a.empty()) {
    return T{1};
  }
  assert(a.size() == a[0].size());
  GaussianElimination(a, static_cast<int>(a[0].size()));
  T d{1};
  for (int i = 0; i < (int)a.size(); i++) {
    d *= a[i][i];
  }
  return d;
}

int main() {
  int cnt = 5;
  vector<vector<Mint>> b(cnt, vector<Mint>(cnt));
  for (int i = 0; i < cnt; i++) {
    for (int j = 0; j < cnt; j++) {
      b[i][j] = i + j + i * j + min(100, i * i) / (j + 1);
      cout << b[i][j] << " ";
    }
    cout << endl;
  }
  cout << Determinant(b);
  return 0;
}
