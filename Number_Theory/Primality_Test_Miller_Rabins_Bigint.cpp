#pragma once
#include <cassert>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

class BigInt {
private:
  int base = 10;
  void str_init(string);
  static AbstractMulter *multer;
  static AbstractDivider *divider;

public:
  int sign = 1;

  vector<int> numbers;

  BigInt() = default;
  BigInt(string);
  BigInt(const BigInt &);
  BigInt(vector<int>::iterator, vector<int>::iterator, int, int);
  BigInt(long long);

  static BigInt add(BigInt, BigInt);
  static BigInt sub(BigInt, BigInt);
  static BigInt binpow(BigInt, int);
  static BigInt binpow(BigInt, int, BigInt);
  static BigInt binpow(BigInt, BigInt);
  static BigInt binpow(BigInt, BigInt, BigInt);
  static BigInt gcd(BigInt, BigInt);
  static AbstractMulter &get_multer();
  static AbstractDivider &get_divider();
  static void set_multer(AbstractMulter *);
  static void set_divider(AbstractDivider *);

  int get_base();
  void set_base(int);

  size_t size();
  void normalize();
  void resize(int);
  vector<int>::iterator begin();
  vector<int>::iterator end();
  BigInt &left_shift(int);
  BigInt &right_shift(int);
  bool even();
  bool odd();

  int &operator[](int i);
  friend ostream &operator<<(ostream &, const BigInt &);
  friend istream &operator>>(istream &, BigInt &);
  friend BigInt operator-(BigInt);
  friend BigInt operator+(BigInt);
  friend bool operator==(BigInt, BigInt);
  friend bool operator!=(BigInt, BigInt);
  friend bool operator>(BigInt, BigInt);
  friend bool operator<(BigInt, BigInt);
  friend bool operator>=(BigInt, BigInt);
  friend bool operator<=(BigInt, BigInt);
  friend BigInt operator+(BigInt, BigInt);
  friend BigInt operator-(BigInt, BigInt);
  friend BigInt operator*(BigInt, int);
  friend BigInt operator*(int, BigInt);
  friend BigInt operator*(BigInt, BigInt);
  friend BigInt operator/(BigInt, int);
  friend BigInt operator/(BigInt, BigInt);
  friend BigInt operator%(BigInt, BigInt);
};

class AbstractMulter {
public:
  virtual BigInt mult(BigInt &, BigInt &) = 0;
};

class AbstractDivider {
public:
  virtual pair<BigInt, BigInt> div(BigInt, BigInt) = 0;
};

class LCG {
private:
  BigInt x;
  BigInt a;
  BigInt c;
  BigInt lower_bound;
  BigInt upper_bound;
  BigInt _upper_bound;

public:
  LCG(BigInt x);
  BigInt next_int();

  void set_x(BigInt x);
  void set_a(BigInt a);
  void set_c(BigInt c);
  void set_lower_bound(BigInt lb);
  void set_upper_bound(BigInt ub);

  BigInt &get_x();
  BigInt &get_a();
  BigInt &get_c();
  BigInt &get_lower_bound();
  BigInt &get_upper_bound();
};

// Miller–Rabin primality test
class MRT {
private:
  static LCG *generator;
  static bool _test(BigInt &, BigInt &);

public:
  static bool is_prime(BigInt n, int k = 50);
};

LCG *MRT::generator = new LCG(42);

bool MRT::_test(BigInt &d, BigInt &n) {
  BigInt a(generator->next_int());
  BigInt x = BigInt::binpow(a, d, n);

  if (x == BigInt(1) || x == n - 1)
    return true;

  while (d != n - 1) {
    x = (x * x) % n;
    d = d * BigInt(2);

    if (x == BigInt(1))
      return false;
    if (x == n - 1)
      return true;
  }
  return false;
}

bool MRT::is_prime(BigInt n, int k) {
  if (n != BigInt(2) && n.even())
    return false;

  if (n <= BigInt(1) || n == BigInt(4))
    return false;

  if (n < BigInt(4))
    return true;

  generator->set_lower_bound(BigInt(2));
  generator->set_upper_bound(BigInt(n - 2));

  BigInt d(n - 1);
  while (d.even()) // d % 2 == 0
    d = d / BigInt(2);

  for (int i = 0; i < k; ++i) {
    if (!_test(d, n))
      return false;
  }

  return true;
}
