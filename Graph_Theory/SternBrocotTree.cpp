#pragma once

#include "exgcd.hpp"

struct SternBrocotTree {
  struct Fraction {
    long long a = 0, b = 1;  // a / b

    explicit Fraction(long long a = 0, long long b = 1) : a(a), b(b) {}

    Fraction& operator^=(const Fraction& other) { a += other.a; b += other.b; return *this; }
    Fraction operator^(const Fraction& other) const { Fraction ret = *this; return ret ^= other; }
    Fraction scaled(long long k) const { return Fraction(a * k, b * k); }

    bool is_inf() const { return a > 0 && b == 0; }

    bool operator<(const Fraction& other) const {
      if (is_inf()) return false;
      if (other.is_inf()) return true;
      return (__int128)a * other.b < (__int128)b * other.a;
    }
    bool operator>(const Fraction& other) const { return other < *this; }
    bool operator==(const Fraction& other) const { return !(*this < other) && !(other < *this); }
  };

  // A fraction (a, b) is valid in SternBrocotTree iff gcd(a, b) == 1.
  static bool validate(long long a, long long b) { return gcd(a, b) == 1; }
  static bool validate(const Fraction& frac) { return validate(frac.a, frac.b); }

  std::array<Fraction, 2> bounds;

  // NOTE: root contains all the valid fractions, except 0/1 and 1/0.
  static SternBrocotTree get_root() { return SternBrocotTree(Fraction(0, 1), Fraction(1, 0)); }

  SternBrocotTree() = default;
  SternBrocotTree(const Fraction& lbound, const Fraction& rbound) : bounds{lbound, rbound} {}

  Fraction eval() const { return bounds[0] ^ bounds[1]; }
  SternBrocotTree go(int side, long long k) const {
    SternBrocotTree ret = *this;
    ret.bounds[side ^ 1] ^= ret.bounds[side].scaled(k);
    return ret;
  }

  // Returns true if the subtree of current node contains frac.
  bool contains(const Fraction& frac) const {
    return bounds[0] < frac && frac < bounds[1];
  }

  // Returns the maximum k after which the target is still inside the subtree of current node.
  // Any reachable fraction (a, b) can be reached in O(log(max(a, b))) non-trivial `approach` calls.
  long long approach(int side, const Fraction& target) const {
    long long l = 0, r = std::max(target.a, target.b);
    while (l <= r) {
      long long mid = (l + r) >> 1;
      if (this->go(side, mid).contains(target)) l = mid + 1;
      else r = mid - 1;
    }
    return l - 1;
  }
};

std::string to_string(const SternBrocotTree::Fraction& frac) {
  std::stringstream ss;
  ss << frac.a << "/" << frac.b;
  return ss.str();
}

using Fraction = SternBrocotTree::Fraction;
//https://github.com/ssstare/icpc/blob/50752b6077df9a1460ea6c38ebdefc869ec9b980/snippets/stern_brocot_tree.hpp
