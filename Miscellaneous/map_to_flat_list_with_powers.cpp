#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>

using MyType = unsigned long long;

// The base-number and the max exponent are the elements of a set
struct Element {
  MyType base{};
  MyType maxExponent{};

  // Actually, this elements does not store one value only, but "maxExponent"
  // elements
  // Do get all of those, we use the well known odometer approach
  // We can always calculate the next value in the list and inform, if we have a
  // rollover
  // In case of a roll over, we can inform other sets about this and the other
  // set can also count up
  MyType currentResult{base};
  MyType currentExponent{};

  // Calculate next value. So, will for {2,4} get 2, 4, 8, 16
  bool next(bool overflowFromPrevious) {

    // Indicates that we will wrap around from previous
    bool overflow{};

    // If a previous odometer roll has overflowed
    if (overflowFromPrevious) {

      // Check, if we are NOT at the end. This will happen, if ther is only 1 as
      // the max exponent.
      if (currentExponent < maxExponent) {

        // Get next exponent
        ++currentExponent;

        // And calculate current data. We could also store those values,
        // but we do not want to waste space
        currentResult *= base;

        // Now check again for overflow (after we incremented the exponent)
        if (currentExponent >= maxExponent) {

          // If overlow, then reset exponent counter and base
          currentExponent = 0;
          currentResult = base;
          overflow = true;
        }
      } else
        overflow = true;
    }
    return overflow;
  }
  // Simple inserter operator override, to be able to do an easier output
  friend std::ostream &operator<<(std::ostream &os, const Element &e) {
    return os << '{' << e.base << ',' << e.maxExponent << '}';
  }
};

// We will use a vetcor and not a set, because we need random access via the
// index operator
using MSet = std::vector<Element>;

void getProducts(MSet &mset) {

  // Selectors for creating the subsets of the power set
  std::vector<bool> selector(mset.size());

  // For all number of elements in the original set
  for (size_t k{}; k < mset.size(); ++k) {

    // Set selecot bits. In each loop one more bit. So, 1 then 11, then 111 and
    // so on
    selector[k] = true;

    // Do all permutations of the above set bits
    do {
      // For debug output
      for (bool b : selector)
        std::cout << b * 1;
      std::cout << "  ";
      std::ostringstream oss{};

      // Here we will store all elements of a resulting subset
      MSet subSet{};

      // Check if the selector bit is set, and if so, then add to subset
      for (size_t n{}; n < mset.size(); ++n) {
        if (selector[n]) {
          subSet.push_back(mset[n]);

          oss << mset[n]; // For debug output
        }
      }
      // Debug output of powerset with number of subsets
      std::cout << "Powerset(" << subSet.size() << "): " << std::left
                << std::setw(22) << oss.str() << ' ';

      // Now, we want to calculate all combinations of subsets, using the
      // odometer approach
      // Here we will store the overall number of combinations. It is the
      // product of all max exponents
      MyType combinations{1};
      for (const Element &element : subSet)
        combinations *= element.maxExponent;

      // Now get the product for all combinations over all subsets
      for (MyType i{}; i < combinations; ++i) {

        // Get the product for one combination
        MyType product{1};
        for (Element &element : subSet)
          product *= element.currentResult;

        std::cout << product << ' '; // For debug output

        // And, using the odometer approach, create the next combination
        bool overflow{true};
        for (Element &element : subSet)
          overflow = element.next(overflow);
      }
      std::cout << '\n'; // For debug output
    } while (std::prev_permutation(selector.begin(), selector.end()));
  }
}
// Test / driver code
int main() {
  MSet mset{{2, 4}, {3, 1}, {5, 3}};

  getProducts(mset);
}
// https://stackoverflow.com/questions/3644858/generating-all-factors-of-a-number-given-its-prime-factorization/3644903#3644903
// https://stackoverflow.com/questions/72220384/generating-all-permuated-products-from-a-map-is-not-working
