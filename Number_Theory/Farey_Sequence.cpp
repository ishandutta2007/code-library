// ULInt numeber_of_fractions(ULInt n_max,
//   ULInt lower_numerator=0, ULInt lower_denominator=1,
//   ULInt upper_numerator=1, ULInt upper_denominator=1)
// /* Using Farey sequence (see: https://en.wikipedia.org/wiki/Farey_sequence).
//  **/
// {
//   ULInt a, b, c, d;
//   a = lower_numerator;
//   b = lower_denominator;
//   c = upper_numerator;
//   d = upper_denominator;

//   ULInt count = 0;
//   ULInt k, temp_a, temp_b;
//   while (!(c == 1 && d == 2)) {
//     k = (n_max + b) / d;
//     temp_a = a;
//     temp_b = b;
//     a = c;
//     b = d;
//     c = k * c - temp_a;
//     d = k * d - temp_b;
//     // std::cout << "\t" << a << "/" << b << std::endl;
//     count += 1;
//   }

//   return count;
// }

/*
  Farey Sequence Generator:

  Author: Zachary Friggstad
  Date: Aug, 2008

  ---------------------------------------------------------------

  The Farey Sequence of order n is the list of
  all reduced fractions between 0 and 1 (inclusive)
  in sorted order.

  e.g. order 6:
  0/1, 1/6, 1/5, 1/4, 2/5, 1/3, 1/2, 2/3, 3/5, 3/4, 4/5, 5/6, 1/1

  Given any positive integer n, this algorithm will generate
  the Farey sequence in order with one term being generated per
  loop iteration.

  ---------------------------------------------------------------

  Complexity:
   - linear time in the size of the output
   - constant space

  Reference:
    "Introduction to the Theory of Numbers" - Hardy & Wright

  Reliablity:
    UVa - 10408
*/

#include <algorithm>
#include <iostream>

using namespace std;

void farey(int n) {
  int h = 0, k = 1, x = 1, y = 0;

  do {
    cout << h << '/' << k << endl;

    int r = (n - y) / k;
    y += r * k;
    x += r * h;

    swap(x, h);
    swap(y, k);
    x = -x;
    y = -y;
  } while (k > 1);
  cout << "1/1" << endl;
}

int main() {
  int n;
  while (cin >> n)
    farey(n);
}