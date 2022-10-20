// C++ program to count
// trailing 0s in n!
#include <iostream>
using namespace std;

int countTrailingZeros(int n) {
  if (n < 0) // Negative Number Edge Case
    return -1;
  int count = 0;
  for (int i = 5; n / i >= 1; i *= 5)
    count += n / i;
  return count;
}

int main() {
  int n = 100;
  cout << "Count of trailing 0s in " << 100 << "! is " << countTrailingZeros(n);
  return 0;
}
