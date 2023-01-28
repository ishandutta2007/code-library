// C++ program to find the N-th number whose
// digital root is X

#include <bits/stdc++.h>
using namespace std;

// Function to find the digital root of
// a number
int findDigitalRoot(int num) {
  int sum = INT_MAX, tempNum = num;

  while (sum >= 10) {
    sum = 0;

    while (tempNum > 0) {
      sum += tempNum % 10;
      tempNum /= 10;
    }

    tempNum = sum;
  }

  return sum;
}

// Function to find the Nth number with
// digital root as X
void findAnswer(int X, int N) {
  // Counter variable to keep the
  // count of valid numbers
  int counter = 0;

  for (int i = 1; counter < N; ++i) {

    // Find digital root
    int digitalRoot = findDigitalRoot(i);

    // Check if is required answer or not
    if (digitalRoot == X) {
      ++counter;
    }

    // Print the answer if you have found it
    // and breakout of the loop
    if (counter == N) {
      cout << i;

      break;
    }
  }
}

// Driver Code
int main() {
  int X = 1, N = 3;

  findAnswer(X, N);

  return 0;
}