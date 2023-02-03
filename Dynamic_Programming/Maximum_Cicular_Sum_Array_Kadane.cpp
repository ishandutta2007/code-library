#include <iostream>
#include <climits>

using namespace std;

int kadane(int array[], int n);

int main() {
  int n;
  cin >> n;

  int array[n];

  for (int i = 0; i < n; i++) {
    cin >> array[i];
  }

  int maxSum = 0;

  int kadanePositive = kadane(array, n);

  int arraySum = 0;

  for (int i = 0; i < n; i++) {
    arraySum += array[i];
  }

  for (int i = 0; i < n; i++) {
    array[i] = -array[i];
  }

  int negativeKadane = arraySum + kadane(array, n);

  maxSum = max(kadanePositive, negativeKadane);

  cout << "Maximum Sum: " << maxSum << endl;
}

int kadane(int array[], int n) {
  int maxKadane, currKadaneSum;
  int kadaneSum[n];

  maxKadane = 0;
  currKadaneSum = 0;

  for (int i = 0; i < n; i++) {

    currKadaneSum = currKadaneSum + array[i];
    if (currKadaneSum < 0) {
      kadaneSum[i] = 0;
      currKadaneSum = 0;
    } else {
      kadaneSum[i] = currKadaneSum;
    }

    maxKadane = max(maxKadane, currKadaneSum);
  }

  return maxKadane;
}