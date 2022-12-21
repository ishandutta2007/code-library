// Function to find the minimum number of Multiplication
//  steps required in multiplying chain of n matrices.
int minimum_Cost_Polygon_Triangulation(int arr[], int low, int high,
                                       vector<vector<int>> &memo) {
  // If we are left with one matrix then
  // we don't need any multiplication steps.
  if (low == high)
    return 0;

  if (memo[low][high] != -1)
    return memo[low][high];

  int minCost = INT_MAX;

  // Iterating from low to high - 1
  for (int k = low; k < high; k++) {

    /*
         Cost = Cost of Multiplying chain on left side +
                Cost of Multiplying chain on right side +
                Cost of Multiplying matrix obtained from left
                and right side.
      */
    int leftCost = solve(arr, low, k, memo);
    int rightCost = solve(arr, k + 1, high, memo);
    int mergeCost = arr[low - 1] * arr[k] * arr[high];

    int tempAns = leftCost + rightCost + mergeCost;

    // If sum of leftCost, rightCost and mergeCost
    // is smaller than minCost then update it.
    minCost = min(minCost, tempAns);
    memo[low][high] = minCost;
  }

  return minCost;
}