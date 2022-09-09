
// 1
// 2 3
// 4 5 6
// 7 8 9 10
// 11 12 13 14 15

// or

// 0
// 1 2
// 3 4 5
// 6 7 8 9
// 10 11 12 13 14

#include <bits/stdc++.h>
using namespace std;
using ll = long long;

pair<int, int> row_col_triagular_arragement(int num, bool zero_indexed = true,
                                            bool zero_starting = true) {
  if (zero_starting == false)
    num--;
  int n_row = (-1 + sqrt(1 + 8 * num)) / 2;
  int n_col = num - (n_row * (n_row + 1)) / 2;
  if (zero_indexed)
    return make_pair(n_row, n_col);
  else
    return make_pair(n_row + 1, n_col + 1);
}

int32_t main() {
  ios_base::sync_with_stdio(0);
  cin.tie(0);
  auto arragement = row_col_triagular_arragement(6);
  cout << arragement.first << " " << arragement.second << endl;
  auto arragement2 = row_col_triagular_arragement(6, false);
  cout << arragement2.first << " " << arragement2.second << endl;
  auto arragement3 = row_col_triagular_arragement(6, true, false);
  cout << arragement3.first << " " << arragement3.second << endl;
  return 0;
}