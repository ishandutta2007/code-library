#include <bits/stdc++.h>
using namespace std;
using ll = long long;

int digs[] = {2, 3, 5, 7};

ll get_next_number(ll num) {
  // cout<<endl<<num<<endl;
  int ld = num % 10;
  bool exists = find(begin(digs), end(digs), ld) != end(digs);
  if (!exists) {
    while (!exists) {
      num++;
      ld = num % 10;
      exists = find(begin(digs), end(digs), ld) != end(digs);
    }
  } else {
    if (ld == 2)
      num++;
    else if (ld == 3 or ld == 5)
      num += 2;
    else if (ld == 7)
      num += 5;
  }
  int pd = (num % 100) / 10;
  exists = find(begin(digs), end(digs), pd) != end(digs);
  if (!exists) {
    return 10 * get_next_number(num / 10) + 2;
  }
  return num;
}

int main() {
  ll nums[] = {3257737, 3257777, 3257787};
  for (ll num : nums) {
    cout << num;
    for (int i = 0; i < 4; i++) {
      num = get_next_number(num);
      cout << "==>" << num;
    }
    cout << "....." << endl;
  }
  return 0;
}

// 3257737 => 3257752 => 3257753 => 3257755 => 3257757 .....
// 3257777 => 3272222 => 3272223 => 3272225 => 3272227 .....
// 3257787 => 3272222 => 3272223 => 3272225 => 3272227 .....
// https://stackoverflow.com/questions/74134479/given-a-number-find-the-next-higher-number-containing-only-certain-digits
