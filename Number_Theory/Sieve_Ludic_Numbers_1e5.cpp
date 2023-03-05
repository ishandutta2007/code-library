#include <bits/stdc++.h>
using namespace std;

vector<int> getLudic(int n) {
  vector<int> ludics;

  for (int i = 1; i <= n; i++)
    ludics.push_back(i);

  for (int index = 1; index < ludics.size(); index++) {
    // Here first item should be included in the list
    // and the deletion is referred by this first item
    // in the loop .
    int first_ludic = ludics[index];

    // Remove_index variable is used to store
    // the next index which we want to delete
    int remove_index = index + first_ludic;

    while (remove_index < ludics.size()) {
      // Removing the next item
      auto it = ludics.begin();
      it = it + remove_index;
      ludics.erase(it);

      // Remove_index is updated so that
      // we get the next index for deletion
      remove_index = remove_index + first_ludic - 1;
    }
  }

  return ludics;
}

int main() {
  int n = 100000;
  vector<int> ans = getLudic(n);
  cout << "[";

  for (int i = ans.size() - 10; i < ans.size() - 1; i++) {
    cout << ans[i] << ", ";
  }

  cout << ans[ans.size() - 1] << "]";
  return 0;
}
