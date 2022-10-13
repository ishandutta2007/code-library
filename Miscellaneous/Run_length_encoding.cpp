#include <bits/stdc++.h>
using namespace std;
string decode(string s) {
  string decoded = "";
  string num = "";
  char lastchar = '\0';
  for (auto i : s) {
    if (isalpha(i)) {
      if (lastchar != '\0' and num != "") {
        decoded = decoded + string(stoi(num), lastchar);
      }
      num = "";
      lastchar = i;
    } else {
      num += i;
    }
  }
  if (lastchar != '\0' and num != "") {
    decoded = decoded + string(stoi(num), lastchar);
  }
  return decoded;
}

pair<char, pair<int, int>> longest_run_length(string s) {
  // string decoded = "";
  string num = "";
  char lastchar = '\0';
  int lrc_si = -1;
  int trl = 0;
  int lrl = 0;
  char lrc = '\0';
  for (auto i : s) {
    if (isalpha(i)) {
      if (lastchar != '\0' and num != "") {
        // decoded = decoded + string(stoi(num), lastchar);
        if (stoi(num) > lrl) {
          lrl = stoi(num);
          lrc = lastchar;
          lrc_si = trl;
        }
        trl += stoi(num);
      }
      num = "";
      lastchar = i;
    } else {
      num += i;
    }
  }
  if (lastchar != '\0' and num != "") {
    // decoded = decoded + string(stoi(num), lastchar);
    if (stoi(num) > lrl) {
      lrl = stoi(num);
      lrc = lastchar;
      lrc_si = trl;
    }
    trl += stoi(num);
  }
  return make_pair(lrc, make_pair(lrc_si, lrl));
}

string encode(string str) {
  string encoded = "";
  int n = str.length();
  for (int i = 0; i < n; i++) {
    int count = 1;
    while (i < n - 1 && str[i] == str[i + 1]) {
      count++;
      i++;
    }
    encoded = encoded + str[i] + to_string(count);
  }
  return encoded;
}
int main() {
  string str = "wwwwaaadexxxxxxywww";
  auto m = encode(str);
  cout << m << endl;

  auto lrlp = longest_run_length(m);
  cout << lrlp.first << " : " << lrlp.second.first
       << "(starting index 0-indexed), " << lrlp.second.second << "(run-length)"
       << endl;

  auto n = decode(m);
  cout << n << endl;

  cout << (n == str) << endl;
  return 0;
}
