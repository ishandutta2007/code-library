#include <bits/stdc++.h>
// Credits: AvinashKartik
using namespace std;

#define GODSPEED                                                               \
  ios::sync_with_stdio(0);                                                     \
  cin.tie(0);                                                                  \
  cout.tie(0);                                                                 \
  cout << fixed;                                                               \
  cout << setprecision(15);
#define f first
#define s second
#define newl cout << "\n";
#define pb push_back
#define mset(a, x) memset(a, x, sizeof(a))
#define debv(a)                                                                \
  for (auto it : a)                                                            \
    cout << it << " ";                                                         \
  newl;
#define deb1(a) cout << a << "\n";
#define deb2(a, b) cout << a << " " << b << "\n";
#define deb3(a, b, c) cout << a << " " << b << " " << c << "\n";
#define deb4(a, b, c, d) cout << a << " " << b << " " << c << " " << d << "\n";
#define uniq(a) a.resize(unique(a.begin(), a.end()) - a.begin());
#define all(a) a.begin(), a.end()

typedef long long ll;
typedef long double ld;

const ll N = 2e5 + 5;
const ll mod = 1e9 + 7;
const ll INF = 0x7f7f7f7f7f7f7f7f;
const int INFi = 0x7f7f7f7f;

template <typename T> struct bitwise_trie {
  int n, b;
  stack<int> s;
  vector<T> trie;

  bitwise_trie(int n, int b) : n(n), b(b) {
    for (int i = n - 1; i >= 1; i--)
      s.push(i);
    trie.resize(n);
  }

  void insert(int x) {
    int v = 0;
    for (int i = b; i >= 0; i--) {
      int k = (x >> i) & 1;
      if (trie[v].go[k] == -1) {
        trie[v].go[k] = s.top();
        s.pop();
      }
      trie[v].cnt++;
      v = trie[v].go[k];
    }
    trie[v].cnt++;
  }

  void erase(int x) {
    int v = 0;
    for (int i = b; i >= 0; i--) {
      int k = (x >> i) & 1;
      trie[v].cnt--;
      int tmp = trie[v].go[k];
      if (trie[tmp].cnt == 1) {
        s.push(tmp);
        trie[v].go[k] = -1;
      }
      v = tmp;
    }
    trie[v].cnt--;
  }

  int min_xor(int x, int f = 0) {
    int v = 0, ans = 0;
    for (int i = b; i >= 0; i--) {
      int k = (x >> i) & 1;
      if (trie[v].go[k ^ f] != -1) {
        ans += (1 << i) * f;
        v = trie[v].go[k ^ f];
      } else {
        ans += (1 << i) * (1 ^ f);
        v = trie[v].go[k ^ 1 ^ f];
      }
    }
    return ans;
  }

  int max_xor(int x) { return min_xor(x, 1); }
};

struct node {
  vector<int> go = {-1, -1};
  int cnt = 0;
};

int l, r, a[N];
bitwise_trie<node> t((1 << 18), 17);

void solve() {
  cin >> l >> r;
  for (int i = 1; i <= r - l + 1; i++) {
    cin >> a[i];
    t.insert(a[i]);
  }

  for (int i = 1; i <= r - l + 1; i++) {
    int x = a[i] ^ l;
    if (t.max_xor(x) != r || t.min_xor(x) != l)
      continue;
    deb1(x) break;
  }

  for (int i = 1; i <= r - l + 1; i++) {
    t.erase(a[i]);
  }
}

int main() {
  GODSPEED;
  int test = 1;
  cin >> test;
  for (int i = 1; i <= test; i++) {
    solve();
  }
}

// https://codeforces.com/contest/1658/problem/D2