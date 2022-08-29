#include <bits/stdc++.h>
using namespace std;

namespace tensor {
template <typename T, int N, int... args>
struct static_tensor : std::array<static_tensor<T, args...>, N> {
  static_tensor() : std::array<static_tensor<T, args...>, N>() {}

  void fill(const T &val) {
    for (int i = 0; i < N; i++) {
      (*this)[i].fill(val);
    }
  }
};

template <typename T, int N> struct static_tensor<T, N> : std::array<T, N> {
  static_tensor() : std::array<T, N>() {}
};

template <typename T, int D>
struct dynamic_tensor : std::vector<dynamic_tensor<T, D - 1>> {
  dynamic_tensor() : std::vector<dynamic_tensor<T, D - 1>>() {}

  template <typename... args>
  dynamic_tensor(size_t n, args... A)
      : std::vector<dynamic_tensor<T, D - 1>>(n,
                                              dynamic_tensor<T, D - 1>(A...)) {}

  dynamic_tensor(const std::initializer_list<dynamic_tensor<T, D - 1>> &list)
      : std::vector<dynamic_tensor<T, D - 1>>(list) {}

  void fill(const T &val) {
    for (auto it = this->begin(); it != this->end(); it++) {
      it->fill(val);
    }
  }
};

template <typename T> struct dynamic_tensor<T, 1> : std::vector<T> {
  dynamic_tensor() : std::vector<T>() {}

  dynamic_tensor(size_t n, const T &val = T()) : std::vector<T>(n, val) {}

  dynamic_tensor(const std::initializer_list<T> &list) : std::vector<T>(list) {}

  void fill(const T &val = T()) { std::fill(this->begin(), this->end(), val); }
};
}

int main() {
  ios::sync_with_stdio(false);
  cin.tie(nullptr);

  int T;
  cin >> T;

  while (T--) {
    int n, k;
    cin >> n >> k;

    vector<int> p(n, -1);
    vector<vector<int>> g(n);
    for (int i = 1; i < n; i++) {
      cin >> p[i], --p[i];
      g[p[i]].push_back(i);
    }

    vector<int> sz(n, 1), d(n, 1);
    auto dfs = [&](int x, const auto &self) -> void {
      for (int i : g[x]) {
        d[i] = 1 + d[x];
        self(i, self);
        sz[x] += sz[i];
      }
    };

    dfs(0, dfs);

    vector<int> path;
    int mx = max_element(d.begin(), d.end()) - d.begin();
    for (int i = mx; i != -1; i = p[i]) {
      path.push_back(i);
    }

    reverse(path.begin(), path.end());
    if (path.size() >= k) {
      cout << k - 1 << '\n';
      for (int i = 0; i < k; i++) {
        cout << path[i] + 1 << ' ';
      }
      cout << '\n';
    } else {
      vector<int> ans;
      auto recurse = [&](int v, const auto &self) -> void {
        ans.push_back(v), k -= 1;
        for (int u : g[v]) {
          if (k > 0) {
            self(u, self);
            ans.push_back(v);
          }
        }
      };

      k -= path.size();
      for (int i = 0; i + 1 < path.size(); i++) {
        ans.push_back(path[i]);
        for (int v : g[path[i]]) {
          if (v != path[i + 1] && k > 0) {
            recurse(v, recurse);
            ans.push_back(path[i]);
          }
        }
      }
      ans.push_back(path.back());
      cout << (int)ans.size() - 1 << '\n';
      for (auto &x : ans)
        cout << x + 1 << ' ';
      cout << '\n';
    }
  }
  return 0;
}

// https://codeforces.com/problemset/problem/1510/G