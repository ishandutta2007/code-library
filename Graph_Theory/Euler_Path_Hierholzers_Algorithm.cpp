// https://open.kattis.com/problems/eulerianpath
// An example that uses Hierholzer’s Algorithm to find an Eulerian Cycle/Path in
// a digraph
#include <bits/extc++.h>

using namespace std;
using namespace __gnu_pbds;

namespace IO {
template <typename T> inline void read(T &t) {
  int n = 0;
  int c = getchar_unlocked();
  t = 0;
  while (!isdigit(c))
    n |= c == '-', c = getchar_unlocked();
  while (isdigit(c))
    t = t * 10 + c - 48, c = getchar_unlocked();
  if (n)
    t = -t;
}

template <typename T, typename... Args> inline void read(T &t, Args &... args) {
  read(t);
  read(args...);
}

template <typename T> inline void write(T x) {
  if (x < 0)
    x = -x, putchar_unlocked('-');
  if (x > 9)
    write(x / 10);
  putchar_unlocked(x % 10 + 48);
}

template <typename T> inline void writeln(T x) {
  write(x);
  putchar_unlocked('\n');
}
}

struct custom_hash {
  static uint64_t splitmix64(uint64_t x) {
    // http://xorshift.di.unimi.it/splitmix64.c
    x += 0x9e3779b97f4a7c15;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
    return x ^ (x >> 31);
  }

  size_t operator()(uint64_t x) const {
    static const uint64_t FIXED_RANDOM =
        chrono::steady_clock::now().time_since_epoch().count();
    return splitmix64(x + FIXED_RANDOM);
  }
};

namespace DSU {
vector<int> father;
vector<int> Size;

void init(int n) {
  vector<int>().swap(father);
  vector<int>().swap(Size);
  father.resize(n + 5);
  iota(father.begin(), father.end(), 0);
  Size.resize(n + 5, 1);
}

int find(int x) {
  if (x != father[x]) {
    father[x] = find(father[x]);
  }
  return father[x];
}

void merge(int x, int y) {
  x = find(x);
  y = find(y);
  if (x == y) {
    return;
  }
  if (Size[x] > Size[y]) {
    swap(x, y);
  }
  Size[y] += Size[x];
  father[x] = y;
}
}

vector<vector<int>> adjList;
vector<int> in_degree;
vector<int> out_degree;
vector<int> answer;
int n, m;

bool has_eulerian_path(const vector<int> &A, int &start_vertex) {
  int a = 0;
  int b = 0;
  int node_with_extra_out_deg = -1;
  int node_with_equal_in_deg_out_deg = -1;
  for (const auto &i : A) {
    if (in_degree.at(i) == out_degree.at(i)) {
      node_with_equal_in_deg_out_deg = i;
      continue;
    } else if (in_degree.at(i) - out_degree.at(i) == 1) {
      ++a;
    } else if (out_degree.at(i) - in_degree.at(i) == 1) {
      ++b;
      node_with_extra_out_deg = i;
    } else {
      return false;
    }
  }
  // Case 1: an Eulerian cycle exists in this CC
  if (a == 0 && b == 0) {
    start_vertex = A.front();
    return true;
  }
  // Case 2: an Eulerian path exists in this CC
  if (a <= 1 && b <= 1) {
    if (node_with_extra_out_deg != -1) {
      start_vertex = node_with_extra_out_deg;
    } else {
      start_vertex = node_with_equal_in_deg_out_deg;
    }
    return true;
  }
  // Case 3: a or b is > 1, no Eulerian cycle/path found
  return false;
}

// The adjList will be changed
vector<int> Hierholzer(int u) {

  vector<int> result;
  stack<int> s;

  s.push(u);
  while (!s.empty()) {
    u = s.top();
    if (!adjList.at(u).empty()) { // still has neighbor
      s.push(adjList.at(u).back());
      adjList.at(u).pop_back();
    } else {
      result.emplace_back(u);
      s.pop();
    }
  }

  reverse(result.begin(), result.end());
  return result;
}

int main() {

  int u, v, start_vertex;
  while (true) {
    IO::read(n, m);
    if (n == 0 && m == 0) {
      break;
    }
    DSU::init(n);
    adjList.resize(n + 5);
    in_degree.resize(n + 5);
    out_degree.resize(n + 5);
    for (int i = 0; i < m; ++i) {
      IO::read(u, v);
      DSU::merge(u, v);
      adjList.at(u).emplace_back(v);
      ++out_degree.at(u);
      ++in_degree.at(v);
    }
    gp_hash_table<int, vector<int>, custom_hash> groups;
    for (int i = 0; i < n; ++i) {
      groups[DSU::find(i)].emplace_back(i);
    }
    bool has_answer = false;
    for (const auto & [ i, j ] : groups) {
      // Skip all nodes that has zero degrees
      if (j.size() == 1 && in_degree.at(j.front()) == 0 &&
          out_degree.at(j.front()) == 0) {
        continue;
      }
      // Stop the iteration if an Eulerian circuit/path found in a CC
      if (has_eulerian_path(j, start_vertex)) {
        answer = Hierholzer(start_vertex);
        has_answer = true;
        break;
      }
    }

    // output
    if (!has_answer) {
      puts("Impossible");
    } else {
      // Eulerian path/cycle in a subgraph of G is not allowed!
      // Since we are trying to find an Eulerian path/cycle for the whole G!
      // Thus the edges of the path must equal to the edges in G
      if ((int)answer.size() != m + 1) {
        puts("Impossible");
      } else {
        for (size_t i = 0; i != answer.size(); ++i) {
          if (i != 0) {
            putchar_unlocked(' ');
          }
          IO::write(answer.at(i));
        }
        putchar_unlocked('\n');
      }
    }
    // clear
    decltype(adjList)().swap(adjList);
    decltype(in_degree)().swap(in_degree);
    decltype(out_degree)().swap(out_degree);
    decltype(answer)().swap(answer);
  }
  return 0;
}
// Kattis/eulerianpath
