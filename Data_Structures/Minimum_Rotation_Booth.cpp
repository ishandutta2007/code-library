#include <bits/extc++.h>

using namespace std;
using namespace __gnu_pbds;
using ll = long long;

namespace IO {
template <typename T> inline bool can_read(T &t) {
  int n = 0;
  int ch = getchar_unlocked();
  while (!isdigit(ch)) {
    if (ch == EOF)
      return false;
    n |= ch == '-';
    ch = getchar_unlocked();
  }
  t = 0;
  while (isdigit(ch)) {
    t = t * 10 + ch - 48;
    ch = getchar_unlocked();
  }
  if (n)
    t = -t;
  return true;
}
template <typename T, typename... Args>
inline bool can_read(T &t, Args &... args) {
  return can_read(t) && can_read(args...);
}
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

// Get the minimal rotation of the vector<T>/string
// Booth's algorithm : O(s.size())
int minimalRotation(vector<int> s) {
  s.resize(s.size() << 1);
  copy(s.begin(), s.begin() + s.size() / 2, s.begin() + s.size() / 2);
  int length = (int)s.size();
  vector<int> f(length, -1); // Failure function
  int k = 0;                 // Least rotation of string found so far
  for (int j = 1, i; j < length; ++j) {
    char sj = s.at(j);
    i = f.at(j - k - 1);
    while (i != -1 && sj != s.at(k + i + 1)) {
      if (sj < s.at(k + i + 1)) {
        k = j - i - 1;
      }
      i = f.at(i);
    }
    if (sj != s.at(k + i + 1)) {
      if (sj < s.at(k)) {
        k = j;
      }
      f.at(j - k) = -1;
    } else {
      f.at(j - k) = i + 1;
    }
  }
  return k;
}

int main() {

  int n;
  IO::read(n);
  vector<int> p(n);
  for (auto &i : p)
    IO::read(i);
  rotate(p.begin(), p.begin() + minimalRotation(p), p.end());
  for (const auto &i : p) {
    IO::write(i);
    putchar_unlocked(' ');
  }

  return 0;
}
// Luogu/P1368
