#include <bits/stdc++.h>
using namespace std;

using Real = double;
using Point = complex<Real>;
const Real EPS = 1e-8, PI = acos(-1);
using Points = vector<Point>;

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_5_A

Real closest_pair(Points ps) {
  if (ps.size() <= 1)
    throw(0);
  sort(begin(ps), end(ps));

  auto compare_y = [&](const Point &a, const Point &b) {
    return imag(a) < imag(b);
  };
  vector<Point> beet(ps.size());
  const Real INF = 1e18;

  function<Real(int, int)> rec = [&](int left, int right) {
    if (right - left <= 1)
      return INF;
    int mid = (left + right) >> 1;
    auto x = real(ps[mid]);
    auto ret = min(rec(left, mid), rec(mid, right));
    inplace_merge(begin(ps) + left, begin(ps) + mid, begin(ps) + right,
                  compare_y);
    int ptr = 0;
    for (int i = left; i < right; i++) {
      if (abs(real(ps[i]) - x) >= ret)
        continue;
      for (int j = 0; j < ptr; j++) {
        auto luz = ps[i] - beet[ptr - j - 1];
        if (imag(luz) >= ret)
          break;
        ret = min(ret, abs(luz));
      }
      beet[ptr++] = ps[i];
    }
    return ret;
  };
  return rec(0, (int)ps.size());
}

// https://github.com/ishandutta2007/codeforces/blob/b90eeeefaeb7a4933a8de26b7b8194dbc4ab7f64/noimi/normal/772/B.cpp
