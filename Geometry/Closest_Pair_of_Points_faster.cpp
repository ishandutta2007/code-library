#include <bits/stdc++.h>
// clang-format off
#define  GODSPEED       ios:: sync_with_stdio(0);cin.tie(0);cout.tie(0);cout<<fixed;cout<<setprecision(15);
#define ASSERT(minval, variable, maxval) {if (!(minval <= variable && variable <= maxval)) {cout << "ASSERT FAILED: " << minval << " " << variable << " " << maxval<< " @ " << __FILE__ << " (" << __LINE__ << ")" << endl; exit(0);  }}
#define ASSERT2(condition) {if(!(condition)){ std::cerr << "ASSERT FAILED: " << #condition << " @ " << __FILE__ << " (" << __LINE__ << ")" << std::endl; } }
// clang-format on

using namespace std;
using ll = long long;
using i64 = long long;
using i128 = __int128;
using u64 = unsigned long long;
using u128 = __uint128_t;
using f64 = double;
using f128 = __float128;

const int MAXN = 1e6;
const int MOD = 1e9 + 7;
const ll INF = 0x7f7f7f7f7f7f7f7f;
const int INFi = 0x7f7f7f7f; // 1e9+7;

typedef complex<double> Point;

struct PointLtX {
  bool operator()(const Point &a, const Point &b) const {
    if (real(a) < real(b)) {
      return true;
    }
    return real(a) == real(b) ? imag(a) < imag(b) : false;
  }
};

struct PointLtY {
  bool operator()(const Point &a, const Point &b) const {
    if (imag(a) < imag(b)) {
      return true;
    }
    return imag(a) == imag(b) ? real(a) < real(b) : false;
  }
};

// pair<Point, Point>
pair<int, int> findClosestPair(vector<Point> &points) {
  sort(points.begin(), points.end(), PointLtX());
  set<Point, PointLtY> candidatePoints;
  pair<Point, Point> closestPair;
  pair<int, int> closestPairIndex;
  double minimumDistance = INFINITY;
  size_t j = 0;
  for (size_t i = 0; i < points.size(); ++i) {
    while (real(points[i]) - real(points[j]) >= minimumDistance && i > j) {
      Point pp = points[j];
      // if (candidatePoints.size()>0 and candidatePoints.find(pp) !=
      // candidatePoints.end())
      candidatePoints.erase(pp);
      ++j;
    }
    double frontYCoordinate = imag(points[i]);
    auto hi = candidatePoints.upper_bound(
        Point(INFINITY, frontYCoordinate + minimumDistance));
    auto lo = candidatePoints.lower_bound(
        Point(INFINITY, frontYCoordinate - minimumDistance));
    for (; lo != hi; ++lo) {
      double dist = abs(*lo - points[i]);
      if (dist < minimumDistance) {
        minimumDistance = dist;
        // closestPair = pair<Point, Point>(*lo, points[i]);
        closestPairIndex =
            pair<int, int>(distance(candidatePoints.begin(), lo), i);
      }
    }
    candidatePoints.insert(points[i]);
  }
  // return closestPair;//You can return values directly if you want or if you
  // have x[], y[]array in solve() you can fetch values from there
  return closestPairIndex;
}

void solve() {
  int N;
  double x[N], y[N];
  vector<Point> points;
  while (scanf("%d", &N), N != 0) {
    for (int i = 0; i < N; ++i) {
      scanf("%lf %lf", &x[i], &y[i]);
      points.emplace_back(x[i], y[i]);
    }
    auto close = findClosestPair(points);
    printf("By Index:%d %d\n", close.first, close.second);
    printf("By Value:%0.2lf %0.2lf %0.2lf %0.2lf\n", x[close.first],
           y[close.first], x[close.second], y[close.second]);
    points.clear();
  }
}

int main(int argc, char **argv) {
  GODSPEED
  // Filenames f = get_filenames(string(argv[0]), "huge");
  // freopen(f.infile_name.c_str(), "r", stdin);
  // freopen(f.outfile_name.c_str(), "w", stdout);
  int test;
  test = 1; //
  for (int i = 1; i <= test; i++) {
    solve();
  }
  return 0;
}
