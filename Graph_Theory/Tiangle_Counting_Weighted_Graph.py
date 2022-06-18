
 
using ll = long long;
void count_triangle(auto edges) {
  ll ret=0;
  std::vector<int> mark(n + 1, -1);
    for (int a = 1; a <= n; ++a) {
      for (auto &&e: edges[a]) mark[e.first] = e.second;
      for (auto &&e: edges[a]) {
        int b = e.first, w1 = e.second;
        for (auto &&ee: edges[b]) {
          int c = ee.first, w2 = ee.second;
          if (mark[c] != -1) {
            ret += w1 * w2 * mark[c];
          }
        }
      }
      for (auto &&e: edges[a]) mark[e.first] = -1;
    }
    return ret;
}
