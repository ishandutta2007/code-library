#include <bits/stdc++.h>
using namespace std;
typedef pair<int, int> ii;
typedef long long LL;

struct Segment {
  int x, y1, y2, id, v;
  Segment(int x, int y1, int y2, int id, int v)
      : x(x), y1(y1), y2(y2), id(id), v(v) {}
};

int getYValue[60010];
map<int, int> getYIdx;

LL S[240010];
int lazy[240010];

void setSum(const int L, const int R, const int idx) {
  if (lazy[idx] > 0)
    S[idx] = getYValue[R + 1] - getYValue[L];
  else
    S[idx] = (L != R) ? S[idx * 2] + S[idx * 2 + 1] : 0;
}

void update(const int l, const int r, const int L, const int R, const int v,
            const int idx = 1) {
  if (l <= L && r >= R)
    lazy[idx] += v;
  else {
    int M = (L + R) / 2;
    if (l <= M)
      update(l, r, L, M, v, idx * 2);
    if (r > M)
      update(l, r, M + 1, R, v, idx * 2 + 1);
  }
  setSum(L, R, idx);
}

int main() {
  int T, N, x1, y1, x2, y2;
  scanf("%d", &T);
  for (int t = 1; t <= T; t++) {
    scanf("%d", &N);
    vector<int> yCoords;
    vector<Segment> segments;
    getYIdx.clear();

    for (int i = 0; i < N; i++) {
      scanf("%d %d %d %d", &x1, &y1, &x2, &y2);
      yCoords.push_back(y1);
      yCoords.push_back(y2);
      segments.push_back(Segment(x1, y1, y2, i, 1));
      segments.push_back(Segment(x2, y1, y2, i, -1));
    }

    sort(segments.begin(), segments.end(),
         [](const Segment &left, const Segment &right) -> bool {
           return make_pair(left.x, ii(left.id, left.v)) <=
                  make_pair(right.x, ii(right.id, right.v));
         });

    sort(yCoords.begin(), yCoords.end());
    yCoords.erase(unique(yCoords.begin(), yCoords.end()), yCoords.end());
    for (int i = 0; i < int(yCoords.size()); i++) {
      getYValue[i] = yCoords[i];
      getYIdx[yCoords[i]] = i;
    }

    memset(S, 0, sizeof(S));
    memset(lazy, 0, sizeof(lazy));

    LL ans = 0;
    int prvX = 0;

    for (auto &segment : segments) {
      ans += (segment.x - prvX) * S[1];
      update(getYIdx[segment.y1], getYIdx[segment.y2] - 1, 0,
             int(yCoords.size()) - 1, segment.v);
      prvX = segment.x;
    }

    printf("Case %d: %lld\n", t, ans);
  }
  return 0;
}

// https://lightoj.com/problem/rectangle-union
