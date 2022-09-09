#include <bits/stdc++.h>
using namespace std;
// https://en.wikipedia.org/wiki/Closest_pair_of_points_problem
// O(nloglogn)
// class Point
// {
// 	public:
// 	int x, y;
// };

bool compareX(const pair<int, int> a, const pair<int, int> b) {
  return (a.first < b.first);
}

bool compareY(const pair<int, int> a, const pair<int, int> b) {
  return (a.second < b.second);
}

float dist(pair<int, int> p1, pair<int, int> p2) {
  return sqrt((p1.first - p2.first) * (p1.first - p2.first) +
              (p1.second - p2.second) * (p1.second - p2.second));
}

float bruteForce(vector<pair<int, int>> P, int st, int end) {
  float min = FLT_MAX;
  for (int i = st; i < end; ++i)
    for (int j = i + 1; j < end; ++j)
      if (dist(P[i], P[j]) < min)
        min = dist(P[i], P[j]);
  return min;
}

float min(float x, float y) { return (x < y) ? x : y; }

float stripClosest(vector<pair<int, int>> strip, int size, float d) {
  float min = d; // Initialize the minimum distance as d
  sort(strip.begin(), strip.end(), compareY);
  for (int i = 0; i < size; ++i)
    for (int j = i + 1; j < size && (strip[j].second - strip[i].second) < min;
         ++j)
      if (dist(strip[i], strip[j]) < min)
        min = dist(strip[i], strip[j]);

  return min;
}

float closestUtil(vector<pair<int, int>> P, int st, int end) {
  // If there are 2 or 3 points, then use brute force
  if (end - st <= 3)
    return bruteForce(P, st, end);

  int mid = (end - st) / 2;
  pair<int, int> midPoint = P[mid];

  float dl = closestUtil(P, st, mid);
  float dr = closestUtil(P, mid, end);

  float d = min(dl, dr);

  vector<pair<int, int>> strip;
  for (int i = st; i < end; i++)
    if (abs(P[i].first - midPoint.first) < d)
      strip.push_back(P[i]);

  return min(d, stripClosest(strip, strip.size(), d));
}

float closest(vector<pair<int, int>> P) {
  sort(P.begin(), P.end(), compareX);
  // for(int i=0;i<P.size();i++)cout<<P[i].first<<" "<<P[i].second<<endl;
  return closestUtil(P, 0, P.size());
}

int main() {
  vector<pair<int, int>> P;
  P.push_back(make_pair(2, 3));
  P.push_back(make_pair(12, 30));
  P.push_back(make_pair(40, 50));
  P.push_back(make_pair(5, 1));
  P.push_back(make_pair(12, 10));
  P.push_back(make_pair(3, 4));
  cout << "The smallest distance is " << closest(P);
  return 0;
}
