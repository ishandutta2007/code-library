#include <bits/stdc++.h>
using namespace std;

vector<pair<int, int>> get_intersecting_rectangle(int xal, int yal, int xar,
                                                  int yar, int xbl, int ybr,
                                                  int x4, int y4) {
  // bottom-left point of intersection rectangle
  int xcl = max(xal, xbl);
  int ycl = max(yal, ybr);

  // top-right point of intersection rectangle
  int xcr = min(xar, x4);
  int ycr = min(yar, y4);

  vector<pair<int, int>> re;
  if (xcl > xcr || ycl > ycr) {
    return re; // no intersection
  }
  re.push_back(make_pair(xcl, ycl));
  re.push_back(make_pair(xcr, ycr));
  return re;
}

int main() {
  // bottom-left and top-right corners of 1st rectangle
  int xal = 0, yal = 0, xar = 10, yar = 8;

  // bottom-left and top-right corners of 2nd rectangle
  int xbl = 2, ybr = 3, x4 = 7, y4 = 9;

  auto inter = get_intersecting_rectangle(xal, yal, xar, yar, xbl, ybr, x4, y4);
  if (inter.size() == 0)
    cout << "No intersection" << endl;
  else if (inter.size() == 2)
    cout << "(" << inter[0].first << "," << inter[0].second << ")"
         << "(" << inter[1].first << "," << inter[1].second << ")" << endl;

  return 0;
}
