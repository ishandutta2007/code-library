#include <bits/stdc++.h>
using namespace std;

const int INF = 1987654321;
const int PHI = -1;
const int NOL = -2;

typedef vector<vector<int>> MatrixInt;
typedef vector<vector<bool>> MatrixBool;

int hungarian_max(MatrixInt &w) {
  int n = w.size();
  int m = w[0].size();
  MatrixBool x(m, vector<bool>(m, false));
  vector<bool> ss(n);
  vector<bool> st(n);
  vector<int> u(m);
  vector<int> v(m);
  vector<int> p(m);
  vector<int> ls(n);
  vector<int> lt(m);
  vector<int> a(n);
  int f = 0;

  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      f = max(f, w[i][j]);

  fill(u.begin(), u.end(), f);
  fill(p.begin(), p.end(), INF);
  fill(lt.begin(), lt.end(), NOL);
  fill(ls.begin(), ls.end(), PHI);
  fill(a.begin(), a.end(), -1);

  while (true) {
    f = -1;
    for (int i = 0; i < n && f == -1; i++) {
      if (ls[i] != NOL && !ss[i])
        f = i;
    }

    if (f != -1) {
      ss[f] = true;
      for (int j = 0; j < m; j++)
        if (!x[f][j] && u[f] + v[j] - w[f][j] < p[j]) {
          lt[j] = f;
          p[j] = u[f] + v[j] - w[f][j];
        }
    } else {
      for (int i = 0; i < m && f == -1; i++)
        if (lt[i] != NOL && !st[i] && p[i] == 0)
          f = i;

      if (f == -1) {
        int d1 = INF, d2 = INF, d;
        for (int i : u)
          d1 = min(d1, i);

        for (int i : p)
          if (i > 0)
            d2 = min(d2, i);

        d = min(d1, d2);

        for (int i = 0; i < n; i++)
          if (ls[i] != NOL)
            u[i] -= d;

        for (int i = 0; i < m; i++) {
          if (p[i] == 0)
            v[i] += d;
          if (p[i] > 0 && lt[i] != NOL)
            p[i] -= d;
        }

        if (d2 >= d1)
          break;
      } else {
        st[f] = true;
        int s = -1;

        for (int i = 0; i < n && s == -1; i++)
          if (x[i][f])
            s = i;

        if (s == -1) {
          for (int l, r;; f = r) {
            r = f;
            l = lt[r];

            if (r >= 0 && l >= 0)
              x[l][r] = !x[l][r];
            else
              break;

            r = ls[l];
            if (r >= 0 && l >= 0)
              x[l][r] = !x[l][r];
            else
              break;
          }

          fill(p.begin(), p.end(), INF);
          fill(lt.begin(), lt.end(), NOL);
          fill(ls.begin(), ls.end(), NOL);
          fill(ss.begin(), ss.end(), false);
          fill(st.begin(), st.end(), false);

          for (int i = 0; i < n; i++) {
            bool ex = true;
            for (int j = 0; j < m && ex; j++)
              ex = !x[i][j];
            if (ex)
              ls[i] = PHI;
          }
        } else
          ls[s] = f;
      }
    }
  }

  int cost = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      if (x[i][j]) {
        a[j] = i;
        cost += w[i][j];
      }
  return cost;
}

int main() {
  int N;
  scanf("%d", &N);
  MatrixInt w(N, vector<int>(N, 0));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      scanf("%d", &w[i][j]);
    }
  cout << hungarian_max(w) << endl;
}
