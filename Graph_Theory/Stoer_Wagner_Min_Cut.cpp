/*
 * Practical Work 2
 * Operational Research 2019/1 - UFMG
 *
 * Francisco Galuppo Azevedo
 * 2017014960
 *
 * Uses the Stoer-Wagner algorithm:
 * A Simple Min-Cut Algorithm, by MECHTHILD STOER and FRANK WAGNER
*/

#include <iostream>
#include <vector>
#include <utility>
#include <numeric>

#define INF 2147483647

using namespace std;

// Global variables
int nv;
vector<bool> deleted;
vector<vector<int>> equivalent;

// st cut
pair<int, int> minimumCutPhase(vector<vector<int>> &g, int a) {
  int n = g.size();

  // Create vectors
  vector<bool> A(n, false);
  A[a] = true;

  vector<int> w(n, 0);

  // Calculate initial weights of vertices
  for (int u = 0; u < n; u++)
    if (!deleted[u])
      w[u] = g[a][u];

  // While A <> V
  int prev = a;
  for (int i = 1; i < nv; i++) {
    // Find the most "fair" vertex
    int z = -1;
    for (int j = 0; j < n; j++)
      if (!deleted[j] && !A[j] && (z < 0 || w[j] > w[z]))
        z = j;

    A[z] = true;

    // If necessary, merge z and prev
    if (i == nv - 1) {
      // Update edge weights
      for (int j = 0; j < n; j++)
        g[j][prev] = g[prev][j] += g[z][j];

      // Mark prev as prev U z
      deleted[z] = true;
      for (auto u : equivalent[z])
        equivalent[prev].push_back(u);

      nv--;
      return make_pair(z, w[z]);
    }

    prev = z;

    // Update weights
    for (int j = 0; j < n; j++)
      if (!A[j] && !deletado[j])
        w[j] += g[z][j];
  }
}

// Mínimo st corte
pair<vector<int>, int> minimumCut(vector<vector<int>> g, int a) {
  //  Values
  int min = INF;
  vector<int> cut;
  pair<int, int> val;

  // Finds the minimum among all st-cuts
  while (nv > 1) {
    val = minimumCutPhase(g, a);
    if (val.second < min) {
      min = val.second;
      cut = equivalente[val.first];
    }
  }

  return make_pair(cut, min);
}

int main() {
  ios_base::sync_with_stdio(false);

  // Read in the dimensions (number of vertices and number of edges) of the
  // graph
  int n, m;
  cin >> n >> m;
  nv = n;

  // Read in the edges of the graph
  vector<vector<int>> g(n, vector<int>(n, 0));

  {
    int u, v, w;
    for (int i = 0; i < m; i++) {
      cin >> u >> v >> w;
      g[u][v] = g[v][u] = w;
    }
  }

  // Vertices remaining
  vector<int> v(n);
  iota(v.begin(), v.end(), 0);

  // Reseta o vetor de deletados
  deletado.resize(n, false);
  equivalente.resize(n);
  for (int i = 0; i < n; i++)
    equivalente[i].push_back(i);

  // Calcula o corte e escreve saída
  pair<vector<int>, int> cut = minimumCut(g, 0);
  cout << cut.first.size() << endl;
  for (auto u : cut.first)
    cout << u << " ";
  cout << endl << cut.second << endl;

  return 0;
}
