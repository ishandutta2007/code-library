// C++ program for above approach
#include <bits/stdc++.h>
using namespace std;

// Graph class represents a undirected graph
// using adjacency list representation
class Graph {
  // No. of vertices
  int V;

  // Pointer to an array containing adjacency lists
  list<int> *adj;

  // A function used by DFS
  void DFSUtil(int v, bool visited[]);

public:
  // Constructor
  Graph(int V);

  void addEdge(int v, int w);
  int NumberOfconnectedComponents();
};

// Function to return the number of
// connected components in an undirected graph
int Graph::NumberOfconnectedComponents() {

  // Mark all the vertices as not visited
  bool *visited = new bool[V];

  // To store the number of connected components
  int count = 0;
  for (int v = 0; v < V; v++)
    visited[v] = false;

  for (int v = 0; v < V; v++) {
    if (visited[v] == false) {
      DFSUtil(v, visited);
      count += 1;
    }
  }

  return count;
}

void Graph::DFSUtil(int v, bool visited[]) {

  // Mark the current node as visited
  visited[v] = true;

  // Recur for all the vertices
  // adjacent to this vertex
  list<int>::iterator i;

  for (i = adj[v].begin(); i != adj[v].end(); ++i)
    if (!visited[*i])
      DFSUtil(*i, visited);
}

Graph::Graph(int V) {
  this->V = V;
  adj = new list<int>[V];
}

// Add an undirected edge
void Graph::addEdge(int v, int w) {
  adj[v].push_back(w);
  adj[w].push_back(v);
}

int32_t main() {
  int n, m;
  cin >> n >> m;
  Graph g(n);
  for (int i = 1; i <= m; i++) {
    int u, v;
    cin >> u >> v;
    g.addEdge(u, v);
  }
  cout << g.NumberOfconnectedComponents() << '\n';
  return 0;
}
