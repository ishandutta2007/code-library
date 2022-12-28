// Credits @Gandhali-Shastri

#include <stdio.h>
#include <iostream>
#include <conio.h>
#include "stdafx.h"

// Basic Definitions

#define WHITE 0
#define GRAY 1
#define BLACK 2
#define MAX_NODES 3000
#define oo 1000000000

// Declarations

int n;                              // number of nodes
int e;                              // number of edges
int capacity[MAX_NODES][MAX_NODES]; // capacity matrix
int flow[MAX_NODES][MAX_NODES];     // flow matrix
int color[MAX_NODES];               // needed for breadth-first search
int pred[MAX_NODES];                // array to store augmenting path
int cut_tree[MAX_NODES][MAX_NODES];

int min(int x, int y) {
  return x < y ? x : y; // returns minimum of x and y
}
// A Queue for Breadth-First Search

int head, tail;
int q[MAX_NODES + 2];

void enqueue(int x) {
  q[tail] = x;
  tail++;
  color[x] = GRAY;
}

int dequeue() {
  int x = q[head];
  head++;
  color[x] = BLACK;
  return x;
}

// Breadth-First Search for an augmenting path

int bfs(int start, int target) {
  int u, v;
  for (u = 0; u < n; u++) {
    color[u] = WHITE;
  }
  head = tail = 0;
  enqueue(start);
  pred[start] = -1;
  while (head != tail) {
    u = dequeue();
    // Search all adjacent white nodes v. If the capacity
    // from u to v in the residual network is positive,
    // enqueue v.
    for (v = 0; v < n; v++) {
      if (color[v] == WHITE && capacity[u][v] - flow[u][v] > 0) {
        enqueue(v);
        pred[v] = u;
      }
    }
  }

  // If the color of the target node is black now,
  // it means that we reached it.
  return color[target] == BLACK;
}

// Ford-Fulkerson Algorithm

int max_flow(int source, int sink) {
  int i, j, u;
  // Initialize empty flow.
  int max_flow = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      flow[i][j] = 0;
    }
  }
  // While there exists an augmenting path,
  // increment the flow along this path.
  while (bfs(source, sink)) {
    // Determine the amount by which we can increment the flow.
    int increment = oo;
    for (u = sink; pred[u] != (-1); u = pred[u]) {
      increment = min(increment, capacity[pred[u]][u] - flow[pred[u]][u]);
    }
    // Now increment the flow.
    for (u = sink; pred[u] != (-1); u = pred[u]) {
      flow[pred[u]][u] += increment;
      flow[u][pred[u]] -= increment; // Reverse in residual
    }
    max_flow += increment;
  }
  // No augmenting path anymore. We are done.
  return max_flow;
}

void gomoryhu() {
  int i, j, s, t, pos, maximumFlow;
  int p[MAX_NODES], f1[MAX_NODES], visited[MAX_NODES];

  for (i = 0; i < n; i++) {
    p[i] = 0;
    f1[i] = 0;
    visited[i] = 0;
    for (j = 0; j < n; j++) {
      cut_tree[i][j] = 0;
    }
  }

  for (s = 1; s < n; s++) {
    t = p[s];
    maximumFlow = max_flow(s, t);
    f1[s] = maximumFlow;

    for (i = 0; i < n; i++) {
      visited[i] = 0;
    }

    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        cut_tree[i][j] = 0;
      }
    }

    for (i = 0; i < n; i++) {
      if (color[i] == BLACK) {
        visited[i] = BLACK;
      }
    }

    for (i = 0; i < n; i++) {
      if (i != s && p[i] == t && visited[i] == BLACK) {
        p[i] = s;
      }
    }

    if (visited[p[t]] == BLACK) {
      p[s] = p[t];
      p[t] = s;
      f1[s] = f1[t];
      f1[t] = maximumFlow;
    }

    for (i = 1; i <= s; i++) {
      cut_tree[i][p[i]] = f1[i];
    }
  }
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      if (cut_tree[i][j] > 0)
        printf("%d %d %d\n", i, j, cut_tree[i][j]);
}

int main() {

  int a, b, c, i, j;
  // FILE* input = fopen("mf.in", "r");
  // read number of nodes and edges
  scanf_s("%d %d", &n, &e);
  // initialize empty capacity matrix
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      capacity[i][j] = 0;
    }
  }
  // read edge capacities
  for (i = 0; i < e; i++) {
    scanf_s("%d %d %d", &a, &b, &c);
    capacity[a][b] = c;
    capacity[b][a] = c; // Could have parallel edges
  }

  // fclose(input);
  gomoryhu();

  /*cin.clear();
  cin.ignore();
  cin.get();*/

  return 0;
}

// https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.flow.gomory_hu_tree.html
// http://www.spoj.com/problems/MCQUERY