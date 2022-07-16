// Dynamic Programming based C++ program to find shortest path with
// exactly k edges
#include <iostream>
#include <climits>
using namespace std;
// O(V3K)


// Define number of vertices in the graph and infinite value
#define V 4
#define INF INT_MAX

// A Dynamic programming based function to find the shortest path from
// u to v with exactly k edges.
int shortestPath(int graph[][V], int u, int v, int k)
{
	// Table to be filled up using DP. The value sp[i][j][e] will store
	// weight of the shortest path from i to j with exactly k edges
	int sp[V][V][k+1];

	// Loop for number of edges from 0 to k
	for (int e = 0; e <= k; e++)
	{
		for (int i = 0; i < V; i++) // for source
		{
			for (int j = 0; j < V; j++) // for destination
			{
				// initialize value
				sp[i][j][e] = INF;

				// from base cases
				if (e == 0 && i == j)
					sp[i][j][e] = 0;
				if (e == 1 && graph[i][j] != INF)
					sp[i][j][e] = graph[i][j];

				//go to adjacent only when number of edges is more than 1
				if (e > 1)
				{
					for (int a = 0; a < V; a++)
					{
						// There should be an edge from i to a and a
						// should not be same as either i or j
						if (graph[i][a] != INF && i != a &&
							j!= a && sp[a][j][e-1] != INF)
						sp[i][j][e] = min(sp[i][j][e], graph[i][a] +
													sp[a][j][e-1]);
					}
				}
		}
		}
	}
	return sp[u][v][k];
}

// driver program to test above function
int main()
{
	/* Let us create the graph shown in above diagram*/
	int graph[V][V] = { {0, 10, 3, 2},
						{INF, 0, INF, 7},
						{INF, INF, 0, 6},
						{INF, INF, INF, 0}
					};
	int u = 0, v = 3, k = 2;
	cout << shortestPath(graph, u, v, k);
	return 0;
}
