#include <iostream>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <functional>
#include <cassert>
using namespace std;

#define fst first
#define snd second

// credits: Jatin Uadav
// 1-indexing has been used.
// vertices are numbered 1, ..., n
// weight_type is the type of edge weights (int, float, long long, double etc.)

template<class weight_type>
struct tree{
	int n; // number of vertices
	vector<vector<pair<int, weight_type>>> con; // adjacency list

	tree(int n): n(n), con(n + 1){}

	void add_edge(int node1, int node2, weight_type W){
		con[node1].push_back({node2, W});
		con[node2].push_back({node1, W});
	}

	/*
	Input :
		C - maximum capacity of any center.
		k - number of centers allowed

	Output :
		A tuple (minimum radius, set of centers opened, assignment of each vertex to some center)

	Algorithm Description :
		First get a sorted list of distances in O(n ^2 log(n))
		Binary search on the index in this list. This works as the answer is going to be
		distance between some pair.

		To find a feasible assignment for a given radius, we use the algorithm in the paper : 
			Technical Note—A Polynomial Algorithm for the Equal
			Capacity p-Center Problem on Trees

	Time Complexity : 
		O(n^2 log(n)) [For sorting the distances]

	Space Complexity : 
		O(n^2) [For storing pairwise distances.]
	*/
	tuple<weight_type, vector<int>, vector<int> > k_centers(int C, int k){
		assert(k * C >= n); // else no solution exists.

		vector<int> centers;
		vector<int> assignment(n + 1);

		vector<weight_type> distances;
		vector<vector<int>> D(n + 1);
		vector<vector<int>> ancestors(n + 1);

		// choose a root of the tree.
		// can be any vertex. Chosen 1 here.
		int root = 1;

		// The answer is the distance between some two vertices.
		// Precompute the distances. O(n^2) in total.
		function<void(int, int, int, int)> dfs = [&](int beg, int s, int p, weight_type d){
			D[beg][s] = d;

			if(beg == root){
				if(p != -1) ancestors[s] = ancestors[p];
				ancestors[s].push_back(s);
			}

			if(beg < s) distances.push_back(d);

			for(auto v : con[s]) if(v.fst != p){
				dfs(beg, v.fst, s, d + v.snd);
			}
		};
			
		for(int i = 1; i <= n; i++){
			D[i].resize(n + 1);
			dfs(i, i, -1, 0);
		}
		sort(distances.begin(), distances.end());

		// check if a given radius can work.
		function <bool(weight_type)> check = [&] (weight_type radius){
			// f(i) = nearest to root vertex, that can cover i
			// U = set of unassigned vertices

			set<int> U;
			vector<int> f(n);
			set<pair<weight_type, int> >st;
			
			for(int i = 1; i <= n; i++){
				U.insert(i);
				weight_type d_i = D[root][i];
				for(auto it : ancestors[i]) if(d_i - D[root][it] <= radius){
					f[i] = it;
					st.insert({-D[root][it], i});
					break;
				}
			}

			centers.clear();

			for(int iter = 0; !U.empty() && iter < k; iter++){
				int center = st.begin()->snd;
				centers.push_back(center);
				st.erase(st.begin());
				vector<int> reachable;
				for(auto it : U){
					if(D[center][it] <= radius) reachable.push_back(it);
				}

				// sort the set reachable in decreasing order of distance from root.
				sort(reachable.begin(), reachable.end(), [&](int i, int j){return D[root][i] > D[root][j];});
				if(reachable.size() > C) reachable.resize(C);

				for(auto it : reachable){
					assignment[it] = center;
					U.erase(it);
				}
			}
			return U.empty();
		};

		int lo = 0, hi = distances.size() - 1;

		// binary search on the index of the answer in the set of distances.
		while(lo < hi){
			int mid = (lo + hi) >> 1;
			if(check(distances[mid])) hi = mid;
			else lo = mid + 1;
		}
		check(distances[lo]);
		return make_tuple(distances[lo], centers, assignment);
	}
};

void usage_example(){
	tree<int> T(4);

//		       	1
//		      / | \
//		   1 / 3|  \ 1
//		    /   |   \
//		   2    3    4

	T.add_edge(1, 2, 1);
	T.add_edge(1, 3, 3);
	T.add_edge(1, 4, 1);

	int radius; vector<int> centers, assignment;

	tie(radius, centers, assignment) = T.k_centers(3, 2);

	// minimum radius
	cout << "radius : "; cout << radius << endl;	

	// centers opened
	cout << "centers : "; for(auto it : centers) cout << it << " "; cout << endl;
	// assignment
	cout << "assignment : ";
	for(int i = 1; i <= T.n; i++) cerr << assignment[i] << " "; cerr << endl;
}

int main(){
	usage_example();	
}
