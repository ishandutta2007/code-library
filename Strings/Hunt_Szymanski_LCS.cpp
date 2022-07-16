#include <bits/stdtr1c++.h>

#define MAX 50010
#define clr(ar) memset(ar, 0, sizeof(ar))
#define read() freopen("lol.txt", "r", stdin)
#define dbg(x) cout << #x << " = " << x << endl

using namespace std;

int ar[MAX];
char A[MAX], B[MAX];

/// Hunt-Szymanski Algorithm for LCS
/// O(R + N) log N, R = numbered of ordered pairs of positions where the two strings match (worst case, R = N^2)

int lcs(char* A, char* B){ /// hash = 935751
    vector <int> adj[256];
    int i, j, l = 0, n = strlen(A), m = strlen(B);
    for (i = 0; i < m; i++) adj[B[i]].push_back(i);

    ar[l++] = -1;
    for (i = 0; i < n; i++){
        for (j = adj[A[i]].size() - 1; j >= 0; j--){
            int x = adj[A[i]][j];
            if (x > ar[l - 1]) ar[l++] = x;
            else ar[lower_bound(ar, ar + l, x) - ar] = x;
        }
    }
    return l - 1;
}

int main(){
    int i, j, k, n, m;
    n = MAX - 10, m = MAX - 10;
    for (i = 0; i < n; i++) A[i] = (rand() % 26) + 'A';
    for (i = 0; i < m; i++) B[i] = (rand() % 26) + 'A';
    A[n] = B[m] = 0;

    clock_t start = clock();
    printf("%d\n", lcs(A, B));

    printf("%0.6f\n", (clock() - start) / (1.0 * CLOCKS_PER_SEC));
    return 0;
}