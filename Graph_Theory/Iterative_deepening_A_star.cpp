class Solution {
public:
    bool dfs(int maxd, int n, int i, int d) {
        if (d == maxd) return !n;
        for (int k = i; ; k++) {
            int sq = k * k;
            if ((maxd - d) * sq > n) break; // A*
            if (dfs(maxd, n - sq, k, d + 1)) return true;
        }
        return false;
    }
    
    int numSquares(int n) {
        // IDA
        for (int d = 1; ; d++) {
            if (dfs(d, n, 1, 0))
                return d;
        }
        return -1;
    }
};
//https://leetcode.com/problems/perfect-squares/