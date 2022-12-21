# https://www.youtube.com/watch?v=h6MxKmDimhg
def minScoreTriangulation(self, values: List[int]):
    @cache
    def dfs(i, j):
        tmp = float("inf")
        for k in range(i + 1, j):
            tmp = min(tmp, values[i] * values[k] * values[j] + dfs(i, k) + dfs(k, j))

        if tmp == float("inf"):
            return 0
        return tmp

    return dfs(0, len(values) - 1)
