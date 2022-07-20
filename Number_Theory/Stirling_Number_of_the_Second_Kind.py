def sterling2(n, k):
    key = str(n) + "," + str(k)
 
    if key in computed.keys():
        return computed[key]
    if n == k == 0:
        return 1
    if (n > 0 and k == 0) or (n == 0 and k > 0):
        return 0
    if n == k:
        return 1
    if k > n:
        return 0
    result = k * sterling2(n - 1, k) + sterling2(n - 1, k - 1)
    computed[key] = result
    return result