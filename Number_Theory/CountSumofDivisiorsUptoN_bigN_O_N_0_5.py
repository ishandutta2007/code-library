import math

def A006218(n):
    return 2*sum(n//k for k in range(1,  math.isqrt(n)+1))- math.isqrt(n)**2

n = 10
print(A006218(n))
