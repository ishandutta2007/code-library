from math import sqrt

def isqrt(n):
  if n <= 0:
    return 0

  x = int(sqrt(n) * (1 + 1e-12))
  while True:
    y = (x + n // x) >> 1
    if y >= x:
      return x
    x = y

def icbrt(n):
  if n <= 0:
    return 0

  x = int(n ** (1. / 3.) * (1 + 1e-12))
  while True:
    y = (2 * x + n // (x * x)) // 3
    if y >= x:
      return x
    x = y

def tabulate_all_primes_sum(N):
  def T(n):
    return n * (n + 1) // 2 - 1

  if N <= 1:
    return ValueError(N)

  v = isqrt(N)
  smalls = [T(i) for i in range(v + 1)]
  larges = [0 if i == 0 else T(N // i) for i in range(v + 1)]

  for p in range(2, v + 1):
    if smalls[p - 1] == smalls[p]:
      continue
    p_sum = smalls[p - 1]
    q = p * p
    end = min(v, N // q)
    for i in range(1, end + 1):
      d = i * p
      if d <= v:
        larges[i] -= (larges[d] - p_sum) * p
      else:
        larges[i] -= (smalls[N // d] - p_sum) * p
    for i in range(v, q - 1, -1):
      smalls[i] -= (smalls[i // p] - p_sum) * p
  return smalls, larges

N = 100

sqrt_N = isqrt(N)
cbrt_N = icbrt(N)

# about O(N^(3/4))
psum_smalls, psum_larges = tabulate_all_primes_sum(N) 
print(psum_smalls)
print(psum_larges)
for i in range(1, sqrt_N+1):
  print("Prime count till {} = {}".format(N//i,psum_larges[i]))


