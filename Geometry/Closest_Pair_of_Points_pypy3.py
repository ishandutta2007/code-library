# This is faster than fastest.cpp in pypy3
from time import perf_counter
import sys
from math import sqrt

# def closest_pair_distance_v1(n=2000000):
#     sprev = 290797
#     pts = []
#     for _ in range(n):
#         s = sprev*sprev%50515093
#         x, sprev = sprev, s
#         s = sprev*sprev%50515093
#         y, sprev = sprev, s
#         pts.append((x,y))
#     pts.sort()
#     x0, y0 = pts[0]
#     x1, y1 = pts[1]
#     h = sqrt((x1-x0)**2+(y1-y0)**2)
#     for i in range(2, n):
#         x, y = pts[i]
#         for j in range(i-1, -1, -1):
#             x0, y0 = pts[j]
#             if x-x0 > h:
#                 break
#             d = sqrt((x-x0)**2+(y-y0)**2)
#             if d < h:
#                 h = d
#     return f"{h:.9f}"

def closest_pair_distance_v2(n=2000000):
    mod = 50515093
    sprev = 290797
    pts = []
    for _ in range(n):
        s = sprev*sprev%mod
        x, sprev = sprev, s
        s = sprev*sprev%mod
        y, sprev = sprev, s
        pts.append(x*mod+y)
    pts.sort()
    xs, ys = [0]*n, [0]*n
    for i, pt in enumerate(pts):
        xs[i] = pt//mod
        ys[i] = pt%mod
    x0, y0 = xs[0], ys[0]
    x1, y1 = xs[1], ys[1]
    h = sqrt((x1-x0)**2+(y1-y0)**2)
    for i in range(2, n):
        x, y = xs[i], ys[i]
        for j in range(i-1, -1, -1):
            x0, y0 = xs[j], ys[j]
            if x-x0 > h:
                break
            d = sqrt((x-x0)**2+(y-y0)**2)
            if d < h:
                h = d
    return f"{h:.9f}"

start_time = perf_counter()
result = closest_pair_distance_v2(n=2000000)
if result is not None: print(result)
end_time = perf_counter()
print(f"Time: {end_time - start_time} secs")

# if __name__ == "__main__":
#     def main():
#         global start_time
#         entrypoint = closest_pair_distance_v2
#         while len(sys.argv) > 1 and sys.argv[1][0] == '-':
#             opt = sys.argv[1].lower()
#             if opt in ('-v1',):
#                 entrypoint = closest_pair_distance_v1
#             elif opt in ('-v2',):
#                 entrypoint = closest_pair_distance_v2
#             else:
#                 print(f"Unknown option {sys.argv[1]}")
#                 sys.exit(1)
#             sys.argv.pop(1)
#         start_time = perf_counter()
#         try:
#             result = entrypoint(*map(int, sys.argv[1:]))
#         except KeyboardInterrupt:
#             result = None
#             print("\nInterrupted")
#         end_time = perf_counter()
#         if result is not None: print(result)
#         print(f"Time: {end_time-start_time} secs")
#     if len(sys.argv) > 1 and sys.argv[1] in ('-profile', '-prof'):
#         import cProfile
#         sys.argv.pop(1)
#         cProfile.run('main()', sort='cumtime')
#     else:
#         main()


