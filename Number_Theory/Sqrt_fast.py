# This is 50x faster than math.sqrt()
# https://stackoverflow.com/questions/53698655/python-looking-for-a-faster-less-accurate-sqrt-function
def sqrt_approx(x):
    delta = 0.1
    runner = x / 2
    while abs(runner - (x / runner)) > delta:
        runner = ((x / runner) + runner) / 2
    return runner
  
