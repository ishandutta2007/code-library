def PPT():
    "Primitive Pythagorean Triples"
    from collections import deque

    fifo = deque([(3, 4, 5)])
    while True:
        (a, b, c) = fifo.popleft()
        yield (a, b, c)
        fifo.append((a - 2 * b + 2 * c, 2 * a - b + 2 * c, 2 * a - 2 * b + 3 * c))
        fifo.append((a + 2 * b + 2 * c, 2 * a + b + 2 * c, 2 * a + 2 * b + 3 * c))
        fifo.append((-a + 2 * b + 2 * c, -2 * a + b + 2 * c, -2 * a + 2 * b + 3 * c))
